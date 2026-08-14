"""Shared, provenance-carrying spectral line tables.

The catalog is calibration-facing.  It provides line positions and lightweight
identification metadata without importing downstream molecular fitting logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Iterable

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback.
    import tomli as tomllib  # type: ignore[no-redef]

from .nist_lamp_calibration import (
    COMMON_NIST_SPECIES,
    default_nist_cache_dir,
    load_nist_asd_exports,
    resolve_cached_line_lists,
)

__all__ = [
    "LINE_FAMILIES",
    "LINE_FAMILY_ISOTOPES",
    "LINE_FAMILY_LABELS",
    "SpectralLine",
    "filter_line_table",
    "load_line_table",
]


LINE_FAMILIES = ("balmer", "fulcher", "thar", "ne", "hg")
LINE_FAMILY_LABELS = {
    "balmer": "Balmer",
    "fulcher": "Fulcher H2",
    "thar": "ThAr",
    "ne": "Ne",
    "hg": "Hg",
}

#: Which hydrogen isotopologues each family can be asked for.  Only the two
#: hydrogen families have an isotopologue at all; a lamp family carries none, so
#: naming one for ThAr, Ne, or Hg is a question the catalog answers with an
#: empty table rather than with the table that was not asked for.
#:
#: The Fulcher entry says ``("H",)`` because the bundled Q-branch anchors are
#: H2 and nothing else.  There is no D2 table on this side; when one arrives it
#: is added here and read by ``_fulcher_lines``, and every caller that already
#: asks for ``isotope="D"`` starts receiving it without changing a line.
LINE_FAMILY_ISOTOPES: dict[str, tuple[str, ...]] = {
    "balmer": ("H", "D"),
    "fulcher": ("H",),
    "thar": (),
    "ne": (),
    "hg": (),
}

#: The isotopologue a family returns when a caller names none.  Every table
#: bundled before deuterium existed was a hydrogen table, so this default is
#: what keeps those callers reading exactly what they always read.
DEFAULT_ISOTOPE = "H"


@dataclass(frozen=True)
class SpectralLine:
    """One identified line in normalized nanometre units."""

    family: str
    label: str
    wavelength_nm: float
    wavelength_medium: str
    species: str
    source_name: str
    source_reference: str
    source_resource: str
    relative_intensity: float | None = None
    notes: str = ""
    blended: bool = False
    """True when another catalog line falls inside the instrumental width.

    A blended row is still a real transition and is still drawn by overlays,
    but a single measured centroid over it belongs to the blend rather than to
    this line, so calibration measurements must skip it.
    """

    isotope: str = ""
    """Hydrogen isotopologue this row belongs to: ``H``, ``D``, or empty.

    Empty means the question does not apply, which is the case for every lamp
    line: a Ne line is a Ne line.
    """

    transition: str = ""
    """Upper and lower level of the transition, when isotopologues are paired.

    The Balmer rows carry ``3-2``/``4-2``/``5-2``/``6-2`` so that H-alpha and
    D-alpha can be recognised as one transition seen through two nuclei.  A
    consumer that must judge a measured centroid against both references pairs
    them on this, never on the text of the label.
    """


def _balmer_lines() -> tuple[SpectralLine, ...]:
    source = "Echelle wavelength-validation convention"
    reference = "docs/line-validation.md"
    values = (
        ("3-2", "H-alpha", 656.2790),
        ("4-2", "H-beta", 486.1350),
        ("5-2", "H-gamma", 434.0470),
        ("6-2", "H-delta", 410.1734),
    )
    return tuple(
        sorted(
            (
                SpectralLine(
                    family="balmer",
                    label=label,
                    wavelength_nm=wavelength,
                    wavelength_medium="air",
                    species="H I",
                    source_name=source,
                    source_reference=reference,
                    source_resource="echelle_spectra.tools.emissiondata",
                    isotope="H",
                    transition=transition,
                )
                for transition, label, wavelength in values
            ),
            key=lambda line: line.wavelength_nm,
        )
    )


def _deuterium_balmer_lines() -> tuple[SpectralLine, ...]:
    """Return the D I counterparts of the Balmer rows above, with their source."""

    resource = files("echelle_spectra.resources").joinpath(
        "line_catalogs/balmer_deuterium.toml"
    )
    payload = tomllib.loads(resource.read_text(encoding="utf-8"))
    notes = f"{payload['adaptation_note']} {payload['citation']}"
    return tuple(
        sorted(
            (
                SpectralLine(
                    family="balmer",
                    label=str(row["label"]),
                    wavelength_nm=float(row["wavelength_nm"]),
                    wavelength_medium=str(payload["wavelength_medium"]),
                    species="D I",
                    source_name=str(payload["source_name"]),
                    source_reference=str(payload["source_reference"]),
                    source_resource=str(payload["source_resource"]),
                    notes=notes,
                    isotope="D",
                    transition=str(row["transition"]),
                )
                for row in payload["lines"]
            ),
            key=lambda line: line.wavelength_nm,
        )
    )


def _fulcher_lines() -> tuple[SpectralLine, ...]:
    resource = files("echelle_spectra.resources").joinpath(
        "line_catalogs/fulcher_h2_q.toml"
    )
    payload = tomllib.loads(resource.read_text(encoding="utf-8"))
    notes = str(payload["adaptation_note"])
    blended = {float(value) for value in payload.get("blended_wavelengths", ())}
    blend_note = str(payload.get("blend_note", ""))
    result = []
    for band, wavelengths in payload["bands"].items():
        for number, wavelength in enumerate(wavelengths, start=1):
            is_blended = float(wavelength) in blended
            result.append(
                SpectralLine(
                    family="fulcher",
                    label=f"Q{number}({band})",
                    wavelength_nm=float(wavelength),
                    wavelength_medium=str(payload["wavelength_medium"]),
                    species="H2",
                    source_name=str(payload["source_name"]),
                    source_reference=str(payload["source_reference"]),
                    source_resource=str(payload["source_resource"]),
                    notes=f"{notes} {blend_note}".strip() if is_blended else notes,
                    blended=is_blended,
                    isotope="H",
                )
            )
    return tuple(sorted(result, key=lambda line: (line.wavelength_nm, line.label)))


def _nist_lines(family: str) -> tuple[SpectralLine, ...]:
    exports = resolve_cached_line_lists(
        lamps=(family,),
        cache_dir=default_nist_cache_dir(),
    )
    table = load_nist_asd_exports(
        exports,
        min_wavelength_nm=0.0,
        max_wavelength_nm=10_000.0,
    )
    table = table.drop_duplicates(subset=["species", "wavelength_nm"])
    result = []
    for row in table.itertuples(index=False):
        species = str(row.species)
        species_label = COMMON_NIST_SPECIES[species].nist_name
        source_file = str(row.source_path).replace("\\", "/").rsplit("/", 1)[-1]
        result.append(
            SpectralLine(
                family=family,
                label=f"{species_label} {float(row.wavelength_nm):.4f}",
                wavelength_nm=float(row.wavelength_nm),
                wavelength_medium="air",
                species=species_label,
                source_name="NIST Atomic Spectra Database (ASD)",
                source_reference="https://physics.nist.gov/asd",
                source_resource=f"echelle_spectra/resources/nist_asd_cache/{source_file}",
                relative_intensity=float(row.weight),
                notes="Package-cached NIST ASD export; observed air wavelength preferred over Ritz.",
            )
        )
    return tuple(sorted(result, key=lambda line: (line.wavelength_nm, line.species)))


@lru_cache(maxsize=None)
def load_line_table(family: str, *, isotope: str | None = None) -> tuple[SpectralLine, ...]:
    """Return one immutable, wavelength-sorted family table.

    Parameters
    ----------
    family:
        One of ``balmer``, ``fulcher``, ``thar``, ``ne``, or ``hg``.
    isotope:
        Optional hydrogen isotopologue facet, ``"H"`` or ``"D"``.  Omit it to
        receive the family's established table, which is the hydrogen one.  A
        family that holds no table for the named isotopologue returns an empty
        tuple: the honest answer for D2 Fulcher anchors we do not have, and a
        far better one than handing back the H2 table under a D label.
    """

    key = family.strip().lower()
    if key not in LINE_FAMILIES:
        known = ", ".join(LINE_FAMILIES)
        raise ValueError(f"unknown line family {family!r}; known families: {known}")
    available = LINE_FAMILY_ISOTOPES[key]
    if not available:
        if isotope is not None:
            raise ValueError(
                f"{key} lines carry no hydrogen isotopologue, so isotope={isotope!r} "
                "has no meaning for them"
            )
        return _nist_lines(key)
    wanted = (isotope if isotope is not None else DEFAULT_ISOTOPE).strip().upper()
    if wanted not in available:
        return ()
    if key == "balmer":
        return _balmer_lines() if wanted == "H" else _deuterium_balmer_lines()
    return _fulcher_lines()


def filter_line_table(
    lines: Iterable[SpectralLine],
    *,
    minimum_nm: float | None = None,
    maximum_nm: float | None = None,
    minimum_relative_intensity: float | None = None,
) -> tuple[SpectralLine, ...]:
    """Filter a line table without discarding its provenance fields."""

    selected = []
    for line in lines:
        if minimum_nm is not None and line.wavelength_nm < minimum_nm:
            continue
        if maximum_nm is not None and line.wavelength_nm > maximum_nm:
            continue
        if minimum_relative_intensity is not None:
            strength = line.relative_intensity
            if strength is None or strength < minimum_relative_intensity:
                continue
        selected.append(line)
    return tuple(selected)
