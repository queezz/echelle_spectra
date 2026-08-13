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


def _balmer_lines() -> tuple[SpectralLine, ...]:
    source = "Echelle wavelength-validation convention"
    reference = "docs/line-validation.md"
    values = (
        ("H-alpha", 656.2790),
        ("H-beta", 486.1350),
        ("H-gamma", 434.0470),
        ("H-delta", 410.1734),
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
                )
                for label, wavelength in values
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
    result = []
    for band, wavelengths in payload["bands"].items():
        for number, wavelength in enumerate(wavelengths, start=1):
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
                    notes=notes,
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
def load_line_table(family: str) -> tuple[SpectralLine, ...]:
    """Return one immutable, wavelength-sorted family table.

    Parameters
    ----------
    family:
        One of ``balmer``, ``fulcher``, ``thar``, ``ne``, or ``hg``.
    """

    key = family.strip().lower()
    if key not in LINE_FAMILIES:
        known = ", ".join(LINE_FAMILIES)
        raise ValueError(f"unknown line family {family!r}; known families: {known}")
    if key == "balmer":
        return _balmer_lines()
    if key == "fulcher":
        return _fulcher_lines()
    return _nist_lines(key)


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
