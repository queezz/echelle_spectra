"""Shared, provenance-carrying spectral line tables.

The catalog is calibration-facing.  It provides line positions and lightweight
identification metadata without importing downstream molecular fitting logic.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from importlib.resources import as_file, files
from typing import Iterable

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback.
    import tomli as tomllib  # type: ignore[no-redef]

from .calibration_alignment import BH_PAPER_WAVELENGTH_TABLE, load_wavelength_table
from .nist_lamp_calibration import (
    COMMON_LAMP_PRESETS,
    COMMON_NIST_SPECIES,
    default_nist_cache_dir,
    load_nist_asd_exports,
    normalize_species_key,
    resolve_cached_line_lists,
)

__all__ = [
    "CATALOG_LINE_FAMILIES",
    "CURATED_LINE_TABLE",
    "CURATED_MATCH_TOLERANCE_NM",
    "LAMP_LINE_FAMILIES",
    "LINE_FAMILIES",
    "LINE_FAMILY_ISOTOPES",
    "LINE_FAMILY_LABELS",
    "SpectralLine",
    "curated_line_count",
    "filter_line_table",
    "load_line_table",
]


#: The families the *overlays* draw: one checkbox and one hue each, in both the
#: 1-D and 2-D views.  It is a closed set on purpose — the palette in
#: :data:`~.line_overlay.LINE_OVERLAY_STYLES` seats exactly these five in the
#: cool arc, and a sixth cannot be added without redesigning it.
LINE_FAMILIES = ("balmer", "fulcher", "thar", "ne", "hg")

#: Every family :func:`load_line_table` will serve, which is a wider set.
#:
#: ``xe`` is here and not in :data:`LINE_FAMILIES` because the two tuples answer
#: different questions.  A xenon lamp on the calibration bench needs its lines —
#: the expected-lines panel, the fit reference, and the identification help all
#: read this catalog — and it needs them whether or not the main viewer has a
#: sixth overlay toggle to spare.  Which lamps the bench can serve is a
#: question about physics and cached data; how many hues fit between the red and
#: yellow-green spectrum curves is a question about one window's palette, and
#: letting the second answer the first is what would have kept xenon out.
CATALOG_LINE_FAMILIES = (*LINE_FAMILIES, "xe")

#: The families whose tables are NIST lamp exports rather than a hand-kept list.
#: These are the ones a curated wavelength table can vouch for, and the ones the
#: overlays filter by strength.
#:
#: ``xe`` is one of them and is vouched for by nothing: the curated table was
#: written from ThAr, Ne, Hg and H2 frames and carries no xenon row at all.
#: That is not a defect in the family — it is the honest state of a lamp whose
#: cache the package ships and whose vetting nobody has done, and every surface
#: that reads :func:`load_line_table` says so rather than inventing a pedigree.
LAMP_LINE_FAMILIES = ("thar", "ne", "hg", "xe")
LINE_FAMILY_LABELS = {
    "balmer": "Balmer",
    "fulcher": "Fulcher H2",
    "thar": "ThAr",
    "ne": "Ne",
    "hg": "Hg",
    "xe": "Xe",
}

#: Which hydrogen isotopologues each family can be asked for.  Only the two
#: hydrogen families have an isotopologue at all; a lamp family carries none, so
#: naming one for ThAr, Ne, Hg, or Xe is a question the catalog answers with an
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
    "xe": (),
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
    """Lamp-context strength in ``(0, 1]``, or ``None`` for a table without one.

    For the NIST-backed lamp families this is the cached relative intensity as
    a fraction of that spectrum's strongest cached line, scaled by the
    ionization stage's
    :data:`~echelle_spectra.tools.nist_lamp_calibration.LAMP_STAGE_WEIGHTS`.
    The stage scaling is what makes the number comparable **between** species
    of one lamp: NIST prints each spectrum on its own source's scale, so
    without it the strongest Ne II line and the strongest Ne I line both read
    1.0 and a selection ranked on strength prefers ions a neon lamp barely
    excites.

    It is a selection aid, not a photometric quantity: do not read a line ratio
    off two of these.
    """

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

    curated: bool = False
    """True when the vetted wavelength table names this line as one of its own.

    The curated table is :data:`CURATED_LINE_TABLE`, whose ``OK`` marks carry
    the BH paper's vetting: somebody found this line on this instrument, fitted
    it, and signed for it.  That is a stronger statement about "is this line
    here" than any database strength number, which is why the overlays never
    let a curated row lose a selection contest to a NIST threshold — see
    :func:`~echelle_spectra.tools.line_overlay.select_overlay_lines`.
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


#: The wavelength table whose rows the packaged catalogs treat as curated.
#: It is the BH-paper table itself rather than one of its aligned descendants:
#: an alignment moves pixels, never wavelengths, so the *identifications* — the
#: only thing read here — are established once, at the head of the lineage.
CURATED_LINE_TABLE = BH_PAPER_WAVELENGTH_TABLE

#: How near a cached NIST row must sit to a curated row to be the same line.
#: The curated wavelengths were transcribed from ASD to five decimals, so a
#: genuine pair agrees far inside this; 0.005 nm is loose enough for a rounded
#: transcription and tight enough to keep the 667.8277 / 667.8331 Ne I pair
#: apart, where the curated row means the strong one.
CURATED_MATCH_TOLERANCE_NM = 0.005


@lru_cache(maxsize=None)
def _curated_wavelengths(family: str) -> tuple[tuple[str, float], ...]:
    """Return ``(species, wavelength_nm)`` for the curated rows of one lamp.

    Only rows the table marks ``OK`` are read — the same predicate
    :func:`~echelle_spectra.tools.calibration_alignment.select_candidate_lines`
    uses — because it is the ``OK`` marks that carry the vetting.  A row left
    at ``?`` is a question somebody wrote down, not an answer, and it earns no
    more standing here than any other database line.  Rows whose species the
    NIST species table does not know (blends, the H2 anchors) are left to the
    families that own them.
    """

    species_keys = COMMON_LAMP_PRESETS[family].species
    resource = files("echelle_spectra.resources").joinpath(
        f"calibration_files/{CURATED_LINE_TABLE}"
    )
    with as_file(resource) as path:
        rows = load_wavelength_table(path)
    curated = []
    for row in rows:
        if row.species not in species_keys:
            continue
        if "ok" not in row.comment.lower():
            continue
        try:
            name = COMMON_NIST_SPECIES[normalize_species_key(row.species)].nist_name
        except KeyError:  # pragma: no cover - defended, not reachable today
            continue
        curated.append((name, float(row.wavelength_nm)))
    return tuple(sorted(set(curated), key=lambda item: (item[1], item[0])))


def _curated_row_line(species: str, wavelength_nm: float, family: str) -> SpectralLine:
    """A curated row the packaged NIST cache holds no counterpart for.

    It is still a line somebody measured on this instrument and signed for, so
    it is carried under the curated table's own provenance rather than dropped.
    It has no ``relative_intensity``: the cache is where that number lives, and
    inventing one would be the database claim this row precisely is not.
    """

    return SpectralLine(
        family=family,
        label=f"{species} {wavelength_nm:.4f}",
        wavelength_nm=wavelength_nm,
        wavelength_medium="air",
        species=species,
        source_name="Echelle curated wavelength table (BH paper vetting)",
        source_reference=f"echelle_spectra/resources/calibration_files/{CURATED_LINE_TABLE}",
        source_resource=f"echelle_spectra/resources/calibration_files/{CURATED_LINE_TABLE}",
        notes="Curated, OK-marked row with no counterpart in the packaged NIST cache.",
        curated=True,
    )


def _with_curated_rows(family: str, table: tuple[SpectralLine, ...]):
    """Flag every cached row a curated row vouches for, and add the rest.

    Each curated wavelength claims the *nearest* cached row of its own species
    inside :data:`CURATED_MATCH_TOLERANCE_NM`.  Two curated rows are allowed to
    claim the same cached row, because the curated table really does name one
    line twice: an echelle overlap exposes it on two orders, and it is measured
    and written down once per order — Ne I 630.47893 and 630.47800 are one
    line, and Hg I 708.19000 and 708.19010 are another.  Refusing the second
    claim would invent a phantom line beside the real one.
    """

    marked = dict(enumerate(table))
    extra: list[SpectralLine] = []
    for species, wavelength in _curated_wavelengths(family):
        candidates = [
            (abs(line.wavelength_nm - wavelength), index)
            for index, line in marked.items()
            if line.species == species
            and abs(line.wavelength_nm - wavelength) <= CURATED_MATCH_TOLERANCE_NM
        ]
        if not candidates:
            extra.append(_curated_row_line(species, wavelength, family))
            continue
        _distance, index = min(candidates)
        marked[index] = replace(marked[index], curated=True)
    result = list(marked.values()) + extra
    return tuple(sorted(result, key=lambda line: (line.wavelength_nm, line.species)))


def curated_line_count(family: str) -> int:
    """How many rows of one family the curated table vouches for, for QA."""

    return sum(1 for line in load_line_table(family) if line.curated)


@lru_cache(maxsize=None)
def load_line_table(family: str, *, isotope: str | None = None) -> tuple[SpectralLine, ...]:
    """Return one immutable, wavelength-sorted family table.

    Parameters
    ----------
    family:
        One of ``balmer``, ``fulcher``, ``thar``, ``ne``, ``hg``, or ``xe``.
    isotope:
        Optional hydrogen isotopologue facet, ``"H"`` or ``"D"``.  Omit it to
        receive the family's established table, which is the hydrogen one.  A
        family that holds no table for the named isotopologue returns an empty
        tuple: the honest answer for D2 Fulcher anchors we do not have, and a
        far better one than handing back the H2 table under a D label.
    """

    key = family.strip().lower()
    if key not in CATALOG_LINE_FAMILIES:
        known = ", ".join(CATALOG_LINE_FAMILIES)
        raise ValueError(f"unknown line family {family!r}; known families: {known}")
    available = LINE_FAMILY_ISOTOPES[key]
    if not available:
        if isotope is not None:
            raise ValueError(
                f"{key} lines carry no hydrogen isotopologue, so isotope={isotope!r} "
                "has no meaning for them"
            )
        return _with_curated_rows(key, _nist_lines(key))
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
