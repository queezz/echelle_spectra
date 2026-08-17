"""UI-independent campaign memory for the live calibration bench.

The live Qt window is deliberately a thin adapter over this module.  Exposure
triage, explicit file-role classification, the data-driven procedure, exposure
guidance, shared-catalog line help, sphere-factor comparison, TOML composition,
and snapshot save/validation can therefore be rehearsed without an event loop.

The bench accepts whatever is dropped on it: lamp names are free text, the
procedure is derived from the roles actually assigned, and no lamp family is
ever a gate.
"""

from __future__ import annotations

import bisect
import json
import os
import re
import shutil
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from enum import Enum
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import ndimage

try:  # pragma: no cover - selected by the running Python version
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

from .calibration_bench import AlignmentState, BenchFrame, CalibrationBenchSession
from .snapshot import (
    Snapshot,
    SnapshotError,
    create_snapshot,
    load_snapshot,
    reference_path,
)
from .tools.calibration_alignment import (
    AlignmentSettings,
    CalibrationTableLine,
    RigidTransform,
    apply_rigid_correction_to_lines,
    load_wavelength_table,
    save_alignment_settings,
    table_vetting,
    write_wavelength_table,
)
from .tools.line_catalog import LINE_FAMILY_LABELS, SpectralLine, load_line_table

__all__ = [
    "ALIGNMENT_SETTINGS_FILENAME",
    "KNOWN_LAMP_NAMES",
    "LAMP_TABLE_SPECIES",
    "PREVIOUS_CAMPAIGN_LAMPS",
    "AbsoluteCalibrationResult",
    "CalibrationCampaignSession",
    "CatalogOrderLine",
    "ChecklistItem",
    "ChecklistState",
    "ComparisonState",
    "CountsHistogram",
    "ExposureGuidance",
    "ExposureState",
    "ExposureTriage",
    "ExpectedLine",
    "FileRoleSuggestion",
    "LampReferenceSet",
    "LoadedFrameRecord",
    "MeasurementRecord",
    "MeasurementRole",
    "ReferenceState",
    "RoleTriageVerdict",
    "SaturationClusters",
    "SaveState",
    "TomlState",
    "WavelengthCorrection",
    "catalog_family_for_lamp",
    "catalog_lines_for_order",
    "catalog_mismatch_warning",
    "compute_absolute_calibration",
    "counts_histogram",
    "default_validity",
    "evaluate_exposure",
    "expected_lines_for_order",
    "lamp_reference_set",
    "measure_saturation_clusters",
    "normalize_lamp_name",
    "suggest_file_roles",
    "triage_exposure",
    "triage_for_role",
    "write_corrected_wavelength_table",
]

ALIGNMENT_SETTINGS_FILENAME = "alignment.toml"

#: Detector full scale for the 16-bit Andor cameras this bench serves.
FULL_SCALE_COUNTS = 65535.0

#: Lamp names the bench offers as ready-made choices.  Any other name is
#: accepted as free text; the list is convenience, never a permitted set.
KNOWN_LAMP_NAMES = ("ThAr", "Ne", "Hg", "H2")

#: What the previous campaign actually measured, offered as a suggestion.
PREVIOUS_CAMPAIGN_LAMPS = ("Ne",)

# Below this the correction would only reformat the table, never move a row.
IDENTITY_SHIFT_PX = 1e-6

#: The export config a bench-composed one inherits from until the next campaign
#: re-measures these values.  The bench never measures experiment timing or the
#: absolute-calibration crop, but a cube exported without them cannot produce
#: the frozen LHD text header, so the previous campaign's values are carried
#: forward explicitly and labelled as inherited rather than left out.
INHERITED_EXPORT_CONFIG_ID = "lhd_cmos_20250926"
INHERITED_TRIGGER_DELAY_S = 2.50
INHERITED_TIME_AXIS_REFERENCE = "LHD discharge time"
INHERITED_FRAME_TIME_FORMULA = "trigger_delay_s + frame * frame_interval_s"
INHERITED_CROP_MEASURED_AT = "2026-06-05"
INHERITED_WAVELENGTH_MIN_NM = 403.0


class MeasurementRole(Enum):
    """Explicit role assigned to one observed SIF."""

    SPHERE = "sphere"
    SPHERE_BACKGROUND = "sphere-background"
    LAMP = "lamp"
    LAMP_BACKGROUND = "lamp-background"
    OTHER = "other"


#: The two lamp roles.  ``MeasurementRole.OTHER`` parks an experiment frame on
#: the bench without giving it any calibration meaning.
LAMP_ROLES = (MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND)

#: The two background roles.  A background is shot with the light off; being
#: dark is the whole of its job, so DIM is its correct state and telling its
#: operator to expose it longer is advice that would spoil it (owner, deadpan:
#: "Background is dim. Wow, genius").
BACKGROUND_ROLES = (MeasurementRole.LAMP_BACKGROUND, MeasurementRole.SPHERE_BACKGROUND)

#: Which signal role a background belongs to, so its counts can be read against
#: the frame it will be subtracted from.
BACKGROUND_PARTNERS = {
    MeasurementRole.LAMP_BACKGROUND: MeasurementRole.LAMP,
    MeasurementRole.SPHERE_BACKGROUND: MeasurementRole.SPHERE,
}

#: The share of its partner's peak at which a background stops being a
#: background.  A dark frame reading a third of the signal it is subtracted
#: from is a light leak, a shutter that stayed open, or the wrong file — and
#: that, not dimness, is what deserves to be loud.
BACKGROUND_LEAK_FRACTION = 0.30

#: The roles exactly one file of a campaign can hold.  Two files proposing the
#: same one is a filename clash the bench must not resolve on its own — which
#: is the whole test of whether a folder's names may be applied unasked.  The
#: two lamp roles are deliberately absent: a lamp family is shot as a
#: bright/dim pair, so several files legitimately propose the same lamp role,
#: and the campaign already reads the last one assigned as the pair's member.
SINGLE_HOLDER_ROLES = (MeasurementRole.SPHERE, MeasurementRole.SPHERE_BACKGROUND)


class ExposureState(Enum):
    """Whole-frame exposure verdict used to name the next safe action."""

    NO_DATA = "no-data"
    SATURATED = "saturated"
    DIM = "dim"
    GOOD = "good"


class ComparisonState(Enum):
    """Lifecycle of candidate-versus-previous absolute factors."""

    NOT_RUN = "not-run"
    READY = "ready"
    INSUFFICIENT_DATA = "insufficient-data"
    FAILED = "failed"


class ChecklistState(Enum):
    """Derived procedure-item state."""

    WAITING = "waiting"
    DONE = "done"
    ATTENTION = "attention"
    SUGGESTION = "suggestion"


class TomlState(Enum):
    """Generated configuration lifecycle."""

    NOT_GENERATED = "not-generated"
    GENERATED = "generated"
    FAILED = "failed"


class ReferenceState(Enum):
    """Lifecycle of the reference rows one assigned lamp earns."""

    UNSCOPED = "unscoped"
    MATCHED = "matched"
    NO_ROWS = "no-rows"
    NO_CATALOG = "no-catalog"


class SaveState(Enum):
    """Snapshot save and validator lifecycle."""

    NOT_READY = "not-ready"
    READY = "ready"
    SAVING = "saving"
    VALIDATED = "validated"
    FAILED = "failed"


@dataclass(frozen=True)
class FileRoleSuggestion:
    """Filename-based help that never confirms a measurement role by itself."""

    roles: tuple[MeasurementRole, ...]
    reason: str
    lamp_name: str = ""

    @property
    def is_unambiguous(self) -> bool:
        return len(self.roles) == 1


@dataclass(frozen=True)
class ExposureGuidance:
    """Raw-count summary and an explicit next acquisition action."""

    state: ExposureState
    peak_value: float | None
    saturated_pixels: int
    saturated_fraction: float
    exposure_s: float | None
    suggested_exposure_s: float | None
    next_action: str
    anomalous_pixels: int = 0


@dataclass(frozen=True)
class CountsHistogram:
    """Raw-count distribution of one frame stack, ready to draw."""

    edges: np.ndarray
    counts: np.ndarray
    full_scale: float


@dataclass(frozen=True)
class SaturationClusters:
    """Full-scale pixels separated into real saturation and lone anomalies.

    A connected cluster of at least two full-scale pixels is real saturation.
    An isolated full-scale pixel is an anomaly — a cosmic ray or a hot pixel —
    and is counted rather than held against the frame.
    """

    saturation_level: float
    cluster_count: int
    cluster_pixels: int
    largest_cluster_pixels: int
    anomalous_pixels: int
    peak_counts: float | None
    clean_peak_counts: float | None
    finite_pixels: int

    @property
    def is_saturated(self) -> bool:
        return self.cluster_count > 0


@dataclass(frozen=True)
class ExposureTriage:
    """One-glance exposure verdict for a freshly loaded frame, before roles."""

    state: ExposureState
    headline: str
    details: tuple[str, ...]
    full_scale: float
    saturation: SaturationClusters
    histogram: CountsHistogram
    top_histogram: CountsHistogram
    headroom_fraction: float | None
    safe_gain: float | None
    noise_floor: float | None
    peak_over_noise: float | None
    #: The peak the verdict was judged on: the brightest non-anomalous pixel.
    peak_counts: float | None = None

    @property
    def is_usable(self) -> bool:
        """Whether lines may be fitted on this frame at all."""

        return self.state in {ExposureState.GOOD, ExposureState.DIM}


@dataclass(frozen=True)
class RoleTriageVerdict:
    """How one frame's exposure reads once its measurement role is known.

    The front-door triage judges a file before any role exists and never
    changes.  This is the second reading, taken once the operator has said what
    the frame is for: a lamp frame is *supposed* to saturate its strong lines —
    that is precisely why bright/dim pairs are shot — so clustered saturation on
    a lamp frame is information, not a failure.  A sphere frame keeps the hard
    verdict, because a saturated sphere frame cannot produce absolute factors.
    """

    role: MeasurementRole | None
    state: ExposureState
    #: Short badge for the files table.
    label: str
    #: The single line the operator reads beside the file.
    headline: str
    #: Whether this verdict refuses the frame the purpose it was given.
    blocking: bool
    #: The longer explanation, for a tooltip or the details dock.
    advice: str
    #: What to do next, when the role changes the answer.  Empty means the
    #: role-blind exposure guidance already says it correctly; a background
    #: frame is the case where it does not, because that guidance would have
    #: the operator expose a dark frame for twenty seconds.
    next_action: str = ""

    @property
    def is_usable(self) -> bool:
        """Whether this frame may still serve the role it was given."""

        return not self.blocking

    @property
    def is_background(self) -> bool:
        """Whether this verdict is reading a frame shot with the light off.

        The surface asks this, not the label: a label is a badge to reword and
        the routing behind colour, marks and advice must not move when it is.
        A frame with no data at all never reached the background reading, so it
        is not one of these however it was classified.
        """

        return self.role in BACKGROUND_ROLES and self.state is not ExposureState.NO_DATA


@dataclass(frozen=True)
class MeasurementRecord:
    """One explicitly classified, read-only source measurement."""

    path: Path
    role: MeasurementRole
    lamp_family: str = ""
    exposure: ExposureGuidance | None = None
    size_bytes: int = 0
    modified_ns: int = 0


@dataclass(frozen=True)
class LoadedFrameRecord:
    """A file the bench has read and triaged, whatever it turns out to be."""

    path: Path
    triage: ExposureTriage
    exposure: ExposureGuidance
    suggestion: FileRoleSuggestion


@dataclass(frozen=True)
class AbsoluteCalibrationResult:
    """Absolute factor curve produced by the established calibration engine."""

    wavelength_nm: np.ndarray
    factors_wmsr: np.ndarray


@dataclass(frozen=True)
class SphereComparison:
    """Candidate/previous factor comparison over their finite shared samples."""

    state: ComparisonState
    reason: str
    candidate: AbsoluteCalibrationResult | None = None
    previous: AbsoluteCalibrationResult | None = None
    sample_count: int = 0
    median_ratio: float | None = None
    p05_ratio: float | None = None
    p95_ratio: float | None = None


@dataclass(frozen=True)
class ChecklistItem:
    """One derived procedure row.

    A row that is not yet possible names what would unblock it, and a row that
    is only advice is marked non-blocking so it can never jail the operator.
    """

    key: str
    label: str
    state: ChecklistState
    detail: str
    blocking: bool = True
    unblocked_by: str = ""


@dataclass(frozen=True)
class WavelengthCorrection:
    """What the saved ``wavelength.txt`` is, and why it is that."""

    applied: bool
    reason: str
    max_shift_px: float
    line_count: int


@dataclass(frozen=True)
class CatalogOrderLine:
    """A shared-catalog line mapped approximately onto one detector order."""

    line: SpectralLine
    order_idx: int
    detector_pixel: float


@dataclass(frozen=True)
class ExpectedLine:
    """One line the assigned lamp is expected to show inside one order.

    This is the bench's *only* expected-line source.  The labelled stick on the
    spectrum, the row in the expected-lines panel, and the per-order count in
    that panel's header are all rendered from a list of these and from nothing
    else, so the three can never disagree again.

    The row is the anchorable curated wavelength-table row — the same object a
    click snaps to — because a line the operator can see but cannot anchor is
    not help.  The packaged NIST entry, when the caches happen to carry that
    wavelength, rides along as annotation (relative intensity, provenance) and
    is never allowed to add or remove a line.
    """

    row: CalibrationTableLine
    order_idx: int
    detector_pixel: float
    label: str
    wavelength_nm: float
    species: str
    relative_intensity: float | None
    source: str
    catalog: SpectralLine | None = None


@dataclass(frozen=True)
class LampReferenceSet:
    """The curated lookup rows one assigned lamp's own catalog contributes.

    Click-to-fit measures every candidate and accepted anchor against these
    rows and no others, so a neon frame is never referenced against a thorium
    line.  A lamp the packaged catalogs have never met contributes no rows at
    all and says why, which is a better answer than another element's lines.
    """

    lamp: str
    catalog_family: str
    catalog_label: str
    species: tuple[str, ...]
    lines: tuple[CalibrationTableLine, ...]
    state: ReferenceState
    message: str

    @property
    def is_referenceable(self) -> bool:
        """Whether an anchor may be looked up against this set at all."""

        return self.state in {ReferenceState.UNSCOPED, ReferenceState.MATCHED}

    @property
    def best_order(self) -> int | None:
        """The order carrying most of this lamp's own rows, if it has any."""

        counts: dict[int, int] = {}
        for line in self.lines:
            counts[line.order_idx] = counts.get(line.order_idx, 0) + 1
        if not counts:
            return None
        return max(counts, key=lambda order: (counts[order], -order))


_LAMP_ALIASES = {
    "thar": "ThAr",
    "th": "ThAr",
    "ne": "Ne",
    "neon": "Ne",
    "hg": "Hg",
    "mercury": "Hg",
    "h2": "H2",
    "fulcher": "H2",
    "deuterium": "D2",
    "d2": "D2",
}

_LAMP_NAME_EXTRA_CHARACTERS = "+-_."

_CATALOG_FAMILIES = {"ThAr": "thar", "Ne": "ne", "Hg": "hg", "H2": "fulcher"}

#: Which curated wavelength-table species each lamp may legitimately emit.
#:
#: ThAr, Ne, and Hg carry their NIST lamp preset's own spectra, spelled as the
#: curated table spells them.  ``ThBlend`` is deliberately absent: it is a
#: weighted compound centroid of two transitions rather than one line, so it is
#: never an anchor.  H2 carries the molecular Fulcher rows together with the
#: atomic Balmer rows a hydrogen discharge shows beside them.
LAMP_TABLE_SPECIES: dict[str, tuple[str, ...]] = {
    "ThAr": ("ThI", "ThII", "ArI", "ArII"),
    "Ne": ("NeI", "NeII"),
    "Hg": ("HgI", "HgII"),
    "H2": ("H2", "H-a", "H-g"),
}


def normalize_lamp_name(value: str) -> str:
    """Return the canonical spelling of a lamp name, accepting any lamp.

    Known spellings collapse onto their canonical form so ``th-ar`` and
    ``ThAr`` are one lamp.  Every other non-empty, path-safe name is kept as
    the operator typed it: the bench never refuses a lamp it has not met.
    """

    text = str(value).strip()
    if not text:
        raise ValueError("a lamp role needs a lamp name")
    compact = text.casefold().replace("-", "").replace(" ", "")
    if compact in _LAMP_ALIASES:
        return _LAMP_ALIASES[compact]
    if not all(
        character.isalnum() or character in _LAMP_NAME_EXTRA_CHARACTERS
        for character in text
    ):
        raise ValueError(
            "a lamp name may hold only letters, digits, '+', '-', '_', and '.'"
        )
    return text


def catalog_family_for_lamp(lamp_name: str) -> str | None:
    """Return the packaged line-catalog family for *lamp_name*, if one exists."""

    try:
        return _CATALOG_FAMILIES.get(normalize_lamp_name(lamp_name))
    except ValueError:
        return None


def lamp_reference_set(
    lamp_name: str, lines: Sequence[CalibrationTableLine]
) -> LampReferenceSet:
    """Return the curated rows *lamp_name* earns from a wavelength table.

    A lamp the catalogs know keeps only its own species' rows, which is what
    stops a nearest-row lookup from ever crossing lamps.  A free-text lamp with
    no catalog keeps none and states that plainly, and a lamp whose catalog
    exists but whose species this particular table never carries says that too.
    """

    rows = tuple(lines)
    text = str(lamp_name).strip()
    if not text:
        return LampReferenceSet(
            "",
            "",
            "",
            (),
            rows,
            ReferenceState.UNSCOPED,
            f"no lamp is assigned, so all {len(rows)} curated rows stay clickable; "
            "assign a lamp role so anchors reference that lamp's own lines",
        )
    try:
        lamp = normalize_lamp_name(text)
    except ValueError:
        lamp = text
    family = catalog_family_for_lamp(lamp)
    if family is None:
        return LampReferenceSet(
            lamp,
            "",
            "",
            (),
            (),
            ReferenceState.NO_CATALOG,
            f"no line catalog for {lamp} — anchors cannot be auto-referenced",
        )
    label = LINE_FAMILY_LABELS.get(family, lamp)
    species = LAMP_TABLE_SPECIES[lamp]
    scoped = tuple(row for row in rows if row.species in species)
    if not scoped:
        return LampReferenceSet(
            lamp,
            family,
            label,
            species,
            (),
            ReferenceState.NO_ROWS,
            f"{lamp} has a packaged {label} catalog, but none of the {len(rows)} "
            f"curated rows carry its species ({', '.join(species)}) — anchors "
            "cannot be auto-referenced from this table",
        )
    return LampReferenceSet(
        lamp,
        family,
        label,
        species,
        scoped,
        ReferenceState.MATCHED,
        f"{len(scoped)} of {len(rows)} curated rows are {label} lines "
        f"({', '.join(species)}); anchors reference {label} only",
    )


def catalog_mismatch_warning(active_catalog: str, assigned_lamp: str) -> str:
    """Return a warning when the shown catalog is not the assigned lamp's own.

    The fit itself always uses the assigned lamp's rows.  This only tells the
    operator that the identification help beside those rows belongs to a
    different lamp, so the blue sticks and the clickable rows disagree on
    purpose rather than by accident.
    """

    lamp = str(assigned_lamp).strip()
    active = str(active_catalog).strip()
    if not lamp or not active:
        return ""
    try:
        lamp = normalize_lamp_name(lamp)
    except ValueError:  # a free-text lamp is compared as the operator typed it
        pass
    try:
        active = normalize_lamp_name(active)
    except ValueError:
        pass
    if active.casefold() == lamp.casefold():
        return ""
    return (
        f"the {active} line help on screen is not {lamp}'s own catalog; "
        f"anchors are still referenced against {lamp} lines only"
    )


def _lamp_name_in_filename(stem: str) -> str:
    """Return the lamp name a filename hints at, or an empty string."""

    tokens = [token for token in re.split(r"[^0-9a-zA-Z]+", stem) if token]
    for token in tokens:
        compact = token.casefold()
        if compact in _LAMP_ALIASES:
            return _LAMP_ALIASES[compact]
    return ""


def suggest_file_roles(path: str | Path) -> FileRoleSuggestion:
    """Suggest a likely role from a filename without ever accepting it.

    The suggestion pre-fills the manual controls and nothing else.  A file
    whose name says nothing still loads, still gets triaged, and can still be
    given any role by hand.
    """

    stem = Path(path).stem
    folded = stem.casefold()
    background = any(token in folded for token in ("background", "_bg", "-bg", "bkg", "dark"))
    sphere = any(token in folded for token in ("sphere", "sphr", "absolute"))
    lamp_name = _lamp_name_in_filename(stem)
    lamp = bool(lamp_name) or "lamp" in folded
    if sphere and background:
        return FileRoleSuggestion(
            (MeasurementRole.SPHERE_BACKGROUND,),
            "filename looks like an integrating-sphere background; confirm explicitly",
        )
    if sphere:
        return FileRoleSuggestion(
            (MeasurementRole.SPHERE,),
            "filename looks like an integrating-sphere signal; confirm explicitly",
        )
    if lamp and background:
        return FileRoleSuggestion(
            (MeasurementRole.LAMP_BACKGROUND,),
            f"filename looks like a {lamp_name or 'lamp'} background; "
            "confirm the role and lamp name",
            lamp_name,
        )
    if lamp:
        return FileRoleSuggestion(
            (MeasurementRole.LAMP,),
            f"filename looks like a {lamp_name or 'lamp'} signal; "
            "confirm the role and lamp name",
            lamp_name,
        )
    if background:
        return FileRoleSuggestion(
            (MeasurementRole.SPHERE_BACKGROUND, MeasurementRole.LAMP_BACKGROUND),
            "background is ambiguous until its source is selected explicitly",
        )
    return FileRoleSuggestion(
        tuple(MeasurementRole),
        "filename does not support a unique role; assign one by hand",
    )


def _metadata_exposure_s(metadata: Mapping[str, object]) -> float | None:
    for key in ("ExposureTime", "exposure_time", "exposure_s"):
        value = metadata.get(key)
        try:
            exposure = float(value)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
        if np.isfinite(exposure) and exposure > 0:
            return exposure
    return None


#: Four-connectivity inside one frame.  Saturation spreads across neighbouring
#: pixels of the same exposure; it never links two frames, so a hot pixel that
#: repeats every frame stays a repeated anomaly instead of becoming a cluster.
_IN_FRAME_CROSS = np.array(
    [[False, True, False], [True, True, True], [False, True, False]]
)

#: Cap on pixels sampled for the robust noise floor; the frames are large and
#: the floor is a property of the background, not of the sample size.
_NOISE_SAMPLE_LIMIT = 250_000


def measure_saturation_clusters(
    images: np.ndarray,
    *,
    saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
) -> SaturationClusters:
    """Separate clustered saturation from isolated full-scale anomalies.

    ``images`` may be ``(rows, cols)`` or ``(frames, rows, cols)``.  Clusters
    are found with four-connectivity inside each frame: two or more touching
    full-scale pixels are real saturation, while a lone full-scale pixel is a
    cosmic ray or hot pixel and is reported as an anomaly instead.
    """

    stack = np.asarray(images, dtype=float)
    if stack.ndim == 2:
        stack = stack[np.newaxis, :, :]
    if stack.ndim != 3:
        raise ValueError("images must have shape (rows, cols) or (frames, rows, cols)")
    level = float(saturation_level)

    finite_pixels = 0
    peak: float | None = None
    clean_peak: float | None = None
    cluster_count = 0
    cluster_pixels = 0
    largest_cluster = 0
    anomalies = 0
    for frame in stack:
        finite = np.isfinite(frame)
        finite_pixels += int(np.count_nonzero(finite))
        if not finite.any():
            continue
        peak = _maximum(peak, float(np.max(frame[finite])))
        saturated = finite & (frame >= level)
        if not saturated.any():
            clean_peak = _maximum(clean_peak, float(np.max(frame[finite])))
            continue
        labels, _found = ndimage.label(saturated, structure=_IN_FRAME_CROSS)
        sizes = np.bincount(labels.reshape(-1))
        sizes[0] = 0
        anomalies += int(np.count_nonzero(sizes == 1))
        clustered = sizes[sizes >= 2]
        cluster_count += int(clustered.size)
        cluster_pixels += int(clustered.sum())
        largest_cluster = max(largest_cluster, int(clustered.max()) if clustered.size else 0)
        lonely = (sizes == 1)[labels]
        candidates = finite & ~lonely
        if candidates.any():
            clean_peak = _maximum(clean_peak, float(np.max(frame[candidates])))
    return SaturationClusters(
        level,
        cluster_count,
        cluster_pixels,
        largest_cluster,
        anomalies,
        peak,
        clean_peak,
        finite_pixels,
    )


def _maximum(current: float | None, candidate: float) -> float:
    return candidate if current is None else max(current, candidate)


def counts_histogram(
    values: np.ndarray,
    *,
    full_scale: float = FULL_SCALE_COUNTS,
    bins: int = 96,
    lowest: float | None = None,
) -> CountsHistogram:
    """Bin finite raw counts between *lowest* and detector full scale."""

    finite = np.asarray(values, dtype=float).reshape(-1)
    finite = finite[np.isfinite(finite)]
    smallest = float(finite.min()) if finite.size else 0.0
    largest = float(finite.max()) if finite.size else float(full_scale)
    low = float(lowest) if lowest is not None else min(0.0, smallest)
    high = max(float(full_scale), largest)
    if high <= low:
        high = low + 1.0
    counts, edges = np.histogram(finite, bins=int(bins), range=(low, high))
    return CountsHistogram(edges, counts, float(full_scale))


def _noise_floor(finite: np.ndarray) -> tuple[float, float]:
    """Return a robust (baseline, noise) pair from a strided pixel sample."""

    stride = max(1, finite.size // _NOISE_SAMPLE_LIMIT)
    sample = finite[::stride]
    baseline = float(np.median(sample))
    noise = float(1.4826 * np.median(np.abs(sample - baseline)))
    if noise <= 0:
        noise = float(np.std(sample))
    if noise <= 0:
        noise = 1.0
    return baseline, noise


def triage_exposure(
    frame: BenchFrame,
    *,
    saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
    full_scale: float = FULL_SCALE_COUNTS,
    dim_fraction: float = 0.20,
    minimum_peak_snr: float = 10.0,
    safe_fraction: float = 0.90,
    top_fraction: float = 0.90,
) -> ExposureTriage:
    """Judge one freshly read frame before it is given any role.

    This is the bench's front door: it needs nothing but a file.  The verdict
    reports clustered saturation, counted anomalies, remaining headroom, and
    whether the frame is bright enough above its own noise floor for lines to
    fit well.
    """

    clusters = measure_saturation_clusters(frame.images, saturation_level=saturation_level)
    finite = np.asarray(frame.images, dtype=float).reshape(-1)
    finite = finite[np.isfinite(finite)]
    histogram = counts_histogram(finite, full_scale=full_scale)
    top_histogram = counts_histogram(
        finite, full_scale=full_scale, bins=48, lowest=top_fraction * full_scale
    )
    # A frame whose only bright pixels are lone spikes still has a real peak:
    # fall back to it rather than pretending the frame holds nothing.
    reference_peak = (
        clusters.clean_peak_counts
        if clusters.clean_peak_counts is not None
        else clusters.peak_counts
    )
    if not finite.size or reference_peak is None:
        return ExposureTriage(
            ExposureState.NO_DATA,
            "NO DATA — the frame holds no finite raw pixels; reacquire.",
            ("Nothing can be judged from this file.",),
            float(full_scale),
            clusters,
            histogram,
            top_histogram,
            None,
            None,
            None,
            None,
        )

    baseline, noise = _noise_floor(finite)
    clean_peak = float(reference_peak)
    headroom = clean_peak / float(full_scale)
    safe_gain = (safe_fraction * float(full_scale)) / max(clean_peak, np.finfo(float).eps)
    peak_over_noise = (clean_peak - baseline) / noise
    state = _triage_state(clusters, headroom, peak_over_noise, dim_fraction, minimum_peak_snr)
    details = _triage_details(clusters, clean_peak, headroom, safe_gain, baseline, peak_over_noise)
    return ExposureTriage(
        state,
        _triage_headline(state, clusters, headroom, safe_gain, peak_over_noise),
        details,
        float(full_scale),
        clusters,
        histogram,
        top_histogram,
        float(headroom),
        float(safe_gain),
        float(baseline),
        float(peak_over_noise),
        clean_peak,
    )


def _triage_state(
    clusters: SaturationClusters,
    headroom: float,
    peak_over_noise: float,
    dim_fraction: float,
    minimum_peak_snr: float,
) -> ExposureState:
    if clusters.is_saturated:
        return ExposureState.SATURATED
    if headroom < dim_fraction or peak_over_noise < minimum_peak_snr:
        return ExposureState.DIM
    return ExposureState.GOOD


def _triage_headline(
    state: ExposureState,
    clusters: SaturationClusters,
    headroom: float,
    safe_gain: float,
    peak_over_noise: float,
) -> str:
    """Compose the single line the operator reads between two exposures."""

    tail = ""
    if clusters.anomalous_pixels:
        tail = (
            f" {clusters.anomalous_pixels} isolated full-scale pixel(s) are anomalies, "
            "not saturation."
        )
    if state is ExposureState.SATURATED:
        return (
            f"SATURATED — {clusters.cluster_pixels} full-scale pixel(s) in "
            f"{clusters.cluster_count} connected cluster(s); lower the exposure "
            f"and shoot again.{tail}"
        )
    if state is ExposureState.DIM:
        return (
            f"TOO DIM FOR LINES — brightest real pixel at {100.0 * headroom:.1f}% of "
            f"full scale, only {peak_over_noise:.0f}x the noise floor: right for a "
            f"background, too weak for a lamp. About {safe_gain:.1f}x brighter is safe."
            f"{tail}"
        )
    return (
        f"HEALTHY — brightest real pixel at {100.0 * headroom:.1f}% of full scale; "
        f"about {safe_gain:.1f}x brighter is still safe.{tail}"
    )


def _triage_details(
    clusters: SaturationClusters,
    clean_peak: float,
    headroom: float,
    safe_gain: float,
    baseline: float,
    peak_over_noise: float,
) -> tuple[str, ...]:
    anomaly_note = (
        f"{clusters.anomalous_pixels} isolated full-scale pixel(s) counted as "
        "anomalies (cosmic rays or hot pixels), not saturation"
        if clusters.anomalous_pixels
        else "no isolated full-scale pixels"
    )
    saturation_note = (
        f"{clusters.cluster_pixels} full-scale pixel(s) in {clusters.cluster_count} "
        f"cluster(s), largest {clusters.largest_cluster_pixels}"
        if clusters.is_saturated
        else "no connected full-scale cluster"
    )
    return (
        f"Brightest non-anomalous pixel {clean_peak:.0f} counts "
        f"({100.0 * headroom:.1f}% of full scale) — about {safe_gain:.1f}x brighter is safe.",
        f"Saturation: {saturation_note}.",
        f"Anomalies: {anomaly_note}.",
        f"Noise floor {baseline:.0f} counts; peak stands {peak_over_noise:.0f}x above it.",
    )


def evaluate_exposure(
    frame: BenchFrame,
    *,
    saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
    dim_fraction: float = 0.20,
    target_fraction: float = 0.70,
) -> ExposureGuidance:
    """State the next safe acquisition action for one triaged frame."""

    triage = triage_exposure(
        frame, saturation_level=saturation_level, dim_fraction=dim_fraction
    )
    exposure_s = _metadata_exposure_s(frame.metadata)
    clusters = triage.saturation
    peak = triage.peak_counts
    if triage.state is ExposureState.NO_DATA or peak is None:
        return ExposureGuidance(
            ExposureState.NO_DATA,
            None,
            0,
            0.0,
            exposure_s,
            None,
            "No finite raw pixels are available; reacquire before continuing.",
            clusters.anomalous_pixels,
        )
    scale = target_fraction * saturation_level / max(peak, np.finfo(float).eps)
    suggested = exposure_s * scale if exposure_s is not None else None
    saturated_fraction = (
        clusters.cluster_pixels / clusters.finite_pixels if clusters.finite_pixels else 0.0
    )
    if triage.state is ExposureState.SATURATED:
        action = "Lower exposure and reacquire; do not accept anchors from this frame."
        if suggested is not None:
            action = f"Lower exposure to about {suggested:.4g} s and reacquire."
    elif triage.state is ExposureState.DIM:
        action = "Increase exposure, then reacquire for stronger unsaturated lines."
        if suggested is not None:
            action = f"Increase exposure toward {suggested:.4g} s, then reacquire."
    else:
        action = "Exposure is usable. Continue with line identification and anchor fitting."
        suggested = exposure_s
    return ExposureGuidance(
        triage.state,
        peak,
        clusters.cluster_pixels,
        saturated_fraction,
        exposure_s,
        suggested,
        action,
        clusters.anomalous_pixels,
    )


#: What a lamp frame is told when its strong lines are clustered-saturated.
LAMP_SATURATION_ADVICE = (
    "A dim-series lamp exposure is meant to saturate its strong lines so the "
    "weak ones emerge; that is what the bright/dim pair is for. Nothing here "
    "is refused: every anchor is saturation-checked on the raw detector "
    "window of its own line, so a click on a saturated line is rejected while "
    "the unsaturated lines beside it fit normally."
)

#: Why a saturated sphere frame is still a hard failure.
SPHERE_SATURATION_ADVICE = (
    "Absolute factors divide the sphere signal by a known radiance, so a "
    "clipped sphere pixel is an unknown number rather than a bright one. "
    "Lower the exposure and reacquire both sphere frames together."
)


def _background_verdict(
    triage: ExposureTriage,
    role: MeasurementRole,
    partner_peak: float | None,
) -> RoleTriageVerdict:
    """Read a background frame by the only standard that applies to it."""

    kind = "lamp" if role is MeasurementRole.LAMP_BACKGROUND else "sphere"
    peak = triage.peak_counts
    partner = BACKGROUND_PARTNERS[role].value
    leaking = (
        peak is not None
        and partner_peak is not None
        and partner_peak > 0
        and peak >= BACKGROUND_LEAK_FRACTION * partner_peak
    )
    if leaking:
        share = 100.0 * peak / partner_peak
        return RoleTriageVerdict(
            role,
            triage.state,
            "background too bright",
            f"reads {share:.0f}% of its {kind} signal — not dark",
            True,
            f"This {kind} background peaks at {peak:.0f} counts against the "
            f"{partner_peak:.0f} counts of the {partner} it will be subtracted "
            "from. A background is shot with the light off, so at that level "
            "something is reaching the detector: a light leak, a shutter that "
            "did not close, or the wrong file in this role. Subtracting it "
            "would remove real signal. Check which file this is and reshoot "
            "the background with the lamp off.",
            f"Reshoot this background with the {kind} off, or correct which "
            "file carries the background role.",
        )
    if triage.state is ExposureState.SATURATED:
        return RoleTriageVerdict(
            role,
            triage.state,
            "background saturated",
            "saturated — it cannot be subtracted",
            True,
            f"A background at full scale is not a background. {triage.headline} "
            "Reshoot it with the light off before using this pair.",
            "Reshoot this background with the light off.",
        )
    return RoleTriageVerdict(
        role,
        triage.state,
        "background",
        "",
        False,
        f"This is a {kind} background, so dark is correct and DIM is the "
        "reading to expect — it is subtracted from its "
        f"{partner}, never measured for lines. Exposure is not raised on a "
        "background; it is matched to the frame it partners. The reading worth "
        "watching is the opposite one: a background approaching its signal's "
        "own counts means light is getting in.",
        "",
    )


def triage_for_role(
    triage: ExposureTriage,
    role: MeasurementRole | None,
    partner_peak: float | None = None,
) -> RoleTriageVerdict:
    """Read one already-triaged frame again, in the light of its role.

    Two roles change what a state means.  A lamp frame is never failed for
    saturating the strong lines it was shot to saturate — that is precisely why
    bright/dim pairs exist.  And a background frame is never scolded for being
    dark: darkness is its entire purpose, so DIM reads as correct and the
    advice never tells anybody to expose it longer.  What a background *can*
    get wrong is being bright, which is why ``partner_peak`` is read against
    it: a dark frame carrying a third of the signal it will be subtracted from
    is a leak or the wrong file, and that gets loud.
    """

    clusters = triage.saturation
    if role in BACKGROUND_ROLES and triage.state is not ExposureState.NO_DATA:
        return _background_verdict(triage, role, partner_peak)
    if triage.state is ExposureState.SATURATED and role in LAMP_ROLES:
        headline = (
            f"saturated in {clusters.cluster_count} cluster(s) — fit "
            "unsaturated lines only"
        )
        return RoleTriageVerdict(
            role,
            triage.state,
            "saturated lines (expected)",
            headline,
            False,
            f"{clusters.cluster_pixels} full-scale pixel(s) in "
            f"{clusters.cluster_count} connected cluster(s), largest "
            f"{clusters.largest_cluster_pixels}. {LAMP_SATURATION_ADVICE}",
        )
    if triage.state is ExposureState.SATURATED:
        return RoleTriageVerdict(
            role,
            triage.state,
            "saturated",
            triage.headline,
            True,
            f"{triage.headline} {SPHERE_SATURATION_ADVICE}"
            if role in {MeasurementRole.SPHERE, MeasurementRole.SPHERE_BACKGROUND}
            else triage.headline,
        )
    return RoleTriageVerdict(
        role,
        triage.state,
        triage.state.value,
        triage.headline,
        triage.state is ExposureState.NO_DATA,
        "\n".join(triage.details),
    )


def catalog_lines_for_order(
    calibration_rows: Sequence[CalibrationTableLine],
    order_idx: int,
    lamp_family: str,
    *,
    maximum_lines: int = 24,
) -> tuple[CatalogOrderLine, ...]:
    """Map shared packaged line knowledge onto an order by table interpolation."""

    rows = sorted(
        (row for row in calibration_rows if row.order_idx == order_idx),
        key=lambda row: row.wavelength_nm,
    )
    if len(rows) < 2:
        return ()
    wavelengths = np.asarray([row.wavelength_nm for row in rows], dtype=float)
    pixels = np.asarray([row.center_pixel for row in rows], dtype=float)
    minimum_nm = float(wavelengths[0])
    maximum_nm = float(wavelengths[-1])
    family = catalog_family_for_lamp(lamp_family)
    if family is None:
        return ()
    catalog = [
        line
        for line in load_line_table(family)
        if minimum_nm <= line.wavelength_nm <= maximum_nm
    ]
    catalog.sort(
        key=lambda line: (
            -(line.relative_intensity if line.relative_intensity is not None else 1.0),
            line.wavelength_nm,
        )
    )
    chosen = sorted(catalog[:maximum_lines], key=lambda line: line.wavelength_nm)
    return tuple(
        CatalogOrderLine(
            line,
            int(order_idx),
            float(np.interp(line.wavelength_nm, wavelengths, pixels)),
        )
        for line in chosen
    )


@lru_cache(maxsize=8)
def _catalog_by_wavelength(family: str) -> tuple[tuple[float, ...], tuple[SpectralLine, ...]]:
    """Return one packaged catalog family sorted by wavelength, read once."""

    if not family:
        return ((), ())
    ordered = sorted(load_line_table(family), key=lambda line: line.wavelength_nm)
    return (
        tuple(float(line.wavelength_nm) for line in ordered),
        tuple(ordered),
    )


def _nearest_catalog_line(
    family: str, wavelength_nm: float, tolerance_nm: float
) -> SpectralLine | None:
    """The packaged entry within *tolerance_nm* of a curated row, if any."""

    wavelengths, lines = _catalog_by_wavelength(family)
    if not wavelengths:
        return None
    index = bisect.bisect_left(wavelengths, wavelength_nm)
    best: SpectralLine | None = None
    best_distance = float(tolerance_nm)
    for candidate in (index - 1, index):
        if not 0 <= candidate < len(lines):
            continue
        distance = abs(wavelengths[candidate] - wavelength_nm)
        if distance <= best_distance:
            best_distance = distance
            best = lines[candidate]
    return best


def expected_lines_for_order(
    reference: LampReferenceSet | None,
    order_idx: int,
    *,
    fallback_lines: Sequence[CalibrationTableLine] = (),
    match_tolerance_nm: float = 0.05,
) -> tuple[ExpectedLine, ...]:
    """Return the single expected-line list one order shows for one lamp.

    The bench used to draw its labelled sticks from the curated wavelength
    table while filling the expected-lines panel from the packaged NIST caches
    interpolated onto the order.  Those are two different lists: the 2025 Ne
    cache stops at 638.3 nm, so order 7 drew a stick for the curated NeI
    640.225 row that the panel had never heard of, and orders 0-6 drew three
    sticks each under a panel reading "0 expected Ne lines in this order".

    There is one list now, and this returns it: the assigned lamp's own
    reference rows falling in this order, in detector order, annotated from the
    packaged catalog where it happens to know the wavelength.  A reference set
    that cannot anchor anything (a free-text lamp, or a lamp whose species this
    table never carries) contributes nothing, and its own message says why.
    """

    if reference is None:
        rows: Sequence[CalibrationTableLine] = tuple(fallback_lines)
        family = ""
    elif reference.is_referenceable:
        rows = reference.lines
        family = reference.catalog_family
    else:
        rows = ()
        family = ""
    order = int(order_idx)
    scoped = sorted(
        (row for row in rows if row.order_idx == order),
        key=lambda row: (row.center_pixel, row.wavelength_nm),
    )
    expected = []
    for row in scoped:
        match = _nearest_catalog_line(family, row.wavelength_nm, match_tolerance_nm)
        source = (
            row.comment.strip() or "the curated wavelength table"
            if match is None
            else match.source_resource.replace("\\", "/").rsplit("/", 1)[-1]
        )
        expected.append(
            ExpectedLine(
                row,
                order,
                float(row.center_pixel),
                f"{row.species} {row.wavelength_nm:.3f}",
                float(row.wavelength_nm),
                row.species,
                None if match is None else match.relative_intensity,
                source,
                match,
            )
        )
    return tuple(expected)


def compute_absolute_calibration(
    *,
    pattern: Path,
    wavelength: Path,
    sphere: Path,
    sphere_background: Path,
    integral: Path,
) -> AbsoluteCalibrationResult:
    """Run the established extraction/absolute-calibration engine for one pair."""

    from .tools.loader import build_calibration

    calibration = build_calibration(
        pattern.parent,
        "CMOS",
        calibration_files={
            "orders": str(pattern.resolve()),
            "wavelength": str(wavelength.resolve()),
            "sphr": str(sphere.resolve()),
            "bkgr": str(sphere_background.resolve()),
            "integral": str(integral.resolve()),
        },
    )
    return AbsoluteCalibrationResult(
        np.asarray(calibration.wavelength, dtype=float).reshape(-1),
        np.asarray(calibration.absolute["wmsr"], dtype=float).reshape(-1),
    )


def write_corrected_wavelength_table(
    base_wavelength: str | Path,
    destination: str | Path,
    *,
    pattern: np.ndarray,
    transform: RigidTransform,
    metadata: Sequence[tuple[str, str]] = (),
    identity_shift_px: float = IDENTITY_SHIFT_PX,
) -> WavelengthCorrection:
    """Write the base lookup table moved by *transform*, or copy it unchanged.

    A transform that moves no curated row would only reformat the base table,
    so its bytes are copied instead and the returned reason states which of the
    two outcomes the caller got.
    """

    source = Path(base_wavelength)
    target = Path(destination)
    rows = load_wavelength_table(source)
    if not rows:
        shutil.copy2(source, target)
        return WavelengthCorrection(
            False,
            f"{source.name} holds no correctable rows, so its bytes were copied unchanged",
            0.0,
            0,
        )
    adjusted = apply_rigid_correction_to_lines(rows, pattern, transform)
    max_shift_px = max(
        abs(new.center_pixel - old.center_pixel) for old, new in zip(rows, adjusted)
    )
    if max_shift_px <= float(identity_shift_px):
        shutil.copy2(source, target)
        return WavelengthCorrection(
            False,
            f"the solved transform moves no line of {source.name} measurably "
            f"({max_shift_px:.3g} px), so its bytes were copied unchanged",
            float(max_shift_px),
            len(rows),
        )
    write_wavelength_table(adjusted, target, metadata=metadata)
    return WavelengthCorrection(
        True,
        f"{len(rows)} rows of {source.name} were moved by the solved transform "
        f"(largest shift {max_shift_px:.3f} px)",
        float(max_shift_px),
        len(rows),
    )


def default_validity(
    snapshot_id: str, valid_from: date | str | None = None
) -> dict[str, str]:
    """Return an open-ended epoch that starts on the acquisition date."""

    if valid_from is None:
        digits = str(snapshot_id).strip()[:8]
        try:
            start = date(int(digits[:4]), int(digits[4:6]), int(digits[6:8]))
        except ValueError as exc:
            raise SnapshotError(
                f"no acquisition date can be read from snapshot id {snapshot_id!r}; "
                "state the epoch start explicitly"
            ) from exc
    elif isinstance(valid_from, date):
        start = valid_from
    else:
        try:
            start = date.fromisoformat(str(valid_from))
        except ValueError as exc:
            raise SnapshotError(
                f"epoch start must be an ISO YYYY-MM-DD date: {valid_from!r}"
            ) from exc
    return {"date_from": start.isoformat()}


def _toml_string(value: object) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _snapshot_reference(snapshot_root: str | Path | None, source: Path) -> str:
    """Name *source* from inside a snapshot folder, the way the binder does.

    With the snapshot folder known the path is computed exactly.  Without it,
    the bench's own layout is the honest assumption: the snapshot lands in
    ``<calibration folder>/calibrations/<id>`` and the frames stay in the
    calibration folder two levels up.
    """

    if snapshot_root is None:
        return f"../../{source.name}"
    return reference_path(snapshot_root, source)


class CalibrationCampaignSession:
    """Campaign-memory state transitions independent of Qt."""

    def __init__(
        self,
        *,
        pattern_source: str | Path,
        wavelength_source: str | Path,
        integral_source: str | Path,
        suggested_lamps: Sequence[str] = PREVIOUS_CAMPAIGN_LAMPS,
        previous_sphere: str | Path | None = None,
        previous_sphere_background: str | Path | None = None,
    ) -> None:
        self.pattern_source = Path(pattern_source)
        self.wavelength_source = Path(wavelength_source)
        self.integral_source = Path(integral_source)
        #: What the previous campaign measured.  Advice for this one, never a gate.
        self.suggested_lamps = tuple(
            dict.fromkeys(normalize_lamp_name(item) for item in suggested_lamps if str(item).strip())
        )
        self.previous_sphere = Path(previous_sphere) if previous_sphere else None
        self.previous_sphere_background = (
            Path(previous_sphere_background) if previous_sphere_background else None
        )
        self.observed: dict[Path, FileRoleSuggestion] = {}
        self.loaded: dict[Path, LoadedFrameRecord] = {}
        self.measurements: dict[Path, MeasurementRecord] = {}
        self.comparison = SphereComparison(ComparisonState.NOT_RUN, "not computed")
        self.toml_state = TomlState.NOT_GENERATED
        self.toml_paths: dict[str, Path] = {}
        self.toml_snapshot_id = ""
        self.save_state = SaveState.NOT_READY
        self.saved_snapshot: Snapshot | None = None
        self.wavelength_correction: WavelengthCorrection | None = None
        self.last_error = ""

    def observe_file(self, path: str | Path) -> FileRoleSuggestion:
        """Remember a file arrival as unclassified and return non-binding help."""

        source = Path(path)
        suggestion = suggest_file_roles(source)
        self.observed[source] = suggestion
        return suggestion

    def record_frame(
        self,
        frame: BenchFrame,
        *,
        saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
    ) -> LoadedFrameRecord:
        """Triage a freshly read frame before it has any role at all.

        The frame itself is not retained — only its verdict — so a bench can
        hold a whole campaign folder without holding every detector image.
        """

        suggestion = self.observe_file(frame.path)
        triage = triage_exposure(frame, saturation_level=saturation_level)
        record = LoadedFrameRecord(
            Path(frame.path),
            triage,
            evaluate_exposure(frame, saturation_level=saturation_level),
            suggestion,
        )
        self.loaded[record.path] = record
        return record

    def forget_file(self, path: str | Path) -> bool:
        """Drop one loaded file and any role it carried."""

        source = Path(path)
        self.observed.pop(source, None)
        removed = self.loaded.pop(source, None) is not None
        return self.remove_classification(source) or removed

    def role_triage(self, path: str | Path) -> RoleTriageVerdict | None:
        """Read one loaded file's exposure in the light of the role it carries."""

        source = Path(path)
        record = self.loaded.get(source)
        if record is None:
            return None
        measurement = self.measurements.get(source)
        role = measurement.role if measurement is not None else None
        return triage_for_role(record.triage, role, self.partner_peak(source, role))

    def partner_peak(
        self, path: str | Path, role: MeasurementRole | None
    ) -> float | None:
        """The peak of the faintest signal frame a background will be subtracted from.

        A lamp background belongs to its own lamp's signal; a sphere background
        to the sphere.  Absent a partner there is nothing to compare against,
        and the background is simply read as a background.

        The faintest partner, not the brightest: a lamp family is shot as a
        bright/dim pair, and the bright frame is deliberately at full scale.
        Judged against that, the leak threshold sits an order of magnitude
        above the dim series the same background is subtracted from, and a
        background carrying most of that dim signal passes as dark.
        """

        partner_role = BACKGROUND_PARTNERS.get(role)
        if partner_role is None:
            return None
        measurement = self.measurements.get(Path(path))
        lamp = measurement.lamp_family if measurement is not None else ""
        peaks = []
        for other in self._records(partner_role, lamp):
            record = self.loaded.get(other.path)
            if record is not None and record.triage.peak_counts is not None:
                peaks.append(float(record.triage.peak_counts))
        return min(peaks) if peaks else None

    def unconfirmed_suggestions(self) -> tuple[tuple[Path, MeasurementRole], ...]:
        """Files whose filename proposes one role that nobody has confirmed.

        A pre-filled control is help, not an assignment.  The bench states which
        files are in that in-between state so the surface can never show a role
        the campaign has not actually been given.
        """

        pending = []
        for path in self.loaded:
            if path in self.measurements:
                continue
            suggestion = self.observed.get(path)
            if suggestion is None or not suggestion.is_unambiguous:
                continue
            pending.append((path, suggestion.roles[0]))
        return tuple(sorted(pending, key=lambda item: item[0].name.casefold()))

    def unanimous_suggestions(
        self, *, declined: Iterable[str | Path] = ()
    ) -> tuple[tuple[Path, MeasurementRole], ...]:
        """Every unassigned file's role when the whole folder leaves nothing to ask.

        A suggestion is help until somebody confirms it, and confirming a
        folder one row at a time is the interaction the operator called
        painful.  When a *set* of filenames says exactly one thing — each file
        naming one role, a lamp role naming its lamp, and no two files claiming
        a role only one file can hold — there is nothing left for the operator
        to decide, and the bench may simply assign it.

        The test is on the set, not the file: one doubtful name anywhere makes
        the whole folder the operator's to sort out, and this returns nothing
        so the confirm flow stays exactly as it was.  Files in *declined* are
        ones the operator has deliberately unassigned; they neither block the
        set nor get their suggestion pushed back onto them.
        """

        ignored = {Path(item) for item in declined}
        claimed: set[tuple[MeasurementRole, str]] = set()
        for record in self.measurements.values():
            # What is already assigned is what a new suggestion has to fit
            # around: a folder holding a sphere already cannot gain a second.
            claimed.add((record.role, record.lamp_family))
        pending: list[tuple[Path, MeasurementRole]] = []
        for path in self.loaded:
            if path in self.measurements or path in ignored:
                continue
            suggestion = self.observed.get(path)
            if suggestion is None or not suggestion.is_unambiguous:
                return ()
            role = suggestion.roles[0]
            lamp = ""
            if role in LAMP_ROLES:
                try:
                    lamp = normalize_lamp_name(suggestion.lamp_name)
                except ValueError:
                    # "lamp" in a filename with no lamp name in it: which lamp
                    # is exactly the question only the operator can answer.
                    return ()
            if role in SINGLE_HOLDER_ROLES and (role, lamp) in claimed:
                return ()
            claimed.add((role, lamp))
            pending.append((path, role))
        return tuple(sorted(pending, key=lambda item: item[0].name.casefold()))

    def apply_unanimous_suggestions(
        self,
        *,
        declined: Iterable[str | Path] = (),
        saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
    ) -> tuple[MeasurementRecord, ...]:
        """Assign a folder's filename roles when they leave nothing to ask.

        Nothing is inferred here that :meth:`unanimous_suggestions` has not
        already found unambiguous, and every assignment goes through
        :meth:`classify_file` — the same door the operator's own combo uses —
        so an applied role is an ordinary assignment they can change.
        """

        return self._assign(
            self.unanimous_suggestions(declined=declined),
            saturation_level=saturation_level,
        )

    def confirm_suggested_roles(
        self, *, saturation_level: float = 0.98 * FULL_SCALE_COUNTS
    ) -> tuple[MeasurementRecord, ...]:
        """Assign every unconfirmed filename suggestion in one explicit step.

        This is still a manual act — the operator asks for it — but it makes
        confirming a whole folder one deliberate press instead of one popup pick
        per row, which is the interaction that silently produced no assignment
        at all.
        """

        return self._assign(
            self.unconfirmed_suggestions(), saturation_level=saturation_level
        )

    def _assign(
        self,
        pending: Sequence[tuple[Path, MeasurementRole]],
        *,
        saturation_level: float,
    ) -> tuple[MeasurementRecord, ...]:
        """Assign each (path, role) pair through the ordinary classification."""

        assigned = []
        for path, role in pending:
            suggestion = self.observed.get(path)
            lamp = suggestion.lamp_name if suggestion is not None else ""
            try:
                assigned.append(
                    self.classify_file(
                        path,
                        role,
                        lamp_family=lamp,
                        saturation_level=saturation_level,
                    )
                )
            except (FileNotFoundError, ValueError):
                # A lamp suggestion with no readable lamp name stays the
                # operator's to make; the rest of the folder still confirms.
                continue
        return tuple(assigned)

    @property
    def assigned_lamps(self) -> tuple[str, ...]:
        """Lamp names the operator actually assigned, in a stable order."""

        return tuple(
            sorted(
                {
                    record.lamp_family
                    for record in self.measurements.values()
                    if record.role in LAMP_ROLES and record.lamp_family
                },
                key=str.casefold,
            )
        )

    def lamp_for_frame(self, alignment: CalibrationBenchSession) -> str:
        """Return the lamp the frame currently open for fitting was assigned.

        A frame that carries no lamp role of its own falls back to the session's
        primary lamp, so a bench holding exactly one lamp keeps referencing it
        even while another file is on screen.
        """

        frame = alignment.frame
        if frame is not None:
            record = self.measurements.get(Path(frame.path))
            if record is not None and record.lamp_family:
                return record.lamp_family
        return self._primary_lamp()

    def lamp_pair(self, lamp: str) -> tuple[Path | None, Path | None]:
        """Return *lamp*'s assigned (signal, background) paths, as far as given.

        Line fitting belongs on the signal, and on the signal minus its own
        background when the pair is complete — the same thing ``echelle-align``
        measures.  Either half may be missing; the caller says what that means.
        """

        if not str(lamp).strip():
            return None, None
        signals = self._records(MeasurementRole.LAMP, lamp)
        backgrounds = self._records(MeasurementRole.LAMP_BACKGROUND, lamp)
        return (
            signals[-1].path if signals else None,
            backgrounds[-1].path if backgrounds else None,
        )

    def scope_alignment_to_lamp(
        self, alignment: CalibrationBenchSession
    ) -> LampReferenceSet:
        """Point the bench's click-to-fit at the open frame's own lamp catalog."""

        lamp = self.lamp_for_frame(alignment)
        current = alignment.reference
        if current is not None and current.lamp == lamp:
            return current
        reference = lamp_reference_set(lamp, alignment.lines)
        alignment.use_lamp_reference(reference)
        return reference

    def classify_file(
        self,
        path: str | Path,
        role: MeasurementRole,
        *,
        lamp_family: str = "",
        frame: BenchFrame | None = None,
        saturation_level: float = 0.98 * FULL_SCALE_COUNTS,
    ) -> MeasurementRecord:
        """Explicitly assign a role; filename suggestions never call this method."""

        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"measurement source not found: {source}")
        stat = source.stat()
        family = ""
        if role in LAMP_ROLES:
            family = normalize_lamp_name(lamp_family)
        loaded = self.loaded.get(source)
        exposure = (
            evaluate_exposure(frame, saturation_level=saturation_level)
            if frame is not None
            else (loaded.exposure if loaded is not None else None)
        )
        action = ""
        if exposure is not None and role in {
            MeasurementRole.SPHERE_BACKGROUND,
            MeasurementRole.LAMP_BACKGROUND,
        }:
            if exposure.state is ExposureState.SATURATED and role in LAMP_ROLES:
                # A lamp background is judged like its lamp: saturation there is
                # reported, never a reason to reshoot the pair.
                action = (
                    f"Saturated in {exposure.saturated_pixels} pixel(s) — expected "
                    "beside a dim-series lamp exposure; keep the paired signal at "
                    "this same exposure."
                )
            elif exposure.state is ExposureState.SATURATED:
                action = "Lower the exposure for both this background and its paired signal."
            else:
                action = (
                    "Keep the paired signal at this same exposure; acquire the signal next."
                )
        elif (
            exposure is not None
            and role is MeasurementRole.LAMP
            and exposure.state is ExposureState.SATURATED
        ):
            # The owner's own words: of course the dim series saturates the
            # strong lines — that is the point of shooting a bright/dim pair.
            action = (
                "Saturated strong lines are expected on a dim-series lamp frame; "
                "fit the unsaturated lines only. Every anchor is saturation-checked "
                "on its own detector window."
            )
        if action:
            assert exposure is not None
            exposure = ExposureGuidance(
                exposure.state,
                exposure.peak_value,
                exposure.saturated_pixels,
                exposure.saturated_fraction,
                exposure.exposure_s,
                exposure.suggested_exposure_s,
                action,
                exposure.anomalous_pixels,
            )
        existing = self.measurements.get(source)
        if (
            existing is not None
            and existing.role is role
            and existing.lamp_family == family
            and existing.size_bytes == int(stat.st_size)
            and existing.modified_ns == int(stat.st_mtime_ns)
        ):
            return existing
        record = MeasurementRecord(
            source,
            role,
            family,
            exposure,
            int(stat.st_size),
            int(stat.st_mtime_ns),
        )
        self.measurements[source] = record
        self.observed.setdefault(source, suggest_file_roles(source))
        self._invalidate_outputs()
        return record

    def remove_classification(self, path: str | Path) -> bool:
        existed = self.measurements.pop(Path(path), None) is not None
        if existed:
            self._invalidate_outputs()
        return existed

    def _invalidate_outputs(self) -> None:
        self.comparison = SphereComparison(ComparisonState.NOT_RUN, "inputs changed")
        self.toml_state = TomlState.NOT_GENERATED
        self.toml_paths = {}
        self.toml_snapshot_id = ""
        self.save_state = SaveState.NOT_READY
        self.saved_snapshot = None
        self.wavelength_correction = None
        self.last_error = ""

    def _records(self, role: MeasurementRole, family: str = "") -> tuple[MeasurementRecord, ...]:
        normalized = normalize_lamp_name(family) if family else ""
        return tuple(
            sorted(
                (
                    record
                    for record in self.measurements.values()
                    if record.role is role
                    and (not normalized or record.lamp_family == normalized)
                ),
                key=lambda record: record.path.name.casefold(),
            )
        )

    def _one(self, role: MeasurementRole) -> MeasurementRecord | None:
        records = self._records(role)
        return records[-1] if records else None

    def compute_sphere_comparison(
        self,
        calculator: Callable[..., AbsoluteCalibrationResult] = compute_absolute_calibration,
    ) -> SphereComparison:
        """Compute candidate factors and compare them, or report insufficient data."""

        sphere = self._one(MeasurementRole.SPHERE)
        background = self._one(MeasurementRole.SPHERE_BACKGROUND)
        if sphere is None or background is None:
            missing = [
                name
                for name, record in (
                    ("sphere signal", sphere),
                    ("sphere background", background),
                )
                if record is None
            ]
            reason = f"classify the {' and '.join(missing)} first"
            pending = [
                path.name
                for path, role in self.unconfirmed_suggestions()
                if role in {MeasurementRole.SPHERE, MeasurementRole.SPHERE_BACKGROUND}
            ]
            if pending:
                reason += (
                    f" — {', '.join(pending)} only shows a filename suggestion; "
                    "confirm the role in the Role column"
                )
            self.comparison = SphereComparison(ComparisonState.FAILED, reason)
            return self.comparison
        try:
            candidate = calculator(
                pattern=self.pattern_source,
                wavelength=self.wavelength_source,
                sphere=sphere.path,
                sphere_background=background.path,
                integral=self.integral_source,
            )
        except Exception as exc:  # numerical/IO boundary becomes recoverable state
            self.comparison = SphereComparison(ComparisonState.FAILED, str(exc))
            return self.comparison

        previous_paths = (self.previous_sphere, self.previous_sphere_background)
        if not all(path is not None and path.is_file() for path in previous_paths):
            self.comparison = SphereComparison(
                ComparisonState.INSUFFICIENT_DATA,
                "previous sphere pair is unavailable; candidate factors were computed",
                candidate=candidate,
            )
            return self.comparison
        try:
            previous = calculator(
                pattern=self.pattern_source,
                wavelength=self.wavelength_source,
                sphere=self.previous_sphere,
                sphere_background=self.previous_sphere_background,
                integral=self.integral_source,
            )
        except Exception as exc:
            self.comparison = SphereComparison(
                ComparisonState.INSUFFICIENT_DATA,
                f"previous factors could not be computed: {exc}",
                candidate=candidate,
            )
            return self.comparison

        count = min(candidate.factors_wmsr.size, previous.factors_wmsr.size)
        candidate_values = candidate.factors_wmsr[:count]
        previous_values = previous.factors_wmsr[:count]
        valid = (
            np.isfinite(candidate_values)
            & np.isfinite(previous_values)
            & (candidate_values > 0)
            & (previous_values > 0)
        )
        ratio = candidate_values[valid] / previous_values[valid]
        if ratio.size < 20:
            self.comparison = SphereComparison(
                ComparisonState.INSUFFICIENT_DATA,
                "fewer than 20 finite positive factor samples overlap",
                candidate,
                previous,
                int(ratio.size),
            )
            return self.comparison
        self.comparison = SphereComparison(
            ComparisonState.READY,
            "candidate and previous factors share a finite comparison range",
            candidate,
            previous,
            int(ratio.size),
            float(np.median(ratio)),
            float(np.percentile(ratio, 5)),
            float(np.percentile(ratio, 95)),
        )
        self._update_save_state(None)
        return self.comparison

    def checklist(self, alignment: CalibrationBenchSession) -> tuple[ChecklistItem, ...]:
        """Derive the procedure from the files and roles that actually exist.

        Nothing here is a fixed lamp list.  The rows follow what was loaded and
        assigned, every row that cannot be done yet names what would unblock
        it, and the previous campaign's lamps appear only as advice.
        """

        return (
            *self._input_items(),
            *self._sphere_items(),
            *self._lamp_items(),
            *self._output_items(alignment),
        )

    def _input_items(self) -> tuple[ChecklistItem, ...]:
        references_ready = all(
            path.is_file()
            for path in (self.pattern_source, self.wavelength_source, self.integral_source)
        )
        loaded = tuple(self.loaded.values())
        saturated = sum(
            1 for record in loaded if record.triage.state is ExposureState.SATURATED
        )
        # Saturation on a lamp frame is expected, so it is counted apart from
        # the saturation that actually refuses a frame its purpose.
        blocking = sum(
            1
            for path in self.loaded
            if (verdict := self.role_triage(path)) is not None
            and verdict.state is ExposureState.SATURATED
            and verdict.blocking
        )
        pending = self.unconfirmed_suggestions()
        on_bench = set(self.loaded) | set(self.measurements)
        detail = "drop SIF files on the bench or use Add files"
        if on_bench:
            expected = saturated - blocking
            saturation_note = f"{blocking} saturated"
            if expected:
                saturation_note += f" ({expected} saturated on purpose, on lamp frames)"
            detail = (
                f"{len(on_bench)} file(s) on the bench, {len(loaded)} file(s) triaged, "
                f"{saturation_note}"
            )
            if pending:
                detail += (
                    f"; {len(pending)} file(s) show only a filename suggestion and "
                    "carry no role yet"
                )
        return (
            ChecklistItem(
                "references",
                "Pattern, wavelength, and sphere reference",
                ChecklistState.DONE if references_ready else ChecklistState.ATTENTION,
                "loaded" if references_ready else "one or more reference files are missing",
                unblocked_by=""
                if references_ready
                else (
                    "these three tables are named when the bench is opened, not "
                    "picked from the files list: reopen the bench pointing at "
                    "tables that exist"
                ),
            ),
            ChecklistItem(
                "files",
                "SIFs loaded and exposure-triaged",
                ChecklistState.DONE if on_bench else ChecklistState.WAITING,
                detail,
                unblocked_by="" if on_bench else "drop any SIF onto the bench window",
            ),
        )

    def _sphere_items(self) -> tuple[ChecklistItem, ...]:
        items = []
        pending = dict(self.unconfirmed_suggestions())
        for role, key, label in (
            (MeasurementRole.SPHERE, "sphere", "Integrating-sphere signal"),
            (
                MeasurementRole.SPHERE_BACKGROUND,
                "sphere-background",
                "Integrating-sphere background",
            ),
        ):
            record = self._one(role)
            # A file the filename already points at is not the same as no file
            # at all, and saying so is what stops a pre-filled control from
            # reading as an assignment nobody made.
            suggested = [path.name for path, proposed in pending.items() if proposed is role]
            if record:
                detail = record.path.name
                unblocked = ""
            elif suggested:
                detail = (
                    f"{', '.join(suggested)} is only suggested by its filename; "
                    "no file carries this role yet"
                )
                unblocked = (
                    f"confirm the {label.casefold()} role on {suggested[0]} — pick it "
                    "in the Role column or press Confirm suggested roles"
                )
            else:
                detail = "no file carries this role yet"
                unblocked = f"assign the {label.casefold()} role to a loaded file"
            items.append(
                ChecklistItem(
                    key,
                    label,
                    ChecklistState.DONE if record else ChecklistState.WAITING,
                    detail,
                    unblocked_by=unblocked,
                )
            )
        comparison_done = self.comparison.state in {
            ComparisonState.READY,
            ComparisonState.INSUFFICIENT_DATA,
        }
        pair_ready = (
            self._one(MeasurementRole.SPHERE) is not None
            and self._one(MeasurementRole.SPHERE_BACKGROUND) is not None
        )
        items.append(
            ChecklistItem(
                "sphere-comparison",
                "Absolute factors compared with previous campaign",
                ChecklistState.DONE
                if comparison_done
                else (
                    ChecklistState.ATTENTION
                    if self.comparison.state is ComparisonState.FAILED
                    else ChecklistState.WAITING
                ),
                self.comparison.reason,
                unblocked_by=""
                if comparison_done
                else (
                    "press Measure sensitivity — the sphere pair is enough"
                    if pair_ready
                    else "assign the sphere signal and its background; no lamp is needed"
                ),
            )
        )
        return tuple(items)

    def _lamp_items(self) -> tuple[ChecklistItem, ...]:
        items = []
        assigned = self.assigned_lamps
        if not assigned:
            items.append(
                ChecklistItem(
                    "lamp-any",
                    "At least one lamp signal",
                    ChecklistState.WAITING,
                    "no lamp has been assigned yet",
                    unblocked_by="give any loaded file a lamp role; any lamp name works",
                )
            )
        for lamp in assigned:
            for role, suffix, label in (
                (MeasurementRole.LAMP, "signal", "signal"),
                (MeasurementRole.LAMP_BACKGROUND, "background", "background"),
            ):
                records = self._records(role, lamp)
                items.append(
                    ChecklistItem(
                        f"lamp-{lamp}-{suffix}",
                        f"{lamp} lamp {label}",
                        ChecklistState.DONE if records else ChecklistState.WAITING,
                        ", ".join(record.path.name for record in records)
                        if records
                        else f"no file carries the {lamp} {label} role yet",
                        unblocked_by=""
                        if records
                        else f"assign the {lamp} lamp {label} role to a loaded file",
                    )
                )
        items.append(self._lamp_suggestion_item(assigned))
        return tuple(items)

    def _lamp_suggestion_item(self, assigned: Sequence[str]) -> ChecklistItem:
        """Advice from the previous campaign, which never blocks anything."""

        folded = {lamp.casefold() for lamp in assigned}
        previous = ", ".join(self.suggested_lamps) or "nothing recorded"
        consider = [
            lamp
            for lamp in (*self.suggested_lamps, *KNOWN_LAMP_NAMES)
            if lamp.casefold() not in folded
        ]
        consider = list(dict.fromkeys(consider))
        detail = f"last time: {previous}"
        if consider:
            detail += f"; consider {', '.join(consider)} if available"
        if assigned:
            detail += f". This session measured {', '.join(assigned)}"
        return ChecklistItem(
            "lamp-suggestions",
            "More lamps are a gift, never a gate",
            ChecklistState.SUGGESTION,
            detail,
            blocking=False,
        )

    def _alignment_item(self, alignment: CalibrationBenchSession) -> ChecklistItem:
        """Report the fit together with the catalog that actually anchored it."""

        reference = alignment.reference
        count = len(alignment.anchors)
        anchors = f"{count} anchor" + ("" if count == 1 else "s")
        rms = "" if alignment.rms_px is None else f", RMS {alignment.rms_px:.2f} px"
        if reference is not None and not reference.is_referenceable:
            return ChecklistItem(
                "alignment",
                "Lamp alignment solved and reviewed",
                ChecklistState.ATTENTION,
                f"{anchors}; {reference.message}",
                unblocked_by=(
                    f"assign one of {', '.join(KNOWN_LAMP_NAMES)} to this frame, or "
                    "add this lamp's rows to the curated wavelength table"
                ),
            )
        if reference is not None and reference.catalog_label:
            detail = f"{anchors} vs {reference.catalog_label} catalog{rms}"
        else:
            detail = (
                f"{anchors} against the whole curated table{rms} — "
                "no lamp catalog is scoping the fit yet"
            )
        done = alignment.alignment_state is AlignmentState.ALIGNED
        if self.assigned_lamps:
            next_step = "load a lamp frame and click two known lines in different orders"
        else:
            next_step = "assign a lamp role to any loaded file, then click two lines"
        if not done:
            return ChecklistItem(
                "alignment",
                "Lamp alignment solved and reviewed",
                ChecklistState.WAITING,
                detail,
                unblocked_by=next_step,
            )
        # A solved fit is not a validated one.  RMS says the anchors agree with
        # each other in pixels; the BH paper's standard was agreement with
        # Fulcher-alpha in nanometres, so the row states which of the two it
        # has (owner, F19 second rider) rather than going green on the easier
        # number and letting the harder one go unasked.
        validation = alignment.validate_science_lines()
        return ChecklistItem(
            "alignment",
            "Lamp alignment solved and reviewed",
            ChecklistState.DONE,
            f"{detail}; {validation.message}",
            unblocked_by=""
            if validation.measured
            else "carry this to first plasma data and validate against Fulcher there",
        )

    def _output_items(self, alignment: CalibrationBenchSession) -> tuple[ChecklistItem, ...]:
        return (
            self._alignment_item(alignment),
            ChecklistItem(
                "tomls",
                "Commented campaign TOMLs generated",
                ChecklistState.DONE
                if self.toml_state is TomlState.GENERATED
                else (
                    ChecklistState.ATTENTION
                    if self.toml_state is TomlState.FAILED
                    else ChecklistState.WAITING
                ),
                ", ".join(path.name for path in self.toml_paths.values())
                if self.toml_paths
                else "generate once the sphere pair, one lamp pair, and the fit are in",
                unblocked_by=""
                if self.toml_state is TomlState.GENERATED
                else "complete the rows above, then press Save alignment settings",
            ),
            ChecklistItem(
                "snapshot",
                "Snapshot saved and validated",
                ChecklistState.DONE
                if self.save_state is SaveState.VALIDATED
                else (
                    ChecklistState.ATTENTION
                    if self.save_state is SaveState.FAILED
                    else ChecklistState.WAITING
                ),
                self._saved_snapshot_detail(),
                unblocked_by=""
                if self.save_state is SaveState.VALIDATED
                else "press Save alignment settings for this snapshot identity, then save and validate",
            ),
        )

    def _saved_snapshot_detail(self) -> str:
        if self.saved_snapshot is None:
            return self.last_error or "save only after every prerequisite is explicit"
        if self.wavelength_correction is None:
            return self.saved_snapshot.snapshot_id
        return f"{self.saved_snapshot.snapshot_id} — {self.wavelength_correction.reason}"

    def _complete_lamps(self) -> tuple[str, ...]:
        """Lamps that carry both a signal and a background this session."""

        return tuple(
            lamp
            for lamp in self.assigned_lamps
            if self._records(MeasurementRole.LAMP, lamp)
            and self._records(MeasurementRole.LAMP_BACKGROUND, lamp)
        )

    def _primary_lamp(self) -> str:
        complete = self._complete_lamps()
        if complete:
            return complete[0]
        assigned = self.assigned_lamps
        return assigned[0] if assigned else ""

    def _measurement_pairs_ready(self) -> bool:
        """One sphere pair and one complete lamp pair are the whole demand."""

        if self._one(MeasurementRole.SPHERE) is None:
            return False
        if self._one(MeasurementRole.SPHERE_BACKGROUND) is None:
            return False
        return bool(self._complete_lamps())

    def _composition_ready(self, alignment: CalibrationBenchSession) -> bool:
        return (
            self._measurement_pairs_ready()
            and self.comparison.state
            in {ComparisonState.READY, ComparisonState.INSUFFICIENT_DATA}
            and alignment.alignment_state is AlignmentState.ALIGNED
            and alignment.transform is not None
        )

    def compose_tomls(
        self,
        snapshot_id: str,
        alignment: CalibrationBenchSession,
        *,
        snapshot_root: str | Path | None = None,
    ) -> dict[str, str]:
        """Compose commented ordinary TOML from the measured session state.

        ``snapshot_root`` is where the snapshot for this identity will be
        written.  The export configuration points at the sphere pair the way the
        snapshot itself does — back out to the calibration folder — so knowing
        the folder makes those two paths exact instead of assumed.
        """

        if not self._composition_ready(alignment):
            raise SnapshotError(
                "sphere/lamp pairs, sphere-factor result, and aligned anchors are required"
            )
        sphere = self._one(MeasurementRole.SPHERE)
        sphere_background = self._one(MeasurementRole.SPHERE_BACKGROUND)
        assert sphere is not None and sphere_background is not None
        lamp_records = tuple(
            record
            for lamp in self.assigned_lamps
            for role in LAMP_ROLES
            for record in self._records(role, lamp)
        )
        transform = alignment.transform
        assert transform is not None

        campaign_lines = [
            "# Generated by echelle-calib from explicitly classified measurements.",
            "# Review freely; this is ordinary authoritative TOML.",
            'schema = "echelle-calibration-campaign/v1"',
            f"snapshot_id = {_toml_string(snapshot_id)}",
            f"comparison_state = {_toml_string(self.comparison.state.value)}",
            f"comparison_reason = {_toml_string(self.comparison.reason)}",
        ]
        if self.comparison.median_ratio is not None:
            campaign_lines.extend(
                [
                    f"factor_median_ratio = {self.comparison.median_ratio:.12g}",
                    f"factor_p05_ratio = {self.comparison.p05_ratio:.12g}",
                    f"factor_p95_ratio = {self.comparison.p95_ratio:.12g}",
                    f"factor_sample_count = {self.comparison.sample_count}",
                ]
            )
        for record in (sphere, sphere_background, *lamp_records):
            campaign_lines.extend(
                [
                    "",
                    "[[measurements]]",
                    f"role = {_toml_string(record.role.value)}",
                    f"source_name = {_toml_string(record.path.name)}",
                ]
            )
            if record.lamp_family:
                campaign_lines.append(f"lamp_family = {_toml_string(record.lamp_family)}")
            if record.exposure and record.exposure.exposure_s is not None:
                campaign_lines.append(f"exposure_s = {record.exposure.exposure_s:.12g}")
            if record.exposure and record.exposure.peak_value is not None:
                campaign_lines.append(f"raw_peak_counts = {record.exposure.peak_value:.12g}")
                campaign_lines.append(
                    f"saturated_pixels = {record.exposure.saturated_pixels}"
                )

        alignment_lines = [
            "# Generated alignment settings; source files are named without machine paths.",
            "# Review anchors and residuals before using this file.",
            'schema = "echelle-calibration-alignment/v1"',
            f"snapshot_id = {_toml_string(snapshot_id)}",
            f"base_pattern_file = {_toml_string(self.pattern_source.name)}",
            f"base_wavelength_file = {_toml_string(self.wavelength_source.name)}",
            f"sphere_file = {_toml_string(sphere.path.name)}",
            f"sphere_background_file = {_toml_string(sphere_background.path.name)}",
            f"lamps = [{', '.join(_toml_string(item) for item in self.assigned_lamps)}]",
            f"n_lines = {len(alignment.anchors)}",
            f"rms_px = {alignment.rms_px:.12g}",
            "",
            "[transform]",
            f"dx_px = {transform.dx_px:.12g}",
            f"dy_px = {transform.dy_px:.12g}",
            f"theta_rad = {transform.theta_rad:.12g}",
        ]
        for anchor in alignment.anchor_rows():
            alignment_lines.extend(
                [
                    "",
                    "[[anchors]]",
                    f"order = {anchor.line.order_idx}",
                    f"wavelength_nm = {anchor.line.wavelength_nm:.12g}",
                    f"expected_pixel = {anchor.line.center_pixel:.12g}",
                    f"measured_pixel = {anchor.fit.center_pixel:.12g}",
                    f"snr = {anchor.fit.snr:.12g}",
                ]
            )

        inherited_note = (
            f"Inherited from {INHERITED_EXPORT_CONFIG_ID}; not measured by this "
            "bench session. Review before the next LHD campaign."
        )
        # The snapshot references the sphere frames rather than copying them, so
        # this configuration names them the same way: out of the snapshot folder
        # and back into the calibration folder that holds the light.
        sphere_reference = _snapshot_reference(snapshot_root, sphere.path)
        sphere_background_reference = _snapshot_reference(
            snapshot_root, sphere_background.path
        )
        export_lines = [
            "# Generated SpectroCube export configuration for this snapshot.",
            "# Paths are relative to the snapshot folder and remain hand-editable.",
            "[metadata]",
            f"config_id = {_toml_string(snapshot_id)}",
            "",
            "# Timing and crop the bench does not measure. Cubes need these to",
            "# write the frozen LHD text header, so the previous campaign's",
            "# values are carried forward and marked as inherited.",
            f"trigger_delay_s = {INHERITED_TRIGGER_DELAY_S:.12g}",
            f"time_axis_reference = {_toml_string(INHERITED_TIME_AXIS_REFERENCE)}",
            f"frame_time_formula = {_toml_string(INHERITED_FRAME_TIME_FORMULA)}",
            f"trigger_delay_note = {_toml_string(inherited_note)}",
            f"crop_measured_at = {_toml_string(INHERITED_CROP_MEASURED_AT)}",
            f"crop_measurement_note = {_toml_string(inherited_note)}",
            "",
            "[calibration]",
            'camera = "CMOS"',
            'instrument_id = "echelle"',
            'wavelength_medium = "air"',
            f"calibration_dir = {_toml_string(snapshot_id)}",
            'order_pattern = "pattern.txt"',
            'wavelength = "wavelength.txt"',
            f"sphere = {_toml_string(sphere_reference)}",
            f"sphere_background = {_toml_string(sphere_background_reference)}",
            'integral = "integral.txt"',
            "",
            "[export]",
            'units = "wmsr"',
            'output_suffix = "_spectrocube_wmsr"',
            "drop_nonfinite_columns = true",
            f"wavelength_min_nm = {INHERITED_WAVELENGTH_MIN_NM:.12g}",
            f"calibration_source = {_toml_string('snapshot ' + snapshot_id)}",
        ]
        return {
            "campaign": "\n".join(campaign_lines) + "\n",
            "alignment": "\n".join(alignment_lines) + "\n",
            "export": "\n".join(export_lines) + "\n",
        }

    def write_tomls(
        self,
        destination_root: str | Path,
        snapshot_id: str,
        alignment: CalibrationBenchSession,
        *,
        overwrite: bool = False,
        snapshot_root: str | Path | None = None,
    ) -> dict[str, Path]:
        """Atomically publish a new identity's generated TOML bundle.

        Refusing to clobber is the default and stays the default.  ``overwrite``
        is the deliberate second press: the new bundle is still staged and
        parsed in full before anything existing is touched, so a failure part
        way through leaves the old files exactly where they were.
        """

        destination_parent = Path(destination_root)
        destination_parent.mkdir(parents=True, exist_ok=True)
        destination = destination_parent / snapshot_id
        if destination.exists() and not overwrite:
            self.toml_state = TomlState.FAILED
            self.last_error = f"configuration identity already exists: {destination}"
            raise SnapshotError(self.last_error)
        superseded_parent = None
        rescued = None
        try:
            texts = self.compose_tomls(snapshot_id, alignment, snapshot_root=snapshot_root)
            staging_parent = Path(
                tempfile.mkdtemp(prefix=f".{snapshot_id}.configs-", dir=destination_parent)
            )
            staging = staging_parent / snapshot_id
            staging.mkdir()
            paths = {}
            for name, text in texts.items():
                path = staging / f"{name}.toml"
                path.write_text(text, encoding="utf-8")
                with path.open("rb") as stream:
                    tomllib.load(stream)
                paths[name] = path
            if destination.exists():
                # Every file is written and parsed by now, so the old bundle is
                # only moved aside once the new one is known to be good — and it
                # is parked in its own directory beside the destination, never
                # inside the staging tree, because the cleanup that follows a
                # failed publish deletes that tree and the parked bundle is at
                # that moment the only copy of the old files.
                superseded_parent = Path(
                    tempfile.mkdtemp(
                        prefix=f".{snapshot_id}.superseded-", dir=destination_parent
                    )
                )
                superseded = superseded_parent / snapshot_id
                os.replace(destination, superseded)
                try:
                    os.replace(staging, destination)
                except Exception:
                    rescued = superseded.resolve()
                    os.replace(superseded, destination)
                    rescued = None
                    raise
            else:
                os.replace(staging, destination)
            shutil.rmtree(staging_parent, ignore_errors=True)
            self._discard_parking(superseded_parent)
        except Exception as exc:
            if "staging_parent" in locals():
                shutil.rmtree(staging_parent, ignore_errors=True)
            self.toml_state = TomlState.FAILED
            if rescued is not None:
                # The old bundle could not be put back, so it stays parked and
                # nothing here may delete it; the operator is told where it is.
                self.last_error = (
                    f"{destination} could not be written and its previous files "
                    f"could not be put back: {exc}. The old bundle is kept at {rescued}"
                )
                raise SnapshotError(self.last_error) from exc
            self._discard_parking(superseded_parent)
            self.last_error = str(exc)
            raise
        self.toml_paths = {
            name: destination / path.name for name, path in paths.items()
        }
        self.toml_snapshot_id = snapshot_id
        self.toml_state = TomlState.GENERATED
        self.last_error = ""
        self._update_save_state(alignment)
        return dict(self.toml_paths)

    @staticmethod
    def _discard_parking(parent: Path | None) -> None:
        """Remove a superseded-bundle parking directory that is no longer needed."""

        if parent is not None:
            shutil.rmtree(parent, ignore_errors=True)

    def _update_save_state(self, alignment: CalibrationBenchSession | None) -> None:
        if self.save_state is SaveState.VALIDATED:
            return
        ready = self.toml_state is TomlState.GENERATED
        if alignment is not None:
            ready = ready and self._composition_ready(alignment)
        self.save_state = SaveState.READY if ready else SaveState.NOT_READY

    def ready_for_snapshot(
        self, snapshot_id: str, alignment: CalibrationBenchSession
    ) -> bool:
        """Return whether measured state and generated TOMLs match this identity."""

        self._update_save_state(alignment)
        return (
            self.save_state is SaveState.READY
            and bool(snapshot_id)
            and self.toml_snapshot_id == snapshot_id
        )

    def alignment_settings(
        self,
        snapshot_id: str,
        detector: str,
        alignment: CalibrationBenchSession,
        *,
        notes: str = "",
        created_at: str = "",
    ) -> AlignmentSettings:
        """Summarise the solved alignment in the established settings shape."""

        if alignment.transform is None:
            raise SnapshotError("no rigid transform has been solved yet")
        sphere = self._one(MeasurementRole.SPHERE)
        sphere_background = self._one(MeasurementRole.SPHERE_BACKGROUND)
        primary = self._primary_lamp()
        lamp = self._records(MeasurementRole.LAMP, primary) if primary else ()
        lamp_background = (
            self._records(MeasurementRole.LAMP_BACKGROUND, primary) if primary else ()
        )
        return AlignmentSettings(
            instrument_id=str(detector).strip().lower(),
            base_wavelength_file=self.wavelength_source.name,
            base_pattern_file=self.pattern_source.name,
            transform=alignment.transform,
            n_lines=len(alignment.anchors),
            rms_px=float(alignment.rms_px or 0.0),
            created_at=created_at,
            alignment_dataset_id=snapshot_id,
            alignment_source_dir=lamp[-1].path.parent.name if lamp else "",
            alignment_lamp=", ".join(self.assigned_lamps),
            signal_file=lamp[-1].path.name if lamp else "",
            background_file=lamp_background[-1].path.name if lamp_background else "",
            sphere_file=sphere.path.name if sphere is not None else "",
            sphere_background_file=(
                sphere_background.path.name if sphere_background is not None else ""
            ),
            output_wavelength_file="wavelength.txt",
            notes=notes or "Rigid detector correction solved on the live calibration bench.",
        )

    @staticmethod
    def _table_provenance(settings: AlignmentSettings) -> list[tuple[str, str]]:
        transform = settings.transform
        return [
            ("Generated", settings.created_at),
            ("Base wavelength file", settings.base_wavelength_file),
            ("Base pattern file", settings.base_pattern_file),
            ("Alignment dataset", settings.alignment_dataset_id),
            ("Alignment source dir", settings.alignment_source_dir),
            ("Lamp", settings.alignment_lamp),
            ("Signal", settings.signal_file),
            ("Background", settings.background_file),
            ("Sphere", settings.sphere_file),
            ("Sphere background", settings.sphere_background_file),
            ("Correction model", "rigid detector transform, dx/dy/theta"),
            (
                "Transform",
                f"dx {transform.dx_px:+.4f} px, dy {transform.dy_px:+.4f} px, "
                f"theta {transform.theta_deg:+.5f} deg",
            ),
            ("RMS", f"{settings.rms_px:.4f} px"),
            ("Fitted lines", str(settings.n_lines)),
            ("Settings file", ALIGNMENT_SETTINGS_FILENAME),
            ("Note", settings.notes),
        ]

    def _vetting_record(self) -> dict[str, object]:
        """Which vetted line set anchored this snapshot, read from its table.

        A snapshot that records only its RMS says how self-consistent the fit
        was.  Which lines were trusted, and on whose authority, is what lets a
        later reader decide whether to believe it — so the manifest carries the
        vetted set by name, and carries its absence by name too.
        """

        vetting = table_vetting(self.wavelength_source)
        return {
            "vetted_set": vetting.vetted_set,
            "vetted_set_source": vetting.vetted_table,
            "vetted_lineage": list(vetting.lineage),
        }

    @staticmethod
    def _validation_record(alignment: CalibrationBenchSession) -> dict[str, object]:
        """Both numbers, never just the flattering one.

        ``rms_px`` beside it is anchor self-consistency.  This says whether the
        solution was ever held against the lines physics knows, and when it was
        not, the manifest says which — so a later reader never mistakes an
        unvalidated snapshot for a validated one.
        """

        validation = alignment.validate_science_lines()
        record: dict[str, object] = {
            "science_validation": validation.state.value,
            "science_validation_note": validation.message,
        }
        if validation.measured:
            record["science_lines_validated"] = validation.line_count
            record["science_residual_rms_nm"] = validation.rms_residual_nm
            record["science_residual_median_nm"] = validation.median_residual_nm
        return record

    def save_snapshot(
        self,
        destination_root: str | Path,
        *,
        snapshot_id: str,
        detector: str,
        alignment: CalibrationBenchSession,
        notes: str = "",
        base_snapshot: str | None = None,
        validity: Mapping[str, object] | None = None,
    ) -> Snapshot:
        """Create through Packet 0's API, then validate through its validator.

        The saved ``wavelength.txt`` is the base table with the solved rigid
        transform applied, so the snapshot carries the calibration the bench
        actually measured rather than the table it started from.
        """

        self._update_save_state(alignment)
        if not self.ready_for_snapshot(snapshot_id, alignment):
            self.save_state = SaveState.FAILED
            if self.toml_snapshot_id and self.toml_snapshot_id != snapshot_id:
                self.last_error = (
                    f"generated TOMLs target {self.toml_snapshot_id!r}, not {snapshot_id!r}"
                )
                raise SnapshotError(self.last_error)
            raise SnapshotError("campaign is not ready: generate TOMLs after all measurements")
        sphere = self._one(MeasurementRole.SPHERE)
        sphere_background = self._one(MeasurementRole.SPHERE_BACKGROUND)
        assert sphere is not None and sphere_background is not None
        lamp_files = [
            (record.lamp_family, record.path)
            for record in self.measurements.values()
            if record.role in LAMP_ROLES
        ]
        self.save_state = SaveState.SAVING
        staging = Path(tempfile.mkdtemp(prefix=f".{snapshot_id}.wavelength-"))
        try:
            epoch = dict(validity) if validity else default_validity(snapshot_id)
            settings = self.alignment_settings(
                snapshot_id,
                detector,
                alignment,
                notes=notes,
                created_at=str(epoch.get("date_from", "")),
            )
            corrected = staging / self.wavelength_source.name
            correction = write_corrected_wavelength_table(
                self.wavelength_source,
                corrected,
                pattern=alignment.pattern,
                transform=alignment.transform,
                metadata=self._table_provenance(settings),
            )
            snapshot = create_snapshot(
                destination_root,
                snapshot_id=snapshot_id,
                detector=detector,
                files={
                    "pattern": self.pattern_source,
                    "wavelength": corrected,
                    "sphere": sphere.path,
                    "sphere_background": sphere_background.path,
                    "integral": self.integral_source,
                },
                lamps=self.assigned_lamps,
                lamp_files=lamp_files,
                # The calibration folder already holds the lamp and sphere
                # frames; the snapshot records where they are and what they
                # hash to, and copies none of their bytes.
                reference_raw=True,
                notes=notes,
                base_snapshot=base_snapshot,
                validity=epoch,
                alignment={
                    "dx_px": alignment.transform.dx_px,
                    "dy_px": alignment.transform.dy_px,
                    "rotation_deg": alignment.transform.theta_deg,
                    "rms_px": alignment.rms_px,
                    "wavelength_correction_applied": correction.applied,
                    "wavelength_max_shift_px": correction.max_shift_px,
                    **self._vetting_record(),
                    **self._validation_record(alignment),
                },
                qc={
                    "lines_used": len(alignment.anchors),
                    "worst_residual_px": max(
                        (item.magnitude_px for item in alignment.residuals), default=0.0
                    ),
                    "sphere_comparison": self.comparison.state.value,
                },
            )
            save_alignment_settings(settings, snapshot.root / ALIGNMENT_SETTINGS_FILENAME)
            validated = load_snapshot(snapshot.root)
        except Exception as exc:
            self.save_state = SaveState.FAILED
            self.last_error = str(exc)
            raise
        finally:
            shutil.rmtree(staging, ignore_errors=True)
        self.saved_snapshot = validated
        self.wavelength_correction = correction
        self.save_state = SaveState.VALIDATED
        self.last_error = ""
        return validated
