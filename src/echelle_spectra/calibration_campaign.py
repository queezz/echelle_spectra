"""UI-independent campaign memory for the live calibration bench.

The live Qt window is deliberately a thin adapter over this module.  Explicit
file-role classification, the self-ticking procedure, exposure guidance,
shared-catalog line help, sphere-factor comparison, TOML composition, and
snapshot save/validation can therefore be rehearsed without an event loop.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from enum import Enum
from pathlib import Path

import numpy as np

try:  # pragma: no cover - selected by the running Python version
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

from .calibration_bench import AlignmentState, BenchFrame, CalibrationBenchSession
from .snapshot import Snapshot, SnapshotError, create_snapshot, load_snapshot
from .tools.calibration_alignment import (
    AlignmentSettings,
    CalibrationTableLine,
    RigidTransform,
    apply_rigid_correction_to_lines,
    load_wavelength_table,
    save_alignment_settings,
    write_wavelength_table,
)
from .tools.line_catalog import SpectralLine, load_line_table

__all__ = [
    "ALIGNMENT_SETTINGS_FILENAME",
    "AbsoluteCalibrationResult",
    "CalibrationCampaignSession",
    "CatalogOrderLine",
    "ChecklistItem",
    "ChecklistState",
    "ComparisonState",
    "ExposureGuidance",
    "ExposureState",
    "FileRoleSuggestion",
    "MeasurementRecord",
    "MeasurementRole",
    "SaveState",
    "TomlState",
    "WavelengthCorrection",
    "catalog_lines_for_order",
    "compute_absolute_calibration",
    "default_validity",
    "evaluate_exposure",
    "suggest_file_roles",
    "write_corrected_wavelength_table",
]

ALIGNMENT_SETTINGS_FILENAME = "alignment.toml"

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


class TomlState(Enum):
    """Generated configuration lifecycle."""

    NOT_GENERATED = "not-generated"
    GENERATED = "generated"
    FAILED = "failed"


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
    """One self-ticking procedure row."""

    key: str
    label: str
    state: ChecklistState
    detail: str


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


def _normalized_lamp_family(value: str) -> str:
    compact = value.strip().lower().replace("-", "").replace(" ", "")
    aliases = {
        "thar": "ThAr",
        "ne": "Ne",
        "hg": "Hg",
        "h2": "H2",
        "fulcher": "H2",
    }
    if compact not in aliases:
        raise ValueError("lamp family must be one of ThAr, Ne, Hg, or H2")
    return aliases[compact]


def _catalog_family(lamp_family: str) -> str:
    return {"ThAr": "thar", "Ne": "ne", "Hg": "hg", "H2": "fulcher"}[
        _normalized_lamp_family(lamp_family)
    ]


def suggest_file_roles(path: str | Path) -> FileRoleSuggestion:
    """Suggest likely roles from a filename without accepting the suggestion."""

    stem = Path(path).stem.casefold()
    background = any(token in stem for token in ("background", "_bg", "-bg", "bkg"))
    sphere = any(token in stem for token in ("sphere", "sphr", "absolute"))
    lamp = any(token in stem for token in ("thar", "th-ar", "neon", "lamp", "h2", "hg"))
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
            "filename looks like a lamp background; confirm the role and lamp family",
        )
    if lamp:
        return FileRoleSuggestion(
            (MeasurementRole.LAMP,),
            "filename looks like a lamp signal; confirm the role and lamp family",
        )
    if background:
        return FileRoleSuggestion(
            (MeasurementRole.SPHERE_BACKGROUND, MeasurementRole.LAMP_BACKGROUND),
            "background is ambiguous until its source is selected explicitly",
        )
    return FileRoleSuggestion(
        tuple(MeasurementRole),
        "filename does not support a unique role; select it explicitly",
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


def evaluate_exposure(
    frame: BenchFrame,
    *,
    saturation_level: float = 0.98 * 65535,
    dim_fraction: float = 0.20,
    target_fraction: float = 0.70,
) -> ExposureGuidance:
    """Inspect raw detector pixels and state the next safe exposure action."""

    values = np.asarray(frame.images, dtype=float)
    finite = values[np.isfinite(values)]
    exposure_s = _metadata_exposure_s(frame.metadata)
    if not finite.size:
        return ExposureGuidance(
            ExposureState.NO_DATA,
            None,
            0,
            0.0,
            exposure_s,
            None,
            "No finite raw pixels are available; reacquire before continuing.",
        )
    peak = float(np.max(finite))
    saturated_pixels = int(np.count_nonzero(finite >= saturation_level))
    saturated_fraction = saturated_pixels / int(finite.size)
    scale = target_fraction * saturation_level / max(peak, np.finfo(float).eps)
    suggested = exposure_s * scale if exposure_s is not None else None
    if saturated_pixels:
        action = "Lower exposure and reacquire; do not accept anchors from this frame."
        if suggested is not None:
            action = f"Lower exposure to about {suggested:.4g} s and reacquire."
        return ExposureGuidance(
            ExposureState.SATURATED,
            peak,
            saturated_pixels,
            saturated_fraction,
            exposure_s,
            suggested,
            action,
        )
    if peak < dim_fraction * saturation_level:
        action = "Increase exposure, then reacquire for stronger unsaturated lines."
        if suggested is not None:
            action = f"Increase exposure toward {suggested:.4g} s, then reacquire."
        return ExposureGuidance(
            ExposureState.DIM,
            peak,
            0,
            0.0,
            exposure_s,
            suggested,
            action,
        )
    return ExposureGuidance(
        ExposureState.GOOD,
        peak,
        0,
        0.0,
        exposure_s,
        exposure_s,
        "Exposure is usable. Continue with line identification and anchor fitting.",
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
    catalog = [
        line
        for line in load_line_table(_catalog_family(lamp_family))
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


class CalibrationCampaignSession:
    """Campaign-memory state transitions independent of Qt."""

    def __init__(
        self,
        *,
        pattern_source: str | Path,
        wavelength_source: str | Path,
        integral_source: str | Path,
        required_lamps: Sequence[str] = ("ThAr",),
        previous_sphere: str | Path | None = None,
        previous_sphere_background: str | Path | None = None,
    ) -> None:
        self.pattern_source = Path(pattern_source)
        self.wavelength_source = Path(wavelength_source)
        self.integral_source = Path(integral_source)
        self.required_lamps = tuple(
            dict.fromkeys(_normalized_lamp_family(item) for item in required_lamps)
        )
        if not self.required_lamps:
            raise ValueError("at least one required lamp family is needed")
        self.previous_sphere = Path(previous_sphere) if previous_sphere else None
        self.previous_sphere_background = (
            Path(previous_sphere_background) if previous_sphere_background else None
        )
        self.observed: dict[Path, FileRoleSuggestion] = {}
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

    def classify_file(
        self,
        path: str | Path,
        role: MeasurementRole,
        *,
        lamp_family: str = "",
        frame: BenchFrame | None = None,
        saturation_level: float = 0.98 * 65535,
    ) -> MeasurementRecord:
        """Explicitly assign a role; filename suggestions never call this method."""

        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"measurement source not found: {source}")
        stat = source.stat()
        family = ""
        if role in {MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND}:
            family = _normalized_lamp_family(lamp_family)
        exposure = (
            evaluate_exposure(frame, saturation_level=saturation_level)
            if frame is not None
            else None
        )
        if exposure is not None and role in {
            MeasurementRole.SPHERE_BACKGROUND,
            MeasurementRole.LAMP_BACKGROUND,
        }:
            if exposure.state is ExposureState.SATURATED:
                action = "Lower the exposure for both this background and its paired signal."
            else:
                action = (
                    "Keep the paired signal at this same exposure; acquire the signal next."
                )
            exposure = ExposureGuidance(
                exposure.state,
                exposure.peak_value,
                exposure.saturated_pixels,
                exposure.saturated_fraction,
                exposure.exposure_s,
                exposure.suggested_exposure_s,
                action,
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
        normalized = _normalized_lamp_family(family) if family else ""
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
            self.comparison = SphereComparison(
                ComparisonState.FAILED,
                "classify both sphere signal and sphere background first",
            )
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
        """Derive the ordered, self-ticking measurement procedure."""

        items = []
        references_ready = all(
            path.is_file()
            for path in (self.pattern_source, self.wavelength_source, self.integral_source)
        )
        items.append(
            ChecklistItem(
                "references",
                "Pattern, wavelength, and sphere reference",
                ChecklistState.DONE if references_ready else ChecklistState.ATTENTION,
                "loaded" if references_ready else "one or more reference files are missing",
            )
        )
        for role, key, label in (
            (MeasurementRole.SPHERE, "sphere", "Integrating-sphere signal"),
            (
                MeasurementRole.SPHERE_BACKGROUND,
                "sphere-background",
                "Integrating-sphere background",
            ),
        ):
            record = self._one(role)
            items.append(
                ChecklistItem(
                    key,
                    label,
                    ChecklistState.DONE if record else ChecklistState.WAITING,
                    record.path.name if record else "classify a measured SIF explicitly",
                )
            )
        comparison_done = self.comparison.state in {
            ComparisonState.READY,
            ComparisonState.INSUFFICIENT_DATA,
        }
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
            )
        )
        for family in self.required_lamps:
            for role, suffix, label in (
                (MeasurementRole.LAMP_BACKGROUND, "background", "lamp background"),
                (MeasurementRole.LAMP, "signal", "lamp signal"),
            ):
                records = self._records(role, family)
                items.append(
                    ChecklistItem(
                        f"lamp-{family}-{suffix}",
                        f"{family} {label}",
                        ChecklistState.DONE if records else ChecklistState.WAITING,
                        records[-1].path.name
                        if records
                        else "classify a measured SIF explicitly",
                    )
                )
        alignment_done = alignment.alignment_state is AlignmentState.ALIGNED
        items.append(
            ChecklistItem(
                "alignment",
                "Lamp alignment solved and reviewed",
                ChecklistState.DONE if alignment_done else ChecklistState.WAITING,
                f"{len(alignment.anchors)} accepted anchor(s)",
            )
        )
        items.append(
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
                else "generate after measurements and alignment are complete",
            )
        )
        items.append(
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
            )
        )
        return tuple(items)

    def _saved_snapshot_detail(self) -> str:
        if self.saved_snapshot is None:
            return self.last_error or "save only after every prerequisite is explicit"
        if self.wavelength_correction is None:
            return self.saved_snapshot.snapshot_id
        return f"{self.saved_snapshot.snapshot_id} — {self.wavelength_correction.reason}"

    def _measurement_pairs_ready(self) -> bool:
        if self._one(MeasurementRole.SPHERE) is None:
            return False
        if self._one(MeasurementRole.SPHERE_BACKGROUND) is None:
            return False
        return all(
            self._records(MeasurementRole.LAMP, family)
            and self._records(MeasurementRole.LAMP_BACKGROUND, family)
            for family in self.required_lamps
        )

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
    ) -> dict[str, str]:
        """Compose commented ordinary TOML from the measured session state."""

        if not self._composition_ready(alignment):
            raise SnapshotError(
                "sphere/lamp pairs, sphere-factor result, and aligned anchors are required"
            )
        sphere = self._one(MeasurementRole.SPHERE)
        sphere_background = self._one(MeasurementRole.SPHERE_BACKGROUND)
        assert sphere is not None and sphere_background is not None
        lamp_records = tuple(
            record
            for family in self.required_lamps
            for role in (MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND)
            for record in self._records(role, family)
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
            f"lamps = [{', '.join(_toml_string(item) for item in self.required_lamps)}]",
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
            'sphere = "sphere.sif"',
            'sphere_background = "sphere_bg.sif"',
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
    ) -> dict[str, Path]:
        """Atomically publish a new identity's generated TOML bundle."""

        destination_parent = Path(destination_root)
        destination_parent.mkdir(parents=True, exist_ok=True)
        destination = destination_parent / snapshot_id
        if destination.exists():
            self.toml_state = TomlState.FAILED
            self.last_error = f"configuration identity already exists: {destination}"
            raise SnapshotError(self.last_error)
        try:
            texts = self.compose_tomls(snapshot_id, alignment)
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
            os.replace(staging, destination)
            shutil.rmtree(staging_parent, ignore_errors=True)
        except Exception as exc:
            if "staging_parent" in locals():
                shutil.rmtree(staging_parent, ignore_errors=True)
            self.toml_state = TomlState.FAILED
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
        lamp = self._records(MeasurementRole.LAMP, self.required_lamps[0])
        lamp_background = self._records(
            MeasurementRole.LAMP_BACKGROUND, self.required_lamps[0]
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
            alignment_lamp=", ".join(self.required_lamps),
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
            if record.role in {MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND}
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
                lamps=self.required_lamps,
                lamp_files=lamp_files,
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
