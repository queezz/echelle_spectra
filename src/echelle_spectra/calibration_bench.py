"""UI-independent state and loading helpers for the live calibration bench.

The Qt bench is intentionally an adapter over this module.  Folder polling,
file-stability decisions, anchor mutation, centroid fitting, saturation
verdicts, rigid alignment, failure, and recovery can therefore be exercised
without an event loop.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from time import time_ns
from typing import TYPE_CHECKING, Callable

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - the campaign owns lamp naming, not this
    from .calibration_campaign import LampReferenceSet

from .tools.calibration_alignment import (
    CalibrationTableLine,
    DetectorWindowSaturation,
    LineCentroidFit,
    RigidTransform,
    TableVetting,
    apply_rigid_correction_to_lines,
    detector_points_from_lines,
    fit_rigid_transform,
    fit_single_gaussian_centroid,
    measure_detector_window_saturation,
    measure_line_centroids,
    select_candidate_lines,
)

__all__ = [
    "BAND_OFFSET_ATTENTION_ROWS",
    "AlignmentState",
    "Anchor",
    "AnchorFitResult",
    "AutoAnchorRejection",
    "AutoAnchorResult",
    "BandCenterOffsets",
    "BenchFrame",
    "CalibrationBenchSession",
    "FileFingerprint",
    "FileLoadState",
    "FrameLoader",
    "OrderBandOffset",
    "Residual",
    "SaturationState",
    "ScienceValidation",
    "ScienceValidationState",
    "StableFileResult",
    "StableFileState",
    "StableSifWatcher",
    "band_center_offsets",
]

#: How far the light may sit from the chosen pattern before the bench says so,
#: in detector rows.  This is well inside the extraction half width, so an
#: offset under it still sums each order over its own band; past it the trace
#: is riding the band's shoulder and the extracted flux starts to leave.
BAND_OFFSET_ATTENTION_ROWS = 1.5


class StableFileState(Enum):
    """Folder-watcher verdict for the newest SIF candidate."""

    EMPTY = "empty"
    CHANGING = "changing"
    STABLE = "stable"
    ALREADY_EMITTED = "already-emitted"
    FAILED = "failed"


class FileLoadState(Enum):
    """Lifecycle of a stable file after the watcher emits it."""

    WAITING = "waiting"
    LOADING = "loading"
    LOADED = "loaded"
    FAILED = "failed"


class AlignmentState(Enum):
    """Alignment lifecycle derived solely from frame and anchor state."""

    WAITING_FOR_FRAME = "waiting-for-frame"
    EMPTY = "empty"
    COLLECTING = "collecting"
    ALIGNED = "aligned"
    FAILED = "failed"


class SaturationState(Enum):
    """Raw-detector saturation verdict for a fitted anchor."""

    CLEAR = "clear"
    SATURATED = "saturated"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class FileFingerprint:
    """Identity used to decide whether a candidate file stopped changing."""

    size: int
    modified_ns: int


@dataclass(frozen=True)
class StableFileResult:
    """One deterministic folder poll."""

    state: StableFileState
    path: Path | None = None
    fingerprint: FileFingerprint | None = None
    unchanged_polls: int = 0
    reason: str = ""
    ready_path: Path | None = None


class StableSifWatcher:
    """Poll a folder and emit the newest SIF after repeated identical stats.

    Stability is intentionally defined by observations, not sleeps: the same
    newest path must have the same byte size and modification timestamp for
    ``required_unchanged_polls`` consecutive polls and meet ``minimum_age_s``.
    The latter defaults to zero so tests and callers can control the policy.
    """

    def __init__(
        self,
        folder: str | Path,
        *,
        required_unchanged_polls: int = 2,
        minimum_age_s: float = 0.0,
    ) -> None:
        if required_unchanged_polls < 1:
            raise ValueError("required_unchanged_polls must be at least 1")
        if minimum_age_s < 0:
            raise ValueError("minimum_age_s must not be negative")
        self.folder = Path(folder)
        self.required_unchanged_polls = int(required_unchanged_polls)
        self.minimum_age_ns = int(float(minimum_age_s) * 1_000_000_000)
        self._candidate: Path | None = None
        self._fingerprint: FileFingerprint | None = None
        self._unchanged_polls = 0
        self._emitted: tuple[Path, FileFingerprint] | None = None

    def poll(self, *, now_ns: int | None = None) -> StableFileResult:
        """Inspect the folder once and return a complete watcher verdict."""

        try:
            candidates = []
            for path in self.folder.iterdir():
                if path.suffix.lower() != ".sif":
                    continue
                try:
                    stat = path.stat()
                except OSError:
                    continue
                if path.is_file():
                    candidates.append(
                        (int(stat.st_mtime_ns), path.name.casefold(), path, stat)
                    )
        except OSError as exc:
            self._reset_candidate()
            return StableFileResult(StableFileState.FAILED, reason=str(exc))

        if not candidates:
            self._reset_candidate()
            return StableFileResult(StableFileState.EMPTY)

        _modified, _name, path, stat = max(candidates, key=lambda item: item[:2])
        fingerprint = FileFingerprint(int(stat.st_size), int(stat.st_mtime_ns))
        if path == self._candidate and fingerprint == self._fingerprint:
            self._unchanged_polls += 1
        else:
            self._candidate = path
            self._fingerprint = fingerprint
            self._unchanged_polls = 1

        age_ns = max(0, int(now_ns if now_ns is not None else time_ns()) - fingerprint.modified_ns)
        stable = (
            self._unchanged_polls >= self.required_unchanged_polls
            and age_ns >= self.minimum_age_ns
        )
        if not stable:
            return StableFileResult(
                StableFileState.CHANGING,
                path,
                fingerprint,
                self._unchanged_polls,
                "waiting for repeated unchanged observations",
            )

        identity = (path, fingerprint)
        if identity == self._emitted:
            return StableFileResult(
                StableFileState.ALREADY_EMITTED,
                path,
                fingerprint,
                self._unchanged_polls,
            )
        self._emitted = identity
        return StableFileResult(
            StableFileState.STABLE,
            path,
            fingerprint,
            self._unchanged_polls,
            ready_path=path,
        )

    def _reset_candidate(self) -> None:
        self._candidate = None
        self._fingerprint = None
        self._unchanged_polls = 0


@dataclass(frozen=True)
class BenchFrame:
    """Detector and extracted-order data consumed by the state machine.

    ``order_spectra`` is the mean over the acquisition's frames, which is what
    a 3-frame acquisition is shot for.  ``frame_order_spectra`` keeps each
    frame's own extraction beside it so a single frame can be fitted instead —
    empty for callers that only ever supply a mean.
    """

    path: Path
    images: np.ndarray
    detector_image: np.ndarray
    order_spectra: tuple[np.ndarray, ...]
    metadata: dict[str, object]
    #: Per-frame extractions, outer index frame, inner index order.
    frame_order_spectra: tuple[tuple[np.ndarray, ...], ...] = ()


class FrameLoader:
    """Load SIF detector/order data while reusing the established extractor."""

    def __init__(self, pattern: np.ndarray, *, half_width_px: int = 8) -> None:
        pattern_array = np.asarray(pattern, dtype=int)
        if pattern_array.ndim != 2 or not pattern_array.size:
            raise ValueError("pattern must have shape (detector columns, orders)")
        self.pattern = pattern_array
        self.half_width_px = int(half_width_px)

    def __call__(self, path: str | Path) -> BenchFrame:
        """Read one SIF and extract every order without absolute calibration."""

        from .tools.echelle import Calibrations, EchelleImage

        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"SIF file not found: {source}")
        calibration = Calibrations(dv=self.half_width_px)
        calibration.pattern = self.pattern
        image = EchelleImage(str(source), clbr=calibration)
        images = np.asarray(image.images, dtype=float)
        if images.ndim != 3:
            raise ValueError("SIF reader returned data outside (frames, rows, columns)")
        calibration.DIMO = images.shape[1]
        calibration.DIMW = images.shape[2]
        if self.pattern.shape[0] != calibration.DIMW:
            raise ValueError(
                "pattern width does not match SIF detector width: "
                f"{self.pattern.shape[0]} != {calibration.DIMW}"
            )
        calibration.make_cutting_masks()
        image.calculate_order_spectra()

        order_spectra = []
        per_frame: list[list[np.ndarray]] = [[] for _ in range(images.shape[0])]
        for order_idx in range(self.pattern.shape[1]):
            frames = [
                np.asarray(image.order_spectra[frame_idx][order_idx], dtype=float).reshape(-1)
                for frame_idx in range(images.shape[0])
            ]
            common_size = min(frame.size for frame in frames)
            for frame_idx, frame in enumerate(frames):
                per_frame[frame_idx].append(frame[:common_size])
            with np.errstate(invalid="ignore"):
                order_spectra.append(
                    np.nanmean(np.vstack([frame[:common_size] for frame in frames]), axis=0)
                )
        with np.errstate(invalid="ignore"):
            detector_image = np.nanmean(images, axis=0)
        return BenchFrame(
            path=source,
            images=images,
            detector_image=detector_image,
            order_spectra=tuple(order_spectra),
            metadata=dict(image.info),
            frame_order_spectra=tuple(tuple(rows) for rows in per_frame),
        )


@dataclass(frozen=True)
class OrderBandOffset:
    """Where one order's light sits relative to its pattern trace, in rows."""

    order_idx: int
    #: Measured centre minus pattern row.  ``None`` when this order could not
    #: be read at all, and ``reason`` then says why.
    offset_rows: float | None
    columns_measured: int = 0
    #: Columns dropped because their window was saturated or unreadable.
    columns_excluded: int = 0
    reason: str = ""


@dataclass(frozen=True)
class BandCenterOffsets:
    """Where a frame's order bands actually sit against a chosen pattern.

    The bench cannot see this any other way.  An anchor is a *dispersion*
    measurement — the centroid supplies the column and the detector row that
    goes with it is read out of the reference pattern, never off the frame — so
    a pattern belonging to another era passes every check the fit can make while
    every order is extracted off the edge of its own band.  That is not
    hypothetical: the 2024 folder calibrated against the packaged 2025 pattern
    reads about five rows out here, and its snapshot saved green.

    A positive offset means the light sits at a *higher* row index than the
    pattern says, which these detector plots draw above the trace; a negative
    one, below it.
    """

    orders: tuple[OrderBandOffset, ...] = ()
    #: Median over the orders that could be read.  ``None`` means none could.
    median_offset_rows: float | None = None
    columns_excluded: int = 0
    #: Why nothing could be measured, when nothing could.
    reason: str = ""

    @property
    def measured(self) -> bool:
        """Whether a real number is in, as opposed to a stated reason."""

        return self.median_offset_rows is not None

    @property
    def order_count(self) -> int:
        """How many orders carried enough light to be read."""

        return sum(1 for item in self.orders if item.offset_rows is not None)

    @property
    def extremes(self) -> tuple[float, float] | None:
        """The smallest and largest per-order offset, or ``None`` if unmeasured."""

        values = [item.offset_rows for item in self.orders if item.offset_rows is not None]
        if not values:
            return None
        return (min(values), max(values))

    def exceeds(self, threshold: float = BAND_OFFSET_ATTENTION_ROWS) -> bool:
        """Whether the light sits far enough off the pattern to be said aloud."""

        if self.median_offset_rows is None:
            return False
        return abs(self.median_offset_rows) > float(threshold)

    def summary(self) -> str:
        """The one sentence an operator reads, in rows and in plain words."""

        if self.median_offset_rows is None:
            return self.reason or "no order band could be read against the pattern"
        side = "below" if self.median_offset_rows < 0 else "above"
        return (
            f"order bands sit {abs(self.median_offset_rows):.1f} rows {side} the "
            f"chosen pattern (median of {self.order_count} order(s))"
        )


def band_center_offsets(
    image: np.ndarray,
    pattern: np.ndarray,
    *,
    dv: int = 8,
    saturation_level: float = 0.98 * 65535,
    column_step: int = 64,
    minimum_columns: int = 5,
    minimum_contrast: float = 1.0,
) -> BandCenterOffsets:
    """Measure a frame's own order-band centres against *pattern*, per order.

    For each sampled column the flux-weighted row centroid is taken over a
    window centred on the pattern's own trace, above the window's own floor.
    The window is deliberately wider than the extraction half width ``dv``: the
    bands of a sphere frame are broader than the strip that is summed out of
    them, so a window the width of ``dv`` sits entirely inside the flat top of
    the band and reports a centre it never saw the edges of.  It is still kept
    inside the gap to the neighbouring order, so no order is measured on its
    neighbour's light.

    A column whose window holds a *cluster* of full-scale pixels is left out,
    the same way exposure triage separates real saturation from a lone hot
    pixel: two adjacent full-scale pixels are saturation and the column goes,
    one on its own is an anomaly and the column stays.

    The reading is only meaningful on a frame that lights every order across
    the detector — the integrating sphere is exactly that frame, which is why
    it is the one the bench reads it on.
    """

    frame = np.asarray(image, dtype=float)
    if frame.ndim == 3:
        with np.errstate(invalid="ignore"):
            frame = np.nanmean(frame, axis=0)
    if frame.ndim != 2:
        raise ValueError("image must be a detector frame or a stack of them")
    rows = np.asarray(pattern, dtype=float)
    if rows.ndim != 2 or not rows.size:
        raise ValueError("pattern must have shape (detector columns, orders)")
    if rows.shape[0] != frame.shape[1]:
        return BandCenterOffsets(
            reason=(
                f"the chosen pattern is {rows.shape[0]} column(s) wide and this "
                f"frame is {frame.shape[1]} — they do not describe the same detector"
            )
        )
    if dv < 1:
        raise ValueError("dv must be at least one row")

    columns = np.arange(0, rows.shape[0], max(1, int(column_step)), dtype=int)
    radius = _band_search_radius(rows, columns, int(dv))
    measured = [
        _order_band_offset(
            frame,
            rows,
            order_idx,
            columns,
            radius=radius,
            saturation_level=float(saturation_level),
            minimum_columns=int(minimum_columns),
            minimum_contrast=float(minimum_contrast),
        )
        for order_idx in range(rows.shape[1])
    ]
    excluded_total = sum(item.columns_excluded for item in measured)
    offsets = [item.offset_rows for item in measured if item.offset_rows is not None]
    if not offsets:
        return BandCenterOffsets(
            tuple(measured),
            None,
            excluded_total,
            "no order carried enough light to locate its band on this frame",
        )
    return BandCenterOffsets(
        tuple(measured), float(np.median(offsets)), excluded_total
    )


def _order_band_offset(
    frame: np.ndarray,
    rows: np.ndarray,
    order_idx: int,
    columns: np.ndarray,
    *,
    radius: int,
    saturation_level: float,
    minimum_columns: int,
    minimum_contrast: float,
) -> OrderBandOffset:
    """One order's own offset: the median of its readable sampled columns."""

    samples: list[float] = []
    excluded = 0
    for column in columns:
        row = int(round(float(rows[column, order_idx])))
        low, high = row - radius, row + radius + 1
        if low < 0 or high > frame.shape[0]:
            excluded += 1
            continue
        window = frame[low:high, column]
        if not np.all(np.isfinite(window)):
            excluded += 1
            continue
        # Two adjacent full-scale pixels are saturation and this column goes;
        # one on its own is a cosmic ray or a hot pixel and it stays.  That is
        # the same separation exposure triage makes over the whole frame.
        full_scale = window >= saturation_level
        if np.any(full_scale[:-1] & full_scale[1:]):
            excluded += 1
            continue
        weights = window - float(np.min(window))
        # One digitizer count.  A window whose brightest pixel does not stand a
        # single count above its own floor holds no band, only arithmetic, and
        # weighting that arithmetic would put a centroid on the window's edge
        # and report it as a twenty-row shift.
        if float(np.max(weights)) < minimum_contrast:
            excluded += 1
            continue
        total = float(np.sum(weights))
        if total <= 0.0:
            excluded += 1
            continue
        index = np.arange(low, high, dtype=float)
        samples.append(float(np.sum(weights * index) / total) - float(row))
    if len(samples) < minimum_columns:
        return OrderBandOffset(
            order_idx,
            None,
            len(samples),
            excluded,
            f"only {len(samples)} of {columns.size} sampled column(s) could be "
            "read; this order carries no usable light here",
        )
    return OrderBandOffset(order_idx, float(np.median(samples)), len(samples), excluded)


def _band_search_radius(pattern: np.ndarray, columns: np.ndarray, dv: int) -> int:
    """A window wide enough to see both edges of a band, and no wider.

    Wide enough is two and a half extraction half widths, which on this
    instrument spans the whole band with room on each side.  No wider is half
    the closest approach of two neighbouring traces, so the window can never
    reach into the next order's light.
    """

    radius = max(int(dv), int(round(2.5 * dv)))
    if pattern.shape[1] < 2:
        return radius
    gaps = np.diff(pattern[columns, :], axis=1)
    closest = float(np.min(np.abs(gaps))) if gaps.size else 0.0
    if closest <= 0.0:
        return radius
    return int(max(dv, min(radius, closest / 2.0 - 1.0)))


@dataclass(frozen=True)
class Anchor:
    """One accepted known-line/centroid pair with raw-detector QC."""

    line: CalibrationTableLine
    fit: LineCentroidFit
    saturation: DetectorWindowSaturation

    @property
    def key(self) -> tuple[int, float, float]:
        return (self.line.order_idx, self.line.center_pixel, self.line.wavelength_nm)


@dataclass(frozen=True)
class Residual:
    """Rigid-fit residual for one accepted anchor."""

    key: tuple[int, float, float]
    order_idx: int
    wavelength_nm: float
    dx_px: float
    dy_px: float
    magnitude_px: float


@dataclass(frozen=True)
class AnchorFitResult:
    """Result of one click-to-fit attempt, accepted or safely rejected."""

    accepted: bool
    reason: str
    anchor: Anchor | None = None
    saturation_state: SaturationState = SaturationState.UNAVAILABLE


class ScienceValidationState(Enum):
    """Whether the solved alignment has been checked against science lines."""

    #: No transform yet, so there is nothing to validate.
    NO_ALIGNMENT = "no-alignment"
    #: Solved, but this bench holds no frame emitting Balmer or Fulcher light.
    NO_FRAME = "no-frame"
    #: A hydrogen frame is here, but no target could be fitted in it.
    NO_LINES = "no-lines"
    #: Science lines were measured against the solved wavelength solution.
    MEASURED = "measured"


@dataclass(frozen=True)
class ScienceValidation:
    """Agreement between the solved solution and the lines physics knows.

    Anchor RMS is self-consistency: it says the anchors agree with each other
    in pixels.  The BH paper's standard was agreement with Fulcher-alpha in
    nanometres, which is a different question and the one that matters.  This
    carries the answer, or states plainly why it cannot be answered here.
    """

    state: ScienceValidationState
    message: str
    results: tuple = ()
    rms_residual_nm: float | None = None
    median_residual_nm: float | None = None
    worst_line: str = ""

    @property
    def measured(self) -> bool:
        """Whether real numbers are in, as opposed to a stated reason."""

        return self.state is ScienceValidationState.MEASURED

    @property
    def line_count(self) -> int:
        return len(self.results)


@dataclass(frozen=True)
class AutoAnchorRejection:
    """One trusted row the auto-anchor pass measured and declined."""

    line: CalibrationTableLine
    reason: str


@dataclass(frozen=True)
class AutoAnchorResult:
    """What one auto-anchor pass considered, accepted, and declined.

    ``ran`` is false when the pass never got as far as measuring anything —
    no frame open, a lamp that cannot be referenced, or a table carrying no
    trusted rows for it — and ``reason`` then says which.
    """

    ran: bool
    reason: str = ""
    accepted: tuple[Anchor, ...] = ()
    rejected: tuple[AutoAnchorRejection, ...] = ()

    @property
    def considered(self) -> int:
        """How many trusted rows the pass actually measured."""

        return len(self.accepted) + len(self.rejected)


class CalibrationBenchSession:
    """Explicit live-bench state transitions with no Qt dependency."""

    def __init__(
        self,
        pattern: np.ndarray,
        lines: Sequence[CalibrationTableLine],
        *,
        saturation_level: float = 0.98 * 65535,
        minimum_snr: float = 5.0,
        fit_window_radius_px: int = 18,
        click_match_radius_px: float = 30.0,
        min_line_width_px: float = 4.0,
        max_line_width_px: float = 40.0,
        vetting: TableVetting | None = None,
    ) -> None:
        pattern_array = np.asarray(pattern, dtype=float)
        if pattern_array.ndim != 2 or not pattern_array.size:
            raise ValueError("pattern must have shape (detector columns, orders)")
        self.pattern = pattern_array
        self.lines = tuple(lines)
        #: The assigned lamp's own rows, once one is known.  ``None`` means no
        #: lamp has been named yet and every curated row is still clickable.
        self.reference: LampReferenceSet | None = None
        self.saturation_level = float(saturation_level)
        self.minimum_snr = float(minimum_snr)
        self.fit_window_radius_px = int(fit_window_radius_px)
        self.click_match_radius_px = float(click_match_radius_px)
        #: Interval widths the auto-anchor pass will consider, matching
        #: ``echelle-align``'s own gates so both paths trust the same rows.
        self.min_line_width_px = float(min_line_width_px)
        self.max_line_width_px = float(max_line_width_px)
        #: Whose vetting the loaded table's ``OK`` marks carry.  ``None`` means
        #: nobody asked the table, which is not the same as it having none.
        self.vetting: TableVetting | None = vetting
        self.file_state = FileLoadState.WAITING
        self.alignment_state = AlignmentState.WAITING_FOR_FRAME
        self.loading_path: Path | None = None
        self.frame: BenchFrame | None = None
        self.selected_order = 0
        #: Which frame of the acquisition is fitted.  ``None`` is the mean of
        #: every frame, which is what a 3-frame acquisition is shot for.
        self.selected_frame: int | None = None
        #: The assigned lamp background subtracted from the fitted spectrum,
        #: which is what ``echelle-align`` measures: signal minus background.
        self.background_path: Path | None = None
        self._background_spectra: tuple[np.ndarray, ...] = ()
        self.anchors: dict[tuple[int, float, float], Anchor] = {}
        self.transform: RigidTransform | None = None
        self.rms_px: float | None = None
        self.residuals: tuple[Residual, ...] = ()
        self.last_error = ""

    def begin_file_load(self, path: str | Path) -> None:
        """Enter the loading state without discarding the last good frame."""

        self.loading_path = Path(path)
        self.file_state = FileLoadState.LOADING
        self.last_error = ""

    def accept_frame(self, frame: BenchFrame) -> None:
        """Commit a successfully loaded frame and reset frame-bound anchors."""

        if frame.detector_image.ndim != 2:
            raise ValueError("detector_image must be two-dimensional")
        if frame.detector_image.shape[1] != self.pattern.shape[0]:
            raise ValueError("loaded frame and pattern have different detector widths")
        if len(frame.order_spectra) != self.pattern.shape[1]:
            raise ValueError("loaded frame and pattern have different order counts")
        self.frame = frame
        self.loading_path = None
        self.file_state = FileLoadState.LOADED
        self.selected_order = min(self.selected_order, len(frame.order_spectra) - 1)
        self.selected_frame = None
        # A new acquisition may belong to another lamp, so the background it is
        # paired with is re-established by the caller rather than carried over.
        self.background_path = None
        self._background_spectra = ()
        self.anchors.clear()
        self.transform = None
        self.rms_px = None
        self.residuals = ()
        self.alignment_state = AlignmentState.EMPTY
        self.last_error = ""

    def fail_file_load(self, path: str | Path, error: BaseException | str) -> None:
        """Record a load failure while preserving the previous usable frame."""

        self.loading_path = Path(path)
        self.file_state = FileLoadState.FAILED
        self.last_error = str(error)

    def load_file(self, path: str | Path, loader: Callable[[Path], BenchFrame]) -> bool:
        """Synchronous convenience transition used by tests and scripts."""

        source = Path(path)
        self.begin_file_load(source)
        try:
            frame = loader(source)
            self.accept_frame(frame)
        except Exception as exc:  # state boundary: surface loader failures uniformly
            self.fail_file_load(source, exc)
            return False
        return True

    def set_selected_order(self, order_idx: int) -> None:
        if order_idx < 0 or order_idx >= self.pattern.shape[1]:
            raise IndexError(f"order {order_idx} is outside the loaded pattern")
        self.selected_order = int(order_idx)

    @property
    def frame_count(self) -> int:
        """How many single frames the open acquisition offers, if any."""

        if self.frame is None:
            return 0
        if self.frame.frame_order_spectra:
            return len(self.frame.frame_order_spectra)
        return int(np.asarray(self.frame.images).shape[0])

    def set_selected_frame(self, frame_idx: int | None) -> None:
        """Fit the mean of the acquisition (``None``) or one single frame.

        Anchors are centroids measured on one particular spectrum, so switching
        which spectrum is on screen drops them rather than carrying numbers
        from the mean into a single-frame solution.
        """

        if frame_idx is not None:
            index = int(frame_idx)
            if index < 0 or index >= self.frame_count:
                raise IndexError(f"frame {frame_idx} is outside the open acquisition")
            frame_idx = index
        if frame_idx == self.selected_frame:
            return
        self.selected_frame = frame_idx
        if self.anchors:
            self.anchors.clear()
            self._recompute_alignment()

    def use_background_frame(self, frame: BenchFrame | None) -> None:
        """Subtract the assigned lamp's background from the fitted spectrum.

        This is what ``echelle-align`` measures — signal minus background — and
        it is the difference between fitting a line and fitting a line sitting
        on the lamp housing's own glow.  Anchors are dropped when the pairing
        changes, because a centroid belongs to the spectrum it was measured on.
        """

        if frame is None:
            changed = bool(self._background_spectra)
            self.background_path = None
            self._background_spectra = ()
        else:
            if len(frame.order_spectra) != self.pattern.shape[1]:
                raise ValueError(
                    "background frame and pattern have different order counts"
                )
            changed = Path(frame.path) != self.background_path
            self.background_path = Path(frame.path)
            self._background_spectra = tuple(
                np.asarray(row, dtype=float) for row in frame.order_spectra
            )
        if changed and self.anchors:
            self.anchors.clear()
            self._recompute_alignment()

    def background_label(self) -> str:
        """Say plainly whether a background is being subtracted, and which."""

        if self.background_path is None:
            return "no background subtracted"
        return f"minus {self.background_path.name}"

    def active_order_spectra(self) -> tuple[np.ndarray, ...]:
        """The extracted spectra a click is currently fitted against."""

        if self.frame is None:
            return ()
        if self.selected_frame is None or not self.frame.frame_order_spectra:
            spectra = tuple(self.frame.order_spectra)
        else:
            spectra = tuple(self.frame.frame_order_spectra[self.selected_frame])
        if not self._background_spectra:
            return spectra
        return tuple(
            self._subtract_background(signal, index)
            for index, signal in enumerate(spectra)
        )

    def _subtract_background(self, signal: np.ndarray, order_idx: int) -> np.ndarray:
        """Subtract one order's background over the samples the two share."""

        if order_idx >= len(self._background_spectra):
            return signal
        background = self._background_spectra[order_idx]
        shared = min(signal.size, background.size)
        if not shared:
            return signal
        corrected = np.array(signal, dtype=float)
        corrected[:shared] = corrected[:shared] - background[:shared]
        return corrected

    def active_images(self) -> np.ndarray:
        """The raw detector data the per-anchor saturation check reads.

        Single-frame fitting is judged on that frame's own raw pixels, so a
        line saturated only in one frame of the acquisition is rejected there
        and still fits on the frames that kept it in range.
        """

        if self.frame is None:
            raise ValueError("no frame is open")
        images = np.asarray(self.frame.images)
        if self.selected_frame is None or images.ndim != 3:
            return images
        return images[self.selected_frame : self.selected_frame + 1]

    def frame_choice_label(self) -> str:
        """Say plainly which spectrum the fit is measuring."""

        count = self.frame_count
        if self.frame is None:
            return "no acquisition open"
        if self.selected_frame is None:
            return f"mean of {count} frame(s)" if count > 1 else "single frame"
        return f"frame {self.selected_frame + 1} of {count}"

    def use_lamp_reference(self, reference: LampReferenceSet | None) -> None:
        """Scope click-to-fit to the rows one assigned lamp's own catalog holds.

        Anchors already accepted against another lamp's rows are dropped: they
        measured this frame against lines it never emitted, and keeping them
        would carry that error straight into the solved transform.
        """

        self.reference = reference
        if reference is None:
            return
        allowed = {
            (line.order_idx, line.center_pixel, line.wavelength_nm)
            for line in reference.lines
        }
        stale = [key for key in self.anchors if key not in allowed]
        for key in stale:
            del self.anchors[key]
        if stale:
            self._recompute_alignment()

    def reference_lines(self) -> tuple[CalibrationTableLine, ...]:
        """Return the rows a click may snap to: the assigned lamp's own."""

        if self.reference is None:
            return self.lines
        return tuple(self.reference.lines)

    def best_reference_order(self) -> int | None:
        """The order carrying most of the assigned lamp's own rows.

        A neon lamp's curated lines live in a handful of orders out of
        twenty-nine.  Opening the fit on order 0 and leaving the operator to
        hunt for them is the empty room again, one tab further along.
        """

        if self.reference is None:
            return None
        return self.reference.best_order

    def corrected_rows(
        self, rows: Sequence[CalibrationTableLine]
    ) -> tuple[CalibrationTableLine, ...]:
        """Place curated rows where *this frame* shows them.

        Until two anchors have solved a transform these are the base table's
        own pixels, which is how the first anchors are found at all: nothing
        yet knows the detector moved, so the operator's eye is the only
        correction there is.  Once a transform exists every row follows it,
        through the same ``apply_rigid_correction_to_lines`` that
        ``echelle-align --save`` writes its adjusted table with, so the sticks,
        the expected-lines panel and the click-to-fit window point at the pixel
        the line is actually on rather than at the one the 2024 table was
        written for.  Re-solving moves them again on the next refresh: there is
        nothing cached here to go stale.
        """

        rows = tuple(rows)
        if self.transform is None or not rows:
            return rows
        return tuple(
            apply_rigid_correction_to_lines(list(rows), self.pattern, self.transform)
        )

    def _order_rows(
        self, order_idx: int | None = None
    ) -> tuple[tuple[CalibrationTableLine, CalibrationTableLine], ...]:
        """One order's curated rows paired with where this frame shows them.

        The base row stays the anchor's identity and the transform's own input
        — solving from corrected positions would fold the correction into
        itself and shrink it to nothing on every pass — while the corrected row
        is what the operator sees and clicks.
        """

        selected = self.selected_order if order_idx is None else int(order_idx)
        base = tuple(
            line for line in self.reference_lines() if line.order_idx == selected
        )
        return tuple(zip(base, self.corrected_rows(base)))

    def lines_for_order(self, order_idx: int | None = None) -> tuple[CalibrationTableLine, ...]:
        """The order's clickable rows, placed where this frame shows them."""

        return tuple(shown for _, shown in self._order_rows(order_idx))

    def fit_anchor_at(self, order_idx: int, clicked_pixel: float) -> AnchorFitResult:
        """Fit the nearest line of the assigned lamp around the clicked pixel."""

        if self.frame is None:
            return AnchorFitResult(False, "no SIF frame loaded")
        reference = self.reference
        if reference is not None and not reference.is_referenceable:
            return AnchorFitResult(False, reference.message)
        candidates = self._order_rows(order_idx)
        if not candidates:
            scope = (
                f"no {reference.catalog_label} rows"
                if reference is not None and reference.catalog_label
                else "no calibration rows"
            )
            return AnchorFitResult(False, f"{scope} for this order")
        line, shown = min(
            candidates, key=lambda pair: abs(pair[1].center_pixel - clicked_pixel)
        )
        if abs(shown.center_pixel - clicked_pixel) > self.click_match_radius_px:
            known = (
                f"{reference.catalog_label} calibration row"
                if reference is not None and reference.catalog_label
                else "known calibration row"
            )
            return AnchorFitResult(False, f"click is not near a {known}")

        # Where to centre the +/-18 px probe.  With no transform solved the
        # click is the only estimate there is, and centring on it is what lets
        # the first anchors be placed on a frame the table is 18 px wrong
        # about.  Once a transform exists the corrected row is the better
        # estimate of the two, so a click landing anywhere within the match
        # radius still opens its window on the line rather than on the line's
        # edge — which is the whole 2019 defect.
        probe_center = float(
            shown.center_pixel if self.transform is not None else clicked_pixel
        )
        probe_line = CalibrationTableLine(
            line.order_idx,
            probe_center - self.fit_window_radius_px,
            probe_center + self.fit_window_radius_px,
            probe_center,
            line.wavelength_nm,
            line.species,
            line.comment,
        )
        detector_qc = measure_detector_window_saturation(
            self.active_images(),
            self.pattern,
            [probe_line],
            x_radius_px=self.fit_window_radius_px,
            saturation_level=self.saturation_level,
        )[0]
        detector_qc = DetectorWindowSaturation(
            line,
            detector_qc.peak_value,
            detector_qc.finite_pixels,
            detector_qc.saturated_pixels,
            detector_qc.saturated_fraction,
            detector_qc.reason,
            detector_qc.saturation_level,
        )
        if detector_qc.is_saturated:
            return AnchorFitResult(
                False,
                f"{line.species} {line.wavelength_nm:.3f} nm is saturated on the raw "
                f"detector ({detector_qc.saturated_pixels} full-scale pixel(s) in its "
                "window) — fit an unsaturated line instead, or the same line on the "
                "bright frame of the pair",
                saturation_state=SaturationState.SATURATED,
            )

        success, center, sigma, amplitude, baseline, snr, reason = fit_single_gaussian_centroid(
            self.active_order_spectra()[order_idx],
            expected_center_px=probe_center,
            window_radius_px=self.fit_window_radius_px,
            min_snr=self.minimum_snr,
        )
        fit = LineCentroidFit(
            line,
            center,
            sigma,
            amplitude,
            baseline,
            snr,
            success,
            reason,
        )
        if not fit.success:
            return AnchorFitResult(
                False,
                fit.reason or "centroid fit failed",
                saturation_state=SaturationState.CLEAR,
            )
        anchor = Anchor(line, fit, detector_qc)
        self.upsert_anchor(anchor)
        return AnchorFitResult(True, "anchor accepted", anchor, SaturationState.CLEAR)

    def trusted_rows(self) -> tuple[CalibrationTableLine, ...]:
        """The vetted rows an auto pass would measure, before any fitting.

        Same gates ``echelle-align`` selects on: the assigned lamp's own rows,
        carrying the curated table's ``OK`` mark, over an interval wide enough
        to be a line and narrow enough to be only one.
        """

        return tuple(
            select_candidate_lines(
                self.reference_lines(),
                species=None,
                require_ok=True,
                min_width_px=self.min_line_width_px,
                max_width_px=self.max_line_width_px,
            )
        )

    def auto_anchor(self) -> AutoAnchorResult:
        """Fit every trusted row of the assigned lamp, across all orders.

        This is the click-free half of :meth:`fit_anchor_at`: the same
        centroid fit, the same raw-detector saturation guard, and the same
        rigid solve, aimed at the rows ``echelle-align`` has always trusted
        — the curated table's own ``OK`` marks — instead of at one pixel the
        operator pointed to.  Anchors already collected survive unless this
        pass measures the same row again, so an auto pass can follow manual
        work without discarding it.  Every declined row comes back with its
        reason: an auto-pass that fails silently is a worse instrument than
        a manual one.
        """

        if self.frame is None:
            return AutoAnchorResult(False, "no SIF frame loaded")
        reference = self.reference
        if reference is not None and not reference.is_referenceable:
            return AutoAnchorResult(False, reference.message)

        spectra = self.active_order_spectra()
        candidates = [
            line for line in self.trusted_rows() if 0 <= line.order_idx < len(spectra)
        ]
        if not candidates:
            return AutoAnchorResult(False, self._no_trusted_rows_message())

        saturation = measure_detector_window_saturation(
            self.active_images(),
            self.pattern,
            candidates,
            x_radius_px=self.fit_window_radius_px,
            saturation_level=self.saturation_level,
        )
        fits = measure_line_centroids(
            spectra,
            candidates,
            window_radius_px=self.fit_window_radius_px,
            saturation_level=None,
            min_snr=self.minimum_snr,
        )

        accepted: list[Anchor] = []
        rejected: list[AutoAnchorRejection] = []
        for line, fit, detector_qc in zip(candidates, fits, saturation):
            if detector_qc.is_saturated:
                rejected.append(
                    AutoAnchorRejection(
                        line,
                        "saturated on the raw detector "
                        f"({detector_qc.saturated_pixels} full-scale pixel(s))",
                    )
                )
                continue
            if not fit.success:
                rejected.append(
                    AutoAnchorRejection(line, fit.reason or "centroid fit failed")
                )
                continue
            anchor = Anchor(line, fit, detector_qc)
            self._accept_anchor(anchor)
            accepted.append(anchor)

        # One solve for the whole pass: recomputing per anchor would fit the
        # transform dozens of times to reach the same answer.
        self._recompute_alignment()
        return AutoAnchorResult(True, "", tuple(accepted), tuple(rejected))

    def _no_trusted_rows_message(self) -> str:
        """Say why an auto pass found nothing it was willing to measure."""

        reference = self.reference
        rows = self.reference_lines()
        catalog = (
            f"{reference.catalog_label} rows"
            if reference is not None and reference.catalog_label
            else "curated rows"
        )
        return (
            f"none of the {len(rows)} {catalog} in this table carry the OK marks "
            "the auto-anchor trusts — click the lines you trust instead"
        )

    def _accept_anchor(self, anchor: Anchor) -> None:
        """Store one validated anchor without resolving the alignment."""

        if not anchor.fit.success:
            raise ValueError("cannot add an unsuccessful centroid fit")
        if anchor.saturation.is_saturated:
            raise ValueError("cannot add a saturated anchor")
        self.anchors[anchor.key] = anchor

    def upsert_anchor(self, anchor: Anchor) -> None:
        """Add or replace an anchor, then deterministically recompute alignment."""

        self._accept_anchor(anchor)
        self._recompute_alignment()

    def remove_anchor(self, key: tuple[int, float, float]) -> bool:
        """Remove an anchor and recompute; return whether it existed."""

        existed = self.anchors.pop(key, None) is not None
        if existed:
            self._recompute_alignment()
        return existed

    def clear_anchors(self) -> None:
        self.anchors.clear()
        self._recompute_alignment()

    def _recompute_alignment(self) -> None:
        self.transform = None
        self.rms_px = None
        self.residuals = ()
        if self.frame is None:
            self.alignment_state = AlignmentState.WAITING_FOR_FRAME
            return
        ordered = sorted(self.anchors.values(), key=lambda item: item.key)
        if not ordered:
            self.alignment_state = AlignmentState.EMPTY
            self.last_error = ""
            return
        if len(ordered) < 2:
            self.alignment_state = AlignmentState.COLLECTING
            self.last_error = ""
            return
        lines = [anchor.line for anchor in ordered]
        centers = [anchor.fit.center_pixel for anchor in ordered]
        try:
            expected = detector_points_from_lines(lines, self.pattern)
            measured = detector_points_from_lines(lines, self.pattern, centers)
            transform, rms_px = fit_rigid_transform(expected, measured)
            predicted = transform.apply(expected)
            delta = predicted - measured
        except (IndexError, ValueError, np.linalg.LinAlgError) as exc:
            self.alignment_state = AlignmentState.FAILED
            self.last_error = str(exc)
            return
        residuals = []
        for anchor, (dx_px, dy_px) in zip(ordered, delta):
            residuals.append(
                Residual(
                    anchor.key,
                    anchor.line.order_idx,
                    anchor.line.wavelength_nm,
                    float(dx_px),
                    float(dy_px),
                    float(np.hypot(dx_px, dy_px)),
                )
            )
        self.transform = transform
        self.rms_px = rms_px
        self.residuals = tuple(residuals)
        self.alignment_state = AlignmentState.ALIGNED
        self.last_error = ""

    def order_wavelengths(self) -> dict[int, np.ndarray]:
        """The solved wavelength of every detector column, per order.

        Built the way ``Calibrations`` builds it — a per-order polynomial
        through the table's own (pixel, wavelength) pairs, straight for an
        order carrying only a couple of rows and quadratic once it carries
        enough to justify the curve — but through the rows *as this frame
        shows them*, with the solved rigid correction already applied.  That
        is what makes the result a wavelength axis for this detector today
        rather than for the detector the table was written on.
        """

        if self.frame is None or self.transform is None:
            return {}
        corrected = apply_rigid_correction_to_lines(
            list(self.lines), self.pattern, self.transform
        )
        columns = self.pattern.shape[0]
        x = np.arange(columns, dtype=float)
        by_order: dict[int, list[CalibrationTableLine]] = {}
        for line in corrected:
            by_order.setdefault(line.order_idx, []).append(line)
        axes: dict[int, np.ndarray] = {}
        for order_idx, rows in by_order.items():
            if order_idx < 0 or order_idx >= self.pattern.shape[1] or len(rows) < 2:
                continue
            pixels = np.asarray([row.center_pixel for row in rows], dtype=float)
            waves = np.asarray([row.wavelength_nm for row in rows], dtype=float)
            degree = 1 if len(rows) < 3 else 2
            try:
                fit = np.polyfit(pixels, waves, degree)
            except (np.linalg.LinAlgError, ValueError):
                continue
            axes[order_idx] = np.poly1d(fit)(x)
        return axes

    def validate_science_lines(
        self, targets: Sequence[object] | None = None
    ) -> ScienceValidation:
        """Check the solved solution against the lines physics already knows.

        The bench's own lamps cannot answer this: neon, thorium and mercury do
        not emit Balmer or Fulcher light, so on a normal calibration evening
        the honest result is that there is no frame to validate against, and
        saying so is the point rather than a shortfall.  When a hydrogen frame
        *is* on the bench, its Balmer and Fulcher targets are fitted against
        the solved wavelength axis and the residuals are reported in
        nanometres, which is the standard the BH paper actually held.
        """

        from .tools.line_validation import (
            balmer_air_targets,
            bundled_fulcher_h2_q_branch_targets,
            summarize_validation,
            validate_lines,
        )

        if self.transform is None:
            return ScienceValidation(
                ScienceValidationState.NO_ALIGNMENT,
                "no transform is solved yet, so there is nothing to validate",
            )
        if not self._carries_science_lines():
            return ScienceValidation(
                ScienceValidationState.NO_FRAME,
                self._no_science_frame_message(),
            )
        axes = self.order_wavelengths()
        spectra = self.active_order_spectra()
        if targets is None:
            targets = [*balmer_air_targets(), *bundled_fulcher_h2_q_branch_targets()]

        results: list[object] = []
        for order_idx, axis in sorted(axes.items()):
            if order_idx >= len(spectra):
                continue
            spectrum = np.asarray(spectra[order_idx], dtype=float)
            shared = min(axis.size, spectrum.size)
            if shared < 8:
                continue
            low, high = float(np.nanmin(axis[:shared])), float(np.nanmax(axis[:shared]))
            in_range = [
                target
                for target in targets
                if low <= float(target.wavelength_nm) <= high
            ]
            if not in_range:
                continue
            results.extend(
                validate_lines(axis[:shared], spectrum[:shared], in_range)
            )
        if not results:
            return ScienceValidation(
                ScienceValidationState.NO_LINES,
                "a hydrogen frame is on the bench, but no Balmer or Fulcher line "
                "could be fitted in it — check the exposure before trusting the "
                "alignment",
            )
        summary = summarize_validation(results)
        worst = max(results, key=lambda item: abs(float(item.residual_nm)))
        return ScienceValidation(
            ScienceValidationState.MEASURED,
            f"{len(results)} science line(s) agree to "
            f"{summary.rms_residual_nm:.4f} nm RMS "
            f"(median {abs(summary.median_residual_nm):.4f} nm)",
            tuple(results),
            float(summary.rms_residual_nm),
            float(summary.median_residual_nm),
            str(worst.line),
        )

    def _carries_science_lines(self) -> bool:
        """Whether the assigned lamp emits the lines validation is made of."""

        reference = self.reference
        return reference is not None and reference.catalog_family == "fulcher"

    def _no_science_frame_message(self) -> str:
        reference = self.reference
        lamp = reference.lamp if reference is not None and reference.lamp else "this lamp"
        return (
            f"{lamp} emits no Balmer or Fulcher light, so the alignment cannot be "
            "checked against science lines here — validate against Fulcher when "
            "the first plasma data exists"
        )

    def anchor_rows(self) -> tuple[Anchor, ...]:
        return tuple(sorted(self.anchors.values(), key=lambda item: item.key))
