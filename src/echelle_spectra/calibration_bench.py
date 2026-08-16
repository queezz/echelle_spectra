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
    detector_points_from_lines,
    fit_rigid_transform,
    fit_single_gaussian_centroid,
    measure_detector_window_saturation,
    measure_line_centroids,
    select_candidate_lines,
)

__all__ = [
    "AlignmentState",
    "Anchor",
    "AnchorFitResult",
    "AutoAnchorRejection",
    "AutoAnchorResult",
    "BenchFrame",
    "CalibrationBenchSession",
    "FileFingerprint",
    "FileLoadState",
    "FrameLoader",
    "Residual",
    "SaturationState",
    "StableFileResult",
    "StableFileState",
    "StableSifWatcher",
]


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

    def lines_for_order(self, order_idx: int | None = None) -> tuple[CalibrationTableLine, ...]:
        selected = self.selected_order if order_idx is None else int(order_idx)
        return tuple(line for line in self.reference_lines() if line.order_idx == selected)

    def fit_anchor_at(self, order_idx: int, clicked_pixel: float) -> AnchorFitResult:
        """Fit the nearest line of the assigned lamp around the clicked pixel."""

        if self.frame is None:
            return AnchorFitResult(False, "no SIF frame loaded")
        reference = self.reference
        if reference is not None and not reference.is_referenceable:
            return AnchorFitResult(False, reference.message)
        candidates = self.lines_for_order(order_idx)
        if not candidates:
            scope = (
                f"no {reference.catalog_label} rows"
                if reference is not None and reference.catalog_label
                else "no calibration rows"
            )
            return AnchorFitResult(False, f"{scope} for this order")
        line = min(candidates, key=lambda item: abs(item.center_pixel - clicked_pixel))
        if abs(line.center_pixel - clicked_pixel) > self.click_match_radius_px:
            known = (
                f"{reference.catalog_label} calibration row"
                if reference is not None and reference.catalog_label
                else "known calibration row"
            )
            return AnchorFitResult(False, f"click is not near a {known}")

        probe_line = CalibrationTableLine(
            line.order_idx,
            clicked_pixel - self.fit_window_radius_px,
            clicked_pixel + self.fit_window_radius_px,
            float(clicked_pixel),
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
            expected_center_px=float(clicked_pixel),
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

    def anchor_rows(self) -> tuple[Anchor, ...]:
        return tuple(sorted(self.anchors.values(), key=lambda item: item.key))
