"""Headless diffraction-order pattern extraction helpers.

This module lifts the order-pattern workflow out of the CMOS calibration
notebook.  It operates on already-loaded detector arrays so notebooks, tests,
or a future CLI can choose how to load SIF/TIFF/FITS data.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, Optional, Sequence

import numpy as np
from scipy.signal import savgol_filter

try:
    import peakutils
except ImportError as exc:  # pragma: no cover - dependency is declared by the package
    raise ImportError("pattern_extraction requires peakutils") from exc

__all__ = [
    "PatternExtractionConfig",
    "PatternColumnDetection",
    "PatternExtractionResult",
    "PatternExtractionTrial",
    "average_detector_frames",
    "subtract_background",
    "sample_columns",
    "amplification_curve",
    "detect_order_peaks_at_column",
    "detect_order_peaks",
    "fit_order_traces",
    "extract_order_pattern",
    "trial_order_pattern_extraction",
]


@dataclass(frozen=True)
class PatternExtractionConfig:
    """Configuration for semi-automated order-pattern extraction."""

    expected_orders: Optional[int] = 29
    sample_step_px: int = 150
    sample_count: int = 10
    smooth_window_px: int = 21
    smooth_polyorder: int = 1
    amplification_rate: float = 3e-3
    peak_threshold: float = 0.13
    peak_min_dist_px: int = 50
    baseline_poly_deg: int = 6
    trace_poly_degree: int = 2
    round_to_int: bool = True


@dataclass(frozen=True)
class PatternColumnDetection:
    """Detected order centers in one detector column."""

    column_px: int
    row_peaks_px: np.ndarray

    @property
    def n_peaks(self) -> int:
        """Number of detected peaks in the column."""
        return int(self.row_peaks_px.size)


@dataclass(frozen=True)
class PatternExtractionResult:
    """Result of fitting an order-pattern table."""

    pattern: np.ndarray
    columns_px: np.ndarray
    detections: List[PatternColumnDetection]
    coefficients: np.ndarray

    @property
    def n_orders(self) -> int:
        """Number of fitted diffraction orders."""
        return int(self.pattern.shape[1])


@dataclass(frozen=True)
class PatternExtractionTrial:
    """One attempted extraction setting and its diagnostics."""

    threshold: float
    columns_px: np.ndarray
    detections: List[PatternColumnDetection]
    result: Optional[PatternExtractionResult] = None
    error: str = ""

    @property
    def peak_counts(self) -> np.ndarray:
        """Detected peak count per sampled column."""
        return np.array([d.n_peaks for d in self.detections], dtype=int)

    @property
    def success(self) -> bool:
        """Return True when this trial produced a fitted pattern."""
        return self.result is not None

    @property
    def expected_count_matches(self) -> int:
        """Number of sampled columns matching the expected order count."""
        if self.result is not None:
            expected = self.result.n_orders
        else:
            expected = int(np.median(self.peak_counts)) if self.peak_counts.size else 0
        return int(np.count_nonzero(self.peak_counts == expected))

    @property
    def peak_count_spread(self) -> int:
        """Peak count range across sampled columns."""
        counts = self.peak_counts
        if counts.size == 0:
            return 0
        return int(counts.max() - counts.min())


def average_detector_frames(images: np.ndarray) -> np.ndarray:
    """Return a 2D mean detector image from a 2D image or 3D frame stack."""
    arr = np.asarray(images, dtype=float)
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3:
        return arr.mean(axis=0)
    raise ValueError(f"images must be 2D or 3D, got shape {arr.shape}")


def subtract_background(signal_images: np.ndarray, background_images: np.ndarray) -> np.ndarray:
    """Average signal/background frames and return the background-subtracted image."""
    signal = average_detector_frames(signal_images)
    background = average_detector_frames(background_images)
    if signal.shape != background.shape:
        raise ValueError(
            f"signal and background image shapes differ: {signal.shape} vs {background.shape}"
        )
    return signal - background


def sample_columns(
    ncols: int,
    step_size: int = 150,
    num_steps: int = 10,
    center_col: Optional[int] = None,
) -> np.ndarray:
    """Return evenly spaced detector columns centered near ``center_col``."""
    if ncols <= 0:
        raise ValueError("ncols must be positive")
    if step_size <= 0:
        raise ValueError("step_size must be positive")
    if num_steps <= 0:
        raise ValueError("num_steps must be positive")

    center = ncols // 2 if center_col is None else int(center_col)
    start = center - (num_steps // 2 * step_size)
    columns = np.arange(start, start + num_steps * step_size, step_size, dtype=int)
    if np.any(columns < 0) or np.any(columns >= ncols):
        raise ValueError(
            f"sample columns must be within [0, {ncols}); got {columns.tolist()}"
        )
    return columns


def amplification_curve(nrows: int, rate: float = 3e-3) -> np.ndarray:
    """Return the row-dependent amplification used for weak high-order tails."""
    if nrows <= 0:
        raise ValueError("nrows must be positive")
    rows = np.arange(nrows, dtype=float)
    return np.exp(float(rate) * rows)


def _validated_savgol_window(window_px: int, nrows: int) -> int:
    if nrows < 3:
        raise ValueError("image must have at least three rows")
    window = int(window_px)
    if window % 2 == 0:
        window += 1
    if window > nrows:
        window = nrows if nrows % 2 == 1 else nrows - 1
    if window < 3:
        window = 3
    return window


def detect_order_peaks_at_column(
    image: np.ndarray,
    column_px: int,
    config: PatternExtractionConfig = PatternExtractionConfig(),
    amplification: Optional[np.ndarray] = None,
) -> PatternColumnDetection:
    """Detect diffraction-order center rows in one detector column."""
    arr = np.asarray(image, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"image must be 2D, got shape {arr.shape}")

    nrows, ncols = arr.shape
    col = int(column_px)
    if not 0 <= col < ncols:
        raise ValueError(f"column_px must be within [0, {ncols}), got {col}")

    amp = amplification
    if amp is None:
        amp = amplification_curve(nrows, config.amplification_rate)
    amp = np.asarray(amp, dtype=float)
    if amp.shape != (nrows,):
        raise ValueError(f"amplification must have shape ({nrows},), got {amp.shape}")

    window = _validated_savgol_window(config.smooth_window_px, nrows)
    smoothed = savgol_filter(arr[:, col], window, config.smooth_polyorder) * amp
    baseline = peakutils.baseline(smoothed, config.baseline_poly_deg)
    peaks = peakutils.indexes(
        smoothed - baseline,
        thres=config.peak_threshold,
        min_dist=config.peak_min_dist_px,
    )
    return PatternColumnDetection(column_px=col, row_peaks_px=np.asarray(peaks, dtype=int))


def detect_order_peaks(
    image: np.ndarray,
    columns_px: Sequence[int],
    config: PatternExtractionConfig = PatternExtractionConfig(),
) -> List[PatternColumnDetection]:
    """Detect order peaks in each requested detector column."""
    arr = np.asarray(image, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"image must be 2D, got shape {arr.shape}")

    amp = amplification_curve(arr.shape[0], config.amplification_rate)
    return [
        detect_order_peaks_at_column(arr, int(column), config=config, amplification=amp)
        for column in columns_px
    ]


def fit_order_traces(
    detections: Sequence[PatternColumnDetection],
    ncols: int,
    config: PatternExtractionConfig = PatternExtractionConfig(),
) -> PatternExtractionResult:
    """Fit polynomial traces through detected order peaks and build a pattern table."""
    if not detections:
        raise ValueError("at least one column detection is required")
    if ncols <= 0:
        raise ValueError("ncols must be positive")

    peak_counts = np.array([d.n_peaks for d in detections], dtype=int)
    if np.any(peak_counts != peak_counts[0]):
        details = ", ".join(f"{d.column_px}: {d.n_peaks}" for d in detections)
        raise ValueError(f"all sampled columns must detect the same order count ({details})")

    n_orders = int(peak_counts[0])
    if config.expected_orders is not None and n_orders != int(config.expected_orders):
        raise ValueError(f"expected {config.expected_orders} orders, detected {n_orders}")
    if len(detections) <= config.trace_poly_degree:
        raise ValueError(
            "number of sampled columns must be greater than trace polynomial degree"
        )

    columns = np.array([d.column_px for d in detections], dtype=float)
    peaks = np.vstack([d.row_peaks_px for d in detections]).astype(float)
    all_columns = np.arange(ncols, dtype=float)

    coefficients = np.vstack(
        [
            np.polyfit(columns, peaks[:, order_idx], config.trace_poly_degree)
            for order_idx in range(n_orders)
        ]
    )
    pattern = np.column_stack([np.poly1d(coef)(all_columns) for coef in coefficients])
    if config.round_to_int:
        pattern = np.rint(pattern).astype(int)

    return PatternExtractionResult(
        pattern=pattern,
        columns_px=columns.astype(int),
        detections=list(detections),
        coefficients=coefficients,
    )


def extract_order_pattern(
    image: np.ndarray,
    config: PatternExtractionConfig = PatternExtractionConfig(),
    columns_px: Optional[Sequence[int]] = None,
) -> PatternExtractionResult:
    """Detect and fit a full order-pattern table from a background-subtracted image."""
    arr = np.asarray(image, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"image must be 2D, got shape {arr.shape}")

    if columns_px is None:
        columns = sample_columns(
            arr.shape[1],
            step_size=config.sample_step_px,
            num_steps=config.sample_count,
        )
    else:
        columns = np.asarray(columns_px, dtype=int)

    detections = detect_order_peaks(arr, columns, config=config)
    return fit_order_traces(detections, arr.shape[1], config=config)


def trial_order_pattern_extraction(
    image: np.ndarray,
    config: PatternExtractionConfig = PatternExtractionConfig(),
    threshold_values: Sequence[float] = (0.10, 0.11, 0.12, 0.13, 0.14),
    column_start_values: Optional[Sequence[int]] = None,
    column_step_px: Optional[int] = None,
    sample_count: Optional[int] = None,
) -> List[PatternExtractionTrial]:
    """Try several threshold/column settings and return ranked diagnostics.

    Successful fits are ranked first. Failed trials are still useful because
    their per-column peak counts identify where detection missed or added an
    order.
    """
    arr = np.asarray(image, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"image must be 2D, got shape {arr.shape}")

    step = config.sample_step_px if column_step_px is None else int(column_step_px)
    count = config.sample_count if sample_count is None else int(sample_count)
    if column_start_values is None:
        centered = sample_columns(arr.shape[1], step_size=step, num_steps=count)
        starts = [int(centered[0]), int(centered[0] + step), int(centered[0] + 2 * step)]
    else:
        starts = [int(start) for start in column_start_values]

    trials: List[PatternExtractionTrial] = []
    for start in starts:
        columns = np.arange(start, start + count * step, step, dtype=int)
        if np.any(columns < 0) or np.any(columns >= arr.shape[1]):
            continue

        for threshold in threshold_values:
            trial_config = replace(config, peak_threshold=float(threshold))
            detections = detect_order_peaks(arr, columns, config=trial_config)
            result: Optional[PatternExtractionResult] = None
            error = ""
            try:
                result = fit_order_traces(detections, arr.shape[1], config=trial_config)
            except ValueError as exc:
                error = str(exc)

            trials.append(
                PatternExtractionTrial(
                    threshold=float(threshold),
                    columns_px=columns,
                    detections=detections,
                    result=result,
                    error=error,
                )
            )

    return sorted(
        trials,
        key=lambda trial: (
            trial.success,
            trial.expected_count_matches,
            -trial.peak_count_spread,
            -abs(trial.threshold - config.peak_threshold),
            -int(trial.columns_px[0]) if trial.columns_px.size else 0,
        ),
        reverse=True,
    )
