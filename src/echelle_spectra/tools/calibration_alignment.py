"""Headless wavelength-calibration alignment helpers.

This module supports the narrow v1 calibration workflow:

1. Read the curated wavelength lookup table already used by ``Calibrations``.
2. Measure single-line centroids from extracted order spectra.
3. Fit a rigid detector correction (translation + rotation, no stretch).
4. Persist the correction and optionally write an adjusted lookup table.

The GUI is deliberately not involved here.  These functions are small enough
to use from notebooks, tests, or a future command-line entry point.
"""

from __future__ import annotations

import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import curve_fit

try:  # pragma: no cover - exercised by Python version
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib

__all__ = [
    "CalibrationTableLine",
    "LineWindowStats",
    "DetectorWindowSaturation",
    "LineCentroidFit",
    "RigidTransform",
    "AlignmentSettings",
    "AlignmentRunConfig",
    "AlignmentRunResult",
    "TableVetting",
    "BH_PAPER_WAVELENGTH_TABLE",
    "load_wavelength_table",
    "read_table_metadata",
    "table_vetting",
    "select_candidate_lines",
    "measure_line_window_stats",
    "measure_detector_window_saturation",
    "rank_candidate_lines",
    "fit_single_gaussian_centroid",
    "measure_line_centroids",
    "fit_rigid_transform",
    "detector_points_from_lines",
    "apply_rigid_correction_to_lines",
    "apply_rigid_correction_to_pattern",
    "shift_lines_in_pixels",
    "write_wavelength_table",
    "PatternCorrection",
    "write_pattern_table",
    "write_corrected_pattern_table",
    "save_alignment_settings",
    "load_alignment_settings",
    "run_calibration_alignment",
]


@dataclass(frozen=True)
class CalibrationTableLine:
    """One active row from an echelle wavelength lookup table."""

    order_idx: int
    pixel_from: float
    pixel_to: float
    center_pixel: float
    wavelength_nm: float
    species: str
    comment: str = ""

    @property
    def width_px(self) -> float:
        """Feature width from the manually selected pixel interval."""
        return abs(self.pixel_to - self.pixel_from)


@dataclass(frozen=True)
class LineWindowStats:
    """Pre-fit diagnostics for a candidate calibration-line window."""

    line: CalibrationTableLine
    peak_pixel: float
    peak_value: float
    baseline: float
    noise: float
    prominence: float
    snr: float
    finite_pixels: int
    saturated_pixels: int
    saturated_fraction: float
    fit_candidate: bool
    reason: str = ""
    saturation_level: Optional[float] = None

    @property
    def is_saturated(self) -> bool:
        """Return True when at least one finite pixel reached the saturation threshold."""
        return self.saturated_pixels > 0


@dataclass(frozen=True)
class DetectorWindowSaturation:
    """Raw-detector saturation diagnostics near one expected line position."""

    line: CalibrationTableLine
    peak_value: float
    finite_pixels: int
    saturated_pixels: int
    saturated_fraction: float
    reason: str = ""
    saturation_level: Optional[float] = None

    @property
    def is_saturated(self) -> bool:
        """Return True when at least one detector pixel reached the saturation threshold."""
        return self.saturated_pixels > 0


@dataclass(frozen=True)
class LineCentroidFit:
    """Measured centroid and quality diagnostics for one calibration line."""

    line: CalibrationTableLine
    center_pixel: float
    sigma_px: float
    amplitude: float
    baseline: float
    snr: float
    success: bool
    reason: str = ""
    diagnostics: Optional[LineWindowStats] = None


@dataclass(frozen=True)
class RigidTransform:
    """Rigid 2D transform mapping expected detector points to measured points."""

    dx_px: float
    dy_px: float
    theta_rad: float

    @property
    def theta_deg(self) -> float:
        """Rotation angle in degrees."""
        return float(np.degrees(self.theta_rad))

    @property
    def matrix(self) -> np.ndarray:
        """2x2 rotation matrix."""
        c = float(np.cos(self.theta_rad))
        s = float(np.sin(self.theta_rad))
        return np.array([[c, -s], [s, c]], dtype=float)

    def apply(self, points_xy: np.ndarray) -> np.ndarray:
        """Apply the rigid transform to an ``(N, 2)`` point array."""
        pts = np.asarray(points_xy, dtype=float)
        return pts @ self.matrix.T + np.array([self.dx_px, self.dy_px], dtype=float)


@dataclass(frozen=True)
class AlignmentSettings:
    """Persisted calibration-alignment summary."""

    instrument_id: str
    base_wavelength_file: str
    transform: RigidTransform
    n_lines: int
    rms_px: float
    created_at: str = ""
    alignment_dataset_id: str = ""
    alignment_source_dir: str = ""
    alignment_lamp: str = ""
    signal_file: str = ""
    background_file: str = ""
    base_pattern_file: str = ""
    sphere_file: str = ""
    sphere_background_file: str = ""
    output_wavelength_file: str = ""
    output_pattern_file: str = ""
    notes: str = ""


@dataclass(frozen=True)
class AlignmentRunConfig:
    """Inputs and thresholds for one headless calibration-alignment run."""

    calibration_dir: Path
    signal_file: Path
    background_file: Path
    sphere_file: Path
    sphere_background_file: Path
    base_wavelength_file: str = "Th_wavelength_CMOS_20240305.txt"
    base_pattern_file: str = "pattern_CMOS_20250926.txt"
    integral_file: str = "integrating_sphere.txt"
    instrument_id: str = "lhd_cmos"
    alignment_dataset_id: str = "20250926"
    alignment_source_dir: str = "local/20250926_calib"
    alignment_lamp: str = "Ne"
    created_at: str = "2026-06-04"
    output_wavelength_file: str = "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    output_pattern_file: str = "pattern_CMOS_20250926_aligned_to_20250926.txt"
    window_radius_px: int = 18
    min_snr: float = 5.0
    saturation_level: float = 0.98 * 65535
    x_radius_px: int = 18
    y_radius_px: int = 4
    species: Tuple[str, ...] = ("NeI",)
    require_ok: bool = True
    min_width_px: float = 4.0
    max_width_px: float = 40.0
    notes: str = "Rigid detector correction from Neon lamp using accepted CMOS order pattern."


@dataclass(frozen=True)
class AlignmentRunResult:
    """Summary and outputs from one headless calibration-alignment run."""

    settings: AlignmentSettings
    adjusted_rows: List[CalibrationTableLine]
    rows: List[CalibrationTableLine]
    candidates: List[CalibrationTableLine]
    fits: List[LineCentroidFit]
    good_fits: List[LineCentroidFit]
    residual_xy: np.ndarray
    residual_px: np.ndarray
    expected_xy: np.ndarray
    measured_xy: np.ndarray
    predicted_xy: np.ndarray
    detector_saturation: List[DetectorWindowSaturation]
    ranked_stats: List[LineWindowStats]

    @property
    def centroid_dx_px(self) -> np.ndarray:
        """Return measured-minus-expected centroid shifts for good fits."""
        return np.asarray(
            [fit.center_pixel - fit.line.center_pixel for fit in self.good_fits],
            dtype=float,
        )


def load_wavelength_table(path: str | Path) -> List[CalibrationTableLine]:
    """Load active calibration rows from an echelle wavelength table.

    Inline comments after ``#`` are preserved in the returned ``comment``
    field so manual ``OK`` / ``?`` annotations remain available for filtering.
    """
    rows: List[CalibrationTableLine] = []
    with Path(path).open() as fh:
        for raw in fh:
            data_part, sep, comment_part = raw.partition("#")
            stripped = data_part.strip()
            if not stripped:
                continue

            tokens = stripped.split()
            if len(tokens) < 6:
                continue

            rows.append(
                CalibrationTableLine(
                    order_idx=int(tokens[0]),
                    pixel_from=float(tokens[1]),
                    pixel_to=float(tokens[2]),
                    center_pixel=float(tokens[3]),
                    wavelength_nm=float(tokens[4]),
                    species=tokens[5],
                    comment=comment_part.strip() if sep else "",
                )
            )
    return rows


#: The curated table whose ``OK`` marks carry the BH paper's own vetting.
#: That alignment was tried and tested for the paper's Balmer and Fulcher
#: analysis, and it is the pedigree that makes "trusted" mean something: an
#: ``OK`` beside a row is worth what the work behind it was worth, never
#: merely what typing those two letters costs.
BH_PAPER_WAVELENGTH_TABLE = "Th_wavelength_CMOS_20240305.txt"

#: Vetted lineages this package can vouch for, by the table that established
#: them.  A table outside this map is not accused of anything — it simply
#: carries no vetting anybody here recorded, and says so rather than
#: borrowing another table's authority.
_VETTED_SETS: dict[str, Tuple[str, str]] = {
    BH_PAPER_WAVELENGTH_TABLE: (
        "BH paper",
        "vetted during the BH-paper calibration and tried against that "
        "paper's own Balmer and Fulcher analysis",
    ),
}

#: How far a lineage is followed before the chain is called a cycle.
_MAX_LINEAGE_DEPTH = 8


@dataclass(frozen=True)
class TableVetting:
    """What one wavelength table's ``OK`` marks are worth, and on whose word.

    An aligned table records the table it was derived from in its own header,
    so a lineage can be walked from the file itself rather than assumed from
    its name.  ``vetted_set`` is empty when nothing in that lineage is a set
    this package can vouch for.
    """

    table: str
    vetted_set: str
    lineage: Tuple[str, ...]
    description: str
    vetted_table: str = ""

    @property
    def is_vetted(self) -> bool:
        """Whether these ``OK`` marks carry a vetting anybody recorded."""

        return bool(self.vetted_set)

    @property
    def message(self) -> str:
        """One sentence an operator can act on, either way."""

        if not self.is_vetted:
            return (
                f"the OK marks in {self.table} carry no recorded vetting — they "
                "are this table's own, not the BH paper's"
            )
        through = (
            f", inherited from {self.vetted_table}"
            if self.vetted_table and self.vetted_table != self.table
            else ""
        )
        return (
            f"the OK marks in {self.table} carry the {self.vetted_set} vetting"
            f"{through}: {self.description}"
        )


def read_table_metadata(path: str | Path) -> dict[str, str]:
    """Return the ``# Key: value`` header a wavelength table opens with.

    Only the leading comment block is read.  Curated tables also carry
    commented-out data rows further down, and those are rows rather than
    metadata however much they look alike.
    """

    metadata: dict[str, str] = {}
    with Path(path).open() as fh:
        for raw in fh:
            if not raw.lstrip().startswith("#"):
                break
            key, sep, value = raw.lstrip().lstrip("#").strip().partition(":")
            if sep:
                metadata.setdefault(key.strip().casefold(), value.strip())
    return metadata


def table_vetting(
    path: str | Path, search_dirs: Optional[Sequence[str | Path]] = None
) -> TableVetting:
    """Report which vetted set, if any, a wavelength table's OK marks carry.

    The lineage is read from the tables themselves: an adjusted table names
    its base in its header, so ``..._aligned_to_20250926.txt`` leads back to
    the 20240305 curated table and inherits its vetting honestly.  A chain
    that cannot be followed to a known set is reported as unvetted rather
    than guessed at.
    """

    path = Path(path)
    lineage = [path.name]
    seen = {path.name}
    current = path
    for _ in range(_MAX_LINEAGE_DEPTH):
        try:
            metadata = read_table_metadata(current)
        except OSError:
            break
        base = metadata.get("base wavelength file")
        if not base or base in seen:
            break
        lineage.append(base)
        seen.add(base)
        roots = [current.parent, current.parent.parent, *(search_dirs or ())]
        nxt = next((Path(root) / base for root in roots if (Path(root) / base).is_file()), None)
        if nxt is None:
            break
        current = nxt
    for name in lineage:
        if name in _VETTED_SETS:
            label, description = _VETTED_SETS[name]
            return TableVetting(path.name, label, tuple(lineage), description, name)
    return TableVetting(path.name, "", tuple(lineage), "")


def select_candidate_lines(
    lines: Sequence[CalibrationTableLine],
    species: Optional[Iterable[str]] = ("NeI",),
    require_ok: bool = True,
    min_width_px: float = 4.0,
    max_width_px: float = 40.0,
) -> List[CalibrationTableLine]:
    """Select reproducible candidate lamp lines for centroid fitting.

    The default keeps manually curated Ne I rows marked as OK/ok and removes
    very narrow or broad intervals that are likely bad picks or blends.
    """
    species_set = set(species) if species is not None else None
    selected: List[CalibrationTableLine] = []
    for line in lines:
        if species_set is not None and line.species not in species_set:
            continue
        if require_ok and "ok" not in line.comment.lower():
            continue
        if not (min_width_px <= line.width_px <= max_width_px):
            continue
        selected.append(line)
    return selected


def _line_window_bounds(
    spectrum_size: int,
    expected_center_px: float,
    window_radius_px: int,
) -> Tuple[int, int]:
    center_i = int(round(expected_center_px))
    lo = max(0, center_i - int(window_radius_px))
    hi = min(spectrum_size, center_i + int(window_radius_px) + 1)
    return lo, hi


def _window_noise(y: np.ndarray) -> Tuple[float, float]:
    edge_n = max(2, min(8, y.size // 5))
    edge = np.concatenate([y[:edge_n], y[-edge_n:]])
    edge_baseline = float(np.median(edge))
    edge_noise = float(np.std(edge - edge_baseline))

    lower_cut = float(np.percentile(y, 60))
    lower = y[y <= lower_cut]
    if lower.size >= 4:
        baseline = float(np.median(lower))
        noise = float(1.4826 * np.median(np.abs(lower - baseline)))
        if noise <= 0:
            noise = float(np.std(lower - baseline))
    else:
        baseline = edge_baseline
        noise = edge_noise

    if noise <= 0:
        noise = edge_noise
    if noise <= 0:
        noise = float(np.std(y - baseline))
    if noise <= 0:
        noise = 1.0
    return baseline, noise


def _measure_one_line_window(
    spectrum: Sequence[float],
    line: CalibrationTableLine,
    window_radius_px: int = 25,
    saturation_level: Optional[float] = None,
    min_snr: float = 5.0,
    max_saturated_fraction: float = 0.0,
) -> LineWindowStats:
    y_all = np.asarray(spectrum, dtype=float)
    if y_all.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    if y_all.size == 0:
        return LineWindowStats(
            line,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.0,
            0,
            0,
            0.0,
            False,
            "empty spectrum",
            saturation_level,
        )

    lo, hi = _line_window_bounds(y_all.size, line.center_pixel, window_radius_px)
    if hi - lo < 7:
        return LineWindowStats(
            line,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.0,
            0,
            0,
            0.0,
            False,
            "window too small",
            saturation_level,
        )

    x = np.arange(lo, hi, dtype=float)
    y = y_all[lo:hi]
    finite = np.isfinite(y)
    finite_count = int(finite.sum())
    if finite_count < 7:
        return LineWindowStats(
            line,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.0,
            finite_count,
            0,
            0.0,
            False,
            "too few finite pixels",
            saturation_level,
        )

    x = x[finite]
    y = y[finite]
    peak_idx = int(np.argmax(y))
    peak_pixel = float(x[peak_idx])
    peak_value = float(y[peak_idx])
    baseline, noise = _window_noise(y)
    prominence = peak_value - baseline
    snr = prominence / noise if noise > 0 else 0.0

    saturated_pixels = 0
    saturated_fraction = 0.0
    if saturation_level is not None:
        saturated_pixels = int(np.count_nonzero(y >= float(saturation_level)))
        saturated_fraction = saturated_pixels / float(finite_count)

    if saturated_fraction > max_saturated_fraction:
        reason = "saturated"
        fit_candidate = False
    elif prominence <= 0 or snr < min_snr:
        reason = "low snr"
        fit_candidate = False
    else:
        reason = ""
        fit_candidate = True

    return LineWindowStats(
        line,
        peak_pixel,
        peak_value,
        baseline,
        noise,
        prominence,
        snr,
        finite_count,
        saturated_pixels,
        saturated_fraction,
        fit_candidate,
        reason,
        saturation_level,
    )


def measure_line_window_stats(
    order_spectra: Sequence[Sequence[float]],
    candidate_lines: Sequence[CalibrationTableLine],
    window_radius_px: int = 25,
    saturation_level: Optional[float] = None,
    min_snr: float = 5.0,
    max_saturated_fraction: float = 0.0,
) -> List[LineWindowStats]:
    """Measure pre-fit line-window diagnostics for candidate rows."""
    results: List[LineWindowStats] = []
    for line in candidate_lines:
        if line.order_idx < 0 or line.order_idx >= len(order_spectra):
            results.append(
                LineWindowStats(
                    line,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.0,
                    0,
                    0,
                    0.0,
                    False,
                    "order missing",
                    saturation_level,
                )
            )
            continue

        results.append(
            _measure_one_line_window(
                order_spectra[line.order_idx],
                line,
                window_radius_px=window_radius_px,
                saturation_level=saturation_level,
                min_snr=min_snr,
                max_saturated_fraction=max_saturated_fraction,
            )
        )
    return results


def measure_detector_window_saturation(
    images: np.ndarray,
    pattern: np.ndarray,
    candidate_lines: Sequence[CalibrationTableLine],
    x_radius_px: int = 18,
    y_radius_px: int = 4,
    saturation_level: Optional[float] = None,
) -> List[DetectorWindowSaturation]:
    """Measure raw-detector saturation near expected calibration-line windows.

    ``images`` may be either ``(rows, cols)`` or ``(frames, rows, cols)``. The
    extracted 1D order spectra can exceed detector full scale because they are
    integrated traces, so saturation should be estimated from these raw image
    pixels whenever possible.
    """
    image_stack = np.asarray(images, dtype=float)
    if image_stack.ndim == 2:
        image_stack = image_stack[np.newaxis, :, :]
    if image_stack.ndim != 3:
        raise ValueError("images must have shape (rows, cols) or (frames, rows, cols)")

    pattern_arr = np.asarray(pattern, dtype=float)
    if pattern_arr.ndim != 2:
        raise ValueError("pattern must have shape (cols, orders)")

    _, n_rows, n_cols = image_stack.shape
    results: List[DetectorWindowSaturation] = []
    for line in candidate_lines:
        if line.order_idx < 0 or line.order_idx >= pattern_arr.shape[1]:
            results.append(
                DetectorWindowSaturation(
                    line,
                    np.nan,
                    0,
                    0,
                    0.0,
                    "order missing",
                    saturation_level,
                )
            )
            continue

        lo_x, hi_x = _line_window_bounds(
            min(n_cols, pattern_arr.shape[0]),
            line.center_pixel,
            x_radius_px,
        )
        window_chunks = []
        for x_i in range(lo_x, hi_x):
            y_i = int(round(pattern_arr[x_i, line.order_idx]))
            lo_y = max(0, y_i - int(y_radius_px))
            hi_y = min(n_rows, y_i + int(y_radius_px) + 1)
            if hi_y > lo_y:
                window_chunks.append(image_stack[:, lo_y:hi_y, x_i])

        if not window_chunks:
            results.append(
                DetectorWindowSaturation(
                    line,
                    np.nan,
                    0,
                    0,
                    0.0,
                    "empty detector window",
                    saturation_level,
                )
            )
            continue

        values = np.concatenate([chunk.reshape(-1) for chunk in window_chunks])
        finite = np.isfinite(values)
        finite_count = int(finite.sum())
        if finite_count == 0:
            results.append(
                DetectorWindowSaturation(
                    line,
                    np.nan,
                    0,
                    0,
                    0.0,
                    "too few finite pixels",
                    saturation_level,
                )
            )
            continue

        finite_values = values[finite]
        peak_value = float(np.nanmax(finite_values))
        saturated_pixels = 0
        saturated_fraction = 0.0
        if saturation_level is not None:
            saturated_pixels = int(np.count_nonzero(finite_values >= float(saturation_level)))
            saturated_fraction = saturated_pixels / float(finite_count)

        reason = "saturated" if saturated_pixels else ""
        results.append(
            DetectorWindowSaturation(
                line,
                peak_value,
                finite_count,
                saturated_pixels,
                saturated_fraction,
                reason,
                saturation_level,
            )
        )
    return results


def rank_candidate_lines(
    order_spectra: Sequence[Sequence[float]],
    candidate_lines: Sequence[CalibrationTableLine],
    window_radius_px: int = 25,
    saturation_level: Optional[float] = None,
    min_snr: float = 5.0,
    max_saturated_fraction: float = 0.0,
) -> List[LineWindowStats]:
    """Return candidate rows sorted by pre-fit quality diagnostics.

    Fit candidates are listed first, then saturated or low-SNR windows. Within
    each group, the strongest estimated SNR and prominence are listed first.
    """
    stats = measure_line_window_stats(
        order_spectra,
        candidate_lines,
        window_radius_px=window_radius_px,
        saturation_level=saturation_level,
        min_snr=min_snr,
        max_saturated_fraction=max_saturated_fraction,
    )
    return sorted(
        stats,
        key=lambda item: (
            item.fit_candidate,
            not item.is_saturated,
            np.nan_to_num(item.prominence, nan=-np.inf),
            np.nan_to_num(item.snr, nan=-np.inf),
        ),
        reverse=True,
    )


def _gaussian_with_linear_baseline(
    x: np.ndarray,
    amplitude: float,
    center: float,
    sigma: float,
    slope: float,
    intercept: float,
) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2) + slope * x + intercept


def fit_single_gaussian_centroid(
    spectrum: Sequence[float],
    expected_center_px: float,
    window_radius_px: int = 25,
    saturation_level: Optional[float] = None,
    min_snr: float = 5.0,
    min_sigma_px: float = 0.35,
    max_sigma_px: float = 8.0,
) -> Tuple[bool, float, float, float, float, float, str]:
    """Fit a single Gaussian centroid near ``expected_center_px``.

    Returns ``(success, center, sigma, amplitude, baseline, snr, reason)``.
    Failure returns finite best-effort numeric fields and a short reason.
    """
    y_all = np.asarray(spectrum, dtype=float)
    if y_all.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    if y_all.size == 0:
        return False, np.nan, np.nan, np.nan, np.nan, 0.0, "empty spectrum"

    lo, hi = _line_window_bounds(y_all.size, expected_center_px, window_radius_px)
    if hi - lo < 7:
        return False, np.nan, np.nan, np.nan, np.nan, 0.0, "window too small"

    x = np.arange(lo, hi, dtype=float)
    y = y_all[lo:hi]
    finite = np.isfinite(y)
    if finite.sum() < 7:
        return False, np.nan, np.nan, np.nan, np.nan, 0.0, "too few finite pixels"
    x = x[finite]
    y = y[finite]

    if saturation_level is not None and float(np.nanmax(y)) >= saturation_level:
        return False, float(x[np.nanargmax(y)]), np.nan, np.nan, np.nan, 0.0, "saturated"

    baseline0, noise = _window_noise(y)

    peak_idx = int(np.argmax(y))
    amp0 = float(y[peak_idx] - baseline0)
    snr0 = amp0 / noise
    if amp0 <= 0 or snr0 < min_snr:
        return False, float(x[peak_idx]), np.nan, amp0, baseline0, snr0, "low snr"

    p0 = [amp0, float(x[peak_idx]), 2.0, 0.0, baseline0]
    bounds = (
        [0.0, x[0], min_sigma_px, -np.inf, -np.inf],
        [np.inf, x[-1], max_sigma_px, np.inf, np.inf],
    )

    try:
        popt, _ = curve_fit(
            _gaussian_with_linear_baseline,
            x,
            y,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
        )
    except (RuntimeError, ValueError) as exc:
        return False, float(x[peak_idx]), np.nan, amp0, baseline0, snr0, str(exc)

    amplitude, center, sigma, slope, intercept = [float(v) for v in popt]
    baseline = slope * center + intercept
    snr = amplitude / noise
    if snr < min_snr:
        return False, center, sigma, amplitude, baseline, snr, "low fitted snr"
    if not (min_sigma_px <= sigma <= max_sigma_px):
        return False, center, sigma, amplitude, baseline, snr, "sigma out of range"

    return True, center, sigma, amplitude, baseline, snr, ""


def measure_line_centroids(
    order_spectra: Sequence[Sequence[float]],
    candidate_lines: Sequence[CalibrationTableLine],
    **fit_kwargs,
) -> List[LineCentroidFit]:
    """Fit centroids for candidate rows against extracted order spectra."""
    diagnostics = measure_line_window_stats(
        order_spectra,
        candidate_lines,
        window_radius_px=int(fit_kwargs.get("window_radius_px", 25)),
        saturation_level=fit_kwargs.get("saturation_level"),
        min_snr=float(fit_kwargs.get("min_snr", 5.0)),
    )
    results: List[LineCentroidFit] = []
    for line, stats in zip(candidate_lines, diagnostics):
        if line.order_idx < 0 or line.order_idx >= len(order_spectra):
            results.append(
                LineCentroidFit(
                    line,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.0,
                    False,
                    "order missing",
                    stats,
                )
            )
            continue

        ok, center, sigma, amp, baseline, snr, reason = fit_single_gaussian_centroid(
            order_spectra[line.order_idx],
            expected_center_px=line.center_pixel,
            **fit_kwargs,
        )
        results.append(LineCentroidFit(line, center, sigma, amp, baseline, snr, ok, reason, stats))
    return results


def _pattern_y_at(pattern: np.ndarray, order_idx: int, x_px: float) -> float:
    if order_idx < 0 or order_idx >= pattern.shape[1]:
        raise IndexError(f"order_idx {order_idx} outside pattern width {pattern.shape[1]}")
    x_i = int(round(x_px))
    x_i = max(0, min(pattern.shape[0] - 1, x_i))
    return float(pattern[x_i, order_idx])


def detector_points_from_lines(
    lines: Sequence[CalibrationTableLine],
    pattern: np.ndarray,
    measured_centers: Optional[Sequence[float]] = None,
) -> np.ndarray:
    """Return detector ``(x, y)`` points from lookup-table rows and pattern."""
    if measured_centers is not None and len(measured_centers) != len(lines):
        raise ValueError("measured_centers length must match lines length")

    points = []
    for i, line in enumerate(lines):
        x = float(measured_centers[i]) if measured_centers is not None else line.center_pixel
        points.append((x, _pattern_y_at(pattern, line.order_idx, x)))
    return np.asarray(points, dtype=float)


def fit_rigid_transform(
    expected_xy: np.ndarray,
    measured_xy: np.ndarray,
) -> Tuple[RigidTransform, float]:
    """Fit translation + rotation mapping ``expected_xy`` to ``measured_xy``.

    Returns the transform and RMS residual in detector pixels.
    """
    expected = np.asarray(expected_xy, dtype=float)
    measured = np.asarray(measured_xy, dtype=float)
    if expected.shape != measured.shape or expected.ndim != 2 or expected.shape[1] != 2:
        raise ValueError("expected_xy and measured_xy must both have shape (N, 2)")
    if expected.shape[0] < 2:
        raise ValueError("at least two points are required for a rigid transform")

    src_centroid = expected.mean(axis=0)
    dst_centroid = measured.mean(axis=0)
    src0 = expected - src_centroid
    dst0 = measured - dst_centroid
    u, _, vt = np.linalg.svd(src0.T @ dst0)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    theta = float(np.arctan2(rotation[1, 0], rotation[0, 0]))
    translation = dst_centroid - src_centroid @ rotation.T
    transform = RigidTransform(float(translation[0]), float(translation[1]), theta)

    residual = transform.apply(expected) - measured
    rms_px = float(np.sqrt(np.mean(np.sum(residual**2, axis=1))))
    return transform, rms_px


def _moved_line(line: CalibrationTableLine, dx_px: float) -> CalibrationTableLine:
    """Return one row slid along the dispersion axis, keeping its wavelength."""
    return CalibrationTableLine(
        order_idx=line.order_idx,
        pixel_from=line.pixel_from + dx_px,
        pixel_to=line.pixel_to + dx_px,
        center_pixel=line.center_pixel + dx_px,
        wavelength_nm=line.wavelength_nm,
        species=line.species,
        comment=line.comment,
    )


def apply_rigid_correction_to_lines(
    lines: Sequence[CalibrationTableLine],
    pattern: np.ndarray,
    transform: RigidTransform,
) -> List[CalibrationTableLine]:
    """Return new lookup rows whose center pixels are moved by ``transform``."""
    expected = detector_points_from_lines(lines, pattern)
    corrected = transform.apply(expected)
    return [
        _moved_line(line, float(point[0] - line.center_pixel))
        for line, point in zip(lines, corrected)
    ]


def apply_rigid_correction_to_pattern(
    pattern: np.ndarray,
    transform: RigidTransform,
) -> np.ndarray:
    """Return the order pattern moved by ``transform``, still one row per column.

    The pattern is the wavelength table's other half: the table says which
    detector column a line sits in, the pattern says which detector row that
    column's order runs along.  ``apply_rigid_correction_to_lines`` moves the
    table's points; this moves the pattern's, so a caller that corrects one
    can correct the other with the same solved transform and end up with two
    files describing the same detector.

    Column ``x`` of order ``o`` carries the point ``(x, pattern[x, o])``.  The
    transform moves that point to ``(x', y')`` — a *point*, not a column — so
    the corrected trace is resampled back onto the integer column grid rather
    than written at the column it started in: that is what keeps the corrected
    pattern consistent with the corrected table, whose rows are likewise
    re-read at their moved ``x``.

    Columns whose corrected trace has no source column — the few at either
    detector edge that the transform slid out of, or in from beyond — hold the
    nearest corrected trace row rather than extrapolating a slope off the end
    of the measurement.  A trace never wraps and never runs away at the edges;
    it goes flat there, which the extraction's own band half-width absorbs.

    Rows are not clipped to the detector height.  A pattern that the transform
    pushes past the top edge stays a truthful description of where the order
    went, and ``Calibrations.make_cutting_masks`` already trims the orders that
    leave the frame; silently pulling the trace back inside would move it off
    its own light instead.
    """

    values = np.asarray(pattern)
    if values.ndim != 2:
        raise ValueError("pattern must be a 2-D (columns, orders) array")
    columns = np.arange(values.shape[0], dtype=float)
    corrected = np.empty(values.shape, dtype=float)
    for order_idx in range(values.shape[1]):
        points = np.column_stack([columns, values[:, order_idx].astype(float)])
        moved = transform.apply(points)
        ordering = np.argsort(moved[:, 0], kind="stable")
        corrected[:, order_idx] = np.interp(
            columns, moved[ordering, 0], moved[ordering, 1]
        )
    if np.issubdtype(values.dtype, np.integer):
        return np.rint(corrected).astype(values.dtype)
    return corrected


def write_pattern_table(
    pattern: np.ndarray,
    path: str | Path,
    metadata: Optional[Sequence[Tuple[str, str]]] = None,
) -> None:
    """Write an order pattern as the integer row-per-column table readers expect.

    ``Calibrations.load_pattern`` reads this back with ``dtype=int``, so the
    rows are rounded here rather than left as the floats the correction
    produced.  The ``#`` header is provenance only; every reader in this
    package parses the file with :func:`numpy.loadtxt`, which skips it.
    """

    values = np.rint(np.asarray(pattern, dtype=float)).astype(int)
    if values.ndim != 2:
        raise ValueError("pattern must be a 2-D (columns, orders) array")
    out = Path(path)
    with out.open("w", newline="\n") as fh:
        fh.write("# Adjusted echelle order pattern\n")
        if metadata:
            for key, value in metadata:
                fh.write(f"# {key}: {value}\n")
        fh.write("# one row per detector column, one column per diffraction order\n")
        for row in values:
            fh.write(" ".join(str(int(item)) for item in row) + "\n")


@dataclass(frozen=True)
class PatternCorrection:
    """What a saved ``pattern.txt`` is, and why it is that."""

    applied: bool
    reason: str
    max_shift_px: float
    order_count: int


#: Below this the correction would only reformat the pattern, never move a trace.
IDENTITY_PATTERN_SHIFT_PX = 1e-6


def write_corrected_pattern_table(
    base_pattern: str | Path,
    destination: str | Path,
    *,
    transform: RigidTransform,
    metadata: Sequence[Tuple[str, str]] = (),
    identity_shift_px: float = IDENTITY_PATTERN_SHIFT_PX,
) -> PatternCorrection:
    """Write the base order pattern moved by *transform*, or copy it unchanged.

    Both writers of a corrected calibration — the bench's snapshot save and
    ``echelle-align --save`` — go through this one function, so the pattern
    beside a corrected wavelength table always speaks the same detector frame
    as that table, whichever surface wrote it.

    A transform that moves no trace by a measurable row would only reformat
    the base pattern, so its bytes are copied instead and the returned reason
    states which of the two outcomes the caller got.
    """

    source = Path(base_pattern)
    target = Path(destination)
    try:
        with warnings.catch_warnings():
            # An empty table is a state to report, never a warning to raise.
            warnings.simplefilter("ignore")
            values = np.loadtxt(source, dtype=int, ndmin=2)
    except (ValueError, IndexError, StopIteration):
        values = np.empty((0, 0), dtype=int)
    if values.size == 0:
        shutil.copy2(source, target)
        return PatternCorrection(
            False,
            f"{source.name} holds no order columns this correction can read, "
            "so its bytes were copied unchanged",
            0.0,
            0,
        )
    corrected = apply_rigid_correction_to_pattern(values, transform)
    max_shift_px = float(np.max(np.abs(corrected.astype(float) - values.astype(float))))
    if max_shift_px <= float(identity_shift_px):
        shutil.copy2(source, target)
        return PatternCorrection(
            False,
            f"the solved transform moves no trace of {source.name} measurably "
            f"({max_shift_px:.3g} px), so its bytes were copied unchanged",
            max_shift_px,
            int(values.shape[1]),
        )
    write_pattern_table(corrected, target, metadata=metadata)
    return PatternCorrection(
        True,
        f"{values.shape[1]} order traces of {source.name} were moved by the solved "
        f"transform (largest shift {max_shift_px:.3f} px)",
        max_shift_px,
        int(values.shape[1]),
    )


def shift_lines_in_pixels(
    lines: Sequence[CalibrationTableLine], dx_px: float
) -> List[CalibrationTableLine]:
    """Return lookup rows translated by ``dx_px`` along the dispersion axis.

    This is the pure-translation case of :func:`apply_rigid_correction_to_lines`
    and moves every row by the same amount, so no order pattern is needed: with
    no rotation, a row's detector ``y`` never enters its corrected ``x``.  A
    caller that measured one rigid detector shift — the drift audit, say — uses
    this instead of building a synthetic pattern for the general routine.
    """
    return [_moved_line(line, float(dx_px)) for line in lines]


def _toml_quote(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def write_wavelength_table(
    lines: Sequence[CalibrationTableLine],
    path: str | Path,
    metadata: Optional[Sequence[Tuple[str, str]]] = None,
) -> None:
    """Write an adjusted wavelength table without mutating the original file."""
    out = Path(path)
    with out.open("w", newline="\n") as fh:
        fh.write("# Adjusted wavelength calibration lookup table\n")
        if metadata:
            for key, value in metadata:
                fh.write(f"# {key}: {value}\n")
        fh.write("# order from to center wavelength band\n")
        for line in lines:
            comment = f"  # {line.comment}" if line.comment else ""
            fh.write(
                f"{line.order_idx:d}\t"
                f"{line.pixel_from:08.3f}\t{line.pixel_to:08.3f}\t"
                f"{line.center_pixel:010.4f}\t{line.wavelength_nm:010.5f}\t"
                f"{line.species}{comment}\n"
            )


def save_alignment_settings(settings: AlignmentSettings, path: str | Path) -> None:
    """Persist alignment settings as small TOML."""
    out = Path(path)
    with out.open("w", newline="\n") as fh:
        fh.write("# Echelle calibration alignment settings\n")
        fh.write(f'instrument_id = "{_toml_quote(settings.instrument_id)}"\n')
        fh.write(f'created_at = "{_toml_quote(settings.created_at)}"\n')
        fh.write(f'alignment_dataset_id = "{_toml_quote(settings.alignment_dataset_id)}"\n')
        fh.write(f'alignment_source_dir = "{_toml_quote(settings.alignment_source_dir)}"\n')
        fh.write(f'alignment_lamp = "{_toml_quote(settings.alignment_lamp)}"\n')
        fh.write(f'signal_file = "{_toml_quote(settings.signal_file)}"\n')
        fh.write(f'background_file = "{_toml_quote(settings.background_file)}"\n')
        fh.write(f'base_wavelength_file = "{_toml_quote(settings.base_wavelength_file)}"\n')
        fh.write(f'base_pattern_file = "{_toml_quote(settings.base_pattern_file)}"\n')
        fh.write(f'sphere_file = "{_toml_quote(settings.sphere_file)}"\n')
        fh.write(f'sphere_background_file = "{_toml_quote(settings.sphere_background_file)}"\n')
        fh.write(f'output_wavelength_file = "{_toml_quote(settings.output_wavelength_file)}"\n')
        fh.write(f'output_pattern_file = "{_toml_quote(settings.output_pattern_file)}"\n')
        fh.write(f"n_lines = {settings.n_lines:d}\n")
        fh.write(f"rms_px = {settings.rms_px:.10g}\n")
        fh.write(f'notes = "{_toml_quote(settings.notes)}"\n')
        fh.write("\n[transform]\n")
        fh.write(f"dx_px = {settings.transform.dx_px:.12g}\n")
        fh.write(f"dy_px = {settings.transform.dy_px:.12g}\n")
        fh.write(f"theta_rad = {settings.transform.theta_rad:.12g}\n")


def load_alignment_settings(path: str | Path) -> AlignmentSettings:
    """Load TOML written by :func:`save_alignment_settings`."""
    with Path(path).open("rb") as fh:
        data = tomllib.load(fh)
    transform_data = data["transform"]
    return AlignmentSettings(
        instrument_id=str(data["instrument_id"]),
        base_wavelength_file=str(data["base_wavelength_file"]),
        n_lines=int(data["n_lines"]),
        rms_px=float(data["rms_px"]),
        created_at=str(data.get("created_at", "")),
        alignment_dataset_id=str(data.get("alignment_dataset_id", "")),
        alignment_source_dir=str(data.get("alignment_source_dir", "")),
        alignment_lamp=str(data.get("alignment_lamp", "")),
        signal_file=str(data.get("signal_file", "")),
        background_file=str(data.get("background_file", "")),
        base_pattern_file=str(data.get("base_pattern_file", "")),
        sphere_file=str(data.get("sphere_file", "")),
        sphere_background_file=str(data.get("sphere_background_file", "")),
        output_wavelength_file=str(data.get("output_wavelength_file", "")),
        output_pattern_file=str(data.get("output_pattern_file", "")),
        notes=str(data.get("notes", "")),
        transform=RigidTransform(
            dx_px=float(transform_data["dx_px"]),
            dy_px=float(transform_data["dy_px"]),
            theta_rad=float(transform_data["theta_rad"]),
        ),
    )


def _line_identity(line: CalibrationTableLine) -> Tuple[int, float, float]:
    return (line.order_idx, line.center_pixel, line.wavelength_nm)


def run_calibration_alignment(config: AlignmentRunConfig) -> AlignmentRunResult:
    """Run the headless Neon wavelength-alignment workflow.

    This packages the current notebook path into a reusable function:

    1. Load calibrations with the selected wavelength table and order pattern.
    2. Extract signal/background order spectra.
    3. Rank curated candidate rows and reject raw-detector saturated windows.
    4. Fit single-Gaussian centroids and a rigid detector transform.
    5. Build adjusted wavelength rows and an :class:`AlignmentSettings` summary.

    The function does not write files. Use :func:`save_alignment_settings` and
    :func:`write_wavelength_table` after reviewing the returned residuals.
    """
    from .echelle import Calibrations, EchelleImage

    calibration_dir = Path(config.calibration_dir)
    files_cmos = {
        "orders": config.base_pattern_file,
        "wavelength": config.base_wavelength_file,
        "sphr": str(config.sphere_file),
        "bkgr": str(config.sphere_background_file),
        "integral": config.integral_file,
    }
    calibration = Calibrations(folder=str(calibration_dir), filenames=files_cmos)
    calibration.start()

    signal = EchelleImage(str(config.signal_file), clbr=calibration)
    signal.calculate_order_spectra()
    signal.correct_order_shapes()

    background = EchelleImage(str(config.background_file), clbr=calibration)
    background.calculate_order_spectra()
    background.correct_order_shapes()

    order_spectra = (
        np.asarray(signal.order_spectra[0], dtype=float)
        - np.asarray(background.order_spectra[0], dtype=float)
    )

    rows = load_wavelength_table(calibration_dir / config.base_wavelength_file)
    candidates = select_candidate_lines(
        rows,
        species=config.species,
        require_ok=config.require_ok,
        min_width_px=config.min_width_px,
        max_width_px=config.max_width_px,
    )
    ranked_stats = rank_candidate_lines(
        order_spectra,
        candidates,
        window_radius_px=config.window_radius_px,
        saturation_level=None,
        min_snr=config.min_snr,
    )
    detector_saturation = measure_detector_window_saturation(
        signal.images,
        calibration.pattern,
        candidates,
        x_radius_px=config.x_radius_px,
        y_radius_px=config.y_radius_px,
        saturation_level=config.saturation_level,
    )
    saturation_by_key = {_line_identity(stat.line): stat for stat in detector_saturation}

    fit_candidates = [stat.line for stat in ranked_stats]
    fits = measure_line_centroids(
        order_spectra,
        fit_candidates,
        window_radius_px=config.window_radius_px,
        saturation_level=None,
        min_snr=config.min_snr,
    )
    good_fits = [
        fit
        for fit in fits
        if fit.success and not saturation_by_key[_line_identity(fit.line)].is_saturated
    ]
    if len(good_fits) < 2:
        raise ValueError(
            f"Need at least two successful non-saturated fits, got {len(good_fits)}"
        )

    good_lines = [fit.line for fit in good_fits]
    measured_centers = [fit.center_pixel for fit in good_fits]
    expected_xy = detector_points_from_lines(good_lines, calibration.pattern)
    measured_xy = detector_points_from_lines(good_lines, calibration.pattern, measured_centers)
    transform, rms_px = fit_rigid_transform(expected_xy, measured_xy)
    predicted_xy = transform.apply(expected_xy)
    residual_xy = predicted_xy - measured_xy
    residual_px = np.sqrt(np.sum(residual_xy**2, axis=1))
    adjusted_rows = apply_rigid_correction_to_lines(rows, calibration.pattern, transform)

    settings = AlignmentSettings(
        instrument_id=config.instrument_id,
        created_at=config.created_at,
        alignment_dataset_id=config.alignment_dataset_id,
        alignment_source_dir=config.alignment_source_dir,
        alignment_lamp=config.alignment_lamp,
        signal_file=Path(config.signal_file).name,
        background_file=Path(config.background_file).name,
        base_wavelength_file=config.base_wavelength_file,
        base_pattern_file=config.base_pattern_file,
        sphere_file=Path(config.sphere_file).name,
        sphere_background_file=Path(config.sphere_background_file).name,
        output_wavelength_file=config.output_wavelength_file,
        output_pattern_file=config.output_pattern_file,
        transform=transform,
        n_lines=len(good_fits),
        rms_px=rms_px,
        notes=config.notes,
    )
    return AlignmentRunResult(
        settings=settings,
        adjusted_rows=adjusted_rows,
        rows=rows,
        candidates=candidates,
        fits=fits,
        good_fits=good_fits,
        residual_xy=residual_xy,
        residual_px=residual_px,
        expected_xy=expected_xy,
        measured_xy=measured_xy,
        predicted_xy=predicted_xy,
        detector_saturation=detector_saturation,
        ranked_stats=ranked_stats,
    )
