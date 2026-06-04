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
    "LineCentroidFit",
    "RigidTransform",
    "AlignmentSettings",
    "load_wavelength_table",
    "select_candidate_lines",
    "fit_single_gaussian_centroid",
    "measure_line_centroids",
    "fit_rigid_transform",
    "detector_points_from_lines",
    "apply_rigid_correction_to_lines",
    "write_wavelength_table",
    "save_alignment_settings",
    "load_alignment_settings",
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
    notes: str = ""


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

    center_i = int(round(expected_center_px))
    lo = max(0, center_i - int(window_radius_px))
    hi = min(y_all.size, center_i + int(window_radius_px) + 1)
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

    edge_n = max(2, min(8, y.size // 5))
    edge = np.concatenate([y[:edge_n], y[-edge_n:]])
    baseline0 = float(np.median(edge))
    noise = float(np.std(edge - baseline0))
    if noise <= 0:
        noise = float(np.std(y - np.median(y)))
    if noise <= 0:
        noise = 1.0

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
    results: List[LineCentroidFit] = []
    for line in candidate_lines:
        if line.order_idx < 0 or line.order_idx >= len(order_spectra):
            results.append(
                LineCentroidFit(line, np.nan, np.nan, np.nan, np.nan, 0.0, False, "order missing")
            )
            continue

        ok, center, sigma, amp, baseline, snr, reason = fit_single_gaussian_centroid(
            order_spectra[line.order_idx],
            expected_center_px=line.center_pixel,
            **fit_kwargs,
        )
        results.append(LineCentroidFit(line, center, sigma, amp, baseline, snr, ok, reason))
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


def apply_rigid_correction_to_lines(
    lines: Sequence[CalibrationTableLine],
    pattern: np.ndarray,
    transform: RigidTransform,
) -> List[CalibrationTableLine]:
    """Return new lookup rows whose center pixels are moved by ``transform``."""
    expected = detector_points_from_lines(lines, pattern)
    corrected = transform.apply(expected)
    adjusted: List[CalibrationTableLine] = []
    for line, point in zip(lines, corrected):
        dx = float(point[0] - line.center_pixel)
        adjusted.append(
            CalibrationTableLine(
                order_idx=line.order_idx,
                pixel_from=line.pixel_from + dx,
                pixel_to=line.pixel_to + dx,
                center_pixel=float(point[0]),
                wavelength_nm=line.wavelength_nm,
                species=line.species,
                comment=line.comment,
            )
        )
    return adjusted


def write_wavelength_table(lines: Sequence[CalibrationTableLine], path: str | Path) -> None:
    """Write an adjusted wavelength table without mutating the original file."""
    out = Path(path)
    with out.open("w", newline="\n") as fh:
        fh.write("# Adjusted wavelength calibration lookup table\n")
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
        fh.write(f'instrument_id = "{settings.instrument_id}"\n')
        fh.write(f'base_wavelength_file = "{settings.base_wavelength_file}"\n')
        fh.write(f"n_lines = {settings.n_lines:d}\n")
        fh.write(f"rms_px = {settings.rms_px:.10g}\n")
        fh.write(f'notes = "{settings.notes.replace(chr(34), chr(39))}"\n')
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
        notes=str(data.get("notes", "")),
        transform=RigidTransform(
            dx_px=float(transform_data["dx_px"]),
            dy_px=float(transform_data["dy_px"]),
            theta_rad=float(transform_data["theta_rad"]),
        ),
    )
