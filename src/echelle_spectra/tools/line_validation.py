"""Known-line wavelength validation helpers.

This module is intentionally narrow: it validates measured peak positions
against supplied line wavelengths. Molecular modeling and intensity analysis
belong in downstream packages such as ``fulcher_analyzer``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import curve_fit

try:
    import pandas as pd
except ModuleNotFoundError:  # pragma: no cover - pandas is a package dependency
    pd = None

__all__ = [
    "LineValidationTarget",
    "LineValidationResult",
    "LineValidationSummary",
    "balmer_air_targets",
    "build_stitched_order_index",
    "fit_validation_line",
    "load_fulcher_h2_q_branch_targets",
    "summarize_validation",
    "validate_lines",
]


@dataclass(frozen=True)
class LineValidationTarget:
    """One expected validation line."""

    label: str
    wavelength_nm: float
    wavelength_medium: str = "air"
    notes: str = ""


@dataclass(frozen=True)
class LineValidationResult:
    """Measured centroid and diagnostics for one validation line."""

    line: str
    expected_nm: float
    measured_nm: float
    residual_nm: float
    order: int | None
    peak_snr: float
    notes: str = ""
    sigma_nm: float | None = None
    frames_fit: int | None = None
    frames_total: int | None = None
    frame_centroid_std_nm: float | None = None


@dataclass(frozen=True)
class LineValidationSummary:
    """Compact residual summary."""

    n_lines: int
    mean_residual_nm: float
    median_residual_nm: float
    rms_residual_nm: float


def _finite_mask(values: np.ndarray) -> np.ndarray:
    if pd is not None:
        return ~pd.isnull(values)
    return np.isfinite(values)


def _gaussian_with_slope(
    x: np.ndarray,
    amplitude: float,
    center: float,
    sigma: float,
    baseline: float,
    slope: float,
) -> np.ndarray:
    return baseline + slope * (x - center) + amplitude * np.exp(
        -0.5 * ((x - center) / sigma) ** 2
    )


def _robust_noise(values: np.ndarray) -> tuple[float, float]:
    baseline = float(np.nanmedian(values))
    noise = float(1.4826 * np.nanmedian(np.abs(values - baseline)))
    if not np.isfinite(noise) or noise <= 0:
        noise = float(np.nanstd(values))
    if not np.isfinite(noise) or noise <= 0:
        noise = 1.0
    return baseline, noise


def balmer_air_targets() -> list[LineValidationTarget]:
    """Return common Balmer validation targets in air wavelengths."""
    return [
        LineValidationTarget("H-alpha", 656.279, "air"),
        LineValidationTarget("H-beta", 486.135, "air"),
        LineValidationTarget("H-gamma", 434.047, "air"),
    ]


def load_fulcher_h2_q_branch_targets(
    path: str | Path,
    *,
    wavelength_medium: str = "air",
) -> list[LineValidationTarget]:
    """Load H2 Fulcher-alpha Q-branch wavelength targets from a table.

    The bundled ``fulcher_analyzer`` table is arranged as Q1..Q11 rows and
    0-0..3-3 vibrational-band columns. Zero entries are omitted.
    """
    table = np.loadtxt(path, comments="#")
    bands = ("0-0", "1-1", "2-2", "3-3")
    targets: list[LineValidationTarget] = []
    for band_index, band in enumerate(bands):
        for q_index, wavelength in enumerate(table[:, band_index], start=1):
            if wavelength > 0:
                targets.append(
                    LineValidationTarget(
                        label=f"Fulcher H2 Q{q_index}({band})",
                        wavelength_nm=float(wavelength),
                        wavelength_medium=wavelength_medium,
                    )
                )
    return targets


def build_stitched_order_index(calibration) -> np.ndarray:
    """Return order indices aligned with ``Spectrum.wavelength``.

    ``Calibrations.wavelength`` includes NaN slots for clipped order tails;
    ``Spectrum.wavelength`` drops those slots and may reverse the stitched
    vector. This helper mirrors those transformations.
    """
    order_grid = np.broadcast_to(
        np.arange(calibration.order_wavel.shape[0])[:, None],
        calibration.order_wavel.shape,
    )
    order_index = order_grid[calibration.order_borders]
    order_index = order_index[_finite_mask(calibration.wavelength)]
    if calibration.direction < 0:
        order_index = np.flip(order_index)
    return order_index.astype(int, copy=False)


def fit_validation_line(
    wavelength_nm: np.ndarray,
    intensity: np.ndarray,
    target: LineValidationTarget,
    *,
    order_index: np.ndarray | None = None,
    window_nm: float = 0.30,
    min_points: int = 8,
    min_snr: float = 5.0,
) -> LineValidationResult | None:
    """Fit one line centroid in a local wavelength window."""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    mask = (
        np.isfinite(wavelength_nm)
        & np.isfinite(intensity)
        & (np.abs(wavelength_nm - target.wavelength_nm) <= window_nm)
    )
    if int(mask.sum()) < min_points:
        return None

    x0 = wavelength_nm[mask]
    y0 = intensity[mask]
    selected_order: int | None = None
    if order_index is not None:
        orders = np.asarray(order_index)[mask]
        best: tuple[float, int, np.ndarray] | None = None
        for order in np.unique(orders):
            keep = orders == order
            if int(keep.sum()) < min_points:
                continue
            score = float(np.nanmax(y0[keep]) - np.nanmedian(y0[keep]))
            if best is None or score > best[0]:
                best = (score, int(order), keep)
        if best is not None:
            selected_order = best[1]
            x0 = x0[best[2]]
            y0 = y0[best[2]]

    if len(x0) < min_points:
        return None

    sort_index = np.argsort(x0)
    x = x0[sort_index]
    y = y0[sort_index]
    edge_count = max(2, min(8, len(y) // 4))
    edge_values = np.r_[y[:edge_count], y[-edge_count:]]
    baseline0, noise = _robust_noise(edge_values)
    peak_index = int(np.nanargmax(y))
    amplitude0 = float(y[peak_index] - baseline0)
    center0 = float(x[peak_index])

    fit_ok = True
    try:
        popt, _ = curve_fit(
            _gaussian_with_slope,
            x,
            y,
            p0=[max(amplitude0, noise), center0, 0.035, baseline0, 0.0],
            bounds=(
                [0.0, target.wavelength_nm - window_nm, 0.004, -np.inf, -np.inf],
                [np.inf, target.wavelength_nm + window_nm, 0.18, np.inf, np.inf],
            ),
            maxfev=20000,
        )
        amplitude, center, sigma, _baseline, _slope = [float(value) for value in popt]
    except Exception:
        amplitude = amplitude0
        center = center0
        sigma = np.nan
        fit_ok = False

    snr = float(amplitude / noise) if noise else float("nan")
    if snr < min_snr:
        return None

    notes = target.notes
    if not fit_ok:
        notes = "; ".join(part for part in (notes, "fallback peak") if part)
    return LineValidationResult(
        line=target.label,
        expected_nm=float(target.wavelength_nm),
        measured_nm=float(center),
        residual_nm=float(center - target.wavelength_nm),
        order=selected_order,
        peak_snr=snr,
        notes=notes,
        sigma_nm=float(sigma) if np.isfinite(sigma) else None,
    )


def validate_lines(
    wavelength_nm: np.ndarray,
    spectra: np.ndarray,
    targets: Sequence[LineValidationTarget],
    *,
    signal_frames: Sequence[int] | None = None,
    order_index: np.ndarray | None = None,
    window_nm: float = 0.30,
    min_snr: float = 5.0,
) -> list[LineValidationResult]:
    """Validate many targets against an aggregate spectrum.

    For multiframe spectra, the aggregate is the mean over *signal_frames*.
    Per-frame fits are attempted afterward and reported as diagnostics.
    """
    spectra = np.asarray(spectra, dtype=float)
    if spectra.ndim == 1:
        spectra = spectra[np.newaxis, :]
    frame_index = (
        np.asarray(signal_frames, dtype=int)
        if signal_frames is not None
        else np.arange(spectra.shape[0], dtype=int)
    )
    aggregate = np.nanmean(spectra[frame_index], axis=0)
    results: list[LineValidationResult] = []
    for target in targets:
        result = fit_validation_line(
            wavelength_nm,
            aggregate,
            target,
            order_index=order_index,
            window_nm=window_nm,
            min_snr=min_snr,
        )
        if result is None:
            continue

        frame_centers: list[float] = []
        for frame in frame_index:
            frame_result = fit_validation_line(
                wavelength_nm,
                spectra[frame],
                target,
                order_index=order_index,
                window_nm=window_nm,
                min_snr=min_snr,
            )
            if frame_result is not None:
                frame_centers.append(frame_result.measured_nm)

        frame_std = float(np.std(frame_centers)) if len(frame_centers) > 1 else None
        frame_note = f"{len(frame_centers)}/{len(frame_index)} signal frames fit"
        if frame_std is not None:
            frame_note += f"; frame sd {frame_std:.4f} nm"
        notes = "; ".join(part for part in (result.notes, frame_note) if part)
        results.append(
            LineValidationResult(
                line=result.line,
                expected_nm=result.expected_nm,
                measured_nm=result.measured_nm,
                residual_nm=result.residual_nm,
                order=result.order,
                peak_snr=result.peak_snr,
                notes=notes,
                sigma_nm=result.sigma_nm,
                frames_fit=len(frame_centers),
                frames_total=len(frame_index),
                frame_centroid_std_nm=frame_std,
            )
        )
    return results


def summarize_validation(results: Iterable[LineValidationResult]) -> LineValidationSummary:
    """Return mean, median, and RMS residuals for validation results."""
    residuals = np.asarray([result.residual_nm for result in results], dtype=float)
    if residuals.size == 0:
        return LineValidationSummary(0, float("nan"), float("nan"), float("nan"))
    return LineValidationSummary(
        n_lines=int(residuals.size),
        mean_residual_nm=float(np.mean(residuals)),
        median_residual_nm=float(np.median(residuals)),
        rms_residual_nm=float(np.sqrt(np.mean(residuals**2))),
    )
