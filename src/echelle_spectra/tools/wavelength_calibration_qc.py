"""QC helpers for Echelle wavelength-calibration lookup tables.

The calibration code fits one polynomial per order from manually identified
line centers. These helpers reproduce that fit, write residual tables, and make
review plots for spotting bad orders or conflicting line anchors.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .calibration_alignment import CalibrationTableLine, load_wavelength_table

__all__ = [
    "WavelengthTable",
    "OrderWavelengthFit",
    "WavelengthQcResult",
    "fit_wavelength_orders",
    "find_focus_orders",
    "run_wavelength_calibration_qc",
    "write_qc_tables",
    "plot_all_orders",
    "plot_focus_orders",
    "plot_order_residuals",
    "plot_focus_line_residuals",
]


@dataclass(frozen=True)
class WavelengthTable:
    """Named wavelength lookup table used in a QC run."""

    name: str
    path: Path
    lines: tuple[CalibrationTableLine, ...]


@dataclass(frozen=True)
class OrderWavelengthFit:
    """Polynomial wavelength fit and residual diagnostics for one order."""

    table: str
    order: int
    degree: int
    n_lines: int
    coeffs: tuple[float, ...]
    rms_residual_nm: float
    max_abs_residual_nm: float
    wavelength_min_nm: float
    wavelength_max_nm: float
    slope_mid_nm_per_px: float

    @property
    def poly(self) -> np.poly1d:
        """Return the fitted polynomial."""
        return np.poly1d(np.asarray(self.coeffs, dtype=float))


@dataclass(frozen=True)
class WavelengthQcResult:
    """Complete in-memory QC result."""

    tables: tuple[WavelengthTable, ...]
    fits: Mapping[str, Mapping[int, OrderWavelengthFit]]
    detector_width_px: int


def _fit_degree(n_lines: int) -> int:
    return 1 if n_lines < 3 else 2


def _line_order(line: CalibrationTableLine) -> int:
    return int(line.order_idx)


def _line_center(line: CalibrationTableLine) -> float:
    return float(line.center_pixel)


def _line_wavelength(line: CalibrationTableLine) -> float:
    return float(line.wavelength_nm)


def _line_species(line: CalibrationTableLine) -> str:
    return str(line.species)


def _short_label(line: CalibrationTableLine) -> str:
    return f"{_line_species(line)} {_line_wavelength(line):.3f}"


def fit_wavelength_orders(
    table_name: str,
    lines: Sequence[CalibrationTableLine],
    *,
    detector_width_px: int = 2560,
) -> dict[int, OrderWavelengthFit]:
    """Fit the per-order wavelength polynomials used by ``Calibrations``."""
    fits: dict[int, OrderWavelengthFit] = {}
    x_eval = np.arange(detector_width_px, dtype=float)
    for order in sorted({_line_order(line) for line in lines}):
        order_lines = [line for line in lines if _line_order(line) == order]
        centers = np.asarray([_line_center(line) for line in order_lines], dtype=float)
        wavelengths = np.asarray([_line_wavelength(line) for line in order_lines], dtype=float)
        degree = _fit_degree(len(order_lines))
        coeffs = np.polyfit(centers, wavelengths, degree)
        poly = np.poly1d(coeffs)
        residual = wavelengths - poly(centers)
        fit_wavelength = poly(x_eval)
        derivative = np.polyder(poly)
        fits[order] = OrderWavelengthFit(
            table=table_name,
            order=order,
            degree=degree,
            n_lines=len(order_lines),
            coeffs=tuple(float(value) for value in coeffs),
            rms_residual_nm=float(np.sqrt(np.mean(residual**2))),
            max_abs_residual_nm=float(np.max(np.abs(residual))),
            wavelength_min_nm=float(np.nanmin(fit_wavelength)),
            wavelength_max_nm=float(np.nanmax(fit_wavelength)),
            slope_mid_nm_per_px=float(derivative(detector_width_px / 2)),
        )
    return fits


def find_focus_orders(
    fits: Mapping[int, OrderWavelengthFit],
    *,
    min_nm: float,
    max_nm: float,
) -> list[int]:
    """Return orders whose fitted wavelength span overlaps a target range."""
    return [
        order
        for order, fit in fits.items()
        if fit.wavelength_min_nm <= max_nm and fit.wavelength_max_nm >= min_nm
    ]


def _lines_for_order(
    lines: Sequence[CalibrationTableLine],
    order: int,
    *,
    min_nm: float | None = None,
    max_nm: float | None = None,
) -> list[CalibrationTableLine]:
    selected = [line for line in lines if _line_order(line) == order]
    if min_nm is not None:
        selected = [line for line in selected if _line_wavelength(line) >= min_nm]
    if max_nm is not None:
        selected = [line for line in selected if _line_wavelength(line) <= max_nm]
    return selected


def _residual_nm(line: CalibrationTableLine, fit: OrderWavelengthFit) -> float:
    return float(_line_wavelength(line) - fit.poly(_line_center(line)))


def _annotate_decluttered(
    ax,
    points: Sequence[tuple[float, float, str, str]],
    *,
    color: str,
    max_labels: int | None = None,
) -> None:
    """Add small labels with alternating offsets and leader lines.

    This is intentionally simple and dependency-free. The goal is to identify a
    few important points without turning dense calibration regions into text.
    """
    if max_labels is not None:
        points = points[:max_labels]
    offsets = [(8, 8), (8, -12), (-8, 10), (-8, -14), (16, 0), (-16, 0)]
    for index, (x, y, label, _key) in enumerate(points):
        dx, dy = offsets[index % len(offsets)]
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
            fontsize=7,
            color=color,
            arrowprops={"arrowstyle": "-", "lw": 0.5, "color": color, "alpha": 0.55},
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.75},
        )


def write_qc_tables(result: WavelengthQcResult, output_dir: str | Path) -> None:
    """Write order summaries and per-line residuals as CSV files."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "order_fit_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "table",
                "order",
                "degree",
                "n_lines",
                "wavelength_min_nm",
                "wavelength_max_nm",
                "slope_mid_nm_per_px",
                "rms_residual_nm",
                "max_abs_residual_nm",
            ]
        )
        for table in result.tables:
            for fit in result.fits[table.name].values():
                writer.writerow(
                    [
                        fit.table,
                        fit.order,
                        fit.degree,
                        fit.n_lines,
                        f"{fit.wavelength_min_nm:.6f}",
                        f"{fit.wavelength_max_nm:.6f}",
                        f"{fit.slope_mid_nm_per_px:.9f}",
                        f"{fit.rms_residual_nm:.9f}",
                        f"{fit.max_abs_residual_nm:.9f}",
                    ]
                )

    with (out / "line_fit_residuals.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "table",
                "order",
                "center_px",
                "wavelength_nm",
                "fit_wavelength_nm",
                "residual_nm",
                "species",
                "comment",
            ]
        )
        for table in result.tables:
            fits = result.fits[table.name]
            for line in table.lines:
                fit = fits[_line_order(line)]
                fit_wavelength = float(fit.poly(_line_center(line)))
                writer.writerow(
                    [
                        table.name,
                        _line_order(line),
                        f"{_line_center(line):.6f}",
                        f"{_line_wavelength(line):.6f}",
                        f"{fit_wavelength:.6f}",
                        f"{_line_wavelength(line) - fit_wavelength:.9f}",
                        _line_species(line),
                        line.comment,
                    ]
                )


def plot_all_orders(result: WavelengthQcResult, output_dir: str | Path) -> None:
    """Plot all order fits and their calibration points."""
    out = Path(output_dir)
    x = np.arange(result.detector_width_px, dtype=float)
    for table in result.tables:
        fits = result.fits[table.name]
        orders = sorted(fits)
        colors = plt.cm.viridis(np.linspace(0, 1, len(orders)))
        fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
        for color, order in zip(colors, orders):
            fit = fits[order]
            order_lines = _lines_for_order(table.lines, order)
            ax.plot(x, fit.poly(x), color=color, lw=1.2)
            ax.scatter(
                [_line_center(line) for line in order_lines],
                [_line_wavelength(line) for line in order_lines],
                color=color,
                s=16,
            )
            ax.text(-30, fit.poly(0), str(order), ha="right", va="center", fontsize=8)
        ax.set_title(f"Wavelength calibration QC: {table.name}")
        ax.set_xlabel("pixel")
        ax.set_ylabel("wavelength, nm")
        ax.grid(True, alpha=0.2)
        fig.savefig(out / f"{table.name}_all_orders.png", dpi=180)
        plt.close(fig)


def plot_focus_orders(
    result: WavelengthQcResult,
    output_dir: str | Path,
    *,
    min_nm: float = 600.0,
    max_nm: float = 640.0,
    label_residual_threshold_nm: float = 0.02,
    max_labels_per_order: int = 4,
) -> None:
    """Plot fitted orders and calibration points around a wavelength range."""
    out = Path(output_dir)
    x = np.arange(result.detector_width_px, dtype=float)
    fig, axes = plt.subplots(
        1,
        len(result.tables),
        figsize=(14, 5),
        sharey=True,
        constrained_layout=True,
    )
    if len(result.tables) == 1:
        axes = [axes]
    for ax, table in zip(axes, result.tables):
        fits = result.fits[table.name]
        focus_orders = find_focus_orders(fits, min_nm=min_nm, max_nm=max_nm)
        colors = plt.cm.plasma(np.linspace(0, 1, len(focus_orders)))
        for color, order in zip(colors, focus_orders):
            fit = fits[order]
            y = fit.poly(x)
            keep = (y >= min_nm - 5.0) & (y <= max_nm + 5.0)
            order_lines = _lines_for_order(
                table.lines,
                order,
                min_nm=min_nm - 5.0,
                max_nm=max_nm + 5.0,
            )
            ax.plot(x[keep], y[keep], color=color, lw=1.5, label=f"order {order}")
            ax.scatter(
                [_line_center(line) for line in order_lines],
                [_line_wavelength(line) for line in order_lines],
                color=color,
                s=26,
                zorder=3,
            )
            label_points = sorted(
                [
                    (
                        _line_center(line),
                        _line_wavelength(line),
                        _short_label(line),
                        f"{abs(_residual_nm(line, fit)):.9f}",
                    )
                    for line in order_lines
                    if abs(_residual_nm(line, fit)) >= label_residual_threshold_nm
                ],
                key=lambda item: float(item[3]),
                reverse=True,
            )
            _annotate_decluttered(
                ax,
                label_points,
                color=color,
                max_labels=max_labels_per_order,
            )
        ax.axhspan(min_nm, max_nm, color="0.8", alpha=0.18)
        ax.set_title(table.name)
        ax.set_xlabel("pixel")
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("wavelength, nm")
    fig.suptitle(f"Wavelength calibration focus: {min_nm:g}-{max_nm:g} nm")
    fig.savefig(out / f"focus_{min_nm:g}_{max_nm:g}_orders.png", dpi=180)
    plt.close(fig)


def plot_order_residuals(result: WavelengthQcResult, output_dir: str | Path) -> None:
    """Plot per-order residual RMS and dispersion smoothness."""
    out = Path(output_dir)
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=False, constrained_layout=True)
    markers = ("o", "s", "^", "D", "v", "P", "X")
    linestyles = ("-", "--", "-.", ":")
    all_orders = sorted(
        {fit.order for fits_by_order in result.fits.values() for fit in fits_by_order.values()}
    )
    for index, table in enumerate(result.tables):
        fits = result.fits[table.name]
        orders = [fit.order for fit in fits.values()]
        axes[0].plot(
            orders,
            [fit.rms_residual_nm for fit in fits.values()],
            marker=markers[index % len(markers)],
            linestyle=linestyles[index % len(linestyles)],
            label=table.name,
        )
        axes[1].plot(
            orders,
            [fit.slope_mid_nm_per_px for fit in fits.values()],
            marker=markers[index % len(markers)],
            linestyle=linestyles[index % len(linestyles)],
            label=table.name,
        )
    if all_orders:
        for ax in axes:
            ax.set_xlim(min(all_orders) - 0.5, max(all_orders) + 0.5)
            ax.set_xticks(all_orders)
            ax.tick_params(axis="x", labelsize=8)
    axes[0].set_title("Per-order polynomial residual RMS")
    axes[0].set_xlabel("order")
    axes[0].set_ylabel("RMS residual, nm")
    axes[0].grid(True, axis="y", alpha=0.25)
    axes[0].grid(True, axis="x", alpha=0.12)
    axes[0].legend()
    axes[1].set_title("Mid-order dispersion smoothness")
    axes[1].set_xlabel("order")
    axes[1].set_ylabel("d wavelength / d pixel, nm/px")
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].grid(True, axis="x", alpha=0.12)
    axes[1].legend()
    fig.savefig(out / "order_residual_and_dispersion_qc.png", dpi=180)
    plt.close(fig)


def plot_focus_line_residuals(
    result: WavelengthQcResult,
    output_dir: str | Path,
    *,
    min_nm: float = 600.0,
    max_nm: float = 640.0,
    label_residual_threshold_nm: float = 0.02,
    max_labels_per_order: int = 4,
) -> None:
    """Plot residuals vs pixel for calibration lines in focused orders."""
    out = Path(output_dir)
    fig, axes = plt.subplots(
        1,
        len(result.tables),
        figsize=(14, 5),
        sharey=True,
        constrained_layout=True,
    )
    if len(result.tables) == 1:
        axes = [axes]
    for ax, table in zip(axes, result.tables):
        fits = result.fits[table.name]
        focus_orders = find_focus_orders(fits, min_nm=min_nm, max_nm=max_nm)
        colors = plt.cm.plasma(np.linspace(0, 1, len(focus_orders)))
        ax.axhline(0.0, color="0.7", lw=0.8)
        for color, order in zip(colors, focus_orders):
            fit = fits[order]
            order_lines = _lines_for_order(
                table.lines,
                order,
                min_nm=min_nm - 5.0,
                max_nm=max_nm + 5.0,
            )
            centers = np.asarray([_line_center(line) for line in order_lines], dtype=float)
            residuals = np.asarray([_residual_nm(line, fit) for line in order_lines], dtype=float)
            ax.scatter(centers, residuals, color=color, s=26, label=f"order {order}")
            if len(centers) >= 2:
                ax.plot(centers, residuals, color=color, lw=0.7, alpha=0.35)
            label_points = sorted(
                [
                    (
                        _line_center(line),
                        _residual_nm(line, fit),
                        _short_label(line),
                        _short_label(line),
                    )
                    for line in order_lines
                    if abs(_residual_nm(line, fit)) >= label_residual_threshold_nm
                ],
                key=lambda item: abs(item[1]),
                reverse=True,
            )
            _annotate_decluttered(
                ax,
                label_points,
                color=color,
                max_labels=max_labels_per_order,
            )
        ax.set_title(table.name)
        ax.set_xlabel("pixel")
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("line wavelength - fitted wavelength, nm")
    fig.suptitle(f"Calibration line residuals in {min_nm:g}-{max_nm:g} nm focus orders")
    fig.savefig(out / f"focus_{min_nm:g}_{max_nm:g}_line_residuals.png", dpi=180)
    plt.close(fig)


def run_wavelength_calibration_qc(
    tables: Sequence[tuple[str, str | Path]],
    output_dir: str | Path,
    *,
    detector_width_px: int = 2560,
    focus_min_nm: float = 600.0,
    focus_max_nm: float = 640.0,
    label_residual_threshold_nm: float = 0.02,
    max_labels_per_order: int = 4,
) -> WavelengthQcResult:
    """Run the standard wavelength-table QC and write plots/tables."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    loaded_tables = tuple(
        WavelengthTable(
            name=name,
            path=Path(path),
            lines=tuple(load_wavelength_table(path)),
        )
        for name, path in tables
    )
    fits = {
        table.name: fit_wavelength_orders(
            table.name,
            table.lines,
            detector_width_px=detector_width_px,
        )
        for table in loaded_tables
    }
    result = WavelengthQcResult(
        tables=loaded_tables,
        fits=fits,
        detector_width_px=detector_width_px,
    )
    write_qc_tables(result, output_path)
    plot_all_orders(result, output_path)
    plot_focus_orders(
        result,
        output_path,
        min_nm=focus_min_nm,
        max_nm=focus_max_nm,
        label_residual_threshold_nm=label_residual_threshold_nm,
        max_labels_per_order=max_labels_per_order,
    )
    plot_order_residuals(result, output_path)
    plot_focus_line_residuals(
        result,
        output_path,
        min_nm=focus_min_nm,
        max_nm=focus_max_nm,
        label_residual_threshold_nm=label_residual_threshold_nm,
        max_labels_per_order=max_labels_per_order,
    )
    return result
