"""Synthetic line-list overlays for lamp wavelength-calibration review."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import percentile_filter
from scipy.signal import find_peaks, peak_widths

from .calibration_alignment import fit_single_gaussian_centroid, load_wavelength_table
from .echelle import Calibrations, EchelleImage
from .nist_lamp_calibration import load_nist_asd_exports

__all__ = [
    "NistOverlayConfig",
    "NistOverlayResult",
    "load_nist_line_csvs",
    "run_nist_synthetic_overlay",
]


@dataclass(frozen=True)
class NistOverlayConfig:
    """Inputs for one synthetic line-list overlay run."""

    calibration_dir: Path
    signal_file: Path
    wavelength_file: str
    pattern_file: str
    sphere_file: str
    sphere_background_file: str
    integral_file: str
    nist_csvs: Mapping[str, Path]
    output_dir: Path
    orders: Sequence[int]
    min_wavelength_nm: float
    max_wavelength_nm: float
    background_file: Path | None = None
    candidate_table_out: Path | None = None
    peak_prominence_sigma: float = 5.0
    peak_prominence_fraction: float = 0.04
    min_nist_weight: float = 0.18
    match_tolerance_nm: float = 0.05
    dominance_threshold: float = 0.45
    centroid_max_shift_px: float = 4.0
    synthetic_sigma_px: float = 2.4


@dataclass(frozen=True)
class NistOverlayResult:
    """Summary of generated overlay artifacts."""

    output_dir: Path
    n_nist_lines: int
    n_measured_peaks: int
    n_candidate_anchors: int
    candidate_table: Path | None


def load_nist_line_csvs(
    nist_csvs: Mapping[str, Path],
    *,
    min_wavelength_nm: float,
    max_wavelength_nm: float,
) -> pd.DataFrame:
    """Load NIST ASD CSV exports into normalized wavelength/weight rows."""
    return load_nist_asd_exports(
        nist_csvs,
        min_wavelength_nm=min_wavelength_nm,
        max_wavelength_nm=max_wavelength_nm,
    )


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _load_order_spectra(config: NistOverlayConfig):
    files = {
        "orders": config.pattern_file,
        "wavelength": config.wavelength_file,
        "sphr": str(config.sphere_file),
        "bkgr": str(config.sphere_background_file),
        "integral": config.integral_file,
    }
    calibration = Calibrations(folder=str(config.calibration_dir), filenames=files)
    calibration.start()
    signal = EchelleImage(str(config.signal_file), clbr=calibration)
    signal.calculate_order_spectra()
    signal.correct_order_shapes()
    spectra = np.nanmean(np.asarray(signal.order_spectra, dtype=float), axis=0)
    if config.background_file is not None:
        background = EchelleImage(str(config.background_file), clbr=calibration)
        background.calculate_order_spectra()
        background.correct_order_shapes()
        spectra = spectra - np.nanmean(np.asarray(background.order_spectra, dtype=float), axis=0)
    return calibration, spectra


def _preprocess(wavelength: np.ndarray, intensity: np.ndarray):
    pixels = np.arange(wavelength.size, dtype=float)
    finite = np.isfinite(wavelength) & np.isfinite(intensity)
    wavelength = wavelength[finite]
    intensity = intensity[finite]
    pixels = pixels[finite]
    order = np.argsort(wavelength)
    wavelength = wavelength[order]
    intensity = intensity[order]
    pixels = pixels[order]
    baseline = percentile_filter(intensity, percentile=12, size=151, mode="nearest")
    corrected = intensity - baseline
    scale = max(float(np.nanpercentile(corrected, 99.5)), 1.0)
    return wavelength, intensity, pixels, baseline, corrected, corrected / scale


def _synthetic_spectrum(wavelength: np.ndarray, lines: pd.DataFrame, sigma_nm: float) -> np.ndarray:
    y = np.zeros_like(wavelength, dtype=float)
    for row in lines.itertuples():
        y += float(row.weight) * np.exp(
            -0.5 * ((wavelength - float(row.wavelength_nm)) / sigma_nm) ** 2
        )
    if np.nanmax(y) > 0:
        y /= float(np.nanmax(y))
    return y


def _detect_peaks(
    wavelength: np.ndarray,
    pixels: np.ndarray,
    corrected: np.ndarray,
    *,
    min_wavelength_nm: float,
    max_wavelength_nm: float,
    peak_prominence_sigma: float,
    peak_prominence_fraction: float,
) -> list[dict]:
    noise = 1.4826 * np.nanmedian(np.abs(corrected - np.nanmedian(corrected)))
    if not np.isfinite(noise) or noise <= 0:
        noise = np.nanstd(corrected)
    peaks, props = find_peaks(
        corrected,
        prominence=max(
            peak_prominence_sigma * float(noise),
            peak_prominence_fraction * float(np.nanmax(corrected)),
        ),
        distance=6,
    )
    widths = peak_widths(corrected, peaks, rel_height=0.5)[0] if len(peaks) else []
    rows = []
    for idx, (peak, width) in enumerate(zip(peaks, widths), start=1):
        if not (min_wavelength_nm <= wavelength[peak] <= max_wavelength_nm):
            continue
        rows.append(
            {
                "peak_id": f"p{idx:03d}",
                "peak_pixel": float(pixels[peak]),
                "peak_nm_current_table": float(wavelength[peak]),
                "height": float(corrected[peak]),
                "prominence": float(props["prominences"][idx - 1]),
                "width_px": float(width),
            }
        )
    return rows


def _local_dominance(
    lines: pd.DataFrame, wavelength_nm: float, sigma_nm: float = 0.035
) -> tuple[float, int]:
    nearby = lines[np.abs(lines["wavelength_nm"] - wavelength_nm) <= 0.12].copy()
    if nearby.empty:
        return 0.0, 0
    kernel = np.exp(-0.5 * ((nearby["wavelength_nm"] - wavelength_nm) / sigma_nm) ** 2)
    weighted = nearby["weight"].to_numpy(dtype=float) * kernel.to_numpy(dtype=float)
    total = float(np.sum(weighted))
    return (float(np.max(weighted)) / total if total > 0 else 0.0), int(len(nearby))


def _match_peaks(
    peaks: list[dict],
    lines: pd.DataFrame,
    raw_intensity: np.ndarray,
    config: NistOverlayConfig,
) -> list[dict]:
    strong = lines[lines["weight"] >= config.min_nist_weight].copy()
    rows = []
    for peak in peaks:
        distances = np.abs(strong["wavelength_nm"] - peak["peak_nm_current_table"])
        if distances.empty:
            continue
        line = strong.loc[distances.idxmin()]
        delta_nm = float(line["wavelength_nm"] - peak["peak_nm_current_table"])
        nearby = lines[np.abs(lines["wavelength_nm"] - float(line["wavelength_nm"])) <= 0.08]
        dominance, nearby_count = _local_dominance(lines, float(line["wavelength_nm"]))
        ok, center, sigma, _amp, _baseline, snr, reason = fit_single_gaussian_centroid(
            raw_intensity,
            peak["peak_pixel"],
            window_radius_px=10,
            min_snr=5.0,
            min_sigma_px=0.5,
            max_sigma_px=7.0,
        )
        rows.append(
            {
                **peak,
                "nearest_species": line["species"],
                "nearest_nist_nm": float(line["wavelength_nm"]),
                "nist_minus_table_nm": delta_nm,
                "nist_weight": float(line["weight"]),
                "nist_intensity": (
                    float(line["intensity"]) if np.isfinite(line["intensity"]) else ""
                ),
                "isolated_nist_line": len(nearby) == 1,
                "local_nist_dominance": dominance,
                "local_nist_lines_0p12nm": nearby_count,
                "within_tolerance": abs(delta_nm) <= config.match_tolerance_nm,
                "gaussian_ok": ok,
                "gaussian_center_px": center,
                "gaussian_center_minus_peak_px": center - peak["peak_pixel"],
                "gaussian_sigma_px": sigma,
                "gaussian_snr": snr,
                "gaussian_reason": reason,
            }
        )
    return rows


def _candidate_rows(matches: list[dict], config: NistOverlayConfig) -> list[dict]:
    return [
        row
        for row in matches
        if row["within_tolerance"]
        and (
            row["isolated_nist_line"]
            or float(row["local_nist_dominance"]) >= config.dominance_threshold
        )
        and row["gaussian_ok"]
        and abs(float(row["gaussian_center_minus_peak_px"])) <= config.centroid_max_shift_px
        and float(row["nist_weight"]) >= config.min_nist_weight
    ]


def _plot_order(
    path: Path,
    order: int,
    wavelength,
    intensity,
    baseline,
    norm,
    synthetic,
    lines,
    peaks,
    candidates,
):
    fig, (ax, ax2) = plt.subplots(
        2,
        1,
        figsize=(15, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.2]},
    )
    ax.plot(wavelength, intensity, color="black", lw=0.7, label=f"lamp order {order}")
    ax.plot(wavelength, baseline, color="0.55", lw=0.8, ls=":", label="baseline")
    colors = {
        "ArI": "#c54e2f",
        "ArII": "#8f2f1e",
        "HI": "#2f6db3",
        "HII": "#1f4d7a",
        "HeI": "#d28a1f",
        "HeII": "#946313",
        "HgI": "#4f7f52",
        "HgII": "#2f5a34",
        "NeI": "#c74386",
        "NeII": "#8e2d5f",
        "ThI": "#5a3da8",
        "ThII": "#2d246d",
    }
    for row in lines.itertuples():
        ax.axvline(
            float(row.wavelength_nm),
            color=colors.get(row.species, "0.5"),
            lw=0.8,
            alpha=0.16 + 0.55 * float(row.weight),
        )
    for row in candidates:
        x = float(row["peak_nm_current_table"])
        ax.annotate(
            f"{row['nearest_species']} {float(row['nearest_nist_nm']):.4f}",
            (x, np.interp(x, wavelength, intensity)),
            xytext=(3, 5),
            textcoords="offset points",
            fontsize=7,
            rotation=35,
            color="0.15",
        )
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylabel("counts")
    ax.set_title(f"Order {order} lamp with synthetic line-list overlay", loc="left")

    ax2.plot(wavelength, norm, color="black", lw=0.8, label="baseline-subtracted lamp")
    ax2.plot(wavelength, synthetic, color="#157a6e", lw=1.0, label="synthetic sticks")
    for peak in peaks:
        x = float(peak["peak_nm_current_table"])
        ax2.scatter(
            x, np.interp(x, wavelength, norm), s=16, color="#f0a51a", edgecolor="black", lw=0.3
        )
    ax2.set_xlabel("Wavelength on current table [nm]")
    ax2.set_ylabel("normalized")
    ax2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_candidate_table(config: NistOverlayConfig, candidates: list[dict]) -> Path | None:
    if config.candidate_table_out is None:
        return None
    base_text = (config.calibration_dir / config.wavelength_file).read_text(encoding="utf-8")
    out = Path(config.candidate_table_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(base_text.rstrip())
        handle.write(
            "\n# Candidate anchors from synthetic line-list overlay; generated diagnostic\n"
        )
        for row in candidates:
            center = float(row["gaussian_center_px"])
            half_width = max(6.0, 2.0 * float(row["gaussian_sigma_px"]))
            handle.write(
                f"{int(row['order']):d}\t"
                f"{center - half_width:08.3f}\t{center + half_width:08.3f}\t"
                f"{center:010.4f}\t{float(row['nearest_nist_nm']):010.5f}\t"
                f"{row['nearest_species']}  # candidate synthetic overlay; "
                f"table_delta={float(row['nist_minus_table_nm']):+.5f} nm\n"
            )
    return out


def run_nist_synthetic_overlay(config: NistOverlayConfig) -> NistOverlayResult:
    """Run synthetic line-list overlay plots and candidate-anchor export."""
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    nist_lines = load_nist_line_csvs(
        config.nist_csvs,
        min_wavelength_nm=config.min_wavelength_nm - 0.5,
        max_wavelength_nm=config.max_wavelength_nm + 0.5,
    )
    nist_lines.to_csv(out_dir / "parsed_line_list.csv", index=False)
    calibration, spectra = _load_order_spectra(config)

    all_peaks: list[dict] = []
    all_matches: list[dict] = []
    all_candidates: list[dict] = []
    for order in config.orders:
        raw_intensity = np.asarray(spectra[int(order)], dtype=float)
        raw_wavelength = np.asarray(calibration.order_wavel[int(order)], dtype=float)
        wavelength, intensity, pixels, baseline, corrected, norm = _preprocess(
            raw_wavelength, raw_intensity
        )
        sigma_nm = float(np.nanmedian(np.abs(np.diff(wavelength)))) * config.synthetic_sigma_px
        in_order_lines = nist_lines[
            (nist_lines["wavelength_nm"] >= float(np.nanmin(wavelength)) - 0.2)
            & (nist_lines["wavelength_nm"] <= float(np.nanmax(wavelength)) + 0.2)
        ].copy()
        synthetic = _synthetic_spectrum(wavelength, in_order_lines, sigma_nm)
        peaks = _detect_peaks(
            wavelength,
            pixels,
            corrected,
            min_wavelength_nm=config.min_wavelength_nm,
            max_wavelength_nm=config.max_wavelength_nm,
            peak_prominence_sigma=config.peak_prominence_sigma,
            peak_prominence_fraction=config.peak_prominence_fraction,
        )
        for peak in peaks:
            peak["order"] = int(order)
        matches = _match_peaks(peaks, in_order_lines, raw_intensity, config)
        candidates = _candidate_rows(matches, config)
        all_peaks.extend(peaks)
        all_matches.extend(matches)
        all_candidates.extend(candidates)
        _plot_order(
            out_dir / f"order{int(order):02d}_synthetic_overlay.png",
            int(order),
            wavelength,
            intensity,
            baseline,
            norm,
            synthetic,
            in_order_lines,
            peaks,
            candidates,
        )

    _write_csv(out_dir / "measured_lamp_peaks.csv", all_peaks)
    _write_csv(out_dir / "peak_to_line_matches.csv", all_matches)
    _write_csv(out_dir / "candidate_anchors.csv", all_candidates)
    candidate_table = _write_candidate_table(config, all_candidates)

    existing_rows = load_wavelength_table(config.calibration_dir / config.wavelength_file)
    readme = [
        "# Synthetic Line-List Overlay",
        "",
        f"- Lamp: `{config.signal_file}`",
        f"- Background: `{config.background_file or ''}`",
        f"- Pattern: `{config.pattern_file}`",
        f"- Wavelength table: `{config.wavelength_file}`",
        f"- Orders: `{', '.join(str(o) for o in config.orders)}`",
        f"- Window: `{config.min_wavelength_nm:g}-{config.max_wavelength_nm:g} nm`",
        "",
        "## Counts",
        "",
        f"- Parsed line-list rows: `{len(nist_lines)}`",
        f"- Measured peaks: `{len(all_peaks)}`",
        f"- Candidate anchors: `{len(all_candidates)}`",
        f"- Existing active table rows: `{len(existing_rows)}`",
    ]
    if candidate_table is not None:
        readme.append(f"- Candidate table: `{candidate_table}`")
    (out_dir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    return NistOverlayResult(
        output_dir=out_dir,
        n_nist_lines=len(nist_lines),
        n_measured_peaks=len(all_peaks),
        n_candidate_anchors=len(all_candidates),
        candidate_table=candidate_table,
    )
