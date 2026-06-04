"""CLI for validating calibrated Echelle wavelengths against known lines."""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

from . import _config
from .tools.echelle import Calibrations, EchelleImage, Spectrum
from .tools.line_validation import (
    balmer_air_targets,
    build_stitched_order_index,
    load_fulcher_h2_q_branch_targets,
    summarize_validation,
    validate_lines,
)


def _default_calibration_dir() -> Path:
    return _config["base_path"] / "resources" / "calibration_files"


def _default_fulcher_table() -> Path:
    return (
        _config["base_path"].parents[2]
        / "fulcheranalyzer"
        / "src"
        / "fulcher_analyzer"
        / "data_molecular"
        / "fulcher-α_band_wavelength.txt"
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-validate-lines",
        description="Validate a calibrated Echelle SIF file against Balmer and optional Fulcher lines.",
    )
    parser.add_argument("sif", help="Experimental .SIF file to validate.")
    parser.add_argument(
        "--calibration-dir",
        default=str(_default_calibration_dir()),
        help="Directory containing calibration resources.",
    )
    parser.add_argument(
        "--pattern",
        default="pattern_CMOS_20250926.txt",
        help="Order-pattern table in --calibration-dir.",
    )
    parser.add_argument(
        "--wavelength",
        default="alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        help="Wavelength lookup table in --calibration-dir.",
    )
    parser.add_argument("--sphere", default="sphere_cmos_20240305.sif")
    parser.add_argument("--sphere-background", default="sphere_cmos_20240305_bkg.sif")
    parser.add_argument("--integral", default="integrating_sphere.txt")
    parser.add_argument(
        "--line-set",
        choices=("balmer", "balmer-fulcher"),
        default="balmer-fulcher",
        help="Known-line set to validate.",
    )
    parser.add_argument(
        "--fulcher-table",
        default=str(_default_fulcher_table()),
        help="Fulcher H2 Q-branch wavelength table.",
    )
    parser.add_argument(
        "--wavelength-medium",
        choices=("air", "vacuum"),
        default="air",
        help="Convention used by the expected line table.",
    )
    parser.add_argument("--window-nm", type=float, default=0.30)
    parser.add_argument("--min-snr", type=float, default=5.0)
    parser.add_argument(
        "--max-abs-residual-nm",
        type=float,
        default=None,
        help="Optional display filter for compact Fulcher tables.",
    )
    return parser


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
        sys.exit(1)


def _print_table(results) -> None:
    print(
        "| line | expected nm | measured nm | residual nm | order | peak SNR | notes |"
    )
    print("| --- | ---: | ---: | ---: | ---: | ---: | --- |")
    for result in results:
        order = "" if result.order is None else str(result.order)
        print(
            f"| {result.line} | "
            f"{result.expected_nm:.5f} | "
            f"{result.measured_nm:.5f} | "
            f"{result.residual_nm:+.5f} | "
            f"{order} | "
            f"{result.peak_snr:.2f} | "
            f"{result.notes} |"
        )


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.wavelength_medium != "air":
        print(
            "ERROR: this CLI currently validates against air wavelength tables only. "
            "Pass --wavelength-medium air or convert the expected table first.",
            file=sys.stderr,
        )
        sys.exit(1)

    sif = Path(args.sif).resolve()
    calibration_dir = Path(args.calibration_dir).resolve()
    fulcher_table = Path(args.fulcher_table).resolve()
    for label, path in [
        ("SIF", sif),
        ("calibration dir", calibration_dir),
        ("pattern", calibration_dir / args.pattern),
        ("wavelength", calibration_dir / args.wavelength),
        ("sphere", calibration_dir / args.sphere),
        ("sphere background", calibration_dir / args.sphere_background),
        ("integral", calibration_dir / args.integral),
    ]:
        if label == "calibration dir":
            if not path.is_dir():
                print(f"ERROR: calibration dir not found: {path}", file=sys.stderr)
                sys.exit(1)
        else:
            _require_file(path, label)

    targets = balmer_air_targets()
    if args.line_set == "balmer-fulcher":
        _require_file(fulcher_table, "Fulcher table")
        targets.extend(load_fulcher_h2_q_branch_targets(fulcher_table))

    filenames = {
        "orders": args.pattern,
        "wavelength": args.wavelength,
        "sphr": args.sphere,
        "bkgr": args.sphere_background,
        "integral": args.integral,
    }
    calibration = Calibrations(folder=str(calibration_dir), filenames=filenames)
    calibration.start()
    image = EchelleImage(str(sif), clbr=calibration)
    image.calculate_order_spectra()
    image.correct_order_shapes()
    image.calculate_spectra()
    spectrum = Spectrum(image)

    background_frames = set(spectrum.info.get("BackgroundFrames", []))
    signal_frames = [
        frame
        for frame in range(spectrum.counts.shape[0])
        if frame not in background_frames
    ]
    order_index = build_stitched_order_index(calibration)
    results = validate_lines(
        spectrum.wavelength,
        spectrum.counts,
        targets,
        signal_frames=signal_frames,
        order_index=order_index,
        window_nm=args.window_nm,
        min_snr=args.min_snr,
    )
    if args.max_abs_residual_nm is not None:
        results = [
            result
            for result in results
            if abs(result.residual_nm) <= args.max_abs_residual_nm
        ]

    print("Validation context:")
    print(f"  sif                = {sif}")
    print(f"  pattern            = {args.pattern}")
    print(f"  wavelength         = {args.wavelength}")
    print("  wavelength medium  = air")
    print(f"  frames             = {spectrum.counts.shape[0]}")
    print(f"  background frames  = {spectrum.info.get('BackgroundFrames', [])}")
    print(f"  signal frames      = {signal_frames[0]}..{signal_frames[-1]} ({len(signal_frames)})")
    print()
    _print_table(results)
    summary = summarize_validation(results)
    print()
    print("Residual summary:")
    for key, value in asdict(summary).items():
        if isinstance(value, float):
            print(f"  {key} = {value:.6f}")
        else:
            print(f"  {key} = {value}")
    sys.exit(0)


if __name__ == "__main__":
    main()
