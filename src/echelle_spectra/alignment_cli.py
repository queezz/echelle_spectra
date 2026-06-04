"""CLI for aligning an Echelle wavelength table to a lamp frame."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from . import _config
from .tools.calibration_alignment import (
    AlignmentRunConfig,
    run_calibration_alignment,
    save_alignment_settings,
    write_wavelength_table,
)


def _default_calibration_dir() -> Path:
    return _config["base_path"] / "resources" / "calibration_files"


def _parse_csv_strings(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-align",
        description=(
            "Align an existing Echelle wavelength table to a lamp frame. "
            "By default this previews diagnostics only; pass --save to write artifacts."
        ),
    )
    parser.add_argument("signal", help="Lamp signal image, usually a Neon .sif file.")
    parser.add_argument("background", help="Matching lamp background image.")
    parser.add_argument("sphere", help="Sphere/flat-field image used for calibration loading.")
    parser.add_argument("sphere_background", help="Matching sphere background image.")
    parser.add_argument(
        "--calibration-dir",
        default=str(_default_calibration_dir()),
        help="Directory containing wavelength, pattern, and integrating-sphere tables.",
    )
    parser.add_argument(
        "--wavelength",
        default="Th_wavelength_CMOS_20240305.txt",
        help="Base wavelength lookup table in --calibration-dir.",
    )
    parser.add_argument(
        "--pattern",
        default="pattern_CMOS_20250926.txt",
        help="Base order-pattern table in --calibration-dir.",
    )
    parser.add_argument(
        "--integral",
        default="integrating_sphere.txt",
        help="Integrating-sphere calibration table in --calibration-dir.",
    )
    parser.add_argument("--dataset-id", default="20250926", help="Alignment dataset id.")
    parser.add_argument(
        "--source-dir",
        default="local/20250926_calib",
        help="Human-readable source directory recorded in settings metadata.",
    )
    parser.add_argument("--instrument-id", default="lhd_cmos", help="Instrument id.")
    parser.add_argument("--lamp", default="Ne", help="Lamp label recorded in settings metadata.")
    parser.add_argument("--created-at", default="2026-06-04", help="Metadata date string.")
    parser.add_argument(
        "--output-name",
        default="Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        help="Adjusted wavelength table filename recorded in settings metadata.",
    )
    parser.add_argument(
        "--settings-out",
        default=None,
        help=(
            "Settings TOML output path. Defaults to "
            "alignments/lhd_cmos_alignment_<dataset>.settings.toml."
        ),
    )
    parser.add_argument(
        "--table-out",
        default=None,
        help="Adjusted wavelength table output path. Defaults to alignments/<output-name>.",
    )
    parser.add_argument("--species", default="NeI", help="Comma-separated species to fit.")
    parser.add_argument("--min-snr", type=float, default=5.0, help="Minimum fitted SNR.")
    parser.add_argument(
        "--window-radius",
        type=int,
        default=18,
        help="1D centroid fit window radius in pixels.",
    )
    parser.add_argument(
        "--saturation-level",
        type=float,
        default=0.98 * 65535,
        help="Raw detector saturation guard level.",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Write settings and adjusted wavelength table after printing diagnostics.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow --save to replace existing output files.",
    )
    return parser


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
        sys.exit(1)


def _print_summary(result) -> None:
    transform = result.settings.transform
    residual = result.residual_px
    centroid_dx = result.centroid_dx_px
    print("Alignment summary:")
    print(f"  candidates = {len(result.candidates)}")
    print(f"  fits       = {len(result.fits)}")
    print(f"  good fits  = {len(result.good_fits)}")
    print(f"  dx         = {transform.dx_px:.4f} px")
    print(f"  dy         = {transform.dy_px:.4f} px")
    print(f"  theta      = {transform.theta_deg:.5f} deg")
    print(f"  RMS        = {result.settings.rms_px:.4f} px")
    print("Residual px:")
    print(f"  mean       = {float(np.mean(residual)):.4f}")
    print(f"  median     = {float(np.median(residual)):.4f}")
    print(f"  max        = {float(np.max(residual)):.4f}")
    print("Centroid dx px:")
    print(f"  mean       = {float(np.mean(centroid_dx)):.4f}")
    print(f"  median     = {float(np.median(centroid_dx)):.4f}")
    print(f"  min/max    = {float(np.min(centroid_dx)):.4f} / {float(np.max(centroid_dx)):.4f}")

    order = np.argsort(result.residual_px)[::-1][:8]
    print("Worst residuals:")
    for idx in order:
        fit = result.good_fits[int(idx)]
        print(
            "  "
            f"order={fit.line.order_idx:02d} "
            f"lambda={fit.line.wavelength_nm:.5f} nm "
            f"residual={result.residual_px[int(idx)]:.4f} px "
            f"dx={fit.center_pixel - fit.line.center_pixel:.4f} px"
        )


def _metadata(result, settings_out: Path) -> list[tuple[str, str]]:
    settings = result.settings
    return [
        ("Generated", settings.created_at),
        ("Base wavelength file", settings.base_wavelength_file),
        ("Base pattern file", settings.base_pattern_file),
        ("Alignment dataset", settings.alignment_dataset_id),
        ("Alignment source dir", settings.alignment_source_dir),
        ("Signal", settings.signal_file),
        ("Background", settings.background_file),
        ("Sphere", settings.sphere_file),
        ("Sphere background", settings.sphere_background_file),
        ("Correction model", "rigid detector transform, dx/dy/theta"),
        ("Settings file", settings_out.name),
        ("Note", settings.notes),
    ]


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    calibration_dir = Path(args.calibration_dir).resolve()
    signal = Path(args.signal).resolve()
    background = Path(args.background).resolve()
    sphere = Path(args.sphere).resolve()
    sphere_background = Path(args.sphere_background).resolve()
    for label, path in [
        ("signal", signal),
        ("background", background),
        ("sphere", sphere),
        ("sphere background", sphere_background),
        ("wavelength", calibration_dir / args.wavelength),
        ("pattern", calibration_dir / args.pattern),
        ("integral", calibration_dir / args.integral),
    ]:
        _require_file(path, label)

    alignment_dir = calibration_dir / "alignments"
    settings_out = (
        Path(args.settings_out)
        if args.settings_out
        else alignment_dir / f"lhd_cmos_alignment_{args.dataset_id}.settings.toml"
    )
    table_out = Path(args.table_out) if args.table_out else alignment_dir / args.output_name
    if args.save and not args.overwrite:
        existing = [path for path in (settings_out, table_out) if path.exists()]
        if existing:
            print(
                "ERROR: output exists; pass --overwrite: "
                + ", ".join(str(path) for path in existing),
                file=sys.stderr,
            )
            sys.exit(1)

    result = run_calibration_alignment(
        AlignmentRunConfig(
            calibration_dir=calibration_dir,
            signal_file=signal,
            background_file=background,
            sphere_file=sphere,
            sphere_background_file=sphere_background,
            base_wavelength_file=args.wavelength,
            base_pattern_file=args.pattern,
            integral_file=args.integral,
            instrument_id=args.instrument_id,
            alignment_dataset_id=args.dataset_id,
            alignment_source_dir=args.source_dir,
            alignment_lamp=args.lamp,
            created_at=args.created_at,
            output_wavelength_file=table_out.name,
            window_radius_px=args.window_radius,
            min_snr=args.min_snr,
            saturation_level=args.saturation_level,
            species=_parse_csv_strings(args.species),
        )
    )
    _print_summary(result)

    if args.save:
        settings_out.parent.mkdir(parents=True, exist_ok=True)
        table_out.parent.mkdir(parents=True, exist_ok=True)
        save_alignment_settings(result.settings, settings_out)
        write_wavelength_table(
            result.adjusted_rows,
            table_out,
            metadata=_metadata(result, settings_out),
        )
        print(f"Saved settings: {settings_out}")
        print(f"Saved adjusted table: {table_out}")
    else:
        print("Preview only; pass --save to write settings and adjusted wavelength table.")

    sys.exit(0)


if __name__ == "__main__":
    main()
