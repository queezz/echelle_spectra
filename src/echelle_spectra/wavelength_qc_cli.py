"""CLI for wavelength-calibration lookup-table QC plots."""

from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

from . import _config
from .tools.wavelength_calibration_qc import run_wavelength_calibration_qc


def _default_calibration_dir() -> Path:
    return _config["base_path"] / "resources" / "calibration_files"


def _default_output_dir() -> Path:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    return Path("local") / "wavelength_calibration_qc" / f"{timestamp}-wavelength-qc"


def _parse_table(value: str) -> tuple[str, Path]:
    if "=" not in value:
        path = Path(value)
        return path.stem, path
    name, path = value.split("=", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError("table name cannot be empty")
    return name, Path(path.strip())


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-wavelength-qc",
        description=(
            "Reproduce Echelle wavelength-calibration polynomial fits and write "
            "per-order QC plots/tables."
        ),
    )
    parser.add_argument(
        "--calibration-dir",
        default=str(_default_calibration_dir()),
        help="Directory containing calibration lookup tables.",
    )
    parser.add_argument(
        "--table",
        action="append",
        type=_parse_table,
        metavar="NAME=PATH",
        help=(
            "Wavelength table to inspect. Relative paths are resolved from "
            "--calibration-dir. May be passed more than once."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=str(_default_output_dir()),
        help="Directory for PNG/CSV outputs.",
    )
    parser.add_argument(
        "--detector-width-px",
        type=int,
        default=2560,
        help="Detector width used to evaluate order polynomials.",
    )
    parser.add_argument(
        "--focus-min-nm",
        type=float,
        default=600.0,
        help="Lower wavelength for focused QC plots.",
    )
    parser.add_argument(
        "--focus-max-nm",
        type=float,
        default=640.0,
        help="Upper wavelength for focused QC plots.",
    )
    parser.add_argument(
        "--label-residual-threshold-nm",
        type=float,
        default=0.02,
        help="Only label focused calibration lines at or above this residual.",
    )
    parser.add_argument(
        "--max-labels-per-order",
        type=int,
        default=4,
        help="Maximum labels to draw per focused order.",
    )
    return parser


def _resolve_tables(calibration_dir: Path, tables: list[tuple[str, Path]] | None):
    if tables is None:
        tables = [
            ("base_20240305", Path("Th_wavelength_CMOS_20240305.txt")),
            (
                "aligned_20250926",
                Path("alignments") / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
            ),
        ]
    resolved = []
    for name, path in tables:
        full_path = path if path.is_absolute() else calibration_dir / path
        if not full_path.is_file():
            raise FileNotFoundError(f"wavelength table not found: {full_path}")
        resolved.append((name, full_path))
    return resolved


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    calibration_dir = Path(args.calibration_dir)
    try:
        tables = _resolve_tables(calibration_dir, args.table)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    result = run_wavelength_calibration_qc(
        tables,
        args.output_dir,
        detector_width_px=args.detector_width_px,
        focus_min_nm=args.focus_min_nm,
        focus_max_nm=args.focus_max_nm,
        label_residual_threshold_nm=args.label_residual_threshold_nm,
        max_labels_per_order=args.max_labels_per_order,
    )

    print(f"Wrote wavelength calibration QC to {Path(args.output_dir).resolve()}")
    for table in result.tables:
        worst = max(
            result.fits[table.name].values(),
            key=lambda fit: fit.rms_residual_nm,
        )
        print(
            f"  {table.name}: worst order {worst.order} "
            f"rms={worst.rms_residual_nm:.6f} nm "
            f"max={worst.max_abs_residual_nm:.6f} nm"
        )
    sys.exit(0)


if __name__ == "__main__":
    main()
