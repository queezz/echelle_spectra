"""CLI for exporting Echelle SIF files to the SpectroCube NetCDF format.

Entry point installed as ``echelle-spectrocube`` (see pyproject.toml).

Usage examples
--------------
Single file::

    echelle-spectrocube shot_042.sif --units wm -o shot_042.nc

Batch folder (all .SIF files)::

    echelle-spectrocube /data/shots/ --units wm -o /data/nc/

Dry run to preview what would happen::

    echelle-spectrocube /data/shots/ --dry-run --verbose

Limitation
----------
Unattended conversion requires the bundled calibration SIF files
(sphere images) in ``resources/calibration_files/``.  If these binary files
are absent from the installed package, use ``--calibration-dir`` to point
to a directory containing them.  The same calibration must be appropriate
for all SIF files in a batch.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .spectrocube_config import export_config_from_toml, export_plan_from_toml

_DEFAULTS = {
    "units": "counts",
    "camera": "CMOS",
    "instrument_id": "echelle",
    "wavelength_medium": "air",
    "drop_nonfinite_columns": True,
}

# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="echelle-spectrocube",
        description=(
            "Export Echelle .sif files to SpectroCube NetCDF (.nc) format.\n\n"
            "Single file:  echelle-spectrocube shot.sif -o shot.nc\n"
            "Batch folder: echelle-spectrocube /data/shots/ --units wm -o /out/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "input",
        nargs="?",
        metavar="INPUT",
        help="Input .sif file or folder containing .sif files (batch mode). Optional with --plan.",
    )
    p.add_argument(
        "--config",
        default=None,
        metavar="TOML",
        help="Calibration/export config TOML with stable camera/spectrometer settings.",
    )
    p.add_argument(
        "--plan",
        default=None,
        metavar="TOML",
        help="SpectroCube generation plan TOML. Can supply input/output and --config.",
    )
    p.add_argument(
        "--units",
        choices=["counts", "wm", "wmsr", "phmsr"],
        default=None,
        help=(
            "Intensity quantity to export.  One of: counts (default), wm "
            "[W m⁻² nm⁻¹], wmsr [W m⁻² sr⁻¹ nm⁻¹], phmsr [ph s⁻¹ m⁻² sr⁻¹ nm⁻¹]."
        ),
    )
    p.add_argument(
        "--frame",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Reserved: all frames are stored in the SpectroCube by default.  "
            "Kept for forward compatibility."
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="PATH",
        help=(
            "Output file (single-file mode) or output directory (batch mode).  "
            "Default: same directory as INPUT."
        ),
    )
    p.add_argument(
        "--camera",
        choices=["CMOS", "CCD"],
        default=None,
        help="Bundled calibration file set to use (default: CMOS).",
    )
    p.add_argument(
        "--calibration-dir",
        default=None,
        metavar="DIR",
        help=(
            "Path to calibration files directory.  "
            "Defaults to the bundled resources/calibration_files/ folder."
        ),
    )
    p.add_argument(
        "--order-pattern",
        default=None,
        metavar="FILE",
        help="Order-pattern table to use instead of the selected camera default.",
    )
    p.add_argument(
        "--wavelength",
        default=None,
        metavar="FILE",
        help="Wavelength lookup table to use instead of the selected camera default.",
    )
    p.add_argument(
        "--sphere",
        default=None,
        metavar="FILE",
        help="Integrating-sphere SIF to use instead of the selected camera default.",
    )
    p.add_argument(
        "--sphere-background",
        default=None,
        metavar="FILE",
        help="Integrating-sphere background SIF to use instead of the selected camera default.",
    )
    p.add_argument(
        "--integral",
        default=None,
        metavar="FILE",
        help="Integrating-sphere spectral table to use instead of the selected camera default.",
    )
    p.add_argument(
        "--instrument-id",
        default=None,
        metavar="ID",
        help="Instrument identifier stored in SpectroCube metadata (default: echelle).",
    )
    p.add_argument(
        "--wavelength-medium",
        choices=["air", "vacuum"],
        default=None,
        help="Wavelength convention stored in SpectroCube metadata (default: air).",
    )
    p.add_argument(
        "--wavelength-min-nm",
        type=float,
        default=None,
        help="Crop exported wavelengths below this inclusive lower bound.",
    )
    p.add_argument(
        "--wavelength-max-nm",
        type=float,
        default=None,
        help="Crop exported wavelengths above this inclusive upper bound.",
    )
    p.add_argument(
        "--calibration-source",
        default=None,
        help="Calibration source metadata for absolute intensity exports.",
    )
    p.add_argument(
        "--no-drop-nonfinite-columns",
        action="store_false",
        dest="drop_nonfinite_columns",
        default=None,
        help="Do not drop wavelength columns containing non-finite intensities.",
    )
    p.add_argument(
        "--output-suffix",
        default=None,
        help="Batch output suffix before .nc (default: _spectrocube).",
    )
    p.add_argument(
        "--pattern",
        default="*.SIF",
        metavar="GLOB",
        help="Glob pattern for batch SIF discovery (default: *.SIF).",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files (default: skip).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be converted without writing any files.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-file progress.",
    )
    return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _output_path_for(sif_path: Path, output_dir: Path) -> Path:
    """Derive the .nc output path for a given SIF file."""
    return output_dir / f"{sif_path.stem}_spectrocube.nc"


def _output_path_with_suffix(sif_path: Path, output_dir: Path, suffix: str | None) -> Path:
    suffix = suffix or "_spectrocube"
    return output_dir / f"{sif_path.stem}{suffix}.nc"


def _settings_from_args(args: argparse.Namespace) -> tuple[argparse.Namespace, dict]:
    """Merge built-in defaults, config TOML, plan TOML, and CLI overrides."""
    plan = export_plan_from_toml(args.plan) if args.plan else {}

    config_path = args.config or plan.get("config")
    settings = dict(_DEFAULTS)
    settings["calibration_files"] = {}
    settings["extra_attrs"] = {}
    if config_path:
        config_settings = export_config_from_toml(config_path)
        for key, value in config_settings.items():
            if value not in (None, {}, ""):
                if key == "extra_attrs":
                    settings["extra_attrs"].update(value)
                else:
                    settings[key] = value

    for key, attr in {
        "units": "units",
        "camera": "camera",
        "instrument_id": "instrument_id",
        "wavelength_medium": "wavelength_medium",
        "wavelength_min_nm": "wavelength_min_nm",
        "wavelength_max_nm": "wavelength_max_nm",
        "calibration_source": "calibration_source",
        "drop_nonfinite_columns": "drop_nonfinite_columns",
        "output_suffix": "output_suffix",
    }.items():
        value = getattr(args, attr)
        if value is not None:
            settings[key] = value

    calibration_files = dict(settings.get("calibration_files") or {})
    for key, value in {
        "orders": args.order_pattern,
        "wavelength": args.wavelength,
        "sphr": args.sphere,
        "bkgr": args.sphere_background,
        "integral": args.integral,
    }.items():
        if value:
            calibration_files[key] = value
    settings["calibration_files"] = calibration_files

    if args.calibration_dir is not None:
        settings["calibration_dir"] = args.calibration_dir

    if args.input is None:
        args.input = plan.get("input") or plan.get("input_dir")
    if args.output is None:
        args.output = plan.get("output") or plan.get("output_dir")
    if args.pattern == "*.SIF" and plan.get("pattern"):
        args.pattern = plan["pattern"]
    args.overwrite = bool(args.overwrite or plan.get("overwrite", False))
    args.dry_run = bool(args.dry_run or plan.get("dry_run", False))
    args.verbose = bool(args.verbose or plan.get("verbose", False))

    return args, settings


def _export_one(
    sif_path: Path,
    output_nc: Path,
    *,
    units: str,
    camera: str,
    calibration_dir: Path | None,
    calibration_files: dict[str, str] | None,
    instrument_id: str,
    wavelength_medium: str,
    wavelength_min_nm: float | None,
    wavelength_max_nm: float | None,
    calibration_source: str | None,
    drop_nonfinite_columns: bool,
    extra_attrs: dict,
    overwrite: bool,
    dry_run: bool,
    verbose: bool,
    calibration: object | None = None,
) -> bool:
    """Export one SIF file.  Returns True on success, False on failure."""
    from .tools.spectrocube_export import export_spectrocube

    if output_nc.exists() and not overwrite:
        if verbose:
            print(f"  SKIP  {sif_path.name}  (output exists; use --overwrite)")
        return True

    if dry_run:
        print(f"  DRY   {sif_path}  ->  {output_nc}")
        return True

    if verbose:
        print(f"  ...   {sif_path.name}")

    try:
        from .tools.loader import load_spectrum

        sp = load_spectrum(
            sif_path,
            calibration_folder=calibration_dir,
            camera=camera,
            calibration=calibration,
            calibration_files=calibration_files,
        )
        export_spectrocube(
            sp,
            str(output_nc),
            units=units,
            instrument_id=instrument_id,
            wavelength_medium=wavelength_medium,
            wavelength_min_nm=wavelength_min_nm,
            wavelength_max_nm=wavelength_max_nm,
            calibration_source=calibration_source,
            drop_nonfinite_columns=drop_nonfinite_columns,
            squeeze_single_frame=False,
            **extra_attrs,
        )
    except Exception as exc:
        print(f"  FAIL  {sif_path.name}: {exc}", file=sys.stderr)
        return False

    if verbose:
        print(f"  OK    {output_nc}")
    return True


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    """Main entry point for the ``echelle-spectrocube`` console script."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    args, settings = _settings_from_args(args)

    # Verify spectrocube is importable before doing any work
    try:
        import spectrocube  # noqa: F401
    except ImportError:
        print(
            "ERROR: The 'spectrocube' package is not installed.\n"
            "Install it with:\n"
            "    pip install spectrocube\n"
            "Or for local development:\n"
            "    pip install -e /path/to/2026-spectrocube",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.input is None:
        parser.error("INPUT is required unless supplied by --plan.")

    input_path = Path(args.input)
    cal_dir = Path(settings["calibration_dir"]) if settings.get("calibration_dir") else None
    calibration_files = settings["calibration_files"]

    # ---- single file -------------------------------------------------------
    if input_path.is_file():
        if args.output:
            out = Path(args.output)
            if out.is_dir():
                out = _output_path_with_suffix(input_path, out, settings.get("output_suffix"))
        else:
            out = _output_path_with_suffix(input_path, input_path.parent, settings.get("output_suffix"))

        ok = _export_one(
            input_path,
            out,
            units=settings["units"],
            camera=settings["camera"],
            calibration_dir=cal_dir,
            calibration_files=calibration_files,
            instrument_id=settings["instrument_id"],
            wavelength_medium=settings["wavelength_medium"],
            wavelength_min_nm=settings.get("wavelength_min_nm"),
            wavelength_max_nm=settings.get("wavelength_max_nm"),
            calibration_source=settings.get("calibration_source"),
            drop_nonfinite_columns=settings["drop_nonfinite_columns"],
            extra_attrs=settings["extra_attrs"],
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            verbose=True,
        )
        sys.exit(0 if ok else 1)

    # ---- batch folder ------------------------------------------------------
    elif input_path.is_dir():
        sif_files = sorted(input_path.glob(args.pattern))
        # If user kept the default *.SIF pattern and nothing was found, also try *.sif
        if not sif_files and args.pattern == "*.SIF":
            sif_files = sorted(input_path.glob("*.sif"))

        if not sif_files:
            print(
                f"No files matching '{args.pattern}' found in {input_path}",
                file=sys.stderr,
            )
            sys.exit(1)

        out_dir = Path(args.output) if args.output else input_path
        if not args.dry_run:
            out_dir.mkdir(parents=True, exist_ok=True)

        print(f"Batch: {len(sif_files)} file(s)  ->  {out_dir}")

        # Load calibration once and reuse across all files (avoids re-reading
        # the sphere SIF images for every input file).
        clbr = None
        if not args.dry_run:
            try:
                from .tools.loader import build_calibration

                if args.verbose:
                    print(f"  Loading {settings['camera']} calibration …")
                clbr = build_calibration(
                    cal_dir,
                    settings["camera"],
                    calibration_files=calibration_files,
                )
                if args.verbose:
                    print("  Calibration ready.")
            except Exception as exc:
                print(f"ERROR: Could not load calibration: {exc}", file=sys.stderr)
                sys.exit(1)

        failed: list[Path] = []
        for sif in sif_files:
            nc_out = _output_path_with_suffix(sif, out_dir, settings.get("output_suffix"))
            ok = _export_one(
                sif,
                nc_out,
                units=settings["units"],
                camera=settings["camera"],
                calibration_dir=cal_dir,
                calibration_files=calibration_files,
                instrument_id=settings["instrument_id"],
                wavelength_medium=settings["wavelength_medium"],
                wavelength_min_nm=settings.get("wavelength_min_nm"),
                wavelength_max_nm=settings.get("wavelength_max_nm"),
                calibration_source=settings.get("calibration_source"),
                drop_nonfinite_columns=settings["drop_nonfinite_columns"],
                extra_attrs=settings["extra_attrs"],
                overwrite=args.overwrite,
                dry_run=args.dry_run,
                verbose=args.verbose,
                calibration=clbr,
            )
            if not ok:
                failed.append(sif)

        n_ok = len(sif_files) - len(failed)
        if failed:
            print(f"\n{len(failed)} failure(s):", file=sys.stderr)
            for f in failed:
                print(f"  {f}", file=sys.stderr)
            print(f"{n_ok}/{len(sif_files)} exported successfully.", file=sys.stderr)
            sys.exit(1)

        if args.dry_run:
            print(f"Dry run complete. {len(sif_files)} file(s) would be converted.")
        else:
            print(f"Done. {n_ok}/{len(sif_files)} exported successfully.")
        sys.exit(0)

    # ---- path not found ----------------------------------------------------
    else:
        print(f"ERROR: Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
