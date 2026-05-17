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
        metavar="INPUT",
        help="Input .sif file or folder containing .sif files (batch mode).",
    )
    p.add_argument(
        "--units",
        choices=["counts", "wm", "wmsr", "phmsr"],
        default="counts",
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
        default="CMOS",
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
        "--instrument-id",
        default="echelle",
        metavar="ID",
        help="Instrument identifier stored in SpectroCube metadata (default: echelle).",
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


def _export_one(
    sif_path: Path,
    output_nc: Path,
    *,
    units: str,
    camera: str,
    calibration_dir: Path | None,
    instrument_id: str,
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
        )
        export_spectrocube(
            sp,
            str(output_nc),
            units=units,
            instrument_id=instrument_id,
            squeeze_single_frame=False,
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

    input_path = Path(args.input)
    cal_dir = Path(args.calibration_dir) if args.calibration_dir else None

    # ---- single file -------------------------------------------------------
    if input_path.is_file():
        if args.output:
            out = Path(args.output)
            if out.is_dir():
                out = _output_path_for(input_path, out)
        else:
            out = _output_path_for(input_path, input_path.parent)

        ok = _export_one(
            input_path,
            out,
            units=args.units,
            camera=args.camera,
            calibration_dir=cal_dir,
            instrument_id=args.instrument_id,
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
                    print(f"  Loading {args.camera} calibration …")
                clbr = build_calibration(cal_dir, args.camera)
                if args.verbose:
                    print("  Calibration ready.")
            except Exception as exc:
                print(f"ERROR: Could not load calibration: {exc}", file=sys.stderr)
                sys.exit(1)

        failed: list[Path] = []
        for sif in sif_files:
            nc_out = _output_path_for(sif, out_dir)
            ok = _export_one(
                sif,
                nc_out,
                units=args.units,
                camera=args.camera,
                calibration_dir=cal_dir,
                instrument_id=args.instrument_id,
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
