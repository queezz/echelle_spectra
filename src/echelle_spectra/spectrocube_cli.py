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
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .campaign_run import (
    RunReceipt,
    default_volume_label,
    find_resumable_run,
    new_run_directory,
    utc_now,
)
from .spectrocube_config import export_config_from_toml, export_plan_from_toml

_DEFAULTS = {
    "units": "counts",
    "camera": "CMOS",
    "instrument_id": "echelle",
    "wavelength_medium": "air",
    "drop_nonfinite_columns": True,
}

_ExportStatus = Literal["exported", "skipped", "dry-run", "failed"]


@dataclass(frozen=True)
class ExportResult:
    status: _ExportStatus
    reason: str = ""


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def _build_parser(*, prog: str = "echelle-spectrocube") -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Export Echelle .sif files to SpectroCube NetCDF (.nc) format.\n\n"
            f"Single file:  {prog} shot.sif -o shot.nc\n"
            f"Batch folder: {prog} /data/shots/ --units wm -o /out/"
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
            "Intensity quantity to export. One of: counts (default), wm "
            "[W m^-2 nm^-1], wmsr [W m^-2 sr^-1 nm^-1], "
            "phmsr [ph s^-1 m^-2 sr^-1 nm^-1]."
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
    p.add_argument(
        "--runs-dir",
        default="local/runs",
        metavar="DIR",
        help="Root for durable batch receipts (default: local/runs).",
    )
    p.add_argument(
        "--run-dir",
        default=None,
        metavar="DIR",
        help="Use or resume one exact run receipt directory.",
    )
    p.add_argument(
        "--new-run",
        action="store_true",
        help="Start a new receipt even when a matching interrupted run exists.",
    )
    p.add_argument(
        "--snapshot-id",
        default=None,
        metavar="ID",
        help="Calibration snapshot ID recorded in every receipt (default: unassigned).",
    )
    p.add_argument(
        "--volume-label",
        default=None,
        metavar="LABEL",
        help="USB volume label recorded in every receipt (default: drive/root identity).",
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


def _color(text: str, code: str, *, stream=sys.stdout) -> str:
    if hasattr(stream, "isatty") and stream.isatty():
        return f"\033[{code}m{text}\033[0m"
    return text


def _batch_header(
    *,
    input_path: Path,
    output_dir: Path,
    n_files: int,
    pattern: str,
    dry_run: bool,
) -> None:
    mode = "DRY RUN" if dry_run else "export"
    print(_color("SpectroCube batch", "1;36"))
    print(f"Source:      {input_path}")
    print(f"Destination: {output_dir}")
    print(f"Pattern:     {pattern}")
    print(f"Files:       {n_files} ({mode})")


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
) -> ExportResult:
    """Export one SIF file."""
    from .tools.spectrocube_export import export_spectrocube

    if output_nc.exists() and not overwrite:
        if verbose:
            print(_color(f"SKIP {sif_path.name} (output exists; use --overwrite)", "33"))
        return ExportResult("skipped", "output exists; use --overwrite to replace it")

    if dry_run:
        if verbose:
            print(_color(f"DRY {sif_path.name}", "36"))
        return ExportResult("dry-run", "dry run")

    temporary_output = output_nc.with_name(f".{output_nc.name}.{os.getpid()}.tmp")
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
            str(temporary_output),
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
        os.replace(temporary_output, output_nc)
    except KeyboardInterrupt:
        temporary_output.unlink(missing_ok=True)
        raise
    except Exception as exc:
        temporary_output.unlink(missing_ok=True)
        print(_color(f"FAIL {sif_path.name}: {exc}", "31", stream=sys.stderr), file=sys.stderr)
        return ExportResult("failed", str(exc))

    return ExportResult("exported")


def _normalize_export_result(result: object) -> ExportResult:
    """Keep older test/integration mocks compatible with the richer result."""
    if isinstance(result, ExportResult):
        return result
    if result is True:
        return ExportResult("exported")
    if result in {"exported", "skipped", "dry-run", "failed"}:
        return ExportResult(result)  # type: ignore[arg-type]
    return ExportResult("failed", f"unexpected exporter result: {result!r}")


def _progress_line(index: int, total: int, started: float, status: str) -> str:
    elapsed = time.monotonic() - started
    if elapsed < 0.05:
        return f"[{index}/{total}] {status} | rate measuring | ETA --"
    rate = index / elapsed
    remaining = max(total - index, 0)
    eta_s = remaining / rate if rate else 0.0
    return f"[{index}/{total}] {status} | {rate:.2f} file/s | ETA {eta_s:.1f}s"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None, *, prog: str = "echelle-spectrocube") -> None:
    """Main entry point for the ``echelle-spectrocube`` console script."""
    parser = _build_parser(prog=prog)
    args = parser.parse_args(argv)
    args, settings = _settings_from_args(args)

    # Dry-run discovery and planning do not need the optional container package.
    if not args.dry_run:
        try:
            import spectrocube  # noqa: F401
        except ImportError:
            print(
                "ERROR: The 'spectrocube' package is not installed.\n"
                "Install it with:\n"
                "    pip install 'echelle_spectra[spectrocube]'\n"
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
            out = _output_path_with_suffix(
                input_path, input_path.parent, settings.get("output_suffix")
            )

        result = _normalize_export_result(
            _export_one(
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
        )
        sys.exit(0 if result.status != "failed" else 1)

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

        _batch_header(
            input_path=input_path,
            output_dir=out_dir,
            n_files=len(sif_files),
            pattern=args.pattern,
            dry_run=args.dry_run,
        )

        receipt = None
        if not args.dry_run:
            runs_root = Path(args.runs_dir)
            if args.run_dir:
                receipt_dir = Path(args.run_dir)
            elif not args.new_run:
                receipt_dir = find_resumable_run(runs_root, input_path, out_dir, args.pattern)
                if receipt_dir is None:
                    receipt_dir = new_run_directory(runs_root, input_path)
            else:
                receipt_dir = new_run_directory(runs_root, input_path)

            resuming = (receipt_dir / "run.toml").is_file()
            if resuming:
                receipt = RunReceipt.load(receipt_dir)
                if receipt.source_root.resolve() != input_path.resolve():
                    print("ERROR: Run receipt source does not match this input.", file=sys.stderr)
                    sys.exit(1)
                if receipt.output_root.resolve() != out_dir.resolve():
                    print(
                        "ERROR: Run receipt destination does not match this output.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if receipt.pattern != args.pattern:
                    print("ERROR: Run receipt pattern does not match --pattern.", file=sys.stderr)
                    sys.exit(1)
            else:
                receipt = RunReceipt.create(
                    receipt_dir,
                    source_root=input_path,
                    output_root=out_dir,
                    pattern=args.pattern,
                    volume_label=args.volume_label or default_volume_label(input_path),
                    snapshot_id=args.snapshot_id or "unassigned",
                    expected_files=len(sif_files),
                )
                print(f"Receipt:     {receipt.directory}")
            if args.snapshot_id and receipt.snapshot_id != args.snapshot_id:
                print("ERROR: Run receipt snapshot does not match --snapshot-id.", file=sys.stderr)
                sys.exit(1)
            if args.volume_label and receipt.volume_label != args.volume_label:
                print("ERROR: Run receipt volume does not match --volume-label.", file=sys.stderr)
                sys.exit(1)
            if resuming:
                receipt.state = "running"
                receipt.expected_files = max(receipt.expected_files, len(sif_files))
                receipt.write_manifest()
                print(f"Resuming:    {receipt.directory}")

        # Load calibration once and reuse across all files (avoids re-reading
        # the sphere SIF images for every input file).
        clbr = None
        if not args.dry_run:
            try:
                from .tools.loader import build_calibration

                if args.verbose:
                    print(f"Loading {settings['camera']} calibration...")
                clbr = build_calibration(
                    cal_dir,
                    settings["camera"],
                    calibration_files=calibration_files,
                )
                if args.verbose:
                    print(_color("Calibration ready.", "32"))
            except KeyboardInterrupt:
                if receipt is not None:
                    receipt.finish("interrupted")
                    print(f"\nInterrupted safely. Resume with --run-dir {receipt.directory}")
                sys.exit(130)
            except Exception as exc:
                print(f"ERROR: Could not load calibration: {exc}", file=sys.stderr)
                if receipt is not None:
                    for sif in sif_files:
                        source = receipt.identity_for(sif)
                        receipt.append(
                            source,
                            _output_path_with_suffix(sif, out_dir, settings.get("output_suffix")),
                            status="failed",
                            started_at=utc_now(),
                            finished_at=utc_now(),
                            duration_s=0.0,
                            reason=f"calibration could not load: {exc}",
                        )
                    receipt.finish("partial")
                sys.exit(1)

        failed: list[Path] = []
        n_exported = 0
        n_skipped = 0
        n_dry_run = 0
        total = len(sif_files)
        batch_started = time.monotonic()
        for index, sif in enumerate(sif_files, start=1):
            nc_out = _output_path_with_suffix(sif, out_dir, settings.get("output_suffix"))
            if args.verbose:
                print(_color(f"[{index}/{total}] {sif.name}", "36"))
            item_started = time.monotonic()
            started_at = utc_now()
            try:
                source = receipt.identity_for(sif) if receipt is not None else None
            except KeyboardInterrupt:
                if receipt is not None:
                    receipt.finish("interrupted")
                    print(f"\nInterrupted safely. Resume with --run-dir {receipt.directory}")
                sys.exit(130)
            if (
                receipt is not None
                and source is not None
                and not args.overwrite
                and receipt.completed_output_is_valid(source, nc_out)
            ):
                result = ExportResult("skipped", "completed output verified from prior receipt")
            elif (
                receipt is not None
                and source is not None
                and not args.overwrite
                and nc_out.exists()
                and receipt.has_export_record(source)
            ):
                result = ExportResult(
                    "failed",
                    "recorded completed output changed; inspect it or use --overwrite",
                )
            else:
                try:
                    result = _normalize_export_result(
                        _export_one(
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
                    )
                except KeyboardInterrupt:
                    if receipt is not None and source is not None:
                        receipt.append(
                            source,
                            nc_out,
                            status="interrupted",
                            started_at=started_at,
                            finished_at=utc_now(),
                            duration_s=time.monotonic() - item_started,
                            reason="keyboard interrupt",
                        )
                        receipt.finish("interrupted")
                        print(f"\nInterrupted safely. Resume with --run-dir {receipt.directory}")
                    sys.exit(130)

            if result.status == "exported" and not nc_out.is_file():
                result = ExportResult("failed", "exporter returned success without an output file")

            status = result.status
            if receipt is not None and source is not None:
                receipt.append(
                    source,
                    nc_out,
                    status=status,
                    started_at=started_at,
                    finished_at=utc_now(),
                    duration_s=time.monotonic() - item_started,
                    reason=result.reason,
                )
            if status == "failed":
                failed.append(sif)
            elif status == "skipped":
                n_skipped += 1
            elif status == "dry-run":
                n_dry_run += 1
            else:
                n_exported += 1
            print(_progress_line(index, total, batch_started, status))

        n_ok = len(sif_files) - len(failed)
        if failed:
            if receipt is not None:
                receipt.finish("partial")
            print(f"\n{len(failed)} failure(s):", file=sys.stderr)
            for f in failed:
                print(f"  {f}", file=sys.stderr)
            print(f"{n_ok}/{len(sif_files)} exported successfully.", file=sys.stderr)
            sys.exit(1)

        if args.dry_run:
            print(
                _color(f"Dry run complete. {n_dry_run}/{total} file(s) would be converted.", "36")
            )
        else:
            if receipt is not None:
                receipt.finish("completed")
            skipped = f", {n_skipped} skipped" if n_skipped else ""
            print(_color(f"Done. {n_exported}/{total} exported successfully{skipped}.", "32"))
        sys.exit(0)

    # ---- path not found ----------------------------------------------------
    else:
        print(f"ERROR: Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
