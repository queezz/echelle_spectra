"""Command adapters for the Packets 9--13 campaign domain modules."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def txt_main(argv: list[str] | None = None, *, prog: str = "echelle-cube2txt") -> int:
    from .lhd_text import write_cube_text

    parser = argparse.ArgumentParser(prog=prog, description="Write canonical LHD text from a saved cube.")
    parser.add_argument("input", help="Saved SpectroCube .nc file.")
    parser.add_argument("output", help="Destination text file.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    path = write_cube_text(args.input, args.output, overwrite=args.overwrite)
    print(path)
    return 0


def catalog_main(argv: list[str] | None = None, *, prog: str = "echelle catalog") -> int:
    from .catalog import build_drive_catalog, merge_catalogs

    parser = argparse.ArgumentParser(prog=prog, description="Build and merge portable cube catalogs.")
    commands = parser.add_subparsers(dest="action", metavar="COMMAND")
    build = commands.add_parser("build", help="Write one catalog beside a drive's cubes.")
    build.add_argument("cubes")
    build.add_argument("--volume-label", required=True)
    build.add_argument("--receipt-dir")
    build.add_argument("--output")
    merge = commands.add_parser("merge", help="Write an all-years index.")
    merge.add_argument("catalogs", nargs="+")
    merge.add_argument("-o", "--output", required=True)
    args = parser.parse_args(argv)
    if args.action is None:
        parser.print_help()
        return 0
    if args.action == "build":
        path = build_drive_catalog(
            args.cubes,
            volume_label=args.volume_label,
            receipt_dir=args.receipt_dir,
            output=args.output,
        )
    else:
        path = merge_catalogs(args.catalogs, args.output)
    print(path)
    return 0


def recal_main(argv: list[str] | None = None, *, prog: str = "echelle recal-cube") -> int:
    from .recalibration import recalibrate_cube

    parser = argparse.ArgumentParser(
        prog=prog,
        description="Revise saved-cube wavelength/flux calibration without reopening raw SIF data.",
    )
    parser.add_argument("input")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--new-snapshot", required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--wavelength-only", action="store_true")
    mode.add_argument("--factor-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    output, manifest = recalibrate_cube(
        args.input,
        args.output,
        new_snapshot_path=args.new_snapshot,
        update_wavelength=not args.factor_only,
        update_factor=not args.wavelength_only,
        overwrite=args.overwrite,
    )
    print(f"cube: {output}\nmanifest: {manifest}")
    return 0


def drift_main(argv: list[str] | None = None, *, prog: str = "echelle drift") -> int:
    from .drift import audit_cubes, create_refinement_snapshot, write_drift_evidence

    parser = argparse.ArgumentParser(prog=prog, description="Audit science-line drift and accept refinements.")
    commands = parser.add_subparsers(dest="action", metavar="COMMAND")
    audit = commands.add_parser(
        "audit", help="Measure sampled cubes and write one immutable verdict file."
    )
    audit.add_argument("cubes", nargs="+")
    audit.add_argument("--every", type=int, default=1)
    audit.add_argument("--shot", action="append", default=[])
    audit.add_argument("-o", "--output", required=True)
    refine = commands.add_parser(
        "refine", help="Accept a shifted verdict and emit an immutable -rN snapshot."
    )
    refine.add_argument("evidence")
    refine.add_argument("--calibrations", required=True)
    refine.add_argument("--accept-shift", required=True, type=float)
    args = parser.parse_args(argv)
    if args.action is None:
        parser.print_help()
        return 0
    if args.action == "audit":
        payload = audit_cubes(args.cubes, every=args.every, shots=set(args.shot))
        path = write_drift_evidence(args.output, payload)
        print(f"{payload['verdict']}: {path}")
        if payload.get("repair_command"):
            print(payload["repair_command"].replace("DRIFT_EVIDENCE.json", str(path)))
        return 0 if payload["verdict"] in {"aligned", "shifted"} else 1
    snapshot, accepted = create_refinement_snapshot(
        args.evidence,
        calibrations_root=args.calibrations,
        accepted_shift_nm=args.accept_shift,
    )
    print(f"snapshot: {snapshot.snapshot_id}\naccepted verdict: {accepted}")
    return 0


def web_main(argv: list[str] | None = None, *, prog: str = "echelle web") -> int:
    from .reading_room import build_reading_room

    parser = argparse.ArgumentParser(prog=prog, description="Build the read-only campaign reading room.")
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--drift", action="append", default=[])
    parser.add_argument("--document", action="append", default=[])
    args = parser.parse_args(argv)
    documents = args.document
    if not documents:
        documents = [
            str(path)
            for path in (
                Path("docs/operator-cheat-sheet.md"),
                Path("docs/calibration-snapshots.md"),
                Path("docs/calibration-epoch-registry.md"),
                Path("docs/harbor-candidate.md"),
            )
            if path.is_file()
        ]
    path = build_reading_room(
        args.catalog,
        args.output,
        drift_paths=args.drift,
        document_paths=documents,
    )
    print(path)
    return 0


def historical_main(argv: list[str] | None = None, *, prog: str = "echelle historical") -> int:
    from .historical import bundled_historical_manifests

    parser = argparse.ArgumentParser(prog=prog, description="Validate thin historical calibration binders.")
    parser.parse_args(argv)
    payload = bundled_historical_manifests()
    print(json.dumps([item.snapshot_id for item in payload]))
    return 0
