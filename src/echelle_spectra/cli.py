"""The campaign-oriented ``echelle`` umbrella command."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .calibration_registry import CalibrationRegistryError, load_calibration_registry
from .campaign_run import GATE_UNRECORDED, latest_run_summaries, list_run_summaries
from .snapshot import SnapshotValidationError, load_snapshot


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle",
        description=(
            "Echelle calibration and processing campaign tools.\n\n"
            "Start with 'echelle status'. Create or inspect a calibration with\n"
            "'echelle snapshot'. Convert SIF files with 'echelle process'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # A campaign report starts with which build produced it, and the answer used
    # to require importing the package by hand from a Python prompt.
    parser.add_argument(
        "--version",
        action="version",
        version=f"echelle {__version__}",
        help="Print the installed echelle_spectra version and exit.",
    )
    commands = parser.add_subparsers(dest="command", metavar="COMMAND")
    status = commands.add_parser("status", help="Summarize calibration snapshots and registry.")
    status.add_argument(
        "--calibrations",
        default="calibrations",
        metavar="DIR",
        help="Snapshot root to inspect (default: calibrations).",
    )
    status.add_argument(
        "--registry",
        default="calibration_registry.toml",
        metavar="FILE",
        help="Epoch registry to report (default: calibration_registry.toml).",
    )
    status.add_argument(
        "--runs",
        default="local/runs",
        metavar="DIR",
        help="Run receipts to inspect (default: local/runs).",
    )
    commands.add_parser("snapshot", help="Create, validate, or inspect a calibration snapshot.")
    commands.add_parser("process", help="Export SIF files through the existing batch processor.")
    commands.add_parser("catalog", help="Build or merge cube catalogs.")
    commands.add_parser("txt", help="Write canonical LHD text from a saved cube.")
    commands.add_parser("recal-cube", help="Revise saved-cube wavelength/flux calibration.")
    commands.add_parser(
        "inventory", help="Say what is on one campaign drive: logbooks, days, calibrations, space."
    )
    commands.add_parser(
        "drift", help="Audit science-line drift, survey pattern geometry, or accept a refinement."
    )
    commands.add_parser(
        "web", help="Build the read-only campaign page: stepper, drives, calibration, reading room."
    )
    commands.add_parser(
        "serve",
        help="Serve the campaign page on loopback so it can browse folders and set a home.",
    )
    commands.add_parser("historical", help="Validate bundled historical manifests.")
    return parser


def _print_combined_run_status(runs_root: Path) -> None:
    targets = latest_run_summaries(runs_root)
    if len(targets) <= 1:
        return
    statuses = {status for target in targets for status in target["counts"]}
    combined_counts = {
        status: sum(int(target["counts"].get(status, 0)) for target in targets)
        for status in statuses
    }
    combined_accounted = sum(combined_counts.values())
    combined_expected = sum(int(target["expected_files"]) for target in targets)
    combined_details = ", ".join(
        f"{count} {status}" for status, count in sorted(combined_counts.items()) if count
    )
    print(f"  targets:   {len(targets)} independent source(s)")
    print(
        f"  combined:  {combined_accounted}/{combined_expected} "
        f"({combined_details or 'no terminal records'})"
    )
    for target in targets:
        target_counts = target["counts"]
        target_accounted = sum(target_counts.values())
        print(
            f"    {target['volume_label']}: {target_accounted}/"
            f"{target['expected_files']} [{target['state']}] {target['id']}"
        )


def _print_registry_status(registry: Path, snapshots_root: Path) -> bool:
    """Print registry epochs and return whether validation failed."""

    if not registry.is_file():
        print(f"  registry:  not found ({registry})")
        return False
    try:
        loaded_registry = load_calibration_registry(registry, snapshots_root=snapshots_root)
    except CalibrationRegistryError as exc:
        print(f"  registry:  INVALID ({exc})")
        return True
    print(f"  registry:  {registry} ({len(loaded_registry.epochs)} epoch(s))")
    for epoch in loaded_registry.epochs:
        shot = (
            f"shot {epoch.shot_from if epoch.shot_from is not None else '-inf'}.."
            f"{epoch.shot_to if epoch.shot_to is not None else '+inf'}"
            if epoch.needs_shot
            else ""
        )
        day = (
            f"date {epoch.date_from.isoformat() if epoch.date_from else '-inf'}.."
            f"{epoch.date_to.isoformat() if epoch.date_to else '+inf'}"
            if epoch.needs_date
            else ""
        )
        bounds = " and ".join(value for value in (shot, day) if value)
        print(f"    {epoch.position}. {epoch.snapshot_id}: {bounds} (inclusive)")
    return False


def _status(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="echelle status")
    parser.add_argument("--calibrations", default="calibrations", metavar="DIR")
    parser.add_argument("--registry", default="calibration_registry.toml", metavar="FILE")
    parser.add_argument("--runs", default="local/runs", metavar="DIR")
    args = parser.parse_args(argv)
    # Absolute, always: a status line that echoes a relative path back tells an
    # operator standing in the wrong folder nothing about where it looked.
    root = Path(args.calibrations).expanduser().resolve()
    valid = []
    invalid: list[tuple[Path, SnapshotValidationError]] = []
    if root.is_dir():
        for manifest in sorted(root.glob("*/snapshot.toml")):
            try:
                valid.append(load_snapshot(manifest.parent))
            except SnapshotValidationError as exc:
                invalid.append((manifest.parent, exc))

    print("Echelle campaign status")
    if valid:
        current = sorted(valid, key=lambda item: item.snapshot_id)[-1]
        print(f"  snapshots: {len(valid)} valid")
        print(f"  current:   {current.snapshot_id} ({', '.join(current.lamps)})")
    else:
        print(f"  snapshots: none found under {root}")
    if invalid:
        print(f"  invalid:   {len(invalid)}")
        for path, exc in invalid:
            print(f"    {path.name}: {exc.errors[0]}")
    registry = Path(args.registry).expanduser().resolve()
    registry_invalid = _print_registry_status(registry, root)
    runs_root = Path(args.runs).expanduser().resolve()
    runs = list_run_summaries(runs_root)
    if runs:
        latest = runs[0]
        counts = latest["counts"]
        accounted = sum(counts.values())
        details = (
            ", ".join(f"{count} {status}" for status, count in counts.items() if count)
            or "no terminal records"
        )
        print(f"  runs:      {len(runs)} receipt(s) under {runs_root}")
        print(f"  latest:    {latest['id']} [{latest['state']}]")
        print(f"  progress:  {accounted}/{latest['expected_files']} ({details})")
        print(f"  snapshot:  {latest['snapshot_id']}")
        gate = str(latest.get("gate") or GATE_UNRECORDED)
        if latest.get("drift_evidence"):
            gate = f"{gate} ({latest.get('drift_verdict')}; {latest['drift_evidence']})"
        elif latest.get("sample"):
            gate = f"{gate} (unverified epoch sample)"
        print(f"  gate:      {gate}")
        _print_combined_run_status(runs_root)
    else:
        print(f"  runs:      none found under {runs_root}")
    return 1 if invalid or registry_invalid else 0


def main(argv: list[str] | None = None) -> int:
    """Dispatch one campaign command while preserving legacy entry points."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()
    if not arguments:
        parser.print_help()
        return 0
    if arguments[0] in {"-h", "--help"}:
        parser.print_help()
        return 0
    # Dispatch is by hand rather than by parse_args, so the umbrella's own flags
    # are answered here; argparse only carries them into --help.
    if arguments[0] == "--version":
        print(f"echelle {__version__}")
        return 0
    command, remainder = arguments[0], arguments[1:]
    if command == "status":
        return _status(remainder)
    if command == "snapshot":
        from .snapshot_cli import main as snapshot_main

        return snapshot_main(remainder, prog="echelle snapshot")
    if command == "process":
        from .spectrocube_cli import _build_parser as build_process_parser
        from .spectrocube_cli import main as process_main

        if not remainder:
            build_process_parser(prog="echelle process").print_help()
            return 0

        result = process_main(remainder, prog="echelle process")
        return 0 if result is None else int(result)
    if command == "serve":
        from .campaign_server import serve_main

        return int(serve_main(remainder, prog="echelle serve"))
    if command in {"catalog", "txt", "recal-cube", "drift", "inventory", "web", "historical"}:
        from . import campaign_tools_cli

        entry = {
            "catalog": campaign_tools_cli.catalog_main,
            "txt": campaign_tools_cli.txt_main,
            "recal-cube": campaign_tools_cli.recal_main,
            "drift": campaign_tools_cli.drift_main,
            "inventory": campaign_tools_cli.inventory_main,
            "web": campaign_tools_cli.web_main,
            "historical": campaign_tools_cli.historical_main,
        }[command]
        return int(entry(remainder, prog=f"echelle {command}"))
    parser.error(f"unknown command: {command}")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
