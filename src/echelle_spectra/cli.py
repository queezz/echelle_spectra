"""The campaign-oriented ``echelle`` umbrella command."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

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
    commands.add_parser("snapshot", help="Create, validate, or inspect a calibration snapshot.")
    commands.add_parser("process", help="Export SIF files through the existing batch processor.")
    return parser


def _status(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="echelle status")
    parser.add_argument("--calibrations", default="calibrations", metavar="DIR")
    parser.add_argument("--registry", default="calibration_registry.toml", metavar="FILE")
    args = parser.parse_args(argv)
    root = Path(args.calibrations)
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
    registry = Path(args.registry)
    if registry.is_file():
        print(f"  registry:  {registry}")
    else:
        print(f"  registry:  not found ({registry})")
    print("  runs:      receipts are not implemented yet")
    return 1 if invalid else 0


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
    command, remainder = arguments[0], arguments[1:]
    if command == "status":
        return _status(remainder)
    if command == "snapshot":
        from .snapshot_cli import main as snapshot_main

        return snapshot_main(remainder, prog="echelle snapshot")
    if command == "process":
        from .spectrocube_cli import main as process_main

        result = process_main(remainder)
        return 0 if result is None else int(result)
    parser.error(f"unknown command: {command}")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
