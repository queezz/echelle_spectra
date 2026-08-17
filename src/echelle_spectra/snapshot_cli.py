"""Command-line interface for calibration snapshots."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .snapshot import (
    Snapshot,
    SnapshotError,
    SnapshotValidationError,
    create_snapshot,
    load_snapshot,
)


def _build_parser(prog: str = "echelle-snapshot") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Create and verify immutable calibration snapshot folders. "
            "A snapshot.toml binder records every artifact role and SHA-256 digest."
        ),
    )
    commands = parser.add_subparsers(dest="snapshot_command", metavar="COMMAND")

    create = commands.add_parser(
        "create",
        help="Copy calibration inputs into a new snapshot and write snapshot.toml.",
    )
    create.add_argument("snapshot_id", help="Folder identity: YYYYMMDD_<detector>[-rev].")
    create.add_argument(
        "--output-root",
        default="calibrations",
        metavar="DIR",
        help="Parent directory for the snapshot (default: calibrations).",
    )
    create.add_argument("--detector", required=True, help="Detector token used in the id.")
    create.add_argument("--pattern", required=True, metavar="FILE")
    create.add_argument("--wavelength", required=True, metavar="FILE")
    create.add_argument("--sphere", required=True, metavar="FILE")
    create.add_argument("--sphere-background", required=True, metavar="FILE")
    create.add_argument("--integral", required=True, metavar="FILE")
    create.add_argument(
        "--lamp-used",
        action="append",
        required=True,
        metavar="NAME",
        help="Lamp species used in this calibration; repeat for several lamps.",
    )
    create.add_argument(
        "--lamp-file",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Source lamp SIF with its species label; repeat as needed.",
    )
    create.add_argument("--base-snapshot", metavar="ID")
    create.add_argument("--notes", default="", help="Free-text instrument notes.")
    create.add_argument("--shot-from", type=int)
    create.add_argument("--shot-to", type=int)
    create.add_argument("--date-from", metavar="YYYY-MM-DD")
    create.add_argument("--date-to", metavar="YYYY-MM-DD")
    create.add_argument("--alignment-dx-px", type=float)
    create.add_argument("--alignment-dy-px", type=float)
    create.add_argument("--alignment-rotation-deg", type=float)
    create.add_argument("--alignment-rms-px", type=float)
    create.add_argument("--qc-lines-used", type=int)
    create.add_argument("--qc-worst-residual-px", type=float)

    historical = commands.add_parser(
        "import-historical",
        help="Convert a bundled historical binder into a registrable snapshot folder.",
    )
    historical.add_argument(
        "binder",
        metavar="ID",
        help="Bundled binder identity (20190529_cmos, 2024-03-05, ...) or a binder TOML path.",
    )
    historical.add_argument(
        "--calibrations",
        default="calibrations",
        metavar="DIR",
        help="Snapshot root the imported epoch is written into (default: calibrations).",
    )
    historical.add_argument(
        "--artifact-root",
        action="append",
        default=[],
        metavar="DIR",
        help=(
            "Extra folder searched for the artifacts a binder names, such as the "
            "campaign folder holding a sphere pair too large to package. Repeat as needed."
        ),
    )
    historical.add_argument(
        "--valid-from",
        metavar="YYYY-MM-DD",
        help="Epoch start (default: the binder's acquired_date).",
    )
    historical.add_argument("--valid-to", metavar="YYYY-MM-DD")
    historical.add_argument("--shot-from", type=int)
    historical.add_argument("--shot-to", type=int)

    validate = commands.add_parser(
        "validate", help="Verify schema, required roles, paths, sizes, and SHA-256 digests."
    )
    validate.add_argument("snapshot", metavar="DIR")

    show = commands.add_parser("show", help="Print a compact snapshot summary.")
    show.add_argument("snapshot", metavar="DIR")
    show.add_argument(
        "--no-digest-check",
        action="store_true",
        help="Read metadata without hashing artifact contents.",
    )
    return parser


def _values(**values):
    selected = {key: value for key, value in values.items() if value is not None}
    return selected or None


def _print_validation_error(exc: SnapshotValidationError) -> None:
    print("Snapshot is invalid:")
    for error in exc.errors:
        print(f"  - {error}")


def _lamp_files(values: list[str]) -> list[tuple[str, str]]:
    parsed = []
    for value in values:
        label, separator, path = value.partition("=")
        if not separator or not label.strip() or not path.strip():
            raise SnapshotError("--lamp-file must be NAME=PATH, for example H2=h2.sif")
        parsed.append((label.strip(), path.strip()))
    return parsed


def _print_summary(snapshot: Snapshot) -> None:
    """Print the compact ``show`` view, references and all."""

    print(f"id:        {snapshot.snapshot_id}")
    print(f"detector:  {snapshot.detector}")
    print(f"lamps:     {', '.join(snapshot.lamps)}")
    print(f"artifacts: {len(snapshot.artifacts)}")
    for artifact in snapshot.artifacts:
        if not artifact.is_reference:
            continue
        # Where the frames really are, said out loud: an operator asking what a
        # snapshot uses should not have to open the binder to learn that the
        # bytes are one folder over.
        label = f"{artifact.role}/{artifact.label}" if artifact.label else artifact.role
        print(f"reference: {label} -> {snapshot.path_for(artifact)}")
    base = snapshot.manifest.get("base_snapshot")
    if base:
        print(f"base:      {base}")
    validity = snapshot.manifest.get("validity")
    if validity:
        values = ", ".join(f"{key}={value}" for key, value in validity.items())
        print(f"validity:  {values}")


def main(argv: list[str] | None = None, *, prog: str = "echelle-snapshot") -> int:
    """Run the snapshot CLI and return a process exit code."""

    parser = _build_parser(prog)
    args = parser.parse_args(argv)
    if args.snapshot_command is None:
        parser.print_help()
        return 0

    try:
        if args.snapshot_command == "create":
            files = {
                "pattern": args.pattern,
                "wavelength": args.wavelength,
                "sphere": args.sphere,
                "sphere_background": args.sphere_background,
                "integral": args.integral,
            }
            snapshot = create_snapshot(
                args.output_root,
                snapshot_id=args.snapshot_id,
                detector=args.detector,
                files=files,
                lamps=args.lamp_used,
                lamp_files=_lamp_files(args.lamp_file),
                notes=args.notes,
                base_snapshot=args.base_snapshot,
                validity=_values(
                    shot_from=args.shot_from,
                    shot_to=args.shot_to,
                    date_from=args.date_from,
                    date_to=args.date_to,
                ),
                alignment=_values(
                    dx_px=args.alignment_dx_px,
                    dy_px=args.alignment_dy_px,
                    rotation_deg=args.alignment_rotation_deg,
                    rms_px=args.alignment_rms_px,
                ),
                qc=_values(
                    lines_used=args.qc_lines_used,
                    worst_residual_px=args.qc_worst_residual_px,
                ),
            )
            print(f"Created calibration snapshot {snapshot.snapshot_id}")
            print(f"  {snapshot.root}")
            print(f"  {len(snapshot.artifacts)} artifact(s), all digests verified")
            return 0

        if args.snapshot_command == "import-historical":
            from .historical import HistoricalError, import_historical_snapshot

            try:
                snapshot = import_historical_snapshot(
                    args.binder,
                    args.calibrations,
                    valid_from=args.valid_from,
                    valid_to=args.valid_to,
                    shot_from=args.shot_from,
                    shot_to=args.shot_to,
                    artifact_roots=tuple(args.artifact_root),
                )
            except HistoricalError as exc:
                print(f"Historical binder was not imported: {exc}", file=sys.stderr)
                return 1
            validity = snapshot.manifest.get("validity") or {}
            print(f"Imported historical calibration {snapshot.snapshot_id}")
            print(f"  {snapshot.root}")
            print(f"  {len(snapshot.artifacts)} artifact(s), all digests verified")
            print(
                "  validity: "
                + ", ".join(f"{key}={value}" for key, value in validity.items())
            )
            print(f"  register it by adding [[epochs]] snapshot_id = {snapshot.snapshot_id!r}")
            return 0

        snapshot = load_snapshot(
            Path(args.snapshot),
            verify_files=not getattr(args, "no_digest_check", False),
        )
        references = [artifact for artifact in snapshot.artifacts if artifact.is_reference]
        if args.snapshot_command == "validate":
            print(
                f"Snapshot {snapshot.snapshot_id} is valid: "
                f"{len(snapshot.artifacts)} artifact(s), all digests verified"
            )
            if references:
                print(
                    f"  {len(references)} referenced source(s) verified where they live"
                )
            return 0

        _print_summary(snapshot)
        return 0
    except SnapshotValidationError as exc:
        _print_validation_error(exc)
        return 1
    except SnapshotError as exc:
        print(f"Snapshot was not created: {exc}")
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
