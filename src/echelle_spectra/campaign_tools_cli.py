"""Command adapters for the Packets 9--13 campaign domain modules."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Bounded refusals
#
# Every command below reads paths an operator typed, and a typed path is wrong
# often enough that the wrong one has to be answered rather than raised.  A
# copied documentation example is the sharpest case: ``/data/all-years.json``
# pasted into PowerShell is a perfectly legal path on the current drive, so the
# process does not fail where the operator is looking -- it fails several
# frames inside pathlib, and prints a traceback about a file nobody named.
#
# One ``CommandError`` is the finished sentence for that operator: what was
# missing, the absolute path this machine actually looked at, the flag that
# supplied it, and the command that writes such a file.  ``_bounded`` prints it
# to stderr and leaves with a nonzero status.  Library modules keep raising
# their own precise exceptions; translating them is this layer's job.
#
# Messages stay ASCII on purpose: these run on institute consoles whose code
# page is cp1252 or cp932, where a decorative dash can turn a diagnosis into a
# UnicodeEncodeError.
# ---------------------------------------------------------------------------


class CommandError(Exception):
    """One operator-facing refusal, printed as a single line without a traceback."""


def _bounded(run: Callable[[], int]) -> int:
    """Run one command body, turning its refusals into one line and exit 1."""

    try:
        return run()
    except CommandError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def _resolved(value: str | Path) -> Path:
    """Return the absolute path *this* machine reads ``value`` as.

    A Unix-shaped absolute path copied onto Windows silently becomes a path on
    the current drive.  The answer therefore has to show where the process
    looked, never merely echo what was typed.
    """

    try:
        return Path(value).expanduser().resolve()
    except (OSError, ValueError, RuntimeError):  # pragma: no cover - exotic paths
        return Path(str(value))


def _from(flag: str) -> str:
    return f"(from {flag})"


def _tail(remedy: str) -> str:
    return f" -- {remedy}" if remedy else ""


def _require_file(value: str | Path, *, flag: str, what: str, remedy: str = "") -> Path:
    """Return the resolved path of an input file, or refuse in one line."""

    path = _resolved(value)
    if path.is_file():
        return path
    if path.is_dir():
        raise CommandError(f"{what} is a folder, not a file: {path} {_from(flag)}")
    raise CommandError(f"{what} not found: {path} {_from(flag)}{_tail(remedy)}")


def _require_dir(value: str | Path, *, flag: str, what: str, remedy: str = "") -> Path:
    """Return the resolved path of an input folder, or refuse in one line."""

    path = _resolved(value)
    if path.is_dir():
        return path
    if path.is_file():
        raise CommandError(f"{what} is a file, not a folder: {path} {_from(flag)}")
    raise CommandError(f"{what} not found: {path} {_from(flag)}{_tail(remedy)}")


def _require_input(value: str | Path, *, flag: str, what: str, remedy: str = "") -> Path:
    """Return a path that may be either a file or a folder, or refuse."""

    path = _resolved(value)
    if path.exists():
        return path
    raise CommandError(f"{what} not found: {path} {_from(flag)}{_tail(remedy)}")


def _require_parent(value: str | Path, *, flag: str, what: str) -> Path:
    """Return the resolved output path whose parent folder already exists.

    The commands here happily create missing parents, which is exactly how a
    mistyped root ends up writing a whole tree somewhere nobody asked for.
    """

    path = _resolved(value)
    parent = path.parent
    if parent.is_dir():
        return path
    if parent.is_file():
        raise CommandError(f"{what} cannot be written under a file: {parent} {_from(flag)}")
    raise CommandError(
        f"{what} has no parent folder: {parent} {_from(flag)} -- create it first, "
        f"or point {flag} inside a folder that exists"
    )


def _parse_json(path: Path, *, flag: str, what: str) -> Any:
    """Read one JSON input, reporting a parse failure as one line."""

    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise CommandError(f"{what} is not valid JSON: {path} {_from(flag)} -- {exc}") from None
    except (OSError, UnicodeDecodeError) as exc:
        raise CommandError(f"{what} could not be read: {path} {_from(flag)} -- {exc}") from None


def _parse_toml(path: Path, *, flag: str, what: str) -> Any:
    """Read one TOML input, reporting a parse failure as one line."""

    import tomllib

    try:
        with path.open("rb") as handle:
            return tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise CommandError(f"{what} is not valid TOML: {path} {_from(flag)} -- {exc}") from None
    except OSError as exc:
        raise CommandError(f"{what} could not be read: {path} {_from(flag)} -- {exc}") from None


def _next_free(path: Path) -> Path:
    """Return the next unused ``-N`` sibling of an immutable output path."""

    base = re.sub(r"-\d+$", "", path.stem) or path.stem
    for index in range(2, 1000):
        candidate = path.with_name(f"{base}-{index}{path.suffix}")
        if not candidate.exists():
            return candidate
    return path.with_name(f"{base}-{path.stat().st_mtime_ns}{path.suffix}")  # pragma: no cover


CATALOG_REMEDY = "build one with `echelle catalog build`"
MERGED_REMEDY = "build the per-drive catalogs with `echelle catalog build`, then `echelle catalog merge`"
CUBE_REMEDY = "write cubes with `echelle process`"
SNAPSHOT_REMEDY = "create one with `echelle snapshot create`, or the live bench `echelle-calib`"
CALIBRATIONS_REMEDY = "snapshots are written there by `echelle snapshot create` and by `echelle-calib`"
REGISTRY_REMEDY = "write the ordered epoch registry by hand; `echelle status` reports what it reads"
DRIFT_REMEDY = "write one with `echelle drift audit`"


def txt_main(argv: list[str] | None = None, *, prog: str = "echelle-cube2txt") -> int:
    from .lhd_text import LhdTextError, write_cube_text

    parser = argparse.ArgumentParser(
        prog=prog, description="Write frozen-header LHD text from a saved cube."
    )
    parser.add_argument("input", help="Saved SpectroCube .nc file.")
    parser.add_argument("output", help="Destination text file.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)

    def run() -> int:
        cube = _require_file(
            args.input, flag="the INPUT argument", what="cube file", remedy=CUBE_REMEDY
        )
        destination = _require_parent(
            args.output, flag="the OUTPUT argument", what="text output"
        )
        try:
            path = write_cube_text(cube, destination, overwrite=args.overwrite)
        except FileExistsError as exc:
            raise CommandError(f"{exc} -- pass --overwrite to replace it") from None
        except (LhdTextError, OSError) as exc:
            raise CommandError(str(exc)) from None
        print(path)
        return 0

    return _bounded(run)


def catalog_main(argv: list[str] | None = None, *, prog: str = "echelle catalog") -> int:
    from .catalog import build_drive_catalog, merge_catalogs

    parser = argparse.ArgumentParser(prog=prog, description="Build and merge portable cube catalogs.")
    commands = parser.add_subparsers(dest="action", metavar="COMMAND")
    build = commands.add_parser("build", help="Write one catalog beside a drive's cubes.")
    build.add_argument("cubes")
    build.add_argument("--volume-label", required=True)
    build.add_argument(
        "--drive-id",
        metavar="ID",
        help=(
            "Stable drive identity for this catalog. Default: the id announced by "
            "echelle-drive-id.toml at or above the cube folder."
        ),
    )
    build.add_argument("--receipt-dir")
    build.add_argument("--output")
    merge = commands.add_parser("merge", help="Write an all-years index.")
    merge.add_argument("catalogs", nargs="+")
    merge.add_argument("-o", "--output", required=True)
    args = parser.parse_args(argv)
    if args.action is None:
        parser.print_help()
        return 0

    def run() -> int:
        if args.action == "build":
            cubes = _require_dir(
                args.cubes,
                flag="the CUBES argument",
                what="cube folder",
                remedy=CUBE_REMEDY,
            )
            receipts = (
                _require_dir(
                    args.receipt_dir,
                    flag="--receipt-dir",
                    what="run receipt folder",
                    remedy="`echelle process --runs-dir DIR` writes receipts there",
                )
                if args.receipt_dir
                else None
            )
            output = (
                _require_parent(args.output, flag="--output", what="catalog")
                if args.output
                else None
            )
            path = build_drive_catalog(
                cubes,
                volume_label=args.volume_label,
                drive_id=args.drive_id,
                receipt_dir=receipts,
                output=output,
            )
        else:
            inputs = []
            for raw in args.catalogs:
                candidate = _require_file(
                    raw,
                    flag="the CATALOGS argument",
                    what="catalog file",
                    remedy=CATALOG_REMEDY,
                )
                _parse_json(candidate, flag="the CATALOGS argument", what="catalog file")
                inputs.append(candidate)
            output = _require_parent(args.output, flag="--output", what="merged catalog")
            try:
                path = merge_catalogs(inputs, output)
            except ValueError as exc:
                raise CommandError(str(exc)) from None
            except KeyError as exc:
                raise CommandError(f"catalog file is missing {exc} and cannot be merged") from None
        print(path)
        return 0

    return _bounded(run)


def recal_main(argv: list[str] | None = None, *, prog: str = "echelle recal-cube") -> int:
    from .recalibration import RecalibrationError, recalibrate_cube
    from .snapshot import SnapshotValidationError

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

    def run() -> int:
        source = _require_file(
            args.input, flag="the INPUT argument", what="cube file", remedy=CUBE_REMEDY
        )
        snapshot = _require_input(
            args.new_snapshot,
            flag="--new-snapshot",
            what="calibration snapshot",
            remedy=SNAPSHOT_REMEDY,
        )
        destination = _require_parent(args.output, flag="--output", what="recalibrated cube")
        try:
            output, manifest = recalibrate_cube(
                source,
                destination,
                new_snapshot_path=snapshot,
                update_wavelength=not args.factor_only,
                update_factor=not args.wavelength_only,
                overwrite=args.overwrite,
            )
        except FileExistsError as exc:
            raise CommandError(f"{exc} -- pass --overwrite to replace it") from None
        except (RecalibrationError, SnapshotValidationError, OSError) as exc:
            raise CommandError(str(exc)) from None
        print(f"cube: {output}\nmanifest: {manifest}")
        return 0

    return _bounded(run)


def drift_main(argv: list[str] | None = None, *, prog: str = "echelle drift") -> int:
    parser = argparse.ArgumentParser(prog=prog, description="Audit science-line drift and accept refinements.")
    commands = parser.add_subparsers(dest="action", metavar="COMMAND")
    audit = commands.add_parser(
        "audit", help="Measure sampled cubes and write one immutable verdict file."
    )
    audit.add_argument(
        "cubes",
        nargs="+",
        metavar="CUBE_OR_DIR",
        help="Saved .nc cubes, or a directory whose .nc cubes are all audited.",
    )
    audit.add_argument("--every", type=int, default=1, help="Audit every Nth selected cube.")
    audit.add_argument(
        "--shot",
        action="append",
        default=[],
        metavar="SHOT",
        help="Also audit this exact shot; 42 matches shot 42 and never 142.",
    )
    audit.add_argument(
        "--from",
        dest="date_from",
        metavar="YYYY-MM-DD",
        help="Audit cubes acquired on or after this date (cube t_start, else the "
        "date in source_file, else created_at).",
    )
    audit.add_argument(
        "--to",
        dest="date_to",
        metavar="YYYY-MM-DD",
        help="Audit cubes acquired on or before this date.",
    )
    audit.add_argument(
        "--catalog",
        metavar="JSON",
        help="Merged catalog used to name the drives holding beyond-repair shots.",
    )
    audit.add_argument(
        "--calibrations",
        default=None,
        metavar="DIR",
        help="Snapshot root named by the composed repair commands (default: calibrations).",
    )
    audit.add_argument("-o", "--output", required=True)
    refine = commands.add_parser(
        "refine", help="Accept a shifted verdict and emit an immutable -rN snapshot."
    )
    refine.add_argument("evidence")
    refine.add_argument("--calibrations", required=True)
    refine.add_argument(
        "--accept-shift",
        required=True,
        type=float,
        metavar="PIXELS",
        help="Exactly acknowledge the sampled median detector shift, in pixels.",
    )
    args = parser.parse_args(argv)
    if args.action is None:
        parser.print_help()
        return 0

    def run() -> int:
        if args.action == "audit":
            return _drift_audit(args)
        return _drift_refine(args)

    return _bounded(run)


def _immutable_evidence(output: Path) -> CommandError:
    """The one refusal an already-taken evidence name earns, with the next name."""

    return CommandError(
        f"drift evidence is immutable and already exists: {output} (from --output) -- "
        f"an audit never overwrites its verdict; write this one to {_next_free(output)}"
    )


def _drift_audit_inputs(args: argparse.Namespace) -> dict[str, Any]:
    """Resolve every path the audit reads or writes before any measuring starts."""

    catalog = None
    if args.catalog:
        catalog = _require_file(
            args.catalog, flag="--catalog", what="catalog file", remedy=MERGED_REMEDY
        )
        _parse_json(catalog, flag="--catalog", what="catalog file")
    calibrations: str | Path = "calibrations"
    if args.calibrations:
        calibrations = _require_dir(
            args.calibrations,
            flag="--calibrations",
            what="calibration snapshot root",
            remedy=CALIBRATIONS_REMEDY,
        )
    output = _require_parent(args.output, flag="--output", what="drift evidence")
    # The audit is minutes of work and its evidence file is immutable by
    # design, so an already-taken name is answered before the measuring starts
    # -- and answered with the name the second audit should use, since an
    # insufficient-data first audit makes a second audit the normal path.
    if output.exists():
        raise _immutable_evidence(output)
    return {
        "cubes": [
            _require_input(raw, flag="the CUBE_OR_DIR argument", what="cube", remedy=CUBE_REMEDY)
            for raw in args.cubes
        ],
        "catalog": catalog,
        "calibrations": calibrations,
        "output": output,
    }


def _report_audit(payload: dict[str, Any], path: Path) -> int:
    """Print one finished verdict and the repair commands it composed."""

    print(f"{payload['verdict']}: {path}")
    if payload.get("interval_warning"):
        print(f"warning: {payload['interval_warning']}")
    if payload.get("data_requirement"):
        print(payload["data_requirement"])
    for step in payload.get("repair_commands", []):
        shell = "" if step["shell"] == "any" else f"[{step['shell']}] "
        print(f"# {shell}{step['purpose']}")
        if step["command"]:
            print(step["command"])
    return 0 if payload["verdict"] in {"aligned", "shifted"} else 1


def _drift_audit(args: argparse.Namespace) -> int:
    from .drift import DriftError, audit_cubes, write_drift_evidence

    resolved = _drift_audit_inputs(args)
    cubes = resolved["cubes"]
    catalog = resolved["catalog"]
    calibrations = resolved["calibrations"]
    output = resolved["output"]
    try:
        payload = audit_cubes(
            cubes,
            every=args.every,
            shots=set(args.shot),
            date_from=args.date_from,
            date_to=args.date_to,
            catalog=catalog,
            evidence_path=output,
            calibrations_root=calibrations,
        )
    except DriftError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        raise CommandError(str(exc)) from None
    try:
        path = write_drift_evidence(output, payload)
    except FileExistsError:
        raise _immutable_evidence(output) from None
    except OSError as exc:
        raise CommandError(str(exc)) from None
    return _report_audit(payload, path)


def _drift_refine(args: argparse.Namespace) -> int:
    from .drift import DriftError, create_refinement_snapshot

    evidence = _require_file(
        args.evidence, flag="the EVIDENCE argument", what="drift evidence", remedy=DRIFT_REMEDY
    )
    _parse_json(evidence, flag="the EVIDENCE argument", what="drift evidence")
    calibrations = _require_dir(
        args.calibrations,
        flag="--calibrations",
        what="calibration snapshot root",
        remedy=CALIBRATIONS_REMEDY,
    )
    try:
        snapshot, accepted = create_refinement_snapshot(
            evidence,
            calibrations_root=calibrations,
            accepted_shift_px=args.accept_shift,
        )
    except (DriftError, FileExistsError, OSError) as exc:
        raise CommandError(str(exc)) from None
    print(f"snapshot: {snapshot.snapshot_id}\naccepted verdict: {accepted}")
    return 0


def web_main(argv: list[str] | None = None, *, prog: str = "echelle web") -> int:
    from .catalog import load_catalog
    from .reading_room import build_reading_room

    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Build the read-only campaign page: the Now stepper, the drives, the "
            "calibration evidence, and the packaged reading room."
        ),
    )
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--drift", action="append", default=[])
    parser.add_argument(
        "--document",
        action="append",
        default=[],
        metavar="MARKDOWN",
        help=(
            "Extra Markdown rendered after the packaged vocabulary, procedure and "
            "provenance documents, which are always included."
        ),
    )
    parser.add_argument(
        "--registry",
        metavar="TOML",
        help="Ordered epoch registry read to pre-fill the composer and name the epochs.",
    )
    parser.add_argument(
        "--calibrations",
        metavar="DIR",
        help="Snapshot root used with --registry (default: calibrations beside the registry).",
    )
    args = parser.parse_args(argv)

    def run() -> int:
        catalog = _require_file(
            args.catalog, flag="--catalog", what="catalog file", remedy=CATALOG_REMEDY
        )
        _parse_json(catalog, flag="--catalog", what="catalog file")
        try:
            load_catalog(catalog)
        except ValueError as exc:
            raise CommandError(f"{exc} (from --catalog)") from None
        registry = None
        if args.registry:
            registry = _require_file(
                args.registry, flag="--registry", what="epoch registry", remedy=REGISTRY_REMEDY
            )
            _parse_toml(registry, flag="--registry", what="epoch registry")
        calibrations = None
        if args.calibrations:
            calibrations = _require_dir(
                args.calibrations,
                flag="--calibrations",
                what="calibration snapshot root",
                remedy=CALIBRATIONS_REMEDY,
            )
        drift = []
        for raw in args.drift:
            evidence = _require_file(
                raw, flag="--drift", what="drift evidence", remedy=DRIFT_REMEDY
            )
            _parse_json(evidence, flag="--drift", what="drift evidence")
            drift.append(evidence)
        documents = [
            _require_file(raw, flag="--document", what="Markdown document")
            for raw in args.document
        ]
        output = _require_parent(args.output, flag="--output", what="campaign page folder")
        try:
            path = build_reading_room(
                catalog,
                output,
                drift_paths=drift,
                document_paths=documents,
                registry_path=registry,
                calibrations_root=calibrations,
            )
        except OSError as exc:
            raise CommandError(str(exc)) from None
        print(path)
        return 0

    return _bounded(run)


def historical_main(argv: list[str] | None = None, *, prog: str = "echelle historical") -> int:
    from .historical import bundled_historical_manifests

    parser = argparse.ArgumentParser(prog=prog, description="Validate thin historical calibration binders.")
    parser.parse_args(argv)
    payload = bundled_historical_manifests()
    print(json.dumps([item.snapshot_id for item in payload]))
    return 0
