"""Command adapters for the Packets 9--13 campaign domain modules."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import webbrowser
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
    audit.add_argument(
        "--every",
        type=int,
        default=None,
        help=(
            "Audit every Nth selected cube. Default: derived so about 20 cubes are "
            "measured, max(1, cubes // 20)."
        ),
    )
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
    audit.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="JSON",
        help=(
            "Evidence file to write. Default: the next free drift-evidence-NNN.json "
            "in the folder the audited cubes share."
        ),
    )
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


EVIDENCE_STEM = "drift-evidence"
AUDIT_TARGET_CUBES = 20


def _evidence_home(cubes: list[Path]) -> Path:
    """Return the folder an unnamed evidence file belongs beside.

    The evidence describes the cubes, so it is written where the cubes are: the
    deepest folder all of them share.  Cubes spread across two drives share
    nothing on Windows, and then the first one's folder is the honest answer.
    """

    folders = [path if path.is_dir() else path.parent for path in cubes]
    if not folders:  # pragma: no cover - the audit always has at least one cube
        return _resolved(Path.cwd())
    try:
        return Path(os.path.commonpath([str(folder) for folder in folders]))
    except ValueError:
        return folders[0]


def _derived_evidence_path(cubes: list[Path]) -> Path:
    """Name the next free ``drift-evidence-NNN.json`` beside the audited cubes.

    Numbered from 001 upward, so a second audit of the same cubes never has to
    argue with the immutable first one: it simply takes the next name.
    """

    home = _evidence_home(cubes)
    for index in range(1, 1000):
        candidate = home / f"{EVIDENCE_STEM}-{index:03d}.json"
        if not candidate.exists():
            return candidate
    raise CommandError(  # pragma: no cover - 999 audits of one folder
        f"every derived evidence name from {EVIDENCE_STEM}-001.json to "
        f"{EVIDENCE_STEM}-999.json is taken in {home} -- name this one with --output"
    )


def _derived_every(cubes: list[Path]) -> int:
    """Derive the sampling interval that measures about twenty of these cubes."""

    from .drift import DriftError, resolve_cube_paths

    try:
        count = len(resolve_cube_paths(list(cubes)))
    except (DriftError, OSError):
        # The audit itself answers an empty or unreadable selection; deriving an
        # interval is not the place to refuse it.
        return 1
    return max(1, count // AUDIT_TARGET_CUBES)


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
    output = None
    if args.output:
        output = _require_parent(args.output, flag="--output", what="drift evidence")
        # The audit is minutes of work and its evidence file is immutable by
        # design, so an already-taken name is answered before the measuring
        # starts -- and answered with the name the second audit should use,
        # since an insufficient-data first audit makes a second audit normal.
        if output.exists():
            raise _immutable_evidence(output)
    cubes = [
        _require_input(raw, flag="the CUBE_OR_DIR argument", what="cube", remedy=CUBE_REMEDY)
        for raw in args.cubes
    ]
    if output is None:
        output = _derived_evidence_path(cubes)
        print(f"evidence: {output} (derived; name one with --output)")
    return {
        "cubes": cubes,
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
    every = args.every
    if every is None:
        every = _derived_every(cubes)
        if every > 1:
            print(f"every: {every} (derived; about {AUDIT_TARGET_CUBES} cubes are measured)")
    try:
        payload = audit_cubes(
            cubes,
            every=every,
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


CAMPAIGN_HOME_NAME = "campaign.toml"
CAMPAIGN_HOME_KEYS = ("catalog", "output", "registry", "calibrations")
CAMPAIGN_HOME_REMEDY = (
    "a campaign home is a TOML file naming catalog, output, registry, calibrations "
    "and drift; write one by hand beside the campaign"
)


def _home_path(base: Path, value: str) -> str:
    """Resolve one path written inside a campaign home against the home's folder.

    The file is meant to be copied beside a campaign and read from anywhere, so
    what it says is relative to *itself*, never to whatever folder the operator
    happened to be standing in.
    """

    candidate = Path(value).expanduser()
    return str(candidate if candidate.is_absolute() else base / candidate)


def _campaign_home(value: str | Path, *, flag: str) -> tuple[Path, dict[str, Any]]:
    """Read one campaign home and return it with the defaults it supplies."""

    candidate = _resolved(value)
    if candidate.is_dir():
        candidate = candidate / CAMPAIGN_HOME_NAME
    path = _require_file(
        candidate, flag=flag, what="campaign home", remedy=CAMPAIGN_HOME_REMEDY
    )
    data = _parse_toml(path, flag=flag, what="campaign home")
    if not isinstance(data, dict):  # pragma: no cover - tomllib always returns a table
        raise CommandError(f"campaign home is not a table of keys: {path} {_from(flag)}")
    base = path.parent
    values: dict[str, Any] = {}
    for key in CAMPAIGN_HOME_KEYS:
        if key not in data:
            continue
        raw = data[key]
        if not isinstance(raw, str) or not raw.strip():
            raise CommandError(
                f"campaign home key '{key}' must be one path written as a string: "
                f'{path} {_from(flag)} -- for example {key} = "all-years.json"'
            )
        values[key] = _home_path(base, raw)
    if "drift" in data:
        raw = data["drift"]
        if not isinstance(raw, list) or not all(isinstance(item, str) for item in raw):
            raise CommandError(
                f"campaign home key 'drift' must be a list of paths: {path} {_from(flag)} "
                '-- for example drift = ["epoch-drift.json"]'
            )
        values["drift"] = [_home_path(base, item) for item in raw]
    # Keys this page does not read are left alone rather than refused: the home
    # is hand-edited, and an operator's own notes are not an error.
    return path, values


def _web_home(args: argparse.Namespace) -> tuple[dict[str, Any], str]:
    """Read the campaign home this build stands in, and say what it supplied.

    Two answers are asked of the operator -- which catalog, and where to write
    -- and a home file beside the campaign answers both, so the second build of
    the day is ``echelle web --open``.  A home is read without --home only when
    something is actually missing, so a fully typed command is never quietly
    steered by a file in the current folder.
    """

    flag = "--home"
    if args.home:
        source, home = _campaign_home(args.home, flag=flag)
    elif not (args.catalog and args.output) and (Path.cwd() / CAMPAIGN_HOME_NAME).is_file():
        flag = f"the {CAMPAIGN_HOME_NAME} in this folder"
        source, home = _campaign_home(Path.cwd() / CAMPAIGN_HOME_NAME, flag=flag)
    else:
        return {}, flag
    print(f"campaign home: {source} (supplies {', '.join(sorted(home)) or 'nothing'})")
    return home, flag


def _web_inputs(args: argparse.Namespace, parser: argparse.ArgumentParser) -> dict[str, Any]:
    """Resolve every path the page is built from, whoever supplied it."""

    from .catalog import load_catalog

    home, home_flag = _web_home(args)

    def chosen(key: str, flag: str) -> tuple[Any, str]:
        typed = getattr(args, key)
        return (typed, flag) if typed else (home.get(key), home_flag)

    catalog_value, catalog_flag = chosen("catalog", "--catalog")
    output_value, output_flag = chosen("output", "--output")
    missing = [
        flag
        for flag, value in (("--catalog", catalog_value), ("--output", output_value))
        if not value
    ]
    if missing:
        parser.error("the following arguments are required: " + ", ".join(missing))

    catalog = _require_file(
        catalog_value, flag=catalog_flag, what="catalog file", remedy=CATALOG_REMEDY
    )
    _parse_json(catalog, flag=catalog_flag, what="catalog file")
    try:
        load_catalog(catalog)
    except ValueError as exc:
        raise CommandError(f"{exc} (from {catalog_flag})") from None

    registry = None
    registry_value, registry_flag = chosen("registry", "--registry")
    if registry_value:
        registry = _require_file(
            registry_value, flag=registry_flag, what="epoch registry", remedy=REGISTRY_REMEDY
        )
        _parse_toml(registry, flag=registry_flag, what="epoch registry")

    calibrations = None
    calibrations_value, calibrations_flag = chosen("calibrations", "--calibrations")
    if calibrations_value:
        calibrations = _require_dir(
            calibrations_value,
            flag=calibrations_flag,
            what="calibration snapshot root",
            remedy=CALIBRATIONS_REMEDY,
        )

    drift = []
    drift_values, drift_flag = chosen("drift", "--drift")
    for raw in drift_values or []:
        evidence = _require_file(raw, flag=drift_flag, what="drift evidence", remedy=DRIFT_REMEDY)
        _parse_json(evidence, flag=drift_flag, what="drift evidence")
        drift.append(evidence)

    return {
        "catalog": catalog,
        "output": _require_parent(output_value, flag=output_flag, what="campaign page folder"),
        "registry": registry,
        "calibrations": calibrations,
        "drift": drift,
        "documents": [
            _require_file(raw, flag="--document", what="Markdown document")
            for raw in args.document
        ],
    }


def web_main(argv: list[str] | None = None, *, prog: str = "echelle web") -> int:
    from .reading_room import build_reading_room

    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Build the read-only campaign page: the Now stepper, the drives, the "
            "calibration evidence, and the packaged reading room."
        ),
    )
    parser.add_argument(
        "--catalog",
        default=None,
        help="Merged catalog the page is built from (default: the campaign home's catalog).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Folder the page is written into (default: the campaign home's output).",
    )
    parser.add_argument(
        "--home",
        metavar="DIR",
        help=(
            "Campaign home supplying defaults for --catalog, --output, --registry, "
            "--calibrations and --drift: a folder holding campaign.toml, or the TOML "
            "file itself. Explicit flags win; paths inside it are relative to it. "
            "Without --home, a campaign.toml in the current folder is read when "
            "--catalog or --output is missing."
        ),
    )
    parser.add_argument(
        "--open",
        action="store_true",
        help="Open the built page in the default browser. The page stays a static file.",
    )
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
        resolved = _web_inputs(args, parser)
        try:
            path = build_reading_room(
                resolved["catalog"],
                resolved["output"],
                drift_paths=resolved["drift"],
                document_paths=resolved["documents"],
                registry_path=resolved["registry"],
                calibrations_root=resolved["calibrations"],
            )
        except OSError as exc:
            raise CommandError(str(exc)) from None
        # The absolute path is the whole delivery: the page is a file, and the
        # operator has to be able to find it, mail it, or open it by hand.
        page = _resolved(path)
        print(page)
        if args.open and not webbrowser.open(page.as_uri()):
            print(
                f"WARNING: no default browser answered; open {page} by hand",
                file=sys.stderr,
            )
        return 0

    return _bounded(run)


def historical_main(argv: list[str] | None = None, *, prog: str = "echelle historical") -> int:
    from .historical import bundled_historical_manifests

    parser = argparse.ArgumentParser(prog=prog, description="Validate thin historical calibration binders.")
    parser.parse_args(argv)
    payload = bundled_historical_manifests()
    print(json.dumps([item.snapshot_id for item in payload]))
    return 0
