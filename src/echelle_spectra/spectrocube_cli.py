r"""CLI for exporting Echelle SIF files to the SpectroCube NetCDF format.

Entry point installed as ``echelle-spectrocube`` (see pyproject.toml).

Usage examples
--------------
Single file::

    echelle-spectrocube shot_042.sif --units wm -o shot_042.nc

Batch folder (all .SIF files)::

    echelle-spectrocube D:\NIFS\shots --units wm -o D:\NIFS\cubes

Dry run to preview what would happen, naming every file it found::

    echelle-spectrocube D:\NIFS\shots --dry-run

The examples name a campaign drive the way the guides do: ``D:\NIFS`` on
Windows, ``/Volumes/NIFS`` on macOS with forward slashes.

Calibration authority
---------------------
Legacy conversion uses one bundled or manually configured calibration.  An
epoch registry instead resolves and verifies one immutable snapshot for every
source, allowing a batch to cross reviewed shot/date boundaries safely.
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .calibration_registry import (
    CalibrationEpochRegistry,
    CalibrationRegistryError,
    load_calibration_registry,
)
from .campaign_run import (
    GATE_SAMPLE,
    GATE_UNGATED,
    GATE_VERDICT,
    RESUME_TRUST_STAT,
    RunReceipt,
    default_volume_label,
    ensure_drive_identity,
    find_resumable_run,
    new_run_directory,
    sha256_file,
    target_runs_root,
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

# Units whose cubes claim an absolute flux scale, and therefore owe a named
# standard. Mirrors the "absolute" rows of tools.spectrocube_export._UNIT_INFO,
# which is imported lazily because it reaches for the optional exporter.
ABSOLUTE_UNITS = frozenset({"wm", "wmsr", "phmsr"})

_ExportStatus = Literal["exported", "skipped", "dry-run", "failed"]
_PRINT_LOCK = threading.Lock()

# Cubes produced by an unverified --sample run carry this attribute. NetCDF has
# no boolean attribute type, so the flag is stored as the legal integer true.
DRIFT_SAMPLE_ATTR = "drift_sample"
DRIFT_SAMPLE_TRUE = 1

# --sample auto: the operator says "take a first look", and the size is derived
# from what the drive turns out to hold rather than guessed before opening it.
SAMPLE_AUTO = "auto"
SAMPLE_AUTO_MIN = 5
SAMPLE_AUTO_MAX = 30
SAMPLE_AUTO_DIVISOR = 25

# A campaign drive is a tree of date-named day folders, and its own lamp frames
# live beside them: a folder called `calibrations`, or any folder holding the
# `snapshot.toml` of a saved calibration snapshot. Those are never science
# shots, so discovery prunes them whole rather than exporting cubes from them.
CALIBRATION_DIR_NAME = "calibrations"
SNAPSHOT_MARKER = "snapshot.toml"

# Pruning is a silent decision unless it is stated, so every run that prunes
# anything names it. A long drive can prune dozens of folders, and a header
# nobody reads to the end is another way of saying nothing: the console lists
# this many and counts the rest, while the receipt keeps the complete list.
PRUNED_LISTING_LIMIT = 5

# Documented example paths follow the convention the guides use for a campaign
# drive. Invented Unix roots (`/data/...`) are the copy-paste trap this repo
# already swept out of the docs: pasted into PowerShell they are legal paths on
# whatever drive the operator happens to be sitting on.
EXAMPLE_SHOTS = r"D:\NIFS\shots"
EXAMPLE_CUBES = r"D:\NIFS\cubes"
EXAMPLE_SECOND_SHOTS = r"E:\NIFS-B\shots"


@dataclass(frozen=True)
class ExportResult:
    status: _ExportStatus
    reason: str = ""


def _sample_size(value: str) -> int | str:
    """Parse ``--sample``: a positive count, or the derived-size word ``auto``."""

    text = str(value).strip()
    if text.lower() == SAMPLE_AUTO:
        return SAMPLE_AUTO
    try:
        return int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"invalid sample size: {value!r} -- pass a whole number of files, or 'auto'"
        ) from None


def _auto_sample_count(n_files: int) -> int:
    """Derive how many files a first look at ``n_files`` should export.

    Small folders are sampled whole -- there is nothing to sample from -- and
    everything larger gets one file per 25, held between 5 and 30 so the first
    run stays minutes rather than hours whatever the drive turns out to hold.
    """

    if n_files < SAMPLE_AUTO_MIN:
        return n_files
    derived = -(-n_files // SAMPLE_AUTO_DIVISOR)  # ceil without importing math
    return max(SAMPLE_AUTO_MIN, min(SAMPLE_AUTO_MAX, derived))


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def _build_parser(*, prog: str = "echelle-spectrocube") -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Export Echelle .sif files to SpectroCube NetCDF (.nc) format.\n\n"
            f"Single file:    {prog} shot.sif -o shot.nc\n"
            f"Batch folder:   {prog} {EXAMPLE_SHOTS} --units wm -o {EXAMPLE_CUBES}\n"
            f"Several drives: {prog} {EXAMPLE_SHOTS} {EXAMPLE_SECOND_SHOTS} "
            f"-o {EXAMPLE_CUBES}\n\n"
            "On macOS the same drives are /Volumes/NIFS/shots and /Volumes/NIFS-B/shots."
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
        "additional_inputs",
        nargs="*",
        metavar="INPUT",
        help="Additional source folders; one sequential worker runs per source.",
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
            "With several inputs, each source gets a named child directory.  "
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
        help=(
            "Glob pattern for batch SIF discovery (default: *.SIF). A plain filename "
            "pattern searches the whole source tree, skipping calibration folders; "
            "every skipped folder is named in the batch header and recorded in the "
            "run receipt. A junction or symbolic link sitting among the day folders "
            "is matched at its own level but is never descended into, so a linked "
            "day contributes nothing below itself. A pattern containing / or ** is "
            "used exactly as typed."
        ),
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
        help="Root for durable batch receipts and per-target trees (default: local/runs).",
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
        action="append",
        metavar="LABEL",
        help=(
            "USB volume label recorded in receipts (default: drive/root identity). "
            "Repeat once per source when processing several targets."
        ),
    )
    p.add_argument(
        "--central-index",
        default=None,
        metavar="JSON",
        help=(
            "Durable all-years index each completed target's per-drive catalog is "
            "merged into automatically. Without it the merge stays manual."
        ),
    )
    p.add_argument(
        "--registry",
        default=None,
        metavar="TOML",
        help=(
            "Ordered calibration epoch registry. Referenced snapshot [validity] "
            "boundaries select exactly one immutable snapshot per source."
        ),
    )
    p.add_argument(
        "--calibrations",
        default=None,
        metavar="DIR",
        help=(
            "Snapshot root used with --registry (default: calibrations beside the registry)."
        ),
    )
    p.add_argument(
        "--drift-verdict",
        default=None,
        metavar="JSON",
        help=(
            "Sampled drift evidence required before any registry-backed processing; "
            "insufficient or unaccepted shifted evidence is refused."
        ),
    )
    p.add_argument(
        "--sample",
        type=_sample_size,
        default=None,
        metavar="N",
        help=(
            "Legal first registry run: process the first N resolved files without a "
            "verdict. The receipt and every produced cube are marked as an unverified "
            "sample, which 'echelle drift audit' then turns into --drift-verdict evidence. "
            "'auto' derives N from the files this drive actually holds: "
            "max(5, min(30, ceil(files / 25))), or every file when the drive holds "
            "fewer than 5."
        ),
    )
    return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _is_calibration_dir(path: Path) -> bool:
    """Say whether a folder holds calibration frames rather than science shots."""

    if path.name.lower() == CALIBRATION_DIR_NAME:
        return True
    try:
        return (path / SNAPSHOT_MARKER).is_file()
    except OSError:  # pragma: no cover - unreadable mount points
        return False


def _science_tree(root: Path) -> tuple[list[Path], list[Path]]:
    """Split ``root``'s tree into the folders discovery searches and the ones it prunes.

    A matched folder prunes its whole subtree: a snapshot's sources sit inside
    it, and none of them are shots.  The named root itself is always searched --
    an operator who points the command at a folder has already said what it is.

    Only the topmost pruned folder of each subtree is returned, which is what a
    reader wants named: everything below it went with it.  ``os.walk`` does not
    follow links, so a junction or symbolic link among the day folders is
    matched where it sits but is never descended into.
    """

    searched = [root]
    pruned: list[Path] = []
    for current, dirnames, _ in os.walk(root):
        here = Path(current)
        keep: list[str] = []
        for name in sorted(dirnames):
            if _is_calibration_dir(here / name):
                pruned.append(here / name)
            else:
                keep.append(name)
        dirnames[:] = keep
        searched.extend(here / name for name in keep)
    return searched, pruned


def _pattern_is_recursive(pattern: str) -> bool:
    """Say whether a pattern is a plain filename pattern this command may walk."""

    return "/" not in pattern and "**" not in pattern


@dataclass(frozen=True)
class Discovery:
    """What one batch target's search found, and what it deliberately skipped."""

    files: tuple[Path, ...] = ()
    pruned: tuple[Path, ...] = ()


def _discover_science_files(input_path: Path, pattern: str) -> Discovery:
    """Find the science files one batch target offers, in a deterministic order.

    Campaign drives are trees of date-named day folders, so a plain filename
    pattern searches the whole tree; a pattern the operator wrote with `/` or
    `**` is theirs and is used exactly as typed.  The pruned folders come back
    with the files because a skip nobody is told about is indistinguishable
    from a drive that never held those shots.
    """

    if not _pattern_is_recursive(pattern):
        return Discovery(files=tuple(sorted(input_path.glob(pattern))))

    directories, pruned = _science_tree(input_path)
    # The lowercase retry is the same fallback the flat search always had: it
    # runs only when the typed case found nothing, so a case-insensitive
    # filesystem -- where `*.SIF` already matched `a.sif` -- never sees a file
    # twice.  Normcase keys keep that promise even if it ever did.
    attempts = [pattern] + (["*.sif"] if pattern == "*.SIF" else [])
    for attempt in attempts:
        found: dict[str, Path] = {}
        for directory in directories:
            for path in directory.glob(attempt):
                if path.is_file():
                    found.setdefault(os.path.normcase(str(path)), path)
        if found:
            return Discovery(files=tuple(sorted(found.values())), pruned=tuple(pruned))
    return Discovery(pruned=tuple(pruned))


def _discover_input_files(input_path: Path, pattern: str) -> list[Path]:
    """Return only the files one batch target offers."""

    return list(_discover_science_files(input_path, pattern).files)


def _display_path(path: Path, root: Path) -> str:
    """Name a folder against the source root the header printed one line above."""

    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def _pruned_lines(pruned: tuple[Path, ...], root: Path) -> list[str]:
    """Say what discovery removed, whenever it removed anything.

    The bench writes a snapshot into whichever folder it was launched at, and
    that can legitimately be a day folder that also holds science SIFs.  The
    whole day then disappears from the walk, and until this said so the console
    header, the receipt, and the campaign page all read as a complete run.
    """

    if not pruned:
        return []
    shown = [_display_path(path, root) for path in pruned[:PRUNED_LISTING_LIMIT]]
    lines = [
        f"Skipped:     {len(pruned)} calibration folder(s) not searched "
        "(named 'calibrations', or holding a snapshot.toml):",
        *(f"               {name}" for name in shown),
    ]
    remaining = len(pruned) - len(shown)
    if remaining:
        lines.append(f"             ... and {remaining} more")
    return lines


def _discovered_lines(files: Sequence[Path], root: Path) -> list[str]:
    """Name every file a dry run would convert.

    Naming them is the whole point of a dry run -- the documented step says to
    "inspect the listed sources before allowing writes" -- yet the listing sat
    behind ``--verbose``, so ``--dry-run`` alone printed a count and nothing to
    inspect.  The paths are relative to the source root because a recursive walk
    finds shots several day folders down, and the folder is what identifies
    them.
    """

    return [
        "Would convert:",
        *(f"  {_display_path(path, root)}" for path in files),
    ]


def _no_files_message(input_path: Path, pattern: str) -> str:
    """Explain an empty discovery, including what it deliberately skipped."""

    if not _pattern_is_recursive(pattern):
        return f"No files matching '{pattern}' found in {input_path}"
    return (
        f"No files matching '{pattern}' found in {input_path} or anywhere below it. "
        "The search was recursive and excluded calibration folders: any folder named "
        "'calibrations', and any folder holding a snapshot.toml. Lamp frames are not "
        "science shots; point --pattern at a path pattern (one containing / or **) to "
        "search exactly where you mean instead. A junction or symbolic link among the "
        "day folders is matched where it sits but is never descended into."
    )


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


def _emit_target(message: str, *, target_label: str | None = None, stream=None) -> None:
    prefix = f"[{target_label}] " if target_label else ""
    selected_stream = sys.stdout if stream is None else stream
    with _PRINT_LOCK:
        for line in message.splitlines() or [""]:
            print(f"{prefix}{line}", file=selected_stream)


def _batch_header(
    *,
    input_path: Path,
    output_dir: Path,
    n_files: int,
    pattern: str,
    dry_run: bool,
    pruned: tuple[Path, ...] = (),
    target_label: str | None = None,
) -> None:
    mode = "DRY RUN" if dry_run else "export"
    _emit_target(
        "\n".join(
            [
                _color("SpectroCube batch", "1;36"),
                f"Source:      {input_path}",
                f"Destination: {output_dir}",
                f"Pattern:     {pattern}",
                f"Files:       {n_files} ({mode})",
                *_pruned_lines(pruned, input_path),
            ]
        ),
        target_label=target_label,
    )


def _warn_unattributed_absolute(settings: dict, target_label: str | None = None) -> None:
    """Say when an absolute export is about to record no flux standard at all.

    ``calibration_source`` is the one attribute naming what the absolute scale
    was traced to, and it is written only when something supplies it.  Nothing
    used to say so, so a whole drive could be exported in W/m2/nm/sr with the
    provenance field simply absent, discovered only by reading a finished cube.
    """

    if str(settings.get("units", "")).lower() not in ABSOLUTE_UNITS:
        return
    if settings.get("calibration_source"):
        return
    _emit_target(
        f"WARNING: units '{settings['units']}' are absolute but no calibration source "
        "is named, so the cubes will not record what the absolute scale was traced "
        "to. Pass --calibration-source, or set [export] calibration_source in the "
        "export config.",
        target_label=target_label,
        stream=sys.stderr,
    )


def _settings_from_args(args: argparse.Namespace) -> tuple[argparse.Namespace, dict]:
    """Merge built-in defaults, config TOML, plan TOML, and CLI overrides."""
    plan = export_plan_from_toml(args.plan) if args.plan else {}

    config_path = args.config or plan.get("config")
    settings = dict(_DEFAULTS)
    settings["calibration_files"] = {}
    settings["extra_attrs"] = {}
    settings["_camera_explicit"] = False
    if config_path:
        config_settings = export_config_from_toml(config_path)
        settings["_camera_explicit"] = config_settings.get("camera") not in (
            None,
            "",
        )
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
            if key == "camera":
                settings["_camera_explicit"] = True

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
    if args.registry is not None:
        settings["registry"] = args.registry
    elif plan.get("registry") and not settings.get("registry"):
        settings["registry"] = plan["registry"]
    if args.calibrations is not None:
        settings["calibrations"] = args.calibrations
    elif plan.get("calibrations") and not settings.get("calibrations"):
        settings["calibrations"] = plan["calibrations"]

    if args.input is None:
        args.input = plan.get("input") or plan.get("input_dir")
    if args.output is None:
        args.output = plan.get("output") or plan.get("output_dir")
    if args.pattern == "*.SIF" and plan.get("pattern"):
        args.pattern = plan["pattern"]
    if args.drift_verdict is None and plan.get("drift_verdict"):
        args.drift_verdict = plan["drift_verdict"]
    if args.central_index is None and plan.get("central_index"):
        args.central_index = plan["central_index"]
    args.overwrite = bool(args.overwrite or plan.get("overwrite", False))
    args.dry_run = bool(args.dry_run or plan.get("dry_run", False))
    args.verbose = bool(args.verbose or plan.get("verbose", False))

    return args, settings


def _load_epoch_registry(
    args: argparse.Namespace,
    settings: dict,
    parser: argparse.ArgumentParser,
) -> CalibrationEpochRegistry | None:
    registry_path = settings.get("registry")
    if not registry_path:
        return None
    manual_files = settings.get("calibration_files") or {}
    if (
        args.snapshot_id
        or settings.get("_camera_explicit")
        or settings.get("calibration_dir")
        or manual_files
    ):
        parser.error(
            "--registry cannot be combined with --snapshot-id, --camera, "
            "--calibration-dir, or manual calibration-file overrides; the selected "
            "snapshot is the calibration authority"
        )
    try:
        return load_calibration_registry(
            registry_path,
            snapshots_root=settings.get("calibrations"),
        )
    except CalibrationRegistryError as exc:
        parser.error(str(exc))
    return None  # pragma: no cover - argparse exits


def _snapshot_camera(snapshot) -> str:
    camera = snapshot.detector.upper()
    if camera not in {"CMOS", "CCD"}:
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} detector {snapshot.detector!r} "
            "does not map to a supported camera (CMOS or CCD)"
        )
    return camera


def _registry_provenance(
    registry: CalibrationEpochRegistry,
    *,
    position: int,
) -> dict[str, object]:
    return {
        "calibration_registry_schema": "echelle-calibration-registry/v1",
        "calibration_registry_sha256": sha256_file(registry.path),
        "calibration_registry_epoch_position": int(position),
    }


@dataclass(frozen=True)
class GateAuthorization:
    """How one target earned the right to process its calibration epoch."""

    gate: str
    files: tuple[Path, ...]
    evidence_path: str = ""
    evidence_sha256: str = ""
    verdict: str = ""
    sample: bool = False


def _date_scan_root(input_path: Path) -> Path:
    """Return the declared source root that bounds a source's date scan.

    A named folder is its own root.  A named file is bounded by the folder
    holding it, so ``2024-03-05/lamp.SIF`` still dates itself while a dated
    mount above it cannot.  The gate and the export must agree on this, or a
    run could be authorized for one epoch and exported under another.
    """

    return input_path if input_path.is_dir() else input_path.parent


def _shell_quote(value: object) -> str:
    text = str(value)
    return f'"{text}"' if " " in text else text


def _gate_refusal(
    reason: str,
    args: argparse.Namespace,
    settings: dict,
    input_path: Path,
) -> str:
    """Refuse a registry run and name the exact commands that produce the evidence."""

    output = args.output or str(input_path)
    sample_command = [
        "echelle process",
        _shell_quote(input_path),
        f"--registry {_shell_quote(settings.get('registry') or 'calibration_registry.toml')}",
    ]
    if settings.get("calibrations"):
        sample_command.append(f"--calibrations {_shell_quote(settings['calibrations'])}")
    sample_command.extend([f"-o {_shell_quote(output)}", "--sample 5"])
    return "\n".join(
        [
            f"ERROR: {reason}",
            "  Take the sampled evidence this gate needs:",
            f"    {' '.join(sample_command)}",
            # The audit takes the cube directory itself, so this line is the
            # same on PowerShell and on a POSIX shell: no glob is expanded.
            f"    echelle drift audit {_shell_quote(output)} -o drift-evidence.json",
            "  Then repeat this command with --drift-verdict drift-evidence.json",
        ]
    )


def _authorize_run(
    args: argparse.Namespace,
    settings: dict,
    registry: CalibrationEpochRegistry | None,
    sif_files: list[Path],
    *,
    input_path: Path,
    target_label: str | None = None,
) -> GateAuthorization | None:
    """Authorize, sample, or refuse one target; refusals are already reported."""

    if registry is None:
        # An explicitly configured calibration is legal, but the receipt must
        # say forever that nothing verified this epoch.
        return GateAuthorization(gate=GATE_UNGATED, files=tuple(sif_files))

    from .drift import DriftError, require_sampled_verdict

    if args.sample is not None:
        size = args.sample
        if size == SAMPLE_AUTO:
            size = _auto_sample_count(len(sif_files))
            _emit_target(
                f"sample auto: {size} of {len(sif_files)} files",
                target_label=target_label,
            )
        selected = tuple(sif_files[:size])
        _emit_target(
            f"Gate:        sample of {len(selected)}/{len(sif_files)} file(s) with no "
            "verdict; every cube is marked drift_sample",
            target_label=target_label,
        )
        return GateAuthorization(gate=GATE_SAMPLE, files=selected, sample=True)

    if not args.drift_verdict:
        _emit_target(
            _gate_refusal(
                "registry-backed processing requires a sampled epoch audit "
                "(--drift-verdict) or an explicit unverified first sample (--sample N).",
                args,
                settings,
                input_path,
            ),
            target_label=target_label,
            stream=sys.stderr,
        )
        return None

    try:
        selected_ids = {
            registry.resolve_source(path, root=_date_scan_root(input_path)).snapshot_id
            for path in sif_files
        }
        payload = require_sampled_verdict(args.drift_verdict, selected_ids)
    except (CalibrationRegistryError, DriftError, OSError, ValueError) as exc:
        _emit_target(
            _gate_refusal(f"drift gate failed: {exc}", args, settings, input_path),
            target_label=target_label,
            stream=sys.stderr,
        )
        return None

    evidence = Path(args.drift_verdict)
    _emit_target(
        f"Gate:        verdict '{payload.get('verdict')}' from {evidence}",
        target_label=target_label,
    )
    return GateAuthorization(
        gate=GATE_VERDICT,
        files=tuple(sif_files),
        evidence_path=str(evidence.resolve()),
        evidence_sha256=sha256_file(evidence),
        verdict=str(payload.get("verdict", "")),
    )


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
    snapshot: object | None = None,
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
            snapshot=snapshot,
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


def _resume_message(receipt: RunReceipt, target_label: str | None) -> str:
    if target_label:
        return f"Interrupted safely. Rerun the same campaign command to resume {receipt.directory}"
    return f"Interrupted safely. Resume with --run-dir {receipt.directory}"


def _write_drive_catalog(
    args: argparse.Namespace,
    out_dir: Path,
    receipt: RunReceipt | None,
    target_label: str | None,
) -> None:
    """Write the per-drive catalog and, when configured, fold it into the index.

    Only a completed target is auto-merged: an interrupted or partial run's
    catalog is written for the operator but is not published into the durable
    all-years index as though the drive were finished.
    """

    if receipt is None:
        return
    from .catalog import build_drive_catalog, merge_into_central_index

    catalog = build_drive_catalog(
        out_dir,
        volume_label=receipt.volume_label,
        # An empty id means this run never saw one; let the catalog look for a
        # drive that announced itself around the cubes instead of asserting none.
        drive_id=receipt.drive_id or None,
        receipt_dir=receipt.directory,
    )
    _emit_target(f"Catalog:     {catalog}", target_label=target_label)
    if not args.central_index or receipt.state != "completed":
        return
    try:
        merged = merge_into_central_index(catalog, args.central_index)
    except (OSError, ValueError) as exc:
        _emit_target(
            f"WARNING: Central index was not updated: {exc}",
            target_label=target_label,
            stream=sys.stderr,
        )
        return
    _emit_target(f"Index:       {merged}", target_label=target_label)


def _run_batch_target(
    args: argparse.Namespace,
    settings: dict,
    input_path: Path,
    out_dir: Path,
    *,
    volume_label: str | None = None,
    runs_root: Path | None = None,
    target_label: str | None = None,
    stop_event: threading.Event | None = None,
) -> int:
    """Process one source sequentially and return its independent exit code."""
    registry: CalibrationEpochRegistry | None = settings.get("epoch_registry")
    discovery = _discover_science_files(input_path, args.pattern)
    sif_files = list(discovery.files)
    if not sif_files:
        _emit_target(
            _no_files_message(input_path, args.pattern),
            target_label=target_label,
            stream=sys.stderr,
        )
        return 1

    authorization = _authorize_run(
        args,
        settings,
        registry,
        sif_files,
        input_path=input_path,
        target_label=target_label,
    )
    if authorization is None:
        return 1
    sif_files = list(authorization.files)

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
    _batch_header(
        input_path=input_path,
        output_dir=out_dir,
        n_files=len(sif_files),
        pattern=args.pattern,
        dry_run=args.dry_run,
        pruned=discovery.pruned,
        target_label=target_label,
    )
    if args.dry_run:
        _emit_target(
            "\n".join(_discovered_lines(sif_files, input_path)), target_label=target_label
        )
    _warn_unattributed_absolute(settings, target_label)

    receipt = None
    if not args.dry_run:
        receipt_root = runs_root or Path(args.runs_dir)
        # First processing of this target announces the drive's stable id at its
        # root, so a later reconnection under another letter or mount point stays
        # one drive in every catalog. A resumed receipt written before the id
        # existed picks it up here rather than staying anonymous forever.
        identity = ensure_drive_identity(input_path, label=volume_label)
        if identity.warning:
            _emit_target(
                f"WARNING: {identity.warning}",
                target_label=target_label,
                stream=sys.stderr,
            )
        # The receipt folder is named for the drive as well as the folder, so
        # two drives whose shots both live in `\data` do not produce two
        # receipts that read alike.
        receipt_name_label = volume_label or identity.label
        if args.run_dir:
            receipt_dir = Path(args.run_dir)
        elif not args.new_run:
            receipt_dir = find_resumable_run(receipt_root, input_path, out_dir, args.pattern)
            if receipt_dir is None:
                receipt_dir = new_run_directory(
                    receipt_root, input_path, volume_label=receipt_name_label
                )
        else:
            receipt_dir = new_run_directory(
                receipt_root, input_path, volume_label=receipt_name_label
            )

        resuming = (receipt_dir / "run.toml").is_file()
        if resuming:
            receipt = RunReceipt.load(receipt_dir)
            if receipt.source_root.resolve() != input_path.resolve():
                _emit_target(
                    "ERROR: Run receipt source does not match this input.",
                    target_label=target_label,
                    stream=sys.stderr,
                )
                return 1
            if receipt.output_root.resolve() != out_dir.resolve():
                _emit_target(
                    "ERROR: Run receipt destination does not match this output.",
                    target_label=target_label,
                    stream=sys.stderr,
                )
                return 1
            if receipt.pattern != args.pattern:
                _emit_target(
                    "ERROR: Run receipt pattern does not match --pattern.",
                    target_label=target_label,
                    stream=sys.stderr,
                )
                return 1
        else:
            receipt_fields = {
                "source_root": input_path,
                "output_root": out_dir,
                "pattern": args.pattern,
                "volume_label": volume_label or identity.label,
                "snapshot_id": args.snapshot_id
                or ("per-source-registry" if registry is not None else "unassigned"),
                "expected_files": len(sif_files),
                "drive_id": identity.drive_id,
                "drive_warning": identity.warning,
            }
            try:
                receipt = RunReceipt.create(receipt_dir, **receipt_fields)
            except FileExistsError:
                receipt_dir = new_run_directory(receipt_root, input_path)
                receipt = RunReceipt.create(receipt_dir, **receipt_fields)
            _emit_target(
                f"Receipt:     {receipt.directory}", target_label=target_label
            )
        expected_receipt_snapshot = args.snapshot_id or (
            "per-source-registry" if registry is not None else None
        )
        if expected_receipt_snapshot and receipt.snapshot_id != expected_receipt_snapshot:
            _emit_target(
                "ERROR: Run receipt snapshot does not match --snapshot-id.",
                target_label=target_label,
                stream=sys.stderr,
            )
            return 1
        if volume_label and receipt.volume_label != volume_label:
            _emit_target(
                "ERROR: Run receipt volume does not match --volume-label.",
                target_label=target_label,
                stream=sys.stderr,
            )
            return 1
        if resuming:
            receipt.state = "running"
            receipt.expected_files = max(receipt.expected_files, len(sif_files))
            receipt.drive_id = receipt.drive_id or identity.drive_id
            receipt.drive_warning = identity.warning
            receipt.write_manifest()
            _emit_target(f"Resuming:    {receipt.directory}", target_label=target_label)
        # The receipt is what an auditor and the campaign page read long after
        # the console has scrolled away, so it carries the same skip the header
        # announced -- in full, not the bounded listing.
        receipt.record_discovery(
            _display_path(path, input_path) for path in discovery.pruned
        )
        # A resumed receipt records how the run in front of it was authorized,
        # not how an earlier one was.
        receipt.record_authorization(
            gate=authorization.gate,
            sample=authorization.sample,
            sample_files=len(sif_files) if authorization.sample else 0,
            evidence_path=authorization.evidence_path,
            evidence_sha256=authorization.evidence_sha256,
            verdict=authorization.verdict,
        )

    if stop_event is not None and stop_event.is_set():
        if receipt is not None:
            receipt.finish("interrupted")
            _write_drive_catalog(args, out_dir, receipt, target_label)
        return 130

    cal_dir = Path(settings["calibration_dir"]) if settings.get("calibration_dir") else None
    calibration_files = settings["calibration_files"]
    clbr = None
    if not args.dry_run and registry is None:
        try:
            from .tools.loader import build_calibration

            if args.verbose:
                _emit_target(
                    f"Loading {settings['camera']} calibration...", target_label=target_label
                )
            clbr = build_calibration(
                cal_dir,
                settings["camera"],
                calibration_files=calibration_files,
            )
            if args.verbose:
                _emit_target(_color("Calibration ready.", "32"), target_label=target_label)
        except KeyboardInterrupt:
            if stop_event is not None:
                stop_event.set()
            if receipt is not None:
                receipt.finish("interrupted")
                _write_drive_catalog(args, out_dir, receipt, target_label)
                _emit_target(
                    _resume_message(receipt, target_label),
                    target_label=target_label,
                )
            return 130
        except Exception as exc:
            _emit_target(
                f"ERROR: Could not load calibration: {exc}",
                target_label=target_label,
                stream=sys.stderr,
            )
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
                _write_drive_catalog(args, out_dir, receipt, target_label)
            return 1

    calibration_cache: dict[str, object] = {}

    failed: list[Path] = []
    n_exported = 0
    n_skipped = 0
    n_resumed = 0
    n_dry_run = 0
    total = len(sif_files)
    batch_started = time.monotonic()
    for index, sif in enumerate(sif_files, start=1):
        if stop_event is not None and stop_event.is_set():
            if receipt is not None:
                receipt.finish("interrupted")
                _write_drive_catalog(args, out_dir, receipt, target_label)
                _emit_target(
                    _resume_message(receipt, target_label),
                    target_label=target_label,
                )
            return 130
        nc_out = _output_path_with_suffix(sif, out_dir, settings.get("output_suffix"))
        if args.verbose:
            _emit_target(
                _color(f"[{index}/{total}] {sif.name}", "36"), target_label=target_label
            )
        item_started = time.monotonic()
        started_at = utc_now()
        try:
            source = receipt.identity_for(sif) if receipt is not None else None
        except KeyboardInterrupt:
            if stop_event is not None:
                stop_event.set()
            if receipt is not None:
                receipt.finish("interrupted")
                _write_drive_catalog(args, out_dir, receipt, target_label)
                _emit_target(
                    _resume_message(receipt, target_label),
                    target_label=target_label,
                )
            return 130
        result: ExportResult | None = None
        selected_snapshot = None
        selected_snapshot_id: str | None = None
        selected_calibration = clbr
        selected_extra_attrs = dict(settings["extra_attrs"])
        if authorization.sample:
            selected_extra_attrs[DRIFT_SAMPLE_ATTR] = DRIFT_SAMPLE_TRUE
        if registry is not None:
            try:
                epoch = registry.resolve_source(sif, root=_date_scan_root(input_path))
                selected_snapshot = epoch.snapshot
                selected_snapshot_id = epoch.snapshot_id
                selected_extra_attrs.update(
                    _registry_provenance(registry, position=epoch.position)
                )
                if not args.dry_run:
                    selected_calibration = calibration_cache.get(epoch.snapshot_id)
                    if selected_calibration is None:
                        from .tools.loader import build_calibration

                        selected_calibration = build_calibration(
                            epoch.snapshot.root,
                            _snapshot_camera(epoch.snapshot),
                            calibration_files=epoch.snapshot.calibration_files(),
                        )
                        calibration_cache[epoch.snapshot_id] = selected_calibration
            except (CalibrationRegistryError, OSError, ValueError) as exc:
                result = ExportResult("failed", f"calibration epoch selection failed: {exc}")
                _emit_target(
                    f"FAIL {sif.name}: {result.reason}",
                    target_label=target_label,
                    stream=sys.stderr,
                )

        resumed_record = (
            receipt.resumable_record(source, nc_out, snapshot_id=selected_snapshot_id)
            if result is None
            and receipt is not None
            and source is not None
            and not args.overwrite
            else None
        )
        if resumed_record is not None:
            # The earlier record proved this source's digest; the resumed record
            # carries that proof forward rather than re-reading the file to
            # restate it.
            source = source.with_recorded_sha256(str(resumed_record.get("source_sha256", "")))
            n_resumed += 1
            receipt.record_resume_trust(RESUME_TRUST_STAT)
            result = ExportResult(
                "skipped", "resumed on recorded identity; digests not re-proved"
            )
        elif result is None and (
            receipt is not None
            and source is not None
            and not args.overwrite
            and nc_out.exists()
            and receipt.has_export_record(source, snapshot_id=selected_snapshot_id)
            # Only the output is in question here. A source that no longer
            # matches its record is a different file, and re-converting it is
            # the answer to that, not a complaint about the cube beside it.
            and receipt.source_matches_record(source, snapshot_id=selected_snapshot_id)
        ):
            result = ExportResult(
                "failed",
                "recorded completed output changed; inspect it or use --overwrite",
            )
        elif result is None:
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
                        extra_attrs=selected_extra_attrs,
                        overwrite=args.overwrite,
                        dry_run=args.dry_run,
                        verbose=args.verbose,
                        calibration=selected_calibration,
                        snapshot=selected_snapshot,
                    )
                )
            except KeyboardInterrupt:
                if stop_event is not None:
                    stop_event.set()
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
                    _write_drive_catalog(args, out_dir, receipt, target_label)
                    _emit_target(
                        _resume_message(receipt, target_label),
                        target_label=target_label,
                    )
                return 130

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
                snapshot_id=selected_snapshot_id,
            )
        if status == "failed":
            failed.append(sif)
        elif status == "skipped":
            n_skipped += 1
        elif status == "dry-run":
            n_dry_run += 1
        else:
            n_exported += 1
        _emit_target(
            _progress_line(index, total, batch_started, status), target_label=target_label
        )

    # One plain line, once, so a fast resume never reads as a verified one.
    if n_resumed:
        _emit_target(
            f"Resumed {n_resumed} file(s) on recorded identity (size+mtime); "
            "'echelle catalog verify' re-proves digests on demand.",
            target_label=target_label,
        )

    n_ok = len(sif_files) - len(failed)
    if failed:
        if receipt is not None:
            receipt.finish("partial")
            _write_drive_catalog(args, out_dir, receipt, target_label)
        failure_lines = [f"{len(failed)} failure(s):", *(f"  {path}" for path in failed)]
        failure_lines.append(f"{n_ok}/{len(sif_files)} exported successfully.")
        _emit_target(
            "\n".join(failure_lines), target_label=target_label, stream=sys.stderr
        )
        return 1

    if args.dry_run:
        _emit_target(
            _color(f"Dry run complete. {n_dry_run}/{total} file(s) would be converted.", "36"),
            target_label=target_label,
        )
    else:
        if receipt is not None:
            receipt.finish("completed")
            _write_drive_catalog(args, out_dir, receipt, target_label)
        skipped = f", {n_skipped} skipped" if n_skipped else ""
        _emit_target(
            _color(f"Done. {n_exported}/{total} exported successfully{skipped}.", "32"),
            target_label=target_label,
        )
    return 0


def _volume_labels(args: argparse.Namespace, count: int, parser: argparse.ArgumentParser) -> list[str | None]:
    labels = list(args.volume_label or [])
    if not labels:
        return [None] * count
    if len(labels) != count:
        parser.error("--volume-label must be repeated once per INPUT when processing several targets.")
    return labels


def _multi_output_dirs(inputs: list[Path], output: str | None) -> list[Path]:
    if output is None:
        return inputs
    root = Path(output)
    names = [source.name or source.anchor.rstrip("\\/").replace(":", "") or "drive" for source in inputs]
    outputs = []
    for source, name in zip(inputs, names):
        target_name = target_runs_root(Path(), source).name if names.count(name) > 1 else name
        outputs.append(root / target_name)
    return outputs


def _run_multi_targets(
    args: argparse.Namespace,
    settings: dict,
    inputs: list[Path],
    labels: list[str | None],
) -> int:
    outputs = _multi_output_dirs(inputs, args.output)
    shared_runs_root = Path(args.runs_dir)
    # The drive each target declares itself to be. It names that target's
    # receipt tree, so two drives whose shots both sit in `\data` no longer get
    # two receipt roots that read alike; the console label may be decorated
    # further below, and must not follow it onto disk.
    drive_labels = [label or default_volume_label(source) for source, label in zip(inputs, labels)]
    display_labels = list(drive_labels)
    if len(set(display_labels)) != len(display_labels):
        display_labels = [
            f"{label}/{target_runs_root(Path(), source).name or index}"
            for index, (label, source) in enumerate(zip(display_labels, inputs), start=1)
        ]

    stop_event = threading.Event()
    executor = ThreadPoolExecutor(max_workers=len(inputs), thread_name_prefix="echelle-drive")
    futures = {
        executor.submit(
            _run_batch_target,
            args,
            settings,
            source,
            output,
            volume_label=label,
            runs_root=target_runs_root(shared_runs_root, source, volume_label=drive),
            target_label=display,
            stop_event=stop_event,
        ): display
        for source, output, label, drive, display in zip(
            inputs, outputs, labels, drive_labels, display_labels
        )
    }
    codes: list[int] = []
    interrupted = False
    try:
        for future in as_completed(futures):
            try:
                code = future.result()
            except Exception as exc:
                _emit_target(
                    f"ERROR: Worker failed unexpectedly: {exc}",
                    target_label=futures[future],
                    stream=sys.stderr,
                )
                code = 1
            codes.append(code)
            if code == 130:
                stop_event.set()
                interrupted = True
    except KeyboardInterrupt:
        interrupted = True
        stop_event.set()
        _emit_target("Interrupt received; finishing active files safely.", stream=sys.stderr)
    finally:
        executor.shutdown(wait=True, cancel_futures=False)

    if interrupted:
        return 130
    succeeded = sum(code == 0 for code in codes)
    failed = len(codes) - succeeded
    _emit_target(f"Campaign complete: {succeeded} target(s) completed, {failed} failed.")
    return 1 if failed else 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def _run_single_file(
    args: argparse.Namespace,
    settings: dict,
    input_path: Path,
) -> int:
    """Run the established one-file export path."""
    if args.output:
        output = Path(args.output)
        if output.is_dir():
            output = _output_path_with_suffix(
                input_path, output, settings.get("output_suffix")
            )
    else:
        output = _output_path_with_suffix(
            input_path, input_path.parent, settings.get("output_suffix")
        )

    _warn_unattributed_absolute(settings)
    cal_dir = Path(settings["calibration_dir"]) if settings.get("calibration_dir") else None
    registry: CalibrationEpochRegistry | None = settings.get("epoch_registry")
    snapshot = None
    calibration = None
    extra_attrs = dict(settings["extra_attrs"])
    if registry is not None:
        authorization = _authorize_run(
            args, settings, registry, [input_path], input_path=input_path
        )
        if authorization is None:
            return 1
        if authorization.sample:
            extra_attrs[DRIFT_SAMPLE_ATTR] = DRIFT_SAMPLE_TRUE
        try:
            epoch = registry.resolve_source(input_path, root=_date_scan_root(input_path))
            snapshot = epoch.snapshot
            extra_attrs.update(_registry_provenance(registry, position=epoch.position))
            if not args.dry_run:
                from .tools.loader import build_calibration

                calibration = build_calibration(
                    snapshot.root,
                    _snapshot_camera(snapshot),
                    calibration_files=snapshot.calibration_files(),
                )
        except (CalibrationRegistryError, OSError, ValueError) as exc:
            print(f"ERROR: Calibration epoch selection failed: {exc}", file=sys.stderr)
            return 1
    result = _normalize_export_result(
        _export_one(
            input_path,
            output,
            units=settings["units"],
            camera=settings["camera"],
            calibration_dir=cal_dir,
            calibration_files=settings["calibration_files"],
            instrument_id=settings["instrument_id"],
            wavelength_medium=settings["wavelength_medium"],
            wavelength_min_nm=settings.get("wavelength_min_nm"),
            wavelength_max_nm=settings.get("wavelength_max_nm"),
            calibration_source=settings.get("calibration_source"),
            drop_nonfinite_columns=settings["drop_nonfinite_columns"],
            extra_attrs=extra_attrs,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            verbose=True,
            calibration=calibration,
            snapshot=snapshot,
        )
    )
    return 0 if result.status != "failed" else 1


def _spectrocube_is_available(*, dry_run: bool) -> bool:
    if dry_run:
        return True
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
        return False
    return True


def main(argv: list[str] | None = None, *, prog: str = "echelle-spectrocube") -> None:
    """Export one file or run isolated sequential workers across source folders."""
    parser = _build_parser(prog=prog)
    args = parser.parse_args(argv)
    args, settings = _settings_from_args(args)
    settings["epoch_registry"] = _load_epoch_registry(args, settings, parser)
    if args.sample is not None:
        if args.drift_verdict:
            parser.error(
                "--sample takes an unverified first sample and --drift-verdict authorizes a "
                "verified run; they cannot be combined. Sample first, audit the sample with "
                "'echelle drift audit', then rerun with --drift-verdict alone."
            )
        if args.sample != SAMPLE_AUTO and args.sample < 1:
            parser.error("--sample must select at least one file.")
        if settings.get("epoch_registry") is None:
            parser.error(
                "--sample marks cubes and receipts as unverified epoch samples and "
                "therefore requires --registry."
            )
    if args.input is None:
        parser.error("INPUT is required unless supplied by --plan.")

    inputs = [Path(args.input), *(Path(value) for value in args.additional_inputs)]
    labels = _volume_labels(args, len(inputs), parser)
    normalized_inputs = [os.path.normcase(str(path.resolve())) for path in inputs]
    if len(set(normalized_inputs)) != len(normalized_inputs):
        parser.error("Each source folder may appear only once in a campaign command.")
    if not _spectrocube_is_available(dry_run=args.dry_run):
        sys.exit(1)
    folder_inputs = [path for path in inputs if path.is_dir()]
    if folder_inputs:
        if len(folder_inputs) != len(inputs):
            parser.error("Several INPUT values must all be source folders; files run one at a time.")
        if len(inputs) > 1 and args.run_dir:
            parser.error("--run-dir names one receipt and cannot be used with several INPUT folders.")
        if len(inputs) > 1 and args.plan:
            parser.error("--plan supplies one target and cannot be combined with several INPUT folders.")
        if args.output and Path(args.output).is_file():
            parser.error("--output must be a directory when processing folders.")
        if len(inputs) > 1:
            sys.exit(_run_multi_targets(args, settings, inputs, labels))
        input_path = inputs[0]
        out_dir = Path(args.output) if args.output else input_path
        sys.exit(
            _run_batch_target(
                args,
                settings,
                input_path,
                out_dir,
                volume_label=labels[0],
            )
        )

    if len(inputs) > 1:
        missing = ", ".join(str(path) for path in inputs if not path.exists())
        if missing:
            print(f"ERROR: Input path not found: {missing}", file=sys.stderr)
            sys.exit(1)
        parser.error("Several INPUT values must be source folders; files run one at a time.")
    input_path = inputs[0]
    if input_path.is_file():
        sys.exit(_run_single_file(args, settings, input_path))
    print(f"ERROR: Input path not found: {input_path}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
