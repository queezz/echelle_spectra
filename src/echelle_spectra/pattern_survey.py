"""Where does this drive's light sit, against which era's order pattern?

Two field questions turned out to be one measurement.  The first is a diary:
sample a campaign through time, read one plasma-bright frame per sampled day,
and record where its order bands sit against each candidate era pattern -- the
geometry between calibrations, which nothing else records.  The second is a
gate: the owner's rule (decision 2026-08-27) that a bulk run over a failing
wavelength audit is permitted only when the extraction *geometry* is
separately verified, inside a 1.5-row alarm, against the epoch's own pattern.

Wavelength and flux stay revisable on a saved cube; extraction geometry does
not, because the wrong pattern sums the wrong rows and the light that was not
summed is gone.  That is why this is the one check that gates conversion, and
why the same sampling answers both questions: the diary is the evidence, and
its last row against the chosen pattern is the verdict.

Nothing here converts anything and nothing here writes to the drive.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np

SURVEY_SCHEMA = "echelle-pattern-survey/v1"

#: Past this many rows the extraction is reading the wrong part of the band,
#: and a bulk run is not authorised.
ALARM_ROWS = 1.5

#: Sample one day in six by default: enough rows to see a drift develop over a
#: cycle, few enough that the survey finishes while an operator waits.
DEFAULT_EVERY = 6

#: The eras the package ships a pattern for, newest last -- the verdict is
#: taken against the last pattern given, so the default verdict is the current
#: era's, which is the one a fresh drive is about to be converted under.
PACKAGED_PATTERNS = (("2024", "pattern_CMOS_20240305.txt"), ("2025", "pattern_CMOS_20250926.txt"))

OK = "GEOMETRY_OK"
ALARM = "GEOMETRY_ALARM"

#: Exit codes the caretaker's PowerShell already branches on.
EXIT_OK = 0
EXIT_ALARM = 2
EXIT_UNMEASURED = 3


class PatternSurveyError(Exception):
    """One refusal about the inputs a survey was handed."""


@dataclass(frozen=True)
class PatternChoice:
    """One candidate era pattern, named the way its column is headed."""

    key: str
    source: str
    rows: np.ndarray


def _load_pattern(path: Path) -> np.ndarray:
    try:
        rows = np.loadtxt(str(path), dtype=int)
    except (OSError, ValueError) as exc:
        raise PatternSurveyError(f"pattern file could not be read: {path} -- {exc}") from None
    rows = np.atleast_2d(rows)
    if rows.ndim != 2 or not rows.size:
        raise PatternSurveyError(
            f"pattern file is not a table of detector columns by orders: {path}"
        )
    return rows


def packaged_patterns() -> list[PatternChoice]:
    """The era patterns this package ships, in the order they were taken."""

    calibrations = files("echelle_spectra.resources").joinpath("calibration_files")
    return [
        PatternChoice(key, f"packaged {name}", _load_pattern(Path(str(calibrations.joinpath(name)))))
        for key, name in PACKAGED_PATTERNS
    ]


def load_patterns(paths: Sequence[str | Path] | None) -> list[PatternChoice]:
    """Return the patterns to measure against: the given ones, else the packaged."""

    if not paths:
        return packaged_patterns()
    chosen = []
    for raw in paths:
        path = Path(raw)
        chosen.append(PatternChoice(path.name, str(path), _load_pattern(path)))
    return chosen


def day_folders(root: Path) -> tuple[Path, list[Path]]:
    """Return the folder the days sit in and the day folders themselves.

    A drive carries its observing days under ``DATA``; a folder handed over on
    its own is already that place.  Either way a day folder is named by its
    date and a calibration folder is not a day.
    """

    root = Path(root)
    base = root / "DATA" if (root / "DATA").is_dir() else root
    days = sorted(
        child
        for child in base.iterdir()
        if child.is_dir()
        and not child.name.startswith(".")
        and child.name[:4].isdigit()
        and "calib" not in child.name.lower()
    )
    if not days and shot_files(base):
        # The root is itself one day's folder, which is how a single day gets
        # checked the moment it comes off the camera.
        return base, [base]
    return base, days


def shot_files(day: Path) -> list[Path]:
    """The frames in one day folder, in shot order, metadata siblings left out."""

    frames = sorted(day.glob("*.SIF")) or sorted(day.glob("*.sif"))
    return [path for path in frames if not path.name.startswith(".")]


def sample_days(days: Sequence[Path], every: int) -> list[Path]:
    """Choose the days to read: every Nth, or the three-frame spread.

    ``every`` of zero -- and a root holding a single day -- ask the older
    question, the one the gate was born as: first day, middle day, last day.
    Three frames spread across a campaign is the cheapest reading that can
    still catch a geometry that moved partway through it.
    """

    days = list(days)
    if not days:
        return []
    if every <= 0 or len(days) == 1:
        spread = [days[0], days[len(days) // 2], days[-1]]
        return list(dict.fromkeys(spread))
    return days[::every] or days


def brightest_frame(images: np.ndarray) -> np.ndarray:
    """The frame carrying the most light, which is the one worth measuring."""

    stack = np.asarray(images, dtype=float)
    if stack.ndim == 2:
        return stack
    return stack[int(np.argmax(np.nansum(stack, axis=(1, 2))))]


class BenchFrameReader:
    """Read one SIF's frames through the bench's own loader."""

    def __init__(self, pattern: np.ndarray) -> None:
        from .calibration_bench import FrameLoader

        self._loader = FrameLoader(np.asarray(pattern, dtype=int))

    def __call__(self, path: Path) -> np.ndarray:
        return np.asarray(self._loader(path).images, dtype=float)


def frame_reader(pattern: np.ndarray) -> Callable[[Path], np.ndarray]:
    """The reader a survey uses when the caller supplies none."""

    return BenchFrameReader(pattern)


def _measure(image: np.ndarray, patterns: Sequence[PatternChoice]) -> tuple[dict, dict]:
    """Median and worst per-order offset of one frame against each pattern."""

    from .calibration_bench import band_center_offsets

    medians: dict[str, float | None] = {}
    worst: dict[str, float | None] = {}
    for choice in patterns:
        reading = band_center_offsets(image, choice.rows)
        offsets = [order.offset_rows for order in reading.orders if order.offset_rows is not None]
        if not reading.measured or not offsets:
            medians[choice.key] = None
            worst[choice.key] = None
            continue
        medians[choice.key] = round(float(np.median(offsets)), 3)
        worst[choice.key] = round(float(max(abs(value) for value in offsets)), 3)
    return medians, worst


def survey_geometry(
    root: Path,
    *,
    patterns: Sequence[PatternChoice],
    every: int = DEFAULT_EVERY,
    reader: Callable[[Path], np.ndarray] | None = None,
    report: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Read one frame per sampled day and record where its bands sit.

    ``report`` receives the table as it is measured rather than at the end: on
    a real drive this runs for minutes, and an operator watching a row appear
    per day can tell a slow survey from a stuck one.
    """

    if not patterns:  # pragma: no cover - the loader always returns at least one
        raise PatternSurveyError("no pattern to measure against")
    say = report or (lambda line: None)
    root = Path(root)
    base, days = day_folders(root)
    sampled = sample_days(days, every)
    # The loader wants a pattern only to check that it describes this detector;
    # the band reading takes each pattern in turn off the raw frame. Hand it the
    # verdict pattern, the last one given: that is the era the drive is about to
    # be converted under, so a stale candidate of the wrong width costs its own
    # column and not every frame on the drive.
    read_frames = reader or frame_reader(patterns[-1].rows)

    widths = {choice.key: max(12, len(choice.key) + 2) for choice in patterns}
    say("day        shot    " + "".join(f"{c.key:>{widths[c.key]}}" for c in patterns))

    rows: list[dict[str, Any]] = []
    for day in sampled:
        shots = shot_files(day)
        if not shots:
            continue
        pick = shots[len(shots) // 2]
        label = pick.stem.split("_")[0]
        try:
            images = read_frames(pick)
        except Exception as exc:  # noqa: BLE001 - one bad frame is a row, not the end
            # A frame the SIF reader will not open is a fact about the drive.
            # Refusing the whole survey over it would lose every day after it.
            say(f"{day.name}  {label}  unreadable: {exc}")
            rows.append({"day": day.name, "shot": pick.name, "error": str(exc)})
            continue
        medians, worst = _measure(brightest_frame(images), patterns)
        rows.append(
            {
                "day": day.name,
                "shot": pick.name,
                "offsets_rows": medians,
                "worst_rows": worst,
            }
        )
        cells = "".join(
            f"{'--':>{widths[c.key]}}"
            if medians[c.key] is None
            else f"{medians[c.key]:+{widths[c.key]}.2f}"
            for c in patterns
        )
        say(f"{day.name}  {label:>6}{cells}")

    verdict_key = patterns[-1].key
    measured = [
        row["worst_rows"][verdict_key]
        for row in rows
        if "worst_rows" in row and row["worst_rows"][verdict_key] is not None
    ]
    worst_offset = max(measured) if measured else None
    return {
        "schema": SURVEY_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "drive_root": str(root),
        "data_root": str(base),
        "every_nth_day": int(every),
        "day_folders": len(days),
        "patterns": {choice.key: choice.source for choice in patterns},
        "verdict_pattern": verdict_key,
        "alarm_rows": ALARM_ROWS,
        "worst_offset_rows": worst_offset,
        "verdict": OK if worst_offset is not None and worst_offset <= ALARM_ROWS else ALARM,
        "rows": rows,
    }


def verdict_lines(payload: dict[str, Any]) -> tuple[list[str], int]:
    """The closing sentences and the exit code the caretaker branches on.

    The verdict word is always the last line, and nothing measured is never
    reported as nothing wrong: an unreadable campaign leaves with its own code
    so a script cannot mistake silence for permission.
    """

    base = payload["data_root"]
    if not payload["day_folders"]:
        return [f"no day folders under {base}", ALARM], EXIT_UNMEASURED
    if not payload["rows"]:
        return [f"no SIF files found under {base}", ALARM], EXIT_UNMEASURED
    worst = payload["worst_offset_rows"]
    if worst is None:
        return [
            "no frame could be measured (dim shots?) -- geometry UNKNOWN",
            ALARM,
        ], EXIT_UNMEASURED
    against = payload["verdict_pattern"]
    if worst <= payload["alarm_rows"]:
        return [
            f"worst offset {worst:.2f} rows vs {against}, "
            f"inside the {payload['alarm_rows']}-row alarm",
            OK,
        ], EXIT_OK
    return [
        f"worst offset {worst:.2f} rows vs {against}, "
        f"PAST the {payload['alarm_rows']}-row alarm",
        ALARM,
    ], EXIT_ALARM


def surveyed_name(root: Path) -> str:
    """The drive this survey belongs to, whichever folder was handed over.

    The file is named after the drive so a campaign home accumulates one
    survey per volume beside its inventory -- and a root typed as the drive's
    ``DATA`` folder must not file itself under the name ``DATA``, which every
    drive shares.
    """

    root = Path(root)
    return root.parent.name if root.name.upper() == "DATA" and root.parent.name else root.name


def survey_path(payload: dict[str, Any], output_root: Path) -> Path:
    """Where this drive's geometry diary belongs under an output folder."""

    name = surveyed_name(Path(payload["drive_root"]))
    return Path(output_root) / "inventory" / f"{name}-pattern-survey.json"


def write_survey(payload: dict[str, Any], output_root: Path) -> Path:
    """Write the survey beside the inventory and return its path."""

    path = survey_path(payload, output_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return path
