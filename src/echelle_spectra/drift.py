"""Science-line drift sampling, verdicts, and immutable refinements.

The judge works in detector space.  A rigid detector shift moves every line by
the same number of *pixels*, never by the same number of nanometres: across the
audited Balmer/Fulcher lines the dispersion runs from about 0.0066 nm/px in the
blue orders to about 0.0108 nm/px at H-alpha, a ~60% spread.  Each measured
wavelength residual is therefore divided by the dispersion that the cube's own
stored per-order polynomial has at that line's detector pixel, one rigid shift
is fitted in pixels, and every verdict boundary is a pixel bound.

The emitted correction is per-order for the same reason.  An accepted shift
slides the wavelength table's anchors along the detector; each order's refitted
polynomial then turns those pixels into that order's own wavelength correction,
which is what ``echelle recal-cube`` later applies.  Correcting every order by
one scalar wavelength instead would leave the low-dispersion orders several
pixels wrong while the next audit called them aligned.

The judge is also isotope-aware, because LHD ran deuterium from 2017 to 2022
and D-alpha sits 0.178 nm blueward of H-alpha -- well inside the +/-0.4 nm
window this audit fits.  Judged against hydrogen alone, a perfectly calibrated
deuterium shot reads as a confident ~16.5 px shift, and the refinement it
invites would be wrong by exactly one isotope.  So each Balmer window is judged
against BOTH references and the centroid is assigned to whichever is nearer in
pixels; the shot is tagged with the isotopologue its lines chose; and on a shot
that shows deuterium the H2 Fulcher anchors are dropped rather than fitted
against the wrong molecule.

The verdict is decided on the *strong* lines only.  The 2026-08-18 rehearsal on
real light settled this: a flawless era match read the Balmer residuals at
+0.38/+0.65/+0.66/+1.21 px on lines of SNR up to 440, while weak Fulcher lines
near the detection floor scattered over -10..+17 px, because a +/-0.4 nm fit
window is +/-37 px at H-alpha and +/-60 px in the blue -- wide enough to reach
a neighbour and drag a faint centroid across it.  Judged over every line at
once that shot read ``misaligned-beyond-repair``, and its advice (reprocess the
raw SIF) could not have repaired what was centroid noise.  So the strong lines
carry the verdict and the weak ones are reported beside it as a data-quality
figure with their own words; ``misaligned-beyond-repair`` is reserved for
strong-line residuals that are both large and mutually inconsistent, which is
what real geometric damage looks like.

The same rehearsal produced the isotope hazard in the flesh.  A 2019 cube
recalibrated onto the *wrong* era had its hydrogen Balmer lines re-identified
as deuterium -- H-gamma to D-gamma is 0.1187 nm, ~16.5 px, degenerate with the
era shift -- and the flipped fit then looked *better* than the correct one.  So
whenever every two-reference Balmer window in a shot chooses deuterium and the
hydrogen reading of the same centroids implies one common shift equal to the
H/D separation, the evidence carries a loud ambiguity finding naming both
readings.  When the calendar excludes deuterium for that shot's date the
finding is not advisory: the interval verdict becomes
``era-misassigned-calibration`` and no bulk run is authorized by it.  The
finding still never moves a residual or an assignment -- the spectroscopy
stands as measured; only the interval's authorization changes.

This is a guard, not a precision instrument (owner's scope, 2026-08-15: the
audit is a rough-alignment epoch gate, and fine wavelength work belongs to
local lines in analysis, outside this pipeline).  One centroid per window, no
two-component deconvolution.  The isotope offset is ~16.5 px in every audited
order -- a rigid pixel offset, since dispersion grows with wavelength just as
the offset does -- so it is exactly degenerate with a real 16.5 px detector
shift, and nearest-assignment cannot separate the two past their ~8 px
midpoint.  That is what the bundled deuterium calendar is for: it says where
the isotope question exists at all, and when the calendar and the measurement
disagree the shot is flagged.  Spectroscopy decides; the calendar tells the
operator where to look twice, and where it says hydrogen outright it withholds
the interval's authorization rather than resolving the physics.
"""

from __future__ import annotations

import json
import re
import tempfile
from datetime import date, datetime, timezone
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback.
    import tomli as tomllib  # type: ignore[no-redef]

from .campaign_run import sha256_file
from .snapshot import REQUIRED_ROLES, Snapshot, create_snapshot, load_snapshot
from .tools.calibration_alignment import (
    load_wavelength_table,
    shift_lines_in_pixels,
    write_wavelength_table,
)
from .tools.line_catalog import SpectralLine, load_line_table

DRIFT_SCHEMA = "echelle-drift-evidence/v2"
DRIFT_SCHEMA_V1 = "echelle-drift-evidence/v1"
#: Schemas the bulk-processing gate still reads.  Only v2 is ever written: it
#: carries the detector-space residuals and per-order corrections that v1 could
#: not express.  v1 records stay readable so evidence taken before this packet
#: keeps authorizing the runs it already authorized.
READABLE_DRIFT_SCHEMAS = (DRIFT_SCHEMA, DRIFT_SCHEMA_V1)

# --- Verdict boundaries, in detector pixels --------------------------------
#
# ALIGNMENT_MAX_PIXEL_RESIDUAL — the finest disagreement a sampled audit may
# call a misalignment.  The bench solves its rigid alignment to about 0.67 px
# RMS on real Neon data (the shipped 20250926 alignment settings), so half a
# pixel is already below what the calibration itself is known to.  It is also
# 0.003-0.005 nm across the audited dispersion range, an order inside the
# ~0.03 nm instrumental line width.  Because the median can never exceed the
# maximum, this one bound states alignment completely.  A cube carrying
# ``wavelength_accuracy_nm`` relaxes the bound per line to whatever that
# accuracy is worth in pixels there: the judge never claims a misalignment
# finer than the calibration's own stated accuracy.
ALIGNMENT_MAX_PIXEL_RESIDUAL = 0.5

# SHIFT_CONSISTENCY_PX — how far one *verdict-carrying* line's pixel residual
# may sit from the fitted rigid shift and still belong to the same rigid shift.
# One pixel is ~1.5x the bench's alignment RMS and covers ordinary centroid
# noise; beyond it the residuals are not one translation of the detector, so no
# single shift repairs them and the interval wants splitting -- or, when the
# disagreement is also large, is beyond repair.  The bound may stay at one
# pixel only because VERDICT_SNR_FLOOR keeps noisy centroids out of the fit:
# every line that reaches it is good to <= 0.25 px, so a pixel of disagreement
# is ~4 sigma of centroid noise and means a real disagreement.
SHIFT_CONSISTENCY_PX = 1.0

# GEOMETRIC_DAMAGE_PX — how far the worst strong-line residual must sit from
# zero before an inconsistent set is called damage rather than quality.  1.5 px
# is 1.5x SHIFT_CONSISTENCY_PX and ~2.2x the bench's own 0.67 px rigid-alignment
# RMS: inside it a slid table still lands every strong line within about a
# pixel, which is a correction worth making, and the honest report is a shift
# whose consistency is noted.  Beyond it, lines that also refuse to agree on one
# translation are not a table that slipped -- the detector geometry itself
# moved, and only re-identifying the anchors against lamp data can fix that.
GEOMETRIC_DAMAGE_PX = 1.5

# REPAIR_LIMIT_PX — the largest rigid shift a table correction may still claim.
# Two independent measurement windows bound it: this audit fits centroids in a
# +/-0.4 nm window, which is +/-37 px at H-alpha's 0.0108 nm/px (+/-60 px in
# the blue), and the bench fits its lamp centroids in a +/-18 px window.  At
# 25 px every audited line is still inside the window that measured it, with
# margin, and the corrected anchors still fall inside the pixel intervals the
# table records for them.  Beyond that the anchors must be re-identified
# against lamp data rather than slid.  25 px is also the honest pixel reading
# of the 0.25 nm limit this rule replaces: that limit was 23 px at H-alpha and
# 38 px in the blue, differing only by which order a line happened to fall in.
REPAIR_LIMIT_PX = 25.0

# --- Evidence sufficiency ---------------------------------------------------
#
# MINIMUM_CENTROIDS counts *resolved* centroids: lines the catalog does not
# flag as sub-resolution blends, measured at or above MINIMUM_SNR.  Three of
# them must come from three distinct Echelle orders and must not all sit in one
# half of the audited wavelength coverage.  A rigid shift is a one-parameter
# claim about the whole detector, and one order in one colour cannot see the
# dispersion spread that separates a real shift from a dispersion change; a
# duplicated or blended catalog row must never manufacture that quorum.
MINIMUM_CENTROIDS = 3
MINIMUM_DISTINCT_ORDERS = 3
MINIMUM_SNR = 4.0

# --- Which lines may decide a verdict ---------------------------------------
#
# INSTRUMENT_LINE_SIGMA_NM is this spectrograph's Gaussian line sigma, the ~0.03
# nm instrumental width the alignment bound is already stated against.
INSTRUMENT_LINE_SIGMA_NM = 0.035

# VERDICT_SNR_FLOOR — the SNR at which one centroid can state the alignment
# bound at all, and therefore the line between evidence and quality.
#
# The physics: a line of sigma 0.035 nm is 3.2 px wide at H-alpha's 0.0108 nm/px
# and 5.1 px wide in the bluest audited order (0.0068 nm/px).  A baseline-
# subtracted centroid on such a line is good to roughly sigma/SNR.  Requiring
# that precision to be at most half of ALIGNMENT_MAX_PIXEL_RESIDUAL in the worst
# (bluest) order gives 5.1 px / 0.25 px = 20.  Below SNR 20 a line's own
# centroid noise is a sizeable fraction of the bound the verdict is stated in,
# and the +/-0.4 nm fit window -- +/-37 px at H-alpha, +/-60 px in the blue --
# is wide enough to reach a neighbour that drags a faint centroid many pixels
# across it.  That is not a hypothesis: on 2026-08-18 the strong Balmer lines
# (SNR 68..440) agreed to under a pixel on a flawless era match while weak
# Fulcher lines near the detection floor scattered over -10..+17 px in the same
# shot.  The strong lines carry the alignment truth; the weak ones describe how
# clean the data is.
#
# It is also 5x MINIMUM_SNR, which was always a *detection* floor -- enough to
# say a line is there -- and never a metrology floor.
VERDICT_SNR_FLOOR = 20.0

# --- The isotope/era degeneracy ---------------------------------------------
#
# ISOTOPE_FLIP_MINIMUM_LINES — how many two-reference Balmer windows a shot must
# hold before their agreement means anything.  One window flipping is an
# assignment; two or more flipping together is a pattern.
ISOTOPE_FLIP_MINIMUM_LINES = 2

# ISOTOPE_DEGENERACY_WINDOW_PX — how near the hydrogen reading's implied common
# shift must fall to the H/D separation before the two readings are declared
# indistinguishable.  Equivalently (the algebra collapses to it) how near the
# deuterium reading's own common residual must be to zero: that is the "better
# looking" fit the alarm exists to refuse.  4 px is half the ~8 px midpoint past
# which nearest-isotope assignment cannot separate the two references at all --
# inside it either hypothesis explains the centroids, outside it the deuterium
# reading is itself several pixels off, so the shot's problem is a shift the
# ordinary verdict already reports and not the isotope question.
ISOTOPE_DEGENERACY_WINDOW_PX = 4.0

#: The verdict a shot earns when its lines read as deuterium, the hydrogen
#: reading of the same centroids is one common H/D-sized shift, and the calendar
#: says deuterium was not running.  Named rather than folded into another word:
#: no correction repairs it, because nothing is broken except which calibration
#: epoch the cube was processed against.
ERA_MISASSIGNED_VERDICT = "era-misassigned-calibration"

# --- Plasma-presence frame selection ---------------------------------------
#
# A frame counts as plasma-bright when its total stands PLASMA_FRAME_SIGMA
# robust deviations above the *background* level, which is the same 5-sigma
# spirit ``Spectrum`` uses to recognise its dark frames.  A median over all
# frames instead would dilute every line with the dark frames beside it.
#
# The background is read off the darker of two clusters, never off the median
# of every frame.  The median is only the background when most frames are dark,
# and long shots break that assumption: on a real 40-frame shot with plasma up
# for 30 of them the median lands *inside* the plasma cluster, the threshold
# tops the brightest frame, and the audit reports insufficient-data on a shot
# that is three-quarters lit.  Splitting the sorted totals at the split that
# separates them best makes no assumption about which side holds the majority,
# so a 2-of-20 shot and a 36-of-40 shot are read the same way.
#
# A split alone would find two clusters in pure noise as well, so the 5-sigma
# rule survives as the guard on it: the dimmest frame above the split must
# still clear the dark cluster's own level by PLASMA_FRAME_SIGMA of that
# cluster's robust deviation.  A shot with no plasma in it fails that guard and
# stays honestly unmeasured rather than being handed a verdict.
PLASMA_FRAME_SIGMA = 5.0

# --- Isotopologues ----------------------------------------------------------
#
# ISOTOPE_FAMILY is the one audited family whose windows hold two references.
# BALMER_ISOTOPES is the order the evidence lists a window's candidates in.
ISOTOPE_FAMILY = "balmer"
BALMER_ISOTOPES = ("H", "D")

# HYDROGEN_ONLY_FAMILIES describe hydrogen and nothing else, so a shot that
# shows deuterium may not be judged by them.  The bundled Fulcher Q-branch
# anchors are H2; there is no D2 table on this side, and a D2 band fitted
# against H2 positions is the same systematic misread as D-alpha fitted against
# H-alpha, wearing a molecular hat.  Dropping the anchors costs quorum, which
# costs at worst an honest insufficient-data.  When a D2 catalog arrives it
# enters through the line catalog's isotope facet and this exclusion lifts for
# it without a change here.
HYDROGEN_ONLY_FAMILIES = ("fulcher",)

#: The bundled, owner-editable LHD deuterium calendar, cited by the evidence.
DEUTERIUM_CALENDAR_RESOURCE = "echelle_spectra/resources/lhd_deuterium_campaign.toml"

#: The cube attributes the audit reads for an acquisition date, in the order it
#: tries them: SpectroCube's ISO acquisition start, then a calendar date inside
#: the recorded source filename, then the export timestamp as a last resort.
DATE_ATTRIBUTES = ("t_start", "source_file", "created_at")

_DATE_PATTERN = re.compile(
    r"(?<!\d)(?P<year>20\d{2})[-_]?(?P<month>0[1-9]|1[0-2])"
    r"[-_]?(?P<day>0[1-9]|[12]\d|3[01])(?!\d)"
)


class DriftError(ValueError):
    """Drift evidence or refinement is unsafe or incomplete."""


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def resolve_cube_paths(inputs: list[str | Path]) -> list[Path]:
    """Expand directories to the cubes under them, shell-neutrally.

    A directory contributes every ``*.nc`` file beneath it in sorted order, so
    one audit command reads the same set on PowerShell and on a POSIX shell
    without either one's glob syntax.
    """

    resolved: list[Path] = []
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            # A Mac writing to exFAT or SMB leaves an AppleDouble metadata
            # sibling (._cube.nc) beside every real cube; it is not a dataset
            # and must never be audited.
            found = sorted(item for item in path.rglob("*.nc") if not item.name.startswith("."))
            if not found:
                raise DriftError(f"no .nc cubes found under {path}")
            resolved.extend(found)
        elif path.is_file():
            resolved.append(path)
        else:
            raise DriftError(f"cube path does not exist: {path}")
    if not resolved:
        raise DriftError("no cubes were selected")
    return sorted(dict.fromkeys(resolved))


def _shot_identity(path: Path) -> tuple[int | None, str]:
    """Return the parsed shot number and the raw stem token for one cube."""

    from .calibration_registry import CalibrationRegistryError, source_identity_from_path

    try:
        identity = source_identity_from_path(path)
    except CalibrationRegistryError:
        return None, path.stem
    return identity.shot_number, path.stem


def _matches_shot(path: Path, shots: set[str]) -> bool:
    """Match a whole shot token, never a substring: 42 is not part of 142."""

    number, stem = _shot_identity(path)
    for shot in shots:
        token = shot.strip()
        if not token:
            continue
        if number is not None and token.isdigit() and int(token) == number:
            return True
        if token == stem:
            return True
    return False


def select_sample_paths(
    paths: list[str | Path], *, every: int = 1, shots: set[str] | None = None
) -> list[Path]:
    """Select interval samples plus explicitly named shots."""

    if every < 1:
        raise DriftError("sample interval must be at least 1")
    ordered = sorted(Path(path) for path in paths)
    chosen = {path for index, path in enumerate(ordered) if index % every == 0}
    if shots:
        chosen.update(path for path in ordered if _matches_shot(path, shots))
    return sorted(chosen)


def _iso_date(value: str) -> date | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        return date.fromisoformat(text[:10])
    except ValueError:
        pass
    match = _DATE_PATTERN.search(text.replace("\\", "/"))
    if match is None:
        return None
    try:
        return date(int(match["year"]), int(match["month"]), int(match["day"]))
    except ValueError:
        return None


def cube_date(attrs: dict[str, Any]) -> tuple[date | None, str]:
    """Return one cube's acquisition date and the attribute it came from."""

    for key in DATE_ATTRIBUTES:
        parsed = _iso_date(str(attrs.get(key, "")))
        if parsed is not None:
            return parsed, key
    return None, ""


def filter_by_date(
    paths: list[Path], *, date_from: str | None, date_to: str | None
) -> list[Path]:
    """Keep cubes whose acquisition date lies inside an inclusive interval."""

    import xarray as xr

    try:
        low = date.fromisoformat(date_from) if date_from else None
        high = date.fromisoformat(date_to) if date_to else None
    except ValueError as exc:
        raise DriftError(f"--from/--to take ISO YYYY-MM-DD dates: {exc}") from exc
    if low is not None and high is not None and low > high:
        raise DriftError("--from must not be after --to")
    selected = []
    for path in paths:
        with xr.open_dataset(path) as opened:
            attrs = dict(opened.attrs)
        when, _source = cube_date(attrs)
        if when is None:
            raise DriftError(
                f"{path.name} carries no acquisition date in "
                + ", ".join(DATE_ATTRIBUTES)
                + "; date selection would silently drop it, so audit it by shot instead"
            )
        if (low is None or when >= low) and (high is None or when <= high):
            selected.append(path)
    if not selected:
        raise DriftError(
            "no sampled cube falls between "
            f"{date_from or '-inf'} and {date_to or '+inf'}"
        )
    return selected


#: Kept so in-package callers that imported the private name keep working.
_filter_by_date = filter_by_date


# ---------------------------------------------------------------------------
# The deuterium calendar, as a prior
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def load_deuterium_calendar() -> dict[str, Any]:
    """Return the bundled LHD deuterium calendar, citations and all."""

    resource = files("echelle_spectra.resources").joinpath("lhd_deuterium_campaign.toml")
    return tomllib.loads(resource.read_text(encoding="utf-8"))


def deuterium_prior(when: date | None, shot_number: str = "") -> dict[str, Any]:
    """Say whether deuterium was even possible for one shot, and on what basis.

    This is a prior and only a prior.  Inside a deuterium window the isotope is
    an open question, because LHD ran hydrogen shots inside its deuterium
    cycles too; outside every window hydrogen is expected.  Nothing here ever
    decides an isotope, moves a residual, or touches a verdict -- the audit
    compares this expectation with what the spectrum actually shows and records
    the disagreement when there is one.
    """

    calendar = load_deuterium_calendar()
    windows = list(calendar.get("windows", ()))
    common = {"calendar": DEUTERIUM_CALENDAR_RESOURCE}
    if when is not None:
        for window in windows:
            low = date.fromisoformat(str(window["date_from"]))
            high = date.fromisoformat(str(window["date_to"]))
            if low <= when <= high:
                return {
                    **common,
                    "expectation": str(window.get("expectation", "deuterium possible")),
                    "window": str(window["name"]),
                    "basis": (
                        f"{when.isoformat()} falls inside {window['name']} "
                        f"({low.isoformat()}..{high.isoformat()})"
                    ),
                }
        return {
            **common,
            "expectation": "hydrogen expected",
            "basis": (
                f"{when.isoformat()} falls outside every deuterium window the "
                "calendar lists"
            ),
        }
    floors = [int(window["shot_from"]) for window in windows if "shot_from" in window]
    token = str(shot_number).strip()
    if floors and token.isdigit():
        first = min(floors)
        if int(token) < first:
            return {
                **common,
                "expectation": "hydrogen expected",
                "basis": f"shot {token} precedes the calendar's first deuterium shot {first}",
            }
        return {
            **common,
            "expectation": "unknown",
            "basis": (
                f"shot {token} is at or beyond the first deuterium shot {first}, but the "
                "calendar states no closing shot number; only an acquisition date places it"
            ),
        }
    return {
        **common,
        "expectation": "unknown",
        "basis": "no acquisition date and no shot number that the calendar can place",
    }


def _isotope_prior_provenance() -> dict[str, Any]:
    """Describe the calendar the evidence's priors came from, inside the evidence."""

    calendar = load_deuterium_calendar()
    return {
        "calendar": DEUTERIUM_CALENDAR_RESOURCE,
        "schema": str(calendar.get("schema", "")),
        "facility": str(calendar.get("facility", "")),
        "role": (
            "prior only: it says where the isotope question exists, never what the "
            "answer is, and never enters the verdict arithmetic"
        ),
        "windows": [
            {
                "name": str(window["name"]),
                "expectation": str(window.get("expectation", "")),
                "date_from": str(window["date_from"]),
                "date_to": str(window["date_to"]),
                **({"shot_from": int(window["shot_from"])} if "shot_from" in window else {}),
            }
            for window in calendar.get("windows", ())
        ],
    }


def _isotope_flag(prior: dict[str, Any], *, shows_deuterium: bool) -> str:
    """Report a calendar/measurement disagreement without resolving it.

    Only one direction is a disagreement.  Hydrogen measured inside a deuterium
    window is ordinary -- the window never predicted deuterium.  Deuterium
    measured where the calendar expects hydrogen is worth a second look, and it
    is worth it precisely because a real blueward shift past ~8 px reads as
    deuterium; the flag names that ambiguity instead of picking a side.
    """

    if not shows_deuterium or prior.get("expectation") != "hydrogen expected":
        return ""
    return (
        "measured deuterium where the bundled LHD deuterium calendar expects hydrogen ("
        + str(prior.get("basis", ""))
        + "); the spectroscopic fit stands and this changes no residual or verdict, but a "
        "real blueward detector shift past ~8 px reads exactly like deuterium, so confirm "
        "which of the two this is before accepting any correction"
    )


#: How the calendar's expectation reads inside the degeneracy finding.
_CALENDAR_WORDS = {
    "hydrogen expected": "H-only",
    "deuterium possible": "D-possible",
}


def isotope_ambiguity(
    lines: list[dict[str, Any]], prior: dict[str, Any]
) -> dict[str, Any] | None:
    """Alarm when a shot's isotope reading and an era shift explain it equally well.

    Every audited Balmer window holds two references ~16.5 px apart, and a rigid
    era shift moves every line by the same pixels.  So the two readings are
    indistinguishable exactly when *all* of a shot's two-reference windows chose
    deuterium **and** the hydrogen reading of the same centroids implies one
    common shift the size of that separation -- which, algebraically, is the
    same statement as "the deuterium fit's own residual is near zero".  That is
    the trap the 2026-08-18 rehearsal walked into: a 2019 cube recalibrated onto
    the wrong era had its hydrogen lines re-identified as deuterium and the
    flipped fit then read *better* than the correct one.

    A single hydrogen assignment anywhere in the shot refutes the era-shift
    hypothesis outright -- one translation cannot move three lines onto
    deuterium and leave a fourth on hydrogen -- so a mixed shot is not
    ambiguous, it is mixed, and this returns nothing for it.

    Returns ``None`` when there is no ambiguity to report.  It never changes an
    assignment or a residual; the caller decides what the finding is worth.
    """

    baseline, flipped_isotope = BALMER_ISOTOPES
    windows = [
        item
        for item in lines
        if item.get("status") == "measured" and item.get("isotope_candidates")
    ]
    if len(windows) < ISOTOPE_FLIP_MINIMUM_LINES:
        return None
    if any(str(item.get("isotope")) != flipped_isotope for item in windows):
        return None
    baseline_px: list[float] = []
    separations: list[float] = []
    for item in windows:
        candidates = {
            str(entry.get("isotope")): entry for entry in item["isotope_candidates"]
        }
        if not {baseline, flipped_isotope} <= set(candidates):
            return None
        anchor = float(candidates[baseline]["pixel_residual_px"])
        moved = float(candidates[flipped_isotope]["pixel_residual_px"])
        baseline_px.append(anchor)
        separations.append(moved - anchor)
    implied_shift = float(np.median(baseline_px))
    separation = float(np.median(separations))
    # The flipped fit's own common residual: how much better it looks.
    degeneracy = abs(implied_shift + separation)
    if not np.isfinite(degeneracy) or degeneracy > ISOTOPE_DEGENERACY_WINDOW_PX:
        return None
    expectation = str(prior.get("expectation", "unknown"))
    calendar = _CALENDAR_WORDS.get(expectation, "unknown")
    excludes = expectation == "hydrogen expected"
    finding = (
        f"isotope flip or era shift, degenerate; calendar says {calendar}"
    )
    return {
        "finding": finding,
        "flipped_lines": len(windows),
        "baseline_isotope": baseline,
        "measured_isotope": flipped_isotope,
        "implied_common_shift_px": implied_shift,
        "isotope_separation_px": separation,
        "flipped_fit_residual_px": degeneracy,
        "degeneracy_window_px": ISOTOPE_DEGENERACY_WINDOW_PX,
        "calendar": calendar,
        "calendar_expectation": expectation,
        "calendar_basis": str(prior.get("basis", "")),
        "excludes_deuterium": bool(excludes),
        "detail": (
            f"all {len(windows)} two-reference Balmer window(s) chose {flipped_isotope}, and "
            f"reading the same centroids as {baseline} implies one common shift of "
            f"{implied_shift:+.2f} px against an H/D separation of {separation:+.2f} px at "
            f"these orders -- the two explanations differ by {degeneracy:.2f} px, inside the "
            f"{ISOTOPE_DEGENERACY_WINDOW_PX:g} px window where neither can be told from the "
            "other. "
            + (
                "The calendar excludes deuterium for this shot, so the likelier reading is a "
                "calibration from the wrong epoch: the fit that looks better is the flipped "
                "one. Re-run the cube against the calibration its own date selects before "
                "trusting any number here."
                if excludes
                else "The calendar allows deuterium here, so the spectroscopic assignment "
                "stands as measured -- but confirm the epoch before accepting a correction, "
                "because a wrong-era calibration would look exactly like this."
            )
        ),
    }


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


def _robust_deviation(values: np.ndarray) -> float:
    """Return a robust standard deviation, or 0.0 when the sample carries none."""

    if values.size == 0:
        return 0.0
    deviation = 1.4826 * float(np.median(np.abs(values - np.median(values))))
    if not np.isfinite(deviation) or deviation <= 0.0:
        deviation = float(np.std(values))
    return deviation if np.isfinite(deviation) and deviation > 0.0 else 0.0


def _dark_cluster_size(ordered: np.ndarray) -> int:
    """Return how many of the ascending totals form the darker of two clusters.

    The split maximises the separation between the two sides, weighted by how
    many frames each holds -- the classic two-class threshold.  It is looked for
    at every boundary, so it is indifferent to whether the dark frames are the
    majority or the minority.  A boundary inside a run of equal totals would not
    separate anything and is skipped; when no boundary separates at all (one
    frame, or every total identical) the answer is 0 and there is no dark
    cluster to measure a background from.
    """

    count = ordered.size
    if count < 2:
        return 0
    cumulative = np.cumsum(ordered)
    best_size, best_score = 0, -np.inf
    for size in range(1, count):
        if ordered[size] <= ordered[size - 1]:
            continue
        dark_mean = cumulative[size - 1] / size
        bright_mean = (cumulative[-1] - cumulative[size - 1]) / (count - size)
        score = size * (count - size) * (bright_mean - dark_mean) ** 2
        if score > best_score:
            best_size, best_score = size, score
    return best_size


def select_plasma_frames(intensity_2d: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    """Return the plasma-bright frame indices and the evidence for the choice."""

    totals = np.nansum(intensity_2d, axis=1)
    if totals.size == 1:
        return np.array([0]), {
            "rule": "single frame",
            "frames": 1,
            "selected_frames": 1,
        }
    finite = np.isfinite(totals)
    ordered = np.sort(totals[finite])
    dark_size = _dark_cluster_size(ordered)
    # With no split to find, every frame reads as background and none is bright.
    dark = ordered[:dark_size] if dark_size else ordered
    background = float(np.median(dark)) if dark.size else 0.0
    deviation = _robust_deviation(dark)
    if deviation <= 0.0:
        # A one-frame dark cluster cannot state its own noise; borrow the scatter
        # of the whole shot, which is never smaller and so never over-selects.
        deviation = _robust_deviation(ordered)
    threshold = background + PLASMA_FRAME_SIGMA * deviation
    selected = np.array([], dtype=int)
    if dark_size > 0 and float(ordered[dark_size]) > threshold:
        selected = np.flatnonzero(finite & (totals >= ordered[dark_size]))
    return selected, {
        "rule": (
            "frame total in the brighter of two clusters, whose dimmest frame "
            f"clears the darker cluster by {PLASMA_FRAME_SIGMA:g} sigma"
        ),
        "frames": int(totals.size),
        "selected_frames": int(selected.size),
        "dark_frames": int(dark.size),
        "background_total": background,
        "sigma_total": float(deviation),
        "threshold_total": float(threshold),
    }


def _spectrum(ds) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Return the wavelength axis and the median over plasma-bright frames."""

    wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
    intensity = ds["intensity"]
    reduce_dims = [dim for dim in intensity.dims if dim != "wavelength"]
    if not reduce_dims:
        return (
            wavelength,
            np.asarray(intensity.values, dtype=float),
            {"rule": "single frame", "frames": 1, "selected_frames": 1},
        )
    stacked = np.asarray(
        intensity.transpose(*reduce_dims, "wavelength").values, dtype=float
    ).reshape(-1, wavelength.size)
    frames, selection = select_plasma_frames(stacked)
    if frames.size == 0:
        return wavelength, np.full(wavelength.size, np.nan), selection
    return (
        wavelength,
        np.nanmedian(stacked[frames], axis=0),
        selection,
    )


def centroid_evidence(
    wavelength: np.ndarray,
    intensity: np.ndarray,
    *,
    expected_nm: float,
    window_nm: float = 0.4,
) -> dict[str, Any]:
    """Measure a baseline-subtracted centroid with explicit sufficiency evidence.

    The baseline and the noise are read off the whole window robustly, never off
    its outer fifth.  An edge estimate is contaminated the moment the line moves
    towards an edge, so the reported SNR fell as the shift grew -- on the
    fixture's four-order cube an 18 px rigid shift dragged H-alpha's SNR from
    2400 to 19, purely because the line's own wing had reached the strip the
    baseline was read from.  That was tolerable while SNR only gated detection;
    now that it decides which lines carry the verdict, an SNR that depends on
    the very shift being measured would let a large honest shift disqualify its
    own evidence.  The window median is the baseline because an emission line
    occupies about a quarter of a +/-0.4 nm window, and the noise is measured
    from the samples at or below that baseline, which no emission line can
    reach.
    """

    selected = np.abs(wavelength - expected_nm) <= window_nm
    x = wavelength[selected]
    y = intensity[selected]
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 5:
        return {"status": "insufficient-data", "reason": "fewer than five samples in window"}
    baseline = float(np.median(y))
    signal = np.clip(y - baseline, 0.0, None)
    below = (y - baseline)[y <= baseline]
    # For Gaussian noise E[x^2 | x <= 0] is exactly the variance, so the RMS of
    # the samples at or under the baseline is an unbiased sigma that no emission
    # line can enter.  The RMS rather than a MAD because it uses every one of
    # those samples: a MAD over half a window is noisy enough that a pure-noise
    # window crosses MINIMUM_SNR measurably more often than it should.
    noise = float(np.sqrt(np.mean(below**2))) if below.size else 0.0
    peak = float(np.max(signal))
    snr = peak / max(noise, np.finfo(float).eps)
    if not np.isfinite(snr) or snr < MINIMUM_SNR or signal.sum() <= 0:
        return {
            "status": "insufficient-data",
            "reason": f"peak SNR {snr:.3g} is below {MINIMUM_SNR:g}",
            "snr": snr,
        }
    centroid = float(np.sum(x * signal) / np.sum(signal))
    return {
        "status": "measured",
        "centroid_nm": centroid,
        "residual_nm": centroid - expected_nm,
        "snr": snr,
        "window_nm": window_nm,
    }


class CubeGeometry:
    """One cube's own detector geometry: per-order polynomials and coverage.

    Packet 8 stores the wavelength polynomials, the pre-flip ``detector_pixel``
    and the ``echelle_order`` of every sample precisely so a later audit can ask
    what a wavelength residual is worth in pixels without reopening raw data.
    """

    def __init__(self, ds):
        raw = ds.attrs.get("wavelength_polynomials_json")
        if not raw:
            raise DriftError(
                "cube carries no wavelength_polynomials_json; the drift judge measures "
                "detector pixels, so re-export it through a calibration snapshot"
            )
        try:
            payload = json.loads(str(raw))
            self.coefficients = {
                int(record["order"]): [float(value) for value in record["coefficients"]]
                for record in payload["orders"]
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise DriftError(f"cube has malformed wavelength_polynomials_json: {exc}") from exc
        missing = [name for name in ("detector_pixel", "echelle_order") if name not in ds.coords]
        if missing:
            raise DriftError(
                "cube lacks the aligned " + ", ".join(missing) + " coordinate(s) the "
                "detector-space judge needs"
            )
        self.wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
        self.pixel = np.asarray(ds.coords["detector_pixel"].values, dtype=float)
        self.order = np.asarray(ds.coords["echelle_order"].values).astype(int)
        self.pixel_range = {
            int(order): (
                float(np.min(self.pixel[self.order == order])),
                float(np.max(self.pixel[self.order == order])),
            )
            for order in np.unique(self.order)
        }

    @property
    def orders(self) -> tuple[int, ...]:
        return tuple(sorted(self.pixel_range))

    def sample_at(self, wavelength_nm: float) -> tuple[int, float]:
        """Return the ``(order, detector_pixel)`` of the nearest stored sample."""

        index = int(np.argmin(np.abs(self.wavelength - float(wavelength_nm))))
        return int(self.order[index]), float(self.pixel[index])

    def _clamped(self, order: int, pixel: float) -> tuple[float, bool]:
        low, high = self.pixel_range[int(order)]
        clamped = not (low <= pixel <= high)
        return float(min(max(pixel, low), high)), clamped

    def dispersion(self, order: int, pixel: float) -> tuple[float, bool]:
        """Return d(lambda)/dx in nm/px, clamped to the order's covered pixels."""

        if int(order) not in self.coefficients:
            raise DriftError(f"cube has no wavelength polynomial for order {order}")
        clamped_pixel, clamped = self._clamped(order, pixel)
        derivative = np.polyder(np.asarray(self.coefficients[int(order)], dtype=float))
        return float(np.polyval(derivative, clamped_pixel)), clamped

    def wavelength_at(self, order: int, pixel: float) -> tuple[float, bool]:
        """Return the order's wavelength at one pixel, clamped to its coverage."""

        clamped_pixel, clamped = self._clamped(order, pixel)
        return (
            float(np.polyval(np.asarray(self.coefficients[int(order)], dtype=float), clamped_pixel)),
            clamped,
        )


def _alignment_tolerance_px(accuracy_nm: float | None, dispersion_nm_per_px: float) -> float:
    """Return the alignment bound for one line, in detector pixels."""

    tolerance = ALIGNMENT_MAX_PIXEL_RESIDUAL
    if accuracy_nm is not None and np.isfinite(accuracy_nm) and accuracy_nm > 0:
        derived = float(accuracy_nm) / abs(dispersion_nm_per_px)
        tolerance = max(tolerance, derived)
    return float(tolerance)


def _in_detector_space(
    geometry: CubeGeometry, measured: dict[str, Any], *, accuracy_nm: float | None
) -> dict[str, Any]:
    """Turn one measured wavelength residual into a detector-pixel residual."""

    order, pixel = geometry.sample_at(measured["centroid_nm"])
    try:
        dispersion, clamped = geometry.dispersion(order, pixel)
    except DriftError as exc:
        return {"status": "insufficient-data", "reason": str(exc)}
    if not np.isfinite(dispersion) or dispersion == 0.0:
        return {
            "status": "insufficient-data",
            "reason": f"order {order} has no usable dispersion at pixel {pixel:.1f}",
        }
    tolerance_px = _alignment_tolerance_px(accuracy_nm, dispersion)
    result = {
        **measured,
        "echelle_order": int(order),
        "detector_pixel": float(pixel),
        "dispersion_nm_per_px": float(dispersion),
        "pixel_residual_px": float(measured["residual_nm"] / dispersion),
        "alignment_tolerance_px": tolerance_px,
        "alignment_tolerance_nm": float(tolerance_px * abs(dispersion)),
    }
    if clamped:
        low, high = geometry.pixel_range[int(order)]
        result["extrapolated_dispersion"] = True
        result["extrapolation_note"] = (
            f"dispersion evaluated at order {order}'s covered pixel range "
            f"[{low:.1f}, {high:.1f}] instead of extrapolating the polynomial"
        )
    return result


def _cube_accuracy_nm(attrs: dict[str, Any]) -> float | None:
    raw = attrs.get("wavelength_accuracy_nm")
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def _isotope_references() -> dict[str, dict[str, SpectralLine]]:
    """Pair each audited Balmer transition with one reference per isotopologue."""

    paired: dict[str, dict[str, SpectralLine]] = {}
    for isotope in BALMER_ISOTOPES:
        for line in load_line_table(ISOTOPE_FAMILY, isotope=isotope):
            paired.setdefault(line.transition, {})[isotope] = line
    return paired


def _measure_line(
    wavelength: np.ndarray,
    intensity: np.ndarray,
    geometry: CubeGeometry,
    *,
    anchor: SpectralLine,
    candidates: dict[str, SpectralLine],
    accuracy_nm: float | None,
) -> dict[str, Any]:
    """Measure one window and, when it holds two references, choose between them.

    The window and the centroid are exactly what a single-isotope audit
    measures: one baseline-subtracted centroid around the hydrogen anchor, in
    the same +/-0.4 nm window.  The isotope question is asked afterwards, of
    that one number, so a hydrogen shot is measured and judged identically to
    how it was before deuterium was considered at all.
    """

    result = centroid_evidence(wavelength, intensity, expected_nm=anchor.wavelength_nm)
    if result.get("status") != "measured":
        return result
    detector = _in_detector_space(geometry, result, accuracy_nm=accuracy_nm)
    if detector.get("status") != "measured" or len(candidates) < 2:
        return detector
    dispersion = float(detector["dispersion_nm_per_px"])
    centroid = float(detector["centroid_nm"])
    evaluated = [
        {
            "isotope": isotope,
            "line": line.label,
            "expected_nm": float(line.wavelength_nm),
            "residual_nm": centroid - float(line.wavelength_nm),
            "pixel_residual_px": (centroid - float(line.wavelength_nm)) / dispersion,
            "source_reference": line.source_reference,
        }
        for isotope, line in sorted(candidates.items())
    ]
    # Nearest in pixels, not in nanometres: the two references are ~16.5 px
    # apart in every audited order, and pixels are the space this judge works
    # in.  Both candidates are kept so the evidence shows what was rejected and
    # by how much -- the reader can always recover the hydrogen-only reading.
    nearest = min(evaluated, key=lambda item: abs(item["pixel_residual_px"]))
    return {**detector, **nearest, "isotope_candidates": evaluated}


def _shows_deuterium(lines: list[dict[str, Any]]) -> bool:
    """True when any measured line in this shot was assigned to deuterium."""

    return any(
        item.get("status") == "measured" and item.get("isotope") == "D" for item in lines
    )


def _majority_isotope(lines: list[dict[str, Any]]) -> str:
    """Return the isotopologue most of a shot's measured lines were assigned to.

    A tie reports ``mixed``.  LHD ran hydrogen and deuterium within the same
    cycles and a shot can genuinely show both, so an even split is stated
    rather than broken by an arbitrary rule.
    """

    counts: dict[str, int] = {}
    for item in lines:
        if item.get("status") != "measured":
            continue
        name = str(item.get("isotope", ""))
        if name:
            counts[name] = counts.get(name, 0) + 1
    if not counts:
        return ""
    highest = max(counts.values())
    leaders = sorted(name for name, count in counts.items() if count == highest)
    return leaders[0] if len(leaders) == 1 else "mixed"


def _isotope_tally(sampled: list[dict[str, Any]]) -> dict[str, int]:
    """Count the sampled shots by the isotopologue each of them read as."""

    counts: dict[str, int] = {}
    for item in sampled:
        name = str(item.get("isotope") or "undetermined")
        counts[name] = counts.get(name, 0) + 1
    return dict(sorted(counts.items()))


def _measure_cube_lines(
    wavelength: np.ndarray,
    intensity: np.ndarray,
    geometry: CubeGeometry,
    *,
    families: tuple[str, ...],
    common: dict[str, Any],
    accuracy_nm: float | None,
) -> tuple[list[dict[str, Any]], str]:
    """Measure one cube's audited lines and return them with its isotope tag.

    The isotope-bearing family is measured first, because whether this shot
    shows deuterium decides whether the hydrogen-only families may be audited
    at all.  Excluded lines are still written into the evidence, saying why
    they were dropped; a silently shorter table would look like a coverage
    problem rather than a decision.
    """

    minimum, maximum = float(np.min(wavelength)), float(np.max(wavelength))
    ordered = [family for family in families if family == ISOTOPE_FAMILY]
    ordered += [family for family in families if family != ISOTOPE_FAMILY]
    references = _isotope_references()
    records: list[dict[str, Any]] = []
    isotope = ""
    excluding = False
    for family in ordered:
        excluded = excluding and family in HYDROGEN_ONLY_FAMILIES
        for line in load_line_table(family):
            if not minimum <= line.wavelength_nm <= maximum:
                continue
            candidates = (
                references.get(line.transition, {}) if family == ISOTOPE_FAMILY else {}
            )
            entry = {
                **common,
                "family": family,
                "line": line.label,
                "expected_nm": line.wavelength_nm,
                "source_reference": line.source_reference,
                "blended": bool(line.blended),
                # Undetermined until a centroid picks between two references;
                # a one-reference family already knows what it describes.
                "isotope": "" if len(candidates) > 1 else line.isotope,
            }
            if excluded:
                records.append(
                    {
                        **entry,
                        "status": "skipped",
                        "isotope_excluded": True,
                        "reason": (
                            "this shot reads deuterium and these anchors are H2; no D2 "
                            "table exists here, so the line is dropped rather than judged "
                            "against the wrong molecule"
                        ),
                    }
                )
                continue
            if line.blended:
                records.append(
                    {
                        **entry,
                        "status": "skipped",
                        "reason": (
                            "sub-resolution blend: the measured centroid belongs to the "
                            "blend, not to this transition"
                        ),
                    }
                )
                continue
            records.append(
                {
                    **entry,
                    **_measure_line(
                        wavelength,
                        intensity,
                        geometry,
                        anchor=line,
                        candidates=candidates,
                        accuracy_nm=accuracy_nm,
                    ),
                }
            )
        if family == ISOTOPE_FAMILY:
            isotope = _majority_isotope(records)
            excluding = _shows_deuterium(records)
    return records, isotope


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


def _resolved(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        item
        for item in lines
        if item.get("status") == "measured"
        and not item.get("blended")
        and np.isfinite(float(item.get("pixel_residual_px", np.nan)))
    ]


def _line_snr(item: dict[str, Any]) -> float:
    """Return one measured line's peak SNR.

    A record without one is treated as strong: every producer in this package
    writes an SNR, so a missing value belongs to hand-composed evidence whose
    author never expressed a doubt for this function to honour.
    """

    try:
        value = float(item.get("snr", np.inf))
    except (TypeError, ValueError):
        return np.inf
    return value if np.isfinite(value) or value == np.inf else np.inf


def _verdict_lines(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return the resolved lines strong enough to decide an alignment verdict."""

    return [item for item in _resolved(lines) if _line_snr(item) >= VERDICT_SNR_FLOOR]


def _quality_lines(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return the resolved lines that describe data quality but never a verdict."""

    return [item for item in _resolved(lines) if _line_snr(item) < VERDICT_SNR_FLOOR]


def _weak_line_scatter(lines: list[dict[str, Any]]) -> dict[str, Any]:
    """Describe the weak lines' scatter in its own words, beside the verdict.

    This is a quality figure and says so.  A faint centroid measured inside a
    +/-0.4 nm window can be dragged whole pixels by a neighbour the window
    reaches, so this scatter reports how crowded and how clean the audited
    spectrum is -- not whether the detector moved.
    """

    weak = _quality_lines(lines)
    strong = _verdict_lines(lines)
    if not weak:
        return {
            "lines": 0,
            "snr_floor": VERDICT_SNR_FLOOR,
            "note": (
                f"every resolved line reached SNR {VERDICT_SNR_FLOOR:g}, so none was set "
                "aside from the verdict"
            ),
        }
    shifts = np.asarray([float(item["pixel_residual_px"]) for item in weak], dtype=float)
    snrs = np.asarray([_line_snr(item) for item in weak], dtype=float)
    median = float(np.median(shifts))
    spread = float(np.max(np.abs(shifts - median)))
    return {
        "lines": len(weak),
        "snr_floor": VERDICT_SNR_FLOOR,
        "snr_range": [float(np.min(snrs)), float(np.max(snrs))],
        "median_shift_px": median,
        "maximum_pixel_deviation_px": spread,
        "pixel_residual_range_px": [float(np.min(shifts)), float(np.max(shifts))],
        "note": (
            f"{len(weak)} resolved line(s) below SNR {VERDICT_SNR_FLOOR:g} scatter over "
            f"{float(np.min(shifts)):+.2f}..{float(np.max(shifts)):+.2f} px "
            f"({spread:.2f} px about their median). At that SNR a centroid alone cannot "
            "state a shift to better than about half a pixel, and the +/-0.4 nm fit window "
            "reaches far enough to catch a neighbour, so this is a figure for how crowded "
            "and how clean the spectrum is. It did not decide the verdict"
            + (
                f", which rests on {len(strong)} line(s) at or above the floor."
                if strong
                else " -- and no line reached the floor, so no verdict was earned."
            )
        ),
    }


def _quorum(
    lines: list[dict[str, Any]], coverage_nm: tuple[float, float] | None
) -> dict[str, Any]:
    """Report whether the verdict-carrying lines can carry a verdict at all."""

    measured = _verdict_lines(lines)
    set_aside = len(_resolved(lines)) - len(measured)
    orders = sorted({int(item["echelle_order"]) for item in measured})
    if coverage_nm is None:
        expected = [float(item["expected_nm"]) for item in lines if "expected_nm" in item]
        coverage_nm = (min(expected), max(expected)) if expected else (0.0, 0.0)
    midpoint = 0.5 * (float(coverage_nm[0]) + float(coverage_nm[1]))
    blue = any(float(item["expected_nm"]) < midpoint for item in measured)
    red = any(float(item["expected_nm"]) > midpoint for item in measured)
    reasons = []
    if len(measured) < MINIMUM_CENTROIDS:
        reasons.append(f"{len(measured)} resolved line(s) below the required {MINIMUM_CENTROIDS}")
    if len(orders) < MINIMUM_DISTINCT_ORDERS:
        reasons.append(
            f"{len(orders)} distinct order(s) below the required {MINIMUM_DISTINCT_ORDERS}"
        )
    if not (blue and red):
        reasons.append(
            "measured lines do not span both halves of the audited coverage "
            f"({coverage_nm[0]:.3f}-{coverage_nm[1]:.3f} nm)"
        )
    if reasons and set_aside:
        # Say which floor emptied the quorum, so an insufficient-data verdict on
        # a well-lit shot is never mistaken for a shot that showed nothing.
        reasons.append(
            f"{set_aside} further resolved line(s) sit below the SNR "
            f"{VERDICT_SNR_FLOOR:g} floor a verdict is decided on"
        )
    return {
        "resolved_lines": len(measured),
        "distinct_orders": len(orders),
        "orders": orders,
        "coverage_nm": [float(coverage_nm[0]), float(coverage_nm[1])],
        "bluer_half": bool(blue),
        "redder_half": bool(red),
        "below_snr_floor_lines": int(set_aside),
        "snr_floor": VERDICT_SNR_FLOOR,
        "satisfied": not reasons,
        "reason": "; ".join(reasons),
    }


def verdict_from_evidence(
    evidence: list[dict[str, Any]],
    *,
    coverage_nm: tuple[float, float] | None = None,
    isotope_ambiguity: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Classify one interval in detector space and report how it was classified.

    Only the lines above ``VERDICT_SNR_FLOOR`` decide the word; the rest are
    summarised beside it under ``weak_line_scatter``.  ``isotope_ambiguity``, when
    the caller found one whose calendar prior excludes deuterium, replaces an
    otherwise ordinary aligned/shifted reading with
    ``ERA_MISASSIGNED_VERDICT``: those residuals are real, but they were
    measured against references the calendar says the plasma cannot have had.
    """

    quorum = _quorum(evidence, coverage_nm)
    scatter = _weak_line_scatter(evidence)
    if not quorum["satisfied"]:
        return "insufficient-data", {"quorum": quorum, "weak_line_scatter": scatter}
    measured = _verdict_lines(evidence)
    shifts = np.asarray([float(item["pixel_residual_px"]) for item in measured], dtype=float)
    dispersions = np.asarray(
        [abs(float(item["dispersion_nm_per_px"])) for item in measured], dtype=float
    )
    median = float(np.median(shifts))
    maximum = float(np.max(np.abs(shifts)))
    spread = float(np.max(np.abs(shifts - median)))
    tolerances = np.asarray(
        [
            float(item.get("alignment_tolerance_px", ALIGNMENT_MAX_PIXEL_RESIDUAL))
            for item in measured
        ],
        dtype=float,
    )
    summary = {
        "median_shift_px": median,
        "maximum_absolute_pixel_residual_px": maximum,
        "maximum_pixel_deviation_px": spread,
        "alignment_tolerance_px": float(np.max(tolerances)),
        "shift_consistency_px": SHIFT_CONSISTENCY_PX,
        # Reporting only: one rigid pixel shift is worth different wavelengths
        # in different orders, which is the whole reason the judge left nm.
        "median_shift_nm_equivalents": [
            float(median * value) for value in (float(np.min(dispersions)), float(np.max(dispersions)))
        ],
        "quorum": quorum,
        "verdict_lines": len(measured),
        "verdict_snr_floor": VERDICT_SNR_FLOOR,
        "verdict_snr_range": [
            float(np.min([_line_snr(item) for item in measured])),
            float(np.max([_line_snr(item) for item in measured])),
        ],
        "verdict_basis": (
            f"decided on the {len(measured)} resolved line(s) at or above SNR "
            f"{VERDICT_SNR_FLOOR:g}; weaker lines are reported under weak_line_scatter and "
            "never enter the arithmetic"
        ),
        "weak_line_scatter": scatter,
    }
    # Two independent facts, because they mean different things.  Inconsistency
    # says no single translation fits these lines; largeness says the misfit is
    # bigger than a centroid could be wrong by.  Only both together are damage.
    inconsistent = spread > SHIFT_CONSISTENCY_PX
    large = maximum > GEOMETRIC_DAMAGE_PX
    summary["strong_lines_inconsistent"] = bool(inconsistent)
    summary["strong_lines_large"] = bool(large)
    summary["geometric_damage_px"] = GEOMETRIC_DAMAGE_PX
    if isotope_ambiguity is not None and isotope_ambiguity.get("excludes_deuterium"):
        summary["verdict_reason"] = str(isotope_ambiguity.get("finding", ""))
        summary["isotope_ambiguity"] = isotope_ambiguity
        return ERA_MISASSIGNED_VERDICT, summary
    if isotope_ambiguity is not None:
        summary["isotope_ambiguity"] = isotope_ambiguity
    if np.all(np.abs(shifts) <= tolerances):
        summary["verdict_reason"] = (
            "every verdict-carrying line sits inside its own alignment tolerance"
        )
        return "aligned", summary
    if inconsistent and large:
        summary["verdict_reason"] = (
            f"the verdict-carrying lines disagree by {spread:.2f} px, past the "
            f"{SHIFT_CONSISTENCY_PX:g} px one translation allows, and the worst of them sits "
            f"{maximum:.2f} px out, past the {GEOMETRIC_DAMAGE_PX:g} px a centroid can be "
            "wrong by: this is the detector geometry, not a table that slipped"
        )
        return "misaligned-beyond-repair", summary
    if abs(median) > REPAIR_LIMIT_PX:
        summary["verdict_reason"] = (
            f"the verdict-carrying lines agree on one shift of {median:+.2f} px, but it is "
            f"past the {REPAIR_LIMIT_PX:g} px a table correction may claim; the anchors must "
            "be re-identified against lamp data rather than slid"
        )
        return "misaligned-beyond-repair", summary
    summary["verdict_reason"] = (
        f"the verdict-carrying lines describe one rigid shift of {median:+.2f} px"
        + (
            f", though they scatter {spread:.2f} px about it -- past the "
            f"{SHIFT_CONSISTENCY_PX:g} px bound, but with the worst line only {maximum:.2f} px "
            f"out, inside the {GEOMETRIC_DAMAGE_PX:g} px that separates quality from damage"
            if inconsistent
            else ""
        )
    )
    return "shifted", summary


def _per_shot_summary(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group the residuals per shot, ordered by date then shot.

    ``median_shift_px`` and ``pixel_spread_px`` describe the lines that would
    decide a verdict, because those are the numbers an operator acts on and the
    ones ``detect_interval_split`` clusters.  A shot with no line above the SNR
    floor still gets a row -- from every resolved line it has, marked
    ``below_snr_floor`` -- so it is visible rather than absent.  The weaker
    lines' own scatter travels alongside as a quality figure.
    """

    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for item in _resolved(lines):
        key = (str(item.get("date", "")), str(item.get("shot_number", "")), str(item["cube"]))
        groups.setdefault(key, []).append(item)
    summary = []
    for (when, shot, cube), items in sorted(groups.items()):
        strong = [entry for entry in items if _line_snr(entry) >= VERDICT_SNR_FLOOR]
        carrying = strong or items
        shifts = np.asarray([float(item["pixel_residual_px"]) for item in carrying], dtype=float)
        median = float(np.median(shifts))
        row = {
            "shot_number": shot,
            "cube": cube,
            "date": when,
            "lines": len(carrying),
            "isotope": _majority_isotope(carrying),
            "median_shift_px": median,
            "pixel_spread_px": float(np.max(np.abs(shifts - median))),
            "resolved_lines": len(items),
            "weak_lines": len(items) - len(strong),
        }
        if not strong:
            row["below_snr_floor"] = True
        weak = [entry for entry in items if _line_snr(entry) < VERDICT_SNR_FLOOR]
        if weak and strong:
            faint = np.asarray(
                [float(entry["pixel_residual_px"]) for entry in weak], dtype=float
            )
            row["weak_line_spread_px"] = float(np.max(np.abs(faint - np.median(faint))))
        summary.append(row)
    return summary


def _split_boundary(per_shot: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], int] | None:
    """Return the shift-ordered shots and the widest gap when it is a real split."""

    if len(per_shot) < 2:
        return None
    ordered = sorted(per_shot, key=lambda item: item["median_shift_px"])
    medians = [item["median_shift_px"] for item in ordered]
    gaps = [medians[index + 1] - medians[index] for index in range(len(medians) - 1)]
    widest = int(np.argmax(gaps))
    if gaps[widest] <= SHIFT_CONSISTENCY_PX:
        return None
    return ordered, widest


def _shot_sort_key(shot: str) -> tuple[int, float, str]:
    return (0, float(shot), shot) if shot.isdigit() else (1, 0.0, shot)


def detect_interval_split(per_shot: list[dict[str, Any]]) -> str:
    """Flag two separated groups of per-shot shifts instead of condemning them.

    An interval that straddles a calibration change holds two consistent rigid
    shifts, not one inconsistent set of residuals.  When the per-shot medians
    fall into two clusters further apart than the consistency bound, the
    interval wants splitting at that boundary; the four verdict words still
    classify the interval, and this warning travels alongside them.
    """

    split = _split_boundary(per_shot)
    if split is None:
        return ""
    ordered, widest = split
    low = sorted((item["shot_number"] for item in ordered[: widest + 1]), key=_shot_sort_key)
    high = sorted((item["shot_number"] for item in ordered[widest + 1 :]), key=_shot_sort_key)
    return (
        f"residuals form two groups (shots {low[0]}..{low[-1]} vs "
        f"{high[0]}..{high[-1]}); consider splitting the interval at the boundary"
    )


def _label_shot_groups(per_shot: list[dict[str, Any]]) -> None:
    """Mark which cluster each shot belongs to when the interval looks split."""

    split = _split_boundary(per_shot)
    if split is None:
        return
    ordered, widest = split
    for position, item in enumerate(ordered):
        item["group"] = 1 if position <= widest else 2


# ---------------------------------------------------------------------------
# Repair composition
# ---------------------------------------------------------------------------


def _next_refinement_id(base_snapshot_id: str, root: Path) -> str:
    stem = re.sub(r"-r\d+$", "", str(base_snapshot_id))
    revision = 1
    while (root / f"{stem}-r{revision}").exists():
        revision += 1
    return f"{stem}-r{revision}"


def _drive_requirement(
    shots: list[str], catalog_path: str | Path | None
) -> str:
    """Name the drives holding the affected shots, or the flag that would."""

    if catalog_path is None:
        return (
            "recalibrating beyond repair needs the raw SIF data: pass "
            "--catalog <merged-catalog.json> to name the drives holding these shots"
        )
    from .catalog import load_catalog

    payload = load_catalog(catalog_path)
    sources = payload.get("sources") or [
        {"volume_label": payload.get("volume_label", "unknown"), "cubes": payload.get("cubes", [])}
    ]
    wanted = {str(shot) for shot in shots}
    drives = sorted(
        {
            str(source.get("volume_label", "unknown"))
            for source in sources
            if any(str(cube.get("shot_number", "")) in wanted for cube in source.get("cubes", []))
        }
    )
    if not drives:
        return (
            f"no catalog row matches the affected shot(s) {', '.join(sorted(wanted))}; "
            "rebuild the catalog or name another one with --catalog"
        )
    return (
        "recalibrating beyond repair needs the raw SIF data on drive(s) "
        + ", ".join(drives)
        + "; reprocess those shots against a re-measured calibration snapshot"
    )


def _repair_steps(
    *,
    evidence_path: Path,
    calibrations_root: Path,
    refined_id: str,
    refined_is_preview: bool,
    shift_px: float,
    cube_directories: list[Path],
) -> list[dict[str, str]]:
    """Compose the real, paste-ready repair sequence for a shifted interval."""

    snapshot_dir = calibrations_root / refined_id
    steps = [
        {
            "shell": "any",
            "purpose": "accept the sampled shift and emit the immutable -rN snapshot",
            # Paths are quoted so a campaign folder with a space survives both
            # shells, and the shift carries 12 significant digits: the
            # acknowledgement is exact, so a rounded command would refuse itself.
            "command": (
                f'echelle drift refine "{evidence_path}" '
                f'--calibrations "{calibrations_root}" --accept-shift {shift_px:.12g}'
            ),
        },
        {
            "shell": "any",
            "purpose": (
                f"repoint the registry [validity] entry at {refined_id}"
                + (
                    " (the id 'echelle drift refine' prints)"
                    if refined_is_preview
                    else ""
                )
                + " before any bulk run; the accepted shift authorizes only the refinement"
            ),
            "command": "",
        },
    ]
    # ``recal-cube`` revises one cube at a time, so the loop is written out in
    # both shells.  The refinement changed only the wavelength role, so
    # --wavelength-only is what it authorizes: the sphere pair and the integral
    # curve are the base snapshot's own bytes.
    for directory in cube_directories:
        purpose = f"revise every cube already exported from {directory}"
        steps.append(
            {
                "shell": "powershell",
                "purpose": purpose,
                "command": (
                    f'Get-ChildItem "{directory}\\*.nc" | ForEach-Object {{ '
                    f'echelle recal-cube $_.FullName -o ($_.FullName -replace "\\.nc$", '
                    f'"-{refined_id}.nc") --new-snapshot "{snapshot_dir}" --wavelength-only }}'
                ),
            }
        )
        steps.append(
            {
                "shell": "posix",
                "purpose": purpose,
                "command": (
                    f'for cube in "{directory.as_posix()}"/*.nc; do echelle recal-cube '
                    f'"$cube" -o "${{cube%.nc}}-{refined_id}.nc" '
                    f'--new-snapshot "{snapshot_dir.as_posix()}" --wavelength-only; done'
                ),
            }
        )
    return steps


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------


def audit_cubes(  # noqa: C901 - one readable pass over cubes, lines, and shots
    cube_paths: list[str | Path],
    *,
    every: int = 1,
    shots: set[str] | None = None,
    families: tuple[str, ...] = ("balmer", "fulcher"),
    date_from: str | None = None,
    date_to: str | None = None,
    catalog: str | Path | None = None,
    evidence_path: str | Path | None = None,
    calibrations_root: str | Path = "calibrations",
) -> dict[str, Any]:
    """Audit sampled saved cubes and return an immutable JSON-ready verdict."""

    import xarray as xr

    resolved = resolve_cube_paths(cube_paths)
    if date_from or date_to:
        resolved = filter_by_date(resolved, date_from=date_from, date_to=date_to)
    selected = select_sample_paths(resolved, every=every, shots=shots)

    per_line: list[dict[str, Any]] = []
    snapshot_ids: set[str] = set()
    sampled: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    ambiguities: list[dict[str, Any]] = []
    coverage: list[float] = []
    geometries: dict[str, CubeGeometry] = {}
    for path in selected:
        with xr.open_dataset(path) as opened:
            ds = opened.load()
        attrs = dict(ds.attrs)
        snapshot_ids.add(str(attrs.get("snapshot_id", "unassigned")))
        shot_number = str(attrs.get("shot_number", path.stem))
        when, date_source = cube_date(attrs)
        wavelength, intensity, frames = _spectrum(ds)
        accuracy_nm = _cube_accuracy_nm(attrs)
        prior = deuterium_prior(when, shot_number)
        record = {
            "cube": path.name,
            "sha256": sha256_file(path),
            "shot_number": shot_number,
            "date": when.isoformat() if when else "",
            "date_attribute": date_source,
            "frame_selection": frames,
            "wavelength_accuracy_nm": accuracy_nm,
            "isotope": "",
            "isotope_prior": prior["expectation"],
            "isotope_prior_basis": prior["basis"],
        }
        sampled.append(record)
        common = {
            "cube": path.name,
            "shot_number": shot_number,
            "date": record["date"],
        }
        if frames["selected_frames"] == 0:
            reason = "no plasma-bright frames"
            skipped.append({**common, "reason": reason})
            per_line.append({**common, "status": "insufficient-data", "reason": reason})
            continue
        try:
            geometry = CubeGeometry(ds)
        except DriftError as exc:
            skipped.append({**common, "reason": str(exc)})
            per_line.append({**common, "status": "insufficient-data", "reason": str(exc)})
            continue
        geometries[path.name] = geometry
        coverage.extend((float(np.min(wavelength)), float(np.max(wavelength))))
        measured, isotope = _measure_cube_lines(
            wavelength,
            intensity,
            geometry,
            families=families,
            common=common,
            accuracy_nm=accuracy_nm,
        )
        per_line.extend(measured)
        record["isotope"] = isotope
        flag = _isotope_flag(prior, shows_deuterium=_shows_deuterium(measured))
        if flag:
            record["isotope_flag"] = flag
        ambiguity = isotope_ambiguity(measured, prior)
        if ambiguity is not None:
            record["isotope_ambiguity"] = ambiguity
            ambiguities.append({"shot_number": shot_number, "cube": path.name, **ambiguity})

    coverage_nm = (min(coverage), max(coverage)) if coverage else None
    # One shot the calendar contradicts is enough to withhold the interval: the
    # cubes were audited as one epoch, so if any of them was measured against
    # references its own date forbids, the epoch's authorization is not sound.
    forcing = next(
        (item for item in ambiguities if item.get("excludes_deuterium")),
        ambiguities[0] if ambiguities else None,
    )
    verdict, summary = verdict_from_evidence(
        per_line, coverage_nm=coverage_nm, isotope_ambiguity=forcing
    )
    per_shot = _per_shot_summary(per_line)
    interval_warning = detect_interval_split(per_shot)
    _label_shot_groups(per_shot)
    summary["sampled_cubes"] = len(sampled)
    summary["skipped_cubes"] = len(skipped)
    summary["blended_lines_skipped"] = sum(
        1
        for item in per_line
        if item.get("status") == "skipped"
        and item.get("blended")
        and not item.get("isotope_excluded")
    )
    summary["isotope_excluded_lines"] = sum(
        1 for item in per_line if item.get("isotope_excluded")
    )
    summary["isotope_tags"] = _isotope_tally(sampled)
    summary["isotope_flagged_shots"] = [
        str(item["shot_number"]) for item in sampled if item.get("isotope_flag")
    ]
    summary["isotope_ambiguous_shots"] = [
        str(item["shot_number"]) for item in ambiguities
    ]

    order_corrections: list[dict[str, Any]] = []
    repair_steps: list[dict[str, str]] = []
    repair = ""
    refined_preview = ""
    if verdict == "shifted":
        shift_px = float(summary["median_shift_px"])
        order_corrections = _order_corrections(geometries, shift_px)
        root = Path(calibrations_root)
        base_ids = sorted(snapshot_ids)
        refined_is_preview = True
        if len(base_ids) == 1 and (root / base_ids[0]).is_dir():
            refined_preview = _next_refinement_id(base_ids[0], root)
            refined_is_preview = False
        else:
            refined_preview = f"{base_ids[0] if base_ids else 'snapshot'}-rN"
        evidence_target = Path(evidence_path) if evidence_path else Path("drift-evidence.json")
        repair_steps = _repair_steps(
            evidence_path=evidence_target,
            calibrations_root=root,
            refined_id=refined_preview,
            refined_is_preview=refined_is_preview,
            shift_px=shift_px,
            cube_directories=sorted({path.parent for path in selected}),
        )
        repair = repair_steps[0]["command"]

    payload: dict[str, Any] = {
        "schema": DRIFT_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample_rule": {
            "every": every,
            "shots": sorted(shots or set()),
            "date_from": date_from or "",
            "date_to": date_to or "",
            "date_attributes": list(DATE_ATTRIBUTES),
        },
        "snapshot_ids": sorted(snapshot_ids),
        "sampled_cubes": sampled,
        "skipped_cubes": skipped,
        "thresholds_px": {
            "alignment_maximum_residual": ALIGNMENT_MAX_PIXEL_RESIDUAL,
            "shift_consistency": SHIFT_CONSISTENCY_PX,
            "repair_limit": REPAIR_LIMIT_PX,
            "minimum_centroids": MINIMUM_CENTROIDS,
            "minimum_distinct_orders": MINIMUM_DISTINCT_ORDERS,
            "minimum_snr": MINIMUM_SNR,
            "verdict_snr_floor": VERDICT_SNR_FLOOR,
            "geometric_damage": GEOMETRIC_DAMAGE_PX,
            "isotope_degeneracy_window": ISOTOPE_DEGENERACY_WINDOW_PX,
            "plasma_frame_sigma": PLASMA_FRAME_SIGMA,
        },
        "isotope_prior": _isotope_prior_provenance(),
        "verdict": verdict,
        "summary": summary,
        "per_shot": per_shot,
        "lines": per_line,
        "order_corrections": order_corrections,
        "repair_command": repair,
        "repair_commands": repair_steps,
    }
    if refined_preview:
        payload["refined_snapshot_id"] = refined_preview
    if interval_warning:
        payload["interval_warning"] = interval_warning
    if ambiguities:
        payload["isotope_ambiguity"] = ambiguities
    if verdict == ERA_MISASSIGNED_VERDICT and forcing is not None:
        # No raw data is needed and no shift repairs this: the cube is fine and
        # the calibration it was processed against is the wrong one.  The fix is
        # a *cross-era* recalibration, so it is the full 'echelle recal-cube'
        # and never the --wavelength-only form the shifted verdict composes:
        # that flag exists for a refinement of the same base snapshot, which
        # keeps its sphere pair and integral curve, and another era's snapshot
        # carries neither of those.
        payload["verdict_advice"] = (
            f"{forcing['finding']} — recalibrate shot(s) "
            + ", ".join(sorted({str(item["shot_number"]) for item in ambiguities}))
            + " onto the calibration snapshot their own acquisition date selects, then audit "
            "again. This is a full 'echelle recal-cube <cube> -o <out> --new-snapshot "
            "<other-era-snapshot>' across eras, NOT the --wavelength-only form: the other "
            "era's sphere pair and integral curve must be re-derived too. Accepting a shift "
            "from this evidence instead would bake one isotope's worth of error "
            f"({forcing['isotope_separation_px']:+.1f} px) into a snapshot."
        )
    if verdict == "misaligned-beyond-repair":
        payload["data_requirement"] = _drive_requirement(
            sorted({str(item["shot_number"]) for item in per_shot}), catalog
        )
        payload["verdict_advice"] = str(summary.get("verdict_reason", ""))
    return payload


def _order_corrections(
    geometries: dict[str, CubeGeometry], shift_px: float
) -> list[dict[str, Any]]:
    """Predict the wavelength correction the refinement will apply per order.

    Moving the table's anchors by ``shift_px`` makes the refitted solution
    ``p_B(x) = p_A(x - shift_px)``, so this is exactly the change
    ``recal-cube`` re-evaluates from the refinement's own polynomials.
    """

    merged: dict[int, tuple[CubeGeometry, tuple[float, float]]] = {}
    for geometry in geometries.values():
        for order, bounds in geometry.pixel_range.items():
            if order not in merged:
                merged[order] = (geometry, bounds)
    corrections = []
    for order in sorted(merged):
        geometry, (low, high) = merged[order]
        reference = 0.5 * (low + high)
        before, _ = geometry.wavelength_at(order, reference)
        after, clamped = geometry.wavelength_at(order, reference - shift_px)
        dispersion, _ = geometry.dispersion(order, reference)
        corrections.append(
            {
                "order": int(order),
                "reference_pixel": float(reference),
                "pixel_range": [float(low), float(high)],
                "dispersion_nm_per_px": float(dispersion),
                "predicted_shift_nm": float(after - before),
                **({"extrapolated": True} if clamped else {}),
            }
        )
    return corrections


def write_drift_evidence(path: str | Path, payload: dict[str, Any]) -> Path:
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(f"drift evidence is immutable and already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return destination


# ---------------------------------------------------------------------------
# Refinement
# ---------------------------------------------------------------------------


def _shift_wavelength_table(
    source: Path,
    destination: Path,
    *,
    shift_px: float,
    metadata: list[tuple[str, str]],
) -> int:
    """Slide every anchor along the detector, keeping its wavelength.

    The anchor that produced wavelength ``lambda`` at pixel ``c`` is now seen at
    ``c + shift_px``, so the refit solution assigns ``lambda`` there and each
    order's own dispersion turns the one pixel shift into that order's
    wavelength correction.
    """

    rows = load_wavelength_table(source)
    if not rows:
        raise DriftError(f"wavelength table {source.name} holds no correctable rows")
    write_wavelength_table(shift_lines_in_pixels(rows, shift_px), destination, metadata=metadata)
    return len(rows)


def create_refinement_snapshot(
    evidence_path: str | Path,
    *,
    calibrations_root: str | Path,
    accepted_shift_px: float,
) -> tuple[Snapshot, Path]:
    """Accept a shifted verdict and emit a new immutable ``-rN`` snapshot."""

    evidence_file = Path(evidence_path)
    accepted_path = evidence_file.with_suffix(".accepted.json")
    if accepted_path.exists():
        raise FileExistsError(f"accepted drift record already exists: {accepted_path}")
    evidence = json.loads(evidence_file.read_text(encoding="utf-8"))
    if evidence.get("schema") == DRIFT_SCHEMA_V1:
        raise DriftError(
            "this evidence states a constant wavelength offset, which no longer describes a "
            "rigid detector shift; re-run 'echelle drift audit' to measure it in pixels"
        )
    if evidence.get("schema") != DRIFT_SCHEMA or evidence.get("verdict") != "shifted":
        raise DriftError("only a shifted sampled verdict can create a refinement snapshot")
    measured = float(evidence.get("summary", {}).get("median_shift_px", np.nan))
    # Exact acknowledgement, with room only for the decimal text the composed
    # command carries: 1e-9 relative on a pixel is a billionth of a pixel.
    if not np.isfinite(measured) or not np.isclose(
        accepted_shift_px, measured, rtol=1e-9, atol=1e-9
    ):
        raise DriftError("--accept-shift must exactly acknowledge the sampled median shift")
    snapshot_ids = list(evidence.get("snapshot_ids", []))
    if len(snapshot_ids) != 1:
        raise DriftError("refinement evidence must refer to exactly one base snapshot")
    root = Path(calibrations_root)
    base = load_snapshot(root / snapshot_ids[0])
    refined_id = _next_refinement_id(base.snapshot_id, root)
    with tempfile.TemporaryDirectory(prefix="echelle-refinement-") as temporary:
        adjusted = Path(temporary) / "wavelength.txt"
        moved = _shift_wavelength_table(
            base.source_path("wavelength"),
            adjusted,
            shift_px=measured,
            metadata=[
                ("base snapshot", base.snapshot_id),
                ("evidence", evidence_file.name),
                ("evidence sha256", sha256_file(evidence_file)),
                ("model", "rigid detector translation; anchors moved, wavelengths unchanged"),
                ("applied pixel shift", f"{measured:+.4f} px"),
            ],
        )
        files = {
            role: (adjusted if role == "wavelength" else base.source_path(role))
            for role in REQUIRED_ROLES
        }
        lamps = [
            (artifact.label, base.path_for(artifact))
            for artifact in base.artifacts
            if artifact.role == "lamp"
        ]
        # A refinement only re-solves the wavelength table.  If the base points
        # at its raw frames rather than holding them, the refinement points at
        # the same ones instead of hoarding a second copy.
        references_raw = any(artifact.is_reference for artifact in base.artifacts)
        snapshot = create_snapshot(
            root,
            snapshot_id=refined_id,
            detector=base.detector,
            files=files,
            lamps=base.lamps,
            lamp_files=lamps,
            notes=f"Accepted science-line refinement from {evidence_file.name}.",
            base_snapshot=base.snapshot_id,
            reference_raw=references_raw,
            validity=base.manifest.get("validity"),
            alignment={
                "method": "sampled Balmer/Fulcher rigid detector shift",
                "observed_shift_px": measured,
                "applied_pixel_shift_px": measured,
                "corrected_rows": moved,
                "evidence_sha256": sha256_file(evidence_file),
            },
            qc={"sampled_verdict": "shifted", "accepted": True},
        )
    accepted = {
        **evidence,
        "accepted": True,
        "accepted_snapshot_id": snapshot.snapshot_id,
        # Accepting a shift condemns the audited snapshot: only the refinement
        # carries the corrected wavelength table, so only the refinement is
        # authorized for bulk processing.
        "snapshot_ids": [snapshot.snapshot_id],
        "base_snapshot_ids": sorted(snapshot_ids),
        "base_evidence_sha256": sha256_file(evidence_file),
    }
    write_drift_evidence(accepted_path, accepted)
    return snapshot, accepted_path


# ---------------------------------------------------------------------------
# Gate
# ---------------------------------------------------------------------------


def _require_authorized_snapshots(selected: set[str], authorized: set[str]) -> None:
    if not selected <= authorized:
        missing = ", ".join(sorted(selected - authorized))
        raise DriftError(f"drift evidence does not cover selected snapshot(s): {missing}")


def require_sampled_verdict(path: str | Path, snapshot_ids: set[str]) -> dict[str, Any]:
    """Validate the bulk-processing prerequisite without treating insufficiency as aligned."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") not in READABLE_DRIFT_SCHEMAS:
        raise DriftError("unsupported drift evidence schema")
    verdict = payload.get("verdict")
    if verdict == "aligned":
        _require_authorized_snapshots(
            snapshot_ids, set(map(str, payload.get("snapshot_ids", [])))
        )
        return payload
    if verdict == "shifted" and payload.get("accepted"):
        refined = str(payload.get("accepted_snapshot_id") or "")
        recorded = set(map(str, payload.get("snapshot_ids", [])))
        # Records written before the refinement narrowed its own authorization
        # list both ids; everything that is not the refinement is superseded.
        superseded = set(map(str, payload.get("base_snapshot_ids") or [])) or (
            recorded - {refined}
        )
        condemned = sorted(snapshot_ids & superseded)
        if condemned:
            raise DriftError(
                f"registry still selects {', '.join(condemned)}; the accepted correction "
                f"produced {refined} — update the registry validity to point at the "
                "refined snapshot"
            )
        _require_authorized_snapshots(snapshot_ids, {refined} if refined else set())
        return payload
    if verdict == "insufficient-data":
        raise DriftError("bulk processing refused: sampled verdict is insufficient-data")
    if verdict == ERA_MISASSIGNED_VERDICT:
        raise DriftError(
            "bulk processing refused: the audited lines read as deuterium where the "
            "calendar expects hydrogen, and the hydrogen reading of the same centroids is "
            "one H/D-sized shift — recalibrate these cubes against the snapshot their own "
            "dates select and audit again"
        )
    raise DriftError(f"bulk processing refused: sampled verdict is {verdict}")
