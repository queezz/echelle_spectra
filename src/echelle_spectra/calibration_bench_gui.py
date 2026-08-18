"""Separate pyqtgraph live calibration bench for ``echelle-calib``."""

from __future__ import annotations

import argparse
import html
import platform
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass, replace
from datetime import date, datetime
from pathlib import Path, PurePath

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets

from . import __version__
from .calibration_bench import (
    AlignmentState,
    BenchFrame,
    CalibrationBenchSession,
    FileLoadState,
    FrameLoader,
    StableFileState,
    StableSifWatcher,
)
from .calibration_campaign import (
    KNOWN_LAMP_NAMES,
    PREVIOUS_CAMPAIGN_LAMPS,
    SELF_COMPARISON_NOTE,
    CalibrationCampaignSession,
    ChecklistState,
    ComparisonState,
    ExposureState,
    ExposureTriage,
    LampReferenceSet,
    MeasurementRole,
    TomlState,
    background_sibling,
    catalog_mismatch_warning,
    default_validity,
    expected_lines_for_order,
    lamp_reference_set,
    triage_for_role,
)
from .snapshot import SnapshotError
from .tools.calibration_alignment import load_wavelength_table, table_vetting
from .tools.pattern_extraction import (
    DEFAULT_SEARCH_RADIUS_PX,
    extract_pattern_from_sphere,
    pattern_row_offsets,
    subtract_background,
)

_PACKAGE_DIR = Path(__file__).parent
_CALIBRATION_DIR = _PACKAGE_DIR / "resources" / "calibration_files"
_DEFAULT_PATTERN = _CALIBRATION_DIR / "pattern_CMOS_20250926.txt"
_DEFAULT_WAVELENGTH = (
    _CALIBRATION_DIR
    / "alignments"
    / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
)
_DEFAULT_INTEGRAL = _CALIBRATION_DIR / "integrating_sphere.txt"
#: The packaged 2024-03-05 sphere pair, used as the previous campaign unless
#: the operator names another.  It is a *copy* of frames from a real folder,
#: which is why the comparison checks content and says so when the folder
#: being calibrated is that same folder.
_DEFAULT_PREVIOUS_SPHERE = _CALIBRATION_DIR / "sphere_cmos_20240305.sif"
_DEFAULT_PREVIOUS_SPHERE_BACKGROUND = (
    _CALIBRATION_DIR / "sphere_cmos_20240305_bkg.sif"
)

_SPHERE_FACTORS_TITLE = "Absolute calibration factors"
_SPHERE_FACTORS_EXPLANATION = (
    "The sphere signal minus its background, divided by the integrating "
    "sphere's known radiance, gives the factor curve that turns counts into "
    "W m⁻² sr⁻¹ nm⁻¹. The median ratio compares this campaign's curve with "
    "the previous one: near 1 means the instrument's response has not moved, "
    "and a large departure is either real ageing or an exposure-normalisation "
    "mismatch worth chasing before the trip. Only the sphere pair is needed — "
    "no lamp."
)

#: The one folder everything the bench generates lives under, wherever its
#: roots end up — owner, 2026-08-17: "It's kinder to keep all generated stuff
#: in one tidy folder."
SNAPSHOT_ROOT_NAME = "calibrations"

#: And its own subfolder for the reviewable settings bundles, so the snapshot
#: root holds snapshot identities and exactly one other child.  A sibling
#: ``calibration-configs\`` beside it was a second generated folder to explain,
#: to find, and to carry, for no gain over a subfolder of the first.
CONFIG_ROOT_NAME = "configs"

#: And its own subfolder for patterns the bench itself extracts, beside the
#: settings bundles rather than loose in the snapshot root: a pattern is an
#: input the campaign was rebased on, not a snapshot.
PATTERN_ROOT_NAME = "patterns"


def packaged_calibration_files(suffix: str = ".txt") -> tuple[Path, ...]:
    """The packaged tables an operator may name by their bare filename."""

    if not _CALIBRATION_DIR.is_dir():  # pragma: no cover - a broken installation
        return ()
    return tuple(
        sorted(
            (path for path in _CALIBRATION_DIR.rglob(f"*{suffix}") if path.is_file()),
            key=lambda path: path.name.casefold(),
        )
    )


def resolve_calibration_file(value: str | PurePath, kind: str) -> Path:
    """Resolve a table argument as a real path first, then as a packaged name.

    The packaged names *are* the documented vocabulary — the README's own
    examples say ``pattern_CMOS_20240305.txt`` — and typing exactly that got the
    owner "pattern file not found", because nothing but a literal path was ever
    tried.  A path on disk still wins, so nothing that worked stops working, and
    the refusal lists the handful of packaged names rather than leaving the
    vocabulary to be guessed.
    """

    source = Path(value)
    if source.is_file():
        return source.resolve()
    packaged = packaged_calibration_files(source.suffix or ".txt")
    if source.name == str(source):
        match = next((path for path in packaged if path.name == source.name), None)
        if match is not None:
            return match.resolve()
    # The vocabulary this argument is about, not every packaged table: a
    # refusal that lists the wavelength tables under --pattern teaches nothing.
    offered = [path for path in packaged if kind in path.name.casefold()] or list(packaged)
    names = ", ".join(path.name for path in offered)
    raise SystemExit(
        f"{kind} file not found: {source}\n"
        f"Packaged {kind} tables you can name directly: {names or 'none'}"
    )


def absolute_root(path: str | PurePath) -> Path:
    """Return *path* as an absolute path without asking the filesystem anything.

    ``resolve()`` is the obvious call and the wrong one here: it asks the
    operating system, and asking about ``\\\\server\\share\\…`` costs a network
    round trip that can hang when the share is away.  A path that is already
    absolute is already the honest answer — only a relative one needs the
    working directory written in front of it.  A pure path is handed back as it
    came, so a Windows UNC path can be reasoned about on any platform.
    """

    base = path if isinstance(path, PurePath) else Path(path)
    if base.is_absolute():
        return base  # type: ignore[return-value]
    return Path.cwd() / base


def default_bench_roots(folder: str | PurePath | None = None) -> tuple[Path, Path]:
    """Where snapshots and generated settings land when no flag says otherwise.

    The bench is launched *at* a calibration folder, and that folder is where
    the operator goes looking for what the bench wrote.  Deriving both roots
    from it is what stopped a run against a NAS folder from scattering its
    TOMLs into whatever directory the shortcut happened to start in.  The
    computed calibration belongs beside the frames it was computed from —
    owner, 2026-08-17: "our metadata inside the folder with real lamp frames.
    Side by side, and that folder holds it all. Raw and calculated stuff."
    With no folder to derive from there is nothing better than the working
    directory, so that is what it falls back to — but the Save tab says so out
    loud either way.

    The settings bundles land *inside* the snapshot root rather than in a
    sibling folder, so the calibration folder gains one generated child and not
    two.  Snapshot enumeration is unaffected: everything that reads a
    calibrations root looks for ``<child>/snapshot.toml`` or resolves a
    snapshot id it was given, so a ``configs`` child is simply not a snapshot.

    Only ``pathlib`` joins are used, which is what keeps a UNC root intact: a
    ``Path`` appends below ``\\\\server\\share`` without ever losing the two
    leading separators that string surgery eats.
    """

    base = absolute_root(Path.cwd() if folder is None else folder)
    snapshots = base / SNAPSHOT_ROOT_NAME
    return snapshots, snapshots / CONFIG_ROOT_NAME


def configure_folder_dialog(dialog: QtWidgets.QFileDialog) -> None:
    """Make a folder picker show the files it is not offering to pick.

    Windows' native folder picker hides files outright: a calibration folder
    holding six SIFs reads "No items match your search", so the only way to
    tell one dated folder from another was to open it, look, close it, and open
    the right one — owner, 2026-08-18: "we should show contents, but make them
    gray. Otherwise I have to open the same folder twice."

    Qt's own dialog does exactly the right thing in ``Directory`` mode with
    ``ShowDirsOnly`` switched *off*: the files are listed and greyed, so the
    folder's contents are evidence without ever becoming the answer.  The
    native dialog cannot be told this, which is why this one asks for Qt's by
    name.  Only folder pickers have the defect — a dialog that picks files
    shows files already — so the pattern and previous-pair pickers keep the
    native ``getOpenFileName`` they have always used.
    """

    dialog.setFileMode(QtWidgets.QFileDialog.Directory)
    dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly, False)
    dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
    dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptOpen)


def choose_calibration_folder(parent, start_dir: str | PurePath) -> str:
    """Ask which calibration folder to open, through one patchable seam.

    A module-level function rather than an inline dialog, for the same reason
    the main GUI's ``choose_snapshot_folder`` is one: a test can answer this
    without a real modal appearing off-screen and hanging the run.
    """

    dialog = QtWidgets.QFileDialog(parent, "Open calibration folder", str(start_dir))
    configure_folder_dialog(dialog)
    if not dialog.exec_():
        return ""
    chosen = dialog.selectedFiles()
    return chosen[0] if chosen else ""


#: A calibration folder's name begins with the day its frames were taken —
#: ``20250926_calib``.  Anchored, and refusing a ninth digit, so a shot number
#: that merely starts with eight digits is never mistaken for a date.
_ACQUISITION_DATE_IN_NAME = re.compile(r"^(\d{4})(\d{2})(\d{2})(?!\d)")


def acquisition_date_from_name(folder: str | PurePath | None) -> date | None:
    """The acquisition date a calibration folder's own name leads with.

    ``20250926_calib`` is the calibration measured on 2025-09-26, and that is
    the date the snapshot identity is built from.  Only the folder's own name
    is read, never the path above it: a campaign tree may carry any number of
    dated components, and only the folder the bench was pointed at is a claim
    about *these* frames.
    """

    if folder is None:
        return None
    name = folder if isinstance(folder, str) else PurePath(folder).name
    match = _ACQUISITION_DATE_IN_NAME.match(str(name).strip())
    if match is None:
        return None
    try:
        return date(int(match.group(1)), int(match.group(2)), int(match.group(3)))
    except ValueError:  # 20259999_calib is not a date, it is a typo
        return None


def acquisition_date_from_metadata(metadata: object) -> date | None:
    """The day an Andor SIF header says its acquisition began.

    ``ExperimentTime`` is the Unix timestamp the acquisition started at, the
    same field the cube exporter records as ``t_start``.  It is read as a
    *local* calendar day on purpose: the identity convention names the day at
    the instrument, which is the same day the folder beside it was named for.
    """

    if not isinstance(metadata, dict):
        return None
    raw = metadata.get("ExperimentTime")
    if raw in (None, ""):
        return None
    try:
        seconds = int(raw)
    except (TypeError, ValueError):
        return None
    if seconds <= 0:
        return None
    try:
        return datetime.fromtimestamp(seconds).date()
    except (OSError, OverflowError, ValueError):
        return None


def snapshot_identity(day: date, detector: str = "cmos") -> str:
    """``YYYYMMDD_<detector>`` — the identity every epoch in the registry uses.

    Dated by acquisition, never by computation: ``20240305_cmos`` and
    ``20250926_cmos`` name the days those calibrations were *measured*, and a
    snapshot computed months later still belongs to its own frames.  (Owner,
    2026-08-17: "calibration IS dated by the day the images were taken. When we
    calculate is our bookkeeping. Not physics.")
    """

    return f"{day:%Y%m%d}_{str(detector).strip() or 'cmos'}"


def open_containing_folder(path: str | Path) -> bool:
    """Show one saved thing's folder in the system's own file manager.

    ``QUrl.fromLocalFile`` is what makes this portable and what makes it work
    for a UNC path — it renders ``\\\\server\\share\\x`` as a ``file:`` URL
    Explorer opens.  A directory is shown as itself; anything else is shown by
    its parent, because the folder is what the operator wanted to see.
    """

    target = Path(path)
    folder = target if target.is_dir() else target.parent
    return bool(QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(folder))))


#: Role labels are deliberately short.  The Role control shows only the role;
#: nothing about a file's *state* is ever written into it, because a combo box
#: elides its text and the clipped half is always the half that mattered.
_ROLE_CHOICES: tuple[tuple[str, MeasurementRole | None], ...] = (
    ("none", None),
    ("Sphere", MeasurementRole.SPHERE),
    ("Sphere bg", MeasurementRole.SPHERE_BACKGROUND),
    ("Lamp", MeasurementRole.LAMP),
    ("Lamp bg", MeasurementRole.LAMP_BACKGROUND),
    ("Other", MeasurementRole.OTHER),
)

_LAMP_ROLES = (MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND)

_TRIAGE_COLORS = {
    ExposureState.GOOD: "#70d6ae",
    ExposureState.DIM: "#ffb86b",
    ExposureState.SATURATED: "#ff8f8f",
    ExposureState.NO_DATA: "#93a8b8",
}

def _record_exposure_s(record) -> float | None:
    """The exposure time one loaded frame's header gave, if it gave one.

    Read in one place because every surface that asks a verdict for its role
    has to hand it the same number: the saturated readings quote an exposure to
    try next, and a suggestion computed from a missing header is no suggestion.
    """

    exposure = getattr(record, "exposure", None)
    return None if exposure is None else exposure.exposure_s


#: Colour for a correct background: a calm slate blue that belongs to no
#: exposure state.  A dark frame that is dark is not an alarm (amber), not a
#: blockage (red) and not the same reading as a healthy signal (the GOOD
#: green) — on the bench's dark ground the two greens were indistinguishable,
#: and "this one is dark on purpose" is exactly what the operator needs to see
#: without being called over to it (owner, 2026-08-18).
_BACKGROUND_COLOR = "#8aa6e0"

#: Colour for a control showing a filename suggestion nobody has confirmed.
_SUGGESTED_COLOR = "#ffb86b"

#: The badge that marks an unconfirmed suggestion.  It lives in the file cell
#: and in the control's colour — never inside the combo's own text.
_SUGGESTED_BADGE = "SUGGESTED"

#: What the top-end panel says when the last 10% of the range holds nothing.
#: In words, because an empty log histogram draws a solid block that reads
#: like a full one (F17 item 3).
_TOP_END_EMPTY = "No pixels within 10% of full scale."

_SUGGESTED_TOOLTIP = (
    "This role is only a guess read from the filename. The bench has NOT been "
    "given it. Pick the role in this list — or press “Confirm suggested "
    "roles” — before the procedure and the factor computation can use the file."
)

# ----------------------------------------------------------------------
# Loudness is a budget: one loud element per view, everything else body size
# ----------------------------------------------------------------------
#: Ordinary bench prose and every secondary reading, in points.  Point sizes
#: rather than pixels so a lower-resolution display scales the whole surface
#: with the platform font instead of shrinking it.
BENCH_BODY_POINT_SIZE = 11.0

#: The one loud element: the verdict headline the next action is read from.
BENCH_HEADLINE_POINT_SIZE = 17.0

#: A tooltip is one short line.  Everything longer belongs in the Why dock,
#: where it can be read calmly instead of floating over the controls.
BENCH_TOOLTIP_LIMIT = 88

#: The raw-counts histogram's floor, in lines of body text.  It is the primary
#: reading of the triage view and never shrinks to its neighbour's size.
BENCH_HISTOGRAM_LINES = 8

#: The near-saturation strip's ceiling, in the same unit.  Deliberately a
#: fraction of the primary histogram's floor: the top end is the number the
#: lamp is adjusted by, and it is read *after* the raw counts, never instead
#: of them (F21 item 4b).
BENCH_TOP_END_LINES = 4

#: The anchor-residual strip's ceiling, in the same unit and for exactly the
#: same reason as the top-end strip above.  The residuals live in the anchors
#: panel now (owner, 2026-08-17: "what about moving the anchor plot into the
#: anchors tab?"), and the table is what that panel is *for*: the strip is the
#: shape of the residuals at a glance, read beside the rows, never instead of
#: them.  Capped, so it can never grow into a second working surface.
BENCH_RESIDUAL_LINES = 5

#: And what it shrinks to when the panel is short.  The subordination here has
#: to be stronger than the triage strip's, because that strip shares its column
#: with a plot and this one shares its panel with a *table* whose workable row
#: count is a law of the rail (``EXPECTED_LINE_ROWS``).  So the strip takes
#: what is left over after the table is workable, between these two bounds:
#: the shortest strip that still carries a row of wavelengths under its bars,
#: and the full strip once the window is maximized or the rail is dragged
#: wider.  A fixed cap at the ceiling would have cost the anchor table four of
#: its six rows at the size the bench opens at.
BENCH_RESIDUAL_FLOOR_LINES = 2

#: The point size the strip's wavelengths are drawn at.  Below the body text on
#: purpose and the only thing on the bench that is: these are tick labels on a
#: strip a rail wide, and the alternative to a smaller font is fewer of them.
#: The operator reads a wavelength off the axis rather than having to click
#: (owner, 2026-08-17: "dropping wavelength hints on the plots? Why? Quite
#: usable. Useful.").
BENCH_RESIDUAL_TICK_POINT_SIZE = 7.0

#: Clear space either side of a wavelength label before the next one may be
#: drawn.  What separates "thinned" from "a smear": labels are dropped until
#: the ones that remain have this much room.
BENCH_RESIDUAL_TICK_GAP = 14

#: How many summary rows the triage table shows before it scrolls.  Enough for
#: a whole campaign folder; past that the reading is a scroll either way.
BENCH_SUMMARY_ROWS = 8

#: The measure prose is set to in the Why dock, in characters.  Typography's
#: own number, not a pixel one: a line much past this is measurably harder to
#: read because the eye loses its place on the way back to the left margin.
#: The dock spans the whole window, so without a cap an explanation wrapped at
#: whatever the window happened to be — two thousand pixels on the owner's
#: screen ("all long and breaky").
BENCH_READING_MEASURE_CHARS = 96

#: How tall the Why dock opens: a title and the five or six wrapped lines a
#: typical verdict costs, measured rather than guessed.  It used to open at
#: four lines of the platform's text and crop almost everything written into
#: it.  This is a wish, not a demand — see ``_details_dock_ceiling_lines``.
BENCH_DETAILS_LINES = 12

#: And the floor it may never be squeezed below, dragged or not.
BENCH_DETAILS_FLOOR_LINES = 3


def bench_point_sizes(base_point_size: float | None = None) -> tuple[float, float]:
    """Return the (body, headline) point sizes for this platform's font.

    Two tiers, not three: the verdict shouts and everything else is body text.
    The floors are absolute, so a small platform font never shrinks the bench
    below legibility and a large one scales both tiers with it together.
    """

    if base_point_size is None or base_point_size <= 0:
        application = QtWidgets.QApplication.instance()
        base_point_size = (
            application.font().pointSizeF() if application is not None else 0.0
        )
    if base_point_size <= 0:  # a pixel-sized platform font reports no points
        base_point_size = 9.0
    body = max(BENCH_BODY_POINT_SIZE, float(base_point_size) * 1.15)
    return body, max(BENCH_HEADLINE_POINT_SIZE, body * 1.55)


#: The bench's geometry, in lines of the platform's own text rather than in
#: pixels measured on one designer's display.  A text line is the one unit that
#: moves with the font, the display's DPI and the desktop's scaling factor at
#: once, which is what "DPI-independent" has to mean here: 1500 px is a
#: comfortable window at 100% scaling and two thirds of a *screen* at 150%.
BENCH_PREFERRED_LINES = (79, 48)

#: The floor the two-rail layout is built to satisfy: both rails working, the
#: readings strip in view, and the plots reduced but present.  Narrower than
#: the preference but not shorter than it, because the two behave differently:
#: taking width away shrinks the plots and the tables lose columns they can
#: scroll to, while taking height away puts a scrollbar on the controls — and a
#: scrolled controls column is the defect this packet exists to end.
BENCH_FLOOR_LINES = (72, 48)

#: What a window costs its screen beyond its own client area — frame, title
#: bar, taskbar — expressed in the same unit rather than in a pixel guess.
_SCREEN_MARGIN_LINES = (2, 3)

#: The bench's two rails and its middle, as shares of whatever width the window
#: has (F20, owner: "TWO RAILS").  Controls on the left, the two working tables
#: on the right, plots between them — deliberately the smallest of the three,
#: because the plots are the element that can be zoomed and the tables are the
#: ones being read.  These are shares and stretch factors rather than pixel
#: cuts precisely so that maximizing widens all three instead of spending the
#: whole gain on the plots and leaving the rails where they were.
_ROOT_SHARES = (0.28, 0.38, 0.34)
_ROOT_STRETCH = (28, 38, 34)

#: The right rail's own cut.  Near enough even; the expected-line table gets
#: the larger half because it carries a wrapping header the anchor table has
#: no equivalent of, so an even split leaves it the fewer readable rows.
_TABLES_SHARES = (0.48, 0.52)
_TABLES_STRETCH = (1, 1)


def bench_layout_unit(font: QtGui.QFont | None = None) -> int:
    """One line of the bench's own body text: the unit every size is quoted in.

    Everything the bench measures itself against — window sizes, screen
    margins, table floors — is a multiple of this.  A larger platform font, a
    higher DPI or a 150% desktop scale all move it together, so the geometry
    stays the same *layout* rather than the same number of pixels.
    """

    if font is None:
        application = QtWidgets.QApplication.instance()
        font = QtGui.QFont(application.font()) if application is not None else QtGui.QFont()
        font.setPointSizeF(bench_point_sizes()[0])
    return max(8, QtGui.QFontMetrics(font).height())


def bench_minimum_size(
    available: QtCore.QSize | None = None, unit: int | None = None
) -> tuple[int, int]:
    """The smallest usable bench — but never larger than the screen it is on.

    F18's law ("if the window is unusable at a size, never open at that size")
    was written as a pixel floor, and a pixel floor is exactly what a scaled
    display breaks: a 1920x1080 screen at 150% offers 1280x720 *logical*
    pixels, so a 1300x880 floor asks for a window taller and wider than the
    whole desktop.  Windows then hands back a window smaller than the floor and
    every panel inside it is crushed — the owner's cramped first paint.  The
    floor therefore yields to the screen: on a display that cannot hold it, the
    honest answer is the display, and the two-rail layout degrades by showing
    fewer table rows rather than by hiding controls behind a scrollbar.
    """

    unit = unit or bench_layout_unit()
    floor = (BENCH_FLOOR_LINES[0] * unit, BENCH_FLOOR_LINES[1] * unit)
    available = _available_screen_size(available)
    if available is None:
        return floor
    return (
        int(min(floor[0], max(unit, available.width() - _SCREEN_MARGIN_LINES[0] * unit))),
        int(min(floor[1], max(unit, available.height() - _SCREEN_MARGIN_LINES[1] * unit))),
    )


def bench_default_geometry(
    available: QtCore.QSize | None = None, unit: int | None = None
) -> tuple[int, int]:
    """Return a default window size that fits the screen it will open on.

    The preferred size is what the bench wants; the screen is what it gets, and
    the floor is what it settles for — in that order, so the result can never
    exceed the desktop the window has to be painted on.
    """

    unit = unit or bench_layout_unit()
    preferred = (BENCH_PREFERRED_LINES[0] * unit, BENCH_PREFERRED_LINES[1] * unit)
    floor = bench_minimum_size(available, unit)
    available = _available_screen_size(available)
    if available is None:
        return (int(max(preferred[0], floor[0])), int(max(preferred[1], floor[1])))
    width = min(preferred[0], max(1, available.width() - _SCREEN_MARGIN_LINES[0] * unit))
    height = min(preferred[1], max(1, available.height() - _SCREEN_MARGIN_LINES[1] * unit))
    return int(max(width, floor[0])), int(max(height, floor[1]))


def _available_screen_size(available: QtCore.QSize | None) -> QtCore.QSize | None:
    """The desktop's usable *logical* size, or ``None`` when there is no screen."""

    if available is None:
        screen = QtWidgets.QApplication.primaryScreen()
        available = screen.availableGeometry().size() if screen is not None else None
    if available is None or available.width() <= 0 or available.height() <= 0:
        return None
    return available

#: The narrowest a reading panel is ever asked to be before the band folds it
#: onto its own row instead.  Prose wraps happily at this width.
_STATUS_PANEL_FLOOR = 210


def one_line(text: str, limit: int = BENCH_TOOLTIP_LIMIT) -> str:
    """Reduce an explanation to the one short line a tooltip may carry."""

    collapsed = " ".join(str(text).split())
    if not collapsed:
        return ""
    sentence = collapsed.split(". ")[0].rstrip(".")
    if len(sentence) + 1 <= limit:
        return f"{sentence}." if sentence else ""
    return sentence[: max(1, limit - 1)].rstrip() + "…"


def role_combo_minimum_width(combo) -> int:
    """The width at which *combo* can show its longest entry unelided.

    A role is state the operator must read exactly, so the control is sized to
    its content rather than left to shrink into an ellipsis at some width
    nobody predicted.
    """

    metrics = combo.fontMetrics()
    widest = max(
        (metrics.horizontalAdvance(combo.itemText(index)) for index in range(combo.count())),
        default=0,
    )
    style = combo.style()
    options = QtWidgets.QStyleOptionComboBox()
    options.initFrom(combo)
    arrow = style.subControlRect(
        QtWidgets.QStyle.CC_ComboBox,
        options,
        QtWidgets.QStyle.SC_ComboBoxArrow,
        combo,
    ).width()
    return int(widest + max(arrow, 18) + 16)


def _emphasise(widget, point_size: float, *, bold: bool = True) -> None:
    """Render one widget's text at *point_size*, keeping its family."""

    font = widget.font()
    font.setPointSizeF(float(point_size))
    font.setBold(bold)
    widget.setFont(font)


#: The bench's own Windows shell identity.  The main GUI has claimed one since
#: it was written (``echelle_spectra-<version>``); the bench never did, which is
#: the whole of the "the icon is gone" report: a Windows process with no
#: explicit AppUserModelID is filed under whatever launched it, so the taskbar
#: button showed the console-script stub's generic icon no matter what
#: ``setWindowIcon`` said.  A *different* id from the main GUI is also what
#: keeps the two windows in two taskbar groups (F14 item 5).
BENCH_APP_USER_MODEL_ID = f"echelle_spectra.calibration_bench-{__version__}"

#: Sizes Windows and X11 actually ask an icon for.  Handing the shell a real
#: pixmap at each one is what stops the 16 px title-bar copy from being a smear
#: of a 128 px drawing — a smear reads, at a glance, as no icon at all.
_ICON_SIZES = (16, 20, 24, 32, 48, 64, 128, 256)


def apply_windows_taskbar_identity(app_id: str = BENCH_APP_USER_MODEL_ID) -> str | None:
    """Claim this process's own Windows taskbar identity.

    Must run *before* the first window exists, which is why it is the first
    thing :func:`main` does.  Returns the id claimed, or ``None`` where there
    is no Windows shell to tell — the call is advisory everywhere else and
    never a reason for the bench not to start.
    """

    if platform.system() != "Windows":
        return None
    try:
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
    except (OSError, AttributeError):  # pragma: no cover - shell refusing is not fatal
        return None
    return app_id


def bench_window_icon(
    source: Path | None = None, *, size: int = 128
) -> QtGui.QIcon:
    """Derive a visibly distinct bench icon from the shared echelle graphic.

    The bench and the main GUI are two windows of one program, so they share
    one drawing rather than gaining a second piece of artwork.  The bench copy
    is tinted toward the bench's own cyan and carries a corner badge, which is
    what makes the two entries in a taskbar tellable apart at a glance.
    """

    origin = Path(source) if source is not None else (
        _PACKAGE_DIR / "resources" / "graphics" / "echelle.png"
    )
    base = QtGui.QPixmap(str(origin))
    if base.isNull():  # pragma: no cover - only when the resource is missing
        base = QtGui.QPixmap(size, size)
        base.fill(QtCore.Qt.transparent)
    canvas = base.scaled(
        size,
        size,
        QtCore.Qt.KeepAspectRatio,
        QtCore.Qt.SmoothTransformation,
    )

    tint = QtGui.QPixmap(canvas.size())
    tint.fill(QtCore.Qt.transparent)
    painter = QtGui.QPainter(tint)
    painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
    painter.drawPixmap(0, 0, canvas)
    # Tint only where the drawing already has pixels, so the transparent
    # background stays transparent.
    painter.setCompositionMode(QtGui.QPainter.CompositionMode_SourceAtop)
    painter.fillRect(tint.rect(), QtGui.QColor(64, 200, 255, 110))
    painter.setCompositionMode(QtGui.QPainter.CompositionMode_SourceOver)

    width = tint.width()
    height = tint.height()
    radius = max(6, int(round(min(width, height) * 0.30)))
    centre = QtCore.QPointF(width - radius - 1, height - radius - 1)
    painter.setPen(QtGui.QPen(QtGui.QColor("#0d1218"), max(1, radius // 6)))
    painter.setBrush(QtGui.QBrush(QtGui.QColor("#ffb020")))
    painter.drawEllipse(centre, float(radius), float(radius))
    # A caliper-like tick pair: the badge says "measuring bench", not "viewer".
    painter.setPen(QtGui.QPen(QtGui.QColor("#0d1218"), max(1, radius // 4)))
    painter.drawLine(
        QtCore.QPointF(centre.x() - radius * 0.45, centre.y() - radius * 0.35),
        QtCore.QPointF(centre.x() - radius * 0.45, centre.y() + radius * 0.35),
    )
    painter.drawLine(
        QtCore.QPointF(centre.x() + radius * 0.45, centre.y() - radius * 0.35),
        QtCore.QPointF(centre.x() + radius * 0.45, centre.y() + radius * 0.35),
    )
    painter.drawLine(
        QtCore.QPointF(centre.x() - radius * 0.45, centre.y()),
        QtCore.QPointF(centre.x() + radius * 0.45, centre.y()),
    )
    painter.end()

    icon = QtGui.QIcon(tint)
    for edge in _ICON_SIZES:
        if edge >= tint.width():
            continue
        icon.addPixmap(
            tint.scaled(
                edge, edge, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation
            )
        )
    return icon


#: Splitter cuts the operator made, kept for as long as this process lives, and
#: kept as *shares of the splitter* rather than as pixel counts.  A cut stored
#: in pixels is a cut tuned at one window size: restore it into a maximized
#: window and the left rail is still the 240 px it was in the small one, which
#: is precisely the collapse the owner had to drag back open before every use.
#: Session state rather than a settings file on purpose: a dragged handle is a
#: working preference for this sitting, not a decision to write to disk.
_SESSION_SPLITTER_SHARES: dict[str, tuple[float, ...]] = {}

#: The same working preference, for the docks — kept in *lines of body text*
#: rather than pixels, because that is the unit every other size here is quoted
#: in and it survives a font or DPI change the way a pixel count does not.
_SESSION_DOCK_LINES: dict[str, float] = {}


def forget_session_layout() -> None:
    """Drop every remembered splitter cut (used by tests and by a fresh run)."""

    _SESSION_SPLITTER_SHARES.clear()
    _SESSION_DOCK_LINES.clear()


class _ElidingLabel(QtWidgets.QLabel):
    """A one-line reading that shortens itself instead of wrapping or clipping.

    ``setText`` keeps the whole string; the widget only ever *draws* as much of
    it as its current width allows, shortened in the middle and offered in full
    as the tooltip.  Nothing here carries state the operator must read exactly
    — those controls are sized to their content instead (F16 item 1).
    """

    def __init__(self, text: str = "", parent=None) -> None:
        super().__init__(parent)
        self._full_text = ""
        self.setText(text)

    def setText(self, text: str) -> None:  # noqa: N802 - Qt naming
        self._full_text = str(text)
        self.setToolTip(self._full_text)
        self._redraw_elided()

    def text(self) -> str:
        """The whole reading, not the shortened one currently drawn.

        Everything that asks a label what it says wants the reading; only the
        painting is shortened.
        """

        return self._full_text

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt naming
        super().resizeEvent(event)
        self._redraw_elided()

    def minimumSizeHint(self) -> QtCore.QSize:  # noqa: N802 - Qt naming
        hint = super().minimumSizeHint()
        return QtCore.QSize(0, hint.height())

    def _redraw_elided(self) -> None:
        width = max(24, self.width() - 2)
        super().setText(
            self.fontMetrics().elidedText(
                self._full_text, QtCore.Qt.ElideMiddle, width
            )
        )


class FrameLoadThread(QtCore.QThread):
    """Read/extract one stable SIF away from the Qt event loop."""

    loaded = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str, str)

    def __init__(self, path: Path, loader: FrameLoader, parent=None) -> None:
        super().__init__(parent)
        self.path = path
        self.loader = loader

    def run(self) -> None:
        try:
            self.loaded.emit(self.loader(self.path))
        except Exception as exc:  # GUI boundary: report domain/IO failures in state
            self.failed.emit(str(self.path), str(exc))


@dataclass(frozen=True)
class PatternRebase:
    """A pattern file the bench is about to stand on, and the light it measured.

    Produced off the event loop — reading a sphere pair and fitting order traces
    over a 2560x2160 frame is seconds of work — and applied on the GUI thread,
    where the session, the reader and the campaign live.
    """

    path: Path
    #: Whether this bench fitted the table, as opposed to an operator naming it.
    extracted: bool = False
    sphere_path: Path | None = None
    #: The sphere signal's own detector image, so the band guard can be
    #: re-measured against the new pattern without opening the file again.
    sphere_image: np.ndarray | None = None
    background_path: Path | None = None
    #: Median row shift of the extracted table from the pattern it was fitted
    #: against, which is what the guard was complaining about.
    median_offset_rows: float | None = None
    caveat: str = ""


def extract_pattern_beside_campaign(
    *,
    loader,
    prior: np.ndarray,
    sphere_path: Path,
    background_path: Path | None,
    output_path: Path,
    search_radius_px: int = DEFAULT_SEARCH_RADIUS_PX,
) -> PatternRebase:
    """Fit a pattern to one sphere pair and write it beside the campaign.

    The fit is ``extract_pattern_from_sphere`` — the very function
    ``echelle-pattern`` calls — with the bench's current pattern as the prior,
    so a bench that fixes its own geometry and an operator who runs the CLI get
    the same table from the same frames.
    """

    sphere_image = np.asarray(loader(sphere_path).detector_image, dtype=float)
    caveats: list[str] = []
    if background_path is None:
        image = sphere_image
        caveats.append(
            "no sphere background carries a role, so the pattern was fitted to "
            "the signal alone; assign the background and extract again for the "
            "subtracted fit echelle-pattern performs"
        )
    else:
        image = subtract_background(
            sphere_image, np.asarray(loader(background_path).detector_image, dtype=float)
        )
    extraction = extract_pattern_from_sphere(
        image, prior, search_radius_px=search_radius_px
    )
    if not extraction.trial_succeeded:
        caveats.append(
            "no unguided trial found every order, so the prior alone guided the "
            "fit — check the traces over the detector image"
        )
    caveat = "; ".join(caveats)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_path, extraction.pattern, fmt="%d")
    return PatternRebase(
        path=output_path,
        extracted=True,
        sphere_path=sphere_path,
        sphere_image=sphere_image,
        background_path=background_path,
        median_offset_rows=float(
            np.median(pattern_row_offsets(extraction.pattern, prior))
        ),
        caveat=caveat,
    )


def read_pattern_beside_campaign(*, loader, chosen_path: Path, sphere_path: Path | None):
    """Adopt a pattern file somebody named, reading the sphere to re-judge it.

    The sphere frame is opened again for one reason: the band guard's verdict is
    a comparison, and a comparison against a pattern that is gone is not a
    reading.  Everything else about the chosen file is already on disk.
    """

    sphere_image = None
    if sphere_path is not None and loader is not None:
        sphere_image = np.asarray(loader(sphere_path).detector_image, dtype=float)
    return PatternRebase(
        path=Path(chosen_path), sphere_path=sphere_path, sphere_image=sphere_image
    )


class CampaignTaskThread(QtCore.QThread):
    """Run one potentially slow campaign operation away from the event loop."""

    completed = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, operation, parent=None) -> None:
        super().__init__(parent)
        self.operation = operation

    def run(self) -> None:
        try:
            self.completed.emit(self.operation())
        except Exception as exc:  # GUI boundary: domain state retains the detail
            self.failed.emit(str(exc))


class CalibrationBenchWindow(QtWidgets.QMainWindow):
    """Thin Qt adapter over :class:`CalibrationBenchSession`."""

    def __init__(
        self,
        session: CalibrationBenchSession,
        *,
        campaign: CalibrationCampaignSession | None = None,
        watcher: StableSifWatcher | None = None,
        loader: FrameLoader | None = None,
        folder: str | Path | None = None,
        output_root: str | Path = SNAPSHOT_ROOT_NAME,
        config_root: str | Path | None = None,
        snapshot_id: str = "",
        snapshot_date: date | None = None,
        snapshot_id_explicit: bool = False,
        detector: str = "cmos",
        base_snapshot: str = "",
        valid_from: date | None = None,
        poll_interval_ms: int = 1000,
        start_timer: bool = True,
    ) -> None:
        super().__init__()
        self.session = session
        self.campaign = campaign
        self.watcher = watcher
        self.loader = loader
        #: The calibration folder this session is about: what the roots hang
        #: off, what the identity is dated by, where the file dialogs open, and
        #: the one fact the Bench state strip could not previously state.
        #: ``None`` is a bench nobody has pointed at a folder yet.
        self.calibration_folder = None if folder is None else absolute_root(folder)
        # Absolute from here on, whatever the caller passed.  A bare name is a
        # promise about the working directory, and a display that repeats the
        # bare name keeps that promise secret — which is exactly how a campaign
        # got saved into the launcher's home folder instead of beside its data.
        self.output_root = absolute_root(output_root)
        # Derived, the settings bundles are a child of the snapshot root: one
        # generated folder per calibration folder, not two.  An explicitly
        # passed root is honoured wherever it points.
        self.config_root = (
            self.output_root / CONFIG_ROOT_NAME
            if config_root is None
            else absolute_root(config_root)
        )
        self.initial_snapshot_id = snapshot_id
        self.initial_detector = detector
        self.initial_base_snapshot = base_snapshot
        self.valid_from = valid_from
        #: The acquisition date the prefilled identity is dated by, and where
        #: it was read.  ``None`` means no acquisition date could be derived
        #: yet, so the identity carries today's date as a placeholder and the
        #: Save tab says so — the one case the operator must check by hand.
        self.snapshot_date = snapshot_date
        self._snapshot_date_source = (
            "the --snapshot-id you passed"
            if snapshot_id_explicit
            else ("the folder name the bench was launched at" if snapshot_date else "")
        )
        #: Whether the identity is somebody's decision rather than the bench's
        #: guess — a flag, or a hand edit.  A decision is never re-derived when
        #: frames land and say something else.
        self._snapshot_id_decided = bool(snapshot_id_explicit)
        self._load_thread: FrameLoadThread | None = None
        self._background_thread: FrameLoadThread | None = None
        self._campaign_thread: CampaignTaskThread | None = None
        #: The verb the next-step button carries.  It is set from the checklist
        #: on every refresh, but the button is built with a caption before any
        #: refresh runs and a campaign-less bench never gets one — so the
        #: caption's own action is what it starts with, never nothing.
        self._next_action = self._pick_files
        #: The snapshot identity an "already exists" refusal was about.  Only
        #: that identity may be overwritten, and only while the id field still
        #: names it.
        self._refused_identity = ""
        #: The folder the last successful snapshot save wrote, so the offer to
        #: open it points at what was actually written rather than at a guess.
        self._saved_snapshot_root: Path | None = None
        self._pattern_items: list[pg.PlotDataItem] = []
        self._pattern_key: tuple[int, int] | None = None
        self._selected_trace: int | None = None
        self._detector_key: tuple[str, int] | None = None
        #: Pooled order-plot decorations.  An order change moves and relabels
        #: these; it never destroys and rebuilds them, which is what made
        #: arrow-key order scrolling crawl on a real 2560x2160 frame.  There is
        #: one pool because there is one expected-line list (F17 item 2).
        self._line_pool: list[tuple[object, object]] = []
        self._anchor_scatter: pg.ScatterPlotItem | None = None
        #: Each histogram plot's own curve and its two threshold lines, built
        #: once and moved thereafter — the same pooling the order plot uses.
        self._histogram_items: dict[int, tuple] = {}
        self._catalog_cache: dict[tuple, tuple] = {}
        #: Which solved correction the cached expected-line lists were placed
        #: with.  A new solve makes it stale, and the next refresh notices.
        self._drawn_correction: tuple[float, float, float] | None = None
        self._catalog_rows: tuple = ()
        self._queue: list[Path] = []
        #: Whether files have ARRIVED and not yet been judged as a set.  A load
        #: is not the same event as an arrival: the bench opens frames of its
        #: own accord too, and only an arrival puts a folder's filenames up for
        #: the whole-drop judgement in ``_apply_unambiguous_suggestions``.
        self._arrivals_pending = False
        self._file_rows: list[Path] = []
        #: The last line the bench wrote about the roles themselves.  While it
        #: is still the line on screen, the bench's own follow-up narration
        #: holds its tongue rather than burying it (see ``_bench_says``).
        self._role_notice = ""
        #: Files the operator has deliberately put back to no role.  A folder
        #: whose names read unambiguously is applied unasked, and that must
        #: never undo a decision somebody made by hand.
        self._declined_suggestions: set[Path] = set()
        self._explainable_widgets: list[object] = []
        self._checklist_labels: list[object] = []
        #: Files the operator opened deliberately; the fit view never argues
        #: with an explicit choice, it only warns about what it is showing.
        self._explicitly_opened: set[Path] = set()
        self._landed_on: tuple[Path, str] | None = None
        self._auto_following = False
        self._family_override = False
        #: Splitters already listening for the operator's drag.
        self._splitter_keys: set[str] = set()
        #: Whether a real layout pass has happened yet.  Splitter cuts laid
        #: down before one are cuts of a size the window does not have.
        self._first_layout_done = False
        #: Re-entrancy guard.  Setting a splitter's sizes resizes its children,
        #: which resizes this window's own layout, which would ask for the cut
        #: again — a recursion Qt answers by overflowing the paint stack.
        self._distributing = False
        #: The Why dock's height as this window last asked for it, and whether
        #: the answer to that ask is still on its way.  A height that arrives
        #: outside an ask is the operator's own drag and is remembered as one.
        self._details_dock_height = 0
        self._sizing_details = False
        self.last_folder = (
            self.calibration_folder
            if self.calibration_folder is not None
            else (Path(watcher.folder) if watcher is not None else Path.cwd())
        )
        self._build_ui()
        self._connect_ui()
        self.setAcceptDrops(True)
        self.refresh()

        self.poll_timer = QtCore.QTimer(self)
        self.poll_timer.setInterval(int(poll_interval_ms))
        self.poll_timer.timeout.connect(self.poll_watch_folder)
        if start_timer and self.watcher is not None:
            self.poll_timer.start()
            QtCore.QTimer.singleShot(0, self.poll_watch_folder)

    def _build_ui(self) -> None:
        self.setWindowTitle("Echelle calibration bench")
        self.setWindowIcon(bench_window_icon())
        self.layout_unit = bench_layout_unit()
        self.setMinimumSize(*bench_minimum_size(unit=self.layout_unit))
        self.resize(*bench_default_geometry(unit=self.layout_unit))
        self.body_pt, self.headline_pt = bench_point_sizes()

        # The readings strip runs across the top of the whole window rather
        # than down the middle column.  The owner offered "top of center or top
        # of a rail"; with two rails flanking it the centre is the narrowest
        # column on screen, and a strip of readings put there folds onto two
        # and three rows and eats exactly the plot height the rails were meant
        # to leave it.  Full width, it is one row at any window size, and it is
        # in view whichever rail, tab or plot is being worked in.
        outer = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        outer.setChildrenCollapsible(False)
        outer.setHandleWidth(8)
        self.readings_splitter = outer
        self.setCentralWidget(outer)

        root = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        root.setChildrenCollapsible(False)
        root.setHandleWidth(8)
        self.root_splitter = root

        controls = QtWidgets.QWidget()
        # No maximum: the left rail is the operator's to widen, and its text is
        # the reading, not a caption.
        controls_layout = QtWidgets.QVBoxLayout(controls)
        controls_layout.setContentsMargins(14, 12, 12, 12)
        controls_layout.setSpacing(8)

        title = QtWidgets.QLabel("LIVE CALIBRATION")
        title.setObjectName("benchTitle")
        subtitle = QtWidgets.QLabel(
            "Drop files → triage → assign roles → fit lines → compare → validate"
        )
        subtitle.setWordWrap(True)
        subtitle.setObjectName("benchSubtitle")
        controls_layout.addWidget(title)
        controls_layout.addWidget(subtitle)

        # LEFT RAIL — nothing but controls.  The tables used to be stacked
        # under them, which is why every control tab had to scroll: two working
        # surfaces were competing for one column's height (F20).
        self.control_tabs = QtWidgets.QTabWidget()
        self._build_files_tab()
        self._build_procedure_tab()
        self._build_lamp_tab()
        self._build_save_tab()
        controls_layout.addWidget(self.control_tabs, 1)

        # The procedure is good but it is behind a tab, so it is only there
        # when you go looking (owner, 2026-08-16: "Do we have a spot for it, so
        # it will be shown?").  This is the spot, and deliberately the smallest
        # one that answers the question the checklist is consulted for: what is
        # the next thing to do, and what is stopping it.  It sits below the
        # tabs, so it is on screen whichever tab is open.  Minimal on purpose —
        # the owner's own warning was that agents would make a big mess of the
        # full summary, so this is a prototype for his verdict, not a design.
        # One line that shortens itself rather than wrapping: a strip below the
        # tabs takes its height from every one of them, and a wrapped "Next:"
        # running to three lines put the Save tab into a scrollbar. The whole
        # sentence is a click away in the Why dock, which is F16's own rule for
        # readings.
        # The step and the button that performs it, side by side.  They used to
        # be on opposite diagonals of the window — the next-step line bottom
        # left, its action top right — and the operator was expected to carry
        # one to the other (owner, 2026-08-16: "Not this scatter of actions
        # here and there!"). Whatever the checklist says is next, the button
        # beside it does that.
        next_panel = QtWidgets.QWidget()
        next_panel.setObjectName("nextStep")
        next_layout = QtWidgets.QHBoxLayout(next_panel)
        next_layout.setContentsMargins(7, 5, 7, 5)
        next_layout.setSpacing(8)
        self.next_step_value = _ElidingLabel("Drop a SIF onto the bench to begin.")
        next_layout.addWidget(self.next_step_value, 1)
        self.next_step_button = QtWidgets.QPushButton("Add SIF files…")
        self.next_step_button.setObjectName("nextStepButton")
        next_layout.addWidget(self.next_step_button)
        controls_layout.addWidget(next_panel)
        self.next_step_panel = next_panel

        self.controls_rail = controls
        root.addWidget(controls)

        # CENTER — the plots, and nothing else.  They are the flexible element:
        # they give width to the rails rather than the other way round, because
        # a plot can be zoomed and a table cannot (owner: "smaller plot view is
        # fine, I can zoom, but mostly I don't need anyway").
        self.view_tabs = QtWidgets.QTabWidget()
        self._build_triage_view()
        self._build_lamp_fit_view()
        self._build_sphere_view()
        root.addWidget(self.view_tabs)

        # RIGHT RAIL — the two working tables, one above the other, each as
        # long as the window allows.
        self._build_tables_rail()
        root.addWidget(self.tables_rail)

        # Stretch factors, not pixel cuts: every one of these is a share of
        # whatever width the window turns out to have, so maximizing widens all
        # three rather than handing the whole gain to the plots.
        for index, stretch in enumerate(_ROOT_STRETCH):
            root.setStretchFactor(index, stretch)
        root.splitterMoved.connect(lambda *_args: self._relayout_wrapped_text())

        self._build_status_band()
        outer.addWidget(self.status_band)
        outer.addWidget(root)
        outer.setStretchFactor(0, 0)
        outer.setStretchFactor(1, 1)

        self._build_details_dock()

        body = self.body_pt
        headline = self.headline_pt
        # Sizes are set through real fonts rather than through the stylesheet:
        # a stylesheet font-size is resolved through the platform's pixel
        # metrics and quietly comes back smaller on some displays, which is
        # exactly the squinting this is meant to end.
        base_font = self.font()
        base_font.setPointSizeF(body)
        self.setFont(base_font)
        # A group box's title sits in its top margin.  Sizing that margin from
        # the real font metrics is what stops "Bench state" from garbling into
        # its own border on a display whose font is larger than the designer's.
        title_height = QtGui.QFontMetrics(base_font).height()
        self.setStyleSheet(
            f"""
            QMainWindow, QWidget {{ background: #151b22; color: #dce8f2; }}
            QGroupBox {{ border: 1px solid #334252; border-radius: 7px;
                        margin-top: {title_height + 6}px;
                        padding-top: {max(6, title_height // 2)}px; font-weight: 600; }}
            QGroupBox::title {{ subcontrol-origin: margin; subcontrol-position: top left;
                        left: 10px; padding: 0 4px; color: #8fd9ff; }}
            #benchTitle {{ color: #80ddff; font-weight: 700; letter-spacing: 2px; }}
            #benchSubtitle, #benchHelp, #mutedText {{ color: #93a8b8; }}
            #triageHeadline {{ font-weight: 700; padding: 6px; }}
            #stateBadge {{ color: #7ee2b8; font-weight: 700; }}
            #warningPanel {{ background: #2a1e13; border-left: 3px solid #ffb86b;
                            padding: 7px; color: #ffd6a3; }}
            #messagePanel {{ background: #0f141a; border-left: 3px solid #49b5df;
                            padding: 7px; color: #bed4e1; }}
            #nextStep {{ background: #14231d; border-left: 3px solid #70d6ae;
                        padding: 7px; color: #cdeadd; }}
            QTableWidget, QTreeWidget, QListWidget, QPlainTextEdit, QTextBrowser {{
                background: #10151b; alternate-background-color: #18212a;
                gridline-color: #2b3946; }}
            QHeaderView::section {{ background: #202b36; color: #b9d5e5; padding: 4px; }}
            QPushButton {{ background: #273746; border: 1px solid #416078; border-radius: 5px;
                          padding: 6px 8px; }}
            QPushButton:hover {{ background: #315069; }}
            QPushButton:disabled {{ color: #657786; border-color: #33404a; }}
            QSpinBox, QComboBox, QLineEdit {{
                background: #0f141a; border: 1px solid #416078; padding: 3px 4px; }}
            QSplitter::handle {{ background: #22303c; }}
            QSplitter::handle:hover {{ background: #3a5f79; }}
            QTabWidget::pane {{ border: 1px solid #334252; }}
            QTabBar::tab {{ background: #202b36; color: #9fb6c6; padding: 6px 10px; }}
            QTabBar::tab:selected {{ background: #294052; color: #8fe3ff; }}
            QToolTip {{ background: #0f141a; color: #dce8f2;
                       border: 1px solid #49b5df; padding: 5px; }}
            """
        )
        # The loudness budget: exactly one element carries the headline size,
        # and it is the verdict the next action is read from.  Everything else
        # — including this window's own title — is body text.
        _emphasise(title, body)
        for widget in (
            self.rms_value,
            self.anchor_count_value,
            self.comparison_value,
            self.file_state_value,
            self.alignment_state_value,
            self.save_state_value,
        ):
            widget.setObjectName("reading")
            _emphasise(widget, body)
        _emphasise(self.triage_headline, headline)
        _emphasise(self.exposure_value, body, bold=False)
        controls.setMinimumWidth(self._controls_minimum_width())
        self.tables_rail.setMinimumWidth(self._tables_minimum_width())
        self.expected_lines_panel.setMinimumHeight(
            self._table_panel_minimum_height(
                self.expected_lines_panel, self.line_help_table, self.line_panel_header
            )
        )
        self.anchors_panel.setMinimumHeight(
            self._table_panel_minimum_height(
                self.anchors_panel,
                self.anchor_table,
                self.anchor_buttons,
                # Only the strip's floor is a demand on the panel: above that it
                # spends what the table has not claimed, so it can never be the
                # reason this panel refuses to shrink.
                extra=BENCH_RESIDUAL_FLOOR_LINES * self.layout_unit,
            )
        )
        self._measure_status_band()
        self._distribute_space()

    #: What "workable" means for either working table, in rows rather than in a
    #: pixel guess that a larger platform font would quietly invalidate.
    #: Below this either table is the corner box the owner was handed, not a
    #: working surface — so it is what the acceptance tests read the rail
    #: against at every window size the bench opens at.
    EXPECTED_LINE_ROWS = 6

    #: And what the rail refuses to shrink either table below, whatever else
    #: has to give.  Lower than the workable count on purpose: a display too
    #: small for the whole layout should lose rows from both tables evenly
    #: rather than crush one of them out of existence to spare the other.
    TABLE_FLOOR_ROWS = 3

    def _table_panel_minimum_height(
        self, panel, table, header_widget, *, extra: int = 0
    ) -> int:
        """The height below which a rail panel stops being a working surface.

        Derived from the table's own metrics, so the floor moves with the
        platform font instead of pinning a row count to one designer's display.
        *extra* is anything else the panel carries beside the table and its
        header — a capped strip, say — which the floor has to include or the
        table pays for it out of its own rows.
        """

        layout = panel.layout()
        margins = layout.contentsMargins()
        return int(
            (extra + layout.spacing() if extra else 0)
            + self.TABLE_FLOOR_ROWS * table.verticalHeader().defaultSectionSize()
            + table.horizontalHeader().sizeHint().height()
            + 2 * table.frameWidth()
            + header_widget.sizeHint().height()
            + margins.top()
            + margins.bottom()
            + layout.spacing()
            # The group box's title sits in its own top margin.
            + QtGui.QFontMetrics(self.font()).height()
            + 12
        )

    def _distribute_space(self) -> None:
        """Cut every splitter as a share of the size the window actually has.

        The old cuts were pixel numbers tuned at one geometry, laid down before
        the first layout pass — which is how a maximized bench still opened
        with a 240 px left rail the owner had to drag wide before he could use
        it.  Shares survive any window size; the operator's own drag is
        remembered as a share too, so it survives a resize as a *proportion*
        rather than as the pixel count it happened to be when the handle
        stopped moving.  The readings strip is the one exception: it is a strip
        and takes exactly what its content costs.
        """

        if self._distributing:
            return
        self._distributing = True
        try:
            self._apply_splitter_shares(self.root_splitter, "root", _ROOT_SHARES)
            self._apply_splitter_shares(
                self.tables_splitter, "tables", self._tables_default_shares()
            )
            self._watch_splitter(self.readings_splitter, "readings")
            self._pin_status_band_height()
            self._fit_residual_strip()
            if getattr(self, "details_dock", None) is not None:
                # The measure is a property of the font and the dock's width;
                # both are only true after a layout pass, so it is re-read on
                # the same passes every other cut in this window is made on.
                self._apply_reading_measure()
                self._fit_details_dock()
        finally:
            self._distributing = False

    def _workable_panel_heights(self):
        """What each rail panel costs with a workable table, measured.

        Measured rather than estimated, for the same reason the strip's own
        height is: the furniture — a group title, a wrapping header, two
        buttons — is whatever this platform's font makes of it, and a formula
        that guesses it is wrong by a row exactly when the rail is tightest.
        ``None`` while there is no geometry to read yet.
        """

        measured = []
        for panel, table, strip in (
            (self.anchors_panel, self.anchor_table, self.residual_widget),
            (self.expected_lines_panel, self.line_help_table, None),
        ):
            if panel.height() <= 0 or table.height() <= 0:
                return None
            furniture = panel.height() - table.height()
            floor = 0
            if strip is not None:
                furniture -= strip.height()
                floor = BENCH_RESIDUAL_FLOOR_LINES * self.layout_unit
            measured.append(
                furniture
                + floor
                + self.EXPECTED_LINE_ROWS * table.verticalHeader().defaultSectionSize()
                + (table.height() - table.viewport().height())
            )
        return tuple(measured)

    def _tables_default_shares(self):
        """Cut the rail where each panel becomes workable, not at a fixed guess.

        The fixed pair was tuned when the anchors panel held nothing but a
        table.  It now holds a labelled residual strip as well, and a share
        that cannot know that spends the difference out of the anchor table's
        own rows — which is the one thing the rail's law forbids.  So the
        default cut is read off what the two panels actually cost with workable
        tables, and only what neither of them needs is split evenly.  A rail
        too short for both falls back to the fixed shares, which lose rows from
        the two tables evenly rather than starving one of them.
        """

        needs = self._workable_panel_heights()
        extent = self.tables_splitter.height() or self.height()
        if needs is None or extent <= 0 or sum(needs) > extent:
            return _TABLES_SHARES
        spare = (extent - sum(needs)) / 2.0
        return tuple((need + spare) / extent for need in needs)

    def _fit_residual_strip(self) -> None:
        """Spend on the residual strip only what the anchor table has not.

        Layout, not data: this runs on the same show, resize and splitter-drag
        passes every other cut in the rail is made on, and never on an anchor
        being added or a transform being solved.  The table keeps its workable
        rows first and the strip takes the remainder, so the panel reads the
        same way on the geometry the bench opens at and on a maximized one —
        the strip is a sparkline on the first and a real strip on the second,
        rather than a fixed block that eats four of the table's six rows.
        """

        table = self.anchor_table
        if self.anchors_panel.height() <= 0 or table.height() <= 0:
            # Before the first layout pass there is no geometry to divide; the
            # show and resize passes both come back through here.
            return
        # Measured, not estimated: what the table and the strip share is
        # whatever the panel's furniture — the group title, the buttons, the
        # margins — is not already using, and the table's claim on it is its
        # workable rows plus its own header and frame.
        shared = table.height() + self.residual_widget.height()
        workable = self.EXPECTED_LINE_ROWS * table.verticalHeader().defaultSectionSize() + (
            table.height() - table.viewport().height()
        )
        spare = shared - workable
        height = int(
            min(
                BENCH_RESIDUAL_LINES * self.layout_unit,
                max(BENCH_RESIDUAL_FLOOR_LINES * self.layout_unit, spare),
            )
        )
        if height != self.residual_widget.maximumHeight():
            self.residual_widget.setMaximumHeight(height)
        # How many wavelengths fit is a question about the width this pass just
        # settled, so it is answered in the same pass.
        self._label_residual_axis()

    def _apply_splitter_shares(self, splitter, key: str, default) -> None:
        """Lay this splitter's cut down as shares of its current extent."""

        shares = _SESSION_SPLITTER_SHARES.get(key, tuple(default))
        extent = (
            splitter.width()
            if splitter.orientation() == QtCore.Qt.Horizontal
            else splitter.height()
        )
        if extent <= 0:
            # Before the first layout pass a splitter has no extent to divide.
            # Hand it the shares against the window instead, and the show/resize
            # passes will re-cut it against the real one.
            extent = (
                self.width()
                if splitter.orientation() == QtCore.Qt.Horizontal
                else self.height()
            )
        total = sum(shares) or 1.0
        splitter.setSizes([max(1, int(extent * share / total)) for share in shares])
        self._watch_splitter(splitter, key)

    def _watch_splitter(self, splitter, key: str) -> None:
        """Remember this splitter's cut from the moment the operator drags it."""

        if key in self._splitter_keys:
            return
        self._splitter_keys.add(key)
        splitter.splitterMoved.connect(
            lambda *_args, _s=splitter, _k=key: self._cut_moved(_s, _k)
        )

    def _cut_moved(self, splitter, key: str) -> None:
        """A dragged handle is remembered, and re-spends the rail's leftovers."""

        self._remember_cut(splitter, key)
        self._fit_residual_strip()

    @staticmethod
    def _remember_cut(splitter, key: str) -> None:
        """Record a dragged handle as shares, which any window size can restore."""

        sizes = splitter.sizes()
        total = sum(sizes)
        if total <= 0:
            return
        _SESSION_SPLITTER_SHARES[key] = tuple(size / total for size in sizes)

    def _controls_minimum_width(self) -> int:
        """Widen the left rail until its own controls stop clipping.

        The minimum is read off the content rather than guessed, which is what
        makes the default geometry legible instead of merely large.
        """

        widest = self.control_tabs.tabBar().sizeHint().width()
        for widget in self.control_tabs.findChildren(
            (QtWidgets.QPushButton, QtWidgets.QComboBox)
        ):
            widest = max(widest, widget.sizeHint().width())
        return int(min(560, max(400, widest + 60)))

    def _tables_minimum_width(self) -> int:
        """The width at which the rail still shows every fixed column it has.

        Read off the anchor table's own columns rather than guessed, so the
        rail's floor moves with the platform font instead of being a number
        measured once on one display at 100% scaling.
        """

        columns = sum(
            self.anchor_table.columnWidth(index)
            for index in range(self.anchor_table.columnCount())
        )
        # Only the anchor table is a demand: its columns are fixed, so a
        # narrower rail would simply hide one.  The expected-line table's own
        # columns resize to their contents around a stretching first column and
        # shrink honestly with the rail.
        return int(columns + 3 * self.layout_unit)

    def _build_tables_rail(self) -> None:
        """The right rail: the anchor table and the expected-line table.

        The owner's law, verbatim: TWO RAILS.  Controls on the left, tables on
        the right, plots in between and deliberately smaller — "smaller plot
        view is fine, I can zoom, but mostly I don't need anyway".  What the
        rail is really for is that these two tables stop competing with the
        controls for one column's height: stacked under the tabs they turned
        every control tab into a scrolling one and left both tables four rows
        tall.  Down their own rail each of them is as long as the window.
        """

        rail = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        rail.setChildrenCollapsible(False)
        rail.setHandleWidth(8)
        self.tables_splitter = rail
        self.tables_rail = rail
        self._build_anchors_panel()
        self._build_expected_lines_panel()
        for index, stretch in enumerate(_TABLES_STRETCH):
            rail.setStretchFactor(index, stretch)

    def _build_anchors_panel(self) -> None:
        """The anchor table, its residual strip, and its two buttons.

        Down the right rail beside the lines it anchors.  The scalar numbers
        the fit produces — RMS, the transform — are readings and live on the
        readings strip, where they are in view from every tab rather than
        behind whichever one happens to be open.  The per-anchor residuals are
        not a scalar: they are the shape of the disagreement, one bar per row
        of this very table, so they belong *here*, under the rows they belong
        to (owner, 2026-08-17: "what about moving the anchor plot into the
        anchors tab?").  Under the spectrum they were a plot the operator had
        to leave the anchor table to read.
        """

        panel = QtWidgets.QGroupBox("Anchors")
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 4)
        # Three stacked things now instead of two, so the gaps between them are
        # tighter than the rail's other panel: every pixel of white space here
        # is one the table or the wavelengths under it did not get, and this
        # panel is the one that has to fit both.
        layout.setSpacing(3)

        self.anchor_table = QtWidgets.QTableWidget(0, 5)
        self.anchor_table.setHorizontalHeaderLabels(
            ["Ord", "λ nm", "Δx px", "Resid", "QC"]
        )
        self.anchor_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Fixed)
        self.anchor_table.horizontalHeader().setFixedHeight(28)
        for column, width in enumerate((44, 78, 68, 68, 48)):
            self.anchor_table.setColumnWidth(column, width)
        self.anchor_table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.anchor_table.verticalHeader().setVisible(False)
        self.anchor_table.verticalHeader().setDefaultSectionSize(24)
        self.anchor_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.anchor_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        layout.addWidget(self.anchor_table, 1)

        # The residuals, directly under the rows they score.  A strip, not a
        # plot: the same subordination as the triage top end (BENCH_TOP_END_
        # LINES), because the table is this panel's working surface and the
        # residuals qualify it.  Capped height, no appetite for growth, and no
        # minimum width at all — the rail's width is the rail's, and this
        # adapts to it rather than widening it (the two-rails law).
        self.residual_widget = pg.GraphicsLayoutWidget()
        self.residual_widget.setBackground("#10151b")
        # Every pixel of this widget's height is a pixel the anchor table did
        # not get, so the graphics layout's own padding is given back: at the
        # height a short panel can spare, pyqtgraph's default margins were most
        # of the strip and the bars were a hairline.
        self.residual_widget.ci.layout.setContentsMargins(0, 0, 0, 0)
        self.residual_widget.ci.layout.setSpacing(0)
        # No plot title either: at that height a title would be the whole of
        # it.  What the strip is, is said by the axis it is drawn against, by
        # the Resid column right above it, and in full by the Why dock when the
        # strip is clicked.
        self.residual_plot = self.residual_widget.addPlot(row=0, col=0)
        self.residual_plot.setContentsMargins(0, 2, 2, 0)
        self.residual_plot.setLabel("left", "resid", units="px")
        self.residual_plot.getAxis("left").enableAutoSIPrefix(False)
        self.residual_plot.getAxis("left").setWidth(38)
        # The wavelengths stay under the bars.  They overlapped even across the
        # full width of the old centre view, and a rail is narrower still — so
        # what adapts is the labelling, not the labels: a smaller font than the
        # bench uses anywhere else, and as many wavelengths as the width can
        # hold, spread evenly across the anchors when it cannot hold them all.
        # Every anchor is still named somewhere: clicking a bar selects its row,
        # which carries the order, the wavelength and the residual in figures.
        self._residual_tick_font = QtGui.QFont(self.font())
        self._residual_tick_font.setPointSizeF(BENCH_RESIDUAL_TICK_POINT_SIZE)
        self.residual_plot.getAxis("bottom").setStyle(
            tickFont=self._residual_tick_font,
            tickLength=3,
            tickTextOffset=1,
        )
        self.residual_plot.addLine(
            y=0, pen=pg.mkPen("#64748b", style=QtCore.Qt.DashLine)
        )
        # A strip is read, not navigated: panning it would only lose the bars,
        # and a stray drag must never be mistaken for the click that selects.
        self.residual_plot.setMouseEnabled(x=False, y=False)
        self.residual_plot.setMenuEnabled(False)
        self.residual_widget.setMaximumHeight(
            BENCH_RESIDUAL_FLOOR_LINES * self.layout_unit
        )
        self.residual_widget.setMinimumWidth(0)
        self.residual_widget.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum
        )
        self._explainable(
            self.residual_widget,
            "How far each anchor sits from the solved transform",
            "One bar per accepted anchor, in the same order as the rows above, "
            "showing the along-dispersion residual in detector pixels. The bars "
            "carry no wavelength labels because the rail is too narrow for one "
            "per anchor: click a bar and the table above selects that anchor, "
            "names its line, and explains it here. Two anchors in different "
            "orders are needed before there is a transform to residual against.",
        )
        layout.addWidget(self.residual_widget, 0)

        #: The bars currently drawn, kept so a table selection can recolour one
        #: of them without rebuilding the item.
        self._residual_bars = None

        # The action that FILLS this table stands on it, always visible, in
        # every tab.  Putting it away in the Lamp fit control tab was the
        # owner's own complaint back again — a bench that does not teach where
        # the next action lives — and an empty table with two greyed buttons
        # says nothing about how to stop being empty.
        buttons = QtWidgets.QWidget()
        button_column = QtWidgets.QVBoxLayout(buttons)
        button_column.setContentsMargins(0, 0, 0, 0)
        button_column.setSpacing(4)
        self.auto_anchor_button = QtWidgets.QPushButton("Auto-anchor lines")
        self._explainable(
            self.auto_anchor_button,
            "Anchoring the lines the calibration already trusts",
            "The curated wavelength table marks the rows vetted during the BH "
            "paper's own calibration — tried and tested for the Balmer and "
            "Fulcher analysis, which is the pedigree that makes them "
            "trustworthy, not their brightness. This measures every one of "
            "them, in every order, with the same centroid fit and the same "
            "raw-detector saturation guard a click uses, and anchors the ones "
            "that pass. It is what echelle-align has always done headlessly. "
            "Anything it declines is listed with its reason, anchors you "
            "placed by hand are kept, and any anchor it places comes off with "
            "a right-click on the spectrum.",
        )
        button_column.addWidget(self.auto_anchor_button)
        # Remove and Clear pair on one row so the primary action costs this
        # panel no height at all: they are short, secondary, and act on a
        # selected row rather than on the frame.  Each carries a minimum width
        # read off its own text, so neither can elide however narrow the rail.
        secondary = QtWidgets.QHBoxLayout()
        secondary.setContentsMargins(0, 0, 0, 0)
        secondary.setSpacing(4)
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.clear_button = QtWidgets.QPushButton("Clear")
        for button in (self.remove_button, self.clear_button):
            button.setMinimumWidth(button.sizeHint().width())
            secondary.addWidget(button)
        button_column.addLayout(secondary)
        layout.addWidget(buttons)

        self.anchor_buttons = buttons
        self.anchors_panel = panel
        self.tables_splitter.addWidget(panel)

    def loud_widgets(self) -> list[object]:
        """Every widget currently drawn at the headline size.

        The budget is one.  This is how the surface proves it, rather than a
        promise in a docstring nobody can check.
        """

        threshold = self.headline_pt - 0.01
        return [
            widget
            for widget in self.findChildren(QtWidgets.QWidget)
            if widget.font().pointSizeF() >= threshold
        ]

    def _build_details_dock(self) -> None:
        """One dock that explains whatever the operator just clicked.

        The whys are learned once and then remembered, so they live behind a
        click instead of standing permanently between the operator and the
        numbers — and behind a click only, so moving the mouse across the
        window never rewrites what is being read.
        """

        dock = QtWidgets.QDockWidget("Why this reading", self)
        dock.setObjectName("benchDetailsDock")
        dock.setAllowedAreas(QtCore.Qt.BottomDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        self.details_view = QtWidgets.QTextBrowser()
        self.details_view.setOpenExternalLinks(False)
        # A dock that explains things on request had claimed 213 px of a 920 px
        # window whether or not anything had been asked.  It reads calmly in
        # a few lines and drags taller the moment it needs to — and its height
        # is quoted in lines of the platform's own text, like every other size
        # in this window.
        self.details_view.setMinimumHeight(BENCH_DETAILS_FLOOR_LINES * self.layout_unit)
        self.details_view.setPlaceholderText(
            "Click any verdict, checklist row, file, or anchor and the whole "
            "explanation is written here. It changes only when asked."
        )
        # An absolute path is one long unbreakable word.  Left to itself it
        # would push the whole block wider than the measure and hand the dock a
        # horizontal scrollbar; broken anywhere, it wraps like everything else
        # and stays selectable text the operator can copy out.
        option = self.details_view.document().defaultTextOption()
        option.setWrapMode(QtGui.QTextOption.WrapAtWordBoundaryOrAnywhere)
        self.details_view.document().setDefaultTextOption(option)
        dock.setWidget(self.details_view)
        self.details_dock = dock
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        # The operator's own drag is remembered from here, by the same law the
        # splitter cuts follow: the filter reports the dock's height, and any
        # height this window did not itself ask for is the operator's.
        dock.installEventFilter(self)
        self._apply_reading_measure()
        self._fit_details_dock()
        self.explain(
            "Exposure triage is the front door",
            "Drop any SIF and the bench judges the exposure before any role "
            "exists: clustered full-scale pixels are real saturation, isolated "
            "ones are cosmic rays, and the histogram's top end is the number "
            "you adjust the lamp by.",
        )

    def _reading_measure(self) -> int:
        """How wide prose may get in the Why dock, in pixels of this font.

        A line stops being readable long before it stops fitting.  The dock is
        as wide as the window, so an explanation was setting itself across two
        thousand pixels and the eye lost its place on the way back to the left
        margin (owner, 2026-08-18: "all long and breaky").  The cap is stated
        in characters — typography's own measure — and turned into pixels by
        the platform's own font metrics, like every other size in this window.
        """

        metrics = QtGui.QFontMetrics(self.details_view.font())
        return int(BENCH_READING_MEASURE_CHARS * max(1, metrics.averageCharWidth()))

    def _apply_reading_measure(self) -> int:
        """Set the dock's wrap to the measure, or to the dock when it is narrower.

        Narrower wins so a dock dragged to the side rail wraps at the width it
        has instead of growing a horizontal scrollbar; the measure is a cap on
        prose, never a demand for space.
        """

        width = self._reading_measure()
        viewport = self.details_view.viewport().width()
        if viewport > 0:
            width = min(width, viewport)
        self.details_view.setLineWrapMode(QtWidgets.QTextEdit.FixedPixelWidth)
        self.details_view.setLineWrapColumnOrWidth(int(width))
        return int(width)

    def explain(self, title: str, text: str) -> None:
        """Show one full explanation in the details dock.

        Blank-line-separated paragraphs are set as paragraphs, with real space
        between them.  Several explanations now carry a path and a SELF-CHECK
        passage after the prose, and joining those with a bare line break ran
        three separate statements together into one grey wall.
        """

        gap = max(6, self.layout_unit // 2)
        paragraphs = [
            block.strip()
            for block in re.split(r"\n[ \t]*\n", str(text).strip())
            if block.strip()
        ]
        body = "".join(
            "<p style='margin-top:0px;margin-bottom:{gap}px'>{para}</p>".format(
                gap=gap, para=paragraph.replace("\n", "<br>")
            )
            for paragraph in paragraphs
        )
        self.details_view.setHtml(
            f"<h3 style='color:#8fd9ff;margin:0 0 6px 0'>{title}</h3>"
            f"<div style='color:#cfe0ec;line-height:145%'>{body}</div>"
        )

    def _spare_height(self) -> float | None:
        """Height no rail is using, measured — ``None`` before there is geometry.

        The dock spans the bottom of the whole window, so every pixel it takes
        comes out of both rails and the plots at once.  Both rails have a
        contract it must not break: six workable rows in each table, and no
        control column that scrolls.  So the surplus is the smaller of the two
        surpluses, and it is measured rather than assumed.
        """

        needs = self._workable_panel_heights()
        extent = self.tables_splitter.height()
        if needs is None or extent <= 0:
            return None
        spare = float(extent - sum(needs))
        for index in range(self.control_tabs.count()):
            tab = self.control_tabs.widget(index)
            inner = tab.widget() if isinstance(tab, QtWidgets.QScrollArea) else None
            if inner is None or tab.viewport().height() <= 0:
                continue
            spare = min(spare, float(tab.viewport().height() - inner.sizeHint().height()))
        return spare

    def _details_dock_ceiling_lines(self) -> float:
        """How tall the dock may open before it starts costing a rail its rows.

        At the bench's default geometry the vertical budget is already spent —
        the two tables are at exactly their six workable rows — so a dock that
        simply demanded the height an explanation costs would buy its comfort
        with the anchor table's rows.  It asks for that height and settles for
        whatever the rails genuinely are not using, keeping a line of slack in
        hand so a font or a scaling that measures slightly differently cannot
        turn the last spare pixels into a scrolling control column.
        """

        current = self.details_dock.height()
        spare = self._spare_height()
        if current <= 0 or self.layout_unit <= 0 or spare is None:
            return float(BENCH_DETAILS_LINES)
        # Signed on purpose: a dock that has already taken more than the rails
        # could spare has to give it back, not merely stop growing.
        return max(
            BENCH_DETAILS_FLOOR_LINES,
            (current + spare - self.layout_unit) / self.layout_unit,
        )

    def _details_dock_lines(self) -> float:
        """How many lines tall the Why dock should be: the drag, or the default.

        A drag is an instruction and is obeyed as given; the default is a wish
        and yields to the rails.
        """

        dragged = _SESSION_DOCK_LINES.get("details")
        if dragged is not None:
            return max(BENCH_DETAILS_FLOOR_LINES, dragged)
        return max(
            BENCH_DETAILS_FLOOR_LINES,
            min(float(BENCH_DETAILS_LINES), self._details_dock_ceiling_lines()),
        )

    def _fit_details_dock(self) -> None:
        """Open the dock at a height that holds a whole ordinary explanation."""

        height = int(self._details_dock_lines() * self.layout_unit)
        self._details_dock_height = height
        if height == self.details_dock.height():
            return
        # Whatever comes back from this is this window's own doing, however far
        # it lands from what was asked, and must never be filed as a drag.
        self._sizing_details = True
        self.resizeDocks([self.details_dock], [height], QtCore.Qt.Vertical)

    def _details_dock_resized(self, height: int) -> None:
        """Remember a height this window did not ask for: the operator's drag."""

        if height <= 0 or self.layout_unit <= 0:
            return
        if self._sizing_details:
            self._sizing_details = False
            # What the layout could actually give, which is now the baseline.
            self._details_dock_height = height
            return
        if abs(height - self._details_dock_height) <= 2:
            return
        _SESSION_DOCK_LINES["details"] = height / self.layout_unit
        self._details_dock_height = height

    def _explainable(self, widget, title: str, text: str, *, hint: str = "") -> None:
        """Route one widget's explanation to the Why dock, not to a tooltip.

        The tooltip keeps a single short line — enough to say that there is
        more — and clicking or focusing the widget writes the whole text into
        the dock, where it can be read without covering the controls.  Hover
        deliberately does nothing: the dock holds still while it is read.
        """

        widget.setToolTip(hint or one_line(text))
        if isinstance(widget, QtWidgets.QLabel):
            # A button already looks clickable; a reading does not, so only the
            # readings advertise that there is more behind them.
            widget.setCursor(QtCore.Qt.WhatsThisCursor)
        widget.setProperty("explainTitle", title)
        widget.setProperty("explainText", text)
        if widget not in self._explainable_widgets:
            self._explainable_widgets.append(widget)
            widget.installEventFilter(self)

    def _forget_explainable(self, widget) -> None:
        """Drop a rebuilt widget's routing so the filter list cannot grow."""

        if widget in self._explainable_widgets:
            self._explainable_widgets.remove(widget)
            widget.removeEventFilter(self)

    #: Asking is a click or a focus, never a hover (F17 item 5).  A dock that
    #: rewrote itself under the pointer yanked the sentence away mid-read
    #: whenever the mouse crossed something on its way somewhere else.
    _EXPLAIN_EVENTS = (
        QtCore.QEvent.MouseButtonPress,
        QtCore.QEvent.FocusIn,
    )

    def eventFilter(self, source, event) -> bool:  # noqa: N802 - Qt naming
        """Send any explained widget's click or focus to the dock."""

        if event.type() == QtCore.QEvent.Resize and source is getattr(
            self, "details_dock", None
        ):
            self._details_dock_resized(event.size().height())
        if event.type() in self._EXPLAIN_EVENTS:
            title = source.property("explainTitle")
            if title:
                self.explain(title, source.property("explainText") or "")
        return super().eventFilter(source, event)

    def showEvent(self, event) -> None:  # noqa: N802 - Qt naming
        """Cut the rails against the size the window really got.

        A splitter sized inside ``_build_ui`` is sized against a window Qt has
        not laid out yet, so the numbers belong to no geometry at all — that is
        half of "opens cramped, drag it before use".  The first show re-cuts
        them, and the queued pass re-cuts them once more after the layout the
        show itself triggers has settled.
        """

        super().showEvent(event)
        self._distribute_space()
        if not self._first_layout_done:
            self._first_layout_done = True
            QtCore.QTimer.singleShot(0, self._distribute_space)

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt naming
        super().resizeEvent(event)
        # Shares, re-applied: a maximized window widens all three columns
        # instead of spending the whole gain on the plots and leaving the rails
        # at the width they had in a small window.
        self._distribute_space()
        self._relayout_wrapped_text()

    def _relayout_wrapped_text(self) -> None:
        """Give every wrapping label the height its own text needs.

        A word-wrapped QLabel in a vertical layout is handed its one-line
        height and silently swallows the rest, which is how the left pane came
        to overlap itself.  Asking each label what its text costs at its real
        width, and reserving that, is what makes the default geometry legible.
        """

        if getattr(self, "checklist_tree", None) is None:
            return
        if getattr(self, "status_band", None) is not None:
            # Ask the column, not the strip: the splitter has its width before
            # the strip inside it does, so reading the strip on the first pass
            # folds the band against a width it is about to grow out of.
            self._reflow_status_band(
                self._status_band_columns(self.readings_splitter.width())
            )
        width = self._checklist_row_width()
        for row in range(self.checklist_tree.count()):
            item = self.checklist_tree.item(row)
            label = self.checklist_tree.itemWidget(item)
            if label is None:
                continue
            self._fit_checklist_row(item, label, width)
        for label in self.findChildren(QtWidgets.QLabel):
            if self._sizes_itself(label):
                continue
            needed = label.heightForWidth(label.width())
            if needed > 0 and label.minimumHeight() != needed:
                label.setMinimumHeight(needed)
        self.file_table.resizeRowsToContents()
        if getattr(self, "status_band", None) is not None:
            # The wrapping readings have just been given their real heights, so
            # this is the first moment the strip's true cost is known.
            self._pin_status_band_height()

    def _sizes_itself(self, label) -> bool:
        """Whether *label* is one this pass must keep its hands off.

        Two surfaces size their own labels: the checklist gives each row the
        height of its own item, and the readings strip works through its size
        policy.  Writing a minimum onto either would feed the next pass its
        own output, and the layout would walk itself taller every resize.
        """

        if not label.wordWrap() or label.width() <= 0:
            return True
        if label.parent() is self.checklist_tree.viewport():
            return True
        return self.status_band.isAncestorOf(label)

    def _build_procedure_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        intro = QtWidgets.QLabel(
            "This list is built from the files and roles you actually assigned. "
            "A row that is not possible yet says what would unblock it, and lamp "
            "advice never blocks anything."
        )
        intro.setWordWrap(True)
        intro.setObjectName("messagePanel")
        layout.addWidget(intro)
        checklist_header = QtWidgets.QLabel("Procedure from your data · measured evidence")
        checklist_header.setStyleSheet("background: #202b36; padding: 8px;")
        layout.addWidget(checklist_header)
        self.checklist_tree = QtWidgets.QListWidget()
        self.checklist_tree.setAlternatingRowColors(True)
        self.checklist_tree.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        # A row is one reading across the whole width.  Selecting it is
        # meaningful — it writes the row's explanation into the Why dock — so
        # selection stays on, but it must never paint a highlight over bare
        # background beside the text (F18 item 4).
        self.checklist_tree.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.checklist_tree.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self._fills_its_share(self.checklist_tree)
        layout.addWidget(self.checklist_tree, 1)
        self.control_tabs.addTab(self._scrollable(tab), "Procedure")

    def _checklist_row_width(self) -> int:
        """The width one checklist row is allowed — the whole viewport.

        The old code pinned every row label to ``max(280, viewport - 18)``,
        computed once while the tab had never been laid out.  The label then
        stayed 280 px wide inside a pane twice that, and the strip of untouched
        list background to its right — complete with the selected row's
        highlight — is the "phantom empty second column" in the owner's
        screenshot.  There is one column, and it is as wide as the viewport.
        """

        tree = self.checklist_tree
        # The viewport already excludes a scrollbar when one is shown, so there
        # is nothing further to reserve — the old ``- 18`` was reserving it a
        # second time and leaving that strip bare.
        return max(1, tree.viewport().width() - 2 * tree.frameWidth())

    def _fit_checklist_row(self, item, label, width: int) -> None:
        """Give one row the full width and exactly the height its text needs.

        ``heightForWidth`` is the answer to the only question worth asking —
        how tall is this text at the width it will actually be drawn at — so
        it is used, alone.  Taking ``max()`` of it and ``sizeHint().height()``
        was what padded every row: a word-wrapped label's size hint is a
        heuristic computed against a width the label does not have, it came
        out taller than the truth on nine rows out of ten, and it won the
        comparison every time.  Rows stood 48 to 114 px taller than their own
        text, which is the "LOT of vertical padding" in the owner's screenshot
        and the reason the list scrolled with five items in it.
        """

        label.setMinimumWidth(width)
        label.setMaximumWidth(width)
        needed = label.heightForWidth(width)
        if needed <= 0:
            needed = label.sizeHint().height()
        item.setSizeHint(QtCore.QSize(width, needed + 4))

    @staticmethod
    def _checklist_row_html(symbol: str, headline: str, unblocked_by: str) -> str:
        """One row as a hanging-indented block.

        A two-cell table is what makes a wrapped line start under the first
        character of the text rather than under the status glyph, and the
        sub-line is indented by a real margin instead of by leading spaces that
        a proportional font renders at some arbitrary width.
        """

        body = html.escape(headline)
        if unblocked_by:
            body += (
                '<div style="margin-top:3px; margin-left:14px;">'
                f"unblocked by: {html.escape(unblocked_by)}</div>"
            )
        return (
            '<table cellspacing="0" cellpadding="0" style="margin:0;">'
            '<tr><td style="padding-right:8px;" valign="top">'
            f"{html.escape(symbol)}</td>"
            f'<td valign="top">{body}</td></tr></table>'
        )

    def _fills_its_share(self, widget) -> None:
        """Let this widget take the room it is given and shrink when it isn't.

        A control page lives inside a QScrollArea, and a resizable scroll area
        sizes its page by the layout's ``heightForWidth`` — which counts every
        child's *preferred* height, not its minimum.  One list or preview with
        a 192 px preferred height therefore put a scrollbar on a page that fit
        its viewport with room to spare, which is the "Files tab scrolls with
        the Confirm button clipped" the owner reported: nothing was too big,
        the page was merely being asked what it would like rather than what it
        needs.  ``Ignored`` is the honest answer for a list that grows into
        whatever it is handed: give it a real minimum, and no preference.
        """

        policy = widget.sizePolicy()
        policy.setVerticalPolicy(QtWidgets.QSizePolicy.Ignored)
        widget.setSizePolicy(policy)
        widget.setMinimumHeight(4 * self.layout_unit)

    @staticmethod
    def _scrollable(widget) -> QtWidgets.QScrollArea:
        """Let a control tab scroll instead of crushing its own contents.

        On the owner's smaller screen the left pane is shorter than the
        controls it holds. Scrolling is honest; overlapping is not.
        """

        area = QtWidgets.QScrollArea()
        area.setWidgetResizable(True)
        area.setFrameShape(QtWidgets.QFrame.NoFrame)
        area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        area.setWidget(widget)
        return area

    @staticmethod
    def _form_layout(parent) -> QtWidgets.QFormLayout:
        """A form whose labels wrap instead of clipping at a narrow width."""

        form = QtWidgets.QFormLayout(parent)
        form.setRowWrapPolicy(QtWidgets.QFormLayout.WrapLongRows)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.ExpandingFieldsGrow)
        form.setLabelAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(6)
        return form

    def _build_files_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        # No dashed "DROP SIF FILES HERE" panel.  It taught a gesture to the one
        # person who does not need teaching — owner, 2026-08-18: "I don't need
        # reminding. I built and am building GUI that accepts drops. And
        # others.. No others to teach. And that's intuitive in any case." — and
        # it spent a hundred pixels of the shortest column doing it.  The empty
        # table is its own invitation; the drop path is unchanged and still
        # takes files anywhere on the window.
        self.open_folder_button = QtWidgets.QPushButton("Open calibration folder…")
        self._explainable(
            self.open_folder_button,
            "Opening another calibration folder without restarting",
            "Picks up the whole bench and puts it down at another folder, "
            "exactly as launching echelle-calib there would: the snapshot and "
            "settings roots are derived inside it, the snapshot identity is "
            "re-dated from its name, and every SIF in it is loaded and triaged "
            "through the ordinary path so the roles suggest themselves as they "
            "always do. The current session is cleared — files, roles, anchors, "
            "factors, generated settings — and the bench asks first when any of "
            "that was not saved. The pattern and wavelength tables you chose "
            "stay chosen: they are your decision, not the folder's, and the "
            "band guard judges them against the new sphere as it does now.",
        )
        layout.addWidget(self.open_folder_button)

        button_row = QtWidgets.QHBoxLayout()
        self.add_files_button = QtWidgets.QPushButton("Add SIF files…")
        self.remove_file_button = QtWidgets.QPushButton("Remove selected")
        button_row.addWidget(self.add_files_button, 2)
        button_row.addWidget(self.remove_file_button, 1)
        layout.addLayout(button_row)

        self.confirm_roles_button = QtWidgets.QPushButton("Confirm suggested roles")
        self._explainable(
            self.confirm_roles_button,
            "What is left to confirm",
            "A folder whose filenames say exactly one thing is assigned the "
            "moment it finishes loading, and the bench says so in one line. "
            "This button is for the rest: a drop holding a name the bench "
            "cannot read, or two files claiming the same role, is left "
            "entirely to you, and every row of it is marked SUGGESTED until "
            "you say what it is. Pressing this assigns each unambiguous one "
            "at once; change any of them in the Role column afterwards.",
        )
        layout.addWidget(self.confirm_roles_button)

        self.file_table = QtWidgets.QTableWidget(0, 3)
        self.file_table.setHorizontalHeaderLabels(["File · triage", "Role", "Lamp"])
        header = self.file_table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        # The role is state the operator must read exactly, so its column is
        # sized by its content and never squeezed into an ellipsis.
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)
        # Mouse tracking is deliberately off: nothing in this window reacts to
        # the pointer merely passing over it (F17 item 5).
        self.file_table.setMouseTracking(False)
        self.file_table.verticalHeader().setVisible(False)
        self.file_table.verticalHeader().setDefaultSectionSize(34)
        self.file_table.verticalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents
        )
        # The file cell carries the verdict; it wraps rather than ellipsising,
        # because the part that gets cut is always the part that mattered.
        self.file_table.setWordWrap(True)
        self.file_table.setTextElideMode(QtCore.Qt.ElideNone)
        self.file_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.file_table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.file_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self._fills_its_share(self.file_table)
        layout.addWidget(self.file_table, 1)

        # "Open" is what you do to a file; this loads the selected acquisition
        # into the fit, which is a different verb and the one the operator is
        # actually choosing (F21 item 3).  It sits directly under the table it
        # takes its selection from — the only place it is ever pressed.
        self.show_frame_button = QtWidgets.QPushButton("Load for fitting")
        self._explainable(
            self.show_frame_button,
            "Loading the selected acquisition into the fit",
            "Takes whichever file is selected above and makes it the frame the "
            "Lamp fit tab measures — its orders, its spectrum, its raw pixels "
            "for the saturation checks. Nothing is opened in another program "
            "and nothing on disk changes. The bench also does this on its own "
            "when a lamp signal is the obvious next thing to fit; this is for "
            "when you want a different file.",
        )
        layout.addWidget(self.show_frame_button)
        self.control_tabs.addTab(self._scrollable(tab), "Files")

    def _build_status_band(self) -> None:
        """The readings that must never be behind a tab or a scrollbar.

        Bench state and the factors line used to be stacked below the file
        table *inside* the scrolling Files tab, which cost the left column some
        three hundred pixels it did not have and still left both invisible from
        every other tab.  They are readings, not controls, so they live across
        the top of the view column — where the space actually was — and they
        are in view whatever the operator is doing (F18 items 2 and 3).
        """

        band = QtWidgets.QWidget()
        band_layout = QtWidgets.QGridLayout(band)
        band_layout.setContentsMargins(10, 6, 10, 6)
        band_layout.setSpacing(10)

        status_group = QtWidgets.QGroupBox("Bench state")
        status_form = self._form_layout(status_group)
        # These two are one-line facts — an input mode and a filename — and a
        # readings strip is the wrong place for either of them to grow to four
        # wrapped lines and push the plots down the window.  They shorten in
        # the middle, where a SIF name is least informative, and the whole name
        # stays one hover away.
        self.watch_value = _ElidingLabel("no folder open — drag and drop or Add files")
        self.file_value = _ElidingLabel("no file open")
        self.file_state_value = QtWidgets.QLabel("WAITING")
        self.file_state_value.setObjectName("stateBadge")
        # Which pattern the bench is wearing, on the strip that is in view
        # whatever tab is open.  Every order the bench extracts, every anchor
        # row it reads and every factor it sums is taken off this file, and
        # until it was said here the only way to know which one was in use was
        # to remember what the launcher passed.
        self.pattern_source_value = _ElidingLabel("—")
        status_form.addRow("Input", self.watch_value)
        # The open frame and the state it is in are one fact in two halves, and
        # pairing them on one row is what buys the pattern its own — the strip
        # is a strip, and it was already as tall as the plots could afford.
        status_form.addRow("Open frame", self._frame_state_row())
        status_form.addRow("Pattern", self.pattern_source_value)
        self.bench_state_group = status_group

        # "Sphere factors" is what the view tab has always called this, and a
        # group box's title cannot wrap or elide: the long form was the single
        # widest thing on the readings strip and folded it for no reason.
        comparison_group = QtWidgets.QGroupBox("Sphere factors")
        comparison_layout = QtWidgets.QVBoxLayout(comparison_group)
        self.comparison_value = QtWidgets.QLabel(
            "NOT RUN — the sphere pair alone unblocks this."
        )
        self.comparison_value.setWordWrap(True)
        self.comparison_value.setObjectName("stateBadge")
        self._explainable(
            self.comparison_value,
            _SPHERE_FACTORS_TITLE,
            _SPHERE_FACTORS_EXPLANATION,
        )
        comparison_layout.addWidget(self.comparison_value)
        # No button on the readings strip.  Bench state, Alignment and Sphere
        # factors are things to READ; the one action they used to carry sat at
        # the opposite corner of the window from the step that calls for it and
        # existed twice over once the next-step panel could run it (owner,
        # 2026-08-16: "Why is your measure sensitivity still at a diagonal and
        # duplicated?").  The verb lives beside the step now, and nowhere else.
        self.sphere_factors_group = comparison_group

        # The alignment numbers are readings too, and they were the last ones
        # still buried in a control tab: RMS decides whether the fit is done,
        # and reading it meant opening the Lamp fit tab first (F20).
        alignment_group = QtWidgets.QGroupBox("Alignment")
        alignment_form = self._form_layout(alignment_group)
        self.alignment_state_value = QtWidgets.QLabel("WAITING FOR FRAME")
        self.alignment_state_value.setObjectName("stateBadge")
        self.anchor_count_value = QtWidgets.QLabel("0")
        self.rms_value = QtWidgets.QLabel("—")
        # One line, middle-elided if the rail is ever that narrow: word wrap
        # put the θ on a second line the row's height budget then hid, so the
        # header promised three numbers and showed two and a comma (owner
        # screenshot, 2026-08-17 — his own law: cut off is the same as absent).
        self.transform_value = _ElidingLabel("—")
        alignment_form.addRow("State", self.alignment_state_value)
        alignment_form.addRow("Anchors / RMS", self._anchor_rms_row())
        alignment_form.addRow("dx / dy / θ", self.transform_value)
        self._explainable(
            self.rms_value,
            "Fit RMS in detector pixels",
            "The root-mean-square distance between where the solved rigid "
            "transform predicts each anchored line and where its centroid "
            "actually is. Sub-pixel is what this instrument gives when the "
            "anchors reference the right lamp's catalog; several pixels means "
            "the anchors are being measured against the wrong element's lines, "
            "or that one bad anchor is dragging the solution.",
        )
        self._explainable(
            self.anchor_count_value,
            "Anchors",
            "Each anchor is one known line whose centroid was fitted on this "
            "frame. Two in different orders solve the rigid transform; more "
            "tighten it. Anchors on saturated lines are refused one by one, "
            "so a dim-series frame contributes its unsaturated lines and "
            "nothing else.",
        )
        self._explainable(
            self.transform_value,
            "Rigid detector transform",
            "How far the detector has moved since the base wavelength table "
            "was measured: a shift along the dispersion (dx), across the "
            "orders (dy), and a rotation. The saved snapshot's wavelength.txt "
            "is the base table moved by exactly this.",
        )
        self.alignment_group = alignment_group

        band.setSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum
        )
        # Wrapping readings ask Qt for the height their own text costs at the
        # width they are given.  Doing it through the size policy rather than
        # through a manual minimum is what keeps a reflowing band stable: a
        # hand-set minimum computed at one width becomes the input to the next
        # reflow, and the strip walks itself taller every pass.
        for label in band.findChildren(QtWidgets.QLabel):
            if not label.wordWrap():
                continue
            policy = label.sizePolicy()
            policy.setHeightForWidth(True)
            label.setSizePolicy(policy)
        self.status_band = band
        self.status_band_layout = band_layout
        #: Three readings.  Each keeps the width its own longest control needs
        #: and the band never shaves them.  Exposure guidance is still not one
        #: of them: it is prose that runs to several lines and belongs with the
        #: triage it explains; a strip that has to hold it is a strip that
        #: pushes the plots off the window.
        self._status_panels = (status_group, alignment_group, comparison_group)
        self._status_columns = 0
        self._reflow_status_band(len(self._status_panels))

    def _frame_state_row(self) -> QtWidgets.QWidget:
        """The open frame's name and its load state on one line."""

        row = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(self.file_value, 1)
        separator = QtWidgets.QLabel("·")
        separator.setObjectName("mutedText")
        layout.addWidget(separator)
        layout.addWidget(self.file_state_value)
        return row

    def _anchor_rms_row(self) -> QtWidgets.QWidget:
        """Anchor count and RMS on one line: two short numbers, one strip row.

        A readings strip is a strip. Two facts that are each three characters
        long do not each earn a row of the window's height.
        """

        row = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(self.anchor_count_value)
        separator = QtWidgets.QLabel("·")
        separator.setObjectName("mutedText")
        layout.addWidget(separator)
        layout.addWidget(self.rms_value, 1)
        return row

    def _measure_status_band(self) -> None:
        """Record what each reading panel needs, once the real font is on it.

        Measured at the end of ``_build_ui`` rather than while building: the
        bench sets its body font last, and a width measured before that is the
        designer's font, not the operator's.
        """

        for panel in self._status_panels:
            # Wrapping prose can be narrow; a button cannot. The demand is the
            # widest thing in the panel that has no other way to fit, never the
            # panel's full sizeHint — that hint is one long unwrapped line and
            # would fold the band into a column for no reason.
            demand = _STATUS_PANEL_FLOOR
            layout = panel.layout()
            margins = layout.contentsMargins() if layout is not None else None
            padding = (margins.left() + margins.right()) if margins is not None else 0
            for child in panel.findChildren(
                (QtWidgets.QPushButton, QtWidgets.QComboBox)
            ):
                demand = max(demand, child.sizeHint().width() + padding + 10)
            if isinstance(panel, QtWidgets.QGroupBox):
                # A group box's own title is drawn in its frame and cannot
                # wrap, shrink or elide: it is a demand like any button's.
                metrics = QtGui.QFontMetrics(panel.font())
                demand = max(demand, metrics.horizontalAdvance(panel.title()) + 40)
            panel.setMinimumWidth(int(demand))

    def _reflow_status_band(self, columns: int) -> None:
        """Re-lay the always-visible readings into *columns* columns.

        The owner's law applied to a row of panels: rather than shave the three
        panels until the compute button no longer fits its own label, the band
        folds onto a second and third row.  Narrower is fewer columns, never a
        squeezed control.
        """

        columns = max(1, min(columns, len(self._status_panels)))
        if columns == self._status_columns:
            return
        self._status_columns = columns
        layout = self.status_band_layout
        for panel in self._status_panels:
            layout.removeWidget(panel)
        for index in range(layout.columnCount()):
            layout.setColumnStretch(index, 0)
        for index, panel in enumerate(self._status_panels):
            row, column = divmod(index, columns)
            span = 1
            # A lone trailing panel takes the whole width rather than leaving a
            # hole beside itself.
            if index == len(self._status_panels) - 1 and column == 0:
                span = columns
            layout.addWidget(panel, row, column, 1, span)
        for index in range(columns):
            layout.setColumnStretch(index, 1)
        self._pin_status_band_height()

    def _pin_status_band_height(self) -> None:
        """Give the readings strip exactly the height it needs, once.

        Left to the splitter's proportions the band grew to a third of the
        window and squeezed the plots — the same hoarding, pointed the other
        way.  It is a strip: it takes what its content costs and the plots keep
        everything else, until the operator drags the handle and takes over.
        """

        band = self.status_band
        layout = band.layout()
        layout.activate()
        needed = max(layout.minimumSize().height(), layout.sizeHint().height())
        if layout.hasHeightForWidth() and band.width() > 0:
            needed = max(needed, layout.heightForWidth(band.width()))
        band.setMinimumHeight(needed)
        if "readings" in _SESSION_SPLITTER_SHARES:
            return  # the operator has cut this column themselves
        sizes = self.readings_splitter.sizes()
        if len(sizes) == 2 and sizes[0] != needed and sum(sizes) > needed:
            self.readings_splitter.setSizes([needed, sum(sizes) - needed])

    def _status_band_columns(self, width: int) -> int:
        """How many columns of readings *width* can hold without squeezing."""

        spacing = self.status_band_layout.spacing()
        margins = self.status_band_layout.contentsMargins()
        usable = width - margins.left() - margins.right()
        demands = [panel.minimumWidth() for panel in self._status_panels]
        for columns in range(len(demands), 1, -1):
            # A grid column is as wide as the widest panel that lands in it.
            per_column = [0] * columns
            for index, demand in enumerate(demands):
                column = index % columns
                per_column[column] = max(per_column[column], demand)
            if sum(per_column) + spacing * (columns - 1) <= usable:
                return columns
        return 1

    def _build_fit_bar(self) -> QtWidgets.QWidget:
        """The controls that steer the spectrum, on one line above it.

        These used to be a tall "Order and frame" group box in the left rail,
        a whole column away from the plot they steer and from the lines being
        clicked (owner, 2026-08-16: "we'd rather have a small box next to the
        action — line selection"). They are three short controls; a rail-wide
        group box was the wrong shape for them in the wrong place.
        """

        bar = QtWidgets.QWidget()
        bar.setObjectName("fitBar")
        # The bench's ordinary padding is generous on purpose — these are
        # controls to be hit, not read. On one line across the centre that
        # generosity is fatal: at the standard padding this row demands 1227 px
        # of minimum width, and a centre column with a floor that high pushes
        # both rails to their own minimums and breaks the rule that the plots
        # are the flexible element (F20). A strip is a different shape from a
        # column and gets its own, tighter, rule.
        bar.setStyleSheet(
            "#fitBar QLabel { padding: 0 2px; }"
            "#fitBar QComboBox { padding: 2px 6px; }"
            "#fitBar QToolButton { padding: 1px 6px; }"
            "#fitBar QSpinBox { padding: 2px 4px; }"
            "#fitBar QCheckBox { padding: 0 2px; }"
        )
        row = QtWidgets.QHBoxLayout(bar)
        row.setContentsMargins(0, 0, 0, 2)
        row.setSpacing(6)

        row.addWidget(QtWidgets.QLabel("Order"))
        self.previous_order_button = QtWidgets.QToolButton()
        self.previous_order_button.setText("◀")
        self.previous_order_button.setToolTip("Previous Echelle order")
        self.order_spin = QtWidgets.QSpinBox()
        self.order_spin.setRange(0, self.session.pattern.shape[1] - 1)
        self.order_spin.setToolTip(
            "Which Echelle order the spectrum plot shows and a click fits."
        )
        self.order_spin.setMaximumWidth(70)
        self.next_order_button = QtWidgets.QToolButton()
        self.next_order_button.setText("▶")
        self.next_order_button.setToolTip("Next Echelle order")
        self.order_total_value = QtWidgets.QLabel(
            f"of {self.session.pattern.shape[1] - 1}"
        )
        row.addWidget(self.previous_order_button)
        row.addWidget(self.order_spin)
        row.addWidget(self.next_order_button)
        row.addWidget(self.order_total_value)

        row.addSpacing(10)
        row.addWidget(QtWidgets.QLabel("Lines"))
        self.line_family_combo = QtWidgets.QComboBox()
        self.line_family_combo.addItems(list(KNOWN_LAMP_NAMES))
        self.line_family_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToMinimumContentsLength
        )
        self.line_family_combo.setMinimumContentsLength(5)
        row.addWidget(self.line_family_combo)
        self._explainable(
            self.line_family_combo,
            "Which catalog fills the expected-lines table",
            "This follows the assigned lamp on its own: assign a Ne lamp and "
            "the expected-lines panel fills with neon. Changing it by hand is "
            "an override for comparing one lamp's catalog against another "
            "frame; the anchors themselves always reference the assigned "
            "lamp's own rows, and the panel warns when the two disagree.",
        )

        row.addSpacing(10)
        # The owner's own name for it, and his own placement: this changes how
        # the detector reads, so it belongs in the strip that steers the view,
        # not parked in a tab a column away. "1:1" costs almost no width, which
        # is what let it come back.
        self.equal_aspect_check = QtWidgets.QToolButton()
        # "1:1" alone was unfindable — the owner went looking for the aspect
        # control and reported it gone while it sat right there. His own words
        # for it, on the button.
        self.equal_aspect_check.setText("Aspect 1:1")
        self.equal_aspect_check.setCheckable(True)
        self._explainable(
            self.equal_aspect_check,
            "Showing the detector at its true shape",
            "The detector view stretches to fill the space it is given, which "
            "is what makes the order traces far enough apart to aim at. Press "
            "this and one detector pixel is drawn square in both directions, "
            "so the frame appears at its real 2560x2160 proportions. Useful "
            "for judging the geometry; unhelpful for clicking lines, which is "
            "why it is off unless you ask.",
        )
        row.addWidget(self.equal_aspect_check)
        row.addStretch(1)
        return bar

    def _build_lamp_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        # The order stepper and the lamp catalog went to the strip above the
        # spectrum, where the clicking happens.  These two did not: choosing
        # which frame to fit is a decision made once per acquisition, and
        # squaring the pixels is a way of looking rather than a way of working.
        frame_row = QtWidgets.QHBoxLayout()
        frame_row.addWidget(QtWidgets.QLabel("Fit on"))
        self.frame_combo = QtWidgets.QComboBox()
        self.frame_combo.addItem("Mean of all frames", None)
        frame_row.addWidget(self.frame_combo, 1)
        layout.addLayout(frame_row)

        self.frame_choice_value = QtWidgets.QLabel("no acquisition open")
        self.frame_choice_value.setWordWrap(True)
        self.frame_choice_value.setObjectName("benchHelp")
        self._explainable(
            self.frame_choice_value,
            "Single frame or the mean of the acquisition",
            "A 3-frame acquisition is shot so the frames can be averaged: the "
            "mean is quieter and is the right thing to fit most of the time. "
            "A single frame is the honest choice when one frame is spoiled by "
            "a cosmic ray, or when a line clips in one frame and not the "
            "others — saturation is checked on the raw pixels of whichever "
            "frame you are fitting. Changing this drops the anchors already "
            "collected, because a centroid belongs to the spectrum it was "
            "measured on.",
        )
        layout.addWidget(self.frame_choice_value)

        # The button itself lives on the anchor panel in the right rail, where
        # the table it fills is and where every tab can see it.  What stays
        # here is the reading behind it: which rows it would measure, and whose
        # vetting they carry.
        self.auto_anchor_value = QtWidgets.QLabel("no acquisition open")
        self.auto_anchor_value.setWordWrap(True)
        self.auto_anchor_value.setObjectName("benchHelp")
        layout.addWidget(self.auto_anchor_value)

        # The anchor table and the alignment readings it produces live in the
        # right rail (F20).  Stacked here they were the single heaviest thing
        # in the controls column, which is what made this tab scroll at every
        # window size the bench has ever opened at.
        self.message_value = QtWidgets.QLabel("Waiting for data.")
        self.message_value.setWordWrap(True)
        self.message_value.setObjectName("messagePanel")
        layout.addWidget(self.message_value)
        layout.addStretch(1)
        self.control_tabs.addTab(self._scrollable(tab), "Lamp fit")

    def _build_save_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        identity_group = QtWidgets.QGroupBox("Snapshot identity")
        form = self._form_layout(identity_group)
        self.snapshot_id_edit = QtWidgets.QLineEdit(self.initial_snapshot_id)
        self.detector_edit = QtWidgets.QLineEdit(self.initial_detector)
        self.base_snapshot_edit = QtWidgets.QLineEdit(self.initial_base_snapshot)
        self.notes_edit = QtWidgets.QLineEdit()
        form.addRow("ID", self.snapshot_id_edit)
        form.addRow("Detector", self.detector_edit)
        form.addRow("Base snapshot", self.base_snapshot_edit)
        form.addRow("Notes", self.notes_edit)
        layout.addWidget(identity_group)
        self._describe_snapshot_identity()

        # Where the files go, in full.  The bare folder names that used to be
        # here read as a promise that the bench wrote beside the data, and the
        # bench had no such promise to make: it wrote under whatever directory
        # the launcher started in.  The whole path is shown, shortened in the
        # middle only when it will not fit, and is one hover — or one click
        # into the Why dock — away from being read exactly.
        self.snapshot_root_value = _ElidingLabel("")
        self.snapshot_root_value.setObjectName("mutedText")
        layout.addWidget(self.snapshot_root_value)
        self.config_root_value = _ElidingLabel("")
        self.config_root_value.setObjectName("mutedText")
        layout.addWidget(self.config_root_value)
        # Both readings, and both explanations, are written in one place: the
        # roots move when another calibration folder is opened, and a label
        # filled in at build time only would have gone on naming the first one.
        self._describe_output_roots()
        # "TOML" names the file format, which is not what the operator is
        # deciding to do (F21 item 6).  What this step saves is the alignment
        # and the campaign around it; that the files happen to be commented
        # TOML is a detail for whoever opens them afterwards.
        self.generate_tomls_button = QtWidgets.QPushButton("Save alignment settings")
        self._explainable(
            self.generate_tomls_button,
            "Saving the alignment and the campaign around it",
            "Writes the solved alignment, the roles you assigned, and the "
            "campaign's own settings as commented text files you can open and "
            "edit by hand. They are the input to the snapshot below, so this "
            "comes first. Existing files are never overwritten silently — if "
            "they are already there the bench says so and offers to redo them.",
        )
        self.save_snapshot_button = QtWidgets.QPushButton("Save and validate snapshot")
        layout.addWidget(self.generate_tomls_button)
        # Refusing to clobber is right; refusing with no way forward is what
        # left the owner stuck at "already exists" (F21 item 11).
        self.regenerate_tomls_button = QtWidgets.QPushButton("Regenerate")
        self.regenerate_tomls_button.setVisible(False)
        self._explainable(
            self.regenerate_tomls_button,
            "Rewriting settings files that already exist",
            "Appears only when the files are already on disk. It writes them "
            "again from the bench's current state, replacing what is there. "
            "Any edits made by hand in those files are lost, which is why it "
            "is a separate, deliberate press rather than something the first "
            "button does quietly.",
        )
        layout.addWidget(self.regenerate_tomls_button)
        layout.addWidget(self.save_snapshot_button)
        self.save_state_value = QtWidgets.QLabel("NOT READY")
        self.save_state_value.setObjectName("stateBadge")
        # A confirmation that names a folder and leaves the operator to go find
        # it is half a confirmation, and on a UNC path it is the half that
        # costs the most typing.  The offer appears only once there is a real
        # saved folder to open.
        confirmation = QtWidgets.QHBoxLayout()
        confirmation.addWidget(self.save_state_value)
        self.open_snapshot_button = QtWidgets.QPushButton("Open folder")
        self.open_snapshot_button.setVisible(False)
        self.open_snapshot_button.clicked.connect(self._open_saved_snapshot_folder)
        confirmation.addWidget(self.open_snapshot_button)
        confirmation.addStretch(1)
        layout.addLayout(confirmation)
        # The bench's message line lives on the Lamp fit tab, and a save is
        # pressed from this one.  A refusal set on a widget the operator cannot
        # see is a silent refusal: the owner pressed Save on a signal-only 2019
        # folder, was told "sphere/lamp pairs … are required" one tab away, and
        # reported that nothing happened at all.  Whatever the save says, it
        # says here too, where the button that caused it is.
        #
        # One line that shortens itself, for the same reason the next-step
        # strip is one: a wrapping label on this tab is charged its full
        # wrapped height, and three lines of refusal put the Save tab into the
        # scrollbar this bench was rebuilt to abolish.  The whole sentence is
        # the tooltip, and the message line on the Lamp fit tab carries it in
        # full.
        self.save_message_value = _ElidingLabel("Nothing has been saved yet.")
        self.save_message_value.setObjectName("benchHelp")
        layout.addWidget(self.save_message_value)
        self.toml_preview = QtWidgets.QPlainTextEdit()
        self.toml_preview.setReadOnly(True)
        self.toml_preview.setPlaceholderText(
            "Generated campaign.toml appears here; all files remain ordinary and editable."
        )
        self._fills_its_share(self.toml_preview)
        layout.addWidget(self.toml_preview, 1)
        self.control_tabs.addTab(self._scrollable(tab), "Save")

    def _describe_snapshot_identity(self) -> None:
        """Say where the identity's date came from — especially when nowhere.

        The identity is the calibration's date, and the date is the day the
        frames were *taken*: today is only ever a placeholder for it. So the
        one state that must be unmissable is the one where nothing could be
        derived, because then — and only then — the date on screen is the
        bench's bookkeeping rather than the measurement's own day.
        """

        if self.snapshot_date is None and self._snapshot_id_decided:
            detail = (
                "This identity carries no leading YYYYMMDD, so nothing here "
                "states an acquisition date and the saved epoch has no day to "
                "start on. Give it the day these frames were taken."
            )
        elif self.snapshot_date is None:
            detail = (
                "No acquisition date could be read — not from the folder name "
                "the bench was launched at, and not from the frames' own SIF "
                "headers — so this is TODAY'S date standing in for it. CHECK "
                "it against the day these frames were taken before saving, and "
                "correct it here if it is wrong."
            )
        else:
            detail = (
                f"Dated {self.snapshot_date:%Y-%m-%d}, read from "
                f"{self._snapshot_date_source}. Loading frames may still fill "
                "this in from their headers while it is untouched; typing here "
                "settles it and the bench never overwrites what you typed."
            )
        self._explainable(
            self.snapshot_id_edit,
            "The snapshot identity, and the day it is dated by",
            "A calibration is dated by the day its images were TAKEN, not by "
            "the day it is computed: YYYYMMDD_<detector>, matching "
            f"20250926_cmos and the other epochs in the registry. {detail}",
        )

    def _snapshot_id_hand_edited(self, _text: str = "") -> None:
        """A typed identity is a decision, and outranks anything derived."""

        self._snapshot_id_decided = True
        self._snapshot_date_source = "the identity you typed here"
        self.snapshot_date = acquisition_date_from_name(
            self.snapshot_id_edit.text().strip()
        )
        self._describe_snapshot_identity()

    def _adopt_acquisition_date(self, frame) -> None:
        """Take the identity's date from the frames when nothing else supplied one.

        Only ever fills a gap: a date already read from the folder name stands,
        a flag stands, and a hand-typed identity stands absolutely. What this
        catches is the folder whose name says nothing — the frames themselves
        still know which day they were shot on.
        """

        if self.snapshot_date is not None or self._snapshot_id_decided:
            return
        if self.snapshot_id_edit.text().strip() != self.initial_snapshot_id:
            return
        day = acquisition_date_from_metadata(getattr(frame, "metadata", None))
        if day is None:
            return
        self.snapshot_date = day
        self._snapshot_date_source = "the loaded frames' own SIF headers"
        self.initial_snapshot_id = snapshot_identity(
            day, self.detector_edit.text().strip() or self.initial_detector
        )
        # setText, not setPlaceholder: this is the identity that will be saved,
        # and it is still the operator's to overwrite.  ``textEdited`` does not
        # fire for a programmatic change, so this cannot mark itself decided.
        self.snapshot_id_edit.setText(self.initial_snapshot_id)
        self._describe_snapshot_identity()

    def _build_triage_view(self) -> None:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        # The folder before the file (F21 item 4a).  Triage used to be a
        # single-file surface: to judge six acquisitions the operator clicked
        # six rows and held five verdicts in their head.  One compact row per
        # file — its role, its verdict in the verdict's own colour, its peak as
        # a share of full scale, its anomaly count — answers "is this folder
        # any good" in one look, and clicking a row opens exactly the detailed
        # reading that was already here, underneath.
        self.triage_summary_table = QtWidgets.QTableWidget(0, 5)
        self.triage_summary_table.setHorizontalHeaderLabels(
            ["File", "Role", "Verdict", "Peak", "Anomalies"]
        )
        summary_header = self.triage_summary_table.horizontalHeader()
        summary_header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        for column in range(1, 5):
            summary_header.setSectionResizeMode(
                column, QtWidgets.QHeaderView.ResizeToContents
            )
        self.triage_summary_table.verticalHeader().setVisible(False)
        self.triage_summary_table.verticalHeader().setDefaultSectionSize(
            self.layout_unit + 6
        )
        # One line per file is the whole point: a wrapped row is not a glance.
        self.triage_summary_table.setWordWrap(False)
        self.triage_summary_table.setTextElideMode(QtCore.Qt.ElideMiddle)
        self.triage_summary_table.setMouseTracking(False)
        self.triage_summary_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectRows
        )
        self.triage_summary_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SingleSelection
        )
        self.triage_summary_table.setEditTriggers(
            QtWidgets.QAbstractItemView.NoEditTriggers
        )
        self.triage_summary_table.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOff
        )
        # It summarises; it never takes room from the reading it introduces.
        self.triage_summary_table.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum
        )
        self.triage_summary_table.setMinimumHeight(3 * self.layout_unit)
        self._fit_triage_summary_height()
        self._explainable(
            self.triage_summary_table,
            "The whole folder at a glance",
            "One row per loaded file: the role it carries, the verdict its "
            "exposure earns in that role, its brightest honest pixel as a "
            "share of full scale, and how many isolated full-scale pixels "
            "(cosmic rays, hot pixels) were counted rather than held against "
            "it. Click a row to read that file in full below.",
        )
        layout.addWidget(self.triage_summary_table, 0)

        self.triage_headline = QtWidgets.QLabel(
            "Drop a SIF onto the bench. Triage needs nothing but a file."
        )
        self.triage_headline.setObjectName("triageHeadline")
        self.triage_headline.setWordWrap(True)
        layout.addWidget(self.triage_headline)
        # The four-line breakdown that used to sit here now folds into the Why
        # dock: the verdict decides the next action, the arithmetic behind it
        # is one hover away.
        self.triage_next_value = QtWidgets.QLabel(
            "Shoot → drop the file → one glance → adjust the lamp → shoot again."
        )
        self.triage_next_value.setWordWrap(True)
        self.triage_next_value.setObjectName("messagePanel")
        layout.addWidget(self.triage_next_value)

        # Exposure guidance explains this view's verdict, so it reads here,
        # beside the histograms it is derived from, rather than in the file
        # list where it used to take four wrapped lines of the controls column.
        self.exposure_value = QtWidgets.QLabel(
            "Exposure guidance for the selected file appears here."
        )
        self.exposure_value.setWordWrap(True)
        self.exposure_value.setObjectName("messagePanel")
        layout.addWidget(self.exposure_value)

        self.histogram_widget = pg.GraphicsLayoutWidget()
        self.histogram_widget.setBackground("#10151b")
        self.histogram_plot = self.histogram_widget.addPlot(
            row=0, col=0, title="Raw counts histogram"
        )
        self.histogram_plot.setLabel("bottom", "counts")
        self.histogram_plot.setLabel("left", "pixels per bin")
        self.histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.histogram_plot.setLogMode(y=True)
        self.histogram_widget.setMinimumHeight(BENCH_HISTOGRAM_LINES * self.layout_unit)
        layout.addWidget(self.histogram_widget, 1)

        # The top end lives in its own widget so that the honest answer — "no
        # pixels are up here" — can take its place in words (F17 item 3). An
        # empty log histogram drew a solid block, which says nothing at all.
        #
        # It is a strip, not a second histogram (F21 item 4b).  Near-saturation
        # distribution is the number the lamp is adjusted by, so it stays — but
        # it was drawn the same size as the raw-counts histogram and split the
        # view's attention in half.  A fixed strip a fraction of the primary's
        # floor makes the subordination structural rather than a matter of
        # which one happens to be looked at first.
        self.top_histogram_widget = pg.GraphicsLayoutWidget()
        self.top_histogram_widget.setBackground("#10151b")
        self.top_histogram_plot = self.top_histogram_widget.addPlot(
            row=0, col=0, title="Top end — the last 10% before full scale"
        )
        self.top_histogram_plot.setLabel("bottom", "counts")
        # No left label: on a strip this thin the axis title is most of the
        # panel, and the count of pixels per bin is read off the primary.
        self.top_histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.top_histogram_plot.setLogMode(y=True)
        self.top_histogram_widget.setMaximumHeight(
            BENCH_TOP_END_LINES * self.layout_unit
        )
        self.top_histogram_widget.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum
        )
        layout.addWidget(self.top_histogram_widget, 0)
        self.top_end_message = QtWidgets.QLabel(_TOP_END_EMPTY)
        self.top_end_message.setWordWrap(True)
        self.top_end_message.setObjectName("messagePanel")
        self.top_end_message.setVisible(False)
        self._explainable(
            self.top_end_message,
            "The top end of the counts histogram",
            "This panel bins only the last 10% of the detector's range, where "
            "clipping happens. When nothing is up there the frame has headroom "
            "and there is no distribution to draw, so it says so rather than "
            "drawing an empty box that looks like a full one.",
        )
        layout.addWidget(self.top_end_message)
        self.view_tabs.addTab(widget, "Exposure triage")

    def _build_lamp_fit_view(self) -> None:
        """The spectrum surface of the line workflow.

        Its expected-line list lives down the tall left column rather than
        crammed under the plot (F17 item 1), but it is the same list: the
        sticks drawn here and the rows listed there come from one call.
        """

        #: Two pens, made once.  An order change re-points them; it never
        #: allocates a pen per stick per refresh.
        self._expected_pen = pg.mkPen("#f5b95f", width=1, style=QtCore.Qt.DashLine)
        self._anchored_pen = pg.mkPen("#6ee7b7", width=1.6)
        widget = QtWidgets.QWidget()
        outer = QtWidgets.QVBoxLayout(widget)
        outer.setContentsMargins(0, 0, 0, 0)

        # "The active filename is stated in the fit view itself" (F16 item 5) —
        # which is here, beside the spectrum it describes, rather than in the
        # control tab where it was four wrapped lines of a narrow rail's height
        # and nowhere near the plot it is talking about.
        self.fit_file_value = QtWidgets.QLabel("no file open for fitting")
        self.fit_file_value.setWordWrap(True)
        self.fit_file_value.setObjectName("messagePanel")
        self._explainable(
            self.fit_file_value,
            "Which file this fit is measuring",
            "The fit view opens the assigned lamp signal by default, and states "
            "here which file it is showing and whether the assigned lamp "
            "background is being subtracted from it. Opening any other file "
            "from the Files tab overrides this and the line above says so.",
        )
        outer.addWidget(self.fit_file_value)
        self.fit_warning_value = QtWidgets.QLabel("")
        self.fit_warning_value.setWordWrap(True)
        self.fit_warning_value.setObjectName("warningPanel")
        self.fit_warning_value.setVisible(False)
        outer.addWidget(self.fit_warning_value)

        # Which lamp catalog is scoping the sticks reads beside the sticks, for
        # the same reason: it is a statement about this spectrum, and in the
        # control rail it was three wrapped lines of a column that has no
        # spare ones.
        self.reference_value = QtWidgets.QLabel(
            "No lamp catalog is scoping the fit yet."
        )
        self.reference_value.setWordWrap(True)
        self.reference_value.setObjectName("messagePanel")
        outer.addWidget(self.reference_value)

        self.fit_bar = self._build_fit_bar()

        split = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        split.setChildrenCollapsible(False)
        self.lamp_fit_splitter = split
        outer.addWidget(split, 1)

        graphics = pg.GraphicsLayoutWidget()
        graphics.setBackground("#10151b")

        self.detector_plot = graphics.addPlot(row=0, col=0, title="Detector + order traces")
        self.detector_plot.setLabel("bottom", "detector column", units="px")
        self.detector_plot.setLabel("left", "detector row", units="px")
        self.detector_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.detector_plot.getAxis("left").enableAutoSIPrefix(False)
        # Detector row increases upward, so the order control and the image
        # agree: order 0 sits at row 56 and order 28 at row 2094, and stepping
        # the control up moves the highlighted trace up (F21 item 9, the owner
        # twice).  Inverting the axis was what put them in opposition, and the
        # order numbers stay exactly as the wavelength table spells them —
        # which is the constraint the owner set on fixing this.  The main GUI
        # never inverted its own view, so the two tools now match as well.
        self.detector_plot.invertY(False)
        # Equal aspect on request: one detector pixel square in both directions
        # (F21 item 8, "not mandatory, but wanted"), off by default because a
        # 2560x2160 frame squared up leaves the order traces too thin to aim at.
        self.detector_plot.setAspectLocked(False)
        self.detector_image = pg.ImageItem(axisOrder="col-major")
        self.detector_plot.addItem(self.detector_image)

        # The spectrum lives in its own graphics widget so the fit strip can
        # sit between it and the detector image — directly above the trace
        # whose lines are being clicked (owner, 2026-08-16: "the order scroll
        # belongs to the bottom, next to the lines I'm supposed to click").
        # One widget for both plots left no seam to put it in.  The residual
        # bars used to sit under this spectrum and are now in the anchors panel
        # of the right rail, beside the rows they score; the whole of the height
        # they cost here goes back to the spectrum, which is the surface the
        # lines are actually clicked on.
        lower_graphics = pg.GraphicsLayoutWidget()
        lower_graphics.setBackground("#10151b")
        self.order_plot = lower_graphics.addPlot(
            row=0, col=0, title="Selected order spectrum"
        )
        self.order_plot.setLabel("bottom", "raw detector column", units="px")
        self.order_plot.setLabel("left", "mean extracted counts")
        self.order_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.order_plot.getAxis("left").enableAutoSIPrefix(False)
        self.order_curve = self.order_plot.plot(
            pen=pg.mkPen("#76d6ff", width=1.2), antialias=False
        )
        self.order_plot.setDownsampling(auto=True, mode="peak")
        self.order_plot.setClipToView(True)
        # The stick the expected-lines panel points at.  One item, moved — not
        # a new item per selection.
        self.line_highlight = pg.InfiniteLine(
            0.0, angle=90, movable=False, pen=pg.mkPen("#ffe08a", width=2.4)
        )
        self.line_highlight.setVisible(False)
        self.order_plot.addItem(self.line_highlight, ignoreBounds=True)
        # One scatter item for every accepted anchor of the shown order, fed
        # with setData rather than rebuilt one item per anchor.
        self._anchor_scatter = pg.ScatterPlotItem(
            size=11, brush=pg.mkBrush("#6ee7b7"), pen=pg.mkPen("#d7fff1", width=1)
        )
        self.order_plot.addItem(self._anchor_scatter)

        split.addWidget(graphics)
        lower = QtWidgets.QWidget()
        lower_layout = QtWidgets.QVBoxLayout(lower)
        lower_layout.setContentsMargins(0, 0, 0, 0)
        lower_layout.setSpacing(2)
        lower_layout.addWidget(self.fit_bar)
        lower_layout.addWidget(lower_graphics, 1)
        split.addWidget(lower)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 2)
        self.view_tabs.addTab(widget, "Lamp fit")

    def _build_expected_lines_panel(self) -> None:
        """The expected-line list, down the right rail (F17 item 1, F20).

        Under the spectrum the table showed four rows of a twenty-row list; in
        the left column below the controls it was still sharing one column's
        height with four tabs of controls.  In the tables rail it is as long as
        the window, beside the anchor table it feeds, and in view whichever
        control tab is open.
        """

        panel = QtWidgets.QGroupBox("Expected lines")
        panel_layout = QtWidgets.QVBoxLayout(panel)
        panel_layout.setContentsMargins(6, 6, 6, 4)
        panel_layout.setSpacing(5)
        self.line_panel_header = QtWidgets.QLabel(
            "Expected lines — assign a lamp and this fills itself."
        )
        self.line_panel_header.setWordWrap(True)
        self.line_panel_header.setObjectName("messagePanel")
        self._explainable(
            self.line_panel_header,
            "The expected-line panel fills itself",
            "These are the assigned lamp's own wavelength-table rows that fall "
            "inside the selected order — the same list, line for line, as the "
            "labelled sticks on the spectrum, and the same rows a click can "
            "anchor. Selecting a row marks its stick; a row already anchored "
            "on this frame is ticked, and clicking its stick takes the anchor "
            "back off. There is nothing to populate by hand.",
        )
        panel_layout.addWidget(self.line_panel_header)
        self.line_help_table = QtWidgets.QTableWidget(0, 5)
        self.line_help_table.setHorizontalHeaderLabels(
            ["Line", "λ nm", "Pixel", "Rel. I", "Anchor"]
        )
        header = self.line_help_table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        for column in (1, 2, 3, 4):
            header.setSectionResizeMode(column, QtWidgets.QHeaderView.ResizeToContents)
        self.line_help_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.line_help_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.line_help_table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.line_help_table.setAlternatingRowColors(True)
        self.line_help_table.verticalHeader().setVisible(False)
        # The same row height as the anchor table beside it: two working lists
        # of the same kind of thing, and the denser row buys readable rows in
        # the rail rather than padding.
        self.line_help_table.verticalHeader().setDefaultSectionSize(24)
        panel_layout.addWidget(self.line_help_table, 1)
        self.expected_lines_panel = panel
        self.tables_splitter.addWidget(panel)

    def _build_sphere_view(self) -> None:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        self.sphere_view_message = QtWidgets.QLabel(
            "Classify the integrating-sphere signal/background pair, then compute factors."
        )
        self.sphere_view_message.setWordWrap(True)
        self.sphere_view_message.setObjectName("messagePanel")
        header = QtWidgets.QHBoxLayout()
        header.addWidget(self.sphere_view_message, 1)
        # The next-step slot can only offer the measurement while its checklist
        # row is the first blocking one: once the comparison is done, or a row
        # above it is not, nothing on the bench could run it again.  It is read
        # here, so the control that reruns it is here — beside the curves it
        # redraws, and not on the readings strip, which stays button-free.
        self.recompute_sphere_button = QtWidgets.QPushButton("Recompute")
        self.recompute_sphere_button.setEnabled(False)
        self._explainable(
            self.recompute_sphere_button,
            "Measuring the sensitivity again",
            "Recomputes this campaign's absolute factors from the sphere pair "
            "currently assigned and compares them with the previous campaign's "
            "again. Use it after correcting which files carry the sphere roles, "
            "or to repeat a comparison that has already been made — the "
            "procedure's own next step offers this only while it is the step "
            "you are on.",
        )
        header.addWidget(self.recompute_sphere_button)
        layout.addLayout(header)

        # The pattern row, directly under the reading that reports the pattern.
        # The band guard is read here — the factor curves are summed over these
        # very traces — so this is where the two things that can answer it live,
        # by the same rule that put Recompute above rather than on the readings
        # strip.  Quiet by default: one line and the picker; the extraction
        # surfaces as the loud verb only when the guard is actually warning, and
        # the procedure's next-step button offers it then too.
        pattern_row = QtWidgets.QHBoxLayout()
        self.pattern_choice_value = _ElidingLabel("Pattern: —")
        self.pattern_choice_value.setObjectName("mutedText")
        pattern_row.addWidget(self.pattern_choice_value, 1)
        self.choose_pattern_button = QtWidgets.QPushButton("Choose pattern file…")
        self._explainable(
            self.choose_pattern_button,
            "Standing the bench on a different pattern",
            "Picks the order-pattern table the whole bench reads: the traces "
            "drawn over the detector image, the rows every anchor's centroid is "
            "paired with, and the bands the factor curves are summed over. The "
            "open frame is extracted again on it, the sphere's band offsets are "
            "measured against it, and anything derived from the old pattern — "
            "factors, settings bundle, save state — is dropped rather than "
            "quietly carried over. Anchors go with it: they were fitted on the "
            "old geometry.",
        )
        pattern_row.addWidget(self.choose_pattern_button)
        self.extract_pattern_button = QtWidgets.QPushButton(
            "Extract pattern from this sphere"
        )
        self._explainable(
            self.extract_pattern_button,
            "Fitting a pattern to the sphere in front of you",
            "Runs the same fit echelle-pattern runs — the sphere signal minus "
            "its background, order bands traced with the current pattern as the "
            "prior — writes the table under the campaign's own patterns folder, "
            "and stands the bench on it without closing. This is the answer to "
            "a band guard that says the light does not sit on the chosen "
            "pattern: the sphere is the frame that lights every order, so it is "
            "the frame the geometry is measured from.",
        )
        pattern_row.addWidget(self.extract_pattern_button)
        layout.addLayout(pattern_row)

        # The other half of the comparison, on the same terms as the pattern:
        # the previous pair arrived as a command line, and a command line is a
        # launch convenience, not a place to change an input from (owner,
        # 2026-08-18: "CLI OR GUI. Not both").  It sits beside the curve it
        # draws and the ratio it is the denominator of.
        previous_row = QtWidgets.QHBoxLayout()
        self.previous_pair_value = _ElidingLabel("Previous pair: —")
        self.previous_pair_value.setObjectName("mutedText")
        previous_row.addWidget(self.previous_pair_value, 1)
        self.choose_previous_button = QtWidgets.QPushButton("Choose previous pair…")
        self._explainable(
            self.choose_previous_button,
            "Comparing against a different previous campaign",
            "Picks the previous campaign's integrating-sphere pair this "
            "campaign's factors are measured against — the amber curve, and "
            "the denominator of the median ratio beside it. Choose the signal "
            "frame; its background is offered from the file naming beside it "
            "(-bg, _bg, -background) and asked for only when the folder does "
            "not answer. The comparison is dropped the moment the pair "
            "changes, because a ratio measured against the old pair is not a "
            "reading about the new one: Recompute is the next press.",
        )
        previous_row.addWidget(self.choose_previous_button)
        layout.addLayout(previous_row)

        self.sphere_plot = pg.PlotWidget(title="Absolute calibration factors")
        self.sphere_plot.setBackground("#10151b")
        self.sphere_plot.setLabel("bottom", "wavelength", units="nm")
        self.sphere_plot.setLabel("left", "factor", units="W m⁻² sr⁻¹ nm⁻¹ count⁻¹")
        self.sphere_plot.addLegend()
        # One PlotDataItem per curve, created once and fed with setData.  A
        # 42k-sample factor curve is only heavy when it is drawn point for
        # point: downsampled to the view, clipped to it, and unantialiased, the
        # same curve pans at interactive speed.
        self.candidate_curve = self._factor_curve("#70d6ae", 1.6, "new measured pair")
        self.previous_curve = self._factor_curve("#f5b95f", 1.2, "previous campaign")
        # Downsampling is set on the plot, not on the curves: PlotItem pushes
        # its own control-panel state onto every item whenever the log mode
        # changes, which silently undoes a per-curve setting.
        self.sphere_plot.setDownsampling(auto=True, mode="peak")
        # No clipToView here, and it is not a matter of taste.  Clipping finds
        # the visible slice with a binary search, which assumes x ascends.  A
        # factor curve is the echelle orders laid end to end, so its wavelength
        # axis walks *backwards* at nearly every sample (42459 of 42739 steps on
        # the real previous-campaign pair).  The search then returns a slice of
        # about one point and the whole curve disappears at any zoom, while the
        # legend still names it — the owner's vanishing orange curve.  The
        # speed came from downsampling anyway: on that data clipping is worth
        # about a tenth of a millisecond per pan frame (F18 item 5).
        self.sphere_plot.setClipToView(False)
        self.sphere_plot.setLogMode(y=True)
        layout.addWidget(self.sphere_plot)
        self.view_tabs.addTab(widget, "Sphere factors")

    def _factor_curve(self, color: str, width: float, name: str) -> pg.PlotDataItem:
        """Add one persistent factor curve to the plot."""

        curve = pg.PlotDataItem(
            pen=pg.mkPen(color, width=width),
            name=name,
            antialias=False,
            connect="finite",
        )
        self.sphere_plot.addItem(curve)
        return curve

    def _connect_ui(self) -> None:
        self.order_spin.valueChanged.connect(self._order_changed)
        self.previous_order_button.clicked.connect(
            lambda: self.order_spin.setValue(self.order_spin.value() - 1)
        )
        self.next_order_button.clicked.connect(
            lambda: self.order_spin.setValue(self.order_spin.value() + 1)
        )
        self.frame_combo.currentIndexChanged.connect(self._frame_choice_changed)
        self.line_family_combo.currentTextChanged.connect(self._line_family_changed)
        self.order_plot.scene().sigMouseClicked.connect(self._order_plot_clicked)
        self.auto_anchor_button.clicked.connect(self._auto_anchor)
        self.equal_aspect_check.toggled.connect(self.detector_plot.setAspectLocked)
        self.remove_button.clicked.connect(self._remove_selected_anchor)
        self.clear_button.clicked.connect(self._clear_anchors)
        self.open_folder_button.clicked.connect(self._pick_calibration_folder)
        self.add_files_button.clicked.connect(self._pick_files)
        self.remove_file_button.clicked.connect(self._remove_selected_file)
        self.confirm_roles_button.clicked.connect(self._confirm_suggested_roles)
        self.show_frame_button.clicked.connect(self._open_selected_file)
        self.file_table.itemSelectionChanged.connect(self._file_selection_changed)
        self.file_table.itemClicked.connect(self._file_row_clicked)
        # A click, not a selection change: the summary drives the files table,
        # and the files table drives the summary's highlight back.
        self.triage_summary_table.cellClicked.connect(self._summary_row_clicked)
        self.checklist_tree.currentRowChanged.connect(self._checklist_row_selected)
        self.anchor_table.itemSelectionChanged.connect(self._anchor_row_selected)
        # The two directions of the anchors panel's own sync: a bar names its
        # row, a row lights its bar.  Both ride signals that already exist.
        self.anchor_table.itemSelectionChanged.connect(self._highlight_selected_residual)
        self.residual_plot.scene().sigMouseClicked.connect(self._residual_plot_clicked)
        self.line_help_table.itemSelectionChanged.connect(self._expected_line_selected)
        self.next_step_button.clicked.connect(self._run_next_action)
        self.generate_tomls_button.clicked.connect(lambda: self._generate_tomls())
        self.regenerate_tomls_button.clicked.connect(self._regenerate_tomls)
        self.recompute_sphere_button.clicked.connect(self._start_sphere_comparison)
        self.choose_pattern_button.clicked.connect(self._choose_pattern_file)
        self.choose_previous_button.clicked.connect(self._choose_previous_pair)
        self.extract_pattern_button.clicked.connect(self._extract_pattern_from_sphere)
        self.save_snapshot_button.clicked.connect(self._save_snapshot)
        self.snapshot_id_edit.textChanged.connect(self.refresh_campaign)
        # ``textEdited`` fires for the operator's keystrokes and never for the
        # bench's own ``setText``, which is exactly the distinction between a
        # decision and a derivation.
        self.snapshot_id_edit.textEdited.connect(self._snapshot_id_hand_edited)
        # Picking the work on the left brings its own view with it: one lamp
        # choice, one view, the whole line workflow in one place.
        self.control_tabs.currentChanged.connect(self._control_tab_changed)

    #: Which view answers which piece of work.  Tabs the operator has not asked
    #: about are left exactly where they were.
    _VIEW_FOR_CONTROL_TAB = {"Files": "Exposure triage", "Lamp fit": "Lamp fit"}

    def _control_tab_changed(self, index: int) -> None:
        wanted = self._VIEW_FOR_CONTROL_TAB.get(self.control_tabs.tabText(index))
        if not wanted:
            return
        for position in range(self.view_tabs.count()):
            if self.view_tabs.tabText(position) == wanted:
                self.view_tabs.setCurrentIndex(position)
                return

    def _frame_choice_changed(self, index: int) -> None:
        """Fit the mean of the acquisition, or one of its single frames."""

        if self.session.frame is None:
            return
        choice = self.frame_combo.itemData(index)
        try:
            self.session.set_selected_frame(choice)
        except IndexError:
            return
        self.message_value.setText(
            f"Fitting {self.session.frame_choice_label()}; anchors were cleared "
            "because a centroid belongs to the spectrum it was measured on."
        )
        self.refresh()

    def _checklist_row_selected(self, row: int) -> None:
        if self.campaign is None or row < 0:
            return
        items = self.campaign.checklist(self.session)
        if row >= len(items):
            return
        item = items[row]
        text = item.detail
        if item.unblocked_by:
            text += f"\n\nWhat unblocks it: {item.unblocked_by}"
        if not item.blocking:
            text += "\n\nThis row is advice. It never blocks anything."
        self.explain(item.label, text)

    def _anchor_row_selected(self) -> None:
        row = self.anchor_table.currentRow()
        anchors = self.session.anchor_rows()
        if not (0 <= row < len(anchors)):
            return
        anchor = anchors[row]
        residual = {item.key: item for item in self.session.residuals}.get(anchor.key)
        residual_text = (
            "no residual yet — two anchors in different orders solve the transform"
            if residual is None
            else f"{residual.magnitude_px:.3f} px from the solved transform"
        )
        self.explain(
            f"{anchor.line.species} {anchor.line.wavelength_nm:.3f} nm, "
            f"order {anchor.line.order_idx}",
            f"Curated row expects detector pixel {anchor.line.center_pixel:.2f}; "
            f"the centroid on this frame sits at {anchor.fit.center_pixel:.2f} "
            f"(SNR {anchor.fit.snr:.1f}). Residual: {residual_text}. Raw window "
            f"peak {anchor.saturation.peak_value:.0f} counts, "
            f"{anchor.saturation.saturated_pixels} saturated pixel(s) — an "
            "anchor is only accepted when its own detector window is clear.",
        )

    def _bench_says(self, text: str) -> None:
        """Narrate an automatic follow-up without burying the operator's notice.

        Opening the assigned lamp signal and subtracting its background are
        things the bench does *because* roles were assigned, and both are
        already stated permanently beside the fit.  Landing them on the message
        line a moment later left the operator watching six roles change with
        the explanation already scrolled away (F21 item 1: one notice).  As
        soon as anything newer is on the line, the narration resumes.
        """

        if self._role_notice and self.message_value.text() == self._role_notice:
            return
        self.message_value.setText(text)

    def _say_roles(self, text: str) -> None:
        """Put one line about the roles on screen, and hold it there."""

        self._role_notice = text
        self.message_value.setText(text)

    def _save_says(self, text: str) -> None:
        """Answer a Save-tab press where the Save tab can be read.

        The bench has one message line and it lives on the Lamp fit tab, so a
        refusal raised by a button on the Save tab landed somewhere the
        operator was not looking — which is how "doesn't save, no errors" was a
        truthful report of a bench that had refused in full sentences.  Both
        lines carry it: the Save tab because that is where the press was, the
        message line because that is where the bench's running account is.  It
        also clears the held role notice, so this can never be suppressed as a
        follow-up narration by :meth:`_bench_says`.
        """

        self._role_notice = ""
        self.save_message_value.setText(text)
        self.message_value.setText(text)

    def _apply_unambiguous_suggestions(self) -> None:
        """Assign a whole drop's roles when its filenames leave nothing to ask.

        The bench used to hold every filename-derived role at arm's length
        until the operator confirmed it row by row, and the owner's folders
        always read the same way: "the suggestions are always like so".  A set
        of names that says exactly one thing is not a question, and asking it
        anyway is the confirmation dance this replaces (F21 item 1).

        Doubt keeps the old behaviour in full.  One ambiguous or clashing name
        anywhere in the drop and nothing is applied at all — every row keeps
        its SUGGESTED badge and its Confirm button — because a bench that
        guesses half a folder is worse than one that asks about all of it.

        The whole drop is judged at once, after the last queued file is read,
        so a six-file folder is one decision and one notice rather than six.
        """

        if self.campaign is None:
            return
        applied = self.campaign.apply_unanimous_suggestions(
            declined=self._declined_suggestions,
            saturation_level=self.session.saturation_level,
        )
        if not applied:
            return
        self.campaign.scope_alignment_to_lamp(self.session)
        # One line, after the refresh that may report opening a lamp signal:
        # what just happened to the roles is what the operator has to read.
        self.refresh()
        self._say_roles(
            f"Roles assigned from filenames: {len(applied)} file(s) — "
            "review them in the table and change any that are wrong."
        )

    def _confirm_suggested_roles(self) -> None:
        """Assign every filename suggestion the operator has not confirmed."""

        if self.campaign is None:
            return
        pending = self.campaign.unconfirmed_suggestions()
        if not pending:
            self.message_value.setText(
                "Every loaded file already carries an explicitly assigned role."
            )
            return
        confirmed = self.campaign.confirm_suggested_roles(
            saturation_level=self.session.saturation_level
        )
        self.campaign.scope_alignment_to_lamp(self.session)
        names = ", ".join(record.path.name for record in confirmed)
        outcome = (
            f"Assigned {len(confirmed)} suggested role(s): {names}. "
            "Change any of them in the Role column."
            if confirmed
            else "No suggestion could be confirmed; assign those roles by hand."
        )
        # The refresh may open the assigned lamp signal and report that; the
        # operator's own action is what they need to read, so it is set last
        # and held against the bench's own follow-up narration.
        self.refresh()
        self._say_roles(outcome)

    def _line_family_changed(self) -> None:
        """A hand-picked catalog is an override; the default follows the lamp."""

        self._family_override = True
        self._catalog_cache.clear()
        self._refresh_reference()
        self.refresh_plots()

    def _follow_assigned_lamp_catalog(self) -> None:
        """Fill the expected-lines panel from the assigned lamp on its own."""

        if self.campaign is None or self._family_override:
            return
        lamp = self.campaign.lamp_for_frame(self.session)
        index = self.line_family_combo.findText(lamp) if lamp else -1
        if index < 0 or index == self.line_family_combo.currentIndex():
            return
        self.line_family_combo.blockSignals(True)
        try:
            self.line_family_combo.setCurrentIndex(index)
        finally:
            self.line_family_combo.blockSignals(False)
        self._catalog_cache.clear()

    def _expected_line_selected(self) -> None:
        """Mark the selected expected line's stick on the spectrum above."""

        row = self.line_help_table.currentRow()
        if not (0 <= row < len(self._catalog_rows)):
            self.line_highlight.setVisible(False)
            return
        entry = self._catalog_rows[row]
        self.line_highlight.setValue(entry.detector_pixel)
        self.line_highlight.setVisible(True)
        self.explain(
            f"{entry.label} — expected at pixel {entry.detector_pixel:.1f}",
            f"{entry.wavelength_nm:.4f} nm, order {entry.order_idx}, from "
            f"{entry.source}. The pixel is where the wavelength table puts the "
            "line, moved by the solved alignment once there is one, so it says "
            "where to look, not where the line is: click the peak itself to fit "
            "its centroid and accept it as an anchor, and click an anchored "
            "stick again to take it back off.",
        )

    def _refresh_reference(self) -> None:
        """State which catalog anchors this fit, and when it is not the lamp's."""

        reference = self.session.reference
        if reference is None:
            self.reference_value.setText(
                "No lamp catalog is scoping the fit yet — assign a lamp role to "
                "the open file so anchors reference that lamp's own lines."
            )
            return
        text = reference.message
        warning = catalog_mismatch_warning(
            self.line_family_combo.currentText(), reference.lamp
        )
        if warning:
            text = f"{text}\nWARNING — {warning}"
        self.reference_value.setText(text)

    # ------------------------------------------------------------------
    # One folder at a time: opening the next campaign without a restart
    # ------------------------------------------------------------------

    def _pick_calibration_folder(self) -> None:
        """Ask for the next calibration folder, and confirm what it would cost.

        Owner, 2026-08-18: "opening a folder for a calibration would be nice.
        So I can start GUI, and go through many calibration folders. I probably
        would have to do that at NIFS."  A trip is a folder after a folder after
        a folder, and closing the whole bench between them was the restart this
        removes.
        """

        if self._campaign_thread is not None:
            self._save_says("A campaign task is running; open a folder when it finishes.")
            return
        # A read already in flight belongs to the folder being left, and it
        # would land in the fresh campaign a moment after it was emptied — the
        # one file of the old folder that followed the bench to the new one.
        # Checked here rather than by greying the button, because the button's
        # enabled state is only as fresh as the last refresh and a queue that
        # applies no roles ends without one.
        if self._load_thread is not None or self._queue:
            self._save_says(
                "The bench is still reading this folder; open the next one when "
                "it finishes."
            )
            return
        start = self.calibration_folder or self.last_folder
        chosen = choose_calibration_folder(self, start)
        if not chosen:
            return
        folder = Path(chosen)
        if not folder.is_dir():
            self._save_says(f"That is not a folder: {folder}")
            return
        warning = self._unsaved_session_warning()
        if warning:
            answer = QtWidgets.QMessageBox.question(
                self,
                "Open another calibration folder",
                warning,
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No,
            )
            if answer != QtWidgets.QMessageBox.Yes:
                self._save_says(f"{folder.name} was not opened; the session is untouched.")
                return
        self.open_calibration_folder(folder)

    def _unsaved_session_warning(self) -> str:
        """The one sentence a session with unsaved work is owed, else ``""``.

        Unsaved means what the operator would mind losing: an alignment solved,
        or files loaded and roles assigned, with no snapshot saved to show for
        it.  The campaign's own ``saved_snapshot`` answers the second half and
        answers it honestly — it is cleared by ``_invalidate_outputs`` the
        moment anything the snapshot was built from changes, so a save followed
        by more work asks again, and a save followed by nothing does not ask at
        all.  An empty bench is never asked anything.
        """

        loaded = bool(self._file_rows) or bool(self._queue)
        solved = bool(self.session.anchors) or self.session.frame is not None
        if not (loaded or solved):
            return ""
        if self.campaign is not None and self.campaign.saved_snapshot is not None:
            return ""
        identity = (
            self.snapshot_id_edit.text().strip()
            or self.initial_snapshot_id
            or "this session"
        )
        return (
            "Opening a new folder clears the current session — the alignment "
            f"for {identity} was not saved."
        )

    def open_calibration_folder(self, folder: str | Path) -> list[Path]:
        """Stand the whole bench at *folder*, as launching there would have.

        This is the launch, minus the process: everything ``main`` derives from
        its folder argument is derived again here from this one — both output
        roots through :func:`default_bench_roots`, the acquisition-dated
        identity through :func:`acquisition_date_from_name`, the file dialogs'
        starting point, and the watch folder if one is running.  The folder's
        SIFs are then handed to :meth:`add_paths`, the same door a drop comes
        through, so triage runs and the roles suggest and auto-apply exactly as
        they do for a drop.

        What the old session held is dropped through the machinery that already
        knows how — :meth:`CalibrationCampaignSession.forget_all_files` and its
        ``_invalidate_outputs``, :meth:`CalibrationBenchSession.forget_frame`
        and its ``clear_anchors`` — rather than through a second reset here that
        would drift from them.

        The chosen pattern and wavelength tables are deliberately *not* reset.
        The operator picked those, in this window, against their own judgement
        of the detector; a folder does not overrule that.  The band guard reads
        the new folder's sphere against the chosen pattern and says so when the
        two disagree, which is the honest way for a stale choice to surface.
        """

        target = absolute_root(Path(folder))
        self.calibration_folder = target
        self.last_folder = target
        self.output_root, self.config_root = default_bench_roots(target)
        self._describe_output_roots()
        self._rewatch(target)
        self._forget_session()
        self._rename_snapshot_identity(target)
        found = self.add_paths([target])
        self.refresh()
        tail = (
            f"{len(found)} SIF file(s) queued for triage."
            if found
            else "It holds no SIF files yet — drop them in and they will load."
        )
        self._save_says(
            f"Opened {target}. {tail} Snapshots and settings will be written "
            f"under {self.output_root}."
        )
        return found

    def _rewatch(self, folder: Path) -> None:
        """Re-aim a running folder watch at the folder now being calibrated."""

        if self.watcher is None:
            return
        self.watcher = StableSifWatcher(
            folder,
            required_unchanged_polls=self.watcher.required_unchanged_polls,
            minimum_age_s=self.watcher.minimum_age_ns / 1_000_000_000,
        )

    def _forget_session(self) -> None:
        """Empty the bench of one folder's work before the next one lands."""

        self._queue.clear()
        self._file_rows.clear()
        self.file_table.setRowCount(0)
        self._declined_suggestions.clear()
        self._explicitly_opened.clear()
        self._arrivals_pending = False
        self._landed_on = None
        self._auto_following = False
        self._family_override = False
        self._role_notice = ""
        self._refused_identity = ""
        self.regenerate_tomls_button.setVisible(False)
        self._saved_snapshot_root = None
        # Drawn from the old frames and the old solve; a new folder redraws it.
        self._catalog_cache.clear()
        self._catalog_rows = ()
        self._drawn_correction = None
        if self.campaign is not None:
            self.campaign.forget_all_files()
        self.session.forget_frame()

    def _rename_snapshot_identity(self, folder: Path) -> None:
        """Re-date the prefilled identity from the folder just opened.

        The same rule the launch argument obeys, applied to the same kind of
        fact: a calibration is dated by the day its frames were taken, the
        folder's own name is the first thing that says so, and the frames' SIF
        headers are the second — which still fill this in through
        ``_adopt_acquisition_date`` when the name says nothing, because opening
        a folder puts the identity back to being the bench's guess rather than
        anybody's decision.
        """

        self.snapshot_date = acquisition_date_from_name(folder)
        self._snapshot_id_decided = False
        self._snapshot_date_source = (
            "the name of the folder you opened" if self.snapshot_date else ""
        )
        self.initial_snapshot_id = snapshot_identity(
            self.snapshot_date or date.today(),
            self.detector_edit.text().strip() or self.initial_detector,
        )
        # ``textEdited`` does not fire for a programmatic change, so writing the
        # derived identity here cannot mark it as somebody's decision.
        self.snapshot_id_edit.setText(self.initial_snapshot_id)
        self._describe_snapshot_identity()

    def _describe_output_roots(self) -> None:
        """Say where this folder's snapshots and settings bundles will land."""

        self.snapshot_root_value.setText(f"Snapshots: {self.output_root}")
        self._explainable(
            self.snapshot_root_value,
            "Where saved snapshots are written",
            "Every snapshot this bench saves is created inside this folder, "
            "one subfolder per snapshot identity. Unless you passed "
            "--output-root, it is the calibrations folder inside the "
            "calibration folder the bench is open at — the computed "
            "calibration sits beside the lamp frames it was computed from, and "
            "that folder holds it all. Everything the bench generates lives "
            "under this one folder, the settings bundles included, in a "
            f"configs subfolder. In full: {self.output_root}",
            hint=str(self.output_root),
        )
        self.config_root_value.setText(f"Configs: {self.config_root}")
        self._explainable(
            self.config_root_value,
            "Where the generated settings files are written",
            "The commented campaign, alignment and export files are written "
            "inside this folder, one subfolder per snapshot identity. Unless "
            "you passed --config-root, it is the configs folder inside the "
            "calibrations folder above — one tidy folder holds everything the "
            "bench generates, beside the lamp frames it was derived from. In "
            f"full: {self.config_root}",
            hint=str(self.config_root),
        )

    # ------------------------------------------------------------------
    # Manual input: drag and drop, and a plain file dialog
    # ------------------------------------------------------------------

    def dragEnterEvent(self, event) -> None:  # noqa: N802 - Qt naming
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802 - Qt naming
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:  # noqa: N802 - Qt naming
        if not event.mimeData().hasUrls():
            event.ignore()
            return
        event.setDropAction(QtCore.Qt.CopyAction)
        event.accept()
        self.add_paths(
            [url.toLocalFile() for url in event.mimeData().urls() if url.toLocalFile()]
        )

    def _pick_files(self) -> None:
        paths, _filter = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Add SIF files",
            str(self.last_folder),
            "Andor SIF (*.sif *.SIF);;All files (*)",
        )
        self.add_paths(paths)

    def add_paths(self, paths: Sequence[str | Path]) -> list[Path]:
        """Queue any dropped or picked files, whatever they are called."""

        accepted: list[Path] = []
        rejected: list[str] = []
        for item in paths:
            source = Path(item)
            if source.is_dir():
                accepted.extend(sorted(source.glob("*.[sS][iI][fF]")))
                continue
            if not source.is_file():
                rejected.append(f"{source.name} is not a file")
                continue
            accepted.append(source)
        for source in accepted:
            self.last_folder = source.parent
            if source not in self._queue:
                self._queue.append(source)
        # Only files that ARRIVE are judged as a set: a drop, an Add files, a
        # new file in the watch folder.  Re-reading a file the bench already
        # holds is a load but not an arrival — opening one for fitting comes
        # through here too, both when the operator asks and when the bench
        # follows the assigned lamp signal itself — and judging the whole folder
        # again when one of those lands would apply suggestions to a folder the
        # operator is halfway through assigning by hand, on no event he caused.
        if any(source not in self._file_rows for source in accepted):
            self._arrivals_pending = True
        if rejected:
            self.message_value.setText("; ".join(rejected))
        elif accepted:
            self.message_value.setText(
                f"Queued {len(accepted)} file(s) for reading and triage."
            )
        self._start_next_load()
        return accepted

    def _start_next_load(self) -> None:
        if self._load_thread is not None or not self._queue:
            return
        self.load_path(self._queue.pop(0))

    # ------------------------------------------------------------------
    # One manual role control per loaded file
    # ------------------------------------------------------------------

    def _selected_file(self) -> Path | None:
        row = self.file_table.currentRow()
        if 0 <= row < len(self._file_rows):
            return self._file_rows[row]
        return None

    def _file_selection_changed(self) -> None:
        selected = self._selected_file()
        # The buttons follow the selection immediately: a control that ignores
        # a click because some other refresh has not run yet is a dead surface.
        self._refresh_file_buttons()
        if selected is not None:
            self._show_triage(selected)

    def _refresh_file_buttons(self) -> None:
        """Enable the per-file actions for whatever row is selected now."""

        selected = self._selected_file()
        self.remove_file_button.setEnabled(selected is not None)
        self.show_frame_button.setEnabled(
            selected is not None and self.loader is not None and self._load_thread is None
        )

    def _open_selected_file(self) -> None:
        selected = self._selected_file()
        if selected is None:
            return
        # An explicit choice is never argued with; the fit view only says what
        # it is showing and warns when that is a background.
        self._explicitly_opened.add(selected)
        if self.session.frame is not None and self.session.frame.path == selected:
            self.message_value.setText(f"{selected.name} is already open for fitting.")
            return
        self.add_paths([selected])

    # ------------------------------------------------------------------
    # The fit tab fits the assigned lamp signal, and says which file that is
    # ------------------------------------------------------------------

    def _assigned_lamp_files(self) -> tuple[Path | None, Path | None]:
        """The primary lamp's (signal, background) paths, as far as assigned."""

        if self.campaign is None:
            return None, None
        return self.campaign.lamp_pair(self.campaign.lamp_for_frame(self.session))

    def _open_frame_role(self) -> MeasurementRole | None:
        frame = self.session.frame
        if frame is None or self.campaign is None:
            return None
        record = self.campaign.measurements.get(Path(frame.path))
        return record.role if record is not None else None

    def _follow_assigned_lamp_signal(self) -> None:
        """Land the fit on the assigned lamp signal rather than whatever is open.

        The frame that happens to be open after loading a folder is the last
        file read, which on the owner's folder was a background: a lineless
        hump the fit view showed without ever naming it.
        """

        if self.campaign is None or self.loader is None or self._auto_following:
            return
        # Loading refreshes the surface, and the surface asks this question
        # again: without the guard the two call each other forever.
        if self._load_thread is not None or self._queue:
            return
        if self.session.file_state is FileLoadState.LOADING:
            return
        signal, _background = self._assigned_lamp_files()
        if signal is None or not signal.is_file():
            return
        open_path = getattr(self.session.frame, "path", None)
        if open_path is not None and Path(open_path) == signal:
            return
        if open_path is not None and Path(open_path) in self._explicitly_opened:
            return
        if self._open_frame_role() is MeasurementRole.LAMP:
            return
        self._auto_following = True
        try:
            self._bench_says(
                f"Opened the assigned lamp signal {signal.name} for fitting."
            )
            self.add_paths([signal])
        finally:
            self._auto_following = False

    def _pair_lamp_background(self) -> None:
        """Subtract the assigned lamp background once the pair is complete."""

        if self.campaign is None or self.loader is None:
            return
        if self._background_thread is not None:
            return
        signal, background = self._assigned_lamp_files()
        open_path = getattr(self.session.frame, "path", None)
        if (
            open_path is None
            or signal is None
            or background is None
            or Path(open_path) != signal
            or not background.is_file()
        ):
            if self.session.background_path is not None and (
                background is None or open_path is None or Path(open_path) != signal
            ):
                self.session.use_background_frame(None)
            return
        if self.session.background_path == background:
            return
        thread = FrameLoadThread(background, self.loader, self)
        thread.loaded.connect(self._background_loaded)
        thread.failed.connect(self._background_failed)
        thread.finished.connect(self._background_finished)
        self._background_thread = thread
        thread.start()

    @QtCore.pyqtSlot(object)
    def _background_loaded(self, frame: BenchFrame) -> None:
        try:
            self.session.use_background_frame(frame)
        except ValueError as exc:
            self.message_value.setText(f"Lamp background was not subtracted: {exc}")
            return
        self._bench_says(
            f"Fitting {frame.path.name} subtracted from the lamp signal, as "
            "echelle-align does."
        )
        self.refresh()

    @QtCore.pyqtSlot(str, str)
    def _background_failed(self, path: str, reason: str) -> None:
        self.message_value.setText(
            f"The assigned lamp background {Path(path).name} could not be read "
            f"({reason}); the fit shows the raw signal."
        )

    @QtCore.pyqtSlot()
    def _background_finished(self) -> None:
        if self._background_thread is not None:
            self._background_thread.deleteLater()
        self._background_thread = None

    def _land_on_the_lamps_own_orders(self) -> None:
        """Open a lamp on an order that actually carries its lines.

        Landing on order 0 with an empty expected-line panel is the empty room
        the operator had to discover and feed, moved one tab along. This lands
        once per frame-and-lamp; every later order change is the operator's.
        """

        frame = self.session.frame
        if frame is None:
            return
        reference = self.session.reference
        key = (Path(frame.path), "" if reference is None else reference.lamp)
        if key == self._landed_on:
            return
        self._landed_on = key
        order = self.session.best_reference_order()
        if order is None or order == self.order_spin.value():
            return
        self.order_spin.setValue(order)

    def _refresh_fit_target(self) -> None:
        """State which file the fit is measuring, and warn when it is wrong."""

        frame = self.session.frame
        if frame is None:
            self.fit_file_value.setText("no file open for fitting")
            self.fit_warning_value.setVisible(False)
            return
        role = self._open_frame_role()
        lamp = "" if self.campaign is None else self.campaign.lamp_for_frame(self.session)
        described = "unassigned file" if role is None else role.value
        if role in _LAMP_ROLES and lamp:
            described = f"{lamp} {described}"
        self.fit_file_value.setText(
            f"Fitting {frame.path.name} · {described} · "
            f"{self.session.frame_choice_label()} · {self.session.background_label()}"
        )
        warning = ""
        if role in {MeasurementRole.LAMP_BACKGROUND, MeasurementRole.SPHERE_BACKGROUND}:
            warning = (
                f"{frame.path.name} is a background frame — no lines are "
                "expected in it. Open the lamp signal from the Files tab to fit "
                "lines."
            )
        elif role is MeasurementRole.SPHERE:
            warning = (
                f"{frame.path.name} is the integrating-sphere signal — a smooth "
                "continuum, not a line spectrum. Open the lamp signal to fit lines."
            )
        self.fit_warning_value.setText(warning)
        self.fit_warning_value.setVisible(bool(warning))

    def _remove_selected_file(self) -> None:
        selected = self._selected_file()
        if selected is None:
            return
        if self.campaign is not None:
            self.campaign.forget_file(selected)
        self._declined_suggestions.discard(selected)
        self._file_rows.remove(selected)
        self.file_table.removeRow(self.file_table.currentRow())
        self.message_value.setText(f"Removed {selected.name} from the bench.")
        self.refresh()

    def _add_file_row(self, path: Path) -> None:
        row = self.file_table.rowCount()
        self.file_table.insertRow(row)
        self._file_rows.append(path)
        self.file_table.setItem(row, 0, QtWidgets.QTableWidgetItem(path.name))
        role_combo = QtWidgets.QComboBox()
        for label, role in _ROLE_CHOICES:
            role_combo.addItem(label, role)
        lamp_combo = QtWidgets.QComboBox()
        lamp_combo.setEditable(True)
        lamp_combo.addItems(list(KNOWN_LAMP_NAMES))
        lamp_combo.setCurrentText("")
        self.file_table.setCellWidget(row, 1, role_combo)
        self.file_table.setCellWidget(row, 2, lamp_combo)
        role_combo.currentIndexChanged.connect(
            lambda _index, source=path: self._role_changed(source)
        )
        # Picking the pre-filled entry again is a deliberate confirmation, and
        # it emits no index change, so the operator's click is heard here.
        role_combo.activated.connect(
            lambda _index, source=path: self._role_changed(source)
        )
        lamp_combo.currentTextChanged.connect(
            lambda _text, source=path: self._role_changed(source)
        )
        self._prefill_role(path, role_combo, lamp_combo)

    def _prefill_role(self, path: Path, role_combo, lamp_combo) -> None:
        """Offer the filename's guess without assigning anything."""

        if self.campaign is None:
            return
        suggestion = self.campaign.observed.get(path)
        if suggestion is None or not suggestion.is_unambiguous:
            return
        role_combo.blockSignals(True)
        lamp_combo.blockSignals(True)
        try:
            index = role_combo.findData(suggestion.roles[0])
            if index >= 0:
                role_combo.setCurrentIndex(index)
            if suggestion.lamp_name:
                lamp_combo.setCurrentText(suggestion.lamp_name)
        finally:
            role_combo.blockSignals(False)
            lamp_combo.blockSignals(False)
        lamp_combo.setEnabled(suggestion.roles[0] in _LAMP_ROLES)

    def _role_changed(self, path: Path) -> None:
        if self.campaign is None or path not in self._file_rows:
            return
        row = self._file_rows.index(path)
        role_combo = self.file_table.cellWidget(row, 1)
        lamp_combo = self.file_table.cellWidget(row, 2)
        role = role_combo.currentData()
        lamp_combo.setEnabled(role in _LAMP_ROLES)
        if role is None:
            # Taking a role off is a decision, and the next drop must not put
            # the filename's guess back on top of it.
            self._declined_suggestions.add(path)
            if self.campaign.remove_classification(path):
                self.campaign.scope_alignment_to_lamp(self.session)
                self.message_value.setText(f"{path.name} is unassigned again.")
                self.refresh()
            return
        self._declined_suggestions.discard(path)
        try:
            record = self.campaign.classify_file(
                path,
                role,
                lamp_family=lamp_combo.currentText(),
                saturation_level=self.session.saturation_level,
            )
        except (FileNotFoundError, ValueError) as exc:
            self.message_value.setText(f"{path.name} keeps its previous role: {exc}")
            return
        label = record.role.value
        if record.lamp_family:
            label = f"{record.lamp_family} {label}"
            if self.line_family_combo.findText(record.lamp_family) >= 0:
                self.line_family_combo.setCurrentText(record.lamp_family)
        self.campaign.scope_alignment_to_lamp(self.session)
        self.message_value.setText(
            f"{record.path.name} is now the {label}. "
            "Dependent comparison/configuration state was reset."
        )
        self.refresh()

    def _start_campaign_task(self, operation) -> None:
        if self._campaign_thread is not None:
            return
        thread = CampaignTaskThread(operation, self)
        thread.completed.connect(self._campaign_task_completed)
        thread.failed.connect(self._campaign_task_failed)
        thread.finished.connect(self._campaign_task_finished)
        self._campaign_thread = thread
        self.refresh_campaign()
        thread.start()

    def _sphere_pair_assigned(self) -> bool:
        """Whether both halves of the sphere pair carry their role."""

        if self.campaign is None:
            return False
        roles = {record.role for record in self.campaign.measurements.values()}
        return {MeasurementRole.SPHERE, MeasurementRole.SPHERE_BACKGROUND} <= roles

    def _start_sphere_comparison(self) -> None:
        if self.campaign is None:
            return
        self.comparison_value.setText("COMPUTING — using the established absolute engine…")
        self._start_campaign_task(self.campaign.compute_sphere_comparison)

    def _pattern_root(self) -> Path:
        """Where patterns this bench extracts are written: one tidy folder."""

        return self.output_root / PATTERN_ROOT_NAME

    def _extracted_pattern_path(self) -> Path:
        """A name that says what it is and when, and never clobbers an earlier one.

        A pattern a snapshot was computed on is not scratch: an extraction made
        an hour later is a second file, not a silent replacement of the first.
        """

        # Dated by the calibration, not by the day somebody pressed the button:
        # the identity in the Save tab is the acquisition's own date whenever
        # anything could say what it was.
        dated = (
            acquisition_date_from_name(self.snapshot_id_edit.text().strip())
            or self.snapshot_date
            or date.today()
        )
        stamp = dated.strftime("%Y%m%d")
        root = self._pattern_root()
        candidate = root / f"pattern_extracted_{stamp}.txt"
        attempt = 2
        while candidate.exists():
            candidate = root / f"pattern_extracted_{stamp}_{attempt}.txt"
            attempt += 1
        return candidate

    def _extraction_blocker(self) -> str:
        """Why the sphere cannot be extracted from, in words, or empty."""

        if self.campaign is None:
            return "campaign memory was not configured for this bench"
        if self.loader is None:
            return "no SIF reader is configured for this bench"
        sphere_path, _background = self.campaign.sphere_pair_paths()
        if sphere_path is None:
            return (
                "no file carries the sphere signal role — assign it in the Files "
                "tab and the pattern can be fitted to the light it holds"
            )
        return ""

    def _extract_pattern_from_sphere(self) -> None:
        """Fit a pattern to the assigned sphere and stand the bench on it."""

        blocker = self._extraction_blocker()
        if blocker:
            self._save_says(f"No pattern was extracted: {blocker}.")
            return
        assert self.campaign is not None
        sphere_path, background_path = self.campaign.sphere_pair_paths()
        assert sphere_path is not None
        prior = np.rint(self.session.pattern).astype(int)
        output_path = self._extracted_pattern_path()
        loader = self.loader
        self._save_says(
            f"Extracting a pattern from {sphere_path.name} with "
            f"{self.campaign.pattern_source.name} as the prior…"
        )
        self._start_campaign_task(
            lambda: extract_pattern_beside_campaign(
                loader=loader,
                prior=prior,
                sphere_path=sphere_path,
                background_path=background_path,
                output_path=output_path,
            )
        )

    def _choose_pattern_file(self) -> None:
        """Name a different pattern file, on the same terms as an extracted one."""

        if self.campaign is None or self._campaign_thread is not None:
            return
        start = self.campaign.pattern_source.parent
        chosen, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Choose the order-pattern table this bench reads",
            str(start if start.is_dir() else self.last_folder),
            "Pattern tables (*.txt);;All files (*)",
        )
        if not chosen:
            return
        chosen_path = Path(chosen)
        if chosen_path == self.campaign.pattern_source:
            self._save_says(
                f"The bench already wears {chosen_path}; nothing was changed."
            )
            return
        sphere_path, _background = self.campaign.sphere_pair_paths()
        loader = self.loader
        self._save_says(f"Standing the bench on {chosen_path}…")
        self._start_campaign_task(
            lambda: read_pattern_beside_campaign(
                loader=loader, chosen_path=chosen_path, sphere_path=sphere_path
            )
        )

    def _choose_previous_pair(self) -> None:
        """Name the previous campaign's sphere pair from inside the bench.

        One dialog where one is enough: the signal is chosen, and its
        background is read off the file naming beside it and offered for
        confirmation.  The second dialog appears only when the folder does not
        answer the question or the operator says that is not the file — which
        is the difference between a suggestion and a decision, exactly as the
        role suggestions work everywhere else on this bench.
        """

        if self.campaign is None or self._campaign_thread is not None:
            return
        start = self.campaign.previous_sphere
        folder = start.parent if start is not None and start.parent.is_dir() else self.last_folder
        chosen, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Choose the previous campaign's integrating-sphere SIGNAL frame",
            str(folder),
            "Andor SIF frames (*.sif);;All files (*)",
        )
        if not chosen:
            return
        signal = Path(chosen)
        background = self._previous_background_for(signal)
        if background is None:
            return
        self._adopt_previous_pair(signal, background)

    def _previous_background_for(self, signal: Path) -> Path | None:
        """The background for a chosen previous signal: offered, then asked for."""

        assert self.campaign is not None
        suggested = background_sibling(signal)
        if suggested is not None and suggested != signal:
            answer = QtWidgets.QMessageBox.question(
                self,
                "The background beside it",
                f"{suggested.name} looks like the background shot with "
                f"{signal.name}.\n\nUse it as the previous background?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.Yes,
            )
            if answer == QtWidgets.QMessageBox.Yes:
                return suggested
        chosen, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Choose the previous campaign's integrating-sphere BACKGROUND frame",
            str(signal.parent),
            "Andor SIF frames (*.sif);;All files (*)",
        )
        return Path(chosen) if chosen else None

    def _adopt_previous_pair(self, signal: Path, background: Path) -> None:
        """Swap the pair, drop what was measured against the old one, and say so."""

        assert self.campaign is not None
        try:
            self.campaign.adopt_previous_pair(signal, background)
        except (FileNotFoundError, ValueError) as error:
            self._save_says(str(error))
            return
        self.refresh()
        self._save_says(
            f"Comparing against {self.campaign.previous_pair_name}. The previous "
            "comparison was measured against another pair and has been dropped; "
            "press Recompute."
        )

    def _refresh_previous_pair(self, busy: bool) -> None:
        """Say which previous pair the next comparison would run against."""

        if self.campaign is None:
            self.previous_pair_value.setText("Previous pair: —")
            self.choose_previous_button.setEnabled(False)
            return
        name = self.campaign.previous_pair_name
        self.previous_pair_value.setText(
            f"Previous pair: {name}" if name else "Previous pair: none chosen"
        )
        self.previous_pair_value.setToolTip(
            f"{self.campaign.previous_sphere}\n{self.campaign.previous_sphere_background}"
            if name
            else "No previous campaign pair has been chosen for the comparison."
        )
        self.choose_previous_button.setEnabled(not busy)

    def _adopt_pattern_rebase(self, result: PatternRebase) -> None:
        """Swap the pattern under the whole bench, and say what that cost.

        Everything derived from the old pattern is invalidated through the paths
        that already exist — the campaign's own ``_invalidate_outputs`` by way of
        :meth:`CalibrationCampaignSession.adopt_pattern`, and the session's
        ``clear_anchors`` by way of :meth:`CalibrationBenchSession.adopt_pattern`
        — rather than through a second reset that would drift from them.
        """

        if self.campaign is None:
            return
        # Read before the swap: a pattern with a different order count releases
        # the open frame, and it is exactly that frame which must be extracted
        # again on the new geometry.
        frame_path = None if self.session.frame is None else self.session.frame.path
        try:
            rows = np.loadtxt(result.path, dtype=int)
            reading = self.campaign.adopt_pattern(
                result.path,
                sphere_image=result.sphere_image,
                saturation_level=self.session.saturation_level,
            )
            cleared = self.session.adopt_pattern(rows)
        except (OSError, ValueError) as exc:
            self._save_says(f"That pattern was not adopted: {exc}")
            self.refresh()
            return
        adopt = getattr(self.loader, "adopt_pattern", None)
        if adopt is not None:
            adopt(rows)
        # The traces are cached on the pattern's SHAPE, which a re-extracted
        # pattern usually shares with the one it replaces.  Without this the new
        # geometry would be drawn as the old one until the order count changed.
        self._pattern_key = None
        self._catalog_cache.clear()

        sentences = []
        if result.extracted:
            offset = result.median_offset_rows
            moved = "" if offset is None else f", moving the traces {offset:+.1f} row(s)"
            sentences.append(
                f"Extracted a pattern from {result.sphere_path.name}{moved}, wrote "
                f"it to {result.path}, and stood the bench on it."
            )
        else:
            sentences.append(f"The bench now reads its orders off {result.path}.")
        if result.caveat:
            sentences.append(f"Note: {result.caveat}.")
        if cleared:
            sentences.append(
                f"{cleared} anchor(s) were cleared — they were fitted on the old "
                "geometry, and auto-anchor re-runs them in seconds."
            )
        if reading is not None and reading.measured:
            sentences.append(
                f"Against the new pattern, {reading.summary()}."
                if self.campaign.pattern_band_warning()
                else f"The band guard now reads quiet: {reading.summary()}."
            )
        sentences.append(
            "Factors, settings bundle and save state were reset: they were "
            "computed over the old pattern."
        )
        # Re-extracting the open frame is the last thing started, so the
        # narration above is on screen while it runs.
        if frame_path is not None:
            sentences.append(f"Re-extracting {frame_path.name} on the new pattern.")
        self.refresh()
        self._save_says(" ".join(sentences))
        if frame_path is not None:
            # Through the queue rather than straight at the reader: a read
            # already in flight would swallow this one, and the frame would go
            # on being fitted off the old geometry with nothing saying so.
            self._queue.append(frame_path)
            self._start_next_load()

    @QtCore.pyqtSlot(object)
    def _campaign_task_completed(self, result) -> None:
        if isinstance(result, PatternRebase):
            # The file was written off the event loop; the swap itself touches
            # the session, the reader and the plots, so it happens here.
            self._adopt_pattern_rebase(result)
            return
        state = getattr(result, "state", None)
        if state is ComparisonState.READY:
            self.message_value.setText(
                "Sphere factors computed, but the previous pair is a copy of "
                "this campaign's own sphere — the ratio proves nothing."
                if getattr(result, "self_comparison", False)
                else "Sphere factors computed and compared."
            )
        elif state is ComparisonState.INSUFFICIENT_DATA:
            self.message_value.setText(
                "Candidate factors computed; previous comparison is insufficient data."
            )
        elif hasattr(result, "snapshot_id"):
            correction = getattr(self.campaign, "wavelength_correction", None)
            detail = "" if correction is None else f" Saved wavelength.txt: {correction.reason}."
            pattern_correction = getattr(self.campaign, "pattern_correction", None)
            if pattern_correction is not None:
                detail += f" Saved pattern.txt: {pattern_correction.reason}."
            # The warning said before the press must survive the press: this
            # line is the last thing the save panel shows, and a snapshot saved
            # over a pattern the light does not fit has to say so here too.
            band_warning = (
                "" if self.campaign is None else self.campaign.pattern_band_warning()
            )
            if band_warning:
                detail += (
                    f" PATTERN DOES NOT FIT THIS SPHERE — {band_warning}; the "
                    "measured offset is in the manifest."
                )
            # The one place the operator learns where the campaign went.  It
            # says the whole path, because "calibrations" said nothing about
            # which machine's calibrations folder it meant.
            saved_root = absolute_root(
                getattr(result, "root", self.output_root / result.snapshot_id)
            )
            self._saved_snapshot_root = saved_root
            self.open_snapshot_button.setVisible(True)
            self.open_snapshot_button.setToolTip(f"Open {saved_root} in the file manager")
            self._save_says(
                f"Snapshot {result.snapshot_id} saved and validated through "
                f"Packet 0 in {saved_root}.{detail}"
            )
        self.refresh()

    def _open_saved_snapshot_folder(self) -> None:
        """Show the folder the last save wrote, in the system file manager."""

        if self._saved_snapshot_root is None:
            return
        if not open_containing_folder(self._saved_snapshot_root):
            self.message_value.setText(
                f"The file manager would not open {self._saved_snapshot_root}; "
                "the snapshot itself is saved and untouched."
            )

    @QtCore.pyqtSlot(str)
    def _campaign_task_failed(self, reason: str) -> None:
        # Every campaign task is started from a button on the Save tab, so its
        # refusal is said there as well as on the message line.
        self._save_says(f"Campaign action failed safely: {reason}")
        self.refresh()

    @QtCore.pyqtSlot()
    def _campaign_task_finished(self) -> None:
        if self._campaign_thread is not None:
            self._campaign_thread.deleteLater()
        self._campaign_thread = None
        self.refresh_campaign()

    def _generate_tomls(self, *, overwrite: bool = False) -> None:
        if self.campaign is None:
            return
        snapshot_id = self.snapshot_id_edit.text().strip()
        try:
            paths = self.campaign.write_tomls(
                self.config_root,
                snapshot_id,
                self.session,
                overwrite=overwrite,
                # Where Save will put the snapshot, so the generated export
                # configuration can name the sphere frames from there exactly
                # as the snapshot's own binder does.
                snapshot_root=self.output_root / snapshot_id,
            )
        except (OSError, SnapshotError, ValueError) as exc:
            outcome = f"Alignment settings were not saved: {exc}"
            # An existing bundle is a choice to offer, not a dead end: the
            # button that can get past it appears beside the one that refused.
            # The offer is recorded against the identity that was refused, so
            # editing the id above cannot carry it to a different bundle.
            refused = "already exists" in str(exc)
            self._refused_identity = snapshot_id if refused else ""
            self.regenerate_tomls_button.setVisible(refused)
        else:
            self._refused_identity = ""
            self.regenerate_tomls_button.setVisible(False)
            outcome = (
                "Saved the alignment settings, the campaign, and the export "
                f"configuration as commented files you can edit: {paths['campaign'].parent}."
                if not overwrite
                else "Rewrote the alignment settings, the campaign, and the export "
                f"configuration, replacing what was in {paths['campaign'].parent}."
            )
            self.toml_preview.setPlainText(
                paths["campaign"].read_text(encoding="utf-8")
            )
        # The refresh below re-reads the roles and may narrate what the bench
        # did about them; the answer to the press is what the operator needs to
        # read, so it is put on screen last and cannot be narrated over.
        self.refresh_campaign()
        self._save_says(outcome)

    def _regenerate_tomls(self) -> None:
        """Rewrite settings files that already exist, on a deliberate press."""

        if self._campaign_thread is not None:
            return
        if self._refused_identity != self.snapshot_id_edit.text().strip():
            # Overwriting is permission for one identity, never for whatever the
            # id field says at the moment of the press.
            self._refused_identity = ""
            self.regenerate_tomls_button.setVisible(False)
            self._save_says(
                "The snapshot identity changed since the bench refused to "
                "overwrite; press Save alignment settings again for this one."
            )
            return
        self._generate_tomls(overwrite=True)

    def _save_snapshot(self) -> None:
        if self.campaign is None:
            return
        # Said before the save, not instead of it.  An operator who knows why
        # this pattern is the right one for this sphere may go ahead; what the
        # bench owes is that nobody saves this without having been told, and
        # that the number goes into the manifest whichever way it goes.
        warning = self.campaign.pattern_band_warning()
        if warning:
            self._save_says(
                f"PATTERN DOES NOT FIT THIS SPHERE — {warning}. Saving anyway: "
                "the measured offset is recorded in the snapshot manifest."
            )
        snapshot_id = self.snapshot_id_edit.text().strip()
        detector = self.detector_edit.text().strip()
        base_snapshot = self.base_snapshot_edit.text().strip() or None
        notes = self.notes_edit.text().strip()
        self._start_campaign_task(
            lambda: self.campaign.save_snapshot(
                self.output_root,
                snapshot_id=snapshot_id,
                detector=detector,
                alignment=self.session,
                notes=notes,
                base_snapshot=base_snapshot,
                validity=default_validity(snapshot_id, self.valid_from),
            )
        )

    def poll_watch_folder(self) -> None:
        """Poll once; loading begins only for a newly emitted stable file."""

        if self.watcher is None or self.loader is None or self._load_thread is not None:
            return
        result = self.watcher.poll()
        if result.path is not None:
            self.file_value.setText(result.path.name)
        if result.state is StableFileState.FAILED:
            self.message_value.setText(f"Watch failed: {result.reason}")
        elif result.state is StableFileState.CHANGING:
            self.message_value.setText(
                f"Waiting for {result.path.name} to stop changing "
                f"({result.unchanged_polls}/{self.watcher.required_unchanged_polls})."
            )
        elif result.ready_path is not None:
            self.add_paths([result.ready_path])

    def load_path(self, path: str | Path) -> None:
        """Start an asynchronous read of one SIF, whatever it is called."""

        if self._load_thread is not None:
            return
        if self.loader is None:
            self.message_value.setText("No SIF reader is configured for this bench.")
            return
        source = Path(path)
        self.session.begin_file_load(source)
        self.refresh()
        thread = FrameLoadThread(source, self.loader, self)
        thread.loaded.connect(self._frame_loaded)
        thread.failed.connect(self._frame_failed)
        thread.finished.connect(self._load_finished)
        self._load_thread = thread
        thread.start()

    @QtCore.pyqtSlot(object)
    def _frame_loaded(self, frame: BenchFrame) -> None:
        self.session.accept_frame(frame)
        # The frames carry the day they were shot on, which is the day the
        # calibration is dated by; a folder that did not say so is filled in
        # from them here, before anything is saved under a wrong identity.
        self._adopt_acquisition_date(frame)
        if self.campaign is not None:
            record = self.campaign.record_frame(
                frame, saturation_level=self.session.saturation_level
            )
            # The newly opened frame may belong to another lamp than the last
            # one, so the reference set follows the frame rather than lagging it.
            self.campaign.scope_alignment_to_lamp(self.session)
            if record.path not in self._file_rows:
                self._add_file_row(record.path)
            self._select_file_row(record.path)
            # The headline is also the loudest thing in the triage view, so it
            # yields to a fresh notice about the roles rather than covering it
            # when the bench re-opens a file on its own.
            self._bench_says(record.triage.headline)
        self.refresh()

    def _select_file_row(self, path: Path) -> None:
        if path not in self._file_rows:
            return
        self.file_table.selectRow(self._file_rows.index(path))
        self._show_triage(path)

    @QtCore.pyqtSlot(str, str)
    def _frame_failed(self, path: str, reason: str) -> None:
        self.session.fail_file_load(path, reason)
        self.refresh()

    @QtCore.pyqtSlot()
    def _load_finished(self) -> None:
        if self._load_thread is not None:
            self._load_thread.deleteLater()
        self._load_thread = None
        self._start_next_load()
        if self._load_thread is None and not self._queue and self._arrivals_pending:
            # The drop is fully read. Filenames are judged as a set — a folder
            # arrives as a folder — so this waits for the last file rather than
            # deciding on each one as it lands.  And it is the *arrival* that is
            # judged: opening a lamp signal for fitting is a load too, and
            # judging the folder again when that one lands would apply
            # suggestions the operator never dropped anything to ask about.
            self._arrivals_pending = False
            self._apply_unambiguous_suggestions()

    def _order_changed(self, order_idx: int) -> None:
        self.session.set_selected_order(order_idx)
        self._refresh_frame_selector()
        self.refresh_plots()

    def anchor_near(self, order_idx: int, pixel: float):
        """The anchor of *order_idx* whose stick a click at *pixel* lands on.

        Either end counts: the curated row's expected pixel *as this frame
        shows it*, where the stick is drawn, and the fitted centroid, where the
        green marker sits.  Reading the base pixel here once the sticks had
        moved would leave a hot zone under bare spectrum and none under the
        stick the operator is aiming at.
        """

        anchors = [
            anchor
            for anchor in self.session.anchor_rows()
            if anchor.line.order_idx == int(order_idx)
        ]
        shown = self.session.corrected_rows([anchor.line for anchor in anchors])
        best = None
        for anchor, row in zip(anchors, shown):
            tolerance = max(3.0, anchor.line.width_px)
            distance = min(
                abs(pixel - row.center_pixel),
                abs(pixel - anchor.fit.center_pixel),
            )
            if distance <= tolerance and (best is None or distance < best[0]):
                best = (distance, anchor)
        return None if best is None else best[1]

    def _order_plot_clicked(self, event) -> None:
        """Click a line to anchor it; click its anchored stick to un-anchor.

        F17 item 4: removal used to live only in the left table, blind to the
        plot the operator is actually looking at. Both buttons remove here, so
        the whole anchor workflow happens on the spectrum and the table follows.
        """

        buttons = (QtCore.Qt.LeftButton, QtCore.Qt.RightButton)
        if event.button() not in buttons or self.session.frame is None:
            return
        if not self.order_plot.sceneBoundingRect().contains(event.scenePos()):
            return
        point = self.order_plot.getViewBox().mapSceneToView(event.scenePos())
        order_idx = self.session.selected_order
        anchor = self.anchor_near(order_idx, point.x())
        if anchor is not None:
            self.session.remove_anchor(anchor.key)
            self.message_value.setText(
                f"{anchor.line.species} {anchor.line.wavelength_nm:.3f} nm "
                "un-anchored from the spectrum; the anchor table is in step and "
                "the alignment was recomputed. Click the line again to re-anchor it."
            )
            self.refresh()
            return
        if event.button() == QtCore.Qt.RightButton:
            self.message_value.setText(
                "No anchor sits at that pixel. Right-click an anchored line to "
                "remove it; left-click a labelled line to anchor it."
            )
            return
        result = self.session.fit_anchor_at(order_idx, point.x())
        self.message_value.setText(result.reason)
        self.refresh()

    def _auto_anchor(self) -> None:
        """Anchor every trusted line of the assigned lamp in one pass."""

        before = set(self.session.anchors)
        result = self.session.auto_anchor()
        if not result.ran:
            self.refresh()
            self.message_value.setText(f"Nothing was auto-anchored: {result.reason}")
            return

        placed = len(result.accepted)
        vetting = self.session.vetting
        against = ""
        if vetting is not None:
            against = (
                f" against the {vetting.vetted_set} vetted set"
                if vetting.is_vetted
                else " against OK marks carrying no recorded vetting"
            )
        sentences = [
            f"Auto-anchor measured {result.considered} trusted "
            f"{self._reference_label()} line(s) across all orders{against} and "
            f"anchored {placed}."
        ]
        if self.session.rms_px is not None:
            sentences.append(f"The solved alignment reads RMS {self.session.rms_px:.3f} px.")
        if result.rejected:
            sentences.append(
                f"Declined {len(result.rejected)}: {self._decline_summary(result)}."
            )
        kept = len(before - {anchor.key for anchor in result.accepted})
        if kept:
            sentences.append(f"Your {kept} hand-placed anchor(s) were kept.")
        # RMS is the anchors agreeing with each other; this is the solution
        # agreeing with physics, and it is the one the BH paper was held to.
        validation = self.session.validate_science_lines()
        sentences.append(
            f"Science-line agreement: {validation.message}."
            if validation.measured
            else f"Not yet validated against science lines — {validation.message}."
        )
        sentences.append(
            "Review the anchor table and right-click any line on the spectrum "
            "to drop it."
        )
        # The refresh reports whatever state it finds; what the operator needs
        # to read is the outcome of the press, so it is written last.
        self.refresh()
        self.message_value.setText(" ".join(sentences))

    def _reference_label(self) -> str:
        reference = self.session.reference
        if reference is not None and reference.catalog_label:
            return reference.catalog_label
        return "curated"

    @staticmethod
    def _decline_summary(result) -> str:
        """Group the declines by reason so the message stays one line long."""

        counts: dict[str, int] = {}
        for rejection in result.rejected:
            head = rejection.reason.split("(")[0].strip()
            counts[head] = counts.get(head, 0) + 1
        return ", ".join(
            f"{count} {reason}" for reason, count in sorted(counts.items())
        )

    def _remove_selected_anchor(self) -> None:
        row = self.anchor_table.currentRow()
        anchors = self.session.anchor_rows()
        if 0 <= row < len(anchors):
            self.session.remove_anchor(anchors[row].key)
            self.message_value.setText("Anchor removed; alignment recomputed.")
            self.refresh()

    def _clear_anchors(self) -> None:
        self.session.clear_anchors()
        self.message_value.setText("All anchors cleared.")
        self.refresh()

    @staticmethod
    def _short_verdict(path: Path, triage: ExposureTriage, verdict) -> str:
        """The verdict, and at most one reading that qualifies it.

        Not the filename: the file table is a few centimetres to the left with
        that row selected, and repeating it here spent a line of the loudest
        text on screen saying something already said (owner, 2026-08-16). Not
        the peak, the anomalies, and the advice either — those have their own
        panels directly underneath. Three lines of large bold red is not a
        headline, it is shouting.
        """

        if verdict.is_background:
            # A correct background says one word and stops. It was reading
            # "BACKGROUND · dark as it should be" in warning amber, which is a
            # bench explaining darkness to the person who shot it (owner: "I
            # know that darks are dark. Did you write it for yourself?"). Only
            # a background that is NOT dark has anything to add.
            return (
                f"{verdict.label.upper()} · {verdict.headline}"
                if verdict.headline
                else verdict.label.upper()
            )
        if triage.state is ExposureState.SATURATED:
            clusters = triage.saturation.cluster_count
            return f"{verdict.label.upper()} · {clusters} saturated cluster(s)"
        return verdict.label.upper()

    def _show_triage(self, path: Path) -> None:
        """Render one file's exposure verdict; roles play no part in it."""

        record = None if self.campaign is None else self.campaign.loaded.get(path)
        if record is None:
            return
        assert self.campaign is not None
        triage = record.triage
        measurement = self.campaign.measurements.get(path)
        role = measurement.role if measurement is not None else None
        verdict = triage_for_role(
            triage,
            role,
            self.campaign.partner_peak(path, role),
            _record_exposure_s(record),
        )
        color = self._reading_colour(verdict, triage)
        headline = triage.headline
        if verdict.is_background:
            headline = verdict.headline or verdict.label
        elif triage.state is ExposureState.SATURATED and not verdict.blocking:
            # The bright/dim pair exists so the dim series saturates its strong
            # lines; saying FAILED here would be the bench misreading physics.
            headline = verdict.headline.upper()
        self.triage_headline.setText(self._short_verdict(path, triage, verdict))
        self.triage_headline.setStyleSheet(f"color: {color};")
        breakdown = (
            "\n".join(triage.details)
            + "\n\nIsolated full-scale pixels are cosmic rays or hot pixels and "
            "are counted, not held against the frame. Only a connected cluster "
            "of full-scale pixels is real saturation."
        )
        whole_verdict = f"{headline}\n\n{verdict.advice}\n\n{breakdown}"
        self._explainable(
            self.triage_headline, f"{path.name} — exposure verdict", whole_verdict
        )
        self._draw_histogram(self.histogram_plot, triage, triage.histogram)
        self._draw_top_histogram(triage)
        guidance = record.exposure
        peak = "—" if guidance.peak_value is None else f"{guidance.peak_value:.0f} counts"
        # The one line that decides the next action stays on screen; the four
        # lines of arithmetic behind it live in the dock.
        # A background must never be told to expose longer: its role decides
        # the action, and the role-blind guidance only speaks when the role has
        # nothing to add (owner, on a lamp background reading "increase toward
        # 19.98 s": "Background is dim. Wow, genius").
        if verdict.is_background:
            # A background that is dark needs no next action, so it gets no
            # panel — not a panel saying there is nothing to do. The role-blind
            # guidance must never be reached here: it would say "increase
            # exposure toward 19.98 s" about a frame shot with the light off.
            action = verdict.next_action
        else:
            action = verdict.next_action or guidance.next_action
        self.triage_next_value.setText(action)
        self.triage_next_value.setVisible(bool(action))
        self._explainable(
            self.triage_next_value,
            f"{path.name} — how the verdict was reached",
            breakdown,
        )
        # Selecting a file writes its whole verdict into the dock at once:
        # what it says, what to do about it, and the numbers behind both.
        self.explain(f"{path.name} — exposure verdict", whole_verdict)
        # The numbers, and only the numbers. This panel used to repeat the
        # filename, the state and the advice — all three already on screen
        # above it — which is how one file's verdict came to be stated three
        # times in one column (owner, 2026-08-16).
        numbers = [f"peak {peak}"]
        if triage.headroom_fraction is not None:
            numbers.append(f"{100.0 * triage.headroom_fraction:.0f}% of full scale")
        if triage.saturation.anomalous_pixels:
            numbers.append(f"{triage.saturation.anomalous_pixels} anomalies")
        self.exposure_value.setText(" · ".join(numbers))
        self._explainable(
            self.exposure_value,
            f"{path.name} — the numbers behind the verdict",
            breakdown,
        )

    def _draw_histogram(self, plot, triage: ExposureTriage, histogram) -> None:
        """Draw one counts histogram as outlined bins on a log floor.

        The plot is in log mode, so the fill level has to be given in log
        space: filling an empty bin down to y=1 is what turned a histogram of
        zeros into a solid painted block.  Bins are floored just below one
        pixel and outlined, so a bin holding a handful of pixels is visible and
        a bin holding none draws nothing.

        The curve and the two threshold lines are built once per plot and then
        moved, never cleared and rebuilt.  ``clear()`` unparents an
        ``InfiniteLine`` while Qt may still ask it to paint, and a removed line
        answering ``boundingRect`` has no view box to measure itself against —
        which raises inside a paint and aborts the whole process rather than
        failing visibly.  Selecting each file in turn is enough to reach it.
        """

        floor = 0.5
        counts = np.maximum(np.asarray(histogram.counts, dtype=float), floor)
        edges = np.asarray(histogram.edges, dtype=float)
        curve, saturation_line, full_scale_line = self._histogram_furniture(plot, floor)
        curve.setData(edges, counts, stepMode="center", fillLevel=float(np.log10(floor)))
        saturation_line.setValue(triage.saturation.saturation_level)
        full_scale_line.setValue(triage.full_scale)

    def _histogram_furniture(self, plot, floor: float) -> tuple:
        """The curve and threshold lines one histogram plot owns for its life."""

        furniture = self._histogram_items.get(id(plot))
        if furniture is None:
            curve = plot.plot(
                np.array([0.0, 1.0]),
                np.array([floor]),
                stepMode="center",
                fillLevel=float(np.log10(floor)),
                brush=pg.mkBrush("#1f4d63"),
                pen=pg.mkPen("#76d6ff", width=1.4),
            )
            furniture = (
                curve,
                plot.addLine(x=0.0, pen=pg.mkPen("#ffb86b", style=QtCore.Qt.DashLine)),
                plot.addLine(x=0.0, pen=pg.mkPen("#ff8f8f", style=QtCore.Qt.DashLine)),
            )
            self._histogram_items[id(plot)] = furniture
        return furniture

    def _draw_top_histogram(self, triage: ExposureTriage) -> None:
        """Draw the top end, or say in words that nothing is up there.

        F17 item 3: an all-zero log histogram rendered as a solid block that
        conveyed nothing.  Emptiness here is a real and good answer — the frame
        has headroom — so it is written out instead of drawn.
        """

        occupied = int(np.count_nonzero(np.asarray(triage.top_histogram.counts)))
        if occupied == 0:
            self.top_histogram_widget.setVisible(False)
            self.top_end_message.setText(_TOP_END_EMPTY)
            self.top_end_message.setVisible(True)
            return
        self.top_end_message.setVisible(False)
        self.top_histogram_widget.setVisible(True)
        self._draw_histogram(self.top_histogram_plot, triage, triage.top_histogram)

    def _describe_input_source(self) -> None:
        """Name WHICH calibration is on the bench, on the Bench state strip.

        This line used to describe the gesture — "manual — drag and drop or Add
        files" — which is the one thing the operator already knew, while the
        folder that every root, every dialog and the identity's own date are
        derived from went unnamed on a strip called Bench state.  The folder
        leads; the watch, being the unusual mode, says so after it.
        """

        if self.calibration_folder is not None:
            watched = " · watching for new SIFs" if self.watcher is not None else ""
            self.watch_value.setText(f"{self.calibration_folder}{watched}")
        elif self.watcher is not None:
            self.watch_value.setText(f"watching {self.watcher.folder}")
        else:
            self.watch_value.setText("no folder open — drag and drop or Add files")

    def refresh(self) -> None:
        """Render the current domain state without changing it."""

        self._describe_input_source()
        if self.session.frame is not None:
            self.file_value.setText(self.session.frame.path.name)
        elif self.session.loading_path is not None:
            self.file_value.setText(self.session.loading_path.name)
        self.file_state_value.setText(self.session.file_state.value.upper())
        self.alignment_state_value.setText(
            self.session.alignment_state.value.replace("-", " ").upper()
        )
        anchors = self.session.anchor_rows()
        self.anchor_count_value.setText(str(len(anchors)))
        self.rms_value.setText(
            "—" if self.session.rms_px is None else f"{self.session.rms_px:.3f} px"
        )
        if self.session.transform is None:
            self.transform_value.setText("—")
        else:
            transform = self.session.transform
            self.transform_value.setText(
                f"{transform.dx_px:+.2f}, {transform.dy_px:+.2f}, "
                f"{transform.theta_deg:+.3f}°"
            )
        if self.session.file_state is FileLoadState.FAILED:
            self.message_value.setText(
                f"Could not read that SIF: {self.session.last_error}. "
                "The last good frame remains visible; drop another file when ready."
            )
        elif self.session.file_state is FileLoadState.LOADING:
            self.message_value.setText("Reading SIF, extracting orders, triaging exposure…")
        elif self.session.alignment_state is AlignmentState.COLLECTING:
            self.message_value.setText("One anchor accepted. Add another order/line to solve.")
        elif self.session.alignment_state is AlignmentState.ALIGNED:
            self.message_value.setText(
                "Live rigid alignment solved. Review RMS and residuals as anchors accumulate."
            )
        elif self.session.alignment_state is AlignmentState.FAILED:
            self.message_value.setText(f"Alignment failed: {self.session.last_error}")
        self._refresh_frame_selector()
        self._follow_assigned_lamp_catalog()
        self._land_on_the_lamps_own_orders()
        self._refresh_reference()
        self._refresh_fit_target()
        self._refresh_anchor_table()
        self.refresh_plots()
        self.refresh_campaign()

    def _refresh_frame_selector(self) -> None:
        """Offer the mean and every single frame of the open acquisition."""

        count = self.session.frame_count
        wanted = [("Mean of all frames" if count > 1 else "Single frame", None)]
        wanted.extend((f"Frame {index + 1} of {count}", index) for index in range(count))
        present = [
            (self.frame_combo.itemText(index), self.frame_combo.itemData(index))
            for index in range(self.frame_combo.count())
        ]
        if present != wanted:
            self.frame_combo.blockSignals(True)
            try:
                self.frame_combo.clear()
                for label, data in wanted:
                    self.frame_combo.addItem(label, data)
                index = self.frame_combo.findData(self.session.selected_frame)
                self.frame_combo.setCurrentIndex(max(0, index))
            finally:
                self.frame_combo.blockSignals(False)
        self.frame_combo.setEnabled(count > 1)
        tail = "" if count > 1 else " — this acquisition holds one frame only"
        self.frame_choice_value.setText(
            f"Fitting {self.session.frame_choice_label()}{tail}"
        )
        orders = self.session.pattern.shape[1]
        self.order_total_value.setText(f"of {orders - 1}")
        self.previous_order_button.setEnabled(self.order_spin.value() > 0)
        self.next_order_button.setEnabled(self.order_spin.value() < orders - 1)

    def refresh_campaign(self) -> None:
        """Render campaign memory without inferring or mutating measurement roles."""

        enabled = self.campaign is not None
        busy = self._campaign_thread is not None
        self.add_files_button.setEnabled(self.loader is not None)
        self.open_folder_button.setEnabled(self.loader is not None and not busy)
        # And a press while a read is in flight is refused in the handler, where
        # the answer is always current — see ``_pick_calibration_folder``.
        self._refresh_file_buttons()
        self.generate_tomls_button.setEnabled(enabled and not busy)
        # These run their handlers on the GUI thread, on the very campaign and
        # session a running task is reading.  They belong in the same gate as
        # the buttons beside them.
        self.next_step_button.setEnabled(not busy)
        self.regenerate_tomls_button.setEnabled(enabled and not busy)
        self.recompute_sphere_button.setEnabled(
            enabled and not busy and self._sphere_pair_assigned()
        )
        self._refresh_pattern_source(busy)
        self._refresh_previous_pair(busy)
        # An offer to overwrite belongs to the identity that was refused.  The
        # id field is edited between the refusal and the press, and this runs on
        # every keystroke in it, so a different identity withdraws the offer
        # rather than inheriting it.
        if self._refused_identity != self.snapshot_id_edit.text().strip():
            self._refused_identity = ""
            self.regenerate_tomls_button.setVisible(False)
        if not enabled:
            self.checklist_tree.clear()
            self.confirm_roles_button.setEnabled(False)
            self.comparison_value.setText("Campaign memory was not configured.")
            self.generate_tomls_button.setEnabled(False)
            self.save_snapshot_button.setEnabled(False)
            return

        assert self.campaign is not None
        pending = self.campaign.unconfirmed_suggestions()
        self.confirm_roles_button.setEnabled(bool(pending) and not busy)
        self.confirm_roles_button.setText(
            f"Confirm {len(pending)} suggested role(s)"
            if pending
            else "Confirm suggested roles"
        )
        self._refresh_file_table()
        self._refresh_triage_summary()
        self._refresh_comparison_summary()
        self._refresh_checklist()
        self.save_state_value.setText(
            self.campaign.save_state.value.replace("-", " ").upper()
        )
        self.save_snapshot_button.setEnabled(
            not busy
            and self.campaign.ready_for_snapshot(
                self.snapshot_id_edit.text().strip(), self.session
            )
        )
        if self.campaign.toml_state is TomlState.GENERATED and not self.toml_preview.toPlainText():
            campaign_path = self.campaign.toml_paths.get("campaign")
            if campaign_path is not None and campaign_path.is_file():
                self.toml_preview.setPlainText(
                    campaign_path.read_text(encoding="utf-8")
                )
        self._refresh_sphere_plot()
        # The fit belongs on the assigned lamp signal, minus its own
        # background — decided here, after the roles are known.
        self._follow_assigned_lamp_signal()
        self._pair_lamp_background()

    def _refresh_pattern_source(self, busy: bool) -> None:
        """Say which pattern the bench wears, and what may be done about it.

        Quiet by default: the path on the readings strip and beside the factor
        curves, one picker, and an extraction that is only emphasised while the
        band guard is actually warning about the pattern it names.
        """

        if self.campaign is None:
            self.pattern_source_value.setText("packaged pattern (no campaign memory)")
            self.pattern_choice_value.setText("Pattern: —")
            self.choose_pattern_button.setEnabled(False)
            self.extract_pattern_button.setEnabled(False)
            return
        source = self.campaign.pattern_source
        self.pattern_source_value.setText(str(source))
        # The whole path on hover, like every other elided reading; the tooltip
        # rule about short lines is about explanations, and a path is a fact.
        self.pattern_choice_value.setText(f"Pattern: {source}")
        self.pattern_choice_value.setToolTip(str(source))
        self._explainable(
            self.pattern_source_value,
            "The pattern this bench is wearing",
            "Every order the bench extracts, every anchor row it reads, and "
            "every factor it sums is taken off this file. It is no longer fixed "
            "at launch: the Sphere factors view can extract one from the "
            f"assigned sphere or open another, in place. In full: {source}",
            hint=one_line(f"{source.name} — click for the whole path and why it matters"),
        )
        self.choose_pattern_button.setEnabled(not busy)
        blocker = self._extraction_blocker()
        self.extract_pattern_button.setEnabled(not busy and not blocker)
        warning = self.campaign.pattern_band_warning()
        self.extract_pattern_button.setToolTip(
            f"Not yet possible: {blocker}."
            if blocker
            else (
                "The sphere's bands do not sit on this pattern; fit one to them."
                if warning
                else "Fit a pattern to the assigned sphere pair and stand on it."
            )
        )
        # Loud only when something is wrong. A pattern the light fits is not a
        # thing to be pressed about.
        _emphasise(self.extract_pattern_button, self.body_pt, bold=bool(warning))

    @staticmethod
    def _verdict_colour(verdict, triage) -> str:
        """The colour a frame's ROLE earns it, not the one its raw state does.

        A background is DIM by construction and DIM is amber, so every correct
        dark frame was painted in the colour that means attend to me.  It then
        earned the GOOD green, which said the right thing but said it in the
        same voice as a healthy signal: on the bench's dark ground the operator
        could no longer tell a background row from a lamp row at a glance.  So
        a correct background now has a colour of its own — calm, cool, and
        nobody's alarm — and a background that is NOT dark still earns red.

        This is the one place that rule is written; the triage panel, the files
        table and the folder summary all read it here, and it is keyed on the
        role the verdict carries so that rewording a badge can never repaint a
        frame.
        """

        if verdict is not None and verdict.is_background:
            if verdict.blocking:
                return _TRIAGE_COLORS[ExposureState.SATURATED]
            return _BACKGROUND_COLOR
        return _TRIAGE_COLORS[triage.state]

    @classmethod
    def _reading_colour(cls, verdict, triage) -> str:
        """The colour one frame's reading is painted, wherever it is read.

        The role decides it, and expected saturation is demoted to a note:
        the dim series of a bright/dim pair is *shot* to clip its strong lines,
        so painting it alarm red is the bench misreading physics.  Written once
        because the triage headline, the files table and the folder summary all
        have to answer this the same way.
        """

        if (
            verdict.state is ExposureState.SATURATED
            and not verdict.blocking
            and not verdict.is_background
        ):
            return _TRIAGE_COLORS[ExposureState.DIM]
        return cls._verdict_colour(verdict, triage)

    @staticmethod
    def _role_text(measurement, is_suggested: bool) -> str:
        """What the summary's own Role column says, in its own column.

        The files table folds the role in beside the verdict and therefore
        drops it when the verdict already names it; a column cannot do that —
        an empty cell under a heading reads as missing, not as redundant.
        """

        if measurement is None:
            return _SUGGESTED_BADGE.lower() if is_suggested else "no role yet"
        if measurement.lamp_family:
            return f"{measurement.lamp_family} {measurement.role.value}"
        return measurement.role.value

    @staticmethod
    def _role_marks(measurement, verdict, is_suggested: bool) -> list[str]:
        """What one row says about its role, without saying it twice."""

        if measurement is None:
            # One spot, not six (owner, 2026-08-16): "SUGGESTED ONLY — no role
            # assigned" spelled out on every row filled the table with one
            # sentence repeated, and it was already said twice over — by the
            # amber border the Role control wears, and by the Confirm button
            # that counts them. A row says what it is; the colour says whether
            # anybody has confirmed it.
            return [_SUGGESTED_BADGE.lower() if is_suggested else "no role yet"]
        if verdict is not None and verdict.is_background:
            # The verdict already said "background"; the role would say it a
            # second time, which is the duplication this table just lost.
            return [measurement.lamp_family] if measurement.lamp_family else []
        if measurement.lamp_family:
            return [f"{measurement.lamp_family} {measurement.role.value}"]
        return [measurement.role.value]

    def _refresh_file_table(self) -> None:
        """Show each loaded file's verdict and the role it currently carries."""

        assert self.campaign is not None
        suggested = dict(self.campaign.unconfirmed_suggestions())
        for row, path in enumerate(self._file_rows):
            record = self.campaign.loaded.get(path)
            measurement = self.campaign.measurements.get(path)
            verdict = None
            marks = []
            if record is not None:
                row_role = measurement.role if measurement is not None else None
                verdict = triage_for_role(
                    record.triage,
                    row_role,
                    self.campaign.partner_peak(path, row_role),
                    _record_exposure_s(record),
                )
                marks.append(verdict.label.upper())
                if record.triage.saturation.anomalous_pixels:
                    marks.append(f"{record.triage.saturation.anomalous_pixels} anomalies")
            marks.extend(self._role_marks(measurement, verdict, path in suggested))
            if path == getattr(self.session.frame, "path", None):
                marks.append("open")
            item = self.file_table.item(row, 0)
            item.setText(f"{path.name}\n{' · '.join(marks)}" if marks else path.name)
            role_combo = self.file_table.cellWidget(row, 1)
            self._render_role_combo(role_combo, path in suggested)
            if record is not None and verdict is not None:
                colour = (
                    _SUGGESTED_COLOR
                    if measurement is None and path in suggested
                    else self._reading_colour(verdict, record.triage)
                )
                item.setForeground(QtGui.QColor(colour))
                title = "exposure verdict"
                explanation = f"{path}\n{verdict.headline}\n\n{verdict.advice}"
                if measurement is None and path in suggested:
                    title = "suggested role only"
                    explanation = f"{path}\n{_SUGGESTED_TOOLTIP}\n\n{verdict.headline}"
                # One short line hovers; the whole thing goes to the dock.
                item.setToolTip(one_line(verdict.headline))
                item.setData(QtCore.Qt.UserRole, title)
                item.setData(QtCore.Qt.UserRole + 1, explanation)
            lamp_combo = self.file_table.cellWidget(row, 2)
            lamp_combo.setEnabled(
                measurement is not None and measurement.role in _LAMP_ROLES
            )

    def _fit_triage_summary_height(self) -> None:
        """Cap the summary at the folder it summarises, never at the view."""

        table = self.triage_summary_table
        rows = max(1, min(table.rowCount(), BENCH_SUMMARY_ROWS))
        table.setMaximumHeight(
            table.horizontalHeader().height()
            + rows * table.verticalHeader().defaultSectionSize()
            + 2 * table.frameWidth()
        )

    def _summary_row(self, path: Path, is_suggested: bool) -> tuple:
        """One file's summary line: its cells, its colour, and its one-liner."""

        assert self.campaign is not None
        record = self.campaign.loaded.get(path)
        measurement = self.campaign.measurements.get(path)
        role_text = self._role_text(measurement, is_suggested)
        if record is None:
            return (path.name, role_text, "—", "—", "—"), "", ""
        role = measurement.role if measurement is not None else None
        verdict = triage_for_role(
            record.triage,
            role,
            self.campaign.partner_peak(path, role),
            _record_exposure_s(record),
        )
        peak = (
            "—"
            if record.triage.headroom_fraction is None
            else f"{100.0 * record.triage.headroom_fraction:.0f}%"
        )
        counted = record.triage.saturation.anomalous_pixels
        return (
            (
                path.name,
                role_text,
                verdict.label.upper(),
                peak,
                str(counted) if counted else "—",
            ),
            self._reading_colour(verdict, record.triage),
            verdict.headline,
        )

    def _refresh_triage_summary(self) -> None:
        """Summarise every loaded file on one line each, in verdict colour.

        The detailed panel below still reads one file — the selected one — and
        this is how it gets selected: the operator sees which row is worth
        opening instead of clicking every row to find out (F21 item 4a).
        """

        assert self.campaign is not None
        table = self.triage_summary_table
        suggested = dict(self.campaign.unconfirmed_suggestions())
        if table.rowCount() != len(self._file_rows):
            table.setRowCount(len(self._file_rows))
            self._fit_triage_summary_height()
        for row, path in enumerate(self._file_rows):
            cells, colour, headline = self._summary_row(path, path in suggested)
            for column, text in enumerate(cells):
                item = table.item(row, column)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    if column:
                        item.setTextAlignment(QtCore.Qt.AlignCenter)
                    table.setItem(row, column, item)
                item.setText(text)
            table.item(row, 0).setToolTip(one_line(path.name))
            if colour:
                # The verdict word carries the colour, not the whole row: a row
                # painted end to end is a highlight, and a highlight on every
                # row is no highlight at all.
                table.item(row, 2).setForeground(QtGui.QColor(colour))
                table.item(row, 2).setToolTip(one_line(headline))
        self._sync_summary_selection()

    def _sync_summary_selection(self) -> None:
        """Highlight the summary row of whatever file is being read below."""

        selected = self._selected_file()
        if selected is None or selected not in self._file_rows:
            return
        index = self._file_rows.index(selected)
        table = self.triage_summary_table
        if table.currentRow() == index:
            return
        table.blockSignals(True)
        try:
            table.selectRow(index)
        finally:
            table.blockSignals(False)

    def _summary_row_clicked(self, row: int, _column: int = 0) -> None:
        """Open the file whose summary row was clicked in the panel below."""

        if 0 <= row < len(self._file_rows):
            self._select_file_row(self._file_rows[row])

    def _render_role_combo(self, role_combo, is_only_suggested: bool) -> None:
        """Make an unconfirmed suggestion impossible to read as an assignment.

        The combo's *text* is never touched: it holds a short role label and
        nothing else, so it cannot elide away the very state the operator has
        to read.  The SUGGESTED state is carried outside the control — by its
        colour here, and by the badge in the file cell beside it.
        """

        if role_combo is None:
            return
        self._explainable(
            role_combo,
            "Suggested is not assigned"
            if is_only_suggested
            else "The role this file carries",
            _SUGGESTED_TOOLTIP
            if is_only_suggested
            else "The measurement role this file carries. Change it freely; "
            "the bench never infers one behind your back.",
        )
        role_combo.setStyleSheet(
            f"border: 1px solid {_SUGGESTED_COLOR}; color: {_SUGGESTED_COLOR};"
            if is_only_suggested
            else ""
        )
        role_combo.setMinimumWidth(role_combo_minimum_width(role_combo))

    def _file_row_clicked(self, item) -> None:
        """Clicking a file row writes its whole verdict into the Why dock.

        This used to fire on hover, which meant the dock rewrote itself
        whenever the pointer crossed the files table on its way elsewhere.
        """

        row = item.row()
        if 0 <= row < len(self._file_rows):
            path = self._file_rows[row]
            title = item.data(QtCore.Qt.UserRole)
            if title:
                self.explain(f"{path.name} — {title}", item.data(QtCore.Qt.UserRole + 1) or "")

    def _refresh_comparison_summary(self) -> None:
        assert self.campaign is not None
        comparison = self.campaign.comparison
        if comparison.state is ComparisonState.READY:
            text = (
                "READY · new/previous median "
                f"{comparison.median_ratio:.3f}; 5–95% "
                f"{comparison.p05_ratio:.3f}–{comparison.p95_ratio:.3f} "
                f"({comparison.sample_count} samples)."
            )
            # The ratio and the pair it was measured against are one reading:
            # a median of 1.000 said alone is unfalsifiable, and the operator
            # who read exactly that off his own folder was right to doubt it.
            if comparison.reference_name:
                text += f"\nvs {comparison.reference_name}"
            if comparison.self_comparison:
                text += f"\nSELF-CHECK · {SELF_COMPARISON_NOTE}."
            self.comparison_value.setText(text)
        else:
            self.comparison_value.setText(
                f"{comparison.state.value.replace('-', ' ').upper()} · {comparison.reason}"
            )
        self._explain_comparison(comparison)

    def _explain_comparison(self, comparison) -> None:
        """Keep the Why dock's sphere-factor text on the real reference files.

        The panel names the pair; the dock carries the absolute paths, which
        is what settles whether "previous" is really previous.
        """

        text = _SPHERE_FACTORS_EXPLANATION
        hint = ""
        if comparison.reference_name:
            text += (
                "\n\nCompared against the previous pair:\n"
                f"{comparison.previous_sphere}\n"
                f"{comparison.previous_sphere_background}"
            )
            hint = f"Compared against {comparison.reference_name}"
        if comparison.self_comparison:
            text += f"\n\nSELF-CHECK: {SELF_COMPARISON_NOTE}."
            hint = f"Self-check — {comparison.reference_name} is a copy of this sphere pair"
        self._explainable(
            self.comparison_value, _SPHERE_FACTORS_TITLE, text, hint=hint
        )

    def _refresh_checklist(self) -> None:
        assert self.campaign is not None
        while self._checklist_labels:
            self._forget_explainable(self._checklist_labels.pop())
        self.checklist_tree.clear()
        symbols = {
            ChecklistState.DONE: "✓",
            ChecklistState.WAITING: "○",
            ChecklistState.ATTENTION: "!",
            ChecklistState.SUGGESTION: "·",
        }
        colors = {
            ChecklistState.DONE: QtGui.QColor("#70d6ae"),
            ChecklistState.WAITING: QtGui.QColor("#8fa5b5"),
            ChecklistState.ATTENTION: QtGui.QColor("#ffb86b"),
            ChecklistState.SUGGESTION: QtGui.QColor("#8fd9ff"),
        }
        width = self._checklist_row_width()
        checklist = self.campaign.checklist(self.session)
        for item in checklist:
            row = QtWidgets.QListWidgetItem()
            explanation = item.detail
            if item.unblocked_by:
                explanation += f"\n\nWhat unblocks it: {item.unblocked_by}"
            if not item.blocking:
                explanation += "\n\nThis row is advice. It never blocks anything."
            row.setToolTip(one_line(item.detail))
            label = QtWidgets.QLabel(
                self._checklist_row_html(
                    symbols[item.state],
                    f"{item.label} — {item.detail}",
                    item.unblocked_by or "",
                )
            )
            label.setTextFormat(QtCore.Qt.RichText)
            label.setWordWrap(True)
            label.setContentsMargins(8, 6, 8, 6)
            label.setStyleSheet(f"color: {colors[item.state].name()};")
            _emphasise(label, self.body_pt, bold=item.state is ChecklistState.ATTENTION)
            self._explainable(label, item.label, explanation)
            self._checklist_labels.append(label)
            self.checklist_tree.addItem(row)
            self.checklist_tree.setItemWidget(row, label)
            self._fit_checklist_row(row, label, width)
        self._refresh_next_step(checklist)

    def _next_action_for(self, step) -> tuple[str, object] | None:
        """The verb that performs one checklist step, and what performs it.

        Every one of these already existed as a button somewhere on the bench.
        What was missing was any relationship between the row that says a step
        is next and the control that does it.
        """

        key = step.key
        if key == "references":
            # This row used to explain itself and offer nothing, because the
            # tables were named when the bench was opened and no file the SIF
            # picker could add would satisfy it.  That is no longer true of the
            # pattern: when the guard says the sphere's bands do not sit on it,
            # the bench holds the sphere, packages the extraction, and can fit
            # and adopt a pattern without closing.  So the row that reports the
            # mismatch carries the verb that answers it.
            if (
                self.campaign is not None
                and self.campaign.pattern_band_warning()
                and not self._extraction_blocker()
            ):
                return (
                    "Extract pattern from this sphere",
                    self._extract_pattern_from_sphere,
                )
            return None
        if key == "files":
            return ("Add SIF files…", self._pick_files)
        pending = () if self.campaign is None else self.campaign.unconfirmed_suggestions()
        if pending:
            # Whatever the row is, nothing can be assigned until the roles are.
            return (f"Confirm {len(pending)} role(s)", self._confirm_suggested_roles)
        if key == "sphere-comparison":
            return ("Measure sensitivity", self._start_sphere_comparison)
        if key == "alignment":
            return ("Auto-anchor lines", self._auto_anchor)
        if key == "tomls":
            return ("Save alignment settings", lambda: self._generate_tomls())
        if key == "snapshot":
            return ("Save and validate", self._save_snapshot)
        return None

    def _set_next_action(self, step) -> None:
        """Put the verb for this step on the button beside it."""

        action = None if step is None else self._next_action_for(step)
        if action is None:
            self.next_step_button.setVisible(False)
            self._next_action = None
            return
        label, handler = action
        self._next_action = handler
        self.next_step_button.setText(label)
        self.next_step_button.setVisible(True)

    def _run_next_action(self) -> None:
        # The same gate the button's enabled state carries: a running campaign
        # task reads the campaign and the session this handler would mutate.
        if self._campaign_thread is not None:
            return
        if self._next_action is not None:
            self._next_action()

    def _refresh_next_step(self, checklist) -> None:
        """The one line the checklist is usually consulted for.

        Not a copy of the procedure — the first row that is blocking and not
        yet done, with whatever the checklist says would unblock it.  That is
        the question an operator crosses the room to ask.  The rows come from
        the caller, which has just built them: deriving the checklist is
        O(files x measurements) and this runs on every keystroke in the
        snapshot-id field.
        """

        if self.campaign is None:
            return
        pending = [
            item
            for item in checklist
            if item.blocking and item.state is not ChecklistState.DONE
        ]
        if not pending:
            self.next_step_value.setText(
                "Next: nothing — every step of the procedure is done."
            )
            self._set_next_action(None)
            self._explainable(
                self.next_step_value,
                "The procedure is complete",
                "Every blocking row on the Procedure tab is done. What remains "
                "is saving the snapshot, and anything the checklist marked as "
                "advice rather than a gate.",
            )
            return
        step = pending[0]
        self.next_step_value.setText(
            f"Next: {step.label} — {step.unblocked_by or step.detail}"
        )
        self._set_next_action(step)
        self._explainable(
            self.next_step_value,
            f"Next step — {step.label}",
            f"{step.detail}\n\nWhat unblocks it: {step.unblocked_by}"
            if step.unblocked_by
            else step.detail,
        )

    def _refresh_sphere_plot(self) -> None:
        """Feed the two persistent factor curves; never rebuild the plot.

        The factor curves hold tens of thousands of samples.  Recreating them
        on every refresh, drawing every sample, and antialiasing the result is
        what turned a pan into a slideshow.
        """

        if self.campaign is None:
            return
        comparison = self.campaign.comparison
        self._set_factor_curve(self.candidate_curve, comparison.candidate)
        self._set_factor_curve(self.previous_curve, comparison.previous)
        message = (
            f"{comparison.state.value.replace('-', ' ').upper()} — {comparison.reason}"
        )
        # Where this sphere's light actually sits on the chosen pattern, read
        # beside the factors it produced: the factor curve is computed by
        # summing each order over that pattern, so a pattern the bands do not
        # fit is a fact about these very curves.
        warning = self.campaign.pattern_band_warning()
        reading = self.campaign.sphere_band_offsets()
        if warning:
            message += f"\nPATTERN DOES NOT FIT — {warning}."
        elif reading is not None and reading.measured:
            message += f"\n{reading.summary()[0].upper()}{reading.summary()[1:]}."
        self.sphere_view_message.setText(message)

    @staticmethod
    def _set_factor_curve(curve: pg.PlotDataItem, result) -> None:
        """Update one factor curve in place, gaps drawn as gaps."""

        if result is None:
            curve.setData([], [])
            curve.setVisible(False)
            return
        size = min(result.wavelength_nm.size, result.factors_wmsr.size)
        curve.setData(
            np.asarray(result.wavelength_nm[:size], dtype=float),
            np.asarray(result.factors_wmsr[:size], dtype=float),
            connect="finite",
        )
        curve.setVisible(True)

    def _refresh_anchor_table(self) -> None:
        anchors = self.session.anchor_rows()
        residual_by_key = {item.key: item for item in self.session.residuals}
        self.anchor_table.setRowCount(len(anchors))
        for row, anchor in enumerate(anchors):
            residual = residual_by_key.get(anchor.key)
            values = (
                str(anchor.line.order_idx),
                f"{anchor.line.wavelength_nm:.3f}",
                f"{anchor.fit.center_pixel - anchor.line.center_pixel:+.2f}",
                "—" if residual is None else f"{residual.magnitude_px:.3f}",
                "clear",
            )
            for column, value in enumerate(values):
                self.anchor_table.setItem(row, column, QtWidgets.QTableWidgetItem(value))
        enabled = bool(anchors)
        self.remove_button.setEnabled(enabled)
        self.clear_button.setEnabled(enabled)
        self._refresh_auto_anchor()

    def _refresh_auto_anchor(self) -> None:
        """Say what a press would measure, before it is pressed."""

        # The auto pass needs a frame, not an anchor: it is what fills an empty
        # table, so gating it on a non-empty one would be exactly backwards.
        session = self.session
        reference = session.reference
        if session.frame is None:
            self.auto_anchor_button.setEnabled(False)
            self.auto_anchor_value.setText("no acquisition open")
            return
        if reference is not None and not reference.is_referenceable:
            self.auto_anchor_button.setEnabled(False)
            self.auto_anchor_value.setText(reference.message)
            return
        rows = session.trusted_rows()
        self.auto_anchor_button.setEnabled(bool(rows))
        if not rows:
            self.auto_anchor_value.setText(
                f"no {self._reference_label()} row in this table carries an OK "
                "mark — click the lines you trust instead"
            )
            return
        orders = len({line.order_idx for line in rows})
        ready = (
            f"{len(rows)} {self._reference_label()} line(s) across "
            f"{orders} order(s) are ready to be measured in one pass"
        )
        # Whose vetting these OK marks carry is the whole reason to trust the
        # pass, so its absence is stated exactly as loudly as its presence.
        vetting = self.session.vetting
        if vetting is None:
            self.auto_anchor_value.setText(ready)
        elif vetting.is_vetted:
            self.auto_anchor_value.setText(
                f"{ready}, carrying the {vetting.vetted_set} vetting"
            )
        else:
            self.auto_anchor_value.setText(f"{ready}, but {vetting.message}")

    def refresh_plots(self) -> None:
        frame = self.session.frame
        if frame is None:
            return
        self._refresh_detector_view(frame)
        self._refresh_pattern_traces()
        self._refresh_order_view()

    def _refresh_detector_view(self, frame: BenchFrame) -> None:
        """Upload the detector image once per frame, never per order.

        Percentiling a 2560x2160 image and re-uploading it costs tens of
        milliseconds; an order change does not alter one pixel of it, so it is
        keyed to the frame's identity and skipped otherwise.
        """

        key = (str(frame.path), id(frame.detector_image))
        if key == self._detector_key:
            return
        finite = frame.detector_image[np.isfinite(frame.detector_image)]
        levels = None
        if finite.size:
            low, high = np.percentile(finite, (1.0, 99.8))
            if high > low:
                levels = (float(low), float(high))
        self.detector_image.setAutoDownsample(True)
        self.detector_image.setImage(
            frame.detector_image.T,
            autoLevels=levels is None,
            levels=levels,
        )
        self._detector_key = key

    def _refresh_order_view(self) -> None:
        """Redraw only what an order change actually changes.

        The main GUI scrolls a 1 GB acquisition instantly because it extracts
        every order once and then re-slices.  This does the same: the detector
        image, its levels, and the order traces belong to the frame, so an
        order switch moves pens and pooled labels and touches nothing else.
        """

        order_idx = self.session.selected_order
        spectra = self.session.active_order_spectra()
        if order_idx >= len(spectra):
            return
        spectrum = spectra[order_idx]
        self.order_curve.setData(np.arange(spectrum.size), spectrum)
        self.order_plot.setTitle(
            f"Order {order_idx} · {self.session.frame_choice_label()}: "
            "click a labeled expected line"
        )
        self._select_pattern_trace(order_idx)
        top = float(np.nanmax(spectrum)) if np.any(np.isfinite(spectrum)) else 1.0

        anchors = [
            anchor
            for anchor in self.session.anchor_rows()
            if anchor.line.order_idx == order_idx
        ]
        anchored_nm = [anchor.line.wavelength_nm for anchor in anchors]

        # One list, one pool: the sticks the operator reads and the rows the
        # panel lists are literally the same objects (F17 item 2).
        expected = self._expected_rows_for_order(order_idx)
        self._catalog_rows = expected
        for index, row in enumerate(expected):
            marker, label = self._pooled_marker(
                self._line_pool, index, self._expected_pen, "#f5c56f"
            )
            is_anchored = any(
                abs(row.wavelength_nm - value) <= 0.05 for value in anchored_nm
            )
            # An anchored stick wears the anchor's colour, because clicking it
            # is now how an anchor is taken back off (F17 item 4).
            marker.setPen(self._anchored_pen if is_anchored else self._expected_pen)
            marker.setValue(row.detector_pixel)
            label.setText(row.label)
            label.setColor("#6ee7b7" if is_anchored else "#f5c56f")
            label.setPos(
                row.detector_pixel, top * (0.90, 0.72, 0.54, 0.36)[index % 4]
            )
        self._hide_pooled(self._line_pool, len(expected))

        self._anchor_scatter.setData(
            [anchor.fit.center_pixel for anchor in anchors],
            [anchor.fit.amplitude + anchor.fit.baseline for anchor in anchors],
        )
        self.line_highlight.setVisible(False)
        self._refresh_line_help_table(expected)
        self._refresh_residual_plot()

    def _display_reference(self) -> LampReferenceSet | None:
        """The one reference set the sticks, the rows, and the counts read.

        Normally it is the session's own — the assigned lamp's scoped rows,
        exactly what a click can snap to.  A hand-picked catalog is an override
        for comparison, and it moves the sticks with the table rather than only
        one of them; the panel says loudly that the anchors still reference the
        assigned lamp.
        """

        if self._family_override:
            return lamp_reference_set(
                self.line_family_combo.currentText(), self.session.lines
            )
        return self.session.reference

    def _correction_key(self) -> tuple[float, float, float] | None:
        """The solved correction the drawn positions currently carry."""

        transform = self.session.transform
        if transform is None:
            return None
        return (transform.dx_px, transform.dy_px, transform.theta_rad)

    def _expected_rows_for_order(self, order_idx: int) -> tuple:
        """Cache the expected-line list per (order, reference) pair.

        Rebuilding it on every arrow-key press is work the answer does not
        change with, and F16's order-scrolling budget has no room for it.

        The rows are placed where this frame shows them — base pixels until a
        transform is solved, corrected ones after — so the sticks, this panel's
        Pixel column and click-to-fit read one set of positions.  A re-solve
        changes the correction, which empties the cache on the very next
        ordinary refresh rather than through a hook of its own.
        """

        correction = self._correction_key()
        if correction != self._drawn_correction:
            self._catalog_cache.clear()
            self._drawn_correction = correction

        reference = self._display_reference()
        key = (
            int(order_idx),
            "" if reference is None else reference.lamp,
            "" if reference is None else reference.state.value,
            0 if reference is None else len(reference.lines),
        )
        rows = self._catalog_cache.get(key)
        if rows is None:
            # The fallback is read only when no reference set exists, so the
            # whole table is corrected only when it is going to be drawn.
            fallback = (
                ()
                if reference is not None
                else self.session.corrected_rows(self.session.lines)
            )
            rows = expected_lines_for_order(
                self._shown_reference(reference),
                order_idx,
                fallback_lines=fallback,
            )
            self._catalog_cache[key] = rows
        return rows

    def _shown_reference(self, reference: LampReferenceSet | None):
        """The display reference with its rows moved to where they are seen."""

        if reference is None or self.session.transform is None:
            return reference
        return replace(reference, lines=self.session.corrected_rows(reference.lines))

    def _pooled_marker(self, pool: list, index: int, pen, color: str):
        """Reuse the index-th stick and label of *pool*, creating it once."""

        while len(pool) <= index:
            marker = pg.InfiniteLine(0.0, angle=90, movable=False, pen=pen)
            label = pg.TextItem("", color=color, anchor=(0.5, 0.5))
            self.order_plot.addItem(marker, ignoreBounds=True)
            self.order_plot.addItem(label, ignoreBounds=True)
            pool.append((marker, label))
        marker, label = pool[index]
        marker.setVisible(True)
        label.setVisible(True)
        return marker, label

    @staticmethod
    def _hide_pooled(pool: list, used: int) -> None:
        for marker, label in pool[used:]:
            marker.setVisible(False)
            label.setVisible(False)

    def _refresh_line_help_table(self, rows) -> None:
        """Fill the expected-lines panel from the list the sticks were drawn from.

        The count in this header, the rows below it, and the labelled sticks on
        the spectrum are three renderings of one tuple, so "0 expected Ne lines
        in this order" can no longer sit under three labelled Ne sticks.
        """

        anchored = {
            round(anchor.line.wavelength_nm, 3)
            for anchor in self.session.anchor_rows()
            if anchor.line.order_idx == self.session.selected_order
        }
        self.line_help_table.setRowCount(len(rows))
        for index, row in enumerate(rows):
            intensity = row.relative_intensity
            is_anchored = any(
                abs(row.wavelength_nm - value) <= 0.05 for value in anchored
            )
            values = (
                row.label,
                f"{row.wavelength_nm:.4f}",
                f"{row.detector_pixel:.1f}",
                # Lamp-context strength spans five decades, so two fixed
                # decimals print every ionized row as a flat "0.00".
                "—" if intensity is None else f"{intensity:.3g}",
                "✓" if is_anchored else "",
            )
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                item.setToolTip(one_line(f"{row.label} from {row.source}"))
                if is_anchored:
                    item.setForeground(QtGui.QColor("#6ee7b7"))
                self.line_help_table.setItem(index, column, item)
        reference = self._display_reference()
        lamp = "" if reference is None else reference.lamp
        assigned = (
            "" if self.campaign is None else self.campaign.lamp_for_frame(self.session)
        )
        order_idx = self.session.selected_order
        scope = f"{lamp} " if lamp else ""
        header = (
            f"{len(rows)} expected {scope}line(s) in order {order_idx} · "
            f"{len(anchored)} anchored here"
        )
        if assigned and lamp and catalog_mismatch_warning(lamp, assigned):
            header += f" · showing {lamp}, the assigned lamp is {assigned}"
        elif not rows:
            # "Empty" has several causes and only some are the operator's to
            # fix; saying which is the whole point of the panel.
            if reference is None:
                header += " — assign a lamp role and this fills itself"
            elif not reference.is_referenceable:
                header += f" — {reference.message}"
            else:
                header += f" — no {lamp} lines fall in this order"
                best = reference.best_order
                if best is not None and best != order_idx:
                    header += f"; {lamp} is strongest in order {best}"
        self.line_panel_header.setText(header)

    def _refresh_pattern_traces(self) -> None:
        """Build one trace per order once; an order change only moves the pen."""

        key = (self.session.pattern.shape[0], self.session.pattern.shape[1])
        if key == self._pattern_key and self._pattern_items:
            return
        while self._pattern_items:
            self.detector_plot.removeItem(self._pattern_items.pop())
        columns = np.arange(self.session.pattern.shape[0])
        for order_idx in range(self.session.pattern.shape[1]):
            item = self.detector_plot.plot(
                columns,
                self.session.pattern[:, order_idx],
                pen=pg.mkPen("#7b91a4", width=0.7),
                antialias=False,
            )
            item.setZValue(10)
            self._pattern_items.append(item)
        self._pattern_key = key
        self._selected_trace = None

    def _select_pattern_trace(self, order_idx: int) -> None:
        """Restyle only the two traces whose selection actually changed."""

        if order_idx == self._selected_trace:
            return
        for index in (self._selected_trace, order_idx):
            if index is None or not (0 <= index < len(self._pattern_items)):
                continue
            selected = index == order_idx
            item = self._pattern_items[index]
            item.setPen(
                pg.mkPen(
                    "#6ee7b7" if selected else "#7b91a4",
                    width=2.2 if selected else 0.7,
                )
            )
            item.setZValue(12 if selected else 10)
        self._selected_trace = order_idx

    #: The residual strip's two colours: every bar, and the one whose anchor is
    #: selected in the table above it.
    _RESIDUAL_BAR = "#49b5df"
    _RESIDUAL_BAR_SELECTED = "#ffe08a"

    def _refresh_residual_plot(self) -> None:
        """One bar per accepted anchor, in the anchor table's own order.

        Driven by exactly the triggers it always was — an anchor added, removed
        or cleared, and the transform re-solved — all of which arrive through
        the same ``refresh`` this hangs off.  Nothing here watches a view range.
        """

        self.residual_plot.clear()
        self._residual_bars = None
        self.residual_plot.addLine(
            y=0,
            pen=pg.mkPen("#64748b", style=QtCore.Qt.DashLine),
        )
        residuals = self.session.residuals
        self._label_residual_axis()
        if not residuals:
            return
        # ``Residual`` and ``anchor_rows`` are both sorted by anchor key, so bar
        # *i* is row *i*; the click still matches on the key rather than trusting
        # that, because the key is what actually identifies an anchor.
        x = np.arange(len(residuals), dtype=float)
        y = np.array([residual.dx_px for residual in residuals])
        self._residual_bars = pg.BarGraphItem(
            x=x, height=y, width=0.68, brush=self._RESIDUAL_BAR
        )
        self.residual_plot.addItem(self._residual_bars)
        self._highlight_selected_residual()

    def residual_tick_labels(self) -> list[str]:
        """The wavelengths currently drawn under the bars, in bar order."""

        levels = self.residual_plot.getAxis("bottom")._tickLevels or [[]]
        return [text for _value, text in levels[0]]

    def _residual_label_budget(self, labels) -> int:
        """How many of these wavelengths this width can carry, clear of each other.

        Measured off the widest label the data actually has rather than a
        template, because a four-figure wavelength is a wider label than a
        three-figure one and the count has to be honest about the labels being
        drawn, not about a typical one.
        """

        width = self.residual_plot.getViewBox().width()
        if width <= 0:
            # Before the first layout pass the strip has no width to divide.
            # The rail's own floor is the honest stand-in, and the show pass
            # comes straight back through here with the real one.
            width = self._tables_minimum_width()
        metrics = QtGui.QFontMetrics(self._residual_tick_font)
        widest = max(metrics.horizontalAdvance(text) for text in labels)
        return max(1, int(width // (widest + BENCH_RESIDUAL_TICK_GAP)))

    def _label_residual_axis(self) -> None:
        """Wavelengths under the bars, thinned to what the rail can show.

        Sparse anchors are all labelled.  Dense ones are sampled evenly across
        the whole set, ends included — the same thinning the spectrum overlays
        do when a family has more lines than a view can label legibly — so the
        axis always reads as a wavelength axis rather than as a smear, and the
        two ends of the anchor set are always named.

        Geometry, not data: this is called both when the anchors change and
        when the rail is resized or dragged, because how many labels fit is a
        question about the width, and the answer changes with it.
        """

        axis = self.residual_plot.getAxis("bottom")
        residuals = self.session.residuals
        if not residuals:
            axis.setTicks([[], []])
            return
        labels = [f"{residual.wavelength_nm:.1f}" for residual in residuals]
        metrics = QtGui.QFontMetrics(self._residual_tick_font)
        budget = self._residual_label_budget(labels)
        last = len(labels) - 1
        if len(labels) <= budget:
            shown = range(len(labels))
        elif budget <= 1:
            shown = (0,)
        else:
            shown = sorted(
                {round(index * last / (budget - 1)) for index in range(budget)}
            )
        shown = list(shown)
        axis.setTicks([[(float(index), labels[index]) for index in shown], []])
        # The first and last wavelengths are centred on the first and last bars,
        # so half of each hangs outside the bars' own span: without a margin cut
        # to their width they are drawn off the end of the strip and clipped.
        # This is what makes the ends legible, which is the whole reason the
        # ends are the labels thinning never drops.
        box = self.residual_plot.getViewBox()
        width = box.width()
        outer = (
            max(
                metrics.horizontalAdvance(labels[shown[0]]),
                metrics.horizontalAdvance(labels[shown[-1]]),
            )
            / 2.0
        )
        margin = 0.5
        if width > 2 * outer:
            margin += outer * max(1.0, float(last)) / (width - 2 * outer)
        box.setXRange(-margin, last + margin, padding=0.0)

    def _selected_anchor_key(self):
        """The key of the anchor whose row is selected, if a row is."""

        table = self.anchor_table
        model = table.selectionModel()
        if model is None or not model.hasSelection():
            return None
        anchors = self.session.anchor_rows()
        row = table.currentRow()
        return anchors[row].key if 0 <= row < len(anchors) else None

    def _highlight_selected_residual(self) -> None:
        """Draw the selected anchor's own bar in the selection colour.

        The cheap half of the two-way sync: one ``setOpts`` on the item that is
        already drawn, no rebuild, no extra refresh path.
        """

        bars = self._residual_bars
        if bars is None:
            return
        selected = self._selected_anchor_key()
        bars.setOpts(
            brushes=[
                pg.mkBrush(
                    self._RESIDUAL_BAR_SELECTED
                    if selected is not None and residual.key == selected
                    else self._RESIDUAL_BAR
                )
                for residual in self.session.residuals
            ]
        )

    def _residual_plot_clicked(self, event) -> None:
        """Click a bar to select its anchor's row in the table above it.

        This is what replaces the per-anchor tick labels: the rail is nowhere
        near wide enough for one wavelength under every bar — they overlapped
        even across the old full-width view — so the identity of a bar is asked
        for rather than always printed, and the answer arrives in the row, in
        the readings it carries, and in the Why dock the row already fills.
        """

        if event.button() != QtCore.Qt.LeftButton:
            return
        residuals = self.session.residuals
        if not residuals:
            return
        if not self.residual_plot.sceneBoundingRect().contains(event.scenePos()):
            return
        point = self.residual_plot.getViewBox().mapSceneToView(event.scenePos())
        index = int(round(point.x()))
        # Half a slot either side: the gap between bars belongs to the nearer
        # bar, so a click never has to land on a one-pixel-wide sliver.
        if not (0 <= index < len(residuals)) or abs(point.x() - index) > 0.5:
            return
        key = residuals[index].key
        for row, anchor in enumerate(self.session.anchor_rows()):
            if anchor.key == key:
                self.anchor_table.selectRow(row)
                item = self.anchor_table.item(row, 0)
                if item is not None:
                    self.anchor_table.scrollToItem(item)
                return


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-calib",
        description=(
            "Triage Andor SIF exposures, assign measurement roles by hand, and "
            "interactively fit a rigid wavelength-calibration alignment."
        ),
    )
    parser.add_argument(
        "folder",
        nargs="?",
        type=Path,
        default=Path.cwd(),
        help="folder the Add files dialog opens in (default: current folder)",
    )
    parser.add_argument(
        "--file",
        type=Path,
        action="append",
        help="load one SIF at start-up; repeat for several",
    )
    parser.add_argument(
        "--watch",
        action="store_true",
        help="also poll the folder for new stable SIFs (optional convenience)",
    )
    parser.add_argument(
        "--pattern",
        type=Path,
        default=_DEFAULT_PATTERN,
        help=(
            "order-pattern table: a path, or the bare filename of a packaged "
            "one such as pattern_CMOS_20240305.txt (default: "
            f"{_DEFAULT_PATTERN.name})"
        ),
    )
    parser.add_argument(
        "--wavelength",
        type=Path,
        default=_DEFAULT_WAVELENGTH,
        help=(
            "wavelength table: a path, or the bare filename of a packaged one "
            f"such as Th_wavelength_CMOS_20240305.txt (default: "
            f"{_DEFAULT_WAVELENGTH.name})"
        ),
    )
    parser.add_argument("--integral", type=Path, default=_DEFAULT_INTEGRAL)
    parser.add_argument(
        "--previous-sphere",
        type=Path,
        default=_DEFAULT_PREVIOUS_SPHERE,
        help=(
            "sphere signal of the campaign to compare against (default: the "
            f"packaged 2024-03-05 pair, {_DEFAULT_PREVIOUS_SPHERE.name}); "
            "point this elsewhere when calibrating the folder those packaged "
            "frames were copied from"
        ),
    )
    parser.add_argument(
        "--previous-sphere-background",
        type=Path,
        default=_DEFAULT_PREVIOUS_SPHERE_BACKGROUND,
        help=(
            "background of that previous sphere (default: the packaged "
            f"2024-03-05 pair, {_DEFAULT_PREVIOUS_SPHERE_BACKGROUND.name})"
        ),
    )
    parser.add_argument(
        "--lamp",
        action="append",
        metavar="NAME",
        help=(
            "lamp to suggest from the previous campaign; any name is accepted "
            "and none is ever required (default: "
            f"{', '.join(PREVIOUS_CAMPAIGN_LAMPS)})"
        ),
    )
    # Unset by default, like the two roots below and for the same kind of
    # reason: the identity is dated by the day the frames were TAKEN, and that
    # day is not knowable until the folder argument — and then the frames
    # themselves — have been read.  Today is the last resort, not the default.
    parser.add_argument(
        "--snapshot-id",
        default=None,
        help=(
            "planned snapshot identity (default: YYYYMMDD_<detector> dated by "
            "the acquisition — the date leading the folder argument's name, "
            "else the date in the loaded frames' SIF headers, else today)"
        ),
    )
    parser.add_argument("--detector", default="cmos")
    parser.add_argument("--base-snapshot", default="20250926_cmos")
    parser.add_argument(
        "--valid-from",
        type=date.fromisoformat,
        default=None,
        help=(
            "ISO date the saved snapshot's open-ended calibration epoch starts "
            "(default: the acquisition date the snapshot identity is named for)"
        ),
    )
    # Left unset on purpose: the default is not knowable until the folder
    # argument has been parsed, and a default computed from the working
    # directory is precisely the trap this replaces.
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "snapshot parent directory "
            f"(default: {SNAPSHOT_ROOT_NAME}/ inside the folder argument)"
        ),
    )
    parser.add_argument(
        "--config-root",
        type=Path,
        default=None,
        help=(
            "parent for generated commented TOML bundles "
            f"(default: {SNAPSHOT_ROOT_NAME}/{CONFIG_ROOT_NAME}/ inside the "
            "folder argument, so one generated folder holds everything)"
        ),
    )
    parser.add_argument("--poll-ms", type=int, default=1000)
    parser.add_argument("--stable-polls", type=int, default=2)
    parser.add_argument("--minimum-age-s", type=float, default=1.0)
    parser.add_argument("--saturation-level", type=float, default=0.98 * 65535)
    parser.add_argument("--minimum-snr", type=float, default=5.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the calibration bench."""

    # Before any window exists, and before Qt asks the shell for one: an
    # explicit identity is what makes Windows show *this* icon on the taskbar
    # button rather than the launcher stub's.
    apply_windows_taskbar_identity()
    args = _build_parser().parse_args(argv)
    if not args.folder.is_dir():
        raise SystemExit(f"folder not found: {args.folder}")
    # Resolved, then said out loud: which pattern the bench opened on is the one
    # fact the fit cannot tell you afterwards, and a bare packaged name resolves
    # to a file the operator never typed the path of.
    args.pattern = resolve_calibration_file(args.pattern, "pattern")
    args.wavelength = resolve_calibration_file(args.wavelength, "wavelength")
    print(f"Pattern:    {args.pattern}")
    print(f"Wavelength: {args.wavelength}")
    if not args.integral.is_file():
        raise SystemExit(f"integrating-sphere reference not found: {args.integral}")
    if args.poll_ms < 50:
        raise SystemExit("--poll-ms must be at least 50")

    # The folder argument is the campaign's own folder, so it — not whatever
    # directory the shortcut started in — is what the roots hang off.  An
    # explicit flag still wins, and is made absolute so the Save tab can be
    # honest about it too.
    default_output, default_config = default_bench_roots(args.folder)
    output_root = default_output if args.output_root is None else absolute_root(args.output_root)
    config_root = default_config if args.config_root is None else absolute_root(args.config_root)

    # The calibration is dated by the day its images were taken.  The folder
    # the bench was launched at is the first place that says so — the owner's
    # ``20250926_calib`` — and the frames' own headers are the second, once
    # they load.  Today's date is only ever the placeholder for a day nothing
    # could tell us, and the Save tab says as much when it is used.
    snapshot_date = acquisition_date_from_name(args.folder)
    snapshot_id = args.snapshot_id or snapshot_identity(
        snapshot_date or date.today(), args.detector
    )

    pattern = np.loadtxt(args.pattern, dtype=int)
    lines = load_wavelength_table(args.wavelength)
    session = CalibrationBenchSession(
        pattern,
        lines,
        saturation_level=args.saturation_level,
        minimum_snr=args.minimum_snr,
        # Read from the table itself, so an adjusted table inherits the vetting
        # of the curated one it was derived from and a stranger inherits none.
        vetting=table_vetting(args.wavelength, [_CALIBRATION_DIR]),
    )
    campaign = CalibrationCampaignSession(
        pattern_source=args.pattern,
        wavelength_source=args.wavelength,
        integral_source=args.integral,
        suggested_lamps=args.lamp or PREVIOUS_CAMPAIGN_LAMPS,
        previous_sphere=args.previous_sphere,
        previous_sphere_background=args.previous_sphere_background,
    )
    watcher = (
        StableSifWatcher(
            args.folder,
            required_unchanged_polls=args.stable_polls,
            minimum_age_s=args.minimum_age_s,
        )
        if args.watch
        else None
    )
    loader = FrameLoader(pattern)
    # Antialiasing a 42k-sample factor curve costs more than it shows: the
    # paint of one such curve dominated every pan the owner tried.
    pg.setConfigOptions(antialias=False, imageAxisOrder="col-major")
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    application.setApplicationName("Echelle calibration bench")
    application.setWindowIcon(bench_window_icon())
    window = CalibrationBenchWindow(
        session,
        campaign=campaign,
        watcher=watcher,
        loader=loader,
        folder=args.folder,
        output_root=output_root,
        config_root=config_root,
        snapshot_id=snapshot_id,
        snapshot_date=(
            acquisition_date_from_name(args.snapshot_id)
            if args.snapshot_id
            else snapshot_date
        ),
        snapshot_id_explicit=args.snapshot_id is not None,
        detector=args.detector,
        base_snapshot=args.base_snapshot,
        valid_from=args.valid_from,
        poll_interval_ms=args.poll_ms,
    )
    window.show()
    if args.file:
        QtCore.QTimer.singleShot(0, lambda: window.add_paths(args.file))
    return int(application.exec_())


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
