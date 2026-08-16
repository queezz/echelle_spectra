"""Separate pyqtgraph live calibration bench for ``echelle-calib``."""

from __future__ import annotations

import argparse
import html
import platform
import sys
from collections.abc import Sequence
from datetime import date
from pathlib import Path

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
    CalibrationCampaignSession,
    ChecklistState,
    ComparisonState,
    ExposureState,
    ExposureTriage,
    LampReferenceSet,
    MeasurementRole,
    TomlState,
    catalog_mismatch_warning,
    default_validity,
    expected_lines_for_order,
    lamp_reference_set,
    triage_for_role,
)
from .snapshot import SnapshotError
from .tools.calibration_alignment import load_wavelength_table, table_vetting

_PACKAGE_DIR = Path(__file__).parent
_CALIBRATION_DIR = _PACKAGE_DIR / "resources" / "calibration_files"
_DEFAULT_PATTERN = _CALIBRATION_DIR / "pattern_CMOS_20250926.txt"
_DEFAULT_WAVELENGTH = (
    _CALIBRATION_DIR
    / "alignments"
    / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
)
_DEFAULT_INTEGRAL = _CALIBRATION_DIR / "integrating_sphere.txt"
_DEFAULT_PREVIOUS_SPHERE = _CALIBRATION_DIR / "sphere_cmos_20240305.sif"
_DEFAULT_PREVIOUS_SPHERE_BACKGROUND = (
    _CALIBRATION_DIR / "sphere_cmos_20240305_bkg.sif"
)

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


def forget_session_layout() -> None:
    """Drop every remembered splitter cut (used by tests and by a fresh run)."""

    _SESSION_SPLITTER_SHARES.clear()


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
        output_root: str | Path = "calibrations",
        config_root: str | Path = "calibration-configs",
        snapshot_id: str = "",
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
        self.output_root = Path(output_root)
        self.config_root = Path(config_root)
        self.initial_snapshot_id = snapshot_id
        self.initial_detector = detector
        self.initial_base_snapshot = base_snapshot
        self.valid_from = valid_from
        self._load_thread: FrameLoadThread | None = None
        self._background_thread: FrameLoadThread | None = None
        self._campaign_thread: CampaignTaskThread | None = None
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
        self._catalog_rows: tuple = ()
        self._queue: list[Path] = []
        self._file_rows: list[Path] = []
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
        self.last_folder = Path(watcher.folder) if watcher is not None else Path.cwd()
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
            #dropTarget {{ border: 2px dashed #49b5df; border-radius: 9px;
                          color: #8fd9ff; font-weight: 700; letter-spacing: 1px;
                          padding: 12px; }}
            #triageHeadline {{ font-weight: 700; padding: 6px; }}
            #stateBadge {{ color: #7ee2b8; font-weight: 700; }}
            #warningPanel {{ background: #2a1e13; border-left: 3px solid #ffb86b;
                            padding: 7px; color: #ffd6a3; }}
            #messagePanel {{ background: #0f141a; border-left: 3px solid #49b5df;
                            padding: 7px; color: #bed4e1; }}
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
        _emphasise(self.drop_hint, body)
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
                self.anchors_panel, self.anchor_table, self.anchor_buttons
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

    def _table_panel_minimum_height(self, panel, table, header_widget) -> int:
        """The height below which a rail panel stops being a working surface.

        Derived from the table's own metrics, so the floor moves with the
        platform font instead of pinning a row count to one designer's display.
        """

        layout = panel.layout()
        margins = layout.contentsMargins()
        return int(
            self.TABLE_FLOOR_ROWS * table.verticalHeader().defaultSectionSize()
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
            self._apply_splitter_shares(self.tables_splitter, "tables", _TABLES_SHARES)
            self._watch_splitter(self.readings_splitter, "readings")
            self._pin_status_band_height()
        finally:
            self._distributing = False

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
            lambda *_args, _s=splitter, _k=key: self._remember_cut(_s, _k)
        )

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
        """The anchor table, down the right rail beside the lines it anchors.

        Only the table and its two buttons: the numbers the fit produces are
        readings and live on the readings strip, where they are in view from
        every tab rather than behind whichever one happens to be open.
        """

        panel = QtWidgets.QGroupBox("Anchors")
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 4)
        layout.setSpacing(5)

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
        self.details_view.setMinimumHeight(3 * self.layout_unit)
        self.details_view.setPlaceholderText(
            "Click any verdict, checklist row, file, or anchor and the whole "
            "explanation is written here. It changes only when asked."
        )
        dock.setWidget(self.details_view)
        self.details_dock = dock
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        self.resizeDocks([dock], [4 * self.layout_unit], QtCore.Qt.Vertical)
        self.explain(
            "Exposure triage is the front door",
            "Drop any SIF and the bench judges the exposure before any role "
            "exists: clustered full-scale pixels are real saturation, isolated "
            "ones are cosmic rays, and the histogram's top end is the number "
            "you adjust the lamp by.",
        )

    def explain(self, title: str, text: str) -> None:
        """Show one full explanation in the details dock."""

        body = str(text).strip().replace("\n", "<br>")
        self.details_view.setHtml(
            f"<h3 style='color:#8fd9ff;margin:0 0 6px 0'>{title}</h3>"
            f"<div style='color:#cfe0ec;line-height:145%'>{body}</div>"
        )

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

        self.drop_hint = QtWidgets.QLabel(
            "DROP SIF FILES HERE\nany names, any order, as many as you like"
        )
        self.drop_hint.setAlignment(QtCore.Qt.AlignCenter)
        self.drop_hint.setWordWrap(True)
        self.drop_hint.setObjectName("dropTarget")
        self.drop_hint.setMinimumHeight(92)
        layout.addWidget(self.drop_hint)

        button_row = QtWidgets.QHBoxLayout()
        self.add_files_button = QtWidgets.QPushButton("Add SIF files…")
        self.remove_file_button = QtWidgets.QPushButton("Remove selected")
        button_row.addWidget(self.add_files_button, 2)
        button_row.addWidget(self.remove_file_button, 1)
        layout.addLayout(button_row)

        self.confirm_roles_button = QtWidgets.QPushButton("Confirm suggested roles")
        self._explainable(
            self.confirm_roles_button,
            "Suggested is not assigned",
            "A filename may pre-fill the Role control, but the bench is not "
            "given that role until somebody says so. Rows in that state are "
            "marked SUGGESTED. This button assigns every one of them at once, "
            "so a whole acquisition folder is one deliberate press rather than "
            "one popup pick per row.",
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
        self.watch_value = _ElidingLabel("manual — drag and drop or Add files")
        self.file_value = _ElidingLabel("no file open")
        self.file_state_value = QtWidgets.QLabel("WAITING")
        self.file_state_value.setObjectName("stateBadge")
        status_form.addRow("Input", self.watch_value)
        status_form.addRow("Open frame", self.file_value)
        status_form.addRow("File state", self.file_state_value)
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
            "Absolute calibration factors",
            "The sphere signal minus its background, divided by the "
            "integrating sphere's known radiance, gives the factor curve that "
            "turns counts into W m⁻² sr⁻¹ nm⁻¹. The median ratio compares this "
            "campaign's curve with the previous one: near 1 means the "
            "instrument's response has not moved, and a large departure is "
            "either real ageing or an exposure-normalisation mismatch worth "
            "chasing before the trip. Only the sphere pair is needed — no lamp.",
        )
        self.compare_button = QtWidgets.QPushButton("Compute factors")
        comparison_layout.addWidget(self.comparison_value)
        comparison_layout.addWidget(self.compare_button)
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
        self.transform_value = QtWidgets.QLabel("—")
        self.transform_value.setWordWrap(True)
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

    def _build_lamp_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        order_group = QtWidgets.QGroupBox("Order and frame")
        order_layout = QtWidgets.QVBoxLayout(order_group)
        order_row = QtWidgets.QHBoxLayout()
        order_row.addWidget(QtWidgets.QLabel("Order"))
        self.previous_order_button = QtWidgets.QToolButton()
        self.previous_order_button.setText("◀")
        self.previous_order_button.setToolTip("Previous Echelle order")
        self.order_spin = QtWidgets.QSpinBox()
        self.order_spin.setRange(0, self.session.pattern.shape[1] - 1)
        self.order_spin.setToolTip(
            "Which Echelle order the spectrum plot shows and a click fits."
        )
        self.next_order_button = QtWidgets.QToolButton()
        self.next_order_button.setText("▶")
        self.next_order_button.setToolTip("Next Echelle order")
        self.order_total_value = QtWidgets.QLabel(
            f"of {self.session.pattern.shape[1] - 1}"
        )
        order_row.addWidget(self.previous_order_button)
        order_row.addWidget(self.order_spin, 1)
        order_row.addWidget(self.next_order_button)
        order_row.addWidget(self.order_total_value)
        order_layout.addLayout(order_row)

        frame_row = QtWidgets.QHBoxLayout()
        frame_row.addWidget(QtWidgets.QLabel("Fit on"))
        self.frame_combo = QtWidgets.QComboBox()
        self.frame_combo.addItem("Mean of all frames", None)
        frame_row.addWidget(self.frame_combo, 1)
        order_layout.addLayout(frame_row)
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
        order_layout.addWidget(self.frame_choice_value)

        family_row = QtWidgets.QHBoxLayout()
        family_row.addWidget(QtWidgets.QLabel("Expected lines"))
        self.line_family_combo = QtWidgets.QComboBox()
        self.line_family_combo.addItems(list(KNOWN_LAMP_NAMES))
        family_row.addWidget(self.line_family_combo, 1)
        order_layout.addLayout(family_row)
        self._explainable(
            self.line_family_combo,
            "Which catalog fills the expected-lines table",
            "This follows the assigned lamp on its own: assign a Ne lamp and "
            "the expected-lines panel fills with neon. Changing it by hand is "
            "an override for comparing one lamp's catalog against another "
            "frame; the anchors themselves always reference the assigned "
            "lamp's own rows, and the panel warns when the two disagree.",
        )
        layout.addWidget(order_group)

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

        destination = QtWidgets.QLabel(
            f"Snapshots: {self.output_root.name}\nConfigs: {self.config_root.name}"
        )
        destination.setObjectName("mutedText")
        destination.setWordWrap(True)
        layout.addWidget(destination)
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
        layout.addWidget(self.save_state_value)
        self.toml_preview = QtWidgets.QPlainTextEdit()
        self.toml_preview.setReadOnly(True)
        self.toml_preview.setPlaceholderText(
            "Generated campaign.toml appears here; all files remain ordinary and editable."
        )
        self._fills_its_share(self.toml_preview)
        layout.addWidget(self.toml_preview, 1)
        self.control_tabs.addTab(self._scrollable(tab), "Save")

    def _build_triage_view(self) -> None:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
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

        graphics = pg.GraphicsLayoutWidget()
        graphics.setBackground("#10151b")
        self.histogram_plot = graphics.addPlot(row=0, col=0, title="Raw counts histogram")
        self.histogram_plot.setLabel("bottom", "counts")
        self.histogram_plot.setLabel("left", "pixels per bin")
        self.histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.histogram_plot.setLogMode(y=True)
        layout.addWidget(graphics, 1)

        # The top end lives in its own widget so that the honest answer — "no
        # pixels are up here" — can take its place in words (F17 item 3). An
        # empty log histogram drew a solid block, which says nothing at all.
        self.top_histogram_widget = pg.GraphicsLayoutWidget()
        self.top_histogram_widget.setBackground("#10151b")
        self.top_histogram_plot = self.top_histogram_widget.addPlot(
            row=0, col=0, title="Top end — the last 10% before full scale"
        )
        self.top_histogram_plot.setLabel("bottom", "counts")
        self.top_histogram_plot.setLabel("left", "pixels per bin")
        self.top_histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.top_histogram_plot.setLogMode(y=True)
        layout.addWidget(self.top_histogram_widget, 1)
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

        self.equal_aspect_check = QtWidgets.QCheckBox("Square detector pixels")
        self._explainable(
            self.equal_aspect_check,
            "Showing the detector at its true shape",
            "The detector view stretches to fill the space it is given, which "
            "is what makes the order traces far enough apart to aim at. Tick "
            "this and one detector pixel is drawn square in both directions, "
            "so the frame appears at its real 2560x2160 proportions. Useful "
            "for judging the geometry; unhelpful for clicking lines, which is "
            "why it is off unless you ask.",
        )
        outer.addWidget(self.equal_aspect_check)

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

        self.order_plot = graphics.addPlot(row=1, col=0, title="Selected order spectrum")
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

        self.residual_plot = graphics.addPlot(row=2, col=0, title="Anchor residuals")
        self.residual_plot.setMaximumHeight(180)
        self.residual_plot.setLabel("bottom", "accepted anchor")
        self.residual_plot.setLabel("left", "fit residual", units="px")
        self.residual_plot.getAxis("bottom").setHeight(58)
        self.residual_plot.getAxis("left").enableAutoSIPrefix(False)
        self.residual_plot.addLine(y=0, pen=pg.mkPen("#64748b", style=QtCore.Qt.DashLine))
        split.addWidget(graphics)
        split.setStretchFactor(0, 1)
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
        layout.addWidget(self.sphere_view_message)
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
        self.add_files_button.clicked.connect(self._pick_files)
        self.remove_file_button.clicked.connect(self._remove_selected_file)
        self.confirm_roles_button.clicked.connect(self._confirm_suggested_roles)
        self.show_frame_button.clicked.connect(self._open_selected_file)
        self.file_table.itemSelectionChanged.connect(self._file_selection_changed)
        self.file_table.itemClicked.connect(self._file_row_clicked)
        self.checklist_tree.currentRowChanged.connect(self._checklist_row_selected)
        self.anchor_table.itemSelectionChanged.connect(self._anchor_row_selected)
        self.line_help_table.itemSelectionChanged.connect(self._expected_line_selected)
        self.compare_button.clicked.connect(self._start_sphere_comparison)
        self.generate_tomls_button.clicked.connect(lambda: self._generate_tomls())
        self.regenerate_tomls_button.clicked.connect(self._regenerate_tomls)
        self.save_snapshot_button.clicked.connect(self._save_snapshot)
        self.snapshot_id_edit.textChanged.connect(self.refresh_campaign)
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
        # operator's own action is what they need to read, so it is set last.
        self.refresh()
        self.message_value.setText(outcome)

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
            f"{entry.source}. The pixel is where the current wavelength table "
            "puts the line, so it says where to look, not where the line is: "
            "click the peak itself to fit its centroid and accept it as an "
            "anchor, and click an anchored stick again to take it back off.",
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
            self.message_value.setText(
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
        self.message_value.setText(
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
            if self.campaign.remove_classification(path):
                self.campaign.scope_alignment_to_lamp(self.session)
                self.message_value.setText(f"{path.name} is unassigned again.")
                self.refresh()
            return
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

    def _start_sphere_comparison(self) -> None:
        if self.campaign is None:
            return
        self.comparison_value.setText("COMPUTING — using the established absolute engine…")
        self._start_campaign_task(self.campaign.compute_sphere_comparison)

    @QtCore.pyqtSlot(object)
    def _campaign_task_completed(self, result) -> None:
        state = getattr(result, "state", None)
        if state is ComparisonState.READY:
            self.message_value.setText("Sphere factors computed and compared.")
        elif state is ComparisonState.INSUFFICIENT_DATA:
            self.message_value.setText(
                "Candidate factors computed; previous comparison is insufficient data."
            )
        elif hasattr(result, "snapshot_id"):
            correction = getattr(self.campaign, "wavelength_correction", None)
            detail = "" if correction is None else f" Saved wavelength.txt: {correction.reason}."
            self.message_value.setText(
                f"Snapshot {result.snapshot_id} saved and validated through Packet 0.{detail}"
            )
        self.refresh()

    @QtCore.pyqtSlot(str)
    def _campaign_task_failed(self, reason: str) -> None:
        self.message_value.setText(f"Campaign action failed safely: {reason}")
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
        try:
            paths = self.campaign.write_tomls(
                self.config_root,
                self.snapshot_id_edit.text().strip(),
                self.session,
                overwrite=overwrite,
            )
        except (OSError, SnapshotError, ValueError) as exc:
            self.message_value.setText(f"Alignment settings were not saved: {exc}")
            # An existing bundle is a choice to offer, not a dead end: the
            # button that can get past it appears beside the one that refused.
            self.regenerate_tomls_button.setVisible("already exists" in str(exc))
        else:
            self.regenerate_tomls_button.setVisible(False)
            self.message_value.setText(
                "Saved the alignment settings, the campaign, and the export "
                f"configuration as commented files you can edit: {paths['campaign'].parent}."
                if not overwrite
                else "Rewrote the alignment settings, the campaign, and the export "
                f"configuration, replacing what was in {paths['campaign'].parent}."
            )
            self.toml_preview.setPlainText(
                paths["campaign"].read_text(encoding="utf-8")
            )
        self.refresh_campaign()

    def _regenerate_tomls(self) -> None:
        """Rewrite settings files that already exist, on a deliberate press."""

        self._generate_tomls(overwrite=True)

    def _save_snapshot(self) -> None:
        if self.campaign is None:
            return
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
            self.message_value.setText(record.triage.headline)
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

    def _order_changed(self, order_idx: int) -> None:
        self.session.set_selected_order(order_idx)
        self._refresh_frame_selector()
        self.refresh_plots()

    def anchor_near(self, order_idx: int, pixel: float):
        """The anchor of *order_idx* whose stick a click at *pixel* lands on.

        Either end counts: the curated row's expected pixel, where the stick is
        drawn, and the fitted centroid, where the green marker sits.
        """

        best = None
        for anchor in self.session.anchor_rows():
            if anchor.line.order_idx != int(order_idx):
                continue
            tolerance = max(3.0, anchor.line.width_px)
            distance = min(
                abs(pixel - anchor.line.center_pixel),
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
        """Two lines: which file, and what it says. Nothing else shouts.

        The full sentence, its advice, and the four-line breakdown are all in
        the Why dock. A headline that runs to fourteen lines is a paragraph
        wearing a headline's font, and it is exactly what "all loud and big"
        produced.
        """

        readings = []
        if triage.state is ExposureState.SATURATED:
            clusters = triage.saturation.cluster_count
            readings.append(f"{clusters} saturated cluster(s)")
        if triage.headroom_fraction is not None:
            readings.append(f"peak {100.0 * triage.headroom_fraction:.0f}% of full scale")
        if triage.saturation.anomalous_pixels:
            readings.append(f"{triage.saturation.anomalous_pixels} anomalies")
        tail = f" · {' · '.join(readings)}" if readings else ""
        return f"{path.name}\n{verdict.label.upper()}{tail}"

    def _show_triage(self, path: Path) -> None:
        """Render one file's exposure verdict; roles play no part in it."""

        record = None if self.campaign is None else self.campaign.loaded.get(path)
        if record is None:
            return
        assert self.campaign is not None
        triage = record.triage
        measurement = self.campaign.measurements.get(path)
        verdict = triage_for_role(
            triage, measurement.role if measurement is not None else None
        )
        color = _TRIAGE_COLORS[triage.state]
        headline = triage.headline
        if triage.state is ExposureState.SATURATED and not verdict.blocking:
            # The bright/dim pair exists so the dim series saturates its strong
            # lines; saying FAILED here would be the bench misreading physics.
            color = _TRIAGE_COLORS[ExposureState.DIM]
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
        self.triage_next_value.setText(guidance.next_action)
        self._explainable(
            self.triage_next_value,
            f"{path.name} — how the verdict was reached",
            breakdown,
        )
        # Selecting a file writes its whole verdict into the dock at once:
        # what it says, what to do about it, and the numbers behind both.
        self.explain(f"{path.name} — exposure verdict", whole_verdict)
        self.exposure_value.setText(
            f"{path.name} · {verdict.label.upper()} · peak {peak}. "
            f"{guidance.next_action}"
        )
        self._explainable(
            self.exposure_value,
            f"{path.name} — next acquisition action",
            f"{guidance.next_action}\n\n{verdict.advice}",
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

    def refresh(self) -> None:
        """Render the current domain state without changing it."""

        if self.watcher is not None:
            self.watch_value.setText(
                f"drag and drop, Add files, or watching {self.watcher.folder.name or '.'}"
            )
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
        self._refresh_file_buttons()
        self.drop_hint.setVisible(not self._file_rows)
        self.compare_button.setEnabled(enabled and not busy)
        self.generate_tomls_button.setEnabled(enabled and not busy)
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
                verdict = triage_for_role(
                    record.triage,
                    measurement.role if measurement is not None else None,
                )
                marks.append(verdict.label.upper())
                if record.triage.saturation.anomalous_pixels:
                    marks.append(f"{record.triage.saturation.anomalous_pixels} anomalies")
            if measurement is None:
                # The badge for an unconfirmed suggestion lives here, in full,
                # never inside the Role control where it would be elided.
                marks.append(
                    f"{_SUGGESTED_BADGE} ONLY — no role assigned"
                    if path in suggested
                    else "no role yet"
                )
            elif measurement.lamp_family:
                marks.append(f"{measurement.lamp_family} {measurement.role.value}")
            else:
                marks.append(measurement.role.value)
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
                    else _TRIAGE_COLORS[record.triage.state]
                )
                if verdict.state is ExposureState.SATURATED and not verdict.blocking:
                    # Expected saturation is a note, not an alarm.
                    colour = _TRIAGE_COLORS[ExposureState.DIM]
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
            self.comparison_value.setText(
                "READY · new/previous median "
                f"{comparison.median_ratio:.3f}; 5–95% "
                f"{comparison.p05_ratio:.3f}–{comparison.p95_ratio:.3f} "
                f"({comparison.sample_count} samples)."
            )
        else:
            self.comparison_value.setText(
                f"{comparison.state.value.replace('-', ' ').upper()} · {comparison.reason}"
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
        for item in self.campaign.checklist(self.session):
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
        self.sphere_view_message.setText(
            f"{comparison.state.value.replace('-', ' ').upper()} — {comparison.reason}"
        )

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

    def _expected_rows_for_order(self, order_idx: int) -> tuple:
        """Cache the expected-line list per (order, reference) pair.

        Rebuilding it on every arrow-key press is work the answer does not
        change with, and F16's order-scrolling budget has no room for it.
        """

        reference = self._display_reference()
        key = (
            int(order_idx),
            "" if reference is None else reference.lamp,
            "" if reference is None else reference.state.value,
            0 if reference is None else len(reference.lines),
        )
        rows = self._catalog_cache.get(key)
        if rows is None:
            rows = expected_lines_for_order(
                reference, order_idx, fallback_lines=self.session.lines
            )
            self._catalog_cache[key] = rows
        return rows

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
                "—" if intensity is None else f"{intensity:.2f}",
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

    def _refresh_residual_plot(self) -> None:
        self.residual_plot.clear()
        self.residual_plot.addLine(
            y=0,
            pen=pg.mkPen("#64748b", style=QtCore.Qt.DashLine),
        )
        if not self.session.residuals:
            return
        x = np.arange(len(self.session.residuals), dtype=float)
        y = np.array([residual.dx_px for residual in self.session.residuals])
        bars = pg.BarGraphItem(x=x, height=y, width=0.68, brush="#49b5df")
        self.residual_plot.addItem(bars)
        labels = [
            f"o{residual.order_idx}\n{residual.wavelength_nm:.2f}"
            for residual in self.session.residuals
        ]
        axis = self.residual_plot.getAxis("bottom")
        axis.setTicks([list(zip(x, labels))])


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
    parser.add_argument("--pattern", type=Path, default=_DEFAULT_PATTERN)
    parser.add_argument("--wavelength", type=Path, default=_DEFAULT_WAVELENGTH)
    parser.add_argument("--integral", type=Path, default=_DEFAULT_INTEGRAL)
    parser.add_argument("--previous-sphere", type=Path, default=_DEFAULT_PREVIOUS_SPHERE)
    parser.add_argument(
        "--previous-sphere-background",
        type=Path,
        default=_DEFAULT_PREVIOUS_SPHERE_BACKGROUND,
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
    parser.add_argument(
        "--snapshot-id",
        default=f"{date.today():%Y%m%d}_cmos",
        help="planned snapshot identity (default: today's CMOS identity)",
    )
    parser.add_argument("--detector", default="cmos")
    parser.add_argument("--base-snapshot", default="20250926_cmos")
    parser.add_argument(
        "--valid-from",
        type=date.fromisoformat,
        default=date.today(),
        help=(
            "ISO date the saved snapshot's open-ended calibration epoch starts "
            "(default: today)"
        ),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path.cwd() / "calibrations",
        help="snapshot parent directory",
    )
    parser.add_argument(
        "--config-root",
        type=Path,
        default=Path.cwd() / "calibration-configs",
        help="parent for generated commented TOML bundles",
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
    if not args.pattern.is_file():
        raise SystemExit(f"pattern file not found: {args.pattern}")
    if not args.wavelength.is_file():
        raise SystemExit(f"wavelength table not found: {args.wavelength}")
    if not args.integral.is_file():
        raise SystemExit(f"integrating-sphere reference not found: {args.integral}")
    if args.poll_ms < 50:
        raise SystemExit("--poll-ms must be at least 50")

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
        output_root=args.output_root,
        config_root=args.config_root,
        snapshot_id=args.snapshot_id,
        detector=args.detector,
        base_snapshot=args.base_snapshot,
        valid_from=args.valid_from,
        poll_interval_ms=args.poll_ms,
    )
    window.last_folder = args.folder
    window.show()
    if args.file:
        QtCore.QTimer.singleShot(0, lambda: window.add_paths(args.file))
    return int(application.exec_())


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
