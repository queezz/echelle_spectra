"""Separate pyqtgraph live calibration bench for ``echelle-calib``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets

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
    MeasurementRole,
    TomlState,
    catalog_lines_for_order,
    catalog_mismatch_warning,
    default_validity,
    triage_for_role,
)
from .snapshot import SnapshotError
from .tools.calibration_alignment import load_wavelength_table

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

_ROLE_CHOICES: tuple[tuple[str, MeasurementRole | None], ...] = (
    ("— unassigned —", None),
    ("Sphere signal", MeasurementRole.SPHERE),
    ("Sphere background", MeasurementRole.SPHERE_BACKGROUND),
    ("Lamp signal", MeasurementRole.LAMP),
    ("Lamp background", MeasurementRole.LAMP_BACKGROUND),
    ("Experiment / other", MeasurementRole.OTHER),
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

#: Suffix that makes an unconfirmed suggestion unmistakable in the Role column.
_SUGGESTED_SUFFIX = " · SUGGESTED"

_SUGGESTED_TOOLTIP = (
    "This role is only a guess read from the filename. The bench has NOT been "
    "given it. Pick the role in this list — or press “Confirm suggested "
    "roles” — before the procedure and the factor computation can use the file."
)

# ----------------------------------------------------------------------
# Typography: text is data here, not decoration
# ----------------------------------------------------------------------
#: Floor for ordinary bench prose, in points.  Point sizes rather than pixels
#: so a lower-resolution display scales the whole surface with the platform
#: font instead of shrinking it.
BENCH_BODY_POINT_SIZE = 11.0

#: Floor for the readings the operator judges an exposure by.
BENCH_READING_POINT_SIZE = 14.0

#: Floor for the one-glance verdict headline.
BENCH_HEADLINE_POINT_SIZE = 17.0


def bench_point_sizes(base_point_size: float | None = None) -> tuple[float, float, float]:
    """Return (body, reading, headline) point sizes for this platform's font.

    The floors are absolute: a small platform font never shrinks the bench
    below legibility, and a large one scales every tier with it together.
    """

    if base_point_size is None or base_point_size <= 0:
        application = QtWidgets.QApplication.instance()
        base_point_size = (
            application.font().pointSizeF() if application is not None else 0.0
        )
    if base_point_size <= 0:  # a pixel-sized platform font reports no points
        base_point_size = 9.0
    body = max(BENCH_BODY_POINT_SIZE, float(base_point_size) * 1.15)
    return (
        body,
        max(BENCH_READING_POINT_SIZE, body * 1.28),
        max(BENCH_HEADLINE_POINT_SIZE, body * 1.55),
    )


def _emphasise(widget, point_size: float, *, bold: bool = True) -> None:
    """Render one widget's text at *point_size*, keeping its family."""

    font = widget.font()
    font.setPointSizeF(float(point_size))
    font.setBold(bold)
    widget.setFont(font)


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
    icon.addPixmap(
        tint.scaled(32, 32, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
    )
    return icon


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
        self._campaign_thread: CampaignTaskThread | None = None
        self._pattern_items: list[pg.PlotDataItem] = []
        self._line_items: list[object] = []
        self._catalog_items: list[object] = []
        self._anchor_items: list[object] = []
        self._queue: list[Path] = []
        self._file_rows: list[Path] = []
        self._explainable_widgets: list[object] = []
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
        self.resize(1540, 940)
        self.setMinimumSize(1080, 700)
        self.body_pt, self.reading_pt, self.headline_pt = bench_point_sizes()

        root = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        root.setChildrenCollapsible(False)
        root.setHandleWidth(8)
        self.root_splitter = root
        self.setCentralWidget(root)

        controls = QtWidgets.QWidget()
        # No maximum: the left pane is the operator's to widen, and its text is
        # the reading, not a caption.
        controls.setMinimumWidth(420)
        controls_layout = QtWidgets.QVBoxLayout(controls)
        controls_layout.setContentsMargins(18, 18, 14, 18)
        controls_layout.setSpacing(10)

        title = QtWidgets.QLabel("LIVE CALIBRATION")
        title.setObjectName("benchTitle")
        subtitle = QtWidgets.QLabel(
            "Drop files → triage exposure → assign roles → identify → compare → validate"
        )
        subtitle.setWordWrap(True)
        subtitle.setObjectName("benchSubtitle")
        controls_layout.addWidget(title)
        controls_layout.addWidget(subtitle)

        self.control_tabs = QtWidgets.QTabWidget()
        self._build_files_tab()
        self._build_procedure_tab()
        self._build_lamp_tab()
        self._build_save_tab()
        controls_layout.addWidget(self.control_tabs, 1)
        root.addWidget(controls)

        self.view_tabs = QtWidgets.QTabWidget()
        self._build_triage_view()
        self._build_alignment_view()
        self._build_line_help_view()
        self._build_sphere_view()
        root.addWidget(self.view_tabs)
        root.setStretchFactor(0, 0)
        root.setStretchFactor(1, 1)
        root.setSizes([560, 980])
        root.splitterMoved.connect(lambda *_args: self._relayout_wrapped_text())

        self._build_details_dock()

        body = self.body_pt
        reading = self.reading_pt
        headline = self.headline_pt
        # Sizes are set through real fonts rather than through the stylesheet:
        # a stylesheet font-size is resolved through the platform's pixel
        # metrics and quietly comes back smaller on some displays, which is
        # exactly the squinting this is meant to end.
        base_font = self.font()
        base_font.setPointSizeF(body)
        self.setFont(base_font)
        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #151b22; color: #dce8f2; }
            QGroupBox { border: 1px solid #334252; border-radius: 7px; margin-top: 9px;
                        padding-top: 8px; font-weight: 600; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; color: #8fd9ff; }
            #benchTitle { color: #80ddff; font-weight: 700; letter-spacing: 1px; }
            #benchSubtitle, #benchHelp, #mutedText { color: #93a8b8; }
            #dropTarget { border: 2px dashed #49b5df; border-radius: 9px;
                          color: #8fd9ff; font-weight: 700; letter-spacing: 1px;
                          padding: 12px; }
            #triageHeadline { font-weight: 700; padding: 8px; }
            #stateBadge { color: #7ee2b8; font-weight: 700; }
            #reading { color: #e8f4ff; font-weight: 700; }
            #messagePanel { background: #0f141a; border-left: 3px solid #49b5df;
                            padding: 9px; color: #bed4e1; }
            QTableWidget, QTreeWidget, QListWidget, QPlainTextEdit, QTextBrowser {
                background: #10151b; alternate-background-color: #18212a;
                gridline-color: #2b3946; }
            QHeaderView::section { background: #202b36; color: #b9d5e5; padding: 5px; }
            QPushButton { background: #273746; border: 1px solid #416078; border-radius: 5px;
                          padding: 7px; }
            QPushButton:hover { background: #315069; }
            QPushButton:disabled { color: #657786; border-color: #33404a; }
            QSpinBox, QComboBox, QLineEdit {
                background: #0f141a; border: 1px solid #416078; padding: 4px; }
            QSplitter::handle { background: #22303c; }
            QSplitter::handle:hover { background: #3a5f79; }
            QTabWidget::pane { border: 1px solid #334252; }
            QTabBar::tab { background: #202b36; color: #9fb6c6; padding: 7px 12px; }
            QTabBar::tab:selected { background: #294052; color: #8fe3ff; }
            QToolTip { background: #0f141a; color: #dce8f2;
                       border: 1px solid #49b5df; padding: 6px; }
            """
        )
        _emphasise(title, headline * 1.15)
        _emphasise(self.drop_hint, reading)
        for widget in (
            self.rms_value,
            self.anchor_count_value,
            self.comparison_value,
            self.file_state_value,
            self.alignment_state_value,
            self.save_state_value,
        ):
            widget.setObjectName("reading")
            _emphasise(widget, reading)
        _emphasise(self.triage_headline, headline)
        _emphasise(self.exposure_value, reading, bold=False)

    def _build_details_dock(self) -> None:
        """One dock that explains whatever the operator just clicked.

        The whys are learned once and then remembered, so they live behind a
        hover or a click instead of standing permanently between the operator
        and the numbers.
        """

        dock = QtWidgets.QDockWidget("Why this reading", self)
        dock.setObjectName("benchDetailsDock")
        dock.setAllowedAreas(QtCore.Qt.BottomDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        self.details_view = QtWidgets.QTextBrowser()
        self.details_view.setOpenExternalLinks(False)
        self.details_view.setMinimumHeight(120)
        self.details_view.setPlaceholderText(
            "Click any verdict, checklist row, file, or anchor to read the full "
            "explanation here. Hovering shows the same text as a tooltip."
        )
        dock.setWidget(self.details_view)
        self.details_dock = dock
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
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

    def _explainable(self, widget, title: str, text: str) -> None:
        """Give one widget its tooltip and make a click fill the details dock."""

        widget.setToolTip(text)
        if isinstance(widget, QtWidgets.QLabel):
            # A button already looks clickable; a reading does not, so only the
            # readings advertise that there is more behind them.
            widget.setCursor(QtCore.Qt.WhatsThisCursor)
        widget.setProperty("explainTitle", title)
        widget.setProperty("explainText", text)
        if widget not in self._explainable_widgets:
            self._explainable_widgets.append(widget)
            widget.installEventFilter(self)

    def eventFilter(self, source, event) -> bool:  # noqa: N802 - Qt naming
        """Turn a click on any explained reading into its dock entry."""

        if event.type() == QtCore.QEvent.MouseButtonPress:
            title = source.property("explainTitle")
            if title:
                self.explain(title, source.property("explainText") or "")
        return super().eventFilter(source, event)

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt naming
        super().resizeEvent(event)
        self._relayout_wrapped_text()

    def _relayout_wrapped_text(self) -> None:
        """Re-wrap the checklist to the pane's real width; nothing is clipped."""

        if getattr(self, "checklist_tree", None) is None:
            return
        width = max(280, self.checklist_tree.viewport().width() - 18)
        for row in range(self.checklist_tree.count()):
            item = self.checklist_tree.item(row)
            label = self.checklist_tree.itemWidget(item)
            if label is None:
                continue
            label.setFixedWidth(width)
            item.setSizeHint(QtCore.QSize(0, label.sizeHint().height() + 8))
        self.file_table.resizeRowsToContents()

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
        layout.addWidget(self.checklist_tree, 1)
        self.control_tabs.addTab(tab, "Procedure")

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
        self.file_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.Stretch
        )
        self.file_table.setColumnWidth(1, 210)
        self.file_table.setColumnWidth(2, 96)
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
        layout.addWidget(self.file_table, 1)

        self.show_frame_button = QtWidgets.QPushButton("Open selected file for lamp fitting")
        layout.addWidget(self.show_frame_button)

        status_group = QtWidgets.QGroupBox("Bench state")
        status_form = QtWidgets.QFormLayout(status_group)
        self.watch_value = QtWidgets.QLabel("manual — drag and drop or Add files")
        self.watch_value.setWordWrap(True)
        self.file_value = QtWidgets.QLabel("no file open")
        self.file_value.setWordWrap(True)
        self.file_state_value = QtWidgets.QLabel("WAITING")
        self.file_state_value.setObjectName("stateBadge")
        status_form.addRow("Input", self.watch_value)
        status_form.addRow("Open frame", self.file_value)
        status_form.addRow("File state", self.file_state_value)
        layout.addWidget(status_group)

        self.exposure_value = QtWidgets.QLabel(
            "Exposure guidance for the selected file appears here."
        )
        self.exposure_value.setWordWrap(True)
        self.exposure_value.setObjectName("messagePanel")
        layout.addWidget(self.exposure_value)

        comparison_group = QtWidgets.QGroupBox("Integrating-sphere factors")
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
        self.compare_button = QtWidgets.QPushButton("Compute and compare factors")
        comparison_layout.addWidget(self.comparison_value)
        comparison_layout.addWidget(self.compare_button)
        layout.addWidget(comparison_group)
        self.control_tabs.addTab(tab, "Files")

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
        family_row.addWidget(QtWidgets.QLabel("Line help"))
        self.line_family_combo = QtWidgets.QComboBox()
        self.line_family_combo.addItems(["ThAr", "Ne", "Hg", "H2"])
        family_row.addWidget(self.line_family_combo)
        order_layout.addLayout(family_row)
        help_text = QtWidgets.QLabel(
            "Blue sticks use shared packaged line knowledge. Gold rows remain the "
            "curated click-to-fit anchors. Saturated raw pixels are refused."
        )
        help_text.setWordWrap(True)
        help_text.setObjectName("benchHelp")
        order_layout.addWidget(help_text)
        layout.addWidget(order_group)

        self.reference_value = QtWidgets.QLabel(
            "No lamp catalog is scoping the fit yet."
        )
        self.reference_value.setWordWrap(True)
        self.reference_value.setObjectName("messagePanel")
        layout.addWidget(self.reference_value)

        fit_group = QtWidgets.QGroupBox("Rigid alignment")
        fit_form = QtWidgets.QFormLayout(fit_group)
        self.alignment_state_value = QtWidgets.QLabel("WAITING FOR FRAME")
        self.alignment_state_value.setObjectName("stateBadge")
        self.anchor_count_value = QtWidgets.QLabel("0")
        self.rms_value = QtWidgets.QLabel("—")
        self.transform_value = QtWidgets.QLabel("—")
        self.transform_value.setWordWrap(True)
        fit_form.addRow("State", self.alignment_state_value)
        fit_form.addRow("Anchors", self.anchor_count_value)
        fit_form.addRow("RMS", self.rms_value)
        fit_form.addRow("dx / dy / θ", self.transform_value)
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
        layout.addWidget(fit_group)

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

        button_row = QtWidgets.QHBoxLayout()
        self.remove_button = QtWidgets.QPushButton("Remove selected")
        self.clear_button = QtWidgets.QPushButton("Clear anchors")
        button_row.addWidget(self.remove_button)
        button_row.addWidget(self.clear_button)
        layout.addLayout(button_row)

        self.message_value = QtWidgets.QLabel("Waiting for data.")
        self.message_value.setWordWrap(True)
        self.message_value.setObjectName("messagePanel")
        layout.addWidget(self.message_value)
        self.control_tabs.addTab(tab, "Lamp fit")

    def _build_save_tab(self) -> None:
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        layout.setContentsMargins(10, 12, 10, 10)

        identity_group = QtWidgets.QGroupBox("Snapshot identity")
        form = QtWidgets.QFormLayout(identity_group)
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
        self.generate_tomls_button = QtWidgets.QPushButton("Generate commented TOMLs")
        self.save_snapshot_button = QtWidgets.QPushButton("Save and validate snapshot")
        layout.addWidget(self.generate_tomls_button)
        layout.addWidget(self.save_snapshot_button)
        self.save_state_value = QtWidgets.QLabel("NOT READY")
        self.save_state_value.setObjectName("stateBadge")
        layout.addWidget(self.save_state_value)
        self.toml_preview = QtWidgets.QPlainTextEdit()
        self.toml_preview.setReadOnly(True)
        self.toml_preview.setPlaceholderText(
            "Generated campaign.toml appears here; all files remain ordinary and editable."
        )
        layout.addWidget(self.toml_preview, 1)
        self.control_tabs.addTab(tab, "Save")

    def _build_triage_view(self) -> None:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        self.triage_headline = QtWidgets.QLabel(
            "Drop a SIF onto the bench. Triage needs nothing but a file."
        )
        self.triage_headline.setObjectName("triageHeadline")
        self.triage_headline.setWordWrap(True)
        layout.addWidget(self.triage_headline)
        self.triage_detail = QtWidgets.QLabel(
            "Shoot → drop the file → one glance → adjust the lamp → shoot again."
        )
        self.triage_detail.setWordWrap(True)
        self.triage_detail.setObjectName("messagePanel")
        layout.addWidget(self.triage_detail)

        graphics = pg.GraphicsLayoutWidget()
        graphics.setBackground("#10151b")
        self.histogram_plot = graphics.addPlot(row=0, col=0, title="Raw counts histogram")
        self.histogram_plot.setLabel("bottom", "counts")
        self.histogram_plot.setLabel("left", "pixels per bin")
        self.histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.histogram_plot.setLogMode(y=True)
        self.top_histogram_plot = graphics.addPlot(
            row=1, col=0, title="Top end — the last 10% before full scale"
        )
        self.top_histogram_plot.setLabel("bottom", "counts")
        self.top_histogram_plot.setLabel("left", "pixels per bin")
        self.top_histogram_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.top_histogram_plot.setLogMode(y=True)
        layout.addWidget(graphics, 1)
        self.view_tabs.addTab(widget, "Exposure triage")

    def _build_alignment_view(self) -> None:
        graphics = pg.GraphicsLayoutWidget()
        graphics.setBackground("#10151b")

        self.detector_plot = graphics.addPlot(row=0, col=0, title="Detector + order traces")
        self.detector_plot.setLabel("bottom", "detector column", units="px")
        self.detector_plot.setLabel("left", "detector row", units="px")
        self.detector_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.detector_plot.getAxis("left").enableAutoSIPrefix(False)
        self.detector_plot.invertY(True)
        self.detector_plot.setAspectLocked(False)
        self.detector_image = pg.ImageItem(axisOrder="col-major")
        self.detector_plot.addItem(self.detector_image)

        self.order_plot = graphics.addPlot(row=1, col=0, title="Selected order spectrum")
        self.order_plot.setLabel("bottom", "raw detector column", units="px")
        self.order_plot.setLabel("left", "mean extracted counts")
        self.order_plot.getAxis("bottom").enableAutoSIPrefix(False)
        self.order_plot.getAxis("left").enableAutoSIPrefix(False)
        self.order_curve = self.order_plot.plot(pen=pg.mkPen("#76d6ff", width=1.4))

        self.residual_plot = graphics.addPlot(row=2, col=0, title="Anchor residuals")
        self.residual_plot.setMaximumHeight(190)
        self.residual_plot.setLabel("bottom", "accepted anchor")
        self.residual_plot.setLabel("left", "fit residual", units="px")
        self.residual_plot.getAxis("bottom").setHeight(62)
        self.residual_plot.getAxis("left").enableAutoSIPrefix(False)
        self.residual_plot.addLine(y=0, pen=pg.mkPen("#64748b", style=QtCore.Qt.DashLine))
        self.view_tabs.addTab(graphics, "Lamp alignment")

    def _build_line_help_view(self) -> None:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        intro = QtWidgets.QLabel(
            "Expected lines come from the same packaged ThAr, Ne, Hg, and Fulcher H2 "
            "catalogs used by the main GUI and validation tools. Pixel positions are "
            "interpolated from the current wavelength table for identification help."
        )
        intro.setWordWrap(True)
        intro.setObjectName("messagePanel")
        layout.addWidget(intro)
        self.line_help_table = QtWidgets.QTableWidget(0, 5)
        self.line_help_table.setHorizontalHeaderLabels(
            ["Label", "λ (nm)", "Pixel", "Rel. I", "Packaged source"]
        )
        self.line_help_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeToContents
        )
        self.line_help_table.horizontalHeader().setSectionResizeMode(
            4, QtWidgets.QHeaderView.Stretch
        )
        self.line_help_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.line_help_table.setAlternatingRowColors(True)
        layout.addWidget(self.line_help_table)
        self.view_tabs.addTab(widget, "Line identification")

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
        layout.addWidget(self.sphere_plot)
        self.view_tabs.addTab(widget, "Sphere factors")

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
        self.remove_button.clicked.connect(self._remove_selected_anchor)
        self.clear_button.clicked.connect(self._clear_anchors)
        self.add_files_button.clicked.connect(self._pick_files)
        self.remove_file_button.clicked.connect(self._remove_selected_file)
        self.confirm_roles_button.clicked.connect(self._confirm_suggested_roles)
        self.show_frame_button.clicked.connect(self._open_selected_file)
        self.file_table.itemSelectionChanged.connect(self._file_selection_changed)
        self.checklist_tree.currentRowChanged.connect(self._checklist_row_selected)
        self.anchor_table.itemSelectionChanged.connect(self._anchor_row_selected)
        self.compare_button.clicked.connect(self._start_sphere_comparison)
        self.generate_tomls_button.clicked.connect(self._generate_tomls)
        self.save_snapshot_button.clicked.connect(self._save_snapshot)
        self.snapshot_id_edit.textChanged.connect(self.refresh_campaign)

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
        if confirmed:
            self.message_value.setText(
                f"Assigned {len(confirmed)} suggested role(s): {names}. "
                "Change any of them in the Role column."
            )
        else:
            self.message_value.setText(
                "No suggestion could be confirmed; assign those roles by hand."
            )
        self.refresh()

    def _line_family_changed(self) -> None:
        self._refresh_reference()
        self.refresh_plots()

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
        if selected is not None:
            self._show_triage(selected)

    def _open_selected_file(self) -> None:
        selected = self._selected_file()
        if selected is None:
            return
        if self.session.frame is not None and self.session.frame.path == selected:
            self.message_value.setText(f"{selected.name} is already open for fitting.")
            return
        self.add_paths([selected])

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

    def _generate_tomls(self) -> None:
        if self.campaign is None:
            return
        try:
            paths = self.campaign.write_tomls(
                self.config_root,
                self.snapshot_id_edit.text().strip(),
                self.session,
            )
        except (OSError, SnapshotError, ValueError) as exc:
            self.message_value.setText(f"TOMLs were not generated: {exc}")
        else:
            self.message_value.setText(
                "Generated commented campaign, alignment, and export TOMLs."
            )
            self.toml_preview.setPlainText(
                paths["campaign"].read_text(encoding="utf-8")
            )
        self.refresh_campaign()

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

    def _order_plot_clicked(self, event) -> None:
        if event.button() != QtCore.Qt.LeftButton or self.session.frame is None:
            return
        if not self.order_plot.sceneBoundingRect().contains(event.scenePos()):
            return
        point = self.order_plot.getViewBox().mapSceneToView(event.scenePos())
        result = self.session.fit_anchor_at(self.session.selected_order, point.x())
        self.message_value.setText(result.reason)
        self.refresh()

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
        self.triage_headline.setText(f"{path.name}\n{headline}")
        self.triage_headline.setStyleSheet(f"color: {color};")
        self._explainable(
            self.triage_headline,
            f"{path.name} — exposure verdict",
            f"{headline}\n\n{verdict.advice}",
        )
        self.triage_detail.setText("\n".join(triage.details))
        self._explainable(
            self.triage_detail,
            f"{path.name} — how the verdict was reached",
            "\n".join(triage.details)
            + "\n\nIsolated full-scale pixels are cosmic rays or hot pixels and "
            "are counted, not held against the frame. Only a connected cluster "
            "of full-scale pixels is real saturation.",
        )
        self._draw_histogram(self.histogram_plot, triage, triage.histogram)
        self._draw_histogram(self.top_histogram_plot, triage, triage.top_histogram)
        guidance = record.exposure
        peak = "—" if guidance.peak_value is None else f"{guidance.peak_value:.0f} counts"
        self.exposure_value.setText(
            f"{path.name} · {verdict.label.upper()} · peak {peak}. "
            f"{guidance.next_action}"
        )
        self._explainable(
            self.exposure_value,
            f"{path.name} — next acquisition action",
            f"{guidance.next_action}\n\n{verdict.advice}",
        )

    @staticmethod
    def _draw_histogram(plot, triage: ExposureTriage, histogram) -> None:
        plot.clear()
        counts = np.maximum(np.asarray(histogram.counts, dtype=float), 0.5)
        plot.plot(
            np.asarray(histogram.edges, dtype=float),
            counts,
            stepMode="center",
            fillLevel=0,
            brush=pg.mkBrush("#1f4d63"),
            pen=pg.mkPen("#76d6ff", width=1.2),
        )
        plot.addLine(
            x=triage.saturation.saturation_level,
            pen=pg.mkPen("#ffb86b", style=QtCore.Qt.DashLine),
        )
        plot.addLine(
            x=triage.full_scale,
            pen=pg.mkPen("#ff8f8f", style=QtCore.Qt.DashLine),
        )

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
        self._refresh_reference()
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
        selected = self._selected_file()
        self.remove_file_button.setEnabled(selected is not None)
        self.show_frame_button.setEnabled(
            selected is not None and self.loader is not None and self._load_thread is None
        )
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
                marks.append(
                    "SUGGESTED ONLY — no role assigned"
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
                explanation = f"{path}\n{verdict.headline}\n\n{verdict.advice}"
                if measurement is None and path in suggested:
                    explanation = f"{path}\n{_SUGGESTED_TOOLTIP}\n\n{verdict.headline}"
                item.setToolTip(explanation)
            lamp_combo = self.file_table.cellWidget(row, 2)
            lamp_combo.setEnabled(
                measurement is not None and measurement.role in _LAMP_ROLES
            )

    @staticmethod
    def _render_role_combo(role_combo, is_only_suggested: bool) -> None:
        """Make an unconfirmed suggestion impossible to read as an assignment.

        The entry's *data* is untouched, so the control still proposes the
        filename's guess; only its text and colour say that the bench has not
        been given that role yet.
        """

        if role_combo is None:
            return
        for index in range(role_combo.count()):
            label = role_combo.itemText(index).removesuffix(_SUGGESTED_SUFFIX)
            wanted = (
                f"{label}{_SUGGESTED_SUFFIX}"
                if is_only_suggested
                and index == role_combo.currentIndex()
                and role_combo.itemData(index) is not None
                else label
            )
            if wanted != role_combo.itemText(index):
                role_combo.setItemText(index, wanted)
        role_combo.setToolTip(
            _SUGGESTED_TOOLTIP
            if is_only_suggested
            else "The measurement role this file carries. Change it freely; "
            "the bench never infers one behind your back."
        )
        role_combo.setStyleSheet(
            f"border: 1px solid {_SUGGESTED_COLOR}; color: {_SUGGESTED_COLOR};"
            if is_only_suggested
            else ""
        )

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
        width = max(280, self.checklist_tree.viewport().width() - 18)
        for item in self.campaign.checklist(self.session):
            text = f"{symbols[item.state]}  {item.label} — {item.detail}"
            if item.unblocked_by:
                text += f"\n     unblocked by: {item.unblocked_by}"
            row = QtWidgets.QListWidgetItem()
            tooltip = f"{item.label}\n{item.detail}"
            if item.unblocked_by:
                tooltip += f"\n\nWhat unblocks it: {item.unblocked_by}"
            row.setToolTip(tooltip)
            label = QtWidgets.QLabel(text)
            label.setWordWrap(True)
            label.setContentsMargins(8, 6, 8, 6)
            label.setToolTip(tooltip)
            label.setStyleSheet(f"color: {colors[item.state].name()};")
            _emphasise(label, self.body_pt, bold=item.state is ChecklistState.ATTENTION)
            label.setFixedWidth(width)
            row.setSizeHint(QtCore.QSize(0, label.sizeHint().height() + 8))
            self.checklist_tree.addItem(row)
            self.checklist_tree.setItemWidget(row, label)

    def _refresh_sphere_plot(self) -> None:
        self.sphere_plot.clear()
        if self.campaign is None:
            return
        comparison = self.campaign.comparison
        if comparison.candidate is not None:
            size = min(
                comparison.candidate.wavelength_nm.size,
                comparison.candidate.factors_wmsr.size,
            )
            self.sphere_plot.plot(
                comparison.candidate.wavelength_nm[:size],
                comparison.candidate.factors_wmsr[:size],
                pen=pg.mkPen("#70d6ae", width=2),
                name="new measured pair",
            )
        if comparison.previous is not None:
            size = min(
                comparison.previous.wavelength_nm.size,
                comparison.previous.factors_wmsr.size,
            )
            self.sphere_plot.plot(
                comparison.previous.wavelength_nm[:size],
                comparison.previous.factors_wmsr[:size],
                pen=pg.mkPen("#f5b95f", width=1.5),
                name="previous campaign",
            )
        self.sphere_plot.setLogMode(y=True)
        self.sphere_view_message.setText(
            f"{comparison.state.value.replace('-', ' ').upper()} — {comparison.reason}"
        )

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

    def refresh_plots(self) -> None:
        frame = self.session.frame
        if frame is None:
            return
        finite = frame.detector_image[np.isfinite(frame.detector_image)]
        levels = None
        if finite.size:
            low, high = np.percentile(finite, (1.0, 99.8))
            if high > low:
                levels = (float(low), float(high))
        self.detector_image.setImage(
            frame.detector_image.T,
            autoLevels=levels is None,
            levels=levels,
        )
        self._refresh_pattern_traces()

        order_idx = self.session.selected_order
        spectrum = self.session.active_order_spectra()[order_idx]
        self.order_curve.setData(np.arange(spectrum.size), spectrum)
        self.order_plot.setTitle(
            f"Order {order_idx} · {self.session.frame_choice_label()}: "
            "click a labeled expected line"
        )
        self._clear_items(self.order_plot, self._line_items)
        self._clear_items(self.order_plot, self._catalog_items)
        self._clear_items(self.order_plot, self._anchor_items)
        top = float(np.nanmax(spectrum)) if np.any(np.isfinite(spectrum)) else 1.0
        line_rows = self.session.lines_for_order(order_idx)
        for index, line in enumerate(line_rows):
            marker = pg.InfiniteLine(
                line.center_pixel,
                angle=90,
                movable=False,
                pen=pg.mkPen("#f5b95f", width=1, style=QtCore.Qt.DashLine),
            )
            label = pg.TextItem(
                f"{line.species} {line.wavelength_nm:.3f}",
                color="#f5c56f",
                anchor=(0.5, 0.5),
            )
            label_level = (0.90, 0.72, 0.54, 0.36)[index % 4]
            label.setPos(line.center_pixel, top * label_level)
            self.order_plot.addItem(marker, ignoreBounds=True)
            self.order_plot.addItem(label, ignoreBounds=True)
            self._line_items.extend([marker, label])
        catalog_rows = catalog_lines_for_order(
            self.session.lines,
            order_idx,
            self.line_family_combo.currentText(),
            maximum_lines=16,
        )
        for index, row in enumerate(catalog_rows):
            marker = pg.InfiniteLine(
                row.detector_pixel,
                angle=90,
                movable=False,
                pen=pg.mkPen("#6aa7ff", width=0.9, style=QtCore.Qt.DotLine),
            )
            self.order_plot.addItem(marker, ignoreBounds=True)
            self._catalog_items.append(marker)
            if index % 3 == 0:
                label = pg.TextItem(
                    row.line.label,
                    color="#8ebcff",
                    anchor=(0.5, 0.5),
                )
                label.setPos(row.detector_pixel, top * (0.22 + 0.09 * (index % 3)))
                self.order_plot.addItem(label, ignoreBounds=True)
                self._catalog_items.append(label)
        for anchor in self.session.anchor_rows():
            if anchor.line.order_idx != order_idx:
                continue
            marker = pg.ScatterPlotItem(
                [anchor.fit.center_pixel],
                [anchor.fit.amplitude + anchor.fit.baseline],
                size=11,
                brush=pg.mkBrush("#6ee7b7"),
                pen=pg.mkPen("#d7fff1", width=1),
            )
            self.order_plot.addItem(marker)
            self._anchor_items.append(marker)
        self._refresh_line_help_table(catalog_rows)
        self._refresh_residual_plot()

    def _refresh_line_help_table(self, rows) -> None:
        self.line_help_table.setRowCount(len(rows))
        for index, row in enumerate(rows):
            intensity = row.line.relative_intensity
            source = row.line.source_resource.replace("\\", "/").rsplit("/", 1)[-1]
            values = (
                row.line.label,
                f"{row.line.wavelength_nm:.4f}",
                f"{row.detector_pixel:.1f}",
                "—" if intensity is None else f"{intensity:.2f}",
                source,
            )
            for column, value in enumerate(values):
                self.line_help_table.setItem(
                    index, column, QtWidgets.QTableWidgetItem(value)
                )

    def _refresh_pattern_traces(self) -> None:
        while self._pattern_items:
            self.detector_plot.removeItem(self._pattern_items.pop())
        columns = np.arange(self.session.pattern.shape[0])
        for order_idx in range(self.session.pattern.shape[1]):
            selected = order_idx == self.session.selected_order
            item = self.detector_plot.plot(
                columns,
                self.session.pattern[:, order_idx],
                pen=pg.mkPen(
                    "#6ee7b7" if selected else "#7b91a4",
                    width=2.2 if selected else 0.7,
                ),
            )
            item.setZValue(12 if selected else 10)
            self._pattern_items.append(item)

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

    @staticmethod
    def _clear_items(plot, items: list[object]) -> None:
        while items:
            plot.removeItem(items.pop())


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
    pg.setConfigOptions(antialias=True, imageAxisOrder="col-major")
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
