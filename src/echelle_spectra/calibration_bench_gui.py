"""Separate pyqtgraph live calibration bench for ``echelle-calib``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
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
from .tools.calibration_alignment import load_wavelength_table

_PACKAGE_DIR = Path(__file__).parent
_CALIBRATION_DIR = _PACKAGE_DIR / "resources" / "calibration_files"
_DEFAULT_PATTERN = _CALIBRATION_DIR / "pattern_CMOS_20250926.txt"
_DEFAULT_WAVELENGTH = (
    _CALIBRATION_DIR
    / "alignments"
    / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
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


class CalibrationBenchWindow(QtWidgets.QMainWindow):
    """Thin Qt adapter over :class:`CalibrationBenchSession`."""

    def __init__(
        self,
        session: CalibrationBenchSession,
        *,
        watcher: StableSifWatcher | None = None,
        loader: FrameLoader | None = None,
        poll_interval_ms: int = 1000,
        start_timer: bool = True,
    ) -> None:
        super().__init__()
        self.session = session
        self.watcher = watcher
        self.loader = loader
        self._load_thread: FrameLoadThread | None = None
        self._pattern_items: list[pg.PlotDataItem] = []
        self._line_items: list[object] = []
        self._anchor_items: list[object] = []
        self._build_ui()
        self._connect_ui()
        self.refresh()

        self.poll_timer = QtCore.QTimer(self)
        self.poll_timer.setInterval(int(poll_interval_ms))
        self.poll_timer.timeout.connect(self.poll_watch_folder)
        if start_timer and self.watcher is not None:
            self.poll_timer.start()
            QtCore.QTimer.singleShot(0, self.poll_watch_folder)

    def _build_ui(self) -> None:
        self.setWindowTitle("Echelle calibration bench")
        self.resize(1440, 900)
        self.setMinimumSize(980, 640)

        root = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        root.setChildrenCollapsible(False)
        self.setCentralWidget(root)

        controls = QtWidgets.QWidget()
        controls.setMinimumWidth(300)
        controls.setMaximumWidth(390)
        controls_layout = QtWidgets.QVBoxLayout(controls)
        controls_layout.setContentsMargins(18, 18, 14, 18)
        controls_layout.setSpacing(12)

        title = QtWidgets.QLabel("LIVE CALIBRATION")
        title.setObjectName("benchTitle")
        subtitle = QtWidgets.QLabel(
            "Newest stable SIF → detector/order review → click known lines → rigid fit"
        )
        subtitle.setWordWrap(True)
        subtitle.setObjectName("benchSubtitle")
        controls_layout.addWidget(title)
        controls_layout.addWidget(subtitle)

        status_group = QtWidgets.QGroupBox("Acquisition")
        status_form = QtWidgets.QFormLayout(status_group)
        self.watch_value = QtWidgets.QLabel("manual")
        self.watch_value.setWordWrap(True)
        self.file_value = QtWidgets.QLabel("Waiting for a stable SIF")
        self.file_value.setWordWrap(True)
        self.file_state_value = QtWidgets.QLabel("WAITING")
        self.file_state_value.setObjectName("stateBadge")
        status_form.addRow("Folder", self.watch_value)
        status_form.addRow("Newest", self.file_value)
        status_form.addRow("File state", self.file_state_value)
        controls_layout.addWidget(status_group)

        order_group = QtWidgets.QGroupBox("Order interaction")
        order_layout = QtWidgets.QVBoxLayout(order_group)
        order_row = QtWidgets.QHBoxLayout()
        order_row.addWidget(QtWidgets.QLabel("Selected order"))
        self.order_spin = QtWidgets.QSpinBox()
        self.order_spin.setRange(0, self.session.pattern.shape[1] - 1)
        order_row.addWidget(self.order_spin)
        order_layout.addLayout(order_row)
        help_text = QtWidgets.QLabel(
            "Click a labeled row. Saturated raw pixels are refused before fitting."
        )
        help_text.setWordWrap(True)
        help_text.setObjectName("benchHelp")
        order_layout.addWidget(help_text)
        controls_layout.addWidget(order_group)

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
        controls_layout.addWidget(fit_group)

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
        controls_layout.addWidget(self.anchor_table, 1)

        button_row = QtWidgets.QHBoxLayout()
        self.remove_button = QtWidgets.QPushButton("Remove selected")
        self.clear_button = QtWidgets.QPushButton("Clear anchors")
        button_row.addWidget(self.remove_button)
        button_row.addWidget(self.clear_button)
        controls_layout.addLayout(button_row)

        self.message_value = QtWidgets.QLabel("Waiting for data.")
        self.message_value.setWordWrap(True)
        self.message_value.setObjectName("messagePanel")
        controls_layout.addWidget(self.message_value)
        root.addWidget(controls)

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
        root.addWidget(graphics)
        root.setStretchFactor(1, 1)

        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #151b22; color: #dce8f2; }
            QGroupBox { border: 1px solid #334252; border-radius: 7px; margin-top: 9px;
                        padding-top: 8px; font-weight: 600; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; color: #8fd9ff; }
            #benchTitle { color: #80ddff; font-size: 21px; font-weight: 700; letter-spacing: 1px; }
            #benchSubtitle, #benchHelp { color: #93a8b8; }
            #stateBadge { color: #7ee2b8; font-weight: 700; }
            #messagePanel { background: #0f141a; border-left: 3px solid #49b5df;
                            padding: 9px; color: #bed4e1; }
            QTableWidget { background: #10151b; alternate-background-color: #18212a;
                           gridline-color: #2b3946; }
            QHeaderView::section { background: #202b36; color: #b9d5e5; padding: 5px; }
            QPushButton { background: #273746; border: 1px solid #416078; border-radius: 5px;
                          padding: 7px; }
            QPushButton:hover { background: #315069; }
            QPushButton:disabled { color: #657786; border-color: #33404a; }
            QSpinBox { background: #0f141a; border: 1px solid #416078; padding: 4px; }
            """
        )

    def _connect_ui(self) -> None:
        self.order_spin.valueChanged.connect(self._order_changed)
        self.order_plot.scene().sigMouseClicked.connect(self._order_plot_clicked)
        self.remove_button.clicked.connect(self._remove_selected_anchor)
        self.clear_button.clicked.connect(self._clear_anchors)

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
            self.load_path(result.ready_path)

    def load_path(self, path: str | Path) -> None:
        """Start an asynchronous load for a stable SIF."""

        if self.loader is None or self._load_thread is not None:
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
        self.refresh()

    @QtCore.pyqtSlot(str, str)
    def _frame_failed(self, path: str, reason: str) -> None:
        self.session.fail_file_load(path, reason)
        self.refresh()

    @QtCore.pyqtSlot()
    def _load_finished(self) -> None:
        if self._load_thread is not None:
            self._load_thread.deleteLater()
        self._load_thread = None

    def _order_changed(self, order_idx: int) -> None:
        self.session.set_selected_order(order_idx)
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

    def refresh(self) -> None:
        """Render the current domain state without changing it."""

        if self.watcher is not None:
            self.watch_value.setText(str(self.watcher.folder))
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
                f"Could not load newest stable SIF: {self.session.last_error}. "
                "The last good frame remains visible; waiting for the next file."
            )
        elif self.session.file_state is FileLoadState.LOADING:
            self.message_value.setText("Loading stable SIF and extracting order spectra…")
        elif self.session.alignment_state is AlignmentState.COLLECTING:
            self.message_value.setText("One anchor accepted. Add another order/line to solve.")
        elif self.session.alignment_state is AlignmentState.ALIGNED:
            self.message_value.setText(
                "Live rigid alignment solved. Review RMS and residuals as anchors accumulate."
            )
        elif self.session.alignment_state is AlignmentState.FAILED:
            self.message_value.setText(f"Alignment failed: {self.session.last_error}")
        self._refresh_anchor_table()
        self.refresh_plots()

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
        spectrum = frame.order_spectra[order_idx]
        self.order_curve.setData(np.arange(spectrum.size), spectrum)
        self.order_plot.setTitle(f"Order {order_idx}: click a labeled expected line")
        self._clear_items(self.order_plot, self._line_items)
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
        self._refresh_residual_plot()

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
            "Watch for stable Andor SIF files and interactively fit a rigid "
            "wavelength-calibration alignment."
        ),
    )
    parser.add_argument(
        "watch_folder",
        nargs="?",
        type=Path,
        default=Path.cwd(),
        help="folder written by the acquisition software (default: current folder)",
    )
    parser.add_argument("--file", type=Path, help="load one SIF immediately and keep watching")
    parser.add_argument("--pattern", type=Path, default=_DEFAULT_PATTERN)
    parser.add_argument("--wavelength", type=Path, default=_DEFAULT_WAVELENGTH)
    parser.add_argument("--poll-ms", type=int, default=1000)
    parser.add_argument("--stable-polls", type=int, default=2)
    parser.add_argument("--minimum-age-s", type=float, default=1.0)
    parser.add_argument("--saturation-level", type=float, default=0.98 * 65535)
    parser.add_argument("--minimum-snr", type=float, default=5.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the calibration bench."""

    args = _build_parser().parse_args(argv)
    if not args.watch_folder.is_dir():
        raise SystemExit(f"watch folder not found: {args.watch_folder}")
    if not args.pattern.is_file():
        raise SystemExit(f"pattern file not found: {args.pattern}")
    if not args.wavelength.is_file():
        raise SystemExit(f"wavelength table not found: {args.wavelength}")
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
    watcher = StableSifWatcher(
        args.watch_folder,
        required_unchanged_polls=args.stable_polls,
        minimum_age_s=args.minimum_age_s,
    )
    loader = FrameLoader(pattern)
    pg.setConfigOptions(antialias=True, imageAxisOrder="col-major")
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    application.setApplicationName("Echelle calibration bench")
    application.setWindowIcon(QtGui.QIcon(str(_PACKAGE_DIR / "resources" / "graphics" / "echelle.png")))
    window = CalibrationBenchWindow(
        session,
        watcher=watcher,
        loader=loader,
        poll_interval_ms=args.poll_ms,
    )
    window.show()
    if args.file is not None:
        QtCore.QTimer.singleShot(0, lambda: window.load_path(args.file))
    return int(application.exec_())


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
