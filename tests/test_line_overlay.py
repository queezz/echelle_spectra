"""Headless GUI tests for known-line overlay toggles and zoom behavior."""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pyqtgraph as pg
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.resources.window_layout import Ui_MainWindow
from echelle_spectra.tools.line_overlay import LineOverlayManager, select_overlay_lines


@pytest.fixture(scope="module")
def qt_app():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def test_gui_exposes_separate_overlay_toggles_off_by_default(qt_app):
    window = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(window)
    assert list(ui.line_overlay_checks) == ["balmer", "fulcher", "thar", "ne", "hg"]
    assert all(not checkbox.isChecked() for checkbox in ui.line_overlay_checks.values())
    window.close()


def test_manager_toggle_zoom_and_disabled_state(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    plot.setXRange(600.0, 616.0, padding=0)
    manager = LineOverlayManager(max_labels=8)
    manager.register_plot("spectrum", plot)

    assert manager.rendered_lines("spectrum", "fulcher") == ()
    manager.set_family_visible("fulcher", True)
    qt_app.processEvents()
    broad = manager.rendered_lines("spectrum", "fulcher")
    assert 1 < len(broad) <= 8
    assert all(600.0 <= line.wavelength_nm <= 616.0 for line in broad)

    plot.setXRange(400.0, 680.0, padding=0)
    qt_app.processEvents()
    manager.refresh()
    assert len(manager.rendered_lines("spectrum", "fulcher")) <= 4

    plot.setXRange(601.7, 603.3, padding=0)
    qt_app.processEvents()
    manager.refresh()
    assert [line.label for line in manager.rendered_lines("spectrum", "fulcher")] == [
        "Q1(0-0)",
        "Q2(0-0)",
        "Q3(0-0)",
    ]

    manager.set_family_visible("fulcher", False)
    assert manager.rendered_lines("spectrum", "fulcher") == ()
    assert not manager.is_family_visible("fulcher")
    widget.close()


def test_atomic_zoom_reveals_weaker_lines_without_overcrowding():
    broad = select_overlay_lines("thar", 578.0, 640.0, max_labels=10)
    narrow = select_overlay_lines("thar", 618.0, 619.0, max_labels=20)
    assert len(broad) <= 10
    assert broad
    assert narrow
    assert min(line.relative_intensity for line in broad) >= 0.25
    assert min(line.relative_intensity for line in narrow) < 0.25
