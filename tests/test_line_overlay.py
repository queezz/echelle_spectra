"""Headless GUI tests for known-line overlay toggles and zoom behavior."""

from __future__ import annotations

import colorsys
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pyqtgraph as pg
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.resources.window_layout import Ui_MainWindow
from echelle_spectra.tools.image_line_overlay import _brighter
from echelle_spectra.tools.line_catalog import LINE_FAMILIES
from echelle_spectra.tools.line_overlay import (
    LINE_OVERLAY_STYLES,
    SPECTRUM_CURVE_COLORS,
    LineOverlayManager,
    select_overlay_lines,
)


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
    assert min(line.relative_intensity for line in broad) >= 0.08
    assert min(line.relative_intensity for line in narrow) < 0.08


def test_a_broad_view_of_a_neon_lamp_is_not_a_view_of_neon_ions():
    """The whole detector, and the marks land where the lamp puts light.

    Pinned on the 401.5--801.9 nm span of the 20250926 CMOS solution, which is
    what the 2-D detector overlay asks for.  Before the cache was widened this
    returned 33 rows spanning 580--638 nm, 16 of them Ne II.
    """

    selected = select_overlay_lines("ne", 401.5, 801.9, max_labels=200)
    assert selected
    ions = [line for line in selected if line.species == "Ne II"]
    assert len(ions) <= 2, [line.label for line in ions]

    # The strong red neon lines the old cache stopped short of, including the
    # brightest blob on the owner's frame at 640.2248 nm.
    for wavelength in (640.2248, 650.6528, 692.9467, 703.2413, 724.5167):
        assert any(
            abs(line.wavelength_nm - wavelength) < 0.01 for line in selected
        ), f"{wavelength} nm missing from the broad Ne selection"


def test_no_family_colour_hides_in_a_spectrum_curve():
    """Overlay hues must contrast with both curves they are drawn over.

    The Ne family shipped a green stick over the green counts curve, which is
    a marker that annotates by disappearing.  Red and yellow-green bracket the
    warm arc, so the palette lives in the cool one; this pins that it stays
    there and that no two families collapse into each other either.
    """

    def hue(color: str) -> float:
        red, green, blue, _alpha = pg.mkColor(color).getRgb()
        return colorsys.rgb_to_hls(red / 255, green / 255, blue / 255)[0] * 360.0

    def separation(first: float, second: float) -> float:
        gap = abs(first - second) % 360.0
        return min(gap, 360.0 - gap)

    curves = {name: hue(color) for name, color in SPECTRUM_CURVE_COLORS.items()}
    families = {name: hue(style.color) for name, style in LINE_OVERLAY_STYLES.items()}
    assert set(families) == set(LINE_FAMILIES)

    for family, family_hue in families.items():
        for curve, curve_hue in curves.items():
            assert separation(family_hue, curve_hue) >= 40.0, (
                f"{family} at {family_hue:.0f} deg is too close to the "
                f"{curve} curve at {curve_hue:.0f} deg"
            )

    names = sorted(families)
    for index, first in enumerate(names):
        for second in names[index + 1 :]:
            assert separation(families[first], families[second]) >= 25.0, (
                f"{first} and {second} are not far enough apart in hue"
            )


def test_the_detector_boxes_use_the_same_family_hue_as_the_sticks():
    """One family, one hue, both views — the boxes only lift it toward white."""

    for family, style in LINE_OVERLAY_STYLES.items():
        lifted = _brighter(style.color)
        base = pg.mkColor(style.color).getRgb()[:3]
        lit = pg.mkColor(lifted).getRgb()[:3]
        assert all(after >= before for before, after in zip(base, lit)), family
        assert lit != base, family
        assert (
            colorsys.rgb_to_hls(*(value / 255 for value in lit))[2] > 0.15
        ), f"{family} washes out to grey on the detector image"
