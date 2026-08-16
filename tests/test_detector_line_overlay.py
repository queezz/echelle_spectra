"""Packet F15 — known-line markers on the 2-D detector image.

The spectrum overlays answer "where is this line in wavelength".  These pin the
geometric half of the same question: with the order pattern and the wavelength
tables loaded, a catalog line has an expected column and an expected row band
on the sensor, and the main GUI's five overlay toggles drive both views.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pyqtgraph as pg
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.tools.image_line_overlay import (
    DetectorGeometry,
    DetectorLineOverlay,
)
from echelle_spectra.tools.line_catalog import SpectralLine

# A three-order sensor with the black Echelle's falling wavelength direction:
# 0.5 nm per column, so order 0 runs 700 nm down to 550 nm and its neighbours
# start 140 nm lower each.  The orders overlap by 10 nm, which is where the
# border mask has to decide who owns a line.
COLUMNS = 301
ORDER_COUNT = 3
DISPERSION_NM = 0.5
ORDER_START_NM = (700.0, 560.0, 420.0)
ORDER_IDS = (21, 22, 23)
HALF_HEIGHT = 8.0


def _order_wavel() -> np.ndarray:
    pixels = np.arange(COLUMNS, dtype=float)
    return np.array([start - DISPERSION_NM * pixels for start in ORDER_START_NM])


def _pattern() -> np.ndarray:
    """Slanted order traces: each order climbs one row every ten columns."""

    climb = np.arange(COLUMNS) // 10
    return np.column_stack(
        [30 + 20 * index + climb for index in range(ORDER_COUNT)]
    ).astype(int)


def _geometry(order_borders=None) -> DetectorGeometry:
    return DetectorGeometry(
        _pattern(),
        _order_wavel(),
        order_ids=ORDER_IDS,
        order_borders=order_borders,
        half_height=HALF_HEIGHT,
    )


def _line(wavelength_nm: float, label: str = "test", family: str = "balmer") -> SpectralLine:
    return SpectralLine(
        family=family,
        label=label,
        wavelength_nm=wavelength_nm,
        wavelength_medium="air",
        species="H I",
        source_name="synthetic",
        source_reference="tests/test_detector_line_overlay.py",
        source_resource="tests",
    )


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


# ---------------------------------------------------------------- geometry ---


def test_line_lands_on_the_pre_flip_column_of_its_order():
    """690 nm is 20 columns into order 0, not 20 columns from the other end.

    ``Spectrum`` flips its stitched 1-D arrays when wavelength falls with
    column, which is why ``detector_pixel`` is documented as the pre-flip raw
    column.  The image is never flipped, so inverting ``order_wavel`` — which
    is indexed by raw column — must give 20 here.  A mirrored mapping would
    put the same line at column 280.
    """

    geometry = _geometry()
    mark = geometry.mark_for(_line(690.0, "H-test"))

    assert mark is not None
    assert mark.order_index == 0
    assert mark.order == 21
    assert mark.column == pytest.approx(20.0)
    assert mark.column != pytest.approx(COLUMNS - 1 - 20.0)
    # The pattern gives the row band of that order at that column.
    assert mark.row == pytest.approx(32.0)
    assert mark.half_height == HALF_HEIGHT


def test_row_follows_the_slanted_trace_across_the_detector():
    """The same order is marked at different rows at different columns."""

    geometry = _geometry()
    blue_end = geometry.mark_for(_line(560.0))
    red_end = geometry.mark_for(_line(695.0))

    assert blue_end is not None and red_end is not None
    assert blue_end.order_index == red_end.order_index == 0
    assert blue_end.column == pytest.approx(280.0)
    assert red_end.column == pytest.approx(10.0)
    assert blue_end.row == pytest.approx(58.0)
    assert red_end.row == pytest.approx(31.0)


def test_a_wavelength_no_order_carries_is_not_marked():
    geometry = _geometry()
    assert geometry.mark_for(_line(900.0)) is None
    assert geometry.mark_for(_line(200.0)) is None
    assert geometry.wavelength_span() == pytest.approx((270.0, 700.0))


def test_order_borders_decide_who_owns_an_overlapping_line():
    """A line inside the order overlap is marked once, on the stitched owner."""

    overlapping = _line(555.0)
    assert _geometry().mark_for(overlapping).order_index == 0

    # Hand order 0 only its red half, exactly as ``calculate_order_borders``
    # trims the overlap, and the same line belongs to order 1.
    columns = np.arange(COLUMNS)
    borders = np.array(
        [
            columns <= 250,
            np.ones(COLUMNS, dtype=bool),
            np.ones(COLUMNS, dtype=bool),
        ]
    )
    bordered = _geometry(order_borders=borders)
    mark = bordered.mark_for(overlapping)
    assert mark.order_index == 1
    assert mark.order == 22
    assert mark.column == pytest.approx(10.0)


def test_partial_orders_and_broken_calibrations_answer_instead_of_raising():
    """A NaN-padded partial order still places lines; a half-built one is None."""

    wavelengths = _order_wavel()
    wavelengths[2][200:] = np.nan  # the CMOS 28th-order situation
    geometry = DetectorGeometry(_pattern(), wavelengths, order_ids=ORDER_IDS)
    assert geometry.mark_for(_line(310.0)) is None  # off the padded end
    assert geometry.mark_for(_line(370.0)).order_index == 2

    assert DetectorGeometry.from_calibration(None) is None
    assert DetectorGeometry.from_calibration(object()) is None

    class Mismatched:
        pattern = np.zeros((COLUMNS, ORDER_COUNT))
        order_wavel = np.zeros((ORDER_COUNT, COLUMNS + 5))

    assert DetectorGeometry.from_calibration(Mismatched()) is None


# ----------------------------------------------------------------- overlay ---


def test_overlay_is_off_by_default_and_creates_nothing(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    before = len(plot.items)

    overlay = DetectorLineOverlay(plot)
    overlay.set_geometry(_geometry())

    assert all(not overlay.is_family_visible(family) for family in ("balmer", "ne"))
    assert overlay.marks("balmer") == ()
    assert overlay.item("balmer") is None
    assert len(plot.items) == before
    widget.close()


def test_toggling_a_family_draws_and_then_hides_its_marks(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    overlay = DetectorLineOverlay(plot)
    overlay.set_geometry(_geometry())

    drawn = overlay.set_family_visible("balmer", True)
    marks = overlay.marks("balmer")
    assert drawn == len(marks) > 1
    item = overlay.item("balmer")
    assert item is not None and item in plot.items and item.isVisible()

    # Two ticks per line, bracketing the order band without covering it.
    columns, rows = item.getData()
    assert len(columns) == 4 * len(marks)
    for index, mark in enumerate(marks):
        band = rows[4 * index : 4 * index + 4]
        assert band[1] < mark.row < band[2]
        assert band[0] < band[1] and band[2] < band[3]
        assert columns[4 * index] == pytest.approx(mark.column + 0.5)

    labels = {mark.label for mark in marks}
    assert {"H-alpha", "H-beta"} <= labels
    alpha = next(mark for mark in marks if mark.label == "H-alpha")
    assert alpha.order_index == 0
    assert alpha.column == pytest.approx((700.0 - 656.2790) / DISPERSION_NM)

    hidden = overlay.set_family_visible("balmer", False)
    assert hidden == 0
    assert overlay.marks("balmer") == ()
    assert not item.isVisible()
    emptied = item.getData()[0]
    assert emptied is None or len(emptied) == 0
    widget.close()


def test_marks_survive_a_calibration_change_and_die_with_it(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    overlay = DetectorLineOverlay(plot)

    # Enabled before any image is loaded: legal, and simply places nothing.
    assert overlay.set_family_visible("balmer", True) == 0
    assert overlay.item("balmer") is None

    assert overlay.set_geometry(_geometry()) is True
    assert overlay.marks("balmer")

    assert overlay.set_geometry(None) is False
    assert overlay.marks("balmer") == ()
    assert not overlay.item("balmer").isVisible()
    widget.close()


def test_unknown_family_is_named_not_ignored(qt_app):
    widget = pg.PlotWidget()
    overlay = DetectorLineOverlay(widget.getPlotItem())
    with pytest.raises(ValueError, match="unknown line family"):
        overlay.set_family_visible("argon", True)
    widget.close()


def test_dense_lamp_families_stay_bounded(qt_app):
    """A whole-detector view of ThAr must not spray thousands of items."""

    widget = pg.PlotWidget()
    overlay = DetectorLineOverlay(widget.getPlotItem(), max_marks=25)
    overlay.set_geometry(_geometry())
    assert 0 < overlay.set_family_visible("thar", True) <= 25
    widget.close()
