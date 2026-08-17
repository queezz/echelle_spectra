"""Packet F15 — known-line markers on the 2-D detector image.

The spectrum overlays answer "where is this line in wavelength".  These pin the
geometric half of the same question: with the order pattern and the wavelength
tables loaded, a catalog line has an expected column and an expected row band
on the sensor, and the main GUI's overlay toggles drive both views.

The field fix these also pin (owner, 2026-08-17, on the bright Ne frame): the
mark used to be two ticks *outside* the extraction band with a hairline pen, so
on a PSF 26 px tall against a 17-px band the ticks were painted over the very
blob they marked and vanished.  It is now one box per line, wide as a rough
line width and tall as the order's band.  And a line an order overlap exposes
twice is boxed twice, because the detector really does show it twice.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pyqtgraph as pg
import pytest
from PyQt5 import QtCore, QtWidgets

from echelle_spectra.tools.image_line_overlay import (
    DEFAULT_LINE_WIDTH_PX,
    CursorLink,
    DetectorGeometry,
    DetectorLineOverlay,
    OrderTraceOverlay,
)
from echelle_spectra.tools.line_catalog import SpectralLine
from echelle_spectra.tools.line_overlay import select_overlay_lines

# A three-order sensor with the black Echelle's falling wavelength direction:
# 0.5 nm per column, so order 0 runs 700 nm down to 550 nm and its neighbours
# start 140 nm lower each.  The orders overlap by 10 nm, which is where the
# border mask has to decide who owns a line.
COLUMNS = 301
ORDER_COUNT = 3
DISPERSION_NM = 0.5
ORDER_START_NM = (700.0, 560.0, 420.0)
ORDER_IDS = (21, 22, 23)
#: The extraction half-width ``dv``: the rows the 1-D spectrum is summed over.
HALF_HEIGHT = 8.0
#: The traces sit 20 rows apart, so the *band* each order owns reaches 10 rows
#: either side of its trace — wider than ``dv``, and stopping at the neighbour.
BAND_HALF_HEIGHT = 10.0


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
    assert mark.half_height == BAND_HALF_HEIGHT
    assert mark.half_width == DEFAULT_LINE_WIDTH_PX / 2.0
    assert mark.primary is True


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


def _bordered_geometry() -> DetectorGeometry:
    """Order 0 trimmed to its red half, as ``calculate_order_borders`` trims it."""

    columns = np.arange(COLUMNS)
    borders = np.array(
        [
            columns <= 250,
            np.ones(COLUMNS, dtype=bool),
            np.ones(COLUMNS, dtype=bool),
        ]
    )
    return _geometry(order_borders=borders)


def test_order_borders_decide_which_of_the_two_marks_is_primary():
    """A line in the overlap is boxed on both orders; the stitch owner is primary.

    The overlap is not a bookkeeping wrinkle — the grating really does put 555
    nm on two adjacent orders, and the detector really does show two blobs.  So
    both are marked.  ``order_borders``, the same mask that stitches the 1-D
    spectrum, only decides which of the two the spectrum reads it from, and
    that one is drawn in the family's solid pen.
    """

    overlapping = _line(555.0)
    marks = _bordered_geometry().marks_for_line(overlapping)

    assert [mark.order_index for mark in marks] == [0, 1]
    assert [mark.order for mark in marks] == [21, 22]
    assert [mark.primary for mark in marks] == [False, True]
    # Order 0 shows it 290 columns in; order 1 shows the same light at column 10.
    assert marks[0].column == pytest.approx(290.0)
    assert marks[1].column == pytest.approx(10.0)

    # The primary — what the spectrum agrees with — is the bordered owner.
    primary = _bordered_geometry().mark_for(overlapping)
    assert primary.order_index == 1
    assert primary.order == 22
    assert primary.column == pytest.approx(10.0)

    # Untrimmed, the first carrying order owns it and its twin is secondary.
    untrimmed = _geometry().marks_for_line(overlapping)
    assert [mark.primary for mark in untrimmed] == [True, False]
    assert _geometry().mark_for(overlapping).order_index == 0


def test_a_line_outside_any_overlap_is_boxed_exactly_once():
    """Only the overlap doubles a line; everywhere else there is one box."""

    marks = _bordered_geometry().marks_for_line(_line(690.0))
    assert len(marks) == 1
    assert marks[0].primary is True
    assert marks[0].order_index == 0


def test_the_band_reaches_the_neighbouring_trace_not_just_the_extraction_rows():
    """The box height is the order's band, floored at the extraction ``dv``.

    The field defect: a Ne line runs ~26 px FWHM down the rows while ``dv=8``
    gives a 17-px extraction band, so anything drawn at ``±dv`` sits inside the
    blob.  Half the distance to the neighbouring trace puts the box edge in the
    dark gutter between orders instead.
    """

    geometry = _geometry()
    # Traces 20 rows apart -> the band reaches 10 rows, wider than dv = 8.
    assert geometry.band_half_height(1, 100.0) == pytest.approx(BAND_HALF_HEIGHT)
    assert geometry.band_half_height(1, 100.0) > geometry.half_height
    # The edge orders have one neighbour and still get the same reach.
    assert geometry.band_half_height(0, 100.0) == pytest.approx(BAND_HALF_HEIGHT)
    assert geometry.band_half_height(2, 100.0) == pytest.approx(BAND_HALF_HEIGHT)

    # Crowded traces never let the box reach into the neighbouring order, and
    # never shrink below the rows the spectrum was extracted from.
    crowded = DetectorGeometry(
        np.column_stack([np.full(COLUMNS, 30 + 6 * index) for index in range(3)]),
        _order_wavel(),
        order_ids=ORDER_IDS,
        half_height=HALF_HEIGHT,
    )
    assert crowded.band_half_height(1, 0.0) == pytest.approx(HALF_HEIGHT)

    # A lone order has nothing to measure against and keeps ``dv``.
    single = DetectorGeometry(
        np.full((COLUMNS, 1), 40),
        _order_wavel()[:1],
        order_ids=(21,),
        half_height=HALF_HEIGHT,
    )
    assert single.band_half_height(0, 10.0) == pytest.approx(HALF_HEIGHT)


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


def _boxes(item) -> list[tuple[float, float, float, float]]:
    """Read an item's ``connect="pairs"`` segments back as (l, r, low, high)."""

    columns, rows = item.getData()
    assert len(columns) % 8 == 0
    out = []
    for index in range(len(columns) // 8):
        xs = columns[8 * index : 8 * index + 8]
        ys = rows[8 * index : 8 * index + 8]
        left, right = float(np.min(xs)), float(np.max(xs))
        low, high = float(np.min(ys)), float(np.max(ys))
        # Four edges, and nothing painted across the middle of the box.
        assert set(np.round(xs, 6)) == {round(left, 6), round(right, 6)}
        assert set(np.round(ys, 6)) == {round(low, 6), round(high, 6)}
        out.append((left, right, low, high))
    return out


def test_toggling_a_family_draws_and_then_hides_its_marks(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    overlay = DetectorLineOverlay(plot)
    overlay.set_geometry(_geometry())

    drawn = overlay.set_family_visible("balmer", True)
    marks = overlay.marks("balmer")
    primary = overlay.primary_marks("balmer")
    # The count the status bar reports is lines, not boxes: an overlap twin is
    # a second box on the same line.
    assert drawn == len(primary) > 1
    assert drawn + overlay.duplicate_count("balmer") == len(marks)
    item = overlay.item("balmer")
    assert item is not None and item in plot.items and item.isVisible()

    # One box per mark, wrapping the trace row and spanning the order's band.
    boxes = _boxes(item)
    assert len(boxes) == len(primary)
    for (left, right, low, high), mark in zip(boxes, primary):
        assert left < mark.column + 0.5 < right
        assert low < mark.row + 0.5 < high
        assert right - left == pytest.approx(DEFAULT_LINE_WIDTH_PX)
        assert high - low == pytest.approx(2 * BAND_HALF_HEIGHT)
        # The regression the owner reported: the old mark left the band
        # unpainted and put its ink outside it, which on this instrument's PSF
        # meant painting over the blob.  The box must reach past ``dv``.
        assert high - low > 2 * HALF_HEIGHT

    labels = {mark.label for mark in marks}
    assert {"H-alpha", "H-beta"} <= labels
    alpha = next(mark for mark in marks if mark.label == "H-alpha")
    assert alpha.order_index == 0
    assert alpha.primary is True
    assert alpha.column == pytest.approx((700.0 - 656.2790) / DISPERSION_NM)

    hidden = overlay.set_family_visible("balmer", False)
    assert hidden == 0
    assert overlay.marks("balmer") == ()
    assert not item.isVisible()
    emptied = item.getData()[0]
    assert emptied is None or len(emptied) == 0
    widget.close()


def test_an_overlap_twin_gets_its_own_pooled_item_in_a_secondary_pen(qt_app):
    """Both boxes are drawn, and the duplicate is unmistakably the duplicate."""

    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    overlay = DetectorLineOverlay(plot)
    overlay.set_geometry(_bordered_geometry())

    lines = overlay.set_family_visible("balmer", True)
    duplicates = overlay.duplicate_marks("balmer")
    assert duplicates, "H-delta lands on two orders in this fixture"
    assert lines == len(overlay.primary_marks("balmer"))

    # The same wavelength, two orders, one of each kind.
    doubled = {mark.wavelength_nm for mark in duplicates}
    for wavelength in doubled:
        placed = [mark for mark in overlay.marks("balmer") if mark.wavelength_nm == wavelength]
        assert len(placed) == 2
        assert sorted(mark.primary for mark in placed) == [False, True]
        assert placed[0].order_index != placed[1].order_index

    # Two pooled items for the family, never one per box.
    primary_item = overlay.item("balmer")
    twin_item = overlay.duplicate_item("balmer")
    assert primary_item is not None and twin_item is not None
    assert primary_item is not twin_item
    assert len([item for item in plot.items if item in (primary_item, twin_item)]) == 2
    assert len(_boxes(twin_item)) == len(duplicates)

    # Same hue, secondary style: the twin is dashed, the owner is solid.
    assert twin_item.opts["pen"].style() != QtCore.Qt.SolidLine
    assert primary_item.opts["pen"].style() == QtCore.Qt.SolidLine

    overlay.set_family_visible("balmer", False)
    assert not primary_item.isVisible() and not twin_item.isVisible()
    assert overlay.duplicate_count("balmer") == 0
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
    """A whole-detector view of ThAr must not spray thousands of items.

    The bound is the strength threshold and the curated table, not ``max_marks``
    — a lamp family's selection is a union of the two and a curated row is
    never dropped to fit a cap.  The 11 000-row ThAr table has to come out in
    the hundreds, and it all still costs two pooled items.
    """

    widget = pg.PlotWidget()
    overlay = DetectorLineOverlay(widget.getPlotItem(), max_marks=25)
    overlay.set_geometry(_geometry())
    drawn = overlay.set_family_visible("thar", True)
    assert 0 < drawn <= 400
    assert len([item for item in (overlay.item("thar"), overlay.duplicate_item("thar"))
                if item is not None]) <= 2
    widget.close()


# ------------------------------------------------------------ order traces ---


def test_order_traces_are_off_by_default_and_cost_nothing(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    before = len(plot.items)

    traces = OrderTraceOverlay(plot)
    traces.set_geometry(_geometry())

    assert traces.is_visible is False
    assert traces.item() is None
    assert traces.order_count() == 0
    assert len(plot.items) == before
    widget.close()


def test_order_traces_draw_every_order_as_one_pooled_curve(qt_app):
    """One item, one row per order, NaN gaps so the pen never joins two orders."""

    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    traces = OrderTraceOverlay(plot)
    traces.set_geometry(_geometry())

    assert traces.set_visible(True) == ORDER_COUNT
    item = traces.item()
    assert item is not None and item in plot.items and item.isVisible()
    # One pooled curve for the whole pattern, however many orders there are.
    assert len([other for other in plot.items if other is item]) == 1

    columns, rows = item.getData()
    assert len(columns) == ORDER_COUNT * (COLUMNS + 1)
    pattern = _pattern()
    for index in range(ORDER_COUNT):
        start = index * (COLUMNS + 1)
        segment_x = columns[start : start + COLUMNS]
        segment_y = rows[start : start + COLUMNS]
        assert np.isnan(columns[start + COLUMNS])
        assert np.isnan(rows[start + COLUMNS])
        assert segment_x == pytest.approx(np.arange(COLUMNS) + 0.5)
        assert segment_y == pytest.approx(pattern[:, index] + 0.5)

    assert traces.set_visible(False) == 0
    assert not item.isVisible()
    emptied = item.getData()[0]
    assert emptied is None or len(emptied) == 0
    widget.close()


def test_order_traces_follow_the_calibration_and_die_with_it(qt_app):
    widget = pg.PlotWidget()
    traces = OrderTraceOverlay(widget.getPlotItem())

    # Enabled before any image loads: legal, and simply draws nothing.
    assert traces.set_visible(True) == 0
    assert traces.item() is None

    assert traces.set_geometry(_geometry()) is True
    assert traces.order_count() == ORDER_COUNT
    assert traces.item().isVisible()

    assert traces.set_geometry(None) is False
    assert traces.order_count() == 0
    assert not traces.item().isVisible()
    widget.close()


# -------------------------------------------------------- boxes meet light ---


#: The Ne PSF measured on ``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif``:
#: thin across columns, fat down rows.  It is the row half that broke the
#: shipped brackets, so a synthetic frame has to carry the same shape.
PSF_COLUMN_FWHM = 6.0
PSF_ROW_FWHM = 26.0
#: A window narrow enough that no neighbouring order's blob can enter it: the
#: real CMOS traces are never closer than 59 rows.
CENTROID_ROW_REACH = 28
CENTROID_COLUMN_REACH = 8


def _window(frame, row, column):
    top, bottom = max(0, int(row) - CENTROID_ROW_REACH), int(row) + CENTROID_ROW_REACH + 1
    left = max(0, int(column) - CENTROID_COLUMN_REACH)
    right = int(column) + CENTROID_COLUMN_REACH + 1
    return frame[top:bottom, left:right], top, left


def _centroid(frame, mark):
    """Centre of this line's *own* blob, or ``None`` when it carries no light.

    A catalog line that the lamp does not excite has no blob, and a line with a
    bright neighbour a few columns away has that neighbour's light in its
    window without having any of its own.  Both are answered the same way: the
    window's brightest pixel has to sit at this line's column, or there is no
    blob here to check the box against.
    """

    background = float(np.median(frame))
    window, top, left = _window(frame, mark.row, mark.column)
    window = np.clip(window - background, 0.0, None)
    if window.max() < 400.0:
        return None
    peak_row, peak_column = np.unravel_index(int(np.argmax(window)), window.shape)
    if abs(left + peak_column - mark.column) > mark.half_width:
        return None
    window, top, left = _window(frame, top + peak_row, left + peak_column)
    window = np.clip(window - background, 0.0, None)
    row_grid, column_grid = np.mgrid[top : top + window.shape[0], left : left + window.shape[1]]
    total = window.sum()
    return (
        float((row_grid * window).sum() / total),
        float((column_grid * window).sum() / total),
    )


def _real_ne_case():
    """The owner's bright Ne frame with its calibration, or ``None`` if absent."""

    root = Path(__file__).resolve().parents[1]
    frame_path = root / "local" / "20250926_calib" / "Ne-0.02s-x3-bright-lines.sif"
    calibration_path = root / "src" / "echelle_spectra" / "resources" / "calibration_files"
    if not frame_path.exists():
        return None
    from echelle_spectra.tools import echelle as ech

    calibration = ech.Calibrations(str(calibration_path))
    calibration.filenames = {
        "orders": "pattern_CMOS_20250926.txt",
        "wavelength": "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        "sphr": "sphere_cmos_20240305.sif",
        "bkgr": "sphere_cmos_20240305_bkg.sif",
        "integral": "integrating_sphere.txt",
    }
    try:
        calibration.start()
        image = ech.EchelleImage(str(frame_path), clbr=calibration, spec="black")
    except (OSError, ValueError, KeyError, IndexError):  # pragma: no cover - env
        return None
    geometry = DetectorGeometry.from_calibration(calibration)
    if geometry is None:  # pragma: no cover - a broken calibration folder
        return None
    span = geometry.wavelength_span()
    lines = select_overlay_lines("ne", span[0], span[1], max_labels=200)
    return geometry, np.asarray(image.images[0], dtype=float), lines


# A three-order stand-in for the CMOS: traces 60 rows apart as the real ones
# are, orders overlapping by 20 nm, and blobs painted at positions computed
# *forward* from the dispersion — never by the inversion under test.
SYNTHETIC_COLUMNS = 400
SYNTHETIC_ROWS = (60.0, 120.0, 180.0)
SYNTHETIC_START_NM = (700.0, 640.0, 580.0)
SYNTHETIC_DISPERSION_NM = 0.2
SYNTHETIC_EXCITED_NM = (690.0, 630.0, 610.0, 570.0)


def _synthetic_case():
    pixels = np.arange(SYNTHETIC_COLUMNS, dtype=float)
    order_wavel = np.array(
        [start - SYNTHETIC_DISPERSION_NM * pixels for start in SYNTHETIC_START_NM]
    )
    pattern = np.column_stack(
        [row + pixels // 100 for row in SYNTHETIC_ROWS]
    )
    geometry = DetectorGeometry(
        pattern, order_wavel, order_ids=(7, 8, 9), half_height=HALF_HEIGHT
    )

    truth = []
    for wavelength in SYNTHETIC_EXCITED_NM:
        for index, start in enumerate(SYNTHETIC_START_NM):
            column = (start - wavelength) / SYNTHETIC_DISPERSION_NM
            if 0.0 <= column <= SYNTHETIC_COLUMNS - 1:
                truth.append((wavelength, index, column, float(pattern[int(column), index])))
    assert len({wavelength for wavelength, *_ in truth}) == len(SYNTHETIC_EXCITED_NM)
    assert len(truth) > len(SYNTHETIC_EXCITED_NM), "two of these fall in an overlap"

    frame = np.full((int(pattern.max()) + 60, SYNTHETIC_COLUMNS), 100.0)
    row_axis = np.arange(frame.shape[0])[:, None]
    column_axis = np.arange(SYNTHETIC_COLUMNS)[None, :]
    for _wavelength, _index, column, row in truth:
        frame += 8000.0 * np.exp(
            -(((column_axis - column) / (PSF_COLUMN_FWHM / 2.355)) ** 2) / 2.0
            - (((row_axis - row) / (PSF_ROW_FWHM / 2.355)) ** 2) / 2.0
        )
    return geometry, frame, truth


def test_every_box_for_an_excited_line_contains_its_blobs_centroid():
    """The half of F15 that was never checked against real light: the rows.

    The shipped smoke verified mark *columns* against the stitched spectrum and
    stopped there, so nothing ever asked whether the row band landed on a blob.
    This asks it on the owner's bright Ne frame when the data folder is present,
    and falls back to a synthetic frame carrying the same measured PSF when it
    is not.
    """

    real = _real_ne_case()
    if real is None:
        geometry, frame, truth = _synthetic_case()
        lines = [_line(wavelength) for wavelength in SYNTHETIC_EXCITED_NM]
    else:
        geometry, frame, lines = real
    marks = geometry.marks_for(lines)
    assert marks

    checked = 0
    for mark in marks:
        found = _centroid(frame, mark)
        if found is None:
            continue
        row, column = found
        checked += 1
        assert abs(column - mark.column) <= mark.half_width, (
            f"{mark.label} order {mark.order}: blob centroid column {column:.1f} "
            f"is outside the box at {mark.column:.1f} +- {mark.half_width}"
        )
        assert abs(row - mark.row) <= mark.half_height, (
            f"{mark.label} order {mark.order}: blob centroid row {row:.1f} "
            f"is outside the band at {mark.row:.1f} +- {mark.half_height}"
        )
    assert checked, "no marked line carries light on this frame"

    if real is None:
        # Every painted blob — overlap twins included — must have been boxed.
        for wavelength, index, column, row in truth:
            placed = [
                mark
                for mark in marks
                if mark.wavelength_nm == wavelength and mark.order_index == index
            ]
            assert placed, f"{wavelength} nm on order {index} was not boxed"
            mark = placed[0]
            assert abs(mark.column - column) <= mark.half_width
            assert abs(mark.row - row) <= mark.half_height


# ------------------------------------------------------------ cursor link ---
#
# The owner's second field request: "can we also add a 2D to 1D cursor link,
# toggleable?"  Same geometry, read at the pointer instead of at a catalog row.


def _cursor_link(readout=None, spectra=None):
    """A link over three throwaway plots, with the synthetic geometry loaded."""

    image = pg.PlotWidget()
    counts = pg.PlotWidget()
    calibrated = pg.PlotWidget()
    plots = {
        "counts": counts.getPlotItem(),
        "calibrated": calibrated.getPlotItem(),
    }
    link = CursorLink(image.getPlotItem(), plots, readout=readout)
    link.set_geometry(_geometry(spectra))
    # The widgets are returned only so the caller can keep them alive; a
    # PlotWidget collected mid-test takes its ViewBox with it.
    return link, (image, counts, calibrated)


def test_the_cursor_link_is_off_by_default_and_connects_nothing(qt_app):
    """Mouse moves are the highest-rate event in the window: off must be free."""

    link, widgets = _cursor_link()

    assert link.is_enabled is False
    assert link.proxy_count() == 0
    assert link.image_marker() is None
    assert link.spectrum_marker("counts") is None
    assert link.label == ""
    # And a move that arrives anyway maps nothing while it is off.
    assert link.show_for_image_point(20.0, 32.0) is None

    for widget in widgets:
        widget.close()


def test_the_cursor_link_maps_a_detector_point_to_its_order_and_wavelength(qt_app):
    """Column 20 of order 0 is 690 nm, and the band is the box's own band.

    The same inversion the marks use, read the other way: order 0 starts at
    700 nm and falls half a nanometre per column, so twenty columns in is
    690 nm, and the trace sits at row 32 there.
    """

    said = []
    link, widgets = _cursor_link(readout=said.append)
    link.set_enabled(True)

    assert link.proxy_count() >= 1
    assert link.show_for_image_point(20.0, 32.0) == (21, pytest.approx(690.0))
    assert link.label == "order 21 · 690.00 nm"
    assert said[-1] == "order 21 · 690.00 nm"

    marker = link.spectrum_marker("calibrated")
    assert marker is not None
    assert marker.isVisible()
    assert marker.value() == pytest.approx(690.0)
    # The pointer is already on the image, so nothing is drawn under it there.
    assert link.image_marker() is None or not link.image_marker().isVisible()

    # Anywhere inside the band names the same order; the band reaches the
    # gutter, not merely the extraction half-width.
    assert link.show_for_image_point(20.0, 32.0 + BAND_HALF_HEIGHT - 0.5)[0] == 21

    for widget in widgets:
        widget.close()


def test_a_cursor_outside_every_order_band_marks_nothing(qt_app):
    """Below the lowest trace there is no order, and no answer is the answer."""

    said = []
    link, widgets = _cursor_link(readout=said.append)
    link.set_enabled(True)
    assert link.show_for_image_point(20.0, 32.0) is not None

    assert link.show_for_image_point(20.0, 0.0) is None
    assert not link.spectrum_marker("calibrated").isVisible()
    assert not link.spectrum_marker("counts").isVisible()
    assert link.label == ""
    assert said[-1] == ""

    # Off the end of the sensor is the same answer.
    assert link.show_for_image_point(-4.0, 32.0) is None

    for widget in widgets:
        widget.close()


def test_hovering_a_spectrum_marks_where_the_stitch_read_it(qt_app):
    """The reverse link takes the order the stitched spectrum actually reads.

    690 nm is carried only by order 0 here; 555 nm sits in the 0/1 overlap and
    the border mask decides, so the mark follows the mask rather than the first
    order that happens to reach.
    """

    columns = np.arange(COLUMNS)
    borders = np.array(
        [columns <= 250, np.ones(COLUMNS, dtype=bool), np.ones(COLUMNS, dtype=bool)]
    )
    link, widgets = _cursor_link(spectra=borders)
    link.set_enabled(True)

    order, column, row = link.show_for_wavelength(690.0)
    assert (order, column) == (21, pytest.approx(20.0))
    assert row == pytest.approx(32.0)
    marker = link.image_marker()
    assert marker is not None and marker.isVisible()
    x, y = marker.getData()
    assert x[0] == pytest.approx(20.5 - 16.0)
    assert y[2] == pytest.approx(32.5 - 16.0)

    # 555 nm is on both order 0 and order 1; the mask gives it to order 1.
    order, _column, _row = link.show_for_wavelength(555.0)
    assert order == 22

    # A wavelength no order carries takes both marks down.
    assert link.show_for_wavelength(200.0) is None
    assert not link.image_marker().isVisible()

    for widget in widgets:
        widget.close()


def test_switching_the_cursor_link_off_disconnects_and_hides(qt_app):
    """Off again is off again: no proxies left, no marks left, nothing said."""

    said = []
    link, widgets = _cursor_link(readout=said.append)
    link.set_enabled(True)
    link.show_for_image_point(20.0, 32.0)
    assert link.spectrum_marker("calibrated").isVisible()

    link.set_enabled(False)

    assert link.is_enabled is False
    assert link.proxy_count() == 0
    assert not link.spectrum_marker("calibrated").isVisible()
    assert link.label == ""
    assert said[-1] == ""
    # The pooled items are kept, so switching it on again allocates nothing.
    before = link.spectrum_marker("calibrated")
    link.set_enabled(True)
    link.show_for_image_point(20.0, 32.0)
    assert link.spectrum_marker("calibrated") is before

    for widget in widgets:
        widget.close()
