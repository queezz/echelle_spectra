"""Packet F12 — the main GUI load path always reaches a terminal state.

Every test drives the shipped ``EchelleSpectraGUI`` off-screen: the real
``prepare_calibration``, the real ``LoadCalibrationsThread`` and
``LoadImageThread``, the real camera flip-retry, and the real display sequence.
Only detector pixels are synthetic — ``read_image`` is replaced so the fixture
needs no SIF files — which is exactly the seam Packet F2 established.

The field defect this pins: the owner's own Ne lamp SIF had become a
cloud-placeholder that the reader could not open, the thread reported that with
the same bare ``None`` it used for a size mismatch, and the handler answered
every ``None`` by flipping the camera and loading again — for ever.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from PyQt5 import QtCore, QtWidgets

from echelle_spectra import calibration_bench_gui as bench_gui
from echelle_spectra import echelle_spectra_gui as gui
from echelle_spectra.tools import echelle as echelle_module

# Two detectors that agree about wavelengths and disagree about pixels, so a
# frame from one is a genuine mismatch for the other.
CMOS_SHAPE = (2000, 170)
CCD_SHAPE = (1800, 150)
STRANGER_SHAPE = (900, 90)
STRANGER_BINNING = (2, 2)

DV = 8
ORDER_COUNT = 6
ORDER_IDS = tuple(range(21, 21 + ORDER_COUNT))
# Order 0 is the reddest and sits lowest on the sensor; each order spans 100 nm
# and overlaps its neighbour by 20 nm, and wavelength falls with pixel as it
# does on the black Echelle.  Together they cover every band the GUI plots.
ORDER_START_NM = tuple(830.0 - 80.0 * index for index in range(ORDER_COUNT))
ORDER_SPAN_NM = 100.0
EXPOSURE_S = 0.25
CYCLE_S = 0.5

UNREADABLE_NAME = "6001_unreadable.SIF"

# The owner's real case: 2019 frames read through a 2025 solution sit a rigid
# eighteen pixels away from where the overlay boxes are drawn.
PIXEL_SHIFT = 18
SNAPSHOT_ID = "20190314_cmos"


def _order_wavelengths(columns: int) -> np.ndarray:
    """The wavelength each order carries at each detector column."""

    pixels = np.arange(columns, dtype=float)
    return np.array(
        [start - ORDER_SPAN_NM * pixels / columns for start in ORDER_START_NM]
    )


def _order_centers(shape: tuple) -> np.ndarray:
    """Row of each order at each column; the bluest one climbs off the sensor."""

    columns, rows = shape
    pixels = np.arange(columns)
    centers = [np.full(columns, 20 + 20 * index, dtype=int) for index in range(ORDER_COUNT)]
    # The last order runs out of sensor nine tenths of the way across, which is
    # the shipped 29-order CMOS situation: a partial order, a NaN-padded
    # wavelength solution, and a float (never object) spectra array.
    top = 20 + 20 * (ORDER_COUNT - 1)
    climb = top + (pixels * (rows - DV - top)) // int(columns * 0.9)
    centers[-1] = np.minimum(climb, rows - DV)
    return np.column_stack(centers)


def _sphere_orders(columns: int) -> np.ndarray:
    """One hump per order, crossing its neighbour inside the shared overlap."""

    wavelengths = _order_wavelengths(columns)
    centers = np.array([start - ORDER_SPAN_NM / 2.0 for start in ORDER_START_NM])
    return np.array(
        [
            800.0 * np.exp(-(((lam - center) / 60.0) ** 2)) + 200.0
            for lam, center in zip(wavelengths, centers)
        ]
    )


def _paint(image: np.ndarray, centers: np.ndarray, values: np.ndarray) -> None:
    rows = image.shape[0]
    for column, center in enumerate(centers):
        band = np.arange(center - DV, center + DV + 1)
        inside = (band >= 0) & (band < rows)
        image[band[inside], column] = values[column] / (2 * DV + 1)


def _frame(shape: tuple, values: np.ndarray) -> np.ndarray:
    columns, rows = shape
    image = np.zeros((rows, columns))
    pattern = _order_centers(shape)
    for order in range(ORDER_COUNT):
        _paint(image, pattern[:, order], values[order])
    return image


def _detector_frames(shape: tuple) -> dict:
    columns = shape[0]
    sphere_values = _sphere_orders(columns)
    return {
        "sphere": _frame(shape, sphere_values)[np.newaxis],
        "background": _frame(shape, np.full_like(sphere_values, 40.0))[np.newaxis],
        "shot": _frame(shape, np.full_like(sphere_values, 600.0))[np.newaxis],
    }


def _info(shape: tuple, binning: tuple = (1, 1)) -> dict:
    return {
        "NumberOfFrames": 1,
        "xbin": binning[0],
        "ybin": binning[1],
        "size": np.array([shape[0] // binning[0], shape[1] // binning[1]]),
        "DetectorDimensions": (shape[0], shape[1]),
        "ExposureTime": EXPOSURE_S,
        "CycleTime": CYCLE_S,
        "DetectorTemperature": -30.0,
    }


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


@pytest.fixture
def detector(monkeypatch: pytest.MonkeyPatch) -> None:
    """Serve synthetic frames wherever the GUI would read a SIF file."""

    frames = {CMOS_SHAPE: _detector_frames(CMOS_SHAPE), CCD_SHAPE: _detector_frames(CCD_SHAPE)}

    def read_image(fpth, spec="black", crop=(0, -1), exptime=1):
        path = Path(str(fpth))
        name = path.name.casefold()
        if name == UNREADABLE_NAME.casefold():
            # Exactly how a Dropbox online-only placeholder answers on Windows.
            raise OSError(22, "Invalid argument", str(fpth))

        if "stranger" in name:
            shape, binning = STRANGER_SHAPE, STRANGER_BINNING
            columns, rows = shape
            images = np.zeros((1, rows // binning[1], columns // binning[0]))
            return images, _info(shape, binning)

        # A snapshot names its files by role — ``sphere.sif`` says nothing about
        # the detector — so the folder that holds them answers for them.
        folder = path.parent.name.casefold()
        shape = CMOS_SHAPE if "cmos" in name or "cmos" in folder else CCD_SHAPE
        if name.endswith("_bkg.sif") or name.endswith("_bg.sif"):
            images = frames[shape]["background"]
        elif "sphere" in name or "absolute" in name:
            images = frames[shape]["sphere"]
        else:
            images = frames[shape]["shot"]
        return images.copy(), _info(shape)

    monkeypatch.setattr(echelle_module, "read_image", read_image)


def _wavelength_table(columns: int, pixel_shift: int = 0) -> str:
    """The identified-line table, optionally moved bodily along the detector.

    A shift moves the pixel each known wavelength was identified at, which is
    what a detector that has slid actually does — the same light, further along
    the sensor — so the fitted solution comes out translated by that many pixels
    and everything read through it moves with it.
    """

    wavelengths = _order_wavelengths(columns)
    rows = []
    for index, order_id in enumerate(ORDER_IDS):
        for pixel in np.linspace(20, columns * 0.8, 5).astype(int):
            rows.append((order_id, int(pixel) + pixel_shift, wavelengths[index][pixel]))
    return (
        "# synthetic lamp identification\n"
        "# order from to center wavelength\n"
        + "".join(
            f"{order:d} {pixel - 2:d} {pixel + 2:d} {pixel:d} {value:.6f}\n"
            for order, pixel, value in rows
        )
    )


def _write_sphere_table(path: Path) -> None:
    """The integrating sphere's spectral response, flat across the fixture."""

    micrometres = np.linspace(0.310, 0.850, 40)
    np.savetxt(
        path,
        np.column_stack([micrometres, np.full(micrometres.size, 1.0e-2)]),
        fmt="%.8f",
    )


def _write_calibration_files(root: Path) -> None:
    """Write the calibration tables under the names the GUI asks for."""

    root.mkdir(parents=True, exist_ok=True)
    (root / "alignments").mkdir(exist_ok=True)

    for shape, pattern_name, wavelength_name in (
        (CCD_SHAPE, "pattern.txt", "Th_wavelength.txt"),
        (
            CMOS_SHAPE,
            "pattern_CMOS_20250926.txt",
            "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        ),
    ):
        np.savetxt(root / pattern_name, _order_centers(shape), fmt="%d")
        (root / wavelength_name).write_text(
            _wavelength_table(shape[0]), encoding="utf-8"
        )

    _write_sphere_table(root / "integrating_sphere.txt")

    for name in (
        "absolute_20170613_b8_0.2_v2.sif",
        "absolute_20170613_b8_0.2_bkg.sif",
        "sphere_cmos_20240305.sif",
        "sphere_cmos_20240305_bkg.sif",
    ):
        (root / name).write_bytes(b"served by the synthetic detector")


def _config(tmp_path: Path) -> dict:
    _write_calibration_files(tmp_path / "resources" / "calibration_files")
    return {
        "base_path": tmp_path,
        "data_path": str(tmp_path / "data"),
        "output_path": str(tmp_path / "out"),
        "debug": False,
        "diag_name": "spec_div1",
    }


class LoadSpy:
    """Count the handler calls so an unbounded retry cannot pass unnoticed."""

    def __init__(self, window):
        self.window = window
        self.outcomes = []
        self._handler = window._on_image_loaded
        window._on_image_loaded = self

    def __call__(self, outcome):
        self.outcomes.append(outcome)
        if len(self.outcomes) > 8:
            raise AssertionError(
                "the load path retried without bound: "
                + ", ".join(item.status for item in self.outcomes)
            )
        return self._handler(outcome)

    @property
    def statuses(self):
        return [outcome.status for outcome in self.outcomes]


@pytest.fixture
def window(qt_app, detector, tmp_path: Path):
    win = gui.EchelleSpectraGUI(_config(tmp_path))
    yield win
    win.close()
    qt_app.processEvents()


def _write_snapshot(root: Path, pixel_shift: int = 0, snapshot_id: str = SNAPSHOT_ID) -> Path:
    """Build a real snapshot folder carrying this fixture's CMOS calibration.

    With ``pixel_shift`` zero it is byte-for-byte the calibration the packaged
    CMOS path already reads, which is what makes a shifted one's effect
    attributable to the shift and to nothing else.
    """

    from echelle_spectra.snapshot import create_snapshot

    staging = root / f"sources_{snapshot_id}"
    staging.mkdir(parents=True, exist_ok=True)
    np.savetxt(staging / "pattern.txt", _order_centers(CMOS_SHAPE), fmt="%d")
    (staging / "wavelength.txt").write_text(
        _wavelength_table(CMOS_SHAPE[0], pixel_shift), encoding="utf-8"
    )
    _write_sphere_table(staging / "integral.txt")
    for name in ("sphere_cmos.sif", "sphere_cmos_bg.sif"):
        (staging / name).write_bytes(b"served by the synthetic detector")

    snapshot = create_snapshot(
        root / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files={
            "pattern": staging / "pattern.txt",
            "wavelength": staging / "wavelength.txt",
            "sphere": staging / "sphere_cmos.sif",
            "sphere_background": staging / "sphere_cmos_bg.sif",
            "integral": staging / "integral.txt",
        },
        lamps=["Th"],
        notes="synthetic viewer fixture",
    )
    return snapshot.root


def _snapshot_window(base: Path, pixel_shift: int = 0):
    """A viewer launched the way ``--calibration <folder>`` launches it."""

    config = _config(base)
    root = _write_snapshot(base, pixel_shift)
    override = gui.load_calibration_override(root)
    return gui.EchelleSpectraGUI(config, calibration_override=override)


def _balmer_columns(app, win, name: str) -> dict:
    """Where the Balmer overlay puts each line on the detector image."""

    spy = LoadSpy(win)
    _load(app, win, name, spy, attempts=1)
    assert spy.statuses == [gui.IMAGE_LOADED]
    win.line_overlay_checks["balmer"].setChecked(True)
    marks = win.detector_line_overlays.primary_marks("balmer")
    assert marks
    return {round(mark.wavelength_nm, 3): mark.column for mark in marks}


def test_window_icon_resolves_and_is_set(qt_app):
    """The taskbar/app icon must survive off the installed package, not just app.setWindowIcon.

    Regression check: the window used to rely solely on the QApplication-level
    icon set in ``start()``, so any window built without going through
    ``start()`` (or before it ran) showed no icon at all.
    """

    from echelle_spectra import _config as installed_config

    icon_path = installed_config["base_path"] / "resources/graphics/echelle.png"
    assert icon_path.is_file(), f"icon file missing at {icon_path}"

    win = gui.EchelleSpectraGUI(installed_config)
    try:
        assert not win.windowIcon().isNull()
    finally:
        win.close()
        qt_app.processEvents()


def _pump(app, ready, timeout: float = 120.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        app.processEvents()
        if ready():
            return True
        time.sleep(0.005)
    return False


def _settle_calibrations(app, win) -> None:
    assert _pump(app, win.calibrations_settled), "calibrations never finished loading"
    assert win.calibration_errors == {"CCD": "", "CMOS": ""}


def _load(app, win, name: str, spy: LoadSpy, attempts: int) -> None:
    win.filename = str(Path(win.data_path) / name)
    win.load_image()
    assert _pump(app, lambda: len(spy.outcomes) >= attempts), (
        f"{name} never reached a terminal state; saw {spy.statuses}"
    )


def test_matching_file_loads_on_the_selected_camera(qt_app, window):
    """The ordinary case: the selected camera fits, one attempt, no flip."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)

    _load(qt_app, window, "5001_cmos.SIF", spy, attempts=1)

    assert spy.statuses == [gui.IMAGE_LOADED]
    assert window.cameras_tried == ["CMOS"]
    assert window.CameraCMOS.isChecked()
    assert window.btn_open.isEnabled()
    assert window.spectra.wavelength.size


def test_other_camera_is_tried_exactly_once(qt_app, window):
    """A CCD frame under a CMOS selection costs one flip and then loads."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)

    _load(qt_app, window, "5002_ccd.SIF", spy, attempts=2)

    assert spy.statuses == [gui.IMAGE_DIMENSION_MISMATCH, gui.IMAGE_LOADED]
    assert window.cameras_tried == ["CMOS", "CCD"]
    assert window.CameraCCD.isChecked()
    assert window.btn_open.isEnabled()


def test_file_matching_neither_camera_stops_and_explains(qt_app, window):
    """Both cameras mismatch: two attempts, then a terminal, honest message."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)

    _load(qt_app, window, "5003_stranger.SIF", spy, attempts=2)
    # Give an unbounded retry the room to prove itself before we conclude.
    _pump(qt_app, lambda: len(spy.outcomes) > 2, timeout=2.0)

    assert spy.statuses == [gui.IMAGE_DIMENSION_MISMATCH, gui.IMAGE_DIMENSION_MISMATCH]
    assert sorted(window.cameras_tried) == ["CCD", "CMOS"]

    message = window.image_info_bw.toPlainText()
    assert "5003_stranger.SIF reads 900x90" in message
    assert "binning 2x2" in message
    assert f"CCD calibration expects {CCD_SHAPE[0]}x{CCD_SHAPE[1]}" in message
    assert f"CMOS calibration expects {CMOS_SHAPE[0]}x{CMOS_SHAPE[1]}" in message

    assert "Loading Image" not in window.coursor_bw.toPlainText()
    assert window.btn_open.isEnabled()
    assert window.acceptDrops()


def test_unreadable_file_stops_without_blaming_the_calibration(qt_app, window):
    """The field regression: a file the reader cannot open is never a mismatch.

    Conflating the two is what let the owner's cloud-placeholder SIF ping-pong
    between the cameras behind a permanent "Loading Image".
    """

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)

    _load(qt_app, window, UNREADABLE_NAME, spy, attempts=1)
    _pump(qt_app, lambda: len(spy.outcomes) > 1, timeout=2.0)

    assert spy.statuses == [gui.IMAGE_UNREADABLE]
    # No camera flip: the other calibration cannot open the file either.
    assert window.cameras_tried == ["CMOS"]
    assert window.CameraCMOS.isChecked()

    message = window.image_info_bw.toPlainText()
    assert f"{UNREADABLE_NAME} could not be read" in message
    assert "OSError" in message
    assert "do not match" not in message

    assert "Loading Image" not in window.coursor_bw.toPlainText()
    assert window.btn_open.isEnabled()


def test_load_before_calibrations_finish_is_queued_visibly(qt_app, window):
    """A load issued during startup waits, says so, and then runs by itself."""

    assert not window.calibrations_settled()
    spy = LoadSpy(window)

    window.filename = str(Path(window.data_path) / "5004_cmos.SIF")
    window.load_image()

    # The load is held, the surface says why, and nothing was handed a
    # half-built DIMW/DIMO.
    assert window.pending_image == window.filename
    assert spy.outcomes == []
    assert "Calibrations loading" in window.coursor_bw.toPlainText()
    assert not window.btn_open.isEnabled()

    _settle_calibrations(qt_app, window)
    assert _pump(qt_app, lambda: len(spy.outcomes) >= 1), "the queued load never ran"

    assert spy.statuses == [gui.IMAGE_LOADED]
    assert window.pending_image is None
    assert window.btn_open.isEnabled()


def test_broken_calibration_reports_instead_of_disabling_the_window(qt_app, detector, tmp_path):
    """A calibration that cannot load must not leave the controls dead."""

    config = _config(tmp_path)
    (tmp_path / "resources" / "calibration_files" / "pattern.txt").write_text(
        "not a pattern table\n", encoding="utf-8"
    )

    win = gui.EchelleSpectraGUI(config)
    try:
        assert _pump(qt_app, win.calibrations_settled), "a failed calibration never reported"
        assert win.calibration_errors["CCD"]
        assert win.calibration_errors["CMOS"] == ""
        assert win.btn_open.isEnabled()
        assert "CCD" in win.statusBar().currentMessage()

        # The camera whose calibration died refuses in one bounded step.
        spy = LoadSpy(win)
        win.CameraCCD.setChecked(True)
        win.filename = str(Path(win.data_path) / "5005_ccd.SIF")
        win.load_image()
        assert _pump(qt_app, lambda: len(spy.outcomes) >= 1)
        _pump(qt_app, lambda: len(spy.outcomes) > 1, timeout=2.0)

        assert spy.statuses == [gui.IMAGE_CALIBRATION_UNAVAILABLE]
        assert "CCD calibration is unusable" in win.image_info_bw.toPlainText()
        assert win.btn_open.isEnabled()
    finally:
        win.close()
        qt_app.processEvents()


def test_overlay_toggle_marks_the_detector_image_where_the_spectrum_says(qt_app, window):
    """Packet F15 — one toggle serves both views, and the image agrees with it.

    The marks are placed by inverting each order's wavelength solution against
    the raw detector column.  The stitched spectrum carries the same answer in
    its ``detector_pixel``/``echelle_order`` coordinates — after ``Spectrum``
    has flipped them, because this fixture's wavelength falls with column as
    the black Echelle's does.  Agreement between the two is the flip
    convention holding end to end.

    The spectrum shows an overlapping line once; the image shows both blobs the
    overlap exposes.  So the coordinate the spectrum reports is the coordinate
    of the **primary** box, and the twin is checked against the pattern only.
    """

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5006_cmos.SIF", spy, attempts=1)
    assert spy.statuses == [gui.IMAGE_LOADED]

    overlays = window.detector_line_overlays
    assert overlays.geometry is not None
    # Off by default, and costing nothing while it is off.
    assert not any(overlays.is_family_visible(family) for family in window.line_overlay_checks)
    assert all(overlays.item(family) is None for family in window.line_overlay_checks)
    assert all(overlays.duplicate_item(family) is None for family in window.line_overlay_checks)

    window.line_overlay_checks["balmer"].setChecked(True)

    # One control, both views.
    assert window.line_overlays.is_family_visible("balmer")
    marks = overlays.marks("balmer")
    assert marks
    assert overlays.item("balmer").isVisible()
    assert "on the detector image" in window.statusBar().currentMessage()

    spectra = window.spectra
    for mark in marks:
        # Every box, twin included, sits on its own order's trace.
        column = int(round(mark.column))
        assert mark.row == pytest.approx(
            _order_centers(CMOS_SHAPE)[column, mark.order_index], abs=0.5
        )
        # And its band reaches past the extraction rows toward the neighbour.
        assert mark.half_height >= DV
        if not mark.primary:
            continue
        index = int(np.argmin(np.abs(spectra.wavelength - mark.wavelength_nm)))
        assert spectra.wavelength[index] == pytest.approx(mark.wavelength_nm, abs=0.5)
        assert int(spectra.echelle_order[index]) == mark.order
        assert float(spectra.detector_pixel[index]) == pytest.approx(mark.column, abs=1.0)

    # This fixture's orders overlap by 20 nm, so some line is boxed twice — the
    # blob the shipped single-ownership rule left unmarked.
    duplicates = overlays.duplicate_marks("balmer")
    assert duplicates
    assert overlays.duplicate_item("balmer").isVisible()
    assert "doubled in order overlaps" in window.statusBar().currentMessage()
    for twin in duplicates:
        siblings = [
            mark for mark in overlays.primary_marks("balmer")
            if mark.wavelength_nm == twin.wavelength_nm
        ]
        assert len(siblings) == 1
        assert siblings[0].order_index != twin.order_index

    window.line_overlay_checks["balmer"].setChecked(False)
    assert overlays.marks("balmer") == ()
    assert not overlays.item("balmer").isVisible()
    assert not overlays.duplicate_item("balmer").isVisible()
    assert not window.line_overlays.is_family_visible("balmer")


def test_order_trace_toggle_draws_the_loaded_pattern(qt_app, window):
    """The bench's order traces, on the main window's image, behind one check."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5008_cmos.SIF", spy, attempts=1)

    traces = window.order_trace_overlay
    assert window.order_trace_check.isChecked() is False
    assert traces.is_visible is False
    assert traces.item() is None

    window.order_trace_check.setChecked(True)
    assert traces.is_visible is True
    assert traces.order_count() == ORDER_COUNT
    assert traces.item().isVisible()
    assert "Order traces shown" in window.statusBar().currentMessage()

    columns, rows = traces.item().getData()
    pattern = _order_centers(CMOS_SHAPE)
    width = CMOS_SHAPE[0]
    for order in range(ORDER_COUNT):
        start = order * (width + 1)
        assert rows[start : start + width] == pytest.approx(pattern[:, order] + 0.5)
        assert np.isnan(rows[start + width])

    window.order_trace_check.setChecked(False)
    assert traces.order_count() == 0
    assert not traces.item().isVisible()


def test_failed_load_clears_the_detector_marks(qt_app, window):
    """Marks must not outlive the frame they were placed for."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5007_cmos.SIF", spy, attempts=1)
    window.line_overlay_checks["balmer"].setChecked(True)
    window.order_trace_check.setChecked(True)
    assert window.detector_line_overlays.marks("balmer")
    assert window.order_trace_overlay.order_count() == ORDER_COUNT

    _load(qt_app, window, UNREADABLE_NAME, spy, attempts=2)
    assert window.detector_line_overlays.geometry is None
    assert window.detector_line_overlays.marks("balmer") == ()
    # The pattern belongs to the frame that is gone, too.
    assert window.order_trace_overlay.geometry is None
    assert window.order_trace_overlay.order_count() == 0
    assert not window.order_trace_overlay.item().isVisible()


def test_the_cursor_link_starts_off_and_costs_nothing_until_switched_on(
    qt_app, window
):
    """Mouse moves over the image are the fastest event the window sees.

    So the link connects nothing at all until it is asked for, and switching it
    off again takes the connection back down rather than leaving a slot on the
    scene doing nothing.
    """

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5011_cmos.SIF", spy, attempts=1)

    link = window.cursor_link
    assert not window.cursor_link_check.isChecked()
    assert link.is_enabled is False
    assert link.proxy_count() == 0
    assert link.geometry is not None

    window.cursor_link_check.setChecked(True)
    assert link.is_enabled is True
    assert link.proxy_count() >= 1
    assert "Cursor link on" in window.statusBar().currentMessage()

    # The loaded calibration answers for a point on the reddest order's trace.
    pattern = _order_centers(CMOS_SHAPE)
    column = 100
    reported = link.show_for_image_point(float(column), float(pattern[column, 0]))
    assert reported is not None
    order, wavelength = reported
    assert order == ORDER_IDS[0]
    assert wavelength == pytest.approx(
        _order_wavelengths(CMOS_SHAPE[0])[0][column], abs=0.05
    )
    assert window.statusBar().currentMessage() == f"order {order} · {wavelength:.2f} nm"
    assert link.spectrum_marker("calibrated").isVisible()

    window.cursor_link_check.setChecked(False)
    assert link.proxy_count() == 0
    assert not link.spectrum_marker("calibrated").isVisible()
    assert "Cursor link off" in window.statusBar().currentMessage()


def test_a_hand_set_display_level_survives_a_redraw_of_the_same_frame(
    qt_app, window, monkeypatch
):
    """The rule sets the levels; it never holds them.

    Auto-levelling runs when the displayed frame changes and at no other time,
    so an operator who has dragged the histogram keeps his choice through every
    other redraw of that frame.  What the levels *are* is pinned on real data
    in ``test_detector_display``; what is pinned here is that the viewer asks
    for them exactly once per frame.
    """

    asked = []

    def levels(image, **kwargs):
        asked.append(np.asarray(image).shape)
        return 7.0, 42.0

    monkeypatch.setattr(gui, "auto_display_levels", levels)

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5012_cmos.SIF", spy, attempts=1)

    assert len(asked) == 1
    assert window.hist.getLevels() == pytest.approx((7.0, 42.0))

    window.hist.setLevels(min=5.0, max=17.0)
    window.show_image_frame()
    assert len(asked) == 1, "the same frame was levelled twice"
    assert window.hist.getLevels() == pytest.approx((5.0, 17.0))

    # A genuinely new frame is levelled again.
    window._levelled_frame = None
    window.show_image_frame()
    assert len(asked) == 2
    assert window.hist.getLevels() == pytest.approx((7.0, 42.0))

    # And the checkbox still decides whether the rule runs at all; with it off
    # a new frame falls back to pyqtgraph's own min/max, as it always did.
    window.check_autoscale.setChecked(False)
    window._levelled_frame = None
    window.show_image_frame()
    assert len(asked) == 2
    assert window.hist.getLevels() != pytest.approx((7.0, 42.0))


# ---------------------------------------------------------------------------
# Opening files through a saved calibration snapshot
# ---------------------------------------------------------------------------


def test_launching_without_the_flag_is_the_window_it_has_always_been(qt_app, window):
    """No ``--calibration`` means the packaged pair, the old title, both cameras."""

    args, qt_args = gui.parse_args([])
    assert args.calibration is None
    assert qt_args == []
    # Qt keeps its own single-dash switches, so an off-screen launch still works.
    assert gui.parse_args(["-platform", "offscreen"])[1] == ["-platform", "offscreen"]

    _settle_calibrations(qt_app, window)
    assert window.calibration_override is None
    assert window.windowTitle() == "Echelle viewer"
    assert window.camera_names == ("CCD", "CMOS")
    assert window.CameraCCD.isEnabled()
    assert window.CameraCMOS.isEnabled()
    # The selector is the flag's in-GUI form, so it starts where the flagless
    # launch starts: on the packaged set the checked camera stands for.
    assert window.calibration_select.currentText() == "Packaged CMOS"

    spy = LoadSpy(window)
    _load(qt_app, window, "5101_cmos.SIF", spy, attempts=1)
    assert spy.statuses == [gui.IMAGE_LOADED]
    assert "calibration:" not in window.image_info_bw.toPlainText()


def test_a_snapshot_calibration_is_worn_and_the_window_says_so(qt_app, detector, tmp_path):
    """The viewer opens files through the snapshot it was pointed at, and names it."""

    win = _snapshot_window(tmp_path)
    try:
        assert win.windowTitle() == f"Echelle viewer — {SNAPSHOT_ID}"
        # One snapshot is one detector's calibration: that camera and no other,
        # with nothing packaged left in the other slot to fall back to.
        assert win.camera_names == ("CMOS",)
        assert win.cb_CCD is None
        assert win.CameraCMOS.isChecked()
        assert not win.CameraCCD.isEnabled()
        assert not win.CameraCMOS.isEnabled()

        assert _pump(qt_app, win.calibrations_settled), "the snapshot never loaded"
        assert win.calibration_errors == {"CMOS": ""}
        assert Path(win.cb_CMOS.folder).name == SNAPSHOT_ID
        assert win.cb_CMOS.filenames["orders"] == "pattern.txt"
        assert win.cb_CMOS.filenames["wavelength"] == "wavelength.txt"

        spy = LoadSpy(win)
        _load(qt_app, win, "5102_cmos.SIF", spy, attempts=1)
        assert spy.statuses == [gui.IMAGE_LOADED]
        assert win.cameras_tried == ["CMOS"]
        assert SNAPSHOT_ID in win.image_info_bw.toPlainText()
        assert win.spectra.wavelength.size
    finally:
        win.close()
        qt_app.processEvents()


def test_a_shifted_snapshot_carries_the_overlays_and_the_axis_with_it(
    qt_app, detector, tmp_path, window
):
    """The point of the flag: 2019 frames marked where 2019 says, not 2025.

    The snapshot's wavelength table identifies the same lines eighteen pixels
    further along the sensor than the packaged one does.  Every reader of the
    loaded calibration must move by that much — the image overlay, the cursor
    link, and the spectrum's wavelength axis — because they all read the one
    calibration the window loaded and nothing else.
    """

    _settle_calibrations(qt_app, window)
    packaged = _balmer_columns(qt_app, window, "5103_cmos.SIF")
    window.cursor_link_check.setChecked(True)
    pattern = _order_centers(CMOS_SHAPE)
    column = 100
    packaged_reading = window.cursor_link.show_for_image_point(
        float(column), float(pattern[column, 0])
    )
    packaged_red_edge = float(np.nanmax(window.spectra.wavelength))

    shifted = _snapshot_window(tmp_path / "shifted", pixel_shift=PIXEL_SHIFT)
    try:
        assert _pump(qt_app, shifted.calibrations_settled), "the snapshot never loaded"
        assert shifted.calibration_errors == {"CMOS": ""}

        moved = _balmer_columns(qt_app, shifted, "5103_cmos.SIF")
        shared = sorted(set(packaged) & set(moved))
        assert shared, "the two calibrations agree about no line at all"
        for wavelength in shared:
            assert moved[wavelength] == pytest.approx(
                packaged[wavelength] + PIXEL_SHIFT, abs=0.5
            )

        # The same eighteen pixels, in nanometres, under the pointer …
        per_pixel_nm = ORDER_SPAN_NM / CMOS_SHAPE[0]
        shifted.cursor_link_check.setChecked(True)
        shifted_reading = shifted.cursor_link.show_for_image_point(
            float(column), float(pattern[column, 0])
        )
        assert shifted_reading[0] == packaged_reading[0]
        assert shifted_reading[1] - packaged_reading[1] == pytest.approx(
            PIXEL_SHIFT * per_pixel_nm, abs=0.02
        )

        # … and on the spectrum's own wavelength axis.
        shifted_red_edge = float(np.nanmax(shifted.spectra.wavelength))
        assert shifted_red_edge - packaged_red_edge == pytest.approx(
            PIXEL_SHIFT * per_pixel_nm, abs=0.1
        )
    finally:
        shifted.close()
        qt_app.processEvents()


def test_a_folder_that_is_not_a_snapshot_is_refused_before_a_window_opens(
    tmp_path, capsys
):
    """Bounded and explanatory: never a wedge, never the packaged tables instead."""

    plain = tmp_path / "just_a_folder"
    plain.mkdir()

    with pytest.raises(gui.CalibrationOverrideError) as absent:
        gui.load_calibration_override(tmp_path / "nowhere")
    assert "is not a directory" in str(absent.value)

    with pytest.raises(gui.CalibrationOverrideError) as bare:
        gui.load_calibration_override(plain)
    assert "is not a calibration snapshot" in str(bare.value)

    # And the entry point says it once and leaves, rather than opening a window
    # wearing the packaged tables the operator did not ask for.
    assert gui.start(["--calibration", str(plain)]) == 2
    assert "is not a calibration snapshot" in capsys.readouterr().err


def test_a_snapshot_missing_a_table_names_the_missing_table(tmp_path):
    """A snapshot that lost a derived file fails validation, saying which."""

    root = _write_snapshot(tmp_path)
    (root / "wavelength.txt").unlink()

    with pytest.raises(gui.CalibrationOverrideError) as err:
        gui.load_calibration_override(root)

    message = str(err.value)
    assert "did not validate" in message
    assert "wavelength.txt" in message


# ---------------------------------------------------------------------------
# The calibration selector: the same thing, chosen in the window
# ---------------------------------------------------------------------------


def _pick(select, label: str) -> None:
    """Choose one calibration entry the way the operator chooses it."""

    index = select.findText(label)
    assert index >= 0, f"no calibration entry named {label!r}"
    assert index != select.currentIndex(), f"{label!r} was already chosen"
    select.setCurrentIndex(index)


def _entry_labels(select) -> list:
    return [select.itemText(index) for index in range(select.count())]


def _balmer_marks(win) -> dict:
    """Where the Balmer overlay currently puts each line on the detector image."""

    marks = win.detector_line_overlays.primary_marks("balmer")
    assert marks
    return {round(mark.wavelength_nm, 3): mark.column for mark in marks}


def test_the_selector_names_the_packaged_calibration_the_window_wears(qt_app, window):
    """It displays what is worn — including when the flip-retry moves it."""

    _settle_calibrations(qt_app, window)
    select = window.calibration_select

    assert _entry_labels(select) == ["Packaged CCD", "Packaged CMOS", gui.BROWSE_LABEL]
    assert select.currentText() == "Packaged CMOS"
    assert select.isEnabled()
    assert window.recent_snapshots == []

    # A CCD frame flips the camera, which changes what the file is being read
    # through; a selector that kept saying CMOS would be lying.
    spy = LoadSpy(window)
    _load(qt_app, window, "5201_ccd.SIF", spy, attempts=2)
    assert window.CameraCCD.isChecked()
    assert select.currentText() == "Packaged CCD"

    # And under the packaged pair the packaged entries mean the radio buttons.
    _pick(select, "Packaged CMOS")
    assert window.CameraCMOS.isChecked()
    assert window.calibration_override is None
    assert window.windowTitle() == "Echelle viewer"


def test_picking_a_snapshot_folder_in_the_window_rebases_the_open_frame(
    qt_app, window, tmp_path, monkeypatch
):
    """The owner's ask: open the data, then change whose eyes it is read with.

    Nothing about the flag's outcome may depend on having used the flag — the
    live pick must move the overlays, the wavelength axis, the camera, and the
    title exactly as a ``--calibration`` launch does, and it must move the frame
    that is already on screen.
    """

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5202_cmos.SIF", spy, attempts=1)
    window.line_overlay_checks["balmer"].setChecked(True)
    packaged = _balmer_marks(window)
    packaged_red_edge = float(np.nanmax(window.spectra.wavelength))

    root = _write_snapshot(tmp_path / "picked", pixel_shift=PIXEL_SHIFT)
    monkeypatch.setattr(gui, "choose_snapshot_folder", lambda parent, start: str(root))

    _pick(window.calibration_select, gui.BROWSE_LABEL)

    # The switch rides the ordinary calibration threads, so the whole F12
    # startup contract applies to it: busy controls, and the open frame held.
    assert not window.calibrations_settled()
    assert not window.calibration_select.isEnabled()
    assert window.pending_image == window.filename

    assert _pump(qt_app, window.calibrations_settled), "the snapshot never loaded"
    assert window.calibration_errors == {"CMOS": ""}
    assert _pump(qt_app, lambda: len(spy.outcomes) >= 2), "the open frame was never re-read"
    assert spy.statuses == [gui.IMAGE_LOADED, gui.IMAGE_LOADED]

    assert window.calibration_override.snapshot_id == SNAPSHOT_ID
    assert window.windowTitle() == f"Echelle viewer — {SNAPSHOT_ID}"
    assert window.calibration_select.currentText() == SNAPSHOT_ID
    assert SNAPSHOT_ID in window.image_info_bw.toPlainText()
    assert Path(window.cb_CMOS.folder).name == SNAPSHOT_ID

    # The detector mode followed the snapshot, and both surfaces say why.
    assert window.camera_names == ("CMOS",)
    assert window.CameraCMOS.isChecked()
    assert not window.CameraCCD.isEnabled()
    assert not window.CameraCMOS.isEnabled()
    assert "fixed to CMOS" in window.calibration_select.toolTip()
    assert "fixed while it is worn" in window.CameraCMOS.toolTip()

    # And every reader of the loaded calibration moved the snapshot's eighteen
    # pixels with it — the image overlay …
    moved = _balmer_marks(window)
    shared = sorted(set(packaged) & set(moved))
    assert shared, "the two calibrations agree about no line at all"
    for wavelength in shared:
        assert moved[wavelength] == pytest.approx(
            packaged[wavelength] + PIXEL_SHIFT, abs=0.5
        )

    # … and the spectrum's own wavelength axis.
    per_pixel_nm = ORDER_SPAN_NM / CMOS_SHAPE[0]
    shifted_red_edge = float(np.nanmax(window.spectra.wavelength))
    assert shifted_red_edge - packaged_red_edge == pytest.approx(
        PIXEL_SHIFT * per_pixel_nm, abs=0.1
    )


def test_switching_back_to_a_packaged_calibration_restores_the_radios(
    qt_app, window, tmp_path, monkeypatch
):
    """Off is as live as on: both packaged sets back, and the flip-retry with them."""

    _settle_calibrations(qt_app, window)
    root = _write_snapshot(tmp_path / "worn")
    monkeypatch.setattr(gui, "choose_snapshot_folder", lambda parent, start: str(root))

    _pick(window.calibration_select, gui.BROWSE_LABEL)
    assert _pump(qt_app, window.calibrations_settled), "the snapshot never loaded"
    assert window.calibration_override is not None
    assert window.cb_CCD is None

    _pick(window.calibration_select, "Packaged CCD")
    assert _pump(qt_app, window.calibrations_settled), "the packaged pair never came back"

    assert window.calibration_override is None
    assert window.camera_names == ("CCD", "CMOS")
    assert window.calibration_errors == {"CCD": "", "CMOS": ""}
    assert window.windowTitle() == "Echelle viewer"
    assert window.CameraCCD.isEnabled()
    assert window.CameraCMOS.isEnabled()
    assert window.CameraCCD.isChecked()
    assert window.CameraCCD.toolTip() == "Calibration for CCD camera selected"
    assert window.calibration_select.currentText() == "Packaged CCD"
    assert "calibration:" not in window.image_info_bw.toPlainText()

    # The snapshot stays one click away for the rest of the session — and only
    # for the session: nothing about it is written down.
    assert _entry_labels(window.calibration_select) == [
        "Packaged CCD",
        "Packaged CMOS",
        SNAPSHOT_ID,
        gui.BROWSE_LABEL,
    ]

    # The radios mean what they always meant: a CMOS frame costs one flip.
    spy = LoadSpy(window)
    _load(qt_app, window, "5203_cmos.SIF", spy, attempts=2)
    assert spy.statuses == [gui.IMAGE_DIMENSION_MISMATCH, gui.IMAGE_LOADED]
    assert window.CameraCMOS.isChecked()
    assert window.calibration_select.currentText() == "Packaged CMOS"


def test_a_refused_folder_leaves_the_viewer_wearing_what_it_wore(
    qt_app, window, tmp_path, monkeypatch
):
    """A failed pick must never strand the window uncalibrated, and must say why.

    The flag refuses before a window exists; the selector refuses with one in
    front of it, so the answer is the same sentence and the calibration already
    on is the calibration still on.
    """

    _settle_calibrations(qt_app, window)

    plain = tmp_path / "not_a_snapshot"
    plain.mkdir()
    monkeypatch.setattr(gui, "choose_snapshot_folder", lambda parent, start: str(plain))
    _pick(window.calibration_select, gui.BROWSE_LABEL)

    with pytest.raises(gui.CalibrationOverrideError) as bare:
        gui.load_calibration_override(plain)
    assert window.statusBar().currentMessage() == str(bare.value)
    assert "is not a calibration snapshot" in window.statusBar().currentMessage()

    # A snapshot that lost a table is refused the same way.
    gutted = _write_snapshot(tmp_path / "gutted", snapshot_id="20190315_cmos")
    (gutted / "wavelength.txt").unlink()
    monkeypatch.setattr(gui, "choose_snapshot_folder", lambda parent, start: str(gutted))
    _pick(window.calibration_select, gui.BROWSE_LABEL)

    with pytest.raises(gui.CalibrationOverrideError) as missing:
        gui.load_calibration_override(gutted)
    assert window.statusBar().currentMessage() == str(missing.value)
    assert "wavelength.txt" in window.statusBar().currentMessage()

    # Nothing was taken off, nothing was half put on, nothing was remembered.
    assert window.calibration_override is None
    assert window.calibrations_settled()
    assert window.calibration_errors == {"CCD": "", "CMOS": ""}
    assert window.windowTitle() == "Echelle viewer"
    assert window.calibration_select.currentText() == "Packaged CMOS"
    assert window.calibration_select.isEnabled()
    assert _entry_labels(window.calibration_select) == [
        "Packaged CCD",
        "Packaged CMOS",
        gui.BROWSE_LABEL,
    ]

    # And the viewer is still a viewer.
    spy = LoadSpy(window)
    _load(qt_app, window, "5204_cmos.SIF", spy, attempts=1)
    assert spy.statuses == [gui.IMAGE_LOADED]


def test_a_load_during_a_calibration_switch_is_queued_never_raced(
    qt_app, window, tmp_path, monkeypatch
):
    """F12's guard holds across a switch, because the switch reuses it."""

    _settle_calibrations(qt_app, window)
    spy = LoadSpy(window)
    _load(qt_app, window, "5205_cmos.SIF", spy, attempts=1)

    root = _write_snapshot(tmp_path / "midswitch")
    monkeypatch.setattr(gui, "choose_snapshot_folder", lambda parent, start: str(root))
    _pick(window.calibration_select, gui.BROWSE_LABEL)

    # Mid-switch the calibrations are in the air again, so everything that
    # could start a load against half-built DIMW/DIMO is visibly dead.
    assert not window.calibrations_settled()
    assert not window.btn_open.isEnabled()
    assert not window.show_btn.isEnabled()
    assert not window.calibration_select.isEnabled()
    assert not window.acceptDrops()

    # A load issued anyway is queued against the calibration coming up.
    window.filename = str(Path(window.data_path) / "5206_cmos.SIF")
    window.load_image()
    assert window.pending_image == window.filename
    assert len(spy.outcomes) == 1

    # And a second switch attempted mid-switch is refused out loud rather than
    # stacked on top of the threads already running.
    _pick(window.calibration_select, "Packaged CCD")
    assert "still loading" in window.statusBar().currentMessage()
    assert window.calibration_override.snapshot_id == SNAPSHOT_ID
    assert window.calibration_select.currentText() == SNAPSHOT_ID

    assert _pump(qt_app, window.calibrations_settled), "the snapshot never loaded"
    assert _pump(qt_app, lambda: len(spy.outcomes) >= 2), "the queued load never ran"
    # Give a ping-pong the room to prove itself before we conclude there is none.
    _pump(qt_app, lambda: len(spy.outcomes) > 2, timeout=2.0)

    assert spy.statuses == [gui.IMAGE_LOADED, gui.IMAGE_LOADED]
    assert window.pending_image is None
    assert Path(window.filename).name == "5206_cmos.SIF"
    assert window.cameras_tried == ["CMOS"]
    assert window.btn_open.isEnabled()
    assert window.calibration_select.isEnabled()


def test_the_viewer_snapshot_dialog_follows_a_pasted_path_like_the_bench_does(
    qt_app, tmp_path, monkeypatch
):
    """Owner, 2026-08-18: "I want a file dialog which allows me to paste in the
    folder for preview. Now I can just paste it for selection, not going near."

    The viewer and the calibration bench ask the same question with what is now
    the same dialog, so this is the viewer's half of one change: a pasted path
    that names a real folder navigates there, and the greyed contents that make
    a snapshot folder recognisable — its ``snapshot.toml`` — are on screen
    before Choose is pressed rather than never.
    """

    first = _write_snapshot(tmp_path / "first")
    second = _write_snapshot(tmp_path / "second")
    opened = []
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "exec_", lambda dialog: (opened.append(dialog) or 0)
    )

    assert gui.choose_snapshot_folder(None, first) == ""
    dialog = opened[0]
    # Qt's own dialog, files listed and greyed: the same configuration the
    # bench gets, because it is now literally the same code.
    assert dialog.fileMode() == QtWidgets.QFileDialog.Directory
    assert not dialog.testOption(QtWidgets.QFileDialog.ShowDirsOnly)
    assert dialog.testOption(QtWidgets.QFileDialog.DontUseNativeDialog)
    assert Path(dialog.directory().absolutePath()) == first

    edit = dialog.findChildren(QtWidgets.QLineEdit)[0]
    edit.setText(str(second))
    qt_app.processEvents()

    assert Path(dialog.directory().absolutePath()) == second
    assert "snapshot.toml" in set(dialog.directory().entryList(QtCore.QDir.Files))
    assert edit.text() == str(second)


def test_the_viewer_snapshot_dialog_reopens_where_it_was_last_left(
    qt_app, tmp_path, monkeypatch
):
    """The viewer walks back to the same snapshot shelf all day.

    The configured data directory is the right answer for the first opening
    and the wrong one for every opening after it. The bench's own pickers keep
    their own places rather than sharing this one.
    """

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    chosen = _write_snapshot(tmp_path / "shelf" / "20250926_cmos")
    opened = []

    def answer(dialog):
        opened.append(Path(dialog.directory().absolutePath()))
        dialog.setDirectory(str(chosen))
        dialog.selectFile(str(chosen))
        return QtWidgets.QDialog.Accepted

    monkeypatch.setattr(QtWidgets.QFileDialog, "exec_", answer)

    assert Path(gui.choose_snapshot_folder(None, data_dir)) == chosen
    gui.choose_snapshot_folder(None, data_dir)
    assert opened == [data_dir, chosen]

    # Each question keeps its own memory: the bench's calibration-folder
    # picker must not be walked to the viewer's snapshot shelf.
    bench_gui.choose_calibration_folder(None, data_dir)
    assert opened[-1] == data_dir
