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
from PyQt5 import QtWidgets

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
        name = Path(str(fpth)).name.casefold()
        if name == UNREADABLE_NAME.casefold():
            # Exactly how a Dropbox online-only placeholder answers on Windows.
            raise OSError(22, "Invalid argument", str(fpth))

        if "stranger" in name:
            shape, binning = STRANGER_SHAPE, STRANGER_BINNING
            columns, rows = shape
            images = np.zeros((1, rows // binning[1], columns // binning[0]))
            return images, _info(shape, binning)

        shape = CMOS_SHAPE if "cmos" in name else CCD_SHAPE
        if name.endswith("_bkg.sif") or name.endswith("_bg.sif"):
            images = frames[shape]["background"]
        elif "sphere" in name or "absolute" in name:
            images = frames[shape]["sphere"]
        else:
            images = frames[shape]["shot"]
        return images.copy(), _info(shape)

    monkeypatch.setattr(echelle_module, "read_image", read_image)


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

        columns = shape[0]
        wavelengths = _order_wavelengths(columns)
        rows = []
        for index, order_id in enumerate(ORDER_IDS):
            for pixel in np.linspace(20, columns * 0.8, 5).astype(int):
                rows.append((order_id, pixel, wavelengths[index][pixel]))
        (root / wavelength_name).write_text(
            "# synthetic lamp identification\n"
            "# order from to center wavelength\n"
            + "".join(
                f"{order:d} {pixel - 2:d} {pixel + 2:d} {pixel:d} {value:.6f}\n"
                for order, pixel, value in rows
            ),
            encoding="utf-8",
        )

    micrometres = np.linspace(0.310, 0.850, 40)
    np.savetxt(
        root / "integrating_sphere.txt",
        np.column_stack([micrometres, np.full(micrometres.size, 1.0e-2)]),
        fmt="%.8f",
    )

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
