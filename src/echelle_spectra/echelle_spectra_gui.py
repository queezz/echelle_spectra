import argparse
import ctypes
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from pyqtgraph.Qt import QtCore, QtGui

from . import __version__, _config
from .lhd_text import render_lhd_header, write_lhd_text
from .resources import window_layout
from .tools import echelle as ech
from .tools import emissionbands as eb
from .tools import emissiondata as ebd
from .tools.detector_display import auto_display_levels
from .tools.image_line_overlay import (
    CursorLink,
    DetectorLineOverlay,
    OrderTraceOverlay,
)
from .tools.line_overlay import SPECTRUM_CURVE_COLORS, LineOverlayManager

# What one camera's attempt at one file ended up seeing.  A single "it did not
# work" sentinel used to stand for all of these, so a file the reader could not
# open was indistinguishable from a file the calibration did not fit — and the
# camera flip-retry chased the wrong cure forever.
IMAGE_LOADED = "loaded"
IMAGE_DIMENSION_MISMATCH = "dimension-mismatch"
IMAGE_UNREADABLE = "unreadable"
IMAGE_EXTRACTION_FAILED = "extraction-failed"
IMAGE_DISPLAY_FAILED = "display-failed"
IMAGE_CALIBRATION_UNAVAILABLE = "calibration-unavailable"

CAMERA_NAMES = ("CCD", "CMOS")

# What one entry of the calibration selector means.  The packaged entries are
# the CCD/CMOS radio buttons said in words; a snapshot entry is a folder already
# opened this session; the browse entry asks for a new one.
CALIBRATION_PACKAGED = "packaged"
CALIBRATION_SNAPSHOT = "snapshot"
CALIBRATION_BROWSE = "browse"

BROWSE_LABEL = "Snapshot folder…"


def packaged_calibration_label(camera):
    """How one packaged camera set is named in the selector."""
    return f"Packaged {camera}"


def choose_snapshot_folder(parent, start_dir):
    """Ask the operator for a calibration snapshot folder.

    Its own function so the in-GUI selector has one seam a test can stand in
    for: an off-screen run must never put a real modal dialog on the screen.
    """
    return QtWidgets.QFileDialog.getExistingDirectory(
        parent,
        "Open files through a saved calibration snapshot",
        str(start_dir),
    )


class CalibrationOverrideError(ValueError):
    """A folder the viewer was pointed at cannot serve as its calibration."""


@dataclass(frozen=True)
class CalibrationOverride:
    """One saved snapshot standing in for the packaged calibration tables.

    ``filenames`` is the ordinary ``Calibrations`` vocabulary the snapshot
    itself hands out, so the override travels through exactly the loader the
    packaged CCD/CMOS sets travel through — only the folder and the names in it
    differ.
    """

    snapshot_id: str
    camera: str
    folder: Path
    filenames: dict


def load_calibration_override(folder):
    """Resolve a snapshot folder into the single calibration the viewer wears.

    Every failure is named and raised here, before a window exists: an operator
    who asked for a particular era's eyes must never be handed the packaged
    tables instead without being told.
    """
    from .snapshot import SnapshotError, load_snapshot

    path = Path(folder).expanduser()
    if not path.is_dir():
        raise CalibrationOverrideError(
            f"--calibration wants a calibration snapshot folder; {path} is not a directory"
        )
    if not (path / "snapshot.toml").is_file():
        raise CalibrationOverrideError(
            f"{path} is not a calibration snapshot: it holds no snapshot.toml"
        )
    try:
        # Structure and presence only, no digest pass: the viewer is a reader,
        # not an auditor, and a thin snapshot's referenced raws can be hundreds
        # of megabytes on a NAS — hashing them before a window opens is a
        # launch-time cost `echelle snapshot validate` already owns.
        snapshot = load_snapshot(path, verify_files=False)
    except SnapshotError as err:
        raise CalibrationOverrideError(
            f"calibration snapshot {path} did not validate: {err}"
        ) from err

    camera = snapshot.detector.strip().upper()
    if camera not in CAMERA_NAMES:
        raise CalibrationOverrideError(
            f"calibration snapshot {snapshot.snapshot_id!r} names detector "
            f"{snapshot.detector!r}, which is not one of {', '.join(CAMERA_NAMES)}"
        )
    return CalibrationOverride(
        snapshot_id=snapshot.snapshot_id,
        camera=camera,
        folder=snapshot.root,
        filenames=snapshot.calibration_files(),
    )


@dataclass(frozen=True)
class ImageLoadOutcome:
    """One camera's answer for one file, carrying everything a message needs."""

    status: str
    camera: str
    image: object = None
    file_dimensions: tuple = ()
    binning: tuple = ()
    expected_dimensions: tuple = ()
    detail: str = ""


@dataclass(frozen=True)
class CalibrationLoadOutcome:
    """One calibration set's answer, so a broken file cannot hang the startup."""

    name: str
    calibration: object
    detail: str = ""


def calibration_dimensions(calibration):
    """Return a calibration's expected (width, height), empty until it loads."""
    width = getattr(calibration, "DIMW", None)
    height = getattr(calibration, "DIMO", None)
    if width is None or height is None:
        return ()
    return (int(width), int(height))


def detector_dimensions(info):
    """Read the detector size a frame reports, the way calibrations derive theirs."""
    reported = info.get("DetectorDimensions")
    if reported is not None:
        return tuple(int(value) for value in reported)
    width, height = info["size"]
    return (int(width) * int(info["xbin"]), int(height) * int(info["ybin"]))


def format_dimensions(dimensions):
    """Render a (width, height) pair for the operator."""
    if not len(dimensions):
        return "unknown size"
    return "{}x{}".format(*dimensions)


class EchelleSpectraGUI(QMainWindow, window_layout.Ui_MainWindow):
    """GUI window for the echelle_spectra app"""

    def __init__(self, config, calibration_override=None):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.setWindowIcon(
            QtGui.QIcon(str(config["base_path"] / "resources/graphics/echelle.png"))
        )

        # Which calibration this window wears.  ``None`` is the packaged
        # CCD/CMOS pair — what every launch without --calibration gets — and a
        # snapshot narrows the window to that snapshot's one detector.
        self.calibration_override = calibration_override
        self.camera_names = tuple(CAMERA_NAMES)
        # The title before any snapshot named itself in it, so a second switch
        # replaces the first snapshot's name instead of queueing behind it.
        self.base_title = self.windowTitle()
        # Snapshot folders opened through the selector this session, most recent
        # last.  Deliberately not persisted: the next launch starts packaged.
        self.recent_snapshots = []
        # True while the selector is being set to agree with what is worn, so
        # the display update is not mistaken for the operator asking to switch.
        self._selector_syncing = False

        # set widget statuses
        self.CameraCCD.setChecked(False)
        self.CameraCMOS.setChecked(True)
        # The overlay palette is chosen to contrast with these two, so they are
        # named once and read from there rather than written twice.
        self.spec_counts = self.p2.plot(pen=SPECTRUM_CURVE_COLORS["counts"])
        self.spec_wm = self.p3.plot(pen=SPECTRUM_CURVE_COLORS["calibrated"])
        self.line_overlays = LineOverlayManager(max_labels=14)
        self.line_overlays.register_plot("counts", self.p2, labels=False)
        self.line_overlays.register_plot("calibrated", self.p3, labels=True)
        # The same five toggles also drive the 2-D image: one control per
        # family, both views.  Nothing is drawn until a family is switched on.
        self.detector_line_overlays = DetectorLineOverlay(self.p1)
        # And its own toggle for the order pattern underneath them.
        self.order_trace_overlay = OrderTraceOverlay(self.p1)
        # The same geometry read at the pointer instead of at a catalog row.
        # Nothing is connected until its checkbox is ticked.
        self.cursor_link = CursorLink(
            self.p1,
            {"counts": self.p2, "calibrated": self.p3},
            readout=self.statusBar().showMessage,
        )
        # Which (file, frame) the image levels were last set for, so a frame
        # the operator has already levelled by hand is not re-levelled under
        # him on an unrelated redraw.
        self._levelled_frame = None

        # define initial class attributes
        self.config = config
        self.frame_spinners = [self.frame_civ, self.frame_he, self.frame_h]
        self.cameras = [self.CameraCCD, self.CameraCMOS]
        # What the camera buttons say when nothing has fixed them, so switching
        # back off a snapshot puts their own words back rather than a blank.
        self.camera_tooltips = [button.toolTip() for button in self.cameras]

        # carry out init actions
        self.connect_actions()
        self.update_paths()
        self.read_last_shot()
        self._build_calibration_selector()
        self.prepare_calibration()
        self.setup_bands()

        # define auxiliary class attributes
        self.shot_dict = {}
        self.shot_range = []
        self.frame_current = 0
        self.in_loop = False

    def change_camera(self):
        """Change the active camera"""
        for cam in self.cameras:
            if not cam.isChecked():
                cam.setChecked(True)
                break

    def prepare_calibration(self):
        """Select calibration files here"""
        files_ccd = {
            "orders": "pattern.txt",
            "wavelength": "Th_wavelength.txt",
            "sphr": "absolute_20170613_b8_0.2_v2.sif",
            "bkgr": "absolute_20170613_b8_0.2_bkg.sif",
            "integral": "integrating_sphere.txt",
        }

        files_cmos = {
            "orders": "pattern_CMOS_20250926.txt",
            "wavelength": "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
            "sphr": "sphere_cmos_20240305.sif",
            "bkgr": "sphere_cmos_20240305_bkg.sif",
            "integral": "integrating_sphere.txt",
        }

        self.cb_CCD = ech.Calibrations(self.path_calibration)
        self.cb_CCD.name = "CCD"
        self.cb_CMOS = ech.Calibrations(self.path_calibration)
        self.cb_CMOS.name = "CMOS"

        self.cb_CCD.filenames = files_ccd
        self.cb_CMOS.filenames = files_cmos

        if self.calibration_override is not None:
            self._wear_calibration_override()

        self.calib_threads = {name: None for name in self.camera_names}
        # One entry per camera once its thread reports: "" when the calibration
        # loaded, otherwise why it did not.  Until both are in, no image load
        # may run against half-built DIMW/DIMO.
        self.calibration_errors = {}
        self.pending_image = None
        self.cameras_tried = []

        for calibration in self.active_calibrations():
            self.load_calibration(calibration)

        self._sync_calibration_selector()

    # -----------------------------------------------------------------------
    #                        The calibration selector
    # -----------------------------------------------------------------------

    def _build_calibration_selector(self):
        """Fill the control column's calibration combo and listen to it.

        The camera buttons are listened to as well, because under the packaged
        pair they *are* the choice: a flip-retry that moves the camera moves
        what the window is reading through, and the selector must say so.
        """
        self._populate_calibration_selector()
        self.calibration_select.currentIndexChanged.connect(
            self._on_calibration_selected
        )
        for button in self.cameras:
            button.toggled.connect(self._sync_calibration_selector)

    def _populate_calibration_selector(self):
        """Rebuild the entries: the packaged sets, this session's snapshots, browse."""
        syncing, self._selector_syncing = self._selector_syncing, True
        try:
            self.calibration_select.clear()
            for camera in CAMERA_NAMES:
                self.calibration_select.addItem(
                    packaged_calibration_label(camera),
                    (CALIBRATION_PACKAGED, camera),
                )
            for override in self.recent_snapshots:
                self.calibration_select.addItem(
                    override.snapshot_id, (CALIBRATION_SNAPSHOT, override.snapshot_id)
                )
            self.calibration_select.addItem(BROWSE_LABEL, (CALIBRATION_BROWSE, None))
        finally:
            self._selector_syncing = syncing

    def _worn_calibration_key(self):
        """The selector entry that names what the window is reading through now"""
        if self.calibration_override is None:
            selected = [name for name in CAMERA_NAMES if self._camera_button(name).isChecked()]
            camera = selected[0] if selected else CAMERA_NAMES[0]
            return (CALIBRATION_PACKAGED, camera)
        return (CALIBRATION_SNAPSHOT, self.calibration_override.snapshot_id)

    def _calibration_entry_index(self, key):
        """Which entry carries this key, or -1.

        Read out by hand rather than through ``findData``: the entries carry
        plain Python tuples, and Qt's own search compares them as variants.
        """
        for index in range(self.calibration_select.count()):
            if self.calibration_select.itemData(index) == key:
                return index
        return -1

    def _sync_calibration_selector(self, *_):
        """Show what is worn, without that display counting as a request."""
        key = self._worn_calibration_key()
        index = self._calibration_entry_index(key)
        syncing, self._selector_syncing = self._selector_syncing, True
        try:
            if index >= 0:
                self.calibration_select.setCurrentIndex(index)
        finally:
            self._selector_syncing = syncing
        self.calibration_select.setToolTip(self._calibration_tooltip())

    def _calibration_tooltip(self):
        """Say, durably, which calibration every file is being read through."""
        override = self.calibration_override
        if override is None:
            return (
                "Which calibration every file is read through: the packaged CCD\n"
                "and CMOS tables, or a saved snapshot folder"
            )
        return (
            f"Reading every file through calibration snapshot {override.snapshot_id}\n"
            f"— a {override.camera} calibration, so the camera is fixed to "
            f"{override.camera}"
        )

    def _camera_button(self, camera):
        return getattr(self, "Camera" + camera)

    def _on_calibration_selected(self, index):
        """Act on the operator's pick: packaged set, known snapshot, or browse."""
        if self._selector_syncing:
            return

        data = self.calibration_select.itemData(index)
        if not data:
            return
        kind, key = data

        if not self.calibrations_settled():
            # The combo is greyed while calibrations load, so this is the
            # belt to that brace: a switch is refused, said, and undone rather
            # than dropped onto half-built DIMW/DIMO.
            self.statusBar().showMessage(
                "Calibrations are still loading — the calibration selector is "
                "not available yet."
            )
            self._sync_calibration_selector()
            return

        if kind == CALIBRATION_BROWSE:
            self._browse_calibration_snapshot()
            return

        if kind == CALIBRATION_PACKAGED:
            self._wear_packaged_calibration(key)
            return

        if kind == CALIBRATION_SNAPSHOT:
            if (
                self.calibration_override is not None
                and self.calibration_override.snapshot_id == key
            ):
                return
            for override in self.recent_snapshots:
                if override.snapshot_id == key:
                    self._rebase_calibration(override)
                    return

    def _browse_calibration_snapshot(self):
        """Ask for a snapshot folder and put it on, or leave everything as it is."""
        override = self.calibration_override
        start = Path(override.folder).parent if override is not None else self.data_path
        folder = choose_snapshot_folder(self, start)
        if not folder:
            # Cancelled: the combo had already moved to the browse entry, so put
            # it back on what the window is still wearing.
            self._sync_calibration_selector()
            return
        self.wear_snapshot_folder(folder)

    def wear_snapshot_folder(self, folder):
        """Read every file through this snapshot from now on; refuse and stay put if it is not one.

        The same resolver the ``--calibration`` flag uses, so a folder refused
        at launch is refused here for the same stated reason — the difference
        being that a window already exists, and it keeps the calibration it has
        rather than being stranded uncalibrated.
        """
        try:
            override = load_calibration_override(folder)
        except CalibrationOverrideError as err:
            self.statusBar().showMessage(str(err))
            self._sync_calibration_selector()
            return False

        self._rebase_calibration(override)
        return True

    def _wear_packaged_calibration(self, camera):
        """Go back to the packaged pair, with this camera selected."""
        if self.calibration_override is None:
            # Already packaged: the entry means exactly what the radio means.
            self._camera_button(camera).setChecked(True)
            return
        self._rebase_calibration(None, camera=camera)

    def _rebase_calibration(self, override, camera=None):
        """Put on a different calibration live, through the ordinary load path.

        Nothing here loads a table itself: ``prepare_calibration`` rebuilds the
        set and starts the same ``LoadCalibrationsThread`` startup starts, so
        ``calibrations_settled`` goes false for the duration and every guard
        F12 built — greyed controls, the queued image, the bounded retry —
        applies to a mid-session switch exactly as it applies to a launch.
        """
        open_file = getattr(self, "filename", None)

        self.calibration_override = override
        self.camera_names = tuple(CAMERA_NAMES)
        self._release_camera_buttons()
        if override is None:
            self.setWindowTitle(self.base_title)

        self.prepare_calibration()

        if override is None and camera is not None:
            self._camera_button(camera).setChecked(True)

        if open_file:
            # The frame on screen was extracted with the calibration just taken
            # off.  Hand it to the same queue a load issued during startup uses,
            # and it is re-read through the new one the moment that one is up.
            self.pending_image = open_file

        self._sync_calibration_selector()

    def _release_camera_buttons(self):
        """Give the camera buttons their ordinary meaning back."""
        for button, tooltip in zip(self.cameras, self.camera_tooltips):
            button.setEnabled(True)
            button.setToolTip(tooltip)

    def _remember_snapshot(self, override):
        """Keep this session's snapshot folders one click away in the combo."""
        known = [item.snapshot_id for item in self.recent_snapshots]
        if override.snapshot_id in known:
            return
        self.recent_snapshots.append(override)
        self._populate_calibration_selector()

    def _wear_calibration_override(self):
        """Put the snapshot's tables where the packaged ones would have gone.

        A snapshot is one detector's calibration, so the window keeps exactly
        that camera and empties the other slot: there is nothing honest for a
        size-mismatch flip to fall back to, and the packaged tables must not be
        what it silently finds.  The folder and the role filenames come straight
        from the snapshot, and the ordinary ``Calibrations`` loader reads them
        the same way it reads ``resources/calibration_files``.
        """
        override = self.calibration_override
        self.camera_names = (override.camera,)
        self._remember_snapshot(override)

        calibration = ech.Calibrations(str(override.folder))
        calibration.name = override.camera
        calibration.filenames = dict(override.filenames)
        setattr(self, "cb_" + override.camera, calibration)
        for name in CAMERA_NAMES:
            if name != override.camera:
                setattr(self, "cb_" + name, None)

        for button, name in zip(self.cameras, CAMERA_NAMES):
            button.setChecked(name == override.camera)
            # Fixed, not merely preselected: the snapshot decided the detector
            # when it was made, and no other calibration is loaded to switch to.
            button.setEnabled(False)
            button.setToolTip(
                f"Calibration snapshot {override.snapshot_id} is a {override.camera} "
                "calibration, so the camera is fixed while it is worn"
            )

        # Built from the bare title, never from the current one: switching
        # snapshots twice must replace the name, not accumulate names.
        self.setWindowTitle("{} — {}".format(self.base_title, override.snapshot_id))

    def active_calibrations(self):
        """The calibrations this window actually wears, in camera order"""
        return [getattr(self, "cb_" + name) for name in self.camera_names]

    def connect_actions(self):
        self.frame.valueChanged.connect(self.show_image_frame)
        self.frame_h.valueChanged.connect(self.show_balmer_frame)
        self.frame_he.valueChanged.connect(self.show_he_frame)
        self.frame_civ.valueChanged.connect(self.show_c_frame)
        self.btn_open.clicked.connect(self.openfile)
        self.shot_range_btn.clicked.connect(self.get_shot_range)
        self.show_btn.clicked.connect(self.load_shot_number)
        self.start_register.clicked.connect(self.start_loop)
        self.abort_register.clicked.connect(self.abort_loop)
        self.btn_save_cube.clicked.connect(self.save_spectrocube)
        for family, checkbox in self.line_overlay_checks.items():
            checkbox.toggled.connect(
                lambda visible, key=family: self.set_line_overlay(key, visible)
            )
        self.order_trace_check.toggled.connect(self.set_order_traces)
        self.cursor_link_check.toggled.connect(self.set_cursor_link)

    def set_line_overlay(self, family, visible):
        """Toggle one known-line family on both the spectra and the image."""
        self.line_overlays.set_family_visible(family, visible)
        lines = self.detector_line_overlays.set_family_visible(family, visible)
        doubled = self.detector_line_overlays.duplicate_count(family)
        label = self.line_overlay_checks[family].text()
        if not visible:
            self.statusBar().showMessage(f"{label} lines hidden.")
        elif self.detector_line_overlays.geometry is None:
            self.statusBar().showMessage(
                f"{label} lines shown; the detector image is marked once an image loads."
            )
        else:
            # An order overlap exposes the same line twice, and the image boxes
            # both — so say how many lines, then how many are doubled.
            twins = f", {doubled} doubled in order overlaps" if doubled else ""
            self.statusBar().showMessage(
                f"{label} lines shown — {lines} on the detector image{twins}."
            )

    def set_order_traces(self, visible):
        """Toggle the calibration's order pattern over the detector image."""
        orders = self.order_trace_overlay.set_visible(visible)
        if not visible:
            self.statusBar().showMessage("Order traces hidden.")
        elif self.order_trace_overlay.geometry is None:
            self.statusBar().showMessage(
                "Order traces shown; the pattern is drawn once an image loads."
            )
        else:
            self.statusBar().showMessage(
                f"Order traces shown — {orders} orders on the detector image."
            )

    def set_cursor_link(self, enabled):
        """Toggle the pointer link between the detector image and the spectra."""
        ready = self.cursor_link.set_enabled(enabled)
        if not enabled:
            self.statusBar().showMessage("Cursor link off.")
        elif not ready:
            self.statusBar().showMessage(
                "Cursor link on; it starts reporting once an image loads."
            )
        else:
            self.statusBar().showMessage(
                "Cursor link on — hover the image or a spectrum."
            )

    def update_paths(self):
        """Set-up file paths for calibration, data folder, output folder etc."""
        for path in ["data_path", "output_path"]:
            setattr(
                self,
                path,
                Path(
                    self.config[path]
                    .replace("{homedir}", str(Path.home()))
                    .replace("{workdir}", str(Path().absolute()))
                ),
            )
            getattr(self, path).mkdir(parents=True, exist_ok=True)

        self.path_calibration = self.config["base_path"] / "resources/calibration_files"
        self.path_last_shot = self.config["base_path"] / ".last_shot"

        if self.config["debug"]:
            print(self.config)

    def write_last_shot(self):
        """Write last loaded shot number to file"""
        try:
            with open(self.path_last_shot, "w") as f:
                f.write("{} # shot number".format(self.shot_number.value()))
        except OSError:
            pass

    def read_last_shot(self):
        """Read last loaded shot number from file"""
        try:
            with open(self.path_last_shot, "r") as f:
                for line in f.readlines():
                    if "shot number" in line:
                        shot = int(line.split("#")[0])
            self.shot_number.setValue(shot)
        except OSError:
            pass

    def load_calibration(self, clbr):
        """Load calibration files and prepare Calibrations thread
        calibration - calibration instance for CCD or CMOS
        """
        self.statusBar().showMessage("Loading calibration files")
        self._enable_controls(False)
        self.coursor_bw.setText(
            """<font size = 6 color = "#d1451b">Loading Calibrations</font>"""
        )

        self.calib_threads[clbr.name] = LoadCalibrationsThread(clbr, self.config)
        self.calib_threads[clbr.name].taskFinished.connect(self._on_calibration_loaded)
        self.calib_threads[clbr.name].start()

    def calibrations_settled(self):
        """True once every calibration thread has reported, loaded or failed"""
        return set(self.camera_names) <= set(self.calibration_errors)

    def _on_calibration_loaded(self, outcome):
        """Record one calibration; release buttons only once both have reported"""
        if outcome.name == "CCD":
            self.cb_CCD = outcome.calibration

        if outcome.name == "CMOS":
            self.cb_CMOS = outcome.calibration

        self.calibration_errors[outcome.name] = outcome.detail
        if self.config["debug"]:
            print(outcome.name, outcome.detail or "loaded")

        if not self.calibrations_settled():
            self.statusBar().showMessage(
                f"{outcome.name} calibration ready. Still loading the other camera."
            )
            return

        self.coursor_bw.setText("")
        self._enable_controls(True)
        broken = [name for name in self.camera_names if self.calibration_errors.get(name)]
        if broken:
            self.statusBar().showMessage(
                "Calibration files failed to load — "
                + "; ".join(f"{name}: {self.calibration_errors[name]}" for name in broken)
            )
        else:
            self.statusBar().showMessage(self._calibration_ready_message())

        self._load_pending_image()

    def _calibration_ready_message(self):
        """Say what is now worn, and — for a snapshot — why the camera is fixed"""
        override = self.calibration_override
        if override is None:
            return "Calibration files loaded. Ready to work."
        return (
            f"Calibration snapshot {override.snapshot_id} loaded — a "
            f"{override.camera} calibration, so the camera is fixed to "
            f"{override.camera}."
        )

    def _load_pending_image(self):
        """Run a load that arrived while the calibrations were still loading"""
        pending = self.pending_image
        if pending is None:
            return

        self.pending_image = None
        self.filename = pending
        self.load_image()

    def setup_bands(self):
        """Create list of band objects from EmissionBand class"""
        self.hbands = [
            band.copy() for band in [ebd.halpha, ebd.hbeta, ebd.hgamma, ebd.hdelta]
        ]
        self.hebands = [
            band.copy()
            for band in [
                ebd.he447,
                ebd.he492,
                ebd.he501,
                ebd.he504,
                ebd.he587,
                ebd.he667,
                ebd.he706,
                ebd.he728,
            ]
        ]
        self.cbands = [
            band.copy()
            for band in [
                ebd.c515,
                ebd.c464,
                ebd.c580,
                ebd.c444,
                ebd.c465,
                ebd.c772,
                ebd.c706,
                ebd.c547,
            ]
        ]
        self.chband = eb.EB(name="CHD-band", bounds=[429, 431.5])

        self.bands = self.hbands + self.hebands + self.cbands + [self.chband]
        self.bandstofit = self.cbands + self.hebands

        # txt = '<font size = 6 color = "#d1451b">{}</font>'.format(txt)
        for band in self.bands:
            self.cmd_bw.append(band.report_html())

    # ===========================================================================
    #                        End of Initialization
    # ===========================================================================

    # ===========================================================================
    #                        Loop through shots
    # ===========================================================================

    """
    CAUTION
    self.advance_if_inloop() - exit points from the loop.
    If placed in wrong places, the loop will behave unpredictably.

    Exit points:

    _onSpecSaved
    no_fit_intesities
    _onFitEnd
    """

    def start_loop(self):
        """Start cycle trough shots in range(self.start_shot, self.end_shot)"""
        self.shot_range = range(self.start_shot.value(), self.end_shot.value() + 1)
        self.shot_range_index = 0
        self.update_shot_dict()

        self.shot_range = [
            shot for shot in self.shot_range if shot in list(self.shot_dict)
        ]
        if not len(self.shot_range):
            self.coursor_bw.setText(
                '<font size = 6 color = "#e5c40d">Empty Range</font>'
            )
            return

        self.progress_range.setRange(0, len(self.shot_range))
        self.shot_range_len = len(self.shot_range)
        self.abort_register.setEnabled(True)
        self.in_loop = True
        self.loop_step()

    def loop_step(self):
        """Update current shot number and run load_shotnumber"""
        self.shot_number.setValue(self.shot_range[0])
        self.shot_range = self.shot_range[1:]
        self.progress_range.setValue(self.shot_range_index)
        self.shot_range_index += 1
        self.load_shot_number()

    def loop_advance(self):
        """If in the loop, go to next shot"""
        if len(self.shot_range):
            self.loop_step()
        else:
            # end of the loop is here
            self.progress_range.setValue(self.shot_range_index)
            txt = '<font size = 6 color = "#d1451b">Loop\'s Done<br>{}/{}</font>'
            txt = txt.format(self.shot_range_index, self.shot_range_len)
            self.coursor_bw.setText(txt)
            self._enable_controls(True)
            self.in_loop = False

    def advance_if_in_loop(self):
        """If in the loop, go to next shot"""
        if self.in_loop:
            self.loop_advance()

    def abort_loop(self):
        """Abort the loop"""
        if len(self.shot_range):
            self.shot_range = []
            self.coursor_bw.append(
                '<font size = 6 color = "#d1451b">Aborting Loop</font>'
            )
            self.abort_register.setEnabled(False)

    # ===========================================================================
    #                   ^^^  Loop through shots ^^^
    # ===========================================================================

    # ===========================================================================
    #                       Open shot image or file
    # ===========================================================================

    def openfile(self):
        """Open file dialog for selection of SIF file to load"""
        self.filename = QtWidgets.QFileDialog.getOpenFileName(
            None, "Open file to plot", str(self.data_path), "*.sif;;*.SIF;;*.*"
        )[0]
        if self.filename:
            self.emitted()

    def update_shot_dict(self):
        """Convert list of *_Echelle.sif files into shot list and make a dictionary of {shot_number: path to shot}"""
        self.shot_dict = {
            int(f.stem.split("_")[0]): f for f in self.data_path.rglob("*_Echelle.SIF")
        }

    def get_shot_range(self):
        """Get first and last shot numbers from the data folder after indexing it"""
        self.update_shot_dict()
        self.start_shot.setValue(min(self.shot_dict))
        self.end_shot.setValue(max(self.shot_dict))
        if self.config["debug"]:
            print("checking shot_dict min max")
            print(min(self.shot_dict))
            print(max(self.shot_dict))

    def load_shot_number(self):
        """Try to load image of the selected shot number from the data folder"""
        self.update_shot_dict()
        try:
            self.filename = self.shot_dict[self.shot_number.value()]
            self.load_image()
        except KeyError:
            self.statusBar().showMessage(
                f"Failed to load shot {self.shot_number.value()}"
            )

    # ===========================================================================
    #                 ^^^   Open shot image or file   ^^^
    # ===========================================================================

    def load_image(self):
        """Load SIF image, giving each camera calibration at most one attempt"""
        if not self.calibrations_settled():
            self._queue_image_load()
            return

        self.cameras_tried = []
        self._start_image_load()

    def _queue_image_load(self):
        """Hold a load issued during startup instead of racing the calibrations"""
        self.pending_image = self.filename
        self._enable_controls(False)
        self.coursor_bw.setText(
            """<font size = 6 color = "#d1451b">Calibrations loading…</font>"""
        )
        self.statusBar().showMessage(
            f"Calibrations are still loading — {Path(self.filename).name} is queued."
        )

    def _start_image_load(self):
        """Run one attempt with the camera calibration currently selected"""
        if not self.in_loop:
            self.progress_range.setRange(0, 1)
            self.clear_fit_traces()

        self.progress_bands.setRange(0, 1)
        self._enable_controls(False)

        self.fres = {band.name: eb.FitResult(band.name) for band in self.bands}
        self.statusBar().showMessage(f"Loading SIF image: {self.filename}")

        # set shot number from the image filename
        try:
            shot = int(Path(self.filename).stem.split("_")[0])
        except ValueError:
            shot = 0
        self.shot_number.setValue(shot)
        self.write_last_shot()
        self._enable_controls(False)
        self.coursor_bw.setText('<font size = 6 color = "#d1451b">Loading Image</font>')

        if self.calibration_override is not None:
            # The snapshot names the detector; the radio buttons only report it.
            cb = getattr(self, "cb_" + self.calibration_override.camera)
        elif self.CameraCMOS.isChecked():
            cb = self.cb_CMOS
            if self.config["debug"]:
                print("CMOS selected")
        elif self.CameraCCD.isChecked():
            cb = self.cb_CCD
            if self.config["debug"]:
                print("CCD selected")
        else:
            cb = self.cb_CCD
            self.CameraCCD.setChecked(True)

        self.cameras_tried.append(cb.name)
        self.image_load_thread = LoadImageThread(self.filename, cb, self.config)
        self.image_load_thread.taskFinished.connect(self._on_image_loaded)
        self.image_load_thread.start()

    def _on_image_loaded(self, outcome):
        """Receive the loaded Echelle Image"""
        if outcome.status != IMAGE_LOADED:
            self._on_image_load_failed(outcome)
            return

        self.em = outcome.image
        try:
            self._show_loaded_image()
        except Exception as err:
            # An unhandled failure here used to escape the slot and leave the
            # window on its orange "Loading Image" with every control dead.
            self._on_image_load_failed(
                ImageLoadOutcome(
                    status=IMAGE_DISPLAY_FAILED,
                    camera=outcome.camera,
                    file_dimensions=outcome.file_dimensions,
                    binning=outcome.binning,
                    expected_dimensions=outcome.expected_dimensions,
                    detail=f"{type(err).__name__}: {err}",
                )
            )
            return

        if not self.in_loop:
            self._enable_controls(True)
        self.statusBar().showMessage(f"SIF image loaded: {self.filename}")
        self.coursor_bw.setText(
            f'<font size = 5 color = "#187031">{Path(self.filename).stem}</font>'
        )

        if self.specsave_bx.isChecked():
            self.spectra.shotnumber = int(self.shot_number.value())
            self.spectra.output_path = self.output_path
            self.spectra.trigdelay = self.trigger_delay.value()

            self.save_spec_thread = SaveSpectraThread(self.spectra)
            self.save_spec_thread.pass_result.connect(self._on_spec_saved)
            self.save_spec_thread.start()
        else:
            self._on_spec_saved(None)

    def _show_loaded_image(self):
        """Build the spectrum and paint every tab from the freshly loaded image"""
        self.spectra = ech.Spectrum(self.em)
        # The calibration this frame was extracted with is the one that says
        # where a wavelength lands on the sensor.
        self.detector_line_overlays.set_geometry(getattr(self.em, "clbr", None))
        self.order_trace_overlay.set_geometry(getattr(self.em, "clbr", None))
        self.cursor_link.set_geometry(getattr(self.em, "clbr", None))

        self._reset_frame()
        self._setup_frame()

        self.show_info()
        self.show_image_frame()
        self.show_c_frame()
        self.show_he_frame()
        self.show_balmer_frame()

    def _on_image_load_failed(self, outcome):
        """Try the other camera once on a size mismatch, then stop and explain

        Only a dimension mismatch is a question the other calibration can
        answer.  An unreadable file, a failed extraction, or a calibration that
        never loaded fail identically on both cameras, so retrying them only
        buys another spin of the same wheel.
        """
        untried = [name for name in self.camera_names if name not in self.cameras_tried]

        if outcome.status == IMAGE_DIMENSION_MISMATCH and untried:
            self.statusBar().showMessage(
                f"{outcome.camera} calibration expects "
                f"{format_dimensions(outcome.expected_dimensions)} — trying {untried[0]}."
            )
            self.change_camera()
            self._start_image_load()
            return

        self.em = None
        # No image, no geometry: stale marks must not outlive the frame they
        # were placed for.
        self.detector_line_overlays.set_geometry(None)
        self.order_trace_overlay.set_geometry(None)
        self.cursor_link.set_geometry(None)
        self._levelled_frame = None
        message = self._load_failure_text(outcome)
        self.image_info_bw.setText(
            '<font size = 4 color = "#d1451b">{}</font>'.format(message)
        )
        self.coursor_bw.setText('<font size = 5 color = "#d1451b">Load failed</font>')
        self.statusBar().showMessage(message)
        if not self.in_loop:
            self._enable_controls(True)

        # A file that cannot be loaded must not strand a running shot loop.
        self.advance_if_in_loop()

    def _load_failure_text(self, outcome):
        """Say what the file read and what each calibration expected"""
        name = Path(self.filename).name

        if outcome.status == IMAGE_DIMENSION_MISMATCH:
            expectations = "; ".join(self._camera_expectations())
            return (
                f"{name} reads {format_dimensions(outcome.file_dimensions)} "
                f"(binning {format_dimensions(outcome.binning)}); {expectations}."
            )

        if outcome.status == IMAGE_CALIBRATION_UNAVAILABLE:
            return f"The {outcome.camera} calibration is unusable: {outcome.detail}"

        if outcome.status == IMAGE_EXTRACTION_FAILED:
            return (
                f"{name} fits the {outcome.camera} calibration but could not be "
                f"extracted: {outcome.detail}"
            )

        if outcome.status == IMAGE_DISPLAY_FAILED:
            return (
                f"{name} loaded against the {outcome.camera} calibration but could "
                f"not be displayed: {outcome.detail}"
            )

        return f"{name} could not be read: {outcome.detail}"

    def _camera_expectations(self):
        """One phrase per camera describing the frame size it can accept"""
        phrases = []
        for calibration in self.active_calibrations():
            dimensions = calibration_dimensions(calibration)
            if not len(dimensions):
                phrases.append(f"the {calibration.name} calibration is unavailable")
            else:
                phrases.append(
                    f"{calibration.name} calibration expects "
                    f"{format_dimensions(dimensions)}"
                )
        return phrases

    def _on_spec_saved(self, result):
        """Continue after spectra saved (or not)"""
        if result is not None:
            self.spectra = result

        if self.spectra.info["NumberOfFrames"] > 1:
            self.no_fit_intesities()
        else:
            # if in loop saving spectra for 1 frame images
            self.advance_if_in_loop()

    def no_fit_intesities(self):
        """Put the spectra data into the FitResult object and integrate experimental spectra"""
        full_frames = np.arange(self.spectra.info["NumberOfFrames"])
        frames = np.arange(self.spectra.info["NumberOfFrames"])

        for band in self.bands:
            # calculate intensities without fit for all bands for a start
            self.fres[band.name].set_frames(frames)
            self.fres[band.name].init_spectra()
            self.fres[band.name].set_info(self.em.info["NumberOfFrames"])

            for frame in frames:
                x, y = eb.banddata(band, frame, self.spectra, doblin=False)
                self.fres[band.name].spectra.append(y)
            self.fres[band.name].x = x

            self.fres[band.name].integrate_spectra()
            # fill intensities_combined with non-fited results
            self.fres[band.name].integrate_fit()
            self.fres[band.name].fill_nans()

        # plot other traces
        for band in self.hbands:
            self.htraces[band.name].setData(
                x=full_frames, y=self.fres[band.name].intensities
            )
        self.chtrace.setData(x=full_frames, y=self.fres[self.chband.name].intensities)

        # also plot traces for bandstofit, without a fit for a start
        self.show_fit_intensities()

        if self.fit_lines_bx.isChecked():
            self.fit_lines()
        else:
            if self.save_lines_bx.isChecked():
                self.save_intensities()

            self.advance_if_in_loop()  # if saving data in the loop, no fitting

    def fit_lines(self):
        """Fit lines"""
        self.statusBar().showMessage("Fitting Lines")
        self._enable_controls(False)
        self.progress_bands.setRange(1, self.spectra.info["NumberOfFrames"] - 1)

        # array of frames for which fit will be attempted
        frames = np.arange(1, self.em.info["NumberOfFrames"])

        for band in self.bandstofit:
            self.fres[band.name].set_frames(frames)
            self.fres[band.name].set_info(self.em.info["NumberOfFrames"])

        self.fit_thread = FitLinesThread(
            self.spectra, self.bandstofit, self.fres, frames
        )
        self.fit_thread.pass_result.connect(self._on_fit_end)
        self.fit_thread.progress.connect(self._on_fit_progress)
        self.fit_thread.start()

    def _on_fit_progress(self, value):
        """Triggered by individual completed fits"""
        self.progress_bands.setValue(value)

    def _on_fit_end(self, result):
        """Triggered once fitting is completed"""
        self.fres = result
        self.statusBar().showMessage("Fit is Done")
        if not self.in_loop:
            self._enable_controls(True)
        self.progress_bands.setRange(0, 1)

        self.show_fit_intensities()
        self.show_c_frame()
        self.show_he_frame()

        if self.save_lines_bx.isChecked():
            self.save_intensities()

        # go to next shot if in loop.
        self.advance_if_in_loop()

    def save_intensities(self):
        """Save intensities"""
        # make numpy array for no-fits
        names = [b.name for b in self.bands]
        frames = np.arange(self.em.info["NumberOfFrames"])

        # self.trigger_delay.setValue(2.5) # to change the default value if needed
        times = self.trigger_delay.value() + frames * self.em.info["CycleTime"]
        data = np.array(
            [times] + [self.fres[b.name].intensities_combined for b in self.bands]
        )
        data = data.T[1:]  # [1:] to remove frame 0, where values are np.nan

        shot = self.shot_number.value()
        pth = self.output_path / f"{shot}_{self.config['diag_name']}.txt"

        if Path(pth).is_file() and not self.overwrite.isChecked():
            return

        header = render_lhd_header(
            diagnostic=self.config["diag_name"],
            shot=shot,
            dimension_name="Time",
            dimension_size=data.shape[0],
            dimension_unit="s",
            value_names=names,
            value_units=["W/m^2"] * len(names),
            trigger_delay_s=self.trigger_delay.value(),
            frame_interval_s=self.em.info["CycleTime"],
        )
        write_lhd_text(
            pth,
            data,
            header=header,
            formats=["%.5f"] + ["%.6e"] * len(names),
            overwrite=self.overwrite.isChecked(),
        )

    def save_spectrocube(self):
        """Save the currently loaded spectrum as a SpectroCube NetCDF (.nc) file.

        Uses the units currently selected in the units combo-box.
        Calls export_spectrocube() — no conversion logic is duplicated here.
        """
        if not hasattr(self, "spectra") or self.spectra is None:
            self.statusBar().showMessage("No spectrum available to export.")
            return

        from .tools.spectrocube_export import export_spectrocube

        units = str(self.spec_units.currentText())
        shot = self.shot_number.value()
        default_name = f"{shot}_spectrocube.nc"

        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save SpectroCube",
            str(self.output_path / default_name),
            "SpectroCube NetCDF (*.nc);;All files (*)",
        )
        if not path:
            return

        self.statusBar().showMessage(f"Saving SpectroCube: {path} …")
        try:
            export_spectrocube(
                self.spectra,
                path,
                units=units,
                instrument_id="echelle",
                shot_number=str(shot),
            )
            self.statusBar().showMessage(f"SpectroCube saved: {path}")
        except ImportError as exc:
            self.statusBar().showMessage(
                "SpectroCube export requires the 'spectrocube' package — see status bar."
            )
            QtWidgets.QMessageBox.warning(
                self,
                "SpectroCube not installed",
                "SpectroCube export requires the optional 'spectrocube' package.\n\n"
                "Install it with:\n    pip install spectrocube\n\n"
                "Or for local development:\n"
                "    pip install -e /path/to/2026-spectrocube\n\n"
                f"Details: {exc}",
            )
        except Exception as exc:
            self.statusBar().showMessage(f"SpectroCube export failed: {exc}")
            QtWidgets.QMessageBox.critical(
                self,
                "Export failed",
                f"SpectroCube export failed:\n{exc}",
            )

    # ===========================================================================
    #                   Display Plots and info
    # ===========================================================================

    def show_fit_intensities(self):
        """ """
        [self.fres[b.name].fill_nans() for b in self.bandstofit]

        for b in self.bandstofit:
            fr = self.fres[b.name].frames
            ints_comb = self.fres[b.name].intensities_combined[fr]
            m = np.invert(np.isnan(ints_comb))
            self.c_he_trace[b.name].setData(x=fr[m], y=ints_comb[m])

    def show_info(self):
        """Show SIF image info"""
        # show image info in browser
        txt = """
              <font size = {fs}>frames: {frames}</font><br>
              <font size = {fs}>background: {background}</font><br>
              <font size = {fs}>T = {temp}</font><br>
              <font size = {fs}>exposure: {exposure}</font><br>
              <font size = {fs}>cycle: {cycle}</font><br>
              <font size = {fs}>xbin: {xbin}</font><br>
              <font size = {fs}>ybin: {ybin}</font><br>
              """
        background_str = (
            "[{}..{}] ({})".format(
                self.em.info["BackgroundFrames"][0],
                self.em.info["BackgroundFrames"][-1],
                len(self.em.info["BackgroundFrames"]),
            )
            if len(self.em.info["BackgroundFrames"]) != 0
            else "None"
        )
        txt = txt.format(
            frames=self.em.info["NumberOfFrames"],
            background=background_str,
            temp=self.em.info["DetectorTemperature"],
            exposure=self.em.info["ExposureTime"],
            cycle=self.em.info["CycleTime"],
            xbin=self.em.info["xbin"],
            ybin=self.em.info["ybin"],
            fs=4,  # font size
        )

        self.image_info_bw.setText(txt + self._calibration_info_html())

    def _calibration_info_html(self):
        """Name the overriding snapshot beside the frame it was read with"""
        if self.calibration_override is None:
            return ""
        return '<font size = 4 color = "#187031">calibration: {}</font><br>'.format(
            self.calibration_override.snapshot_id
        )

    def show_image_frame(self):
        """Display current frame of the loaded SIF file"""
        frame = int(self.frame.value())
        self.frame_current = frame

        self.show_detector_frame(frame)
        self.show_spectrum()

    def show_detector_frame(self, frame):
        """Put one frame on the image, levelling it only when it is new.

        The rule is the echelle-aware one in
        :mod:`~echelle_spectra.tools.detector_display`: black at the background
        floor, white at the 99.9th percentile, so the strong lamp lines clip
        and the faint ones stay visible.  It runs when the displayed frame
        changes and never again — ``setImage`` is told not to level for itself
        either — so a histogram the operator has dragged survives every redraw
        of the same frame.  The levels are set here; they are not held here.
        """
        key = (str(getattr(self, "filename", "")), int(frame))
        fresh = key != self._levelled_frame
        levels = None
        if fresh and self.check_autoscale.isChecked():
            levels = auto_display_levels(self.em.images[frame])
        # Only a frame nobody has levelled yet falls back to pyqtgraph's own
        # min/max, and only so that switching the rule off still shows an image.
        self.img.setImage(self.em.images[frame], autoLevels=fresh and levels is None)
        if levels is not None:
            self.hist.setLevels(min=levels[0], max=levels[1])
        self._levelled_frame = key

    def show_spectrum(self):
        """Show the current frame spectrum in counts and calibrated"""
        frame = int(self.frame.value())
        counts = self.spectra.counts[frame]
        self.spec_counts.setData(x=self.spectra.wavelength, y=counts)

        # show calibrated spectrum
        # str() - for PyQt4, to convert QString into py string
        units = str(self.spec_units.currentText())
        self.p3.setLabel(
            "left",
            text=f'<font color = "#d62300" size="6">{self.spec_labels[units]}</font>',
        )
        calibrated = self.spectra.spectra_to_save[units][frame]
        self.spec_wm.setData(x=self.spectra.wavelength, y=calibrated)
        # The overlay ranks its labels on the trace that is actually on screen,
        # so a NIST-weak line sitting on a bright blob is named before a
        # NIST-strong line sitting on nothing.
        self.line_overlays.set_measured_spectrum(
            "counts", self.spectra.wavelength, counts
        )
        self.line_overlays.set_measured_spectrum(
            "calibrated", self.spectra.wavelength, calibrated
        )
        self.line_overlays.refresh()

    def show_balmer_frame(self):
        """Update the plots for H-balmer series experimental spectra"""
        frame = int(self.frame_h.value())
        self.frame_current = frame

        if frame >= self.em.info["NumberOfFrames"]:
            return

        for b, k in zip(
            self.hbands + [self.chband], [b.name for b in self.hbands] + ["CH"]
        ):
            if self.fres[b.name].x is None:
                x, y = eb.banddata(b, frame, self.spectra, doblin=False)
                self.h_data[k].setData(x=x, y=y)
            else:
                x = self.fres[b.name].x
                y = self.fres[b.name].spectra[frame - 1]
                self.h_data[k].setData(x=x, y=y)

    def show_c_frame(self):
        """Update the plots for C experimental spectra"""
        frame = int(self.frame_civ.value())
        self.frame_current = frame

        if frame >= self.em.info["NumberOfFrames"]:
            return

        for b in self.cbands:
            if self.fres[b.name].x is None:
                x, y = eb.banddata(b, frame, self.spectra)
                self.civ_data[b.name].setData(x=x, y=y)
            else:
                x = self.fres[b.name].x
                y = self.fres[b.name].spectra[frame - 1]
                self.civ_data[b.name].setData(x=x, y=y)
                if self.fres[b.name].spectra_fit[frame] is not None:
                    xx = self.fres[b.name].x_fit
                    ys = self.fres[b.name].spectra_fit[frame]
                    self.civ_fit[b.name].setData(x=xx, y=ys.sum(axis=0))
                else:
                    self.civ_fit[b.name].setData(x=[], y=[])

    def show_he_frame(self):
        """Update the plots for He experimental spectra"""
        frame = int(self.frame_he.value())
        self.frame_current = frame

        if frame >= self.em.info["NumberOfFrames"]:
            return

        for b in self.hebands:
            if self.fres[b.name].x is None:
                x, y = eb.banddata(b, frame, self.spectra)
                self.he_data[b.name].setData(x=x, y=y)
            else:
                x = self.fres[b.name].x
                y = self.fres[b.name].spectra[frame - 1]
                self.he_data[b.name].setData(x=x, y=y)
                if self.fres[b.name].spectra_fit[frame] is not None:
                    xx = self.fres[b.name].x_fit
                    ys = self.fres[b.name].spectra_fit[frame]
                    self.he_fit[b.name].setData(x=xx, y=ys.sum(axis=0))
                else:
                    self.he_fit[b.name].setData(x=[], y=[])

    def clear_fit_traces(self):
        """Clear traces of fitted curves"""
        for trace in self.c_he_trace.values():
            trace.setData(x=[], y=[])

    # ===========================================================================
    #             ^^^   Display Plots and info   ^^^
    # ===========================================================================

    def _enable_controls(self, enable=True):
        """Enable or disable controls"""
        for ctrl in self.controls_open + self.controls_reg:
            ctrl.setEnabled(enable)
        self.setAcceptDrops(enable)

    def _reset_frame(self):
        """Reset frames on all tabs"""
        frame = self.frame_current
        for gi in self.frame_spinners:
            gi.setMaximum(1)
        window_layout.gui_setup_spinbox(self.frame, 0, 0, 0)
        self.frame_current = frame

    def _setup_frame(self):
        """Setup frames on all tabs"""
        window_layout.gui_setup_spinbox(
            self.frame, 0, 0, self.em.info["NumberOfFrames"] - 1
        )
        for gi in self.frame_spinners:
            gi.setMaximum(self.em.info["NumberOfFrames"] - 1)

        if self.frame_current > int(self.em.info["NumberOfFrames"]) - 1:
            self.frame.setValue(0)
        else:
            for f in [self.frame] + self.frame_spinners:
                f.setValue(self.frame_current)

    def emitted(self):
        self.load_image()

    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.accept()
        else:
            ev.ignore()

    def dropEvent(self, ev):
        if ev.mimeData().hasUrls:
            ev.setDropAction(QtCore.Qt.CopyAction)
            ev.accept()
            filename = ev.mimeData().urls()[0]
            self.filename = filename.toLocalFile()
            self.emitted()
        else:
            ev.ignore()


class LoadCalibrationsThread(QtCore.QThread):
    """Thread for loading calibration files"""

    taskFinished = QtCore.pyqtSignal(object)

    def __init__(self, calibration, config):
        QtCore.QThread.__init__(self)
        self.calibration = calibration
        self.config = config

    def run(self):
        cb = self.calibration
        try:
            cb.start()
        except Exception as err:
            # Reporting the failure keeps the window usable; a dead thread would
            # leave the controls disabled and the startup message forever on.
            self.taskFinished.emit(
                CalibrationLoadOutcome(cb.name, cb, f"{type(err).__name__}: {err}")
            )
            return

        if self.config["debug"]:
            print("\n" + cb.name, cb.DIMO, cb.DIMW)
            [print(a, " " * (12 - len(a)), b) for a, b in cb.filenames.items()]

        self.taskFinished.emit(CalibrationLoadOutcome(cb.name, cb))


class LoadImageThread(QtCore.QThread):
    """Thread for loading SIF image"""

    taskFinished = QtCore.pyqtSignal(object)

    def __init__(self, filename, cb, config):
        QtCore.QThread.__init__(self)
        self.filename = filename
        self.cb = cb
        self.config = config

    def run(self):
        expected = calibration_dimensions(self.cb)
        if not len(expected):
            self.taskFinished.emit(
                ImageLoadOutcome(
                    status=IMAGE_CALIBRATION_UNAVAILABLE,
                    camera=self.cb.name,
                    detail="its calibration files did not load",
                )
            )
            return

        try:
            em = ech.EchelleImage(self.filename, clbr=self.cb)
        except Exception as err:
            if self.config["debug"]:
                print(f"LoadImageThread Error: {err}")
            self.taskFinished.emit(
                ImageLoadOutcome(
                    status=IMAGE_UNREADABLE,
                    camera=self.cb.name,
                    expected_dimensions=expected,
                    detail=f"{type(err).__name__}: {err}",
                )
            )
            return

        dimensions = detector_dimensions(em.info)
        binning = (int(em.info["xbin"]), int(em.info["ybin"]))

        if self.config["debug"]:
            print(dimensions, self.cb.DIMO, self.cb.DIMW, self.cb.name)

        if dimensions != expected:
            self.taskFinished.emit(
                ImageLoadOutcome(
                    status=IMAGE_DIMENSION_MISMATCH,
                    camera=self.cb.name,
                    file_dimensions=dimensions,
                    binning=binning,
                    expected_dimensions=expected,
                )
            )
            return

        try:
            em.calculate_order_spectra()  # image -> order spectra
            em.correct_order_shapes()  # remove out of bounds boundaries
            em.calculate_spectra()  # order spectra -> fullwidth spectra
        except Exception as err:
            # Extraction used to run outside any guard, so a failure here killed
            # the thread without ever emitting and the window waited forever.
            self.taskFinished.emit(
                ImageLoadOutcome(
                    status=IMAGE_EXTRACTION_FAILED,
                    camera=self.cb.name,
                    file_dimensions=dimensions,
                    binning=binning,
                    expected_dimensions=expected,
                    detail=f"{type(err).__name__}: {err}",
                )
            )
            return

        self.taskFinished.emit(
            ImageLoadOutcome(
                status=IMAGE_LOADED,
                camera=self.cb.name,
                image=em,
                file_dimensions=dimensions,
                binning=binning,
                expected_dimensions=expected,
            )
        )


class FitLinesThread(QtCore.QThread):
    """Thread for fitting lines (CIV, He, H)"""

    pass_result = QtCore.pyqtSignal(object)
    progress = QtCore.pyqtSignal(int)

    def __init__(self, spectra, bands, fres, frames):
        QtCore.QThread.__init__(self)
        self.spectra = spectra
        self.bands = bands
        self.fres = fres
        self.frames = frames

    def run(self):
        for band in self.bands:
            self.fres[band.name].set_frames(self.frames)
            self.fres[band.name].init_spectra()

        for frame in self.frames:
            for b in self.bands:
                res = b.fitb(frame, self.spectra, detect=True, doblin=False)
                self.fres[b.name].spectra.append(res["exp"][1])
                self.fres[b.name].x = res["exp"][0]
                try:
                    self.fres[b.name].spectra_fit[frame] = res["fit"][1]
                    self.fres[b.name].x_fit = res["fit"][0]
                    self.fres[b.name].fit_out[frame] = res["out"]
                except TypeError:
                    pass

            self.progress.emit(frame + 1)

        for band in self.bands:
            self.fres[band.name].integrate_spectra()
            self.fres[band.name].integrate_fit()

        self.pass_result.emit(self.fres)


class SaveSpectraThread(QtCore.QThread):
    """Thread for saving spectra"""

    pass_result = QtCore.pyqtSignal(object)

    def __init__(self, arg):
        QtCore.QThread.__init__(self)
        self.spectra = arg

    def run(self):
        self.spectra.save()
        self.pass_result.emit(self.spectra)


def parse_args(argv=None):
    """Read the viewer's own switches and hand the rest to Qt.

    Qt spells its own switches with a single dash (``-platform offscreen``), so
    everything left over that starts with ``--`` is a misspelling of ours rather
    than a message for Qt, and is refused instead of quietly ignored.  With no
    arguments at all nothing is parsed away and Qt receives what it always did.
    """
    parser = argparse.ArgumentParser(
        prog="echelle_spectra",
        description="Echelle viewer: open one SIF at a time.",
    )
    parser.add_argument(
        "--calibration",
        metavar="SNAPSHOT",
        default=None,
        # Plain ASCII, deliberately: --help has to print on a Japanese console
        # too, where an em dash is not encodable.
        help=(
            "Open files through a saved calibration snapshot folder (the one "
            "holding snapshot.toml) instead of the packaged CCD/CMOS tables."
        ),
    )
    args, extra = parser.parse_known_args(argv)
    unknown = [item for item in extra if item.startswith("--")]
    if unknown:
        parser.error("unrecognized arguments: " + " ".join(unknown))
    return args, extra


def start(argv=None):
    import platform

    args, qt_args = parse_args(argv)
    override = None
    if args.calibration is not None:
        try:
            override = load_calibration_override(args.calibration)
        except CalibrationOverrideError as err:
            # Said once, on the way out.  A window wearing the packaged tables
            # while the operator believes he asked for 2019 is the one outcome
            # this flag exists to prevent.
            print(f"echelle_spectra: {err}", file=sys.stderr)
            return 2

    if platform.system() == "Windows":
        try:
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                f"echelle_spectra-{__version__}"
            )
        except OSError:
            pass

    app = QApplication(sys.argv[:1] + list(qt_args))
    app.setWindowIcon(
        QtGui.QIcon(str(_config["base_path"] / "resources/graphics/echelle.png"))
    )
    win = EchelleSpectraGUI(_config, calibration_override=override)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    start()
