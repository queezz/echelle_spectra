import ctypes
import numpy as np
import sys
from datetime import datetime
from pathlib import Path
from PyQt5.QtWidgets import QApplication
from pyqtgraph.Qt import QtCore, QtGui

import tools.echelle as ech
import tools.emissionbands as eb
import tools.emissiondata as ebd
from __init__ import __version__
from __init__ import _config
from resources import window_layout


class EchelleSpectraGUI(QtGui.QMainWindow, window_layout.Ui_MainWindow):
    """GUI window for the echelle_spectra app"""
    def __init__(self, config):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        # set widget statuses
        self.CameraCCD.setChecked(False)
        self.CameraCMOS.setChecked(True)
        self.spec_counts = self.p2.plot(pen="#6ac600")
        self.spec_wm = self.p3.plot(pen="r")

        # define initial class attributes
        self.config = config
        self.frame_spinners = [self.frame_civ, self.frame_he, self.frame_h]
        self.cameras = [self.CameraCCD, self.CameraCMOS]

        # carry out init actions
        self.connect_actions()
        self.update_paths()
        self.read_last_shot()
        self.prepare_calibration()
        self.setup_bands()

        # define auxiliary class attributes
        self.header = self.get_header_template()
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
            "orders": "pattern_cmos.txt",
            "wavelength": "Th_wavelength_CMOS.txt",
            "sphr": "sphere_CMOS.sif",
            "bkgr": "sphere_CMOS_bkg.sif",
            "integral": "integrating_sphere.txt",
        }

        self.cb_CCD = ech.Calibrations(self.path_calibration)
        self.cb_CCD.name = "CCD"
        self.cb_CMOS = ech.Calibrations(self.path_calibration)
        self.cb_CMOS.name = "CMOS"

        self.cb_CCD.filenames = files_ccd
        self.cb_CMOS.filenames = files_cmos

        self.calib_thread_ccd = None
        self.calib_files_cmos = None

        self.calib_threads = {
            "CCD": self.calib_thread_ccd,
            "CMOS": self.calib_files_cmos,
        }

        self.load_calibration(self.cb_CCD)
        self.load_calibration(self.cb_CMOS)

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

    def update_paths(self):
        """Set-up file paths for calibration, data folder, output folder etc."""
        for path in ["data_path", "output_path"]:
            setattr(self, path, Path(self.config[path]
                                     .replace("{homedir}", str(Path.home()))
                                     .replace("{workdir}", str(Path().absolute()))))
            getattr(self, path).mkdir(parents=True, exist_ok=True)

        self.path_calibration = self.config["base_path"] / "resources/calibration_files"
        self.path_header_template = self.config["base_path"] / "resources/header_template.txt"
        self.path_last_shot = self.config["base_path"] / ".last_shot"

    def get_header_template(self):
        """Get header template used when saving data

        Usage:

        header = self.header.format(
            diag_name = 'alfa',
            shot = 10,
            date = 2019,
            dimno = 1,
            dimname = 'time',
            dimsize = 1,
            dimunits = 's',
            nval = 20,
            vnames = 'a b c',
            vunit = 'W/m',
            trigdelay = 2.75,
            cycletime = 0.25
            )
        """
        with open(self.path_header_template, "r") as fl:
            return fl.read()

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
        self.coursor_bw.setText("""<font size = 6 color = "#d1451b">Loading Calibrations</font>""")

        self.calib_threads[clbr.name] = LoadCalibrationsThread(clbr, self.config)
        self.calib_threads[clbr.name].taskFinished.connect(self._on_calibration_loaded)
        self.calib_threads[clbr.name].start()

    def _on_calibration_loaded(self, result):
        """When calibration loaded, release buttons and get ready to work"""
        if result.name == "CCD":
            self.cb_CCD = result
            if self.config['debug']:
                print("CCD")

        if result.name == "CMOS":
            self.cb_CMOS = result
            if self.config['debug']:
                print("CMOS")

        self.statusBar().showMessage("Calibration files loaded. Ready to work.")
        self._enable_controls(True)
        self.coursor_bw.setText("")

    def setup_bands(self):
        """Create list of band objects from EmissionBand class"""
        self.hbands = [band.copy() for band in [ebd.halpha, ebd.hbeta, ebd.hgamma, ebd.hdelta]]
        self.hebands = [band.copy() for band in [ebd.he447, ebd.he492, ebd.he501, ebd.he504, ebd.he587, ebd.he667,
                                                 ebd.he706, ebd.he728]]
        self.cbands = [band.copy() for band in [ebd.c515, ebd.c464, ebd.c580, ebd.c444, ebd.c465, ebd.c772, ebd.c706,
                                                ebd.c547]]
        self.chband = eb.EB(name="CHD-band", bounds=[429, 431.5])

        self.bands = self.hbands + self.hebands + self.cbands + [self.chband]
        self.bandstofit = self.cbands + self.hebands

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

        self.shot_range = [shot for shot in self.shot_range if shot in list(self.shot_dict)]
        if not len(self.shot_range):
            self.coursor_bw.setText('<font size = 6 color = "#e5c40d">Empty Range</font>')
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
            self.coursor_bw.append('<font size = 6 color = "#d1451b">Aborting Loop</font>')
            self.abort_register.setEnabled(False)

    # ===========================================================================
    #                   ^^^  Loop through shots ^^^
    # ===========================================================================

    # ===========================================================================
    #                       Open shot image or file
    # ===========================================================================

    def openfile(self):
        """Open file dialog for selection of SIF file to load"""
        self.filename = QtGui.QFileDialog.getOpenFileName(None, "Open file to plot", str(self.data_path),
                                                          "*.sif;;*.SIF;;*.*")[0]
        if self.filename:
            self.emitted()

    def update_shot_dict(self):
        """Convert list of *_Echelle.sif files into shot list and make a dictionary of {shot_number: path to shot}"""
        self.shot_dict = {int(f.stem.split("_")[0]): f for f in self.data_path.rglob("*_Echelle.SIF")}

    def get_shot_range(self):
        """Get first and last shot numbers from the data folder after indexing it"""
        self.update_shot_dict()
        self.start_shot.setValue(min(self.shot_dict))
        self.end_shot.setValue(max(self.shot_dict))

    def load_shot_number(self):
        """Try to load image of the selected shot number from the data folder"""
        self.update_shot_dict()
        try:
            self.filename = self.shot_dict[self.shot_number.value()]
            self.load_image()
        except KeyError:
            self.statusBar().showMessage(f"Failed to load shot {self.shot_number.value()}")

    # ===========================================================================
    #                 ^^^   Open shot image or file   ^^^
    # ===========================================================================

    def load_image(self):
        """Load SIF image"""
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

        if self.CameraCMOS.isChecked():
            cb = self.cb_CMOS
            if self.config['debug']:
                print("CMOS selected")
        elif self.CameraCCD.isChecked():
            cb = self.cb_CCD
            if self.config['debug']:
                print("CCD selected")
        else:
            cb = self.cb_CCD
            self.CameraCCD.setChecked(True)

        self.image_load_thread = LoadImageThread(self.filename, cb, self.config)
        self.image_load_thread.taskFinished.connect(self._on_image_loaded)
        self.image_load_thread.start()

    def _on_image_loaded(self, result):
        """Receive the loaded Echelle Image"""
        self.em = result

        if result is None:
            self._enable_controls(True)
            txt = "DIMENSIONS of the Image and Calibration do not match"
            txt = '<font size = 6 color = "#d1451b">{}</font>'.format(txt)
            self.image_info_bw.setText(txt)

            self.change_camera()
            self.load_image()
            return

        self.spectra = ech.Spectrum(self.em)

        self._reset_frame()
        self._setup_frame()

        self.show_info()
        self.show_image_frame()
        self.show_c_frame()
        self.show_he_frame()
        self.show_balmer_frame()

        if not self.in_loop:
            self._enable_controls(True)
        self.statusBar().showMessage(f"SIF image loaded: {self.filename}")
        self.coursor_bw.setText(f'<font size = 5 color = "#187031">{Path(self.filename).stem}</font>')

        if self.specsave_bx.isChecked():
            self.spectra.shotnumber = int(self.shot_number.value())
            self.spectra.output_path = self.output_path
            self.spectra.trigdelay = self.trigger_delay.value()

            self.save_spec_thread = SaveSpectraThread(self.spectra)
            self.save_spec_thread.pass_result.connect(self._on_spec_saved)
            self.save_spec_thread.start()
        else:
            self._on_spec_saved(None)

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
            self.htraces[band.name].setData(x=full_frames, y=self.fres[band.name].intensities)
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

        self.fit_thread = FitLinesThread(self.spectra, self.bandstofit, self.fres, frames)
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
        data = np.array([times] + [self.fres[b.name].intensities_combined for b in self.bands])
        data = data.T[1:]  # [1:] to remove frame 0, where values are np.nan

        shot = self.shot_number.value()
        pth = self.output_path / f"{shot}_{self.config['diag_name']}.txt"

        if Path(pth).is_file() and not self.overwrite.isChecked():
            return

        header = self.header.format(
            diag_name=self.config["diag_name"],
            shot=shot,
            date=datetime.now().strftime("%m/%d/%Y %H:%M"),
            dimno=1,
            dimname="Time",
            dimsize=data.shape[0],
            dimunits="s",
            nval=data.shape[1] - 1,
            vnames=",".join([f"'{name}'" for name in names]),
            vunit=",".join(["'W/m^2'"] * len(names)),
            trigdelay=self.trigger_delay.value(),
            cycletime=self.em.info["CycleTime"],
        ).strip()
        np.savetxt(pth, data, delimiter=", ", header=header, comments="", fmt=["%.5f"] + ["%.6e"] * len(names))

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
        background_str = ("[{}..{}] ({})".format(self.em.info["BackgroundFrames"][0],
                                                 self.em.info["BackgroundFrames"][-1],
                                                 len(self.em.info["BackgroundFrames"]))
                          if len(self.em.info["BackgroundFrames"]) != 0 else "None")
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

        self.image_info_bw.setText(txt)

    def show_image_frame(self):
        """Display current frame of the loaded SIF file"""
        frame = int(self.frame.value())
        self.frame_current = frame

        self.img.setImage(self.em.images[frame])
        if self.check_autoscale.isChecked():
            self.hist.setLevels(min=np.percentile(self.em.images[frame], 1),
                                max=np.percentile(self.em.images[frame], 99))
        self.show_spectrum()

    def show_spectrum(self):
        """Show the current frame spectrum in counts and calibrated"""
        frame = int(self.frame.value())
        self.spec_counts.setData(x=self.spectra.wavelength, y=self.spectra.counts[frame])

        # show calibrated spectrum
        # str() - for PyQt4, to convert QString into py string
        units = str(self.spec_units.currentText())
        self.p3.setLabel("left", text=f'<font color = "#d62300" size="6">{self.spec_labels[units]}</font>')
        self.spec_wm.setData(x=self.spectra.wavelength, y=self.spectra.spectra_to_save[units][frame])

    def show_balmer_frame(self):
        """Update the plots for H-balmer series experimental spectra"""
        frame = int(self.frame_h.value())
        self.frame_current = frame

        if frame >= self.em.info["NumberOfFrames"]:
            return

        for b, k in zip(self.hbands + [self.chband], [b.name for b in self.hbands] + ["CH"]):
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
        window_layout.gui_setup_spinbox(self.frame, 0, 0, self.em.info["NumberOfFrames"] - 1)
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
        self.calibration.start()
        cb = self.calibration

        if self.config['debug']:
            print("\n" + cb.name, cb.DIMO, cb.DIMW)
            [print(a, " " * (12 - len(a)), b) for a, b in cb.filenames.items()]

        self.taskFinished.emit(self.calibration)


class LoadImageThread(QtCore.QThread):
    """Thread for loading SIF image"""

    taskFinished = QtCore.pyqtSignal(object)

    def __init__(self, filename, cb, config):
        QtCore.QThread.__init__(self)
        self.filename = filename
        self.cb = cb
        self.config = config

    def run(self):
        try:
            em = ech.EchelleImage(self.filename, clbr=self.cb)
        except Exception as err:
            if self.config["debug"]:
                print(f"LoadImageThread Error: {err}")
            self.taskFinished.emit(None)
            return

        if self.config['debug']:
            print(em.info["DetectorDimensions"], self.cb.DIMO, self.cb.DIMW, self.cb.name)

        if not em.info["DetectorDimensions"] == (self.cb.DIMW, self.cb.DIMO):
            self.taskFinished.emit(None)
            return

        em.calculate_order_spectra()  # image -> order spectra
        em.correct_order_shapes()  # remove out of bounds boundaries
        em.calculate_spectra()  # order spectra -> fullwidth spectra

        self.taskFinished.emit(em)


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


def start():
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(f"echelle_spectra-{__version__}")
    except OSError:
        pass

    app = QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(str(_config["base_path"] / "resources/graphics/echelle.png")))
    win = EchelleSpectraGUI(_config)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    start()
