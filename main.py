from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg

from gui import mainwindow
import tools.echelle as ech
import tools.emissionbands as eb
import tools.emissiondata as ebd
import numpy as np
import os
from os.path import join, basename

import matplotlib.pylab as plt

DEBUG = False

NOF = "NumberOfFrames"  # Just to shorten this key for the number of frames


class EchelleAnalizer(QtGui.QMainWindow, mainwindow.Ui_MainWindow):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.frame_spinners = [self.frame_civ, self.frame_he, self.frame_h]

        self.connect_actions()
        self.update_path()
        self.read_header_template()
        self.shotdict = None
        self.shot_range = []
        self.read_startup_settings()

        self.prepare_calibration()

        self.spec_counts = self.p2.plot(pen="#6ac600")
        self.spec_wm = self.p3.plot(pen="r")

        self.setup_bands()

        self.frame_current = 0
        self.inloop = False

        self.CameraCCD.setChecked(False)
        self.CameraCMOS.setChecked(True)
        self.cameras = [self.CameraCCD, self.CameraCMOS]

    def change_camera(self):
        """ Change the active Camera
        """
        for cam in self.cameras:
            if not cam.isChecked():
                cam.setChecked(True)
                break

    def prepare_calibration(self):
        """ Select calibration files here
        """
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
        self.frame.valueChanged.connect(self.show_imageframe)
        self.frame_h.valueChanged.connect(self.show_balmerframe)
        self.frame_he.valueChanged.connect(self.show_heframe)
        self.frame_civ.valueChanged.connect(self.show_cframe)
        self.btn_open.clicked.connect(self.openfile)
        self.shot_range_btn.clicked.connect(self.get_shot_range)
        self.start_pth = None

        self.show_btn.clicked.connect(self.load_shotnumber)

        self.start_register.clicked.connect(self.start_loop)
        self.abort_register.clicked.connect(self.abort_loop)

    def update_path(self):
        """ Set-up file paths for calibraiton, data folder, output folder etc.
        """
        self.start_path = None
        self.path_calibration = "calibration_files"
        self.path_startup_options = "startup_options.txt"
        try:
            with open("path_datafolder.txt", "r") as f:
                self.path_data = f.readline().strip()
        except:
            self.path_data = ".."

        try:
            with open("path_output.txt", "r") as f:
                self.path_output = f.readline().strip()
        except:
            self.path_output = "../Echelle_register"

        if not os.path.exists(self.path_output):
            os.makedirs(self.path_output)

        print(self.path_output)

        self.path_header_template = "header_template.txt"

    def read_header_template(self):
        """ read hader template to use when saving data
        
        usage:
        
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
            self.header = fl.read()

    def write_startup_settings(self):
        """ Write startup settings to file,
        such as a shotnumber of the last image loaded
        """
        try:
            with open(self.path_startup_options, "w") as f:
                f.write("{} # shot number".format(self.shot_number.value()))
        except:
            pass

    def read_startup_settings(self):
        """ """
        try:
            with open(self.path_startup_options, "r") as f:
                for line in f.readlines():
                    if "shot number" in line:
                        shot = int(line.split("#")[0])
            self.shot_number.setValue(shot)
        except:
            pass

    def load_calibration(self, clbr):
        """ Load calibration files and prepare Calibrations class
        calibration - calibration instance for CCD or CMOS
        """
        msg = "Loading calibration files"
        self.statusBar().showMessage(msg)
        self._enablecontrols(False)
        txt = """<font size = 6 color = "#d1451b">Loading Calibrations</font>"""
        self.coursor_bw.setText(txt)

        self.calib_threads[clbr.name] = LoadCalibrations(clbr)
        self.calib_threads[clbr.name].taskFinished.connect(self.onCalibrationLoaded)
        self.calib_threads[clbr.name].start()

        # self.calib_thread = LoadCalibrations(clbr)
        # self.calib_thread.taskFinished.connect(self.onCalibrationLoaded)
        # self.calib_thread.start()

    def onCalibrationLoaded(self, result):
        """ When calibration loaded, release buttons and get ready to work
        """
        if result.name == "CCD":
            self.cb_CCD = result
            if DEBUG:
                print("CCD")

        if result.name == "CMOS":
            self.cb_CMOS = result
            if DEBUG:
                print("CMOS")

        msg = "Calibration files loaded. Ready to work."
        self.statusBar().showMessage(msg)
        self._enablecontrols(True)
        self.coursor_bw.setText("")

        # self.load_calibration(self.cb_CMOS)

    def setup_bands(self):
        """ create list of bands objects from EmissionBand class
        """
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
        """ Start cycle trhough shots in range(self.start_shot,self.end_shot)
        """
        self.shot_range = range(self.start_shot.value(), self.end_shot.value() + 1)
        self.shot_range_index = 0
        if self.shotdict is None:
            self.make_shotdict()

        shots_in_df = list(self.shotdict)
        self.shot_range = [i for i in self.shot_range if i in shots_in_df]
        if not len(self.shot_range):
            txt = '<font size = 6 color = "#e5c40d">Empyt Range</font>'
            self.coursor_bw.setText(txt)
            return

        self.progress_range.setRange(0, len(self.shot_range))

        self.shot_range_len = len(self.shot_range)
        self.abort_register.setEnabled(True)
        self.inloop = True
        self.loop_step()

    def loop_step(self):
        """ update curent shotnuber and run loadimage
        """
        self.shot_number.setValue(self.shot_range[0])
        self.shot_range = self.shot_range[1:]
        self.progress_range.setValue(self.shot_range_index)
        self.shot_range_index += 1
        self.load_shotnumber()

    def loop_advance(self):
        """ if in the loop go to next shot
        """
        if len(self.shot_range):
            self.loop_step()
        else:
            # End of the loop is here
            self.progress_range.setValue(self.shot_range_index)
            txt = '<font size = 6 color = "#d1451b">Loop\'s Done<br>{}/{}</font>'
            txt = txt.format(self.shot_range_index, self.shot_range_len)
            self.coursor_bw.setText(txt)
            self._enablecontrols(True)
            self.inloop = False

    def advance_if_inloop(self):
        """ if in the loop go to next shot """
        if self.inloop:
            self.loop_advance()

    def abort_loop(self):
        """ Abort the loop
        """
        if len(self.shot_range):
            self.shot_range = []
            txt = '<font size = 6 color = "#d1451b">Aborting Loop</font>'
            self.coursor_bw.append(txt)
            self.abort_register.setEnabled(False)

    # ===========================================================================
    #                   ^^^  Loop through shots ^^^
    # ===========================================================================

    # ===========================================================================
    #                       Open shot image or file
    # ===========================================================================

    def openfile(self):
        """ Get data path from settings.txt
        For now - directly specify here
        """
        if self.start_path is None:
            self.start_pth = self.path_data
        self.filename = QtGui.QFileDialog.getOpenFileName(
            None, "Open file to plot", self.start_pth, "*.sif;;*.SIF;;*.*"
        )[0]
        self.start_pth = self.filename
        if self.filename:
            self.emitted()

    def index_datafolder(self, **kws):
        """ get list of *_Echelle.sif files in the
        Echelle data subfolders and compile a shotlist
        Python3 is preferable, Python2 support is included
        """
        suffix = kws.get("suffix", "_Echelle.SIF")  # _Echelle.sif
        pth = self.path_data
        matches = []
        if sys.version_info[0] < 3:
            import fnmatch

            for root, dirnames, filenames in os.walk(pth):
                for filename in fnmatch.filter(filenames, "*{}".format(suffix)):
                    matches.append(join(root, filename))
        else:
            import glob

            pth += "/**/" + "*{}".format(suffix)
            matches = (filename for filename in glob.iglob(pth, recursive=True))
        return matches

    def make_shotdict(self, **kws):
        """ Convert list of *_Echelle.sif files into shot list and
        make a dictionary of {shot_number : path to shot}
        """
        suffix = kws.get("suffix", "_Ech")
        datalist = self.index_datafolder()
        self.shotdict = {int(os.path.basename(i).split(suffix)[0]): i for i in datalist}

    def get_shot_range(self):
        """ get first and last shotnumbers from the datafolder after indexing it
        """
        self.make_shotdict()
        shots = sorted(list(self.shotdict))
        a, b = shots[0], shots[-1]
        self.start_shot.setValue(a)
        self.end_shot.setValue(b)

    def load_shotnumber(self):
        """ Try to load image of the selected shotnuber from the datafolder
        """
        if self.shotdict is None:
            self.make_shotdict()
        try:
            self.filename = self.shotdict[self.shot_number.value()]
            self.loadimage()
        except KeyError:
            msg = "Failed to loaded shot {}, reindexing"
            msg = msg.format(self.shot_number.value())
            self.statusBar().showMessage(msg)
            self.update_path()
            self.make_shotdict()
            try:
                self.filename = self.shotdict[self.shot_number.value()]
                self.loadimage()
            except:
                msg = "Failed to load shot {}".format(self.shot_number.value())
                self.statusBar().showMessage(msg)

    # ===========================================================================
    #                 ^^^   Open shot image or file   ^^^
    # ===========================================================================

    def loadimage(self):
        """ Load Sif Image """
        if not self.inloop:
            self.progress_range.setRange(0, 1)
            self.clear_fit_traces()

        self.progress_bands.setRange(0, 1)

        self._enablecontrols(False)

        self.fres = {band.name: eb.FitResult(band.name) for band in self.bands}
        msg = "Loading SIF image: {}".format(self.filename)
        self.statusBar().showMessage(msg)
        # setup shotnumber from the image filename
        bname = basename(self.filename)
        try:
            shot = int(bname.split("_Echelle")[0])
        except:
            shot = 0
        self.shot_number.setValue(shot)
        self.write_startup_settings()
        self._enablecontrols(False)
        txt = '<font size = 6 color = "#d1451b">Loading Image</font>'
        self.coursor_bw.setText(txt)

        if self.CameraCMOS.isChecked():
            cb = self.cb_CMOS
            if DEBUG:
                print("CMOS selected")
        elif self.CameraCCD.isChecked():
            cb = self.cb_CCD
            if DEBUG:
                print("CCD selected")
        else:
            cb = self.cb_CCD
            self.CameraCCD.setChecked(True)

        self.image_thread = LoadImage([self.filename, cb])
        self.image_thread.taskFinished.connect(self.onImageLoaded)
        self.image_thread.start()

    def onImageLoaded(self, result):
        """ receive the loaded Echelle Image
        """
        self.em = result

        if result == None:
            self._enablecontrols(True)
            txt = "DIMENSIONS of the Image and Calibration do not match"
            txt = '<font size = 6 color = "#d1451b">{}</font>'.format(txt)
            self.image_info_bw.setText(txt)

            self.change_camera()

            self.loadimage()

            return

        self.spectra = ech.Spectrum(self.em)

        self._reset_frame()
        self._setup_frame()

        self.show_info()
        self.show_imageframe()
        if not self.inloop:
            self._enablecontrols(True)
        msg = "SIF image loaded: {}".format(self.filename)
        self.statusBar().showMessage(msg)

        txt = '<font size = 5 color = "#187031">{}</font>'
        txt = txt.format(basename(self.filename))
        self.coursor_bw.setText(txt)

        if self.specsave_bx.isChecked():
            self.spectra.shotnumber = int(self.shot_number.value())
            self.spectra.path_output = self.path_output
            self.spectra.trigdelay = self.trigger_delay.value()

            self.save_spec_thread = SaveSpectra(self.spectra)
            self.save_spec_thread.passresult.connect(self._onSpecSaved)
            self.save_spec_thread.start()
        else:
            self._onSpecSaved(None)

    def _onSpecSaved(self, result):
        """ Continue after spectra saved (or not) """

        if result is not None:
            self.spectra = result

        if self.spectra.info[NOF] > 1:
            self.no_fit_intesities()
        else:
            # if in loop saving spectra for 1 frame images
            self.advance_if_inloop()

    def no_fit_intesities(self):
        """ Put the spectra data into the FitResult object
        and integrate experimental spectra
        """
        full_frames = np.arange(self.spectra.info[NOF])
        frames = np.arange(self.spectra.info[NOF])
        # bands = self.hbands+[self.chband]
        # calculate intensities without fit for all bands for a start
        bands = self.bands
        [self.fres[b.name].set_frames(frames) for b in bands]
        [self.fres[b.name].init_spectra() for b in bands]
        [self.fres[b.name].set_info(self.em.info[NOF]) for b in bands]

        for band in bands:
            for frame in frames:
                x, y = eb.banddata(band, frame, self.spectra, doblin=False)
                self.fres[band.name].spectra.append(y)

            self.fres[band.name].x = x

        [self.fres[b.name].integrate_spectra() for b in bands]
        # fill intensities_combined with non-fited results
        [self.fres[b.name].integrate_fit() for b in bands]
        [self.fres[b.name].fill_nans() for b in bands]

        # Plot traces

        [
            self.htraces[b.name].setData(x=full_frames, y=self.fres[b.name].intensities)
            for b in self.hbands
        ]
        self.chtrace.setData(x=full_frames, y=self.fres[self.chband.name].intensities)

        # Also plot traces for bandstofit, without a fit for a start
        self.show_fit_intesities()

        if self.fit_lines_bx.isChecked():
            self.fit_lines()
        else:
            if self.save_lines_bx.isChecked():
                self.save_intensities()

            self.advance_if_inloop()  # if saving data in the loop, no fitting

    def fit_lines(self):
        """ Fit lines
        """
        msg = "Fitting Lines"
        self.statusBar().showMessage(msg)
        self._enablecontrols(False)
        self.progress_bands.setRange(1, self.spectra.info[NOF] - 1)

        # array of frames for wich fit will be tryed.
        frames = np.arange(1, self.em.info[NOF])

        [self.fres[b.name].set_frames(frames) for b in self.bandstofit]
        [self.fres[b.name].set_info(self.em.info[NOF]) for b in self.bandstofit]

        args = [self.spectra, self.bandstofit, self.fres, frames]
        self.fit_thread = FitLines(args)
        self.fit_thread.passresult.connect(self._onFitEnd)
        self.fit_thread.progress.connect(self._onFitProgress)
        self.fit_thread.start()

    def _onFitProgress(self, value):
        self.progress_bands.setValue(value)

    def _onFitEnd(self, result):
        """
        """
        self.fres = result
        msg = "Fit is Done"
        self.statusBar().showMessage(msg)
        if not self.inloop:
            self._enablecontrols(True)
        self.progress_bands.setRange(0, 1)

        self.show_fit_intesities()

        if self.save_lines_bx.isChecked():
            self.save_intensities()

        # Go to next shot if in loop.
        self.advance_if_inloop()

    def save_intensities(self):
        """ Save intesities """
        import datetime
        import time

        DIAGNAME = "spec_div1"

        # make numpy array for no-fits
        names = [b.name for b in self.bands]
        frames = np.arange(self.em.info[NOF])
        # self.trigger_delay.setValue(2.5) # to change the default value
        times = self.trigger_delay.value() + frames * self.em.info["CycleTime"]
        data = np.array(
            [times] + [self.fres[b.name].intensities_combined for b in self.bands]
        )
        data = data.T[1:]  # [1:] to remove frame 0, where values are np.nan
        self.path_output
        shot = self.shot_number.value()
        pth = join(self.path_output, "{}_{}.txt".format(shot, DIAGNAME))

        if os.path.exists(pth) and not self.overwrite.isChecked():
            return

        vnames = "'" + "','".join(name for name in names) + "'"
        vunit = "'" + "','".join("W/m^2" for name in names) + "'"
        fmt = ["%.5f"] + ["%.6e" for name in names]
        date = datetime.datetime.fromtimestamp(time.time())
        date = date.strftime("%m/%d/%Y %H:%M")
        header = self.header  # load template from separate file
        header = header.format(
            diag_name=DIAGNAME,
            shot=shot,
            date=date,
            dimno=1,
            dimname="Time",
            dimsize=data.shape[0],
            dimunits="s",
            nval=data.shape[1] - 1,
            vnames=vnames,
            vunit=vunit,
            trigdelay=self.trigger_delay.value(),
            cycletime=self.em.info["CycleTime"],
        )
        header = header.strip()
        np.savetxt(pth, data, delimiter=", ", header=header, comments="", fmt=fmt)

    # ===========================================================================
    #                   Display Plots and info
    # ===========================================================================

    def show_fit_intesities(self):
        """ """
        [self.fres[b.name].fill_nans() for b in self.bandstofit]

        for b in self.bandstofit:
            fr = self.fres[b.name].frames
            ints_comb = self.fres[b.name].intensities_combined[fr]
            m = np.invert(np.isnan(ints_comb))
            self.c_he_trace[b.name].setData(x=fr[m], y=ints_comb[m])

    def show_info(self):
        """ Show SIF image info """
        # show image info in browser
        txt = """
              <font size = {fs}>frames: {frames}</font><br>
              <font size = {fs}>T = {temp}</font><br>
              <font size = {fs}>exposure: {exposure}</font><br>
              <font size = {fs}>cycle: {cycle}</font><br>
              <font size = {fs}>xbin: {xbin}</font><br>
              <font size = {fs}>ybin: {ybin}</font><br>
              """
        txt = txt.format(
            frames=self.em.info[NOF],
            temp=self.em.info["DetectorTemperature"],
            exposure=self.em.info["ExposureTime"],
            cycle=self.em.info["CycleTime"],
            xbin=self.em.info["xbin"],
            ybin=self.em.info["ybin"],
            fs=5,  # font size
        )

        self.image_info_bw.setText(txt)

    def show_imageframe(self):
        """ Display current frame of the loaded SIF file """
        frame = int(self.frame.value())
        self.frame_current = frame
        self.img.setImage(self.em.images[frame])
        # Adjust hystogram
        mean = self.em.images[frame].mean()
        min = self.em.images[frame].min()
        self.hist.setLevels(mean * 1.5, min)
        self.hist.setHistogramRange(mean * 1.5, min)
        self.hist.regionChanged()

        self.show_spectrum()

    def show_spectrum(self):
        """ Show the current frame spectrum in counts and calibrated """
        frame = int(self.frame.value())
        self.spec_counts.setData(
            x=self.spectra.wavelength, y=self.spectra.counts[frame]
        )

        # show calibrated spectrum
        # str() - for PyQt4, to convert Qtstring int py string
        units = str(self.spec_units.currentText())

        txt = u'<font color = "#d62300" size="6">{}</font>'
        txt = txt.format(self.spec_labels[units])
        self.p3.setLabel("left", text=txt)
        spec = self.spectra.spectra_to_save[units]
        self.spec_wm.setData(x=self.spectra.wavelength, y=spec[frame])

    def show_balmerframe(self):
        """ Update the plots for H-balmer series experimental spectra
        """
        frame = int(self.frame_h.value())
        self.frame_current = frame

        if frame >= self.em.info["NumberOfFrames"]:
            return

        bands = self.hbands + [self.chband]
        kw = [b.name for b in self.hbands] + ["CH"]
        for b, k in zip(bands, kw):
            if self.fres[b.name].x is None:
                x, y = eb.banddata(b, frame, self.spectra, doblin=False)
                self.h_data[k].setData(x=x, y=y)
            else:
                x = self.fres[b.name].x
                y = self.fres[b.name].spectra[frame - 1]
                self.h_data[k].setData(x=x, y=y)

    def show_cframe(self):
        """ Update the plots for Carbon experimental spectra
        """
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

    def show_heframe(self):
        """ Update the plots for Carbon experimental spectra
        """
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

    def clear_traces(self):
        """ Clear trace plots """
        [self.htraces[b.name].setData(x=[], y=[]) for b in self.hbands]
        self.chtrace.setData(x=[], y=[])
        [self.c_he_trace[b.name].setData(x=[], y=[]) for b in self.bandstofit]

    def clear_fit_traces(self):
        """ clear traces of the fitted lines """
        [self.c_he_trace[b.name].setData(x=[], y=[]) for b in self.bandstofit]

    # ===========================================================================
    #             ^^^   Display Plots and info   ^^^
    # ===========================================================================

    def _enablecontrols(self, enable=True):
        """ Enable or disable controls
        """
        [i.setEnabled(enable) for i in self.controls_open]
        [i.setEnabled(enable) for i in self.controls_reg]
        self.setAcceptDrops(enable)

    def _reset_frame(self):
        """ Reset frame on all tabs
        """
        frame = self.frame_current
        [gi.setMaximum(1) for gi in self.frame_spinners]
        mainwindow.gui_setup_spinbox(self.frame, 0, 0, 0)
        self.frame_current = frame

    def _setup_frame(self):
        """ Setup frames on all tabs
        """
        mainwindow.gui_setup_spinbox(self.frame, 0, 0, self.em.info[NOF] - 1)
        [gi.setMaximum(self.em.info[NOF] - 1) for gi in self.frame_spinners]

        if self.frame_current > int(self.em.info[NOF]) - 1:
            self.frame.setValue(0)
        else:
            [f.setValue(self.frame_current) for f in [self.frame] + self.frame_spinners]

    def emitted(self):
        self.loadimage()

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


class LoadCalibrations(QtCore.QThread):
    """ Load Calibration files in a separate thread to release GUI """

    taskFinished = QtCore.pyqtSignal(object)

    def __init__(self, calibration):
        QtCore.QThread.__init__(self)
        self.calibration = calibration

    def run(self):
        """ convert Echelle image intol spectra """
        self.calibration.start()
        cb = self.calibration

        # plt.plot(cb.wavelength,cb.absolute['wm'])
        if DEBUG:
            print("\n" + cb.name, cb.DIMO, cb.DIMW)
            [print(a, " " * (12 - len(a)), b) for a, b in cb.filenames.items()]

        self.taskFinished.emit(self.calibration)


class LoadImage(QtCore.QThread):
    """ Load Sif Image in a separate thread to release GUI
    """

    taskFinished = QtCore.pyqtSignal(object)

    def __init__(self, arg):
        QtCore.QThread.__init__(self)
        self.arg = arg

    def run(self):
        """ convert Echelle image into spectra """
        filename, cb = self.arg
        try:
            em = ech.EchelleImage(filename, clbr=cb)
        except Exception as err:
            msg = "loadimage ERROR: {}".format(err)
            self.taskFinished.emit(None)
            # if DEBUG: print(msg)
            print(msg)
            return

        if DEBUG:
            print(em.info["DetectorDimensions"], cb.DIMO, cb.DIMW, cb.name)

        if not em.info["DetectorDimensions"] == (cb.DIMW, cb.DIMO):
            self.taskFinished.emit(None)
            return

        em.calculate_order_spectra()  # image -> order spectra
        em.correct_order_shapes()  # remove out of bounds boundaries
        em.calculate_spectra()  # order spectra -> fullwidth spectra

        self.taskFinished.emit(em)


class FitLines(QtCore.QThread):
    """ Fit lines (CIV, He, H) in a separate thread to release GUI """

    passresult = QtCore.pyqtSignal(object)
    progress = QtCore.pyqtSignal(int)

    def __init__(self, arg):
        QtCore.QThread.__init__(self)
        self.arg = arg

    def run(self):
        """ convert Echelle image into spectra
        """
        spectra, bands, fres, frames = self.arg
        [fres[b.name].set_frames(frames) for b in bands]
        [fres[b.name].init_spectra() for b in bands]

        for frame in frames:
            for b in bands:
                res = b.fitb(frame, spectra, detect=True, doblin=False)
                fres[b.name].spectra.append(res["exp"][1])
                fres[b.name].x = res["exp"][0]
                try:
                    fres[b.name].spectra_fit[frame] = res["fit"][1]
                    fres[b.name].x_fit = res["fit"][0]
                    fres[b.name].fit_out[frame] = res["out"]
                except TypeError:
                    pass

            self.progress.emit(frame + 1)

        [fres[band.name].integrate_spectra() for band in bands]
        [fres[band.name].integrate_fit() for band in bands]

        self.passresult.emit(fres)


class SaveSpectra(QtCore.QThread):
    """ Save spectra in separate thread to keep timing in line """

    passresult = QtCore.pyqtSignal(object)

    def __init__(self, arg):
        QtCore.QThread.__init__(self)
        self.spectra = arg

    def run(self):
        """ save spectra """

        self.spectra.save()
        self.passresult.emit(self.spectra)


if __name__ == "__main__":
    import ctypes

    myappid = "mycompany.myproduct.subproduct.version"  # arbitrary string
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    except:
        pass
    import sys

    app = QtGui.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon("icons/echelle.png"))
    win = EchelleAnalizer()
    win.show()
    sys.exit(app.exec_())
