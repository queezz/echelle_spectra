#! /usr/bin/python

"""

"""
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.dockarea import DockArea, Dock
from PyQt5.QtWidgets import QMainWindow, QTabWidget, QPushButton
from PyQt5 import QtWidgets


def gui_setup_spinbox(spbox, val, min, max):
    """set up spinbox value, min and max """
    spbox.setMinimum(min)
    spbox.setMaximum(max)
    spbox.setValue(val)


def html_snc(text, color="red", size=6):
    """ change text color for html output """
    txt = u'<font color = "{c}" size="{s}">{t}</font>'
    return txt.format(c=color, t=text, s=size)


class Ui_MainWindow(object):
    def setupUi(self, QMainWindow):
        QMainWindow.setObjectName("MainWindow")
        # Interpret image data as row-major instead of col-major
        pg.setConfigOptions(imageAxisOrder="row-major")
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.centralwidget = QMainWindow
        self.tabwidg = QTabWidget()
        self.area = DockArea()
        self.area_civ = DockArea()
        self.area_he = DockArea()
        self.area_h = DockArea()
        self.area_cmd = DockArea()
        self.tabwidg.addTab(self.area, "Image")
        self.tabwidg.addTab(self.area_civ, "CIV")
        self.tabwidg.addTab(self.area_he, "He")
        self.tabwidg.addTab(self.area_h, "H2")
        self.tabwidg.addTab(self.area_cmd, "CMD")
        QMainWindow.setCentralWidget(self.tabwidg)

        QMainWindow.setWindowTitle(u"Echelle viewer")
        QMainWindow.resize(1000, 900)
        QMainWindow.statusBar().showMessage(" ")
        QMainWindow.setAcceptDrops(True)
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # DOCs
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.d1 = Dock("Image", size=(300, 400))
        self.d2 = Dock("Spectrum", size=(200, 300))
        self.d3 = Dock("Controls", size=(50, 20))
        self.d4 = Dock("Info", size=(60, 20))
        self.d5 = Dock("Register", size=(60, 20))
        self.d6 = Dock("Traces", size=(300, 400))

        [dd.setAcceptDrops(True) for dd in [self.d1, self.d2, self.d3]]
        self.area.addDock(self.d6, "top")
        self.area.addDock(self.d1, "below", self.d6)
        self.area.addDock(self.d4, "left")
        self.area.addDock(self.d3, "top", self.d4)
        self.area.addDock(self.d5, "top", self.d4)
        self.area.addDock(self.d2, "bottom", self.d1)
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.w1 = pg.GraphicsLayoutWidget()
        self.w2 = pg.GraphicsLayoutWidget()
        self.w_tr = pg.GraphicsLayoutWidget()
        self.d1.addWidget(self.w1)
        self.d2.addWidget(self.w2)
        self.d6.addWidget(self.w_tr)
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # Controls
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.btn_open = QPushButton("Manual SIF load")
        self.lbl_show = QtWidgets.QLabel("Selected shot:")
        self.btn_open.setToolTip("Manually open single SIF image")
        self.show_btn = QtWidgets.QPushButton("Load selected")
        self.show_btn.setToolTip("Load image from currently selected shot number")
        self.shot_number = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.shot_number, 139529, 1, 599999)
        self.shot_number.setToolTip("Shot number")
        self.last_shot_btn = QtWidgets.QPushButton("Select last")
        self.last_shot_btn.setToolTip("Select last shot from data folder")
        self.check_autoscale = QtWidgets.QCheckBox("Autoscale")
        msg = "Autoscale spectrum when hovering over image with Ctrl depressed"
        self.check_autoscale.setToolTip(msg)
        self.check_autoscale.setChecked(True)

        self.CameraCCD = QtWidgets.QRadioButton()
        self.CameraCCD.setText("CCD")
        self.CameraCCD.setToolTip("Calibration for CCD camera selected")
        self.CameraCMOS = QtWidgets.QRadioButton()
        self.CameraCMOS.setText("CMOS")
        self.CameraCMOS.setToolTip("Calibration for CMOS camera selected")

        self.lbl_units = QtWidgets.QLabel("Spectra units:")
        self.spec_units = QtWidgets.QComboBox()
        [self.spec_units.addItem(i) for i in ["counts", "wm", "wmsr", "phmsr"]]
        self.spec_units.setToolTip("Units of the spectra")
        self.spec_units.setCurrentIndex(1)

        self.specsave_bx = QtWidgets.QCheckBox("Save spec")
        msg = "Save spectra to output folder, always overwrite"
        self.specsave_bx.setToolTip(msg)
        self.fit_lines_bx = QtWidgets.QCheckBox("Fit spec")
        self.fit_lines_bx.setToolTip("Fit spec")
        self.save_lines_bx = QtWidgets.QCheckBox("Save lines")
        self.save_lines_bx.setToolTip("Save lines to output folder")
        self.overwrite = QtWidgets.QCheckBox("Overwrite all existing output")
        self.overwrite.setToolTip("Overwrite all existing spectra and intensity data")
        self.progress_bands = QtWidgets.QProgressBar()
        self.progress_bands.setRange(0, 1)
        self.progress_bands.setToolTip("Fitting lines")

        self.lbl_frame = QtWidgets.QLabel("Selected frame:")
        self.frame = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.frame, 0, 0, 0)
        self.frame.setToolTip("Selected frame number of loaded image")

        self.frame_civ = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.frame_civ, 1, 1, 1)

        self.coursor_bw = QtWidgets.QTextBrowser()
        self.coursor_bw.setMaximumHeight(65)
        self.image_info_bw = QtWidgets.QTextBrowser()
        # == == == == == == == == == == == == == == == == == == == == == == ====

        self.w3 = pg.LayoutWidget()
        self.w3.addWidget(self.btn_open, 0, 0, 1, 2)
        self.w3.addWidget(self.lbl_show, 1, 0)
        self.w3.addWidget(self.shot_number, 1, 1)
        self.w3.addWidget(self.last_shot_btn, 2, 0)
        self.w3.addWidget(self.show_btn, 2, 1)
        self.w3.addWidget(self.lbl_frame, 3, 0)
        self.w3.addWidget(self.frame, 3, 1)
        self.w3.addWidget(self.lbl_units, 4, 0)
        self.w3.addWidget(self.spec_units, 4, 1)
        self.w3.addWidget(self.check_autoscale, 5, 0)
        self.w3.addWidget(self.specsave_bx, 5, 1)
        self.w3.addWidget(self.fit_lines_bx, 6, 0)
        self.w3.addWidget(self.save_lines_bx, 6, 1)
        self.w3.addWidget(self.overwrite, 7, 0, 1, 2)
        self.w3.addWidget(self.CameraCCD, 8, 0)
        self.w3.addWidget(self.CameraCMOS, 8, 1)
        self.w3.addWidget(self.progress_bands, 9, 0, 1, 2)
        self.w3.addWidget(self.coursor_bw, 10, 0, 1, 2)

        self.controls_open = [
            self.btn_open,
            self.shot_number,
            self.last_shot_btn,
            self.show_btn,
            self.overwrite
        ]

        self.d3.addWidget(self.w3)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.w4 = pg.LayoutWidget()
        self.w4.addWidget(self.image_info_bw)
        self.d4.addWidget(self.w4)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        # FIT Gaussian or Lorentz in the Spectrum plot
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.wfit = pg.LayoutWidget()

        self.region_show = QtWidgets.QPushButton("fit tool")
        self.btn_region_hide = QtWidgets.QPushButton("hide tool")
        self.fit_model = QtWidgets.QComboBox()
        self.fit_model.addItem("gauss")
        self.fit_model.addItem("gauss2")
        self.fit_model.addItem("lorentz")
        self.btn_fit = QtWidgets.QPushButton("fit")
        self.fit_bw = QtWidgets.QTextBrowser()

        self.region_fit = pg.LinearRegionItem(values=[0, 1])
        self.region_fit.setBrush(QtGui.QColor(249, 247, 232, 50))

        self.wfit.addWidget(self.region_show, 0, 0)
        self.wfit.addWidget(self.btn_region_hide, 0, 1)
        self.wfit.addWidget(self.btn_fit, 1, 0)
        self.wfit.addWidget(self.fit_model, 1, 1)
        self.wfit.addWidget(self.fit_bw, 11, 0, 1, 2)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        # IMAGE and plots
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # A plot area (ViewBox + axes) for displaying the image
        self.p1 = self.w1.addPlot()  # image plot
        self.p2 = self.w2.addPlot()  # spectrum plot
        self.p3 = self.w2.addPlot(row=1, col=0)

        # Plots for intensities
        # - - - - - - -
        self.h_tr = self.w_tr.addPlot(row=0, col=0)
        self.ch_tr = self.w_tr.addPlot(row=1, col=0)
        self.c_tr = self.w_tr.addPlot(row=2, col=0)
        # - - - - - - -
        self.h_tr.setLabel("left", text=html_snc("H balmer", "#a86ff7", size=6))
        self.ch_tr.setLabel("left", text=html_snc("CH", "#2b71e2", size=6))
        self.c_tr.setLabel("left", text=html_snc("C and He", "#2ae286", size=6))
        # - - - - - - -
        self.p2.setLabel("left", text=html_snc("counts", "#6ac600", size=6))

        slabels = [
            "counts",
            "W/(m<sup>2</sup> nm)",
            "W/(m<sup>2</sup> sr nm)",
            "N<sub>ph.</sub>/(m<sup>2</sup> sr nm)",
        ]
        skeys = ["counts", "wm", "wmsr", "phmsr"]
        self.spec_labels = {i: j for i, j in zip(skeys, slabels)}

        txt = html_snc("W/(m<sup>2</sup> nm)", "#d62300", size=6)
        self.p3.setLabel("left", text=txt)
        self.p3.setXLink(self.p2)
        # Item for displaying image data
        self.img = pg.ImageItem()
        self.p1.addItem(self.img)
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # Contrast/color control
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img)
        self.w1.addItem(self.hist)
        # hist.vb.setMouseEnabled(y=False)#makes user interaction a little easier

        self.fit_model.setEnabled(False)
        self.btn_fit.setEnabled(False)
        pen = pg.mkPen("#11d7ff", width=2, style=QtCore.Qt.DashLine)
        self.cL = pg.InfiniteLine(angle=90, movable=False, pen=pen)
        self.cL1 = pg.InfiniteLine(angle=90, movable=True, pen="b")
        self.cL2 = pg.InfiniteLine(angle=90, movable=True, pen="y")

        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # REGISTER DATA
        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.wreg = pg.LayoutWidget()

        self.shot_range_btn = QtWidgets.QPushButton("Get available shot range")
        msg = "Get the highest and lowest available shot numbers from the data folder"
        self.shot_range_btn.setToolTip(msg)
        self.start_shot = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.start_shot, 100000, 1, 599999)
        self.start_shot.setToolTip("Start shot")
        self.end_shot = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.end_shot, 599999, 1, 599999)
        self.end_shot.setToolTip("End shot")
        self.trigger_delay_lbl = QtWidgets.QLabel("Trigger delay:")
        self.trigger_delay = QtWidgets.QDoubleSpinBox()
        gui_setup_spinbox(self.trigger_delay, 2.75, 0.0, 1000.0)
        self.trigger_delay.setSuffix(" s")
        self.trigger_delay.setSingleStep(0.05)
        self.trigger_delay.setToolTip("Trigger delay (s)")
        self.start_register = QtWidgets.QPushButton("Start")
        msg = "Calculate spectra and line intensities for the shot range"
        self.start_register.setToolTip(msg)
        self.abort_register = QtWidgets.QPushButton("Abort")
        msg = "Abort the loop after current image is analyzed"
        self.abort_register.setToolTip(msg)
        self.progress_range = QtWidgets.QProgressBar()
        self.progress_range.setToolTip("Loop trough available shot numbers")

        self.wreg.addWidget(self.shot_range_btn, 0, 0, 1, 2)
        self.wreg.addWidget(self.start_shot, 1, 0)
        self.wreg.addWidget(self.end_shot, 1, 1)
        self.wreg.addWidget(self.trigger_delay_lbl, 2, 0)
        self.wreg.addWidget(self.trigger_delay, 2, 1)
        self.wreg.addWidget(self.start_register, 3, 0)
        self.wreg.addWidget(self.abort_register, 3, 1)
        self.wreg.addWidget(self.progress_range, 4, 0, 1, 2)

        self.controls_reg = [
            self.shot_range_btn,
            self.start_shot,
            self.end_shot,
            self.trigger_delay,
            self.start_register,
        ]

        self.d5.addWidget(self.wreg)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # CIV tab setup
        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.d_civ_1 = Dock("spectra", size=(300, 400))
        self.d_civ_2 = Dock("controls", size=(100, 400))

        self.area_civ.addDock(self.d_civ_1)
        self.area_civ.addDock(self.d_civ_2, "left")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.w_civ_controls = pg.LayoutWidget()

        self.frame_civ = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.frame_civ, 1, 1, 1)
        self.frame_civ.setSuffix(" frame")
        self.bw_civ = QtWidgets.QTextBrowser()

        self.w_civ_controls.addWidget(self.frame_civ, 0, 0)
        self.w_civ_controls.addWidget(self.bw_civ, 10, 0, 2, 1)

        self.d_civ_2.addWidget(self.w_civ_controls)
        # - - - - - - -
        self.w_civ_plots = pg.GraphicsLayoutWidget()
        self.w_civ_plots.setBackground(background="k")

        self.d_civ_1.addWidget(self.w_civ_plots)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # He tab setup
        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.d_he_1 = Dock("spectra", size=(300, 400))
        self.d_he_2 = Dock("controls", size=(100, 400))

        self.area_he.addDock(self.d_he_1)
        self.area_he.addDock(self.d_he_2, "left")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.w_he_controls = pg.LayoutWidget()

        self.frame_he = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.frame_he, 1, 1, 1)
        self.frame_he.setSuffix(" frame")
        self.bw_he = QtWidgets.QTextBrowser()

        self.w_he_controls.addWidget(self.frame_he, 0, 0)
        self.w_he_controls.addWidget(self.bw_he, 10, 0, 2, 1)

        self.d_he_2.addWidget(self.w_he_controls)
        # - - - - - - -
        self.w_he_plots = pg.GraphicsLayoutWidget()
        self.w_he_plots.setBackground(background="k")

        self.d_he_1.addWidget(self.w_he_plots)

        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        # Hydrogen tab setup
        # == == == == == == == == == == == == == == == == == == == == == == ====
        #
        # == == == == == == == == == == == == == == == == == == == == == == ====
        self.d_h_1 = Dock("spectra", size=(300, 400))
        self.d_h_2 = Dock("controls", size=(100, 400))

        self.area_h.addDock(self.d_h_1)
        self.area_h.addDock(self.d_h_2, "left")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.w_h_controls = pg.LayoutWidget()

        self.frame_h = QtWidgets.QSpinBox()
        gui_setup_spinbox(self.frame_h, 1, 1, 1)
        self.frame_h.setSuffix(" frame")
        self.bw_h = QtWidgets.QTextBrowser()

        self.w_h_controls.addWidget(self.frame_h, 0, 0)
        self.w_h_controls.addWidget(self.bw_h, 10, 0, 2, 1)

        self.d_h_2.addWidget(self.w_h_controls)
        # - - - - - - -
        self.w_h_plots = pg.GraphicsLayoutWidget()
        self.w_h_plots.setBackground(background="k")

        self.d_h_1.addWidget(self.w_h_plots)
        # == == == == == == == == == == == == == == == == == == == == == == == =
        # CMD tab setup
        # == == == == == == == == == == == == == == == == == == == == == == == =
        #
        # == == == == == == == == == == == == == == == == == == == == == == == =
        self.d_cmd = Dock("CMD out")
        self.area_cmd.addDock(self.d_cmd)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.w_cmd = pg.LayoutWidget()
        self.cmd_bw = QtWidgets.QTextBrowser()
        self.w_cmd.addWidget(self.cmd_bw)
        self.d_cmd.addWidget(self.w_cmd)
        # == == == == == == == == == == == == == == == == == == == == == == == =
        # END of GUI
        # == == == == == == == == == == == == == == == == == == == == == == == =

        civ_names = [
            "CII-515",
            "CIII-464",
            "CIV-580",
            "CIV-444",
            "CIV-465",
            "CIV-772",
            "CIV-706",
            "CIV-547",
        ]
        he_names = [
            "He-447",
            "He-492",
            "He-501",
            "He-504",
            "He-587",
            "He-667",
            "He-706",
            "He-728",
        ]
        h_names = ["H-alpha", "H-beta", "H-gamma", "H-delta"]
        self.h_names = h_names

        self.civ_p, self.civ_data, self.civ_fit = self.create_plots(
            civ_names, self.w_civ_plots
        )
        self.he_p, self.he_data, self.he_fit = self.create_plots(
            he_names, self.w_he_plots
        )
        self.h_p, self.h_data, self.h_fit = self.create_plots(h_names, self.w_h_plots)

        clrs = ["#ef0000", "#00efff", "#2500f0", "#6700b5"]
        clrs = {name: clr for name, clr in zip(h_names, clrs)}
        smbls = ["o", "s", "d", "x"]
        smbls = {name: smbl for name, smbl in zip(h_names, smbls)}
        sptext = "&nbsp;&nbsp;&nbsp;&nbsp;"
        kws = {
            name: {
                "pen": clrs[name],
                "symbol": smbls[name],
                "symbolPen": clrs[name],
                "symbolBrush": clrs[name],
                "name": sptext + name,
                "symbolSize": 8,
            }
            for name, clr in zip(h_names, clrs)
        }

        self.htraces = {name: self.h_tr.plot(**kws[name]) for name in h_names}

        clr = "#d9f442"
        kws = {
            "pen": clr,
            "symbol": "o",
            "symbolPen": clr,
            "symbolBrush": clr,
            "name": sptext + "CH",
            "symbolSize": 8,
        }
        self.chtrace = self.ch_tr.plot(**kws)

        self.c_he_trace = {}
        for i, name in enumerate(he_names + civ_names):
            if "He" in name:
                clr = "#d9f442"
                kws = {
                    "pen": clr,
                    "symbol": "s",
                    "symbolPen": clr,
                    "symbolBrush": clr,
                    "name": sptext + "CH",
                    "symbolSize": 8,
                }
            else:
                clr = "#6700b5"
                kws = {
                    "pen": clr,
                    "symbol": "o",
                    "symbolPen": clr,
                    "symbolBrush": clr,
                    "name": sptext + "CH",
                    "symbolSize": 8,
                }

            self.c_he_trace[name] = self.c_tr.plot(**kws)

    def create_plots(self, names, plotwidg):
        """ Create plots for Carbon, He or Hydrogen bands
        """
        bandplots = [
            plotwidg.addPlot(row=ii // 2, col=ii % 2) for ii, name in enumerate(names)
        ]
        bandplots = {name: p for name, p in zip(names, bandplots)}
        clrs = {name: "#43f967" for name in names}
        col = "r"
        symb = "s"
        pen = pg.mkPen("#11d7ff", width=2, style=QtCore.Qt.DashLine)
        tnames = names
        if "C" in names[0]:
            col = "#edc740"
            symb = "o"
            pen = pg.mkPen("#11d7ff", width=2, style=QtCore.Qt.DashLine)
            clrs["CII-515"] = "#e2d409"
            clrs["CIII-464"] = "#e2d409"

        if "He" in names[0]:
            col = "#edc740"
            symb = "x"
            pen = pg.mkPen("#11d7ff", width=2, style=QtCore.Qt.DashLine)

        if "H-" in names[0]:
            col = "#edc740"
            symb = "s"
            pen = pg.mkPen("#11d7ff", width=2, style=QtCore.Qt.DashLine)
            tnames = ["H&alpha;", "H&beta;", "H&gamma;", "H&delta;"]
            clrs["H-alpha"] = "#ef0000"
            clrs["H-beta"] = "#00efff"
            clrs["H-gamma"] = "#2500f0"
            clrs["H-delta"] = "#6700b5"
            bandplots["CH"] = plotwidg.addPlot(len(names) // 2 + 1, 0, 1, 2)
            # bandplots['trace'] = plotwidg.addPlot(len(names)//2+2,0,1,2)

        [
            bandplots[name].setTitle(html_snc(t, clrs[name], size=6))
            for name, t in zip(names, tnames)
        ]

        kws = {
            "pen": None,
            "symbol": symb,
            "name": "&nbsp;&nbsp;&nbsp;&nbsp;data",
            "symbolSize": 8,
            "symbolPen": col,
            "symbolBrush": col,
        }

        band_data = {name: bandplots[name].plot(**kws) for name in names}
        band_fit = {name: bandplots[name].plot(pen=pen) for name in names}
        if "H-" in names[0]:
            band_data["CH"] = bandplots["CH"].plot(**kws)
            bandplots["CH"].setTitle(html_snc("CH", "#e2d409", size=6))
        return bandplots, band_data, band_fit


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

