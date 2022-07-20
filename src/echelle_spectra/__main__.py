"""Graphical tool for the extraction and analysis of calibrated spectra from 2D Echelle spectrometer images"""

import ctypes
import shutil
import sys
import tomli

from PyQt5.QtWidgets import QApplication
from pyqtgraph.Qt import QtGui

from __init__ import __version__
from echelle_spectra import EchelleSpectra


def load_config():
    with open("./config.toml", "rb") as cf:
        return tomli.load(cf)


try:
    config = load_config()
except (tomli.TOMLDecodeError, OSError):
    print("Warning: restoring default config file")
    shutil.copy("./resources/defaults.toml", "./config.toml")
    config = load_config()

try:
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(f"echelle_spectra-{__version__}")
except OSError:
    pass

app = QApplication(sys.argv)
app.setWindowIcon(QtGui.QIcon("./resources/graphics/echelle.png"))
win = EchelleSpectra(config)
win.show()
sys.exit(app.exec_())
