"""Graphical tool for the extraction and analysis of calibrated spectra from 2D Echelle spectrometer images"""

import ctypes
import shutil
import sys
import tomli

from pathlib import Path
from PyQt5.QtWidgets import QApplication
from pyqtgraph.Qt import QtGui

try:
    # script-level imports
    from __init__ import __version__
    from echelle_spectra import EchelleSpectra
except ModuleNotFoundError:
    # module-level imports
    from .__init__ import __version__
    from .echelle_spectra import EchelleSpectra


base_path = Path(__file__).parent.absolute()


def load_config(parent):
    with open(parent / "config.toml", "rb") as cf:
        return tomli.load(cf)


try:
    config = load_config(base_path)
except (tomli.TOMLDecodeError, OSError):
    print("Warning: restoring default config file")
    shutil.copy(base_path / "resources/defaults.toml", base_path / "config.toml")
    config = load_config(base_path)

try:
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(f"echelle_spectra-{__version__}")
except OSError:
    pass

app = QApplication(sys.argv)
app.setWindowIcon(QtGui.QIcon(base_path / "resources/graphics/echelle.png"))
win = EchelleSpectra(config)
win.show()
sys.exit(app.exec_())
