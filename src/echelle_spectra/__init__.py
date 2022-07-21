import sys
import pathlib

# temporarily add this module's directory to PATH
_echelle_base = pathlib.Path(__file__).parent.absolute()
sys.path.append(str(_echelle_base))

# module-level imports
import echelle_spectra_gui as gui
from tools import echelle
from tools import emissionbands
from tools import emissiondata
from tools.config_loader import load_config

# echelle_spectra version
__version__ = "0.0.1"

# other useful module-level variables
_config = load_config(_echelle_base)
