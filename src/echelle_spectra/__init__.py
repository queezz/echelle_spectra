import pathlib

from .tools.config_loader import load_config
from .tools import echelle
from .tools import emissionbands
from .tools import emissiondata

__version__ = "0.0.2"

_echelle_base = pathlib.Path(__file__).parent.absolute()
_config = load_config(_echelle_base)

del pathlib
