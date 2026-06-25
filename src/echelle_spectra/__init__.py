import pathlib

from .tools.config_loader import load_config
from .tools import echelle
from .tools import emissionbands
from .tools import emissiondata
from .tools.spectrocube_export import export_spectrocube, to_spectrocube

__version__ = "0.2.0"

_echelle_base = pathlib.Path(__file__).parent.absolute()
_config = load_config(_echelle_base)

del pathlib
