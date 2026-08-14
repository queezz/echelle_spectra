import pathlib

from .tools import echelle as echelle
from .tools import emissionbands as emissionbands
from .tools import emissiondata as emissiondata
from .tools.config_loader import load_config
from .tools.line_catalog import SpectralLine, load_line_table
from .tools.spectrocube_export import export_spectrocube as export_spectrocube
from .tools.spectrocube_export import to_spectrocube as to_spectrocube

__version__ = "1.5.0"

__all__ = [
    "echelle",
    "emissionbands",
    "emissiondata",
    "export_spectrocube",
    "load_line_table",
    "SpectralLine",
    "to_spectrocube",
]

_echelle_base = pathlib.Path(__file__).parent.absolute()
_config = load_config(_echelle_base)

del pathlib
