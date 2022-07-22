<h1 align="center">
  <br>
  <img src="src/echelle_spectra/resources/graphics/echelle.png" alt="Echelle Spectra Logo" width="200">
  <br>
  Echelle Spectra
  <br>
</h1>

<h4 align="center">Graphical tool for the extraction and analysis of calibrated spectra from 2D Echelle spectrometer images</h4>

<p align="center">
  <a href="https://www.python.org/downloads/release/python-395">
    <img src="https://img.shields.io/badge/python-3.9-brigtgreen.svg" alt="Python 3.7">
  </a>
  
  <a href="https://github.com/queezz/echelle_spectra/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/queezz/echelle_spectra" alt="MIT License">
  </a>
  
  <a href="https://github.com/queezz/echelle_spectra/releases/latest">
    <img src="https://img.shields.io/github/v/release/queezz/echelle_spectra?include_prereleases&sort=semver"
         alt="Latest release">
  </a>
  
  <a href="https://github.com/ajulik1997/queezz/echelle_spectra/latest">
    <img src="https://img.shields.io/github/release-date-pre/queezz/echelle_spectra" alt="Latest release date">
  </a>
  
  <a href="https://github.com/queezz/echelle_spectra/commits">
    <img src="https://img.shields.io/github/commits-since/queezz/echelle_spectra/latest" alt="Commits since latest release">
  </a>
</p>

<p align="center">
  <a href="#quick-start">Quick start</a> •
  <a href="#usage">Usage</a> •
  <a href="#configuration">Configuration</a> •
  <a href="#authors">Authors</a> •
  <a href="#license">License</a>
</p>

<p float="left">
  <img src="images/gui_image.png" width="500" />
  <img src="images/gui_civ.png" width="500" />
</p>

## Quick start

The echelle_spectra module can be installed directly from PyPI via `pip`:

`pip install echelle_spectra --upgrade`

The echelle_spectra app GUI can be invoked directly as a Python module:

`python -m echelle_spectra`

Alternatively, the module can imported and the GUI started from within a Python script:

```python
import echelle_spectra
echelle_spectra.gui.start()
```

If you instead wish to run the tool from its source, clone this repository using `git`, install the project requirements, and execute the `__main__.py` script file as follows:

```shell
git clone https://github.com/queezz/echelle_spectra.git
cd echelle_spectra
pip install -r requirements.txt
python ./src/echelle_spectra/__main__.py
```

## Usage

With the Echelle Spectra GUI open, you will be presented with a set of blank image and spectra plots, as well as a control panel of settings and information on the left.

Firstly, an image file containing data from an Echelle spectrometer needs to be loaded. During the loading process, each frame of the file is converted to a calibrated wavelength-intensity spectrum, and fitting to regions of interest in the spectrum is carried out if requested. Optionally, all spectral and fitting data is saved to disk, which can be imported by another application.

After a file is successfully loaded, the first frame will be graphically displayed in the image viewer, and the corresponding computed spectrum will be shown below it. The control panel on the left can be used to select which frame from the file is being visualised.

Other tabs are also available in the GUI, which contain plots populated by snippets of spectral data in common regions of interest. If fitting was enabled during loading, a convolution of one or more Gaussian curves is plotted that best approximates the emissive intensity in those wavelength regions.

### Data loading

(Info box)

### Calibration

### Data export

### Bulk processing

# GUI for Echelle image processing
To convert Echelle images to spectra two options are avaliable: use this GUI, or import `EchelleImage` class to read, calibrate, and produce a `Spectrum` class, from wicht data could be exported. See examples for [CCD](src/echelle_spectra/examples/testtool-CCD.ipynb) and [CMOS](src/echelle_spectra/examples/testtool-CMOS.ipynb) sensors.
| CCD image                           | CMOS image                           |
| --------------------------------    | --------------------------------     |
| ![UI](src/echelle_spectra/examples/CCD_cut.png) | ![UI](src/echelle_spectra/examples/CMOS_cut.png) |


 ![UI](images/gui.png)