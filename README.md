
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

<h2 align="center">
  <br>
  <img src="src/echelle_spectra/resources/graphics/echelle.png" alt="Echelle Spectra Logo" width="60">
  <br>
  Echelle Spectra
  <br>
</h2>

<h4 align="center">Graphical tool for the extraction and analysis of calibrated spectra from 2D Echelle spectrometer images</h4>


<p align="center">
  <a href="#start-with-venv">Start with Venv</a> •
  <a href="#quick-start">Quick start</a> •
  <a href="#usage">Usage</a> •
  <a href="#configuration">Configuration</a> •
  <a href="#authors">Authors</a> •
  <a href="#license">License</a>
</p>


## Install with venv
There are plentiful ways to install a package, but for now we recommend [venv](https://docs.python.org/3/library/venv.html).
Clone the repo into your folder of choice (use `.` to put it inside current dir)
```bash
git clone https://github.com/queezz/echelle_spectra.git .
```

Install a venv locally in the project folder

```
python -m venv ./.venv/echelle
```

Activate it
```powershell
./.venv/echelle/Scripts/Activate.ps1
```

Update pip
```
python -m pip install -U pip
```

Install `echelle_spectra` in editable mode with dev tools:
```
python -m pip install -e ".[dev]"
```
Run the GUI program
```bash
python ./src/echelle_spectra/__main__.py
```

### External venv
If you use your code with some cloud, it may be better to keep venv oustide the cloud storage. Then put venvs somewhere central. On Windows:
`$env:USERPROFILE\.venvs\` is a goood place.
```shell
python -m venv $env:USERPROFILE\.venvs\echelle
```
Activate
```shell
& $env:USERPROFILE\.venvs\echelle\Scripts\Activate.ps1
```

## Quick start
Run GUI from python:
```python
import echelle_spectra
echelle_spectra.gui.start()
```
Or from the command prompt:
```shell
python ./src/echelle_spectra/__main__.py
```

## Usage

With the Echelle Spectra GUI open, you will be presented with a set of blank image and spectra plots, as well as a control panel of settings and information on the left.

Firstly, an image file containing data from an Echelle spectrometer needs to be loaded. During the loading process, each frame of the file is converted to a calibrated wavelength-intensity spectrum, and fitting to regions of interest in the spectrum is carried out if requested. Optionally, all spectral and fitting data is saved to disk, which can be imported by another application.

After a file is successfully loaded, the first frame will be graphically displayed in the image viewer, and the corresponding computed spectrum will be shown below it. The control panel on the left can be used to select which frame from the file is being visualised.

Other tabs are also available in the GUI, which contain plots populated by snippets of spectral data in common regions of interest. If fitting was enabled during loading, a convolution of one or more Gaussian curves is plotted that best approximates the emissive intensity in those wavelength regions.


# GUI

 ![UI](images/gui.png)
