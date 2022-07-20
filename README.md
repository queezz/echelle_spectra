<h1 align="center">
    <br>
    <img src="src/echelle_spectra/resources/graphics/echelle.svg" alt="echelle_spectra logo" width="120">
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

# GUI for Echelle image processing
To convert Echelle images to spectra two options are avaliable: use this GUI, or import `EchelleImage` class to read, calibrate, and produce a `Spectrum` class, from wicht data could be exported. See examples for [CCD](src/echelle_spectra/examples/testtool-CCD.ipynb) and [CMOS](src/echelle_spectra/examples/testtool-CMOS.ipynb) sensors.
| CCD image                           | CMOS image                           |
| --------------------------------    | --------------------------------     |
| ![UI](src/echelle_spectra/examples/CCD_cut.png) | ![UI](src/echelle_spectra/examples/CMOS_cut.png) |


 ![UI](images/gui.png)