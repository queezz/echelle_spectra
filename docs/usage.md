# Basic Usage

!!! note
    Detailed usage examples and API reference will be added in a future documentation pass. For now, see the example notebooks in `examples/workflows/black_cmos/`.

---

## Intended workflow

The general processing pipeline for an echelle spectrum is:

1. **Load a 2D echelle image**  
   Read the raw detector image from file (FITS, TIFF, or camera-native format).

2. **Define or load the order pattern**  
   Either load a previously saved pattern calibration, or run order detection on a suitable image (flat-field or bright continuum source). This step is only needed when the optical setup changes.

3. **Extract order spectra**  
   Apply the pattern to slice each diffraction order from the 2D image into a 1D intensity array.

4. **Calibrate wavelength**  
   Map pixel positions to wavelengths using a calibration solution (polynomial fit derived from known emission lines). Again, this is only redone when the setup changes.

5. **Inspect or analyze spectral windows**  
   Use the defined wavelength windows (see [Band data](band-data.md)) to locate and integrate lines of interest such as H-alpha, He-587, CII-515, etc.

---

## GUI

Launch the graphical interface with:

```bash
echelle_spectra
```

The GUI wraps the same extraction and calibration pipeline in an interactive window.

---

## Notebooks

For scripted or batch workflows, see the example notebooks:

```
examples/workflows/black_cmos/
├── 01_load_image.ipynb
├── 02_pattern_calibration.ipynb
├── 03_wavelength_calibration.ipynb
└── 04_extract_spectrum.ipynb    ← routine daily use
```
