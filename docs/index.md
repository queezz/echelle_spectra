# Echelle Spectra

**Echelle Spectra** is a tool for extracting and analyzing calibrated spectra from 2D echelle spectrometer images.

For copy-paste commands, start with the [operator cheat
sheet](operator-cheat-sheet.md). It distinguishes portable-kit commands from a
development checkout and lists every installed GUI and CLI surface.

It provides a graphical interface for interactive use, but the underlying workflow can also be driven from Jupyter notebooks or scripts — useful for batch processing or automated pipelines.

---

## What it does

- Loads 2D echelle spectrometer images (e.g. from Andor CMOS cameras)
- Detects and fits diffraction order patterns on the detector
- Extracts per-order spectra
- Applies wavelength calibration
- Allows inspection and analysis of spectral windows for lines of interest (H, He, C, CH bands, etc.)

---

## Quick start

```bash
git clone https://github.com/queezz/echelle_spectra.git
cd echelle_spectra
pip install -e .
echelle_spectra          # launch GUI
```

---

## CMOS workflow notebooks

Notebooks for the Andor CMOS (2560×2160) setup are in `examples/workflows/black_cmos/`:

| Notebook | Purpose | Frequency |
|---|---|---|
| `01_load_image.ipynb` | Sanity check: load image, verify dimensions | as needed |
| `02_automated_pattern_extraction.ipynb` | Run packaged pattern extraction and compare traces | rare — only if optics moved |
| `02_pattern_calibration.ipynb` | Manual/tuning reference for pattern extraction | rare — setup/debugging |
| `03_wavelength_calibration.ipynb` | Manual line ID, polynomial fit | rare — new setup only |
| `04_extract_spectrum.ipynb` | Load → calibrate → extract → save | **routine** |
| `05_calibration_alignment.ipynb` | Align existing wavelength table with new Ne lamp data | session recalibration |

For daily use only `04_extract_spectrum.ipynb` is needed.

---

## Documentation history

The project previously had a minimal static HTML documentation site on the `html-docs` branch. That branch is preserved for reference. This MkDocs site is the new home for documentation going forward.
