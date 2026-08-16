# Echelle Spectra

**Echelle Spectra** is a tool for extracting and analyzing calibrated spectra from 2D echelle spectrometer images.

Three surfaces share one calibration and extraction pipeline: the **live
calibration bench** you use at the instrument, the **`echelle` campaign
commands** that convert whole drives, and the **single-SIF viewer** for looking
at one shot.

[Where to start](usage-overview.md) routes by task in one table. For copy-paste
commands, go to the [operator cheat sheet](operator-cheat-sheet.md): it
distinguishes portable-kit commands from a development checkout and lists every
installed surface.

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
echelle status           # what calibrations and runs exist
echelle_spectra          # open the single-SIF viewer
```

The primary workflow is **snapshot → registry → process**. See
[Where to start](usage-overview.md).

---

## CMOS workflow notebooks

The notebooks in `examples/workflows/black_cmos/` are a **manual and tuning
reference** for the Andor CMOS (2560×2160) setup — kept for a new optical setup
or a fit worth stepping through, not for routine work.

| Notebook | Purpose | When |
|---|---|---|
| `01_load_image.ipynb` | Sanity check: load image, verify dimensions | as needed |
| `02_automated_pattern_extraction.ipynb` | Run packaged pattern extraction and compare traces | rare — only if optics moved |
| `02_pattern_calibration.ipynb` | Manual pattern tuning and debugging | rare — setup/debugging |
| `03_wavelength_calibration.ipynb` | Manual line ID, polynomial fit | rare — new setup only |
| `04_extract_spectrum.ipynb` | One-file walkthrough: load → calibrate → extract → save | reference |
| `05_calibration_alignment.ipynb` | Align existing wavelength table with new Ne lamp data | superseded by the bench |

For routine calibration use the [live bench](calibration-bench.md); for routine
conversion use [`echelle process`](calibration-to-cube.md). Both record
provenance the notebooks do not.

---

## Documentation history

The project previously had a minimal static HTML documentation site on the `html-docs` branch. That branch is preserved for reference. This MkDocs site is the new home for documentation going forward.
