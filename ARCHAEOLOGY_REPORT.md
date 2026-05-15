# Echelle Spectra — Repository Archaeology Report

**Date:** 2026-05-15
**Scope:** Full codebase analysis, workflow tracing, calibration model documentation

---

## Table of Contents

1. [Repository Map](#1-repository-map)
2. [Architecture Overview](#2-architecture-overview)
3. [Core Module Analysis](#3-core-module-analysis)
4. [Calibration Model](#4-calibration-model)
5. [Workflow Map](#5-workflow-map)
6. [GUI Analysis](#6-gui-analysis)
7. [Example Classification](#7-example-classification)
8. [New Spectrometer Initialization](#8-new-spectrometer-initialization)
9. [Where Manual Recalibration Happens](#9-where-manual-recalibration-happens)
10. [Session Alignment Proposal](#10-session-alignment-proposal)
11. [God Objects and Oversized Files](#11-god-objects-and-oversized-files)
12. [Duplicated and Dead Code](#12-duplicated-and-dead-code)
13. [Proposed MkDocs Structure](#13-proposed-mkdocs-structure)
14. [Safe Incremental Improvements](#14-safe-incremental-improvements)
15. [Minimal-Risk Modernization Plan](#15-minimal-risk-modernization-plan)

---

## 1. Repository Map

```
echelle_spectra/
├── pyproject.toml              # Package config (setuptools, v0.0.2)
├── requirements.txt            # Flat deps (slightly out of sync with pyproject.toml)
├── MANIFEST.in                 # Package data inclusion rules
├── README.md                   # Install + quick start instructions
├── LICENSE                     # MIT
├── .gitignore
│
├── images/
│   └── gui.png                 # Screenshot for README
│
├── examples/                   # Jupyter notebooks + output PDFs + pattern data
│   ├── 1.1_Cutting_Pattern_Synthetic.ipynb/.pdf
│   ├── 1.2_Cutting_Pattern_Masks.ipynb/.pdf
│   ├── 2.1_Cutting_Pattern_FE.ipynb
│   ├── 2.2_Cutting_Pattern_NihamaPlus.ipynb
│   ├── 3.1_loading_images.ipynb
│   ├── 3.2_try_load_new_echelle.ipynb
│   ├── 3.3_Wavelength_Calibration.ipynb
│   ├── 3.4_Wavelength_echelle_spectra copy.ipynb
│   ├── 3.5_Load_Black_Echelle.ipynb
│   ├── 4.1_Load_bkack_echelle copy.ipynb
│   ├── 4.2_Pattern_LHD.ipynb
│   ├── 4.3_Wavelength_Calibration_LHD.ipynb
│   ├── 5.0_Bands.ipynb
│   ├── testtool-CCD.ipynb/.pdf
│   ├── testtool-CMOS.ipynb/.pdf
│   ├── pattern.txt                         # CCD pattern (exported)
│   ├── pattern_CMOS_20240305.txt           # CMOS pattern (exported)
│   ├── pattern_fujii.txt / pattern_fujii_rot.txt
│   ├── pattern_image_synthetic.txt         # 100 MB synthetic pattern image
│   ├── CCD_cut.png / CMOS_cut.png
│   └── *.pdf                               # Exported figures
│
└── src/echelle_spectra/
    ├── __init__.py             # Package init, sys.path hack, exports
    ├── __main__.py             # Entry point → gui.start()
    ├── echelle_spectra_gui.py  # Main GUI application (912 lines)
    │
    ├── tools/
    │   ├── echelle.py          # Core: EchelleImage, Calibrations, Spectrum (900 lines)
    │   ├── emissionbands.py    # EmissionBand fitting model (369 lines)
    │   ├── emissiondata.py     # Predefined emission band data (387 lines)
    │   ├── config_loader.py    # TOML config loader (24 lines)
    │   └── note.txt            # "Not cleaned up" notice
    │
    └── resources/
        ├── defaults.toml       # Default config (data_path, output_path, diag_name)
        ├── header_template.txt # Output file header format
        ├── wavelength_fujii.txt # Partial: fujii wavelength cal (H2 Fulcher lines)
        ├── window_layout.py    # Qt UI layout definition (546 lines)
        ├── graphics/
        │   ├── echelle.png     # App icon
        │   └── echelle.svg
        ├── test_data/
        │   └── CCD_Example.SIF # ~15 MB test image
        └── calibration_files/
            ├── pattern.txt                     # CCD order pattern (1024×28)
            ├── pattern_cmos.txt                # CMOS order pattern (older)
            ├── pattern_CMOS_20240305.txt        # CMOS order pattern (2560×29, current)
            ├── Th_wavelength.txt               # CCD wavelength table (28 orders, Th/Ar/Ne)
            ├── Th_wavelength_CMOS.txt          # CMOS wavelength table (older)
            ├── Th_wavelength_CMOS_20240305.txt  # CMOS wavelength table (current)
            ├── Th_wavelength_CMOS_only_Th-Ar.txt # CMOS, Th/Ar subset only
            ├── wavelength_cmos_2024.txt        # Near-identical to Th_wavelength_CMOS.txt
            ├── wavelength_fujii.txt            # Fujii spectrometer (H2 Fulcher only)
            ├── integrating_sphere.txt          # Sphere spectral irradiance table
            ├── absolute_20170613_b8.sif        # CCD sphere image
            ├── absolute_20170613_b8_0.2_bkg.sif # CCD sphere background
            ├── absolute_20170613_b8_0.2_v2.sif  # CCD sphere (v2, used)
            ├── sphere_CMOS.sif / sphere_CMOS_bkg.sif         # CMOS sphere (older)
            ├── sphere_cmos_20240305.sif / _bkg.sif            # CMOS sphere (current)
            └── ThAr_10.0s_16bit.sif / ThAr-0.3s-x3_20240305.sif  # Lamp images
```

### File Size Notes

| File | Size | Notes |
|------|------|-------|
| `pattern_image_synthetic.txt` | **100 MB** | In `examples/`; should probably be in `.gitignore` or LFS |
| `sphere_*.sif` | ~63 MB each | 4 files, ~252 MB total; necessary calibration data |
| `ThAr*.sif` | 22–66 MB | Lamp reference images |
| `CCD_Example.SIF` | 15 MB | Test data |

---

## 2. Architecture Overview

### Class Hierarchy

```
Calibrations              ← Owns calibration state for one detector config
  ├── .pattern            ← Order positions (numpy array from pattern.txt)
  ├── .cutting_masks      ← Pixel masks for each diffraction order
  ├── .order_wavel        ← Polynomial wavelength solution per order
  ├── .order_borders      ← Where adjacent orders overlap/join
  ├── .absolute           ← Sensitivity curve from integrating sphere
  ├── .sphr → EchelleImage  ← Integrating sphere image
  └── .bkgr → EchelleImage  ← Background image

EchelleImage              ← Single loaded detector image
  ├── .images             ← Raw pixel data (frames × rows × cols)
  ├── .info               ← Metadata dict (exposure, frames, binning, etc.)
  ├── .clbr → Calibrations  ← Back-reference to calibration
  ├── .order_spectra      ← Extracted 1D spectra per order
  └── .spectra            ← Stitched full-range spectrum

Spectrum                  ← Calibrated spectrum with physical units
  ├── .wavelength         ← Wavelength array
  ├── .counts             ← Background-subtracted counts
  ├── .wm / .wmsr / .phmsr  ← Physical unit spectra
  └── .absolute           ← Sensitivity from calibration

EmissionBand              ← Spectral line/band with lmfit Gaussian model
  ├── .wavelengths        ← Line center wavelengths (from NIST etc.)
  ├── .Ak / .stat_weights ← Transition probabilities
  ├── .model / .pars      ← lmfit composite Gaussian model
  └── .fitb()             ← Fit to experimental data

EB                        ← Simplified band (bounds + name, no atomic data)

FitResult                 ← Container for fit results across frames

EchelleSpectraGUI         ← PyQt5 main window
  ├── .cb_CCD / .cb_CMOS  ← Two Calibrations instances
  ├── .bands / .bandstofit ← Emission band definitions
  └── Threads: LoadCalibrations, LoadImage, FitLines, SaveSpectra
```

### Data Flow

```
SIF/TIFF file
    │
    ▼
read_image()  ←────── spec="black" (SIF/Andor) or spec="fujii" (TIFF/Hamamatsu)
    │
    ▼
EchelleImage(fpth, clbr=calibrations)
    │
    ├── calculate_order_spectra()  ←── uses clbr.cutting_masks
    ├── correct_order_shapes()     ←── handles partial orders at sensor edge
    └── calculate_spectra()        ←── uses clbr.order_borders to stitch
    │
    ▼
Spectrum(echelle_image)
    │
    ├── Background subtraction (auto-detected from multi-frame images)
    ├── Wavelength assignment from clbr
    ├── Absolute calibration (counts → W/m²/nm, W/m²/sr/nm, Nph/m²/sr)
    └── .save() → text file with header
```

---

## 3. Core Module Analysis

### `tools/echelle.py` (900 lines) — THE core module

| Component | Lines | Role |
|-----------|-------|------|
| `read_image()` | 27–81 | Dispatch: SIF→sif_parser, TIFF→PIL, handles binning/crop/transpose |
| `EchelleImage` | 84–377 | Image container + order extraction + plotting |
| `Calibrations` | 379–769 | **Central calibration class**: pattern loading, mask generation, wavelength calibration, absolute calibration |
| `Spectrum` | 772–876 | Calibrated spectrum: background subtraction, unit conversion, saving |
| `header_template` | 878–896 | Output file header string (duplicated with `resources/header_template.txt`) |

### `tools/emissionbands.py` (369 lines)

| Component | Role |
|-----------|------|
| `baseline_als()` | Asymmetric least squares baseline removal |
| `banddata()` | Extract wavelength/intensity slice for a band from a Spectrum |
| `EB` | Simple band container (bounds + name) |
| `EmissionBand` | Full atomic transition data + lmfit Gaussian model builder |
| `FitResult` | Multi-frame fit result container with integration |

### `tools/emissiondata.py` (387 lines)

Pure data file defining all emission bands used in the application:
- **CIV lines**: c444, c465, c547, c580, c706, c772
- **CII/CIII**: c515, c464
- **Xenon**: xe461, xe462, xe467
- **Helium**: he447, he492, he501, he504, he587, he667, he706, he728
- **Hydrogen Balmer**: halpha, hbeta, hgamma, hdelta
- **set_bounds()**: Adjusts wavelength windows for all bands; called at import time

### `tools/config_loader.py` (24 lines)

Loads `config.toml` from package directory; falls back to `defaults.toml` on error.

### `echelle_spectra_gui.py` (912 lines)

| Section | Lines | Role |
|---------|-------|------|
| `EchelleSpectraGUI.__init__` | 19–49 | Setup UI, connect signals, load calibrations, define bands |
| `prepare_calibration()` | 58–93 | Hardcoded CCD/CMOS filename dicts, creates two `Calibrations` |
| Loop system | 261–316 | Batch-process shot ranges |
| Image loading | 325–444 | File open, thread dispatch, dimension mismatch auto-retry |
| `no_fit_intensities()` | 456–493 | Non-fit intensity extraction for all bands |
| `fit_lines()` | 495–535 | Threaded Gaussian fitting |
| Display methods | 583–727 | Plot updates for image, spectrum, C/He/H tabs |
| `LoadCalibrationsThread` | 780–798 | Background calibration loading |
| `LoadImageThread` | 801–834 | Background image loading + extraction |
| `FitLinesThread` | 837–873 | Background line fitting |
| `SaveSpectraThread` | 876–887 | Background spectrum saving |

### `resources/window_layout.py` (546 lines)

Hand-coded PyQt5 layout using `pyqtgraph.dockarea`. Defines all widgets, tabs (Image, CIV, He, H2, CMD), plot areas, and trace configurations. Not auto-generated — purely manual.

---

## 4. Calibration Model

### 4.1 Order Pattern

**What it is:** A 2D integer array mapping pixel columns to the vertical center of each diffraction order on the detector.

**File format:** Space-separated text, no header.

| Detector | File | Shape | Orders |
|----------|------|-------|--------|
| CCD (Andor, 1024×1024) | `pattern.txt` | 1024 × 28 | 28 |
| CMOS (Andor, 2560×2160) | `pattern_CMOS_20240305.txt` | 2560 × 29 | 29 |
| CMOS (older) | `pattern_cmos.txt` | 2560 × 29 | 29 |

**Interpretation:** Row `i` = pixel column `i` in the wavelength direction. Column `j` = diffraction order `j`. Value = vertical pixel position of order center at that column.

**How loaded:** `Calibrations.load_pattern()` → `np.loadtxt(pth, dtype='int')` → `self.pattern`

**How used:** `Calibrations.make_mask(ordind)` creates a 2D pixel mask for extracting order `ordind`, using `self.dv` (half-width, default 8 pixels → 17 pixel band).

**Key insight:** The pattern is empirically determined, currently by semi-manual procedures in the notebooks (e.g., 2.1, 2.2, 4.2). The pattern is **mostly stable** — it traces the physical optics and does not significantly stretch between sessions.

### 4.2 Wavelength Calibration

**What it is:** A lookup table of identified lamp spectral lines, with their pixel positions and known wavelengths, per order.

**File format:**

```
#order  from    to      center      wavelength  band
0       0632    0657    0644.5921   794.81764   ArI
0       0670    0695    0681.7851   794.31814   NeI
```

- `order`: Diffraction order index
- `from`, `to`: Pixel range where the line appears
- `center`: Sub-pixel center (from peak detection)
- `wavelength`: Known wavelength in nm
- `band`: Species identifier (ArI, NeI, ThI, HgI, XeI, H2, etc.)

**Calibration files:**

| File | Detector | Date | Lamps |
|------|----------|------|-------|
| `Th_wavelength.txt` | CCD | Original | ThAr, Ne |
| `Th_wavelength_CMOS.txt` | CMOS | 2019-05-29 | ThAr, Ne, Hg, Xe |
| `Th_wavelength_CMOS_20240305.txt` | CMOS | 2024-03-05 | ThAr, Ne, Hg, H2 |
| `Th_wavelength_CMOS_only_Th-Ar.txt` | CMOS | ? | ThAr only subset |
| `wavelength_cmos_2024.txt` | CMOS | 2024-03-14 | ThAr, Ne, Hg, Xe (near-duplicate of `Th_wavelength_CMOS.txt`) |
| `wavelength_fujii.txt` (resources/) | Fujii | 2022-09-29 | H2 Fulcher lines only (partial, 1 order) |
| `wavelength_fujii.txt` (calibration_files/) | Fujii | ? | H2 Fulcher d3Pu bands (multiple orders) |

**How used:** `Calibrations.wavelength_calibration()`:
1. Reads the table, groups lines by order
2. Fits polynomial (degree 1 if ≤2 lines, degree 2 if ≥3 lines) per order: `pixel_center → wavelength`
3. Evaluates polynomial over full pixel range → `self.order_wavel[order, pixel]`
4. Handles partial orders (28th CMOS order extends beyond sensor edge) by filling with `np.nan`

### 4.3 Order Border Calculation

**What it is:** Determines where to cut and stitch adjacent orders to form a continuous spectrum.

**How it works:** `Calibrations.calculate_order_borders()`:
1. Uses the integrating sphere spectrum (flat illumination)
2. For each pair of adjacent orders, finds the wavelength where their signals cross
3. Creates boolean masks defining which pixels belong to which order
4. Result: `self.order_borders` — boolean array, `self.wavelength` — continuous wavelength vector

### 4.4 Absolute Calibration

**What it is:** Converts detector counts to physical radiometric units.

**Files:**
- `integrating_sphere.txt`: Known spectral radiance of the integrating sphere (µm vs mW/cm²·sr·µm)
- `sphere_*.sif` / `absolute_*.sif`: Detector images of the integrating sphere
- `sphere_*_bkg.sif` / `absolute_*_bkg.sif`: Background images

**How it works:** `Calibrations.absolute_calibration()`:
1. Interpolates known sphere irradiance at measured wavelengths
2. Extracts sphere spectrum through same pipeline (order extraction → stitching)
3. Computes sensitivity: `sphere_irradiance / (detector_counts / exposure_time)`
4. Produces three conversion factors:
   - `wmsr`: W/(m² sr nm)
   - `wm`: W/(m² nm) — assumes 4π steradians
   - `phmsr`: photons/(m² sr) — via Planck relation

### 4.5 Spectrometer Configurations

| Name | Detector | Image Format | Camera | Pattern | Wavelength Range |
|------|----------|-------------|--------|---------|-----------------|
| `"black"` | Andor CCD (1024×1024) | `.sif` via `sif_parser` | Andor | 28 orders | ~410–795 nm |
| `"black"` | Andor CMOS (2560×2160) | `.sif` via `sif_parser` | Andor | 29 orders | ~410–795 nm |
| `"fujii"` | Hamamatsu Orca QUEST | `.tiff` via PIL | Hamamatsu | (separate) | Fulcher band region |

The `"black"` spectrometer has been upgraded from CCD to CMOS; both calibrations are loaded simultaneously by the GUI and auto-selected based on image dimensions.

---

## 5. Workflow Map

### Workflow 1: Full Calibration (one-time per detector config)

```
1. Take integrating sphere image + background    → sphere.sif, sphere_bkg.sif
2. Take ThAr/Ne/Hg/Xe lamp image                → ThAr.sif
3. Determine order pattern                        → pattern.txt
     (notebooks 1.x, 2.x, 4.2)
     Method: peak detection along vertical slices of sphere/lamp image
4. Identify lamp lines per order                  → Th_wavelength_*.txt
     (notebooks 3.3, 4.3)
     Method: manual identification of known lines at pixel positions
5. Record sphere spectral data                    → integrating_sphere.txt
     (from manufacturer specification)
6. Bundle into calibration_files/
```

### Workflow 2: GUI Data Processing (routine)

```
1. Launch GUI                                     → loads both CCD and CMOS calibrations
2. Open SIF image (single file or shot number)
3. Auto-detect camera from image dimensions
4. Extract order spectra using calibration masks
5. Stitch into continuous spectrum
6. Apply absolute calibration
7. Display image + spectrum
8. Optionally: fit emission bands (CIV, He, H)
9. Optionally: save spectra and/or intensities
```

### Workflow 3: Batch Processing (routine)

```
1. Set shot number range in GUI
2. Click "Start" → loops through all shots
3. For each shot: load → extract → calibrate → fit → save
4. Progress bar tracks completion
```

### Workflow 4: Notebook-Based Analysis (exploratory)

```
1. Create Calibrations instance with file paths
2. Call clbr.start() → loads everything
3. Create EchelleImage with clbr
4. Call em.calibrate() → order spectra → full spectrum
5. Create Spectrum(em) → physical units
6. Analyze emission bands with emissionbands module
```

---

## 6. GUI Analysis

The GUI (`echelle_spectra_gui.py` + `window_layout.py`) is a single-window PyQt5 application with 5 tabs:

| Tab | Purpose |
|-----|---------|
| **Image** | 2D detector image viewer + spectrum plots + intensity time traces |
| **CIV** | Carbon ion emission band spectra (8 bands in 4×2 grid) |
| **He** | Helium emission band spectra (8 bands in 4×2 grid) |
| **H2** | Hydrogen Balmer + CH band spectra (5 bands) |
| **CMD** | Text output showing emission band reports |

### GUI → Data Coupling

The GUI directly creates `Calibrations`, `EchelleImage`, `Spectrum`, and `EmissionBand` objects. There is **no intermediate API layer** — the GUI is the primary orchestrator.

### Threading Model

4 QThread subclasses handle background work:
- `LoadCalibrationsThread` — calls `clbr.start()`
- `LoadImageThread` — loads image, runs extraction pipeline
- `FitLinesThread` — runs Gaussian fits across all frames/bands
- `SaveSpectraThread` — saves to disk

---

## 7. Example Classification

### Essential Workflow Examples (preserve and document)

| Notebook | Topic | Status |
|----------|-------|--------|
| `1.1_Cutting_Pattern_Synthetic.ipynb` | Pattern generation from synthetic data | **Essential** — demonstrates pattern creation concept |
| `1.2_Cutting_Pattern_Masks.ipynb` | Visualizing cutting masks | **Essential** — shows mask validation |
| `2.1_Cutting_Pattern_FE.ipynb` | Pattern for real Fe-lamp data | **Essential** — real-world pattern workflow |
| `3.3_Wavelength_Calibration.ipynb` | Wavelength calibration procedure | **Essential** — core calibration workflow |
| `4.2_Pattern_LHD.ipynb` | Pattern creation for LHD CMOS | **Essential** — latest detector pattern |
| `4.3_Wavelength_Calibration_LHD.ipynb` | Wavelength cal for LHD CMOS | **Essential** — latest detector wavelength cal |
| `5.0_Bands.ipynb` | Emission band analysis | **Essential** — demonstrates band processing |

### Exploratory / Experimental (preserve but flag as experimental)

| Notebook | Topic | Status |
|----------|-------|--------|
| `2.2_Cutting_Pattern_NihamaPlus.ipynb` | Pattern for Nihama+ spectrometer | Exploratory — different spectrometer |
| `3.1_loading_images.ipynb` | Basic image loading | Exploratory — superseded by GUI |
| `3.2_try_load_new_echelle.ipynb` | Testing new echelle loading | Exploratory |
| `3.5_Load_Black_Echelle.ipynb` | Loading black echelle images | Exploratory |

### Likely Obsolete / Duplicate (candidates for cleanup)

| Notebook | Issue |
|----------|-------|
| `3.4_Wavelength_echelle_spectra copy.ipynb` | **"copy"** in filename — likely accidental duplicate |
| `4.1_Load_bkack_echelle copy.ipynb` | **Typo + "copy"** in filename — likely obsolete |
| `testtool-CCD.ipynb` | Test tool — may be obsolete with CMOS upgrade |
| `testtool-CMOS.ipynb` | Test tool — still useful? |

### Data Files in Examples (review)

| File | Size | Status |
|------|------|--------|
| `pattern_image_synthetic.txt` | **100 MB** | Should be generated, not stored; add to `.gitignore` or LFS |
| `pattern.txt` | 94 KB | Duplicate of `calibration_files/pattern.txt` |
| `pattern_CMOS_20240305.txt` | 333 KB | Duplicate of `calibration_files/pattern_CMOS_20240305.txt` |
| `pattern_fujii.txt` / `_rot.txt` | 134 KB each | Fujii spectrometer patterns, not in calibration_files |

---

## 8. New Spectrometer Initialization

**Current process** (reconstructed from code + notebooks):

1. **Physical setup**: Mount spectrometer, connect camera, install lamps
2. **Take calibration images**:
   - Integrating sphere image + background (matching exposure)
   - Lamp image(s): ThAr, Ne, optionally Hg/Xe/H2
3. **Determine order pattern** (notebook workflow):
   - Load sphere or lamp image
   - Detect vertical intensity peaks for each pixel column
   - Extract and fit peak positions → order centers
   - Save as `pattern_NEWSPEC.txt`
4. **Wavelength calibration** (notebook workflow):
   - Load lamp image
   - For each order, identify known spectral lines
   - Record pixel position ↔ wavelength pairs
   - Save as `Th_wavelength_NEWSPEC.txt`
5. **Register calibration files**:
   - Place files in `resources/calibration_files/`
   - **Hardcode filename dict** in `echelle_spectra_gui.py:prepare_calibration()` or construct one manually in notebook
6. **Optionally**: Hardcode new `spec="newname"` branch in `read_image()`

**Pain points:**
- No formal "add a spectrometer" procedure
- Filename selection is hardcoded in the GUI
- `read_image()` requires code changes for new camera formats
- No validation that calibration files are consistent

---

## 9. Where Manual Recalibration Happens

### Pattern Recalibration
- **Notebooks 2.1, 2.2, 4.2**: Peak detection along vertical slices
- **Fully manual**: User must run notebook, inspect results, iterate
- **Output**: New `pattern_*.txt` file

### Wavelength Recalibration
- **Notebooks 3.3, 4.3**: Line identification in lamp spectra
- **Fully manual**: User examines each order, finds known lines, records pixel positions
- **Output**: New `Th_wavelength_*.txt` file
- **Code path**: `Calibrations.wavelength_calibration()` reads this file

### The Manual Bottleneck
The wavelength calibration file (`Th_wavelength_*.txt`) is the critical manually-created artifact:
- Each entry is a manually identified line: `order from to center wavelength species`
- `center` is determined by peak-finding (semi-automated in notebooks)
- `wavelength` is looked up from reference databases (NIST)
- The matching of peak-to-species is done by the human operator

Evidence of recalibration sessions:
- `Th_wavelength.txt` → original CCD
- `Th_wavelength_CMOS.txt` → first CMOS (2019-05-29)
- `Th_wavelength_CMOS_20240305.txt` → latest CMOS (2024-03-05)
- Comments in files like `# OK`, `# ok`, `# ?`, `# not clear` show manual inspection

---

## 10. Session Alignment Proposal

### Conceptual Framework

The order pattern is physically stable (optics don't move). What drifts between sessions is:
- Small global pixel shifts (thermal expansion, camera remounting)
- Small wavelength shifts (same physical cause)

**Key insight**: A few known strong lines (e.g., Balmer series from an H lamp) can constrain a global shift model, avoiding full recalibration.

### Where Alignment Logic Should Live

```
tools/echelle.py
  └── Calibrations
       ├── .start()           ← existing: full calibration from files
       ├── .start_cut()       ← existing: pattern-only calibration
       │
       ├── .align_session()   ← NEW: session alignment entry point
       │     1. Load a session lamp image (H lamp or similar)
       │     2. Extract order spectra using existing pattern
       │     3. Find known lines (Balmer: Hα, Hβ, Hγ, Hδ)
       │     4. Measure pixel positions of found lines
       │     5. Compare to expected positions from current wavelength table
       │     6. Solve for global shift (dx pixels, or Δλ per order)
       │     7. Update self.order_wavel with corrections
       │     8. Recompute order_borders and wavelength
       │
       └── .save_alignment()  ← NEW: persist alignment corrections
```

### Existing Abstractions to Reuse

1. **`Calibrations.wavelength_calibration()`** — already does polynomial fitting per order; extend it to accept correction offsets
2. **`emissionbands.banddata()`** — already extracts spectral windows around known lines; reuse for alignment line detection
3. **`EmissionBand.fitb()`** — already fits Gaussians to known lines; reuse to measure center positions
4. **`Calibrations.calculate_order_borders()`** — needs to be re-run after alignment

### Proposed Minimal Implementation

```python
# In Calibrations class:

def align_session(self, lamp_image_path, reference_lines=None, spec="black"):
    """Align wavelength calibration using known spectral lines.
    
    Parameters
    ----------
    lamp_image_path : str
        Path to a lamp image (e.g., H lamp with Balmer lines)
    reference_lines : list of dict
        Known lines: [{"order": 23, "wavelength": 656.279, "name": "H-alpha"}, ...]
    """
    # 1. Load lamp image using existing pipeline
    lamp = EchelleImage(lamp_image_path, clbr=self, spec=spec)
    lamp.calculate_order_spectra()
    lamp.correct_order_shapes()
    
    # 2. For each reference line, find its actual pixel position
    shifts = []
    for ref in reference_lines:
        order = ref["order"]
        expected_wl = ref["wavelength"]
        # Find expected pixel from current calibration
        expected_pixel = np.interp(expected_wl, self.order_wavel[order], np.arange(self.DIMW))
        # Find actual peak near expected position
        actual_pixel = find_peak_near(lamp.order_spectra[0][order], int(expected_pixel))
        shifts.append(actual_pixel - expected_pixel)
    
    # 3. Solve global shift (median, or per-order polynomial)
    global_shift = np.median(shifts)
    
    # 4. Apply shift to wavelength solution
    x_shifted = np.arange(self.DIMW) - global_shift
    for o in range(self.order_wavel.shape[0]):
        self.order_wavel[o] = np.interp(np.arange(self.DIMW), x_shifted, self.order_wavel[o])
    
    # 5. Recompute downstream
    self.calculate_order_borders()
    self.absolute_calibration()
```

### Reference Lines for H Lamp Alignment

Using Balmer lines already defined in `emissiondata.py`:

| Line | Wavelength (nm) | Approx. Order (CMOS) | Strength |
|------|-----------------|----------------------|----------|
| Hα | 656.279 | ~6–7 | Very strong |
| Hβ | 486.135 | ~17–18 | Strong |
| Hγ | 434.047 | ~24–25 | Moderate |
| Hδ | 410.173 | ~27 | Weaker |

These span the full wavelength/order range and provide excellent constraints for a global shift or low-order polynomial correction.

### Integration Strategy (Minimal Disruption)

1. **Phase 1**: Add `align_session()` method to `Calibrations` in `echelle.py`. No GUI changes needed — usable from notebooks first.
2. **Phase 2**: Add "Align" button to GUI that takes a lamp image path and runs alignment.
3. **Phase 3**: Save alignment corrections to a session file alongside the base calibration.

---

## 11. God Objects and Oversized Files

### `Calibrations` class (~390 lines, moderate)

**Not a god object yet**, but trending toward one. It handles:
- Pattern loading and mask generation
- Wavelength calibration (polynomial fitting)
- Order border calculation
- Absolute calibration
- Integrating sphere image management

**Recommendation**: Keep together for now. The responsibilities are naturally coupled. If/when alignment is added, consider extracting a `WavelengthCalibration` helper class.

### `EchelleSpectraGUI` class (~760 lines) — **the real god object**

This class handles:
- UI initialization and signal connection
- Calibration setup
- Emission band setup
- Shot loop management
- Image loading orchestration
- Spectrum display
- Emission band display (3 separate tabs)
- Intensity display
- Fit orchestration
- Data saving
- Drag-and-drop

**Recommendation**: Consider extracting:
- A `ShotProcessor` class that handles the load→extract→fit→save pipeline
- Band display logic could be templated (C/He/H tabs are near-identical)

### `window_layout.py` (546 lines, `Ui_MainWindow`) — **generated-style but hand-written**

Huge procedural `setupUi()` method. This is typical for Qt UI code, but:
- Could be partially replaced by `.ui` file from Qt Designer
- Band plot creation (`create_plots()`) already shows good factoring

### `echelle.py` (900 lines) — **well-structured but monolithic**

Contains 3 major classes + 1 function + header template. Natural file; the concern is **growth** as features are added.

**Recommendation**: When adding alignment logic, consider splitting into:
- `echelle_image.py` — `read_image()` + `EchelleImage`
- `calibrations.py` — `Calibrations`
- `spectrum.py` — `Spectrum`
- Keep the old `echelle.py` as a re-export shim for backwards compatibility

---

## 12. Duplicated and Dead Code

### Duplicated

| What | Where | Issue |
|------|-------|-------|
| Header template | `echelle.py:878` AND `resources/header_template.txt` | Two different templates! The one in `echelle.py` is simpler (used by `Spectrum.save()`); the one in `header_template.txt` is more detailed (used by GUI's `save_intensities()`) |
| Calibration filenames | `echelle.py:Calibrations.filenames` default AND `gui.py:prepare_calibration()` | GUI always overrides the defaults; defaults point to CCD files |
| Pattern files | `examples/pattern.txt` AND `calibration_files/pattern.txt` | 94 KB duplicate |
| Pattern files | `examples/pattern_CMOS_20240305.txt` AND `calibration_files/pattern_CMOS_20240305.txt` | 333 KB duplicate |
| Wavelength files | `Th_wavelength_CMOS.txt` ≈ `wavelength_cmos_2024.txt` | Near-identical content, different names |
| `wavelength_fujii.txt` | `resources/wavelength_fujii.txt` (1 line) AND `calibration_files/wavelength_fujii.txt` (full) | Different content! Root-level one has only 1 order |

### Dead / Obsolete Code

| What | Where | Issue |
|------|-------|-------|
| `DIMW = 1024` / `DIMO = 1024` | `echelle.py:13–14` | Module-level constants, **never used** — dimensions come from sphere image |
| Commented-out `read_image` method | `echelle.py:120–126` | Dead code with TODO: Remove |
| `show_masks` method | `echelle.py:544–547` | Commented out entirely |
| `print_wavelength` method | `echelle.py:621–636` | Commented out entirely |
| `__init__.py` sys.path hack | `__init__.py:5–6` | Appends package dir to sys.path — fragile, may cause import shadowing |
| Old imports in `__init__.py` | Lines 9–13 | Uses bare `import echelle_spectra_gui as gui` relying on sys.path hack |

### Import System Issues

The `__init__.py` uses a `sys.path` hack to enable bare imports like `import tools.echelle`. This works but is fragile and non-standard. The GUI uses the same pattern:

```python
import tools.echelle as ech          # not: from echelle_spectra.tools import echelle
from __init__ import __version__     # not: from echelle_spectra import __version__
```

This means the package cannot be imported normally (e.g., `from echelle_spectra.tools import echelle` would fail from outside the package directory).

---

## 13. Proposed MkDocs Structure

```yaml
# mkdocs.yml
site_name: Echelle Spectra
theme:
  name: material
  palette:
    scheme: slate

nav:
  - Home: index.md
  - Getting Started:
    - Installation: getting-started/installation.md
    - Quick Start (GUI): getting-started/gui-quickstart.md
    - Quick Start (Python): getting-started/python-quickstart.md

  - Concepts:
    - Echelle Spectrometer Overview: concepts/echelle-overview.md
    - Order Pattern: concepts/order-pattern.md
    - Wavelength Calibration: concepts/wavelength-calibration.md
    - Absolute Calibration: concepts/absolute-calibration.md
    - Emission Bands: concepts/emission-bands.md
    - Spectrometer Configurations: concepts/spectrometers.md

  - Workflows:
    - Full Calibration Workflow: workflows/full-calibration.md
    - Adding a New Spectrometer: workflows/new-spectrometer.md
    - Session Alignment: workflows/session-alignment.md
    - Batch Processing: workflows/batch-processing.md
    - Emission Band Fitting: workflows/band-fitting.md
    - Saving and Exporting Data: workflows/data-export.md

  - Calibration Reference:
    - Calibration Files: reference/calibration-files.md
    - Pattern File Format: reference/pattern-format.md
    - Wavelength Table Format: reference/wavelength-format.md
    - Integrating Sphere Data: reference/sphere-data.md
    - Detector Configurations: reference/detectors.md

  - API Reference:
    - echelle (core): api/echelle.md
    - emissionbands: api/emissionbands.md
    - emissiondata: api/emissiondata.md
    - GUI: api/gui.md

  - Examples:
    - Pattern Creation: examples/pattern-creation.md
    - Wavelength Calibration: examples/wavelength-calibration.md
    - Image Processing: examples/image-processing.md

  - Development:
    - Architecture: dev/architecture.md
    - Contributing: dev/contributing.md
    - Changelog: dev/changelog.md

plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          paths: [src]
```

### Key Documentation Priorities (ordered)

1. **Workflow: Full Calibration** — The most critical undocumented knowledge. Currently lives only in notebooks and the author's memory.
2. **Concepts: Order Pattern** — What a pattern file means, how it's created, why it's stable.
3. **Concepts: Wavelength Calibration** — The manual line identification process, the polynomial model, the file format.
4. **Workflow: Adding a New Spectrometer** — Step-by-step guide from hardware to working calibration.
5. **Workflow: Session Alignment** — Future workflow, document the concept now.
6. **Calibration Reference** — Formal description of every file in `calibration_files/`.

---

## 14. Safe Incremental Improvements

### Priority 1: Low-Risk, High-Value

| # | Improvement | Risk | Impact |
|---|-------------|------|--------|
| 1 | **Fix imports**: Replace `sys.path` hack with proper relative imports | Low | Fixes packaging; allows normal `from echelle_spectra.tools import echelle` |
| 2 | **Remove dead code**: Delete commented-out methods, unused constants `DIMW`/`DIMO` | None | Cleaner codebase |
| 3 | **Consolidate duplicate header templates** | Low | Single source of truth |
| 4 | **Add `.gitignore` entry for `pattern_image_synthetic.txt`** | None | -100 MB from repo |
| 5 | **Remove duplicate pattern files from `examples/`** | Low | Use symlinks or document the path to originals |
| 6 | **Rename misspelled notebooks** | None | `4.1_Load_bkack_echelle copy.ipynb` → archive or rename |
| 7 | **Clarify duplicate wavelength files** | Low | Delete `wavelength_cmos_2024.txt` if truly a duplicate |

### Priority 2: Moderate Effort, High Value

| # | Improvement | Risk | Impact |
|---|-------------|------|--------|
| 8 | **Extract spectrometer config into data files** | Medium | Replace hardcoded filename dicts with TOML/JSON per spectrometer |
| 9 | **Add `read_image()` dispatcher registry** | Low | Make it extensible without code changes |
| 10 | **Create `mkdocs.yml` skeleton** | Low | Start docs migration |
| 11 | **Write calibration workflow documentation** | None | Capture implicit knowledge |
| 12 | **Add type hints to core classes** | Low | Better IDE support and documentation |

### Priority 3: Larger but Safe

| # | Improvement | Risk | Impact |
|---|-------------|------|--------|
| 13 | **Split `echelle.py`** into `echelle_image.py`, `calibrations.py`, `spectrum.py` | Medium | Cleaner imports, easier maintenance; keep re-export shim |
| 14 | **Implement `align_session()`** | Medium | Key workflow improvement |
| 15 | **Extract shot processing from GUI** into reusable class | Medium | Enables headless batch processing |

---

## 15. Minimal-Risk Modernization Plan

### Phase 1: Cleanup (no behavior changes)

1. Fix `__init__.py` imports to use proper relative imports
2. Fix `echelle_spectra_gui.py` and `__main__.py` imports
3. Remove dead/commented-out code in `echelle.py`
4. Remove duplicate files (patterns in examples, near-duplicate wavelength files)
5. Add `pattern_image_synthetic.txt` to `.gitignore`
6. Clean up example naming (remove "copy" duplicates)

### Phase 2: Documentation (no code changes)

1. Create `mkdocs.yml` with proposed structure
2. Write "Order Pattern" concept page
3. Write "Wavelength Calibration" concept page
4. Write "Full Calibration Workflow" page (from notebooks)
5. Write "Calibration Files Reference" page
6. Add docstrings to undocumented methods

### Phase 3: Spectrometer Config Externalization

1. Create `spectrometers/` directory with TOML files per config:
   ```toml
   # spectrometers/black_ccd.toml
   name = "CCD"
   camera = "Andor CCD"
   image_format = "sif"
   pattern = "pattern.txt"
   wavelength = "Th_wavelength.txt"
   sphere = "absolute_20170613_b8_0.2_v2.sif"
   background = "absolute_20170613_b8_0.2_bkg.sif"
   integral = "integrating_sphere.txt"
   dv = 8
   ```
2. Update `Calibrations.__init__()` to optionally load from TOML
3. Update GUI to discover spectrometer configs

### Phase 4: Session Alignment

1. Add `align_session()` to `Calibrations`
2. Test in notebooks with H lamp images
3. Add GUI button for alignment
4. Document the alignment workflow

### Phase 5: Optional Improvements

1. Split `echelle.py` into separate files
2. Extract shot processing from GUI
3. Add proper logging (replace `print()` statements)
4. Add basic tests for calibration pipeline
5. Consider `pandas` for wavelength tables (currently raw `np.loadtxt`)

---

## Appendix: Key Code Entry Points

| Use Case | Entry Point |
|----------|-------------|
| Launch GUI | `python -m echelle_spectra` or `echelle_spectra` CLI |
| Programmatic calibration | `from echelle_spectra.tools import echelle; clbr = echelle.Calibrations(folder, filenames); clbr.start()` |
| Load image | `em = echelle.EchelleImage(path, clbr=clbr); em.calibrate()` |
| Get spectrum | `s = echelle.Spectrum(em)` |
| Emission band fit | `from echelle_spectra.tools import emissiondata as ebd; ebd.halpha.fitb(frame, spectrum)` |
| Config defaults | `resources/defaults.toml` |
| Add spectrometer | Manually create pattern + wavelength files, hardcode filenames in GUI |
