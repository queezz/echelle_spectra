# Calibration Alignment

This page describes the headless v1 workflow for correcting an existing wavelength
lookup table after a small detector or optics shift.  It is intended for notebooks and
future CLI work, not the GUI.

---

## Goal

The baseline CMOS calibration table (`Th_wavelength_CMOS_20240305.txt`) is preserved as
historical truth.  A session alignment measures where strong lamp lines actually appear
in a new lamp frame, fits a rigid detector correction, and writes a separate settings
file or adjusted lookup table.

The v1 correction model is deliberately small:

- translation in detector x/y
- small rotation
- no stretch, shear, or warp

---

## Candidate lines

Use the curated local wavelength table first.  The default helper selection keeps
active `NeI` rows marked `ok`/`OK` in the inline comments and rejects very narrow or
very broad manual pixel intervals.

```python
from pathlib import Path

from echelle_spectra.tools.calibration_alignment import (
    load_wavelength_table,
    select_candidate_lines,
)

cal_dir = Path("src/echelle_spectra/resources/calibration_files")
table_path = cal_dir / "Th_wavelength_CMOS_20240305.txt"

rows = load_wavelength_table(table_path)
candidates = select_candidate_lines(rows)
```

NIST matching is intentionally not part of v1.  The local curated table is reproducible
and keeps the manual comments that already flag doubtful lines.

For calibration review, cached NIST ASD line exports can be overlaid on a lamp
order spectrum with `echelle-nist-overlay`. This is a diagnostic tool for
finding candidate anchors; it does not download from NIST or change resource
tables by itself.

Example using cached Th I and Ar I CSV exports:

```bash
echelle-nist-overlay \
  src/echelle_spectra/resources/calibration_files/ThAr-0.3s-x3_20240305.sif \
  --line-list ThI=local/thar_nist_synthetic_overlay/nist_th_i_600_640.csv \
  --line-list ArI=local/thar_nist_synthetic_overlay/nist_ar_i_600_640.csv \
  --orders 8 \
  --min-nm 608 \
  --max-nm 634.5 \
  --output-dir local/thar_nist_synthetic_overlay/order8-review \
  --candidate-table-out local/thar_nist_synthetic_overlay/order8-review/Th_wavelength_CMOS_20240305_plus_order8_nist_thar_candidates.txt
```

Common lamp presets can also be resolved from cached NIST ASD exports. Presets
include `thar`, `hg`, `ne`, `he`, `ar`, `th`, `h`, and `h2`; ion stages I and II
are included where useful. The package bundles ThAr, Hg, and Ne caches covering
380-810 nm — the whole range the instrument reaches, with margin — so those
presets can be used without an explicit cache directory anywhere on the
detector. Regenerate them with `python -m echelle_spectra.tools.nist_cache_refresh`, which
carries the exact NIST ASD query:

```bash
echelle-nist-overlay \
  src/echelle_spectra/resources/calibration_files/ThAr-0.3s-x3_20240305.sif \
  --lamp thar \
  --orders 8-10 \
  --min-nm 578 \
  --max-nm 635 \
  --output-dir local/runs/2026-06-11_thar-review
```

Use `--line-list-dir` for lamps or wavelength ranges not present in the bundled
cache. The CLI still requires cached NIST ASD exports; it does not query NIST
during a calibration run — refreshing the cache is a deliberate act performed by
`echelle_spectra.tools.nist_cache_refresh`, not a side effect of calibrating. Use `echelle-nist-overlay --list-lamps` to inspect the
built-in presets.

The overlay writes:

- order plots with measured lamp spectra, line-list sticks, and broadened synthetic spectra;
- measured peak and nearest-line CSVs;
- a candidate-anchor CSV;
- optionally, a generated candidate wavelength table.

Review the overlay plot before promoting candidate rows. Dense Th I regions can
have many nearby NIST lines, so the candidate filter requires a local dominant
line, a successful Gaussian centroid, and a wavelength match within the chosen
tolerance. Treat the synthetic trace as line-position and local-structure
evidence, not as a calibrated lamp-intensity model; NIST tabular relative
intensities may not match the measured lamp, detector, and extraction response.

For Fulcher-alpha order-8 review, inspect neighboring overlapping orders as well.
Some wavelengths near the edge of order 8 are better centered in order 7, so
order 7 can verify that a local lamp structure is real even when the order-8
polynomial is poorly constrained.

`echelle-wavelength-qc` writes `order_residual_and_dispersion_qc.png`. The upper
panel shows per-order wavelength-table RMS. The lower panel shows mid-order
dispersion smoothness, defined as the derivative of each fitted wavelength
polynomial at detector center, in nm/pixel. A kink in this curve can indicate a
physically implausible order solution even when an individual order has acceptable
statistics.

---

## Fit centroids

After loading a lamp frame with the existing `Calibrations` and `EchelleImage` path,
rank candidate windows, check saturation on raw 2D detector pixels, and then fit
single-Gaussian centroids in the extracted order spectra:

```python
from echelle_spectra.tools.calibration_alignment import (
    measure_detector_window_saturation,
    measure_line_centroids,
    rank_candidate_lines,
)

ranked = rank_candidate_lines(
    lamp_image.order_spectra[0] - lamp_background.order_spectra[0],
    candidates,
    min_snr=5.0,
)
saturation = measure_detector_window_saturation(
    lamp_image.images,
    calibration.pattern,
    candidates,
    saturation_level=0.98 * 65535,
)

fits = measure_line_centroids(
    lamp_image.order_spectra[0] - lamp_background.order_spectra[0],
    [stat.line for stat in ranked],
    min_snr=5.0,
)

good = [fit for fit in fits if fit.success]
```

Rejected rows keep a short `reason`, for example `saturated`, `low snr`, or
`sigma out of range`.

---

## Fit the rigid correction

Use the order pattern to convert lookup rows and measured centroids into detector
points, then fit translation plus rotation:

```python
import numpy as np

from echelle_spectra.tools.calibration_alignment import (
    AlignmentSettings,
    detector_points_from_lines,
    fit_rigid_transform,
    save_alignment_settings,
)

lines = [fit.line for fit in good]
measured_centers = [fit.center_pixel for fit in good]

expected_xy = detector_points_from_lines(lines, calibration.pattern)
measured_xy = detector_points_from_lines(lines, calibration.pattern, measured_centers)

transform, rms_px = fit_rigid_transform(expected_xy, measured_xy)

settings = AlignmentSettings(
    instrument_id="lhd_cmos",
    created_at="2026-06-04",
    alignment_dataset_id="20250926",
    alignment_source_dir="local/20250926_calib",
    alignment_lamp="Ne",
    signal_file="Ne-0.02s-x3-bright-lines.sif",
    background_file="Ne-0.02s-x3-bright-lines_bg.sif",
    base_wavelength_file="Th_wavelength_CMOS_20240305.txt",
    base_pattern_file="pattern_CMOS_20240305.txt",
    sphere_file="sphere-0.1s-x3.sif",
    sphere_background_file="sphere-0.1s-x3-bg.sif",
    output_wavelength_file="Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
    transform=transform,
    n_lines=len(lines),
    rms_px=rms_px,
    notes="Ne I rigid alignment",
)
save_alignment_settings(
    settings,
    cal_dir / "alignments" / "lhd_cmos_alignment_20250926.settings.toml",
)
```

---

## Write an adjusted lookup table

The original lookup table should not be overwritten.  Generate a new file when the
correction should be consumed by `Calibrations` or by `echelle_optics`:

```python
from echelle_spectra.tools.calibration_alignment import (
    apply_rigid_correction_to_lines,
    write_wavelength_table,
)

adjusted = apply_rigid_correction_to_lines(rows, calibration.pattern, transform)
write_wavelength_table(
    adjusted,
    cal_dir / "alignments" / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
    metadata=[
        ("Generated", "2026-06-04"),
        ("Base wavelength file", "Th_wavelength_CMOS_20240305.txt"),
        ("Base pattern file", "pattern_CMOS_20240305.txt"),
        ("Alignment dataset", "20250926"),
        ("Alignment source dir", "local/20250926_calib"),
        ("Signal", "Ne-0.02s-x3-bright-lines.sif"),
        ("Background", "Ne-0.02s-x3-bright-lines_bg.sif"),
        ("Correction model", "rigid detector transform, dx/dy/theta"),
    ],
)
```

Point `Calibrations.filenames["wavelength"]` at the adjusted table for a session that
needs the correction.

---

## CLI Workflow

Use `echelle-align` when the inputs and thresholds are known and the task is to
reproduce the alignment without Jupyter. Without `--save`, the command prints
diagnostics and does not write artifacts:

```bash
echelle-align \
  local/20250926_calib/Ne-0.02s-x3-bright-lines.sif \
  local/20250926_calib/Ne-0.02s-x3-bright-lines_bg.sif \
  local/20250926_calib/sphere-0.1s-x3.sif \
  local/20250926_calib/sphere-0.1s-x3-bg.sif \
  --pattern pattern_CMOS_20250926.txt
```

To write the reviewed settings and adjusted wavelength table:

```bash
echelle-align \
  local/20250926_calib/Ne-0.02s-x3-bright-lines.sif \
  local/20250926_calib/Ne-0.02s-x3-bright-lines_bg.sif \
  local/20250926_calib/sphere-0.1s-x3.sif \
  local/20250926_calib/sphere-0.1s-x3-bg.sif \
  --pattern pattern_CMOS_20250926.txt \
  --save
```

The default output paths are:

- `src/echelle_spectra/resources/calibration_files/alignments/lhd_cmos_alignment_20250926.settings.toml`
- `src/echelle_spectra/resources/calibration_files/alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt`

The notebook remains useful for plotting candidate windows and inspecting failed
fits, but the CLI is the preferred reproducible workflow once settings are known.

---

## Validation

Recommended checks:

1. Compare lamp-line residuals before and after correction.
2. Run the [Wavelength line validation](line-validation.md) gate on LHD Balmer
   peaks (`H-alpha`, `H-beta`, `H-gamma`).
3. Check identified Fulcher lines as supporting evidence when blends are
   understood.
4. Report systematic offsets in nm, especially near the previously suspected
   approximately 0.1 nm error.
