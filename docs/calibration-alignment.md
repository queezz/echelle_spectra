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

---

## Fit centroids

After loading a lamp frame with the existing `Calibrations` and `EchelleImage` path,
fit single-Gaussian centroids in the extracted order spectra:

```python
from echelle_spectra.tools.calibration_alignment import measure_line_centroids

fits = measure_line_centroids(
    lamp_image.order_spectra[0],
    candidates,
    saturation_level=60000,
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
    base_wavelength_file="Th_wavelength_CMOS_20240305.txt",
    transform=transform,
    n_lines=len(lines),
    rms_px=rms_px,
    notes="Ne I rigid alignment",
)
save_alignment_settings(settings, cal_dir / "lhd_cmos_20240305.settings.toml")
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
write_wavelength_table(adjusted, cal_dir / "Th_wavelength_CMOS_20240305_aligned.txt")
```

Point `Calibrations.filenames["wavelength"]` at the adjusted table for a session that
needs the correction.

---

## Validation

Recommended checks:

1. Compare lamp-line residuals before and after correction.
2. Check LHD Balmer peaks (`H-alpha`, `H-beta`, `H-gamma`) as rough validation.
3. Check identified Fulcher lines for a stricter wavelength residual estimate.
4. Report systematic offsets in nm, especially near the previously suspected
   approximately 0.1 nm error.
