# Pattern Extraction

The CMOS order-pattern notebook now has a headless implementation in
`echelle_spectra.tools.pattern_extraction`.  The helpers keep the same working
model as `examples/workflows/black_cmos/02_pattern_calibration.ipynb`:

1. Average sphere and background frames.
2. Subtract the background.
3. Sample several detector columns.
4. Detect vertical order peaks in each sampled column.
5. Fit a smooth polynomial trace for each order.
6. Write the resulting `(columns, orders)` pattern table after visual review.

```python
import numpy as np

from echelle_spectra.tools.echelle import Calibrations
from echelle_spectra.tools.pattern_extraction import (
    PatternExtractionConfig,
    extract_order_pattern,
    extract_order_pattern_near_prior,
    subtract_background,
    trial_order_pattern_extraction,
)

cb = Calibrations(folder=str(CALIB_DIR), filenames=files_cmos)
cb.load_sphere()

image = subtract_background(cb.sphr.images, cb.bkgr.images)

config = PatternExtractionConfig(
    expected_orders=29,
    sample_step_px=150,
    sample_count=10,
    peak_threshold=0.13,
    peak_min_dist_px=50,
)
result = extract_order_pattern(image, config=config)

np.savetxt(CALIB_DIR / "pattern_CMOS_NEWDATE.txt", result.pattern, fmt="%d")
```

The fitted pattern should still be overlaid on the sphere frame before saving.
The automation validates peak counts and array shape, but it does not replace
the physical sanity check against the detector image.

If a fixed setting misses or adds an order in one sampled column, scan a small
grid of thresholds and column starts first:

```python
trials = trial_order_pattern_extraction(
    image,
    config=config,
    threshold_values=[0.10, 0.11, 0.12, 0.13, 0.14, 0.15],
    column_start_values=[530, 680, 750],
)
best = next(trial for trial in trials if trial.success)
result = best.result
```

Failed trials keep their peak counts, which makes it clear which sampled
columns need manual inspection.

For session-to-session recalibration, an existing pattern can be used as a
prior. This searches near each reference trace and returns one peak per order:

```python
reference_pattern = np.loadtxt("pattern_CMOS_20240305.txt", dtype=int)
result = extract_order_pattern_near_prior(
    image,
    reference_pattern,
    config=config,
    columns_px=best.columns_px,
    search_radius_px=20,
)
```

This is not a full trace-tracking algorithm, but it is useful when the new
geometry is mostly a shifted version of the previous calibration.

## Notebook Migration

The notebook can now focus on review controls and plotting:

- Load the calibration files with `Calibrations`.
- Call `subtract_background()` and `extract_order_pattern()`.
- Plot `result.detections` and `result.pattern` over the sphere image.
- Save only after all sampled columns detect the expected order count and the
  traces look physically reasonable.
