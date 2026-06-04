# Wavelength Line Validation

Use line validation as a gate between wavelength calibration/alignment and
SpectroCube export. It checks known emission-line peak positions in an
experimental SIF file against the calibrated wavelength axis.

This is not a new calibration-generation step. It does not write wavelength
tables, fit molecular populations, or analyze Fulcher intensities.

---

## Current calibration context

For the 20250926 CMOS/LHD data, validate with the accepted calibration files:

- `pattern_CMOS_20250926.txt`
- `alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt`
- `alignments/lhd_cmos_alignment_20250926.settings.toml`

The CLI defaults to the accepted pattern and aligned wavelength table.

---

## Wavelength convention

The current validation pass treats expected Balmer and Fulcher wavelengths as
air wavelengths.

The built-in Balmer targets are:

| line | air wavelength nm |
| --- | ---: |
| H-alpha | 656.279 |
| H-beta | 486.135 |
| H-gamma | 434.047 |

Do not compare these results against vacuum wavelengths unless the expected
line table has been explicitly converted first. The air-to-vacuum difference is
far larger than the current measured Balmer residuals.

---

## Validate Balmer lines

Run the Balmer-only gate first:

```bash
echelle-validate-lines local/193778_Echelle.SIF --line-set balmer
```

Equivalent direct module invocation inside the repository venv:

```bash
/Users/queezz/.venvs/echelle-spectra/bin/python -m echelle_spectra.line_validation_cli \
  local/193778_Echelle.SIF \
  --line-set balmer
```

For multiframe files, `Spectrum` keeps its existing background-frame detection.
The validation aggregate is the mean over non-background signal frames, and the
report includes how many signal frames could be fit for each line.

For `local/193778_Echelle.SIF`, frames `40..49` were detected as background and
frames `0..39` were used as signal frames.

| line | expected nm | measured nm | residual nm | order | peak SNR | notes |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| H-alpha | 656.27900 | 656.27802 | -0.00098 | 6 | 208.41 | 32/40 signal frames fit; frame sd 0.0006 nm |
| H-beta | 486.13500 | 486.13279 | -0.00221 | 19 | 318.81 | 32/40 signal frames fit; frame sd 0.0004 nm |
| H-gamma | 434.04700 | 434.04832 | +0.00132 | 24 | 87.65 | 30/40 signal frames fit; frame sd 0.0011 nm |

Balmer residual RMS was `0.001592 nm`, which is good enough to proceed to
SpectroCube export for this dataset.

---

## Fulcher checks

Fulcher-alpha line checks can be included as supporting evidence:

```bash
echelle-validate-lines \
  local/193778_Echelle.SIF \
  --line-set balmer-fulcher \
  --min-snr 10 \
  --max-abs-residual-nm 0.08
```

The Fulcher targets are loaded from the `fulcheranalyzer` H2 Q-branch wavelength
table:

```text
fulcheranalyzer/src/fulcher_analyzer/data_molecular/fulcher-α_band_wavelength.txt
```

Treat this as a visual and positional validation aid. In the LHD data, many
Fulcher features are blended at the instrumental width, and some table entries
share or overlap nearby measured peaks. Use Balmer lines as the primary gate;
use Fulcher lines to confirm that the molecular band structure lands in the
right detector/wavelength neighborhood.

Downstream Fulcher extraction, deblending, population analysis, and temperature
fits stay in `fulcheranalyzer` or later SpectroCube consumers.

---

## Python API

The core helper operates on extracted spectra and caller-supplied line tables:

```python
from echelle_spectra.tools.line_validation import (
    balmer_air_targets,
    build_stitched_order_index,
    validate_lines,
)

targets = balmer_air_targets()
order_index = build_stitched_order_index(calibration)
results = validate_lines(
    spectrum.wavelength,
    spectrum.counts,
    targets,
    signal_frames=range(40),
    order_index=order_index,
)
```

Each result contains:

- `line`
- `expected_nm`
- `measured_nm`
- `residual_nm`
- `order`
- `peak_snr`
- `notes`

The API deliberately accepts external line lists so molecule-specific knowledge
does not have to live inside `echelle_spectra`.
