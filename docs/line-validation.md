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

### 2026-06-10 order-8 Fulcher-alpha check

A local QC run compared the accepted 20250926 aligned wavelength table against a
candidate table with the order-8 `H2 629.662` calibration row removed:

```text
local/fulcher_alpha_calibration_qc/2026-06-10_22-04-fulcher-alpha-calibration-qc
```

The run selected downstream "good" Fulcher-alpha H2 Q-branch lines from the
`2026-fulcher-extractor` policy set:

- overview QC lines;
- `matrix_action=keep`;
- not Boltzmann-excluded;
- not `suspicious` or `accept_with_warning`;
- including trusted decontaminated overview lines without explicit policy rows.

On `local/193778_Echelle.SIF`, the comparison was:

| calibration | selected | fit | missed | order-8 RMS nm | all-line RMS nm | max abs nm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `aligned_20250926` | 15 | 11 | 4 | 0.079666 | 0.054069 | 0.159501 |
| `drop_order8_h2_629662` | 15 | 10 | 5 | 0.044075 | 0.028627 | 0.087834 |

This supports the earlier lamp-line QC result: the order-8 `H2 629.662` row is
not a safe calibration anchor for the 20250926 wavelength table. Removing it
improves Fulcher-alpha reproduction, but does not fully solve the local
623-630 nm problem. In particular, `H2_Q9_1-1` at `623.7457 nm` remains a large
positive residual (`+0.087834 nm`) in the candidate table, and `H2_Q7_2-2` at
`629.6622 nm` drops below the current fit criteria after the row is removed.

Treat this as evidence for a local order-aware correction pass rather than as a
final replacement calibration.

### 2026-06-10 Th-Ar/NIST order-8 candidate check

A follow-up local run used the old 2024 Th-Ar lamp with cached NIST ASD Th I and
Ar I line exports to generate independent order-8 candidate anchors:

```text
local/thar_nist_synthetic_overlay/2026-06-11_00-05-package-smoke-pixel-fix
```

The valid package-generated run fixed detector-pixel bookkeeping for decreasing
wavelength orders and produced:

- `Th_wavelength_CMOS_20240305_plus_order8_nist_thar_candidates.txt`
- `Th_wavelength_CMOS_20240305_plus_order8_nist_thar_drop_h2_629662_aligned_to_20250926.txt`
- `2026-06-11_00-03-fulcher-alpha-qc-corrected-thar-drop629/`

The most useful candidate disabled the old order-8 `H2 629.6622` row, added
NIST-backed Th/Ar rows to the 2024 table, then applied the normal 2025 Ne-lamp
rigid alignment. Order-8 wavelength-table residuals improved substantially:

| calibration | order-8 rows | order-8 RMS nm | max abs nm |
| --- | ---: | ---: | ---: |
| `aligned_20250926` | 9 | 0.069843 | 0.161373 |
| `aligned_plus_order8_thar_drop629` | 18 | 0.007298 | 0.026176 |

Fulcher-alpha reproduction also improved:

| calibration | selected | fit | missed | order-8 RMS nm | all-line RMS nm | max abs nm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `aligned_20250926` | 15 | 11 | 4 | 0.079666 | 0.054069 | 0.159501 |
| `aligned_plus_order8_thar_drop629` | 15 | 10 | 5 | 0.045576 | 0.029552 | 0.091117 |

The extra missed line is `H2_Q7_2-2` at `629.6622 nm`, which is expected when
the calibration no longer forces that endpoint feature. The mid-order Q branch
improves: for example, `H2_Q2_2-2` moves from `-0.006696 nm` to `-0.001229 nm`,
and `H2_Q3_2-2` moves from `-0.015923 nm` to `-0.002203 nm`. `H2_Q9_1-1`
remains a likely blend/identity problem rather than a clean calibration anchor.

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
