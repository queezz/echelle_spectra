# Batch SpectroCube Workflow

This workflow turns a validated folder of LHD CMOS Echelle `.SIF` files into
cropped, absolutely calibrated SpectroCube `.nc` files.

Use this after the wavelength validation gate has passed for the calibration
context. For the current 20250926 LHD CMOS setup, the accepted context is:

- `pattern_CMOS_20250926.txt`
- `alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt`
- `wavelength_medium = "air"`
- local sphere: `local/20250926_calib/sphere-0.1s-x3.sif`
- local sphere background: `local/20250926_calib/sphere-0.1s-x3-bg.sif`

## Products

The counts export is useful as an IO smoke test, but the routine science product
should be absolute spectral radiance:

| Field | Value |
| --- | --- |
| `units` | `wmsr` |
| `calibration_type` | `absolute` |
| `intensity_units` | `W/m2/nm/sr` |
| `wavelength_medium` | `air` |

The current config crops the unstable low-wavelength edge below `403.0 nm`. This
crop is a measured property of the camera/spectrometer/calibration setup, so it
lives in the config and should not be remeasured for every shot.

Some absolute-calibration columns can be non-finite before export, especially at
weak or unstable wavelength edges. `Spectrum` construction allows those known
columns without printing `0 * inf` runtime warnings, and SpectroCube export drops
non-finite wavelength/intensity columns by default. The dropped count is recorded
in metadata. Unexpected errors or failed exports should still be investigated.

## Calibration config

The reusable calibration/export config is:

```text
src/echelle_spectra/resources/calibration_files/export_configs/lhd_cmos_20250926.toml
```

Key export settings:

```toml
[export]
units = "wmsr"
output_suffix = "_spectrocube_wmsr_403nm"
wavelength_min_nm = 403.0
drop_nonfinite_columns = true
calibration_source = "local 20250926 integrating sphere (echelle_spectra)"
```

The export records crop provenance in SpectroCube metadata:

- `wavelength_crop_min_nm`
- `original_wavelength_min_nm`
- `original_wavelength_max_nm`
- `original_wavelength_points`
- `dropped_wavelength_crop_columns`
- `dropped_nonfinite_wavelength_columns`

## Single-shot plan

The first validated single-shot plan is:

```text
src/echelle_spectra/resources/spectrocube_plans/lhd_193778_wmsr.toml
```

Run it from the repository root:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli --plan src\echelle_spectra\resources\spectrocube_plans\lhd_193778_wmsr.toml
```

Expected output:

```text
local/193778_Echelle_spectrocube_wmsr_403nm.nc
```

## Dataset batch plan

The packaged dataset batch plan template is:

```text
src/echelle_spectra/resources/spectrocube_plans/lhd_20250926_wmsr_batch.toml
```

It uses portable `local/...` paths so package resources do not contain
machine-specific absolute paths. For a real lab dataset, either copy the plan to
an ignored local path and edit the dataset locations there, or keep the packaged
plan and override the input/output paths from the CLI.

Do not commit personal absolute dataset paths into `src/echelle_spectra/resources`.
Keep those in `local/` or pass them at runtime.

Template-only run from the repository root:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli --plan src\echelle_spectra\resources\spectrocube_plans\lhd_20250926_wmsr_batch.toml
```

Real dataset run with CLI path overrides:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli `
  --plan src\echelle_spectra\resources\spectrocube_plans\lhd_20250926_wmsr_batch.toml `
  C:\path\to\20250926 `
  -o C:\path\to\20250926\spectrocubes_wmsr_403nm
```

The batch terminal output starts with the source and destination, then reports
the receipt location and live count/rate/ETA. Add `--registry`,
`--calibrations`, and `--volume-label` so every source resolves one verified
snapshot and every record carries that identity:

```text
SpectroCube batch
Source:      C:\path\to\20250926
Destination: C:\path\to\20250926\spectrocubes_wmsr_403nm
Pattern:     *.SIF
Files:       53 (export)
Receipt:     local/runs/2026-08-13_12-00-00-20250926
Loading CMOS calibration...
Calibration ready.
[1/53] 193778_Echelle.SIF
[1/53] exported | 0.20 file/s | ETA 260.0s
[2/53] 193779_Echelle.SIF
Done. 53/53 exported successfully.
```

Ctrl-C marks the receipt interrupted. Running the same source/destination again
resumes the newest unfinished receipt and skips a file whose recorded identity
still matches on disk: the source by size and modification time, the output by
existence and size. Nothing is re-read, so resuming costs what remains rather
than what is already done, and the run says so in one line. See
[Durable campaign runs](campaign-runs.md) for the receipt schema and explicit
`--run-dir` workflow.

For several USB drives, pass one folder per drive and repeat the volume label:

```powershell
echelle process D:\2024\shots E:\2025\shots `
  -o F:\spectrocubes `
  --runs-dir F:\runs `
  --volume-label NIFS-2024 `
  --volume-label NIFS-2025 `
  --config local\campaign-export.toml `
  --registry F:\calibration\calibration_registry.toml `
  --calibrations F:\calibration\calibrations
```

The output root receives one child directory per source. Processing is
sequential within each source and concurrent across sources; receipts and
failures remain independent, while `echelle status --runs F:\runs` reconciles
the latest state of both targets.

The registry's ordered snapshot IDs point to the immutable snapshot
`[validity]` tables that define inclusive shot/date bounds. See
[Calibration epoch registry and cube provenance](calibration-epoch-registry.md)
for boundary, precedence, optional-field, and recovery details.

To preview a config-driven batch without a plan:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli C:\path\to\20250926 `
  --config src\echelle_spectra\resources\calibration_files\export_configs\lhd_cmos_20250926.toml `
  --pattern "*.SIF" `
  --dry-run `
  --verbose
```

## Quality check

After generation, open a representative cube in SpectroView and check:

- the viewer reports `calibration_type: absolute`
- the viewer reports `intensity_units: W/m2/nm/sr`
- the wavelength range starts near `403 nm`
- the low-wavelength edge no longer dominates the plot
- metadata includes dropped-column counts when non-finite/cropped columns were
  removed

For the first validated product, `spectroview --info` reported:

| Field | Value |
| --- | --- |
| Wavelength range | `403.00-801.86 nm` |
| Wavelength points | `42416` |
| Intensity shape | `frame=50, wavelength=42416` |
| Export config | `lhd_cmos_20250926` |

Later, SpectroView should add a selector for counts versus calibrated intensity.
Until then, treat these batch products as the calibrated `wmsr` product.
