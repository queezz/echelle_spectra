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

The active dataset batch plan is:

```text
src/echelle_spectra/resources/spectrocube_plans/lhd_20250926_wmsr_batch.toml
```

It points at the 20250926 LHD dataset folder and writes calibrated cubes to a
dedicated `spectrocubes_wmsr_403nm` output folder.

Run it from the repository root:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli --plan src\echelle_spectra\resources\spectrocube_plans\lhd_20250926_wmsr_batch.toml
```

To preview a config-driven batch without a plan:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m echelle_spectra.spectrocube_cli C:\Users\queez\Dropbox\Experiments\2025-LHD-BH\Echelle\20250926 `
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

For the first validated product, `spectroview --info` reported:

| Field | Value |
| --- | --- |
| Wavelength range | `403.00-801.86 nm` |
| Wavelength points | `42416` |
| Intensity shape | `frame=50, wavelength=42416` |
| Export config | `lhd_cmos_20250926` |

Later, SpectroView should add a selector for counts versus calibrated intensity.
Until then, treat these batch products as the calibrated `wmsr` product.
