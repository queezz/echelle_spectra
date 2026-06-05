# Exporting to SpectroCube

## What is SpectroCube?

[SpectroCube](https://queezz.github.io/spectrocube/) is a neutral, self-describing
container for calibrated spectral data.  It stores intensity, a wavelength axis, and
structured metadata in a single [NetCDF-4](https://www.unidata.ucar.edu/software/netcdf/)
`.nc` file via [xarray](https://xarray.pydata.org/).

**Architectural split:**

| Package | Role |
|---------|------|
| `echelle_spectra` | Reads raw 2-D Echelle images, extracts order spectra, applies wavelength and absolute calibration — produces `Spectrum` objects. |
| `spectrocube` | Defines the neutral container/schema, validates it, and handles save/load. Knows nothing about Echelle optics. |

The conversion from `Spectrum` → `SpectroCube` lives here, in `echelle_spectra`,
because only this package knows what the calibration quantities mean.

---

## Quick start

```python
from echelle_spectra.tools.spectrocube_export import to_spectrocube, export_spectrocube

# --- assume `sp` is a calibrated Spectrum object ---

# Option A: convert to SpectroCube in memory
sc = to_spectrocube(sp, instrument_id="BlackEchelle")
print(sc)
# SpectroCube(instrument='BlackEchelle', dims=(frame=3, wavelength=3547), ...)

# Option B: convert and save to NetCDF in one call
sc = export_spectrocube(sp, "output/shot_042.nc", instrument_id="BlackEchelle")
```

Both functions are also importable directly from the package root:

```python
from echelle_spectra import to_spectrocube, export_spectrocube
```

---

## Choosing intensity units

The `units` keyword controls which calibrated quantity is stored:

| `units` | Attribute on `Spectrum` | `calibration_type` | `intensity_units` |
|---------|-------------------------|--------------------|-------------------|
| `"counts"` (default) | `spectrum.counts` | `"counts"` | `"counts"` |
| `"wm"` | `spectrum.wm` | `"absolute"` | `"W/m2/nm"` |
| `"wmsr"` | `spectrum.wmsr` | `"absolute"` | `"W/m2/nm/sr"` |
| `"phmsr"` | `spectrum.phmsr` | `"absolute"` | `"ph/s/nm/sr"` |

```python
# Export absolute spectral radiance
sc = export_spectrocube(sp, "shot_042_abs.nc", units="wmsr", instrument_id="BlackEchelle")
```

---

## Preserved metadata

The following fields from `Spectrum` are automatically written to SpectroCube attrs:

| Attr in `.nc` | Source |
|---------------|--------|
| `source_package` | always `"echelle_spectra"` |
| `created_at` | UTC timestamp at export time |
| `exposure_s` | `spectrum.info["ExposureTime"]` |
| `frame_interval_s` | `spectrum.info["CycleTime"]` |
| `source_file` | `spectrum.fpth` |
| `shot_number` | `spectrum.shotnumber` |
| `background_frames` | `spectrum.info["BackgroundFrames"]` (when non-empty) |
| `calibration_folder` | calibration resource folder used to load the spectrum |
| `calibration_order_pattern_file` | order-pattern file used for extraction |
| `wavelength_calibration_file` | wavelength lookup table used for extraction |
| `dropped_nonfinite_wavelength_columns` | count of invalid wavelength columns dropped before export |
| `calibration_source` | `"integrating sphere (echelle_spectra)"` (absolute modes only) |

!!! note "Non-finite absolute-calibration columns"
    Absolute calibration can produce non-finite columns where the sphere response
    is invalid or the low-wavelength edge is outside the useful calibration
    region. These columns are tolerated during `Spectrum` construction and are
    dropped during SpectroCube export by default. Keep
    `drop_nonfinite_columns = true` for routine calibrated exports unless you are
    deliberately debugging the calibration arrays.

Additional arbitrary metadata can be passed as keyword arguments:

```python
sc = export_spectrocube(
    sp,
    "shot_042.nc",
    instrument_id="BlackEchelle",
    notes="He-CD4 mixture, 5 Pa",
    grating="316 l/mm",
    slit_um=50,
)
```

---

## Reloading a saved SpectroCube

```python
from spectrocube import SpectroCube

sc = SpectroCube.load("output/shot_042.nc")
print(sc.wavelength)     # 1-D nm array
print(sc.intensity)      # (frame, wavelength) or (wavelength,) array
print(sc.ds.attrs)       # all metadata
```

---

## Multi-frame vs single-frame

When `spectrum.counts` has shape `(n_frames, n_wavelengths)`:

* **Multiple frames** (`n_frames > 1`) → stored as `(frame, wavelength)` 2-D.  
  The `frame` coordinate is an integer index `[0, 1, …, n_frames-1]`.
* **Single frame** (`n_frames == 1`) → squeezed to `(wavelength,)` 1-D by
  default (`squeeze_single_frame=True`).  Pass `squeeze_single_frame=False`
  to keep the 2-D shape.

Frame *time* is not stored as a separate coordinate because `echelle_spectra`
computes it externally as `trigger_delay + frame * CycleTime`.  Use the
`frame_interval_s` and `exposure_s` attrs to reconstruct it if needed.

---

## GUI — Save SpectroCube button

The viewer has a **Save SpectroCube** button in the Controls dock (below
the existing "Save spec" / "Save lines" checkboxes).

For CMOS/LHD data, the GUI defaults to the accepted 20250926 calibration context:

- `pattern_CMOS_20250926.txt`
- `alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt`
- `wavelength_medium = "air"` in exported SpectroCube metadata

**Workflow:**

1. Open a SIF file with the **Manual SIF load** button or by entering a shot
   number and clicking **Load selected**.
2. Wait for the image and spectrum to finish loading.
3. Select the desired **Spectra units** from the combo-box
   (`counts`, `wm`, `wmsr`, or `phmsr`).
4. Click **Save SpectroCube**.
5. A file-save dialog opens, pre-filled with a name like
   `193777_spectrocube.nc` in the configured output folder.
6. Choose or confirm the path and click **Save**.

The status bar shows `SpectroCube saved: /path/to/file.nc` on success, or a
dialog with an actionable error message on failure (including instructions if
the `spectrocube` package is not installed).

!!! note "spectrocube is optional"
    If `spectrocube` is not installed, a dialog explains what to do rather
    than crashing.  See [Local development installation](#local-development-installation).

---

## CLI — `echelle-spectrocube`

After installation the `echelle-spectrocube` command is available on the PATH.

For the current LHD CMOS dataset workflow, see
[Batch SpectroCube workflow](spectrocube-batch-workflow.md).

### Single-file export

```bash
echelle-spectrocube shot_193777_Echelle.SIF --units wm --wavelength-medium air -o output/shot_193777.nc
```

### Config-driven export

Stable camera/spectrometer/calibration choices can be stored in a calibration
config TOML, including wavelength crop bounds measured once for that calibration
context:

```bash
echelle-spectrocube local/193778_Echelle.SIF \
  --config src/echelle_spectra/resources/calibration_files/export_configs/lhd_cmos_20250926.toml
```

The 20250926 CMOS/LHD config crops the unstable low-wavelength edge below
`403.0 nm`, exports `wmsr` (`W/m2/nm/sr`), and records the original wavelength
range plus dropped-column counts in SpectroCube metadata.

### Plan-driven export

For repeatable single-file or batch generation, use a plan TOML that points at
the calibration config and supplies input/output paths:

```bash
echelle-spectrocube --plan src/echelle_spectra/resources/spectrocube_plans/lhd_193778_wmsr.toml
```

The current plan writes:

```text
local/193778_Echelle_spectrocube_wmsr_403nm.nc
```

For the 20250926 CMOS/LHD calibration pass, use the local sphere/background
pair and export absolute spectral radiance:

```powershell
echelle-spectrocube local\193778_Echelle.SIF `
  --units wmsr `
  --wavelength-medium air `
  --calibration-dir local\20250926_calib `
  --order-pattern pattern_CMOS_20250926.txt `
  --wavelength C:\path\to\calibration_files\alignments\Th_wavelength_CMOS_20240305_aligned_to_20250926.txt `
  --sphere sphere-0.1s-x3.sif `
  --sphere-background sphere-0.1s-x3-bg.sif `
  --integral C:\path\to\calibration_files\integrating_sphere.txt `
  -o local\193778_Echelle_spectrocube_wmsr.nc
```

### Batch folder export

```bash
# Export all .SIF files in a folder, save .nc files alongside them
echelle-spectrocube /data/shots/ --units wm

# Save to a separate output directory
echelle-spectrocube /data/shots/ --units wmsr -o /data/nc/

# Preview what would happen (no files written)
echelle-spectrocube /data/shots/ --dry-run --verbose
```

### Full option reference

| Option | Default | Description |
|--------|---------|-------------|
| `INPUT` | *(required)* | `.sif` file or folder |
| `--config` | *(none)* | Calibration/export config TOML |
| `--plan` | *(none)* | SpectroCube generation plan TOML |
| `--units` | `counts` | `counts`, `wm`, `wmsr`, or `phmsr` |
| `-o / --output` | same dir as INPUT | Output file (single) or directory (batch) |
| `--camera` | `CMOS` | Which bundled calibration to use: `CMOS` or `CCD` |
| `--calibration-dir` | bundled resources | Path to calibration files folder |
| `--order-pattern` | selected camera default | Override order-pattern file |
| `--wavelength` | selected camera default | Override wavelength lookup table |
| `--sphere` | selected camera default | Override integrating-sphere SIF |
| `--sphere-background` | selected camera default | Override integrating-sphere background SIF |
| `--integral` | selected camera default | Override integrating-sphere spectral table |
| `--instrument-id` | `echelle` | Stored in SpectroCube metadata |
| `--wavelength-medium` | `air` | Wavelength convention stored in SpectroCube metadata: `air` or `vacuum` |
| `--wavelength-min-nm` | config/default | Inclusive low-wavelength crop bound |
| `--wavelength-max-nm` | config/default | Inclusive high-wavelength crop bound |
| `--calibration-source` | config/default | Absolute calibration source metadata |
| `--no-drop-nonfinite-columns` | drop enabled | Keep non-finite wavelength columns instead of dropping them |
| `--output-suffix` | `_spectrocube` | Batch output suffix before `.nc` |
| `--pattern` | `*.SIF` | Glob for batch discovery (also tries `*.sif` as fallback) |
| `--overwrite` | *(skip existing)* | Replace existing output files |
| `--dry-run` | — | Print plan without writing |
| `--verbose` | — | Per-file current-progress output |

### Output filename convention

For each input file `shot_042_Echelle.SIF` the output is named
`shot_042_Echelle_spectrocube.nc` in the output directory.

### Calibration limitation

!!! warning "Calibration files required"
    The CLI uses the same Echelle order-pattern, wavelength calibration, and
    integrating-sphere files as the GUI.  These files must be present in the
    installed package under `resources/calibration_files/` (or supplied via
    `--calibration-dir`).

    For batch conversion the calibration is loaded **once** and reused for all
    files, so the cost is equivalent to opening the GUI once.

    If the `.sif` sphere files are not bundled with your installation (they are
    large binaries), point to a local copy:

    ```bash
    echelle-spectrocube /data/shots/ --calibration-dir /lab/calibration_files/
    ```

    The same `--camera` flag selects either the CCD or CMOS calibration file
    set (matching the radio button in the GUI).

    The default CMOS file set is the accepted 20250926/LHD context:
    `pattern_CMOS_20250926.txt` plus
    `alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt`.

---

## Local development installation

Both packages use a `src/` layout and are installed as editable packages.
In a fresh virtual environment:

=== "Unix / macOS"

    ```bash
    # 1. Install spectrocube (the neutral container)
    pip install -e /path/to/2026-spectrocube

    # 2. Install echelle_spectra with the spectrocube optional extra
    pip install -e /path/to/echelle_spectra[spectrocube]

    # Or, if spectrocube is already installed separately:
    pip install -e /path/to/echelle_spectra
    ```

=== "Windows PowerShell"

    ```powershell
    # 1. Install spectrocube (the neutral container)
    pip install -e C:\path\to\2026-spectrocube

    # 2. Install echelle_spectra with the spectrocube optional extra
    pip install -e C:\path\to\echelle_spectra[spectrocube]

    # Or, if spectrocube is already installed separately:
    pip install -e C:\path\to\echelle_spectra
    ```

Using the project venvs (if you already use `~/.venvs/echelle-spectra`):

=== "Unix / macOS"

    ```bash
    source ~/.venvs/echelle-spectra/bin/activate
    pip install -e /path/to/2026-spectrocube   # makes spectrocube importable
    ```

=== "Windows PowerShell"

    ```powershell
    & "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\Activate.ps1"
    pip install -e C:\path\to\2026-spectrocube   # makes spectrocube importable
    ```

---

## API reference

Full parameter lists and return types:
[SpectroCube export API reference](spectrocube-export-api.md).
