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
| `calibration_source` | `"integrating sphere (echelle_spectra)"` (absolute modes only) |

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
