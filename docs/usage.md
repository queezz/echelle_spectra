# The single-SIF viewer

```bash
echelle_spectra
```

`echelle_spectra` is the **viewer**: it opens one SIF at a time and shows you
what is in it — the detector image, the extracted orders, the calibrated
spectrum, and known-line overlays on top. Use it to look at a shot.

It is not where calibrations are made and not where drives are converted:

- **Making a calibration** → the live bench, `echelle-calib`.
  See [Live calibration bench](calibration-bench.md).
- **Converting SIFs to cubes and LHD text** → the `echelle` umbrella command.
  See [Where to start](usage-overview.md) and the
  [operator cheat sheet](operator-cheat-sheet.md).

---

## What the viewer shows

- The raw 2D echelle image, with the fitted order pattern over it.
- Each diffraction order extracted into a 1D spectrum.
- That spectrum on a wavelength axis, from the loaded calibration.
- Optional Balmer, Fulcher H₂, ThAr, Ne, and Hg line markers, each separately
  toggleable. All overlays start off; labels thin out in wide views and reveal
  local lines as you zoom. The viewer and `echelle-validate-lines` read the same
  provenance-carrying tables — see
  [Known-line overlays](known-line-overlays.md).
- The defined spectral windows for lines of interest — H-alpha, He-587,
  CII-515, and the rest. See [Band data](band-data.md).

---

## How a spectrum is built

The same five steps happen whether you are in the viewer, at the bench, or
running `echelle process`:

1. **Load the 2D image.** The raw detector frame from the SIF.
2. **Apply the order pattern.** Where each diffraction order sits on the
   detector. Only redetermined when the optics move.
3. **Extract the orders.** Slice each order out of the image into a 1D
   intensity array.
4. **Apply the wavelength solution.** Map detector column to wavelength through
   the per-order polynomial fit.
5. **Read the windows.** Locate and integrate the lines of interest.

Steps 2 and 4 are the calibration. Where the files behind them come from, and
what a finished cube records about them, is
[From calibration to cube](calibration-to-cube.md).

---

## Opening files through a saved calibration

By default the viewer reads every file through the packaged CCD and CMOS
tables — today's calibration. Frames from an earlier era were taken on an
instrument that has since moved, so today's tables put the order traces and the
line boxes beside the blobs instead of on them.

Point the viewer at that era's snapshot and it sees the files through those
eyes instead:

```powershell
echelle_spectra --calibration "<calibrations>\20190314_cmos"
```

Replace `<calibrations>` with your own snapshot root — `D:\NIFS\calibrations`
on Windows, `/Volumes/NIFS/calibrations` on macOS.

The folder is a calibration snapshot — the one holding `snapshot.toml` next to
`pattern.txt` and `wavelength.txt`. Snapshots are made by the bench and listed
by `echelle status`; see [Calibration snapshots](calibration-snapshots.md).

The **Calibration** selector in the control column is the in-GUI form of this —
it switches between the packaged CCD and CMOS tables and any snapshot folder you
browse to, live, without reopening the file — and the flag simply preconfigures
it at launch.

What changes:

- The title bar names the snapshot: **Echelle viewer — 20190314_cmos**, and the
  info panel repeats it beside each loaded frame. Without the flag, nothing
  about the window changes.
- The order pattern, the line overlays, the cursor link, and the spectrum's
  wavelength axis all come from that snapshot, because they all read the one
  calibration the window loaded.
- The camera is fixed to the detector the snapshot was made for, and the CCD /
  CMOS buttons are greyed out: a snapshot is one detector's calibration, and
  there is no second one loaded to flip to.

If the folder is not a snapshot, or the snapshot is missing a table or fails
its digest check, the viewer says which problem it hit and exits. It never
falls back to the packaged tables behind your back.

---

## When the instrument has moved a little

For a small detector or optics shift, fit a rigid correction rather than
recalibrating from scratch. The live bench does this interactively; the same fit
is available headlessly through
[Calibration alignment](calibration-alignment.md), which writes a **separate**
adjusted wavelength table. Never overwrite the historical calibration files.

Before exporting a dataset, run the
[Wavelength line validation](line-validation.md) gate on real plasma lines:
Balmer first, then Fulcher-alpha positions as supporting checks where the
features are visible and not dominated by blends.

---

## Notebooks

The notebooks in `examples/workflows/black_cmos/` are a **manual and tuning
reference**, kept for the cases the packaged commands do not cover — a new
optical setup, a pattern that needs to be talked into place by hand, a fit worth
inspecting step by step.

| Notebook | Use it for |
| --- | --- |
| `01_load_image.ipynb` | Sanity check: load an image, verify dimensions |
| `02_automated_pattern_extraction.ipynb` | Run packaged pattern extraction and compare traces |
| `02_pattern_calibration.ipynb` | Manual pattern tuning and debugging |
| `03_wavelength_calibration.ipynb` | Manual line identification and polynomial fit — new setup only |
| `04_extract_spectrum.ipynb` | Load → calibrate → extract → save, one file at a time |
| `05_calibration_alignment.ipynb` | Align an existing wavelength table against fresh lamp data |

For routine work, prefer the bench for calibration and `echelle process` for
conversion: both record provenance the notebooks do not. Historical notebooks
are in `examples/obsolete/`.
