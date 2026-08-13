# Live Calibration Bench

`echelle-calib` is a separate, lean pyqtgraph window for lamp work at the
instrument. It does not enlarge the established five-tab analysis GUI and does
not write calibration snapshots yet; snapshot save and procedure memory belong
to the next campaign stage.

## Start watching

Point the bench at the folder written by the Andor acquisition software:

```bash
echelle-calib path/to/acquisition-folder
```

The accepted 2025 CMOS order pattern and aligned wavelength table are packaged
defaults. Override either explicitly when rehearsing another detector or
historical calibration:

```bash
echelle-calib path/to/acquisition-folder \
  --pattern path/to/pattern.txt \
  --wavelength path/to/wavelength.txt
```

Load an existing fixture immediately while keeping the folder watch active:

```bash
echelle-calib path/to/2025-fixture-campaign --file path/to/lamp.sif
```

Run `echelle-calib --help` for polling, stability, saturation, and SNR controls.

## Safe file arrival

The watcher chooses the newest case-insensitive `.sif` file by modification
time, with filename as a deterministic tie-breaker. A candidate is loadable
only when:

1. its byte size and nanosecond modification time are identical across two
   consecutive polls by default; and
2. its modification time is at least one second old by default.

Growth or a timestamp change restarts the unchanged-poll count. One unchanged
fingerprint is emitted only once. A missing/unreadable folder and a SIF load
failure become explicit states instead of crashing the event loop. If a newer
file fails to load, the last good detector frame and its anchors stay visible;
the next stable file can recover the acquisition state. A successful new frame
starts a fresh anchor set because anchors belong to the measured frame.

These rules make polling repeatable and testable. They do not constitute a
measurement of real USB throughput or prove that a particular camera-writing
program closes files in a certain way.

## Watch-to-fit workflow

1. Wait for **File state: LOADED**. The upper plot shows the mean detector image
   and every accepted pattern trace; the selected order is highlighted.
2. Choose an order. The middle plot shows its frame-averaged extracted counts
   against raw detector column, with the packaged wavelength-table rows labeled.
3. Click near the measured peak corresponding to a label. The click chooses the
   nearest known row, and the existing single-Gaussian fitter searches around
   the clicked pixel.
4. Review the verdict. Raw 2-D detector pixels are checked before acceptance.
   Saturated windows are refused with lower-exposure guidance; a low-SNR or bad
   Gaussian is also refused and does not mutate the anchor set.
5. Repeat on another line/order. One accepted anchor is **COLLECTING**; two or
   more solve the established translation-plus-rotation model and enter
   **ALIGNED**.
6. Review dx, dy, rotation, detector-pixel RMS, and the per-anchor residual
   strip. Removing or replacing an anchor recomputes immediately. Dropping
   below two anchors returns to **COLLECTING**; clearing returns to **EMPTY**.

The model is deliberately rigid: translation in detector x/y plus small
rotation, with no stretch, shear, or warp. A low RMS only describes the selected
anchors. It does not replace wavelength-table QC or the Balmer/Fulcher validation
gate on plasma spectra.

## Tested state contract

The workflow logic lives in `echelle_spectra.calibration_bench`, independent of
Qt. Automated tests pin:

- empty, changing, stable, already-emitted, and failed watch states;
- growth reset and newest-file selection;
- good-load, failed-load preservation, and next-file recovery;
- anchor add, replace, remove, clear, saturated rejection, and failed-fit paths;
- waiting, empty, collecting, aligned, failed, and recovered alignment states;
- packaged detector/order loading with the accepted 2025 pattern and a
  historical bundled ThAr SIF fixture;
- off-screen detector, order, anchor, RMS, and residual GUI presentation.

The packaged fixture proves a repeatable software path over repository data. It
does not prove live hardware timing, acquisition-software behavior, lamp
visibility at NIFS, or sustained USB performance.
