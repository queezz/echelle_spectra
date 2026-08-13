# Calibration bench campaign rehearsal

Rehearse the complete watch-to-validated-snapshot interaction on a folder
containing copied calibration fixtures (never modify campaign originals):

```bash
echelle-calib path/to/fixture-folder \
  --file path/to/fixture-folder/lamp.sif \
  --snapshot-id 20250926_cmos-fixture \
  --output-root path/to/rehearsal/calibrations \
  --config-root path/to/rehearsal/configs
```

The packaged 2025 CMOS pattern and wavelength alignment are the defaults. Use
`--pattern`, `--wavelength`, `--integral`, and previous-sphere options for
another campaign. The bench is read-only over source SIFs.

Expected rehearsal:

1. file state reaches `LOADED`, but no checklist role is credited until it is
   explicitly confirmed;
2. sphere signal/background and every required lamp signal/background tick their
   own procedure rows;
3. exposure guidance names the next safe measurement action;
4. shared packaged line help appears for the chosen lamp/order;
5. two accepted anchors produce the live rigid solution, RMS, and residuals;
6. sphere factors compare with the previous pair or report `insufficient data`;
7. commented campaign/alignment/export TOMLs are generated;
8. Packet 0 creates the snapshot and its validator reports `VALIDATED`.

The repository has no packaged 2025 Ne or lamp-background SIF. The automated
complete rehearsal therefore uses historical packaged ThAr/sphere files and an
explicitly named synthetic lamp-background placeholder. It proves state and file
contracts, not a 2025 lamp-response or hardware measurement.
