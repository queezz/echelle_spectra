# Calibration bench rehearsal

Rehearse the watch-to-fit interaction on a folder containing copied or linked
calibration fixtures (never modify the campaign originals):

```bash
echelle-calib path/to/fixture-folder --file path/to/fixture-folder/lamp.sif
```

The packaged 2025 CMOS pattern and wavelength alignment are the defaults. Use
`--pattern` and `--wavelength` for another campaign. The bench is read-only over
SIF inputs; this v1.3.0 interaction core does not save a snapshot.

Expected rehearsal:

1. file state reaches `LOADED`;
2. detector data and 29 order traces appear;
3. selecting an order updates the extracted spectrum and labeled rows;
4. click-to-fit either refuses a saturated/weak window or adds a clear anchor;
5. the second accepted anchor produces a live rigid solution, RMS, and residuals;
6. removing that anchor returns the state to `COLLECTING` without stale fit data.
