# Live Calibration Bench

`echelle-calib` is a separate pyqtgraph window for calibration work at the
instrument. It keeps acquisition, procedure, line-identification, sphere-factor,
configuration, and snapshot state outside the established five-tab analysis GUI.
The Qt callbacks are adapters over UI-independent campaign and alignment state.

## Start a campaign

Point the bench at the folder written by the acquisition software:

```bash
echelle-calib path/to/acquisition-folder
```

The accepted 2025 CMOS pattern and aligned wavelength table, packaged
integrating-sphere reference, historical 2024 sphere pair, and ThAr checklist
are defaults. Override them explicitly when rehearsing another detector or
campaign:

```bash
echelle-calib path/to/acquisition-folder \
  --pattern path/to/pattern.txt \
  --wavelength path/to/wavelength.txt \
  --integral path/to/integrating-sphere.txt \
  --previous-sphere path/to/previous-sphere.sif \
  --previous-sphere-background path/to/previous-sphere-bg.sif \
  --lamp ThAr --lamp Ne \
  --snapshot-id 20260901_cmos \
  --output-root path/to/calibrations \
  --config-root path/to/calibration-configs
```

Load one fixture immediately while keeping the folder watch active:

```bash
echelle-calib path/to/fixture-campaign --file path/to/lamp.sif
```

Run `echelle-calib --help` for polling, stability, saturation, SNR, campaign,
configuration, and snapshot controls.

## Safe file arrival and explicit roles

The watcher chooses the newest case-insensitive `.sif` file by modification
time, with filename as a deterministic tie-breaker. A candidate is loadable
only when:

1. its byte size and nanosecond modification time are identical across two
   consecutive polls by default; and
2. its modification time is at least one second old by default.

Growth or a timestamp change restarts the unchanged-poll count. One unchanged
fingerprint is emitted only once. A missing/unreadable folder and a SIF load
failure become explicit states instead of crashing the event loop. A failed
newer load preserves the last good detector frame and anchors; a later good
file recovers. A successfully loaded new frame starts a fresh anchor set because
anchors belong to that measurement.

Filename matching is help, not evidence. On the **Acquire** tab the bench may
suggest sphere signal, sphere background, lamp signal, or lamp background, but
the operator must press **Confirm role for loaded SIF**. Ambiguous names offer
several choices and tick nothing. Lamp roles also require an explicit ThAr, Ne,
Hg, or H2 family. This prevents a convenient but unsupported filename guess from
completing the procedure.

These rules make polling repeatable and testable. They do not measure real USB
throughput or prove how a particular acquisition writer closes files.

## On-screen procedure

The **Procedure** tab derives its check marks from campaign state in this order:

1. load the pattern, wavelength table, and integrating-sphere reference;
2. measure and explicitly classify the sphere signal;
3. measure and explicitly classify the matching sphere background;
4. compute candidate absolute factors and compare them with the previous pair;
5. for every requested lamp, classify a matching background and signal;
6. solve and review the lamp alignment with at least two accepted anchors;
7. generate the commented campaign, alignment, and export TOMLs;
8. save and validate the snapshot.

For practical lamp work, acquire/classify the background before the signal so
the newest frame visible for click fitting is the lamp signal. A checked item
always names its supporting filename or numerical result. Changing any
classification invalidates comparison/configuration/save readiness rather than
leaving stale green checks.

## Saturation and exposure guidance

Role confirmation inspects finite raw detector pixels, not the extracted curve.
The guidance panel reports the raw peak and one of:

- **SATURATED** — lower exposure and reacquire; fitted saturated windows remain
  refused before they can enter the anchor set;
- **DIM** — increase exposure toward a useful unsaturated peak;
- **GOOD** — continue to line identification and anchor fitting;
- **NO DATA** — reacquire because no finite raw pixels are available.

When SIF metadata carries an exposure time, the advice includes an approximate
next exposure targeting 70% of the configured saturation level. Background
frames instead tell the operator to keep the paired signal at the same exposure.
This is acquisition guidance, not a camera-timing validation.

## Line identification and live alignment

The lamp-fit view retains the Packet 4 interaction:

1. choose an order;
2. read the blue shared-catalog sticks and the dedicated **Line identification**
   table for ThAr, Ne, Hg, or Fulcher H2 provenance and approximate pixel
   positions;
3. click a gold curated calibration row near the measured peak;
4. review the raw-pixel saturation verdict and Gaussian centroid;
5. repeat on another line/order to solve translation plus rotation;
6. review dx, dy, rotation, RMS, and per-anchor residuals.

Shared-catalog positions are interpolated through the current wavelength-table
rows for identification help. They do not silently become fit anchors. Accepted
anchors still use the curated calibration table and established fitting,
saturation, detector-point, and rigid-transform functions. A low selected-anchor
RMS does not replace wavelength-table QC or Balmer/Fulcher plasma validation.

## Sphere factors and `insufficient data`

After the sphere pair is explicit, **Compute and compare factors** runs the
existing extraction and absolute-calibration engine twice: once for the measured
pair and once for the previous pair. The **Sphere factors** view overlays both
curves and reports the median and 5–95% range of their finite positive ratio.

If the previous pair is absent, fails to load, or provides fewer than 20 shared
finite positive samples, the state is **INSUFFICIENT DATA**. Candidate factors
remain visible and the procedure records that the comparison was attempted; the
bench never calls this aligned or invents a ratio. A candidate-factor failure is
**FAILED** and must be corrected before configuration or save.

The packaged default comparison uses the historical 2024 CMOS sphere pair
because no packaged 2025 sphere pair is available. It is software rehearsal
evidence, not a claim about a new lamp response.

## Generate TOMLs and save the snapshot

The **Save** tab writes a new identity directory containing:

```text
calibration-configs/20260901_cmos/
├── campaign.toml
├── alignment.toml
└── export.toml
```

All three files begin with explanatory comments, contain source filenames rather
than machine paths, parse as ordinary TOML, and remain freely hand-editable. The
campaign file records explicit roles, lamp families, exposure evidence, and the
sphere-comparison result. Alignment settings record the measured transform and
anchors. The export configuration uses the existing SpectroCube configuration
shape and snapshot role filenames. Existing configuration identities are not
silently overwritten; choose a corrected snapshot identity if a generated bundle
must change after review.

**Save and validate snapshot** calls `create_snapshot()` from the established
snapshot API. That API copies read-only source inputs into an atomic staging
directory, digests them, validates them, and refuses to replace an existing
identity. The bench then calls `load_snapshot()` explicitly for the completed
snapshot and reports **VALIDATED** only after that validator succeeds. Lamp signal
and background SIFs are both retained under `lamps/` with their explicit family.
The bench does not build a competing manifest or duplicate validation rules.

A save failure preserves classifications, comparison, anchors, and generated
TOMLs. Correct the identity or missing input and retry. The generated-TOML
identity must match the snapshot identity, preventing a corrected save from
quietly using stale configuration.

## Tested evidence boundary

Automated coverage pins watcher, classification, checklist, exposure, shared-line,
comparison, TOML, save, validation, failure, and recovery states without Qt;
focused off-screen tests cover every Packet 5 view. The complete fixture rehearsal
uses the accepted 2025 CMOS geometry/wavelength resources, packaged historical
ThAr and sphere SIFs, and an explicitly named synthetic lamp-background
placeholder because no 2025 lamp-background SIF is packaged.

That rehearsal proves the software path to a validated snapshot. It does not
prove live camera behavior, acquisition-writer timing, real USB throughput, NIFS
field operation, lamp response, or hardware performance.
