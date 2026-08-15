# Live Calibration Bench

`echelle-calib` is a separate pyqtgraph window for calibration work at the
instrument. It keeps file triage, procedure, line-identification, sphere-factor,
configuration, and snapshot state outside the established five-tab analysis GUI.
The Qt callbacks are adapters over UI-independent campaign and alignment state.

The bench takes whatever you throw at it. It has no required lamp, no required
filename, and no required folder: files arrive by hand, roles are assigned by
hand, and the procedure is derived from what you actually loaded.

## Start a campaign

```bash
echelle-calib
```

The optional folder argument only chooses where the **Add SIF files…** dialog
opens first:

```bash
echelle-calib path/to/todays-calibration-folder
```

The accepted 2025 CMOS pattern and aligned wavelength table, packaged
integrating-sphere reference, and historical 2024 sphere pair are defaults.
Override them explicitly when rehearsing another detector or campaign:

```bash
echelle-calib path/to/calibration-folder \
  --pattern path/to/pattern.txt \
  --wavelength path/to/wavelength.txt \
  --integral path/to/integrating-sphere.txt \
  --previous-sphere path/to/previous-sphere.sif \
  --previous-sphere-background path/to/previous-sphere-bg.sif \
  --lamp Ne --lamp ThAr \
  --snapshot-id 20260901_cmos \
  --valid-from 2026-09-01 \
  --output-root path/to/calibrations \
  --config-root path/to/calibration-configs
```

`--lamp` names lamps to *suggest*, not to demand; any name is accepted and none
is ever required. Load files at start-up with a repeatable `--file`, and enable
the optional folder watch with `--watch`:

```bash
echelle-calib path/to/calibration-folder --file path/to/lamp.sif --watch
```

Run `echelle-calib --help` for polling, stability, saturation, SNR, campaign,
configuration, and snapshot controls.

## Reading the bench

Text here is data, not decoration, and the bench is built to be read across a
full-screen window without leaning in:

- readings are sized in **points relative to the platform font**, with floors,
  so a lower-resolution display scales the surface rather than shrinking it.
  Verdict headlines, RMS, anchor counts, save/alignment states, and the factor
  result are rendered larger and bolder than the prose around them;
- the pane divider is a **real splitter with no width ceiling**: drag the
  controls as wide as the reading needs. Triage lines, advice, checklist rows,
  and file cells **wrap**; nothing important is ellipsized;
- every verdict and parameter carries a **one-line tooltip**, and clicking or
  focusing one — or any checklist row, file row, or anchor — writes the full
  explanation into the **Why this reading** dock at the bottom. The dock never
  reacts to the pointer merely passing over a widget: it changes when it is
  asked to, so a sentence cannot be yanked away mid-read.

The bench window carries its own icon, derived at run time from the same
`echelle.png` the main GUI uses — tinted and badged rather than redrawn — so the
two windows are tellable apart in a taskbar without a second artwork file.

## Getting files onto the bench

Primary input is manual, because at acquisition time files are misnamed,
renamed later, appear after the fact, and get retried:

- **drag and drop** one or many SIFs onto the window — the same drop machinery
  the main GUI uses; dropping a folder queues the SIFs inside it;
- **Add SIF files…** opens an ordinary file dialog in the current folder;
- `--watch` additionally polls a folder and loads each newly stable SIF. This is
  a convenience, never a requirement.

Every file is read on a worker thread, one at a time, and lands as a row in the
**Files** table. The frames themselves are not retained: only the verdict is, so
a whole campaign folder can sit on the bench at once. **Open selected file for
lamp fitting** re-reads one row when you want to click lines on it.

A SIF that cannot be read becomes an explicit failure message and leaves the
last good frame in place; the queue continues with the next file.

When the watcher is enabled it chooses the newest case-insensitive `.sif` file
by modification time, with filename as a deterministic tie-breaker, and emits it
only after its byte size and nanosecond modification time repeat across two
consecutive polls and it is at least one second old. These rules make polling
repeatable and testable. They do not measure real USB throughput or prove how a
particular acquisition writer closes files.

## Exposure triage is the front door

Triage needs nothing but a file — no role, no lamp, no folder. As soon as a
frame is read, the **Exposure triage** view shows one verdict line, the raw
counts histogram, and a second histogram of the last 10% before full scale, so
the top end stays legible where the background peak would otherwise dwarf it.
When no pixel reaches that last 10% there is no distribution to draw, and the
panel says so in words — *No pixels within 10% of full scale.* — instead of
painting an empty log histogram as a solid block. This replaces squinting at
the acquisition software: shoot, drop the file, one glance, adjust the lamp,
shoot again.

Saturation is judged by clustering alone:

- a **connected cluster of two or more full-scale pixels** (four-connectivity
  inside one frame) is real saturation;
- an **isolated full-scale pixel is an anomaly** — a cosmic ray or a hot pixel —
  counted and reported, never held against the exposure. A pixel that saturates
  in every frame of a stack stays a repeated anomaly, because clusters are never
  linked across frames.

The verdict states:

- **SATURATED** — how many full-scale pixels in how many clusters; lower the
  exposure and shoot again;
- **TOO DIM FOR LINES** — the brightest real pixel as a percentage of full scale
  and as a multiple of the frame's own noise floor: right for a background, too
  weak for a lamp;
- **HEALTHY** — headroom from the brightest **non-anomalous** pixel, as
  "X% of full scale — about N× brighter is still safe";
- **NO DATA** — no finite raw pixels; reacquire.

The noise floor is a robust median/MAD estimate over a strided pixel sample of
the frame. The verdict is acquisition guidance, not a camera-timing validation.

### Saturation reads differently once a frame has a role

Triage is the front door and knows nothing about roles. Once a file carries
one, the same measurement is read again in the light of it, and saturation is
the only verdict a role can change:

- a **lamp signal or lamp background** is never failed for saturating. A
  bright/dim pair is shot precisely so the strong lines clip in the dim series
  while the weak ones emerge, so the frame reads
  `saturated in N cluster(s) — fit unsaturated lines only` and stays usable.
  The fit is guarded per anchor instead: a click on a line whose own raw
  detector window holds full-scale pixels is refused by name, and the
  unsaturated lines beside it fit normally;
- a **sphere signal or sphere background** keeps the hard verdict, because a
  clipped sphere pixel is an unknown number rather than a bright one and no
  absolute factor can be computed from it.

## Roles are assigned by hand, one control per file

Each row of the **Files** table carries its own role control and its own lamp
name: sphere signal, sphere background, lamp signal, lamp background, or
experiment/other for a frame that is simply parked on the bench. The lamp name
offers ThAr, Ne, Hg, and H2 as ready-made choices and accepts any other name you
type — the known list is convenience, never a permitted set.

Filename patterns only **pre-fill** these controls. `Ne-0.02s-x3-bright-lines_bg.sif`
pre-selects a Ne lamp background; `IMG_0042.sif` pre-selects nothing and still
takes any role you pick. A filename can never complete the procedure by itself.

Because a pre-filled control is a guess and not an assignment, it says so. An
unconfirmed row reads `Sphere signal · SUGGESTED` in amber, its file cell reads
`SUGGESTED ONLY — no role assigned`, and the Procedure tab names the file rather
than claiming nothing exists: "sphere-0.1s-x3.sif is only suggested by its
filename". Confirm one row by picking its role in the control — picking the
entry already shown counts — or press **Confirm N suggested role(s)** to assign
every unconfirmed suggestion of the folder in one deliberate step. Any of them
can then be changed freely.

## On-screen procedure

The **Procedure** tab is built from the files and roles that exist, not from a
fixed lamp list:

1. pattern, wavelength table, and integrating-sphere reference are loaded;
2. SIFs are on the bench and triaged;
3. the sphere signal and its background carry their roles;
4. candidate absolute factors are computed and compared with the previous pair —
   the sphere pair alone unblocks this, no lamp needed;
5. for **every lamp you actually assigned**, its signal and background rows;
6. the lamp alignment is solved with at least two accepted anchors;
7. the commented campaign, alignment, and export TOMLs are generated;
8. the snapshot is saved and validated.

A Ne-only campaign therefore reaches "snapshot saved and validated" with every
tick green and no mention of any lamp it did not measure. The previous
campaign's lamps appear as one non-blocking advisory row — "last time: Ne;
consider ThAr, Hg, H2 if available" — because every extra lamp measured at the
bench is a future refinement, and one calibration is enough to ship.

Every row that is not yet possible names what would unblock it, so no surface is
ever dead: the sphere factors need only the sphere pair, the lamp rows need only
any one classified lamp, and triage needs only a file. With nothing loaded, the
drop target and the **Add SIF files…** button are the primary surface.

A checked item always names its supporting filename or numerical result.
Changing any classification invalidates comparison/configuration/save readiness
rather than leaving stale green checks.

## Exposure guidance for a classified frame

Once a file carries a role, its guidance panel restates the triage verdict as a
next acquisition action. When SIF metadata carries an exposure time, the advice
includes an approximate next exposure targeting 70% of the configured saturation
level. Background frames instead tell the operator to keep the paired signal at
the same exposure. A saturated lamp frame is never told to lower its exposure:
it is told that the clipped strong lines are expected and that the unsaturated
lines are the ones to fit.

## Line identification and live alignment

The lamp-fit tab carries an explicit **Order** control — a spin box with ◀/▶
step buttons and the highest order beside it — and a **Fit on** selector that
chooses between the mean of the acquisition and any one of its frames. A 3-frame
acquisition therefore fits on the quieter mean by default and drops to a single
frame when one frame is spoiled or when a line clips in one frame and not the
others; saturation is checked on the raw pixels of whichever frame is being
fitted. Changing the selection clears the collected anchors, because a centroid
belongs to the spectrum it was measured on.

The lamp-fit view retains the Packet 4 interaction:

1. choose an order and which frame to fit;
2. read the labelled sticks on the spectrum, or the same lines listed in the
   **Expected lines** panel down the left column — selecting a row marks its
   stick;
3. click a labelled line near the measured peak to anchor it;
4. review the raw-pixel saturation verdict and Gaussian centroid;
5. click an anchored stick again — with either mouse button — to take the
   anchor back off; the anchor table follows immediately;
6. repeat on another line/order to solve translation plus rotation;
7. review dx, dy, rotation, RMS, and per-anchor residuals.

**There is one expected-line list.** The labelled sticks, the rows of the
Expected lines panel, and the count in that panel's header are three renderings
of one tuple: the assigned lamp's own wavelength-table rows falling in the
selected order. The packaged NIST caches annotate those rows — relative
intensity and provenance — where they happen to carry the wavelength, and they
never add or remove a line. Earlier rounds drew sticks from the curated table
while filling the panel from the caches interpolated onto the order, which is
how order 7 could label a NeI 640.225 the panel had never heard of and order 6
could show three labelled Ne sticks under "0 expected Ne lines in this order":
the packaged Ne cache spans 580.4–638.3 nm and the curated table does not.

Accepted anchors use the curated calibration table and the established fitting,
saturation, detector-point, and rigid-transform functions. A low selected-anchor
RMS does not replace wavelength-table QC or Balmer/Fulcher plasma validation.

## Anchors reference the lamp you assigned

The gold rows a click may snap to are **the assigned lamp's own lines and
nothing else**. The panel above the alignment box states which catalog is
scoping the fit, for example:

```text
39 of 184 curated rows are Ne lines (NeI, NeII); anchors reference Ne only
```

The lamp comes from the file currently open for fitting; if that file carries no
lamp role, the bench falls back to the session's single lamp. Assigning a role
or opening another file re-scopes the fit immediately, and any anchor that was
accepted against the previous lamp's rows is dropped rather than carried into
the new transform.

This matters because the packaged wavelength table is a **ThAr** table: its 184
curated rows are mostly thorium and argon, and every neon line has a thorium row
a few pixels away. Measuring neon centroids against those neighbours is how the
same real 2025 folder produced RMS 9.37 px from 125 anchors, against the
0.54 px that its 39 curated Ne rows give.

Three honesty states replace the old silence:

- **Scoped.** The panel names the catalog and the species it selected on.
- **Mismatched line help.** When the **Line help** combo shows a different
  catalog from the assigned lamp, the panel adds
  `WARNING — the ThAr line help on screen is not Ne's own catalog; anchors are
  still referenced against Ne lines only`. The blue sticks and the clickable
  rows disagree on purpose; the fit is unaffected.
- **No catalog.** A free-text lamp the packaged catalogs have never met states
  `no line catalog for Kr — anchors cannot be auto-referenced`, refuses clicks
  with that same reason, and turns the procedure's alignment row into an
  attention item. A lamp whose catalog exists but whose species this particular
  table never carries says that instead.

The procedure row names the catalog that anchored the fit, so the evidence
travels with the checklist:

```text
✓  Lamp alignment solved and reviewed — 39 anchors vs Ne catalog, RMS 0.54 px
```

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

## The saved snapshot carries the measured calibration

```text
calibrations/20260901_cmos/
├── snapshot.toml
├── pattern.txt
├── wavelength.txt   # base table moved by the solved transform
├── alignment.toml   # transform, RMS, fitted-line count
├── integral.txt
├── sphere.sif
├── sphere_bg.sif
└── lamps/
```

The snapshot's `wavelength.txt` is not a copy of the base table. Before
`create_snapshot()` runs, the bench applies the solved rigid transform to every
curated row of the base table with the same `apply_rigid_correction_to_lines()`
and `write_wavelength_table()` functions `echelle-align --save` uses, and hands
that corrected table to the snapshot API. The written header names the base
wavelength and pattern files, the alignment dataset, the lamp, signal,
background, sphere and sphere-background filenames, the transform, the RMS, and
the number of fitted lines.

A transform that moves no row measurably would only reformat the table, so the
base bytes are copied instead. The message panel and the procedure checklist say
which of the two happened and how far the largest row moved, and the manifest's
`[alignment]` table records `wavelength_correction_applied` alongside the solved
`dx_px`, `dy_px`, and rotation.

The snapshot folder also gains an `alignment.toml` in the established alignment
settings shape, so `load_alignment_settings()` reads the transform back from the
snapshot itself. The reviewable configuration bundle under `--config-root` keeps
its own richer `alignment.toml`, which additionally lists every accepted anchor.

Every saved snapshot carries a `[validity]` table so it is usable in a
calibration epoch registry without hand editing. The epoch is open-ended and
starts on the acquisition date; `--valid-from YYYY-MM-DD` states another start,
and the default is today.

A save failure preserves classifications, comparison, anchors, and generated
TOMLs. Correct the identity or missing input and retry. The generated-TOML
identity must match the snapshot identity, preventing a corrected save from
quietly using stale configuration.

## Tested evidence boundary

Automated coverage pins triage, watcher, classification, checklist, exposure,
shared-line, comparison, TOML, save, validation, failure, and recovery states
without Qt; focused off-screen tests cover the bench views, the drag-and-drop
path, and per-file role assignment. Triage verdicts are pinned against synthetic
frames for clustered saturation, lone cosmic-ray pixels, a too-dim frame, and a
healthy one, and a fixture folder mirroring the real 2025 Ne-only campaign
reaches a saved, validated, registrable snapshot with no ThAr anywhere. Tests
also pin that a nonzero transform changes the saved table's digest and reproduces the solved
shift when the table is read back, that a transform moving nothing copies the
base bytes, that the saved snapshot resolves through a calibration registry
unedited, and that the snapshot's `alignment.toml` round-trips. The complete fixture rehearsal
uses the accepted 2025 CMOS geometry/wavelength resources, packaged historical
ThAr and sphere SIFs, and an explicitly named synthetic lamp-background
placeholder because no 2025 lamp-background SIF is packaged.

That rehearsal proves the software path to a validated snapshot. It does not
prove live camera behavior, acquisition-writer timing, real USB throughput, NIFS
field operation, lamp response, or hardware performance.

The bench has additionally been driven head-less over the real 2025-09-26
Ne-only folder — six SIFs dropped through the real input path, triaged, given
roles by hand, fitted, compared, saved, validated, and resolved through a
registry. That run exercises the software against real detector data; it is
still not a live-instrument or field-operation claim.

Repeating that run with lamp-scoped anchors is what measures the scoping: the
same folder, the same bright Ne frame, every curated Ne row clicked, gives 39
accepted anchors at RMS 0.536 px with all 39 within 2 px of their expected
pixel, where the lamp-blind lookup gave 125 anchors at RMS 9.372 px with 47
within 2 px. The solved shift is +0.024 px, which says the packaged table was
already aligned to this data and the old 9 px scatter was an artefact of
measuring neon against thorium.

The same folder has since been driven again through the Files-table role
controls themselves rather than through the domain functions. Every dropped file
pre-fills correctly and every control reports `SUGGESTED` until confirmed;
confirming the sphere pair in the Role column computes real absolute factors,
`READY · new/previous median 0.422; 5–95% 0.242–0.484 (42471 samples)`. On the
real `Ne-0.1s-x3-dimm-lines.sif`, which saturates in 33 clusters by design, the
frame stays usable and the per-anchor guard refuses 11 clicks by line name while
accepting 28 unsaturated ones, solving at RMS 0.565 px; fitting frame 2 of that
3-frame acquisition alone accepts 29 and solves at RMS 0.559 px.
