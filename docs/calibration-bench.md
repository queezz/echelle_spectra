# Live Calibration Bench

`echelle-calib` is a separate pyqtgraph window for calibration work at the
instrument. It keeps file triage, procedure, line-identification, sphere-factor,
configuration, and snapshot state in their own window, separate from the
single-SIF viewer. The Qt callbacks are adapters over UI-independent campaign
and alignment state.

The bench takes whatever you throw at it. It has no required lamp, no required
filename, and no required folder: files arrive by hand, roles are assigned by
hand, and the procedure is derived from what you actually loaded.

## Start a campaign

```bash
echelle-calib
```

Name today's calibration folder and the bench uses it twice: it is where the
**Add SIF files…** dialog opens first, and it is what both output roots are
derived from (see [Where the bench writes](#where-the-bench-writes)).

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
  --config-root path/to/calibrations/configs
```

`--lamp` names lamps to *suggest*, not to demand; any name is accepted and none
is ever required. Load files at start-up with a repeatable `--file`, and enable
the optional folder watch with `--watch`:

```bash
echelle-calib path/to/calibration-folder --file path/to/lamp.sif --watch
```

Run `echelle-calib --help` for polling, stability, saturation, SNR, campaign,
configuration, and snapshot controls.

## Where the bench writes

Both output roots are **derived from the folder argument**. Launched at
`T:\2025-LHD-BH\Echelle\20250926_calib`, the bench writes:

```text
…\20250926_calib\calibrations\<snapshot-id>\          the immutable snapshot
…\20250926_calib\calibrations\configs\<snapshot-id>\  the settings bundle
```

The folder argument is the calibration's own folder: the lamp and sphere frames
that were measured, and the calibration computed from them, side by side in the
one folder that holds it all. It is also where you go looking for what the bench
wrote — so that is what the roots hang off, and nothing is scattered into
whatever directory a shortcut happened to start in. With no folder argument at
all there is nothing better than the working directory, and the bench falls back
to it.

Because that folder already holds the raw frames, the snapshot **references
them where they lie and never copies them**: it records the path back to each
sphere and lamp SIF beside the SHA-256 and size that say the frames are still
the ones it measured.

Everything the bench generates lives under that single `calibrations\` folder:
the snapshot identity folders, and one `configs\` subfolder holding the settings
bundles. A calibration folder therefore gains exactly one generated child, not
two. Nothing that reads a calibrations root is confused by `configs\`: snapshot
enumeration keys on a child's `snapshot.toml`, and the epoch registry resolves
the snapshot IDs it was given by name.

`--output-root` and `--config-root` still override each root independently, and
an overriding path is made absolute so the window can be honest about it. A
network path stays intact: only path joins are used, so a `\\server\share`
root keeps its two leading separators.

Either way the **Save** tab states both roots in full — shortened in the middle
only when they will not fit, one hover or one click into the *Why this reading*
dock away from being read exactly. Once a snapshot is saved, an **Open folder**
button appears beside the confirmation and shows the saved folder in the
system's own file manager, so no path has to be retyped.

## The identity is dated by acquisition, not by computation

A calibration belongs to the day its images were **taken**. When it is computed
is bookkeeping, not physics — and `20240305_cmos`, `20250926_cmos` and every
other epoch in the registry are already named that way. So the bench prefills
the snapshot identity as `YYYYMMDD_<detector>` from the acquisition date, and
looks for it in this order:

1. a leading `YYYYMMDD` in the launch folder's own name — `20250926_calib`
   gives `20250926_cmos`;
2. failing that, the acquisition date the loaded SIF headers carry
   (`ExperimentTime`), filled in as soon as the first frame lands;
3. failing both, today's date as a placeholder, with the Save tab's *Why this
   reading* text saying no acquisition date could be derived and that it must
   be checked before saving.

`--snapshot-id` always wins, and so does typing in the **ID** field: once you
have typed there, a later frame never overwrites it. The epoch a saved snapshot
declares follows the identity's own date, so a snapshot computed months later
still starts its epoch on the day it was measured; `--valid-from YYYY-MM-DD`
overrides that.

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
two windows are tellable apart in a taskbar without a second artwork file. It is
built at every size a shell asks for, down to the 16 px the title bar uses, and
on Windows `echelle-calib` claims its own **AppUserModelID** before opening the
window. Without that claim the shell files the process under whatever launched
it and shows the launcher's icon on the taskbar however the window is dressed —
which is why the bench icon kept "disappearing" while the main GUI, which has
always made the claim, kept its own.

## How the window spends its space

The bench opens at the smallest size it is actually usable at, never smaller.
At the default geometry the **file table**, the **Bench state** readings, the
**Sphere factors** line and at least six rows of the **Expected lines** table
are all on screen at once, and the controls column reaches every one of them
without a scrollbar. Three splitters divide the surface and every one of them
is draggable; whatever you drag is what the rest of the session opens with:

- **controls | views** across the window;
- **control tabs | Expected lines** down the left column, so the line table
  keeps a real share of the height instead of the leftovers;
- **readings | plots** down the right column. Bench state and the factors line
  live in a strip across the top of the view column rather than stacked below
  the file table, so they are in view whichever control tab is open. On a
  narrow window the strip folds onto a second row rather than shaving the
  controls inside it.

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
7. the alignment settings and the campaign around them are saved;
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

Once a file carries a role, the guidance line under the **Exposure triage**
histograms restates the triage verdict as a next acquisition action. When SIF metadata carries an exposure time, the advice
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
the packaged Ne cache then spanned only 580.4–638.3 nm and the curated table
does not. The caches now cover 380–810 nm, so they annotate nearly every
curated row — but the panel is still the curated rows, and would be even if the
caches were empty.

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

After the sphere pair is explicit, **Compute factors** runs the
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

## Save the alignment settings, then save the snapshot

The **Save alignment settings** button comes first: it writes the solved
alignment, the roles you assigned, and the campaign's own settings into a new
identity directory. Existing files are never overwritten silently — the bench
says they are already there and offers a separate, deliberate **Regenerate**
press, because that press discards any hand edits.

The directory it writes contains:

```text
calibrations/configs/20260901_cmos/
├── campaign.toml
├── alignment.toml
└── export.toml
```

All three files begin with explanatory comments, contain source filenames rather
than machine paths, parse as ordinary TOML, and remain freely hand-editable. The
campaign file records explicit roles, lamp families, exposure evidence, and the
sphere-comparison result. Alignment settings record the measured transform and
anchors. The export configuration uses the existing SpectroCube configuration
shape and snapshot role filenames.

The export file also carries the LHD timing and crop the bench does not measure
— `trigger_delay_s` above all — carried forward from the previous campaign and
marked as inherited. Review those values before the next LHD campaign; see
[From calibration to cube](calibration-to-cube.md).

**Save and validate snapshot** calls `create_snapshot()` from the established
snapshot API. That API writes the computed files into an atomic staging
directory, digests every input, validates the result, and refuses to replace an
existing identity. The bench then calls `load_snapshot()` explicitly for the
completed snapshot and reports **VALIDATED** only after that validator succeeds.
Every raw frame — the sphere pair and both the signal and background SIF of each
lamp — is digested where it lies and recorded as a reference with its explicit
family. The bench does not build a competing manifest or duplicate validation
rules.

## The saved snapshot carries the measured calibration

```text
20260901_calib/                        the calibration folder you launched on
├── sphere-0.1s-x3.sif                 the frames stay here, and only here
├── sphere-0.1s-x3-bg.sif
├── ThAr-0.3s-x3.sif
├── ThAr-0.3s-x3-bg.sif
└── calibrations/20260901_cmos/
    ├── snapshot.toml   # …and this points back at them, with their digests
    ├── pattern.txt     # base pattern moved by the same solved transform
    ├── wavelength.txt  # base table moved by the solved transform
    ├── alignment.toml  # transform, RMS, fitted-line count
    └── integral.txt
```

The snapshot holds what the bench computed and points at what it read. Moving
or copying the whole calibration folder keeps every reference true, because a
reference to a frame in the same folder is stored relative to the snapshot
(`../../sphere-0.1s-x3.sif`). Moving the snapshot folder out on its own does
not, and validation says so by the full path it looked in.

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

The snapshot's `pattern.txt` is corrected with that same transform, through
`apply_rigid_correction_to_pattern()` and the shared
`write_corrected_pattern_table()` that `echelle-align --save` also writes its
pattern with, and the manifest records `pattern_correction_applied` and
`pattern_max_shift_px` beside the wavelength pair. **The bench's own live view
and the snapshot it saves disagree by design: the bench measures the raw frames
against the *reference* pattern and table, which is what makes a shift
measurable at all, while the snapshot describes the *data* — both files moved
onto the detector the bench just measured.** A transform that moves no trace by
a measurable row copies the base pattern's bytes instead, byte for byte.

The snapshot folder also gains an `alignment.toml` in the established alignment
settings shape, so `load_alignment_settings()` reads the transform back from the
snapshot itself. The reviewable configuration bundle under `--config-root` keeps
its own richer `alignment.toml`, which additionally lists every accepted anchor.

Every saved snapshot carries a `[validity]` table so it is usable in a
calibration epoch registry without hand editing. The epoch is open-ended and
starts on the acquisition date the snapshot identity is named for — see
[The identity is dated by acquisition](#the-identity-is-dated-by-acquisition-not-by-computation)
— so the day of computation never leaks into it; `--valid-from YYYY-MM-DD`
states another start.

A save failure preserves classifications, comparison, anchors, and the settings
already written. Correct the identity or missing input and retry. Those settings
must carry the same identity as the snapshot, which prevents a corrected save
from quietly using the previous attempt's configuration.

## What happens to these files next

Everything the bench saved is now input to processing — but not all of it, and
not equally. Which file supplies which part of a finished cube, and what a cube
records about the calibration behind it, is
[From calibration to cube](calibration-to-cube.md).

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
