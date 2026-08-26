# Echelle operator cheat sheet

This is the short manual for the installed surfaces in Echelle Spectra 1.6.
It separates the portable NIFS kit from a source checkout because their command
prefixes are different. The portable kit does **not** add commands to the
machine-wide `PATH`.

If you are not sure which of the three surfaces your job belongs to, read
[Where to start](usage-overview.md) first; it routes by task in one table.

Lab is queezz's own private development helper for switching environments on
his own machines. It is not an Echelle dependency, it is not required on an
institute workstation, and nothing here needs it. For the trip, use the
portable-kit section below. Administrator access may be available there, but
the kit does not need to use it or change the machine-wide Python.

The portable 1.6 kit is the current released travel kit. Verify the installed
version before a campaign and keep campaign data outside the replaceable kit.

## First identify how Echelle is installed

### Portable kit on Windows

Open PowerShell in the unpacked kit folder. Install once, then use the kit
wrapper for campaign commands:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\install.ps1 offline
.\echelle.ps1 status
```

The viewer, the live bench, and the specialist commands live inside the
kit-local environment:

```powershell
& ".\.venv\Scripts\echelle_spectra.exe"
& ".\.venv\Scripts\echelle-calib.exe" "D:\NIFS\incoming"
& ".\.venv\Scripts\echelle-pattern.exe" --help
```

### Portable kit on macOS

Open Terminal in the unpacked kit folder:

```bash
chmod +x install.sh echelle bin/uv
./install.sh offline
./echelle status
```

The other installed commands are under `.venv/bin`:

```bash
./.venv/bin/echelle_spectra
./.venv/bin/echelle-calib "/Volumes/NIFS/incoming"
./.venv/bin/echelle-pattern --help
```

Use `online` instead of `offline` only when intentionally testing the
network-backed installer. Normal use after installation needs neither command.

### Owner's Lab-managed development checkout

Use Lab to select the external environment:

```powershell
lab activate echelle-spectra
```

After activation, the prompt is using the right Python and ordinary commands
such as `echelle` should resolve directly. Lab manages environment selection;
it does not reinstall the project.

Refresh the editable install once after pulling a version that declares a new
console command:

```powershell
uv pip install --no-deps -e .
echelle --help
echelle-calib --help
```

For a one-command check without changing the current shell:

```powershell
lab run echelle-spectra echelle status
```

Direct activation and absolute `.venvs` paths are troubleshooting fallbacks,
not the normal workflow. If Lab has not yet installed its shell integration,
`lab activate echelle-spectra` prints the PowerShell dot-source command to run;
`lab shell-init --help` explains the permanent setup.

### Ordinary source checkout without Lab

Lab is never required to build or run the source distribution. In a normal
Python environment, install the checkout with its development dependencies and
then use the same bare commands:

```bash
python -m pip install -e ".[dev]"
echelle status
echelle_spectra
```

This route assumes a suitable Python installation and is meant for development.
The portable kit above is the institute-workstation route when Python and
developer tools cannot be assumed.

## What each surface does

| Command | Use it for | Normal trip use |
| --- | --- | --- |
| `echelle_spectra [--calibration DIR]` | Single-SIF viewer: open one shot, see the image, the extracted orders, the calibrated spectrum, and the known-line overlays (`--calibration` points it at a saved snapshot folder, so old files are read through their own era's tables) | Yes, for looking at a shot |
| `echelle-calib [FOLDER]` | Separate live calibration bench; drag SIFs onto it (FOLDER sets where the file dialog opens, where the two output folders are derived, and the acquisition date the snapshot identity is prefilled from; `--watch` adds folder polling), triages every exposure, takes hand-assigned roles for any lamp, fits alignment, compares sphere response, and builds configuration/snapshot evidence | Yes, during calibration |
| `echelle --version` | Prints the installed version, so a report says which build produced it | Yes; harmless anywhere |
| `echelle status` | Summarizes snapshots, registry presence, and durable processing receipts | Yes; safest first command |
| `echelle snapshot create` | Copies role-named calibration inputs into one immutable, digested snapshot | Yes, usually through the bench |
| `echelle snapshot validate DIR` | Rechecks snapshot schema, paths, sizes, and SHA-256 digests | Yes |
| `echelle snapshot show DIR` | Prints a compact snapshot summary | Yes |
| `echelle snapshot import-historical ID --calibrations DIR` | Converts a bundled 2019/2024/2025 binder into a registrable snapshot; `--artifact-root` supplies campaign folders too large to package | Once per historical epoch |
| `echelle process INPUT -o OUTPUT` | Converts one SIF, a folder, or several drives to SpectroCube NetCDF and records resumable receipts; a plain filename pattern (default `*.SIF`) walks the whole source tree, pruning `calibrations/` folders and any folder holding `snapshot.toml` and naming every pruned folder in the header and the receipt; `--dry-run` names every file it would convert | Yes; a registry run needs `--sample N`/`--sample auto` or `--drift-verdict` |
| `echelle catalog build/merge` | Writes per-drive catalogs keyed on stable drive ids and merges them by recency into an all-years index; `--receipt-dir` takes the runs root or one receipt folder, and is what supplies the drive id and the gate | Candidate; audit/catalog work |
| `echelle txt CUBE OUTPUT` / `echelle-cube2txt` | Writes LHD text at the frozen legacy header; refuses a cube missing `trigger_delay_s`, `frame_interval_s` or `exposure_s`, so the run that wrote the cube must have carried timing (see [Getting to `echelle txt`](#getting-to-echelle-txt)) | Candidate; no raw SIF needed |
| `echelle recal-cube CUBE --new-snapshot DIR` | Applies safe wavelength/factor snapshot deltas and refuses geometry changes | Candidate; reviewed repair only |
| `echelle drift audit/refine` | Samples Balmer/Fulcher centroids, solves one rigid detector shift in pixels, emits a four-state verdict, and accepts immutable `-rN` refinements. `-o` is optional — omitted, the evidence lands beside the audited cubes as the next free `drift-evidence-NNN.json`; `--every` is optional too, derived so roughly 20 cubes get measured; with `--from`/`--to` the derivation counts only the cubes inside the window, because the audit filters by date before it samples. Each Balmer window is judged against both the H and D references and assigned to the nearer one, so a deuterium shot is tagged `D` and read as aligned instead of condemned by the 0.178 nm isotope offset; a `D` shot drops the H2 Fulcher anchors, and a tag that disagrees with the bundled LHD deuterium calendar is flagged in the evidence, never silently resolved | Candidate; required before any registry run |
| `echelle web [--open]` | Builds the static read-only campaign page and always prints its absolute `index.html` path; `--open` also opens it in the default browser. Four tabs: Now (a per-drive stepper saying what to do first), Drives (catalog with a local Find), Calibration (one record per saved snapshot — where its folder is, which lamps and spheres went into it, and the bench numbers: dx/dy/theta/RMS, lines used and worst residual, the vetted lineage, and the sphere comparison against its named reference — then epochs and drift verdicts) and the packaged reading room. `--home DIR-or-campaign.toml`, or a `campaign.toml` in the current directory, supplies `--catalog`/`--output`/`--registry`/`--calibrations`/`--drift` defaults that explicit flags always override; `--registry` and `--drift` are what let the stepper place the campaign | Candidate; never controls workers |
| `echelle historical` | Validates the three thin historical calibration binders | Candidate; inspection only |
| `echelle-pattern SPHERE BACKGROUND` | Previews or writes a detector order-pattern fit | Specialist recalibration |
| `echelle-align SIGNAL BACKGROUND SPHERE SPHERE_BG` | Previews or saves a rigid wavelength-table alignment | Specialist recalibration |
| `echelle-validate-lines SIF` | Checks a calibrated experimental spectrum against Balmer and optional Fulcher lines | Calibration QC |
| `echelle-wavelength-qc` | Reproduces per-order wavelength fits and writes QC plots/tables | Calibration QC |
| `echelle-nist-overlay SIF` | Draws cached NIST lamp sticks over extracted order spectra | Lamp identification/QC |
| `echelle-spectrocube` | Compatibility name for the exporter behind `echelle process` | Prefer `echelle process` |

Run `COMMAND --help` for the complete option list. In a portable kit, replace
`COMMAND` with the kit-local form shown above.

## Recommended NIFS trip loop

Keep acquired SIF files, generated calibrations, configuration, cubes, and run
receipts outside the kit folder. Then the kit can be replaced without risking
campaign data.

### A drive this machine can read but not write

A Windows-formatted (NTFS) drive mounts read-only on a Mac. That is fine for
the data — raw SIFs are only ever read — but the campaign home cannot live on
that drive there. Point the served page (`echelle serve`, or `echelle web
--open` from no campaign folder) at the drive anyway: it detects the read-only
volume, explains it, and offers to create a writable home (under
`~/Echelle-campaigns/` by default, or any folder you pick) whose
`campaign.toml` remembers the drive as its `data` folder. Cubes and drift
evidence then derive into that home instead of beside the data; the Advanced
fold can re-point them. On a machine that can write the drive — the same drive
back on Windows — nothing changes: everything derives beside the data as
before.

### Windows PowerShell

```powershell
$Data = "D:\NIFS"

# 1. Confirm the installation and current evidence.
.\echelle.ps1 status `
  --calibrations "$Data\calibrations" `
  --registry "$Data\calibration_registry.toml" `
  --runs "$Data\runs"

# 2. Open the live bench on the calibration folder, then drag SIFs onto it.
#    The two roots below are only needed because this loop keeps them one
#    level up; without them the bench writes everything it generates into
#    calibrations\ inside the folder argument itself.
& ".\.venv\Scripts\echelle-calib.exe" "$Data\20250926_calib" `
  --output-root "$Data\calibrations" `
  --config-root "$Data\calibrations\configs"

# 3. Recheck the snapshot produced by the completed bench procedure.
.\echelle.ps1 snapshot validate "$Data\calibrations\20260814_cmos"

# 4. Preview the first sample without writing files or receipts.
#    --sample auto derives its size from what the folder holds
#    (max(5, min(30, ceil(files/25))), or every file when fewer than 5).
.\echelle.ps1 process "$Data\shots" `
  -o "$Data\cubes" `
  --runs-dir "$Data\runs" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --volume-label NIFS-A `
  --sample auto `
  --dry-run

# 5. Take the sample for real by repeating step 4 without --dry-run.
#    Nothing else may run yet: these cubes are the audit's input.

# 6. Audit the sample into one immutable verdict file. The audit takes the
#    cube folder itself, so no shell glob is involved. -o is optional: left
#    out, the evidence is written beside the cubes as the next free
#    drift-evidence-NNN.json, numbered from 001. --every is optional too,
#    derived as max(1, cubes // 20) so roughly 20 cubes get measured.
.\echelle.ps1 drift audit "$Data\cubes" `
  --calibrations "$Data\calibrations"

# 7. Process the whole drive under that verdict. Name the evidence file the
#    audit just printed, e.g. drift-evidence-001.json beside the cubes.
#    --config carries the trigger timing no calibration measures; without it
#    step 9 cannot write LHD text. See "Getting to echelle txt" below.
.\echelle.ps1 process "$Data\shots" `
  -o "$Data\cubes" `
  --runs-dir "$Data\runs" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --config "$Data\export-timing.toml" `
  --volume-label NIFS-A `
  --drift-verdict "$Data\cubes\drift-evidence-001.json"

# 8. Index what this drive now holds, and fold it into the all-years index.
#    --receipt-dir takes the runs root itself: the receipt for these cubes is
#    found inside it, and with it the drive's stable id and the run's gate.
.\echelle.ps1 catalog build "$Data\cubes" `
  --volume-label NIFS-A `
  --receipt-dir "$Data\runs" `
  --output "$Data\cubes\catalog.json"

.\echelle.ps1 catalog merge "$Data\cubes\catalog.json" `
  -o "$Data\all-years.json"

# 9. Write the LHD deliverable for one shot. Repeat per cube as needed.
.\echelle.ps1 txt "$Data\cubes\193778_Echelle_spectrocube.nc" `
  "$Data\text\193778.txt"

# 10. Build and open the campaign page in one step.
.\echelle.ps1 web `
  --catalog "$Data\all-years.json" `
  --output "$Data\campaign-page" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --drift "$Data\cubes\drift-evidence-001.json" `
  --open
```

### macOS Terminal

```bash
data="/Volumes/NIFS"

# 1. Confirm the installation and current evidence.
./echelle status \
  --calibrations "$data/calibrations" \
  --registry "$data/calibration_registry.toml" \
  --runs "$data/runs"

# 2. Open the live bench on the calibration folder, then drag SIFs onto it.
#    The two roots below are only needed because this loop keeps them one
#    level up; without them the bench writes everything it generates into
#    calibrations/ inside the folder argument itself.
./.venv/bin/echelle-calib "$data/20250926_calib" \
  --output-root "$data/calibrations" \
  --config-root "$data/calibrations/configs"

# 3. Recheck the completed snapshot.
./echelle snapshot validate "$data/calibrations/20260814_cmos"

# 4. Preview the first sample before writing cubes or receipts.
#    --sample auto derives its size from what the folder holds
#    (max(5, min(30, ceil(files/25))), or every file when fewer than 5).
./echelle process "$data/shots" \
  -o "$data/cubes" \
  --runs-dir "$data/runs" \
  --registry "$data/calibration_registry.toml" \
  --calibrations "$data/calibrations" \
  --volume-label NIFS-A \
  --sample auto \
  --dry-run

# 5. Take the sample for real by repeating step 4 without --dry-run.
#    Nothing else may run yet: these cubes are the audit's input.

# 6. Audit the sample into one immutable verdict file. The audit takes the
#    cube folder itself, so no shell glob is involved. -o is optional: left
#    out, the evidence is written beside the cubes as the next free
#    drift-evidence-NNN.json, numbered from 001. --every is optional too,
#    derived as max(1, cubes // 20) so roughly 20 cubes get measured.
./echelle drift audit "$data/cubes" \
  --calibrations "$data/calibrations"

# 7. Process the whole drive under that verdict, in ABSOLUTE units. The LHD
#    deliverable must be absolutely calibrated, so --units wmsr applies the
#    snapshot's sphere factor (radiance, W/m2/nm/sr) instead of raw counts, and
#    --calibration-source names the sphere standard this drive traces to (each
#    campaign has its own sphere; name that drive's). Name the evidence file the
#    audit just printed, e.g. drift-evidence-001.json beside the cubes. --config
#    carries the trigger timing no calibration measures; without it step 9
#    cannot write LHD text. See "Getting to echelle txt" below.
./echelle process "$data/shots" \
  -o "$data/cubes" \
  --runs-dir "$data/runs" \
  --registry "$data/calibration_registry.toml" \
  --calibrations "$data/calibrations" \
  --config "$data/export-timing.toml" \
  --units wmsr \
  --calibration-source "snapshot 20250926_cmos integrating sphere" \
  --volume-label NIFS-A \
  --drift-verdict "$data/cubes/drift-evidence-001.json"

# 8. Index what this drive now holds, and fold it into the all-years index.
#    --receipt-dir takes the runs root itself: the receipt for these cubes is
#    found inside it, and with it the drive's stable id and the run's gate.
./echelle catalog build "$data/cubes" \
  --volume-label NIFS-A \
  --receipt-dir "$data/runs" \
  --output "$data/cubes/catalog.json"

./echelle catalog merge "$data/cubes/catalog.json" \
  -o "$data/all-years.json"

# 9. Write the LHD deliverable for one shot. Repeat per cube as needed.
./echelle txt "$data/cubes/193778_Echelle_spectrocube.nc" \
  "$data/text/193778.txt"

# 10. Build and open the campaign page in one step.
./echelle web \
  --catalog "$data/all-years.json" \
  --output "$data/campaign-page" \
  --registry "$data/calibration_registry.toml" \
  --calibrations "$data/calibrations" \
  --drift "$data/cubes/drift-evidence-001.json" \
  --open
```

Put every accepted snapshot ID into the ordered registry and verify its
inclusive bounds with `echelle status` before processing. Use one
`--volume-label` for each input folder when processing several drives.

## Getting to `echelle txt`

`echelle txt` writes the frozen LHD header, which needs three timing numbers.
Two of them — `frame_interval_s` and `exposure_s` — come out of the detector
header and are recorded for you. The third, `trigger_delay_s`, is measured by
nothing: no calibration produces it, no snapshot carries it, and no processing
flag supplies it. It reaches a cube only through an export config's `[metadata]`
block. Without that config, step 9 refuses every cube the loop just wrote.

So the loop's step 7 names one small file, and this is the whole of it. Save it
as `D:\NIFS\export-timing.toml` (`/Volumes/NIFS/export-timing.toml`):

```toml
# Minimal export config: timing only, safe to use with --registry.
# It deliberately names no camera, calibration folder or calibration file --
# the registry owns all of those, and a config that also named them would be
# refused as two authorities for one run.

[metadata]
# Seconds between the LHD discharge trigger and the first frame. This is the
# one number the calibration bench cannot measure; take it from the machine's
# timing, or carry forward the value the previous campaign used.
trigger_delay_s = 2.5
time_axis_reference = "LHD discharge time"
frame_time_formula = "trigger_delay_s + frame * frame_interval_s"
```

The live bench already writes these three lines: they are the `[metadata]` block
at the top of `calibrations\configs\<snapshot-id>\export.toml`, carried forward
from the previous campaign and marked as inherited. That file cannot be passed
to a registry-backed run as it stands, because it also names the camera and the
calibration files — copy its `[metadata]` block into `export-timing.toml` and
edit the delay if this campaign's timing differs. A run **without** `--registry`
can use the bench's `export.toml` whole, with `--config`.

### `--calibration-source` on absolute runs

The trip loop above already runs absolute (`--units wmsr`), because the LHD
deliverable must be — `echelle txt` refuses a counts cube, so a drive processed
in counts cannot become LHD text at all. The absolute factor itself comes from
the registry-resolved snapshot's own sphere; `--units wmsr` only tells the run
to apply it and write radiance (`W/m2/nm/sr`) instead of raw counts.

`--calibration-source` says what that scale was traced to:

```powershell
--units wmsr --calibration-source "snapshot 20250926_cmos integrating sphere"
```

It is free text and it becomes the cube's `calibration_source` attribute — the
only record of which sphere standard stands behind the numbers, and it is
carried into the LHD text's own comments. It is one value for the whole run, so
name the sphere for the drive being processed; a single `echelle process` call
whose shots span two calibration eras would label them all alike even though
the factors are still applied per snapshot. An absolute run that omits it warns
and then writes cubes carrying no such record at all. The bench's generated
`export.toml` sets it under `[export]` alongside `units`, so a non-registry run
driven by that file already has it.

The absolute scale stays revisable after the fact: `echelle recal-cube
--factor-only` re-applies a different era's sphere factor onto a saved cube
without reopening the raw SIF, so a cube can be re-scaled if a calibration is
later refined.

Steps 4–7 are the whole gate: a registry-backed run needs either `--sample N`
or `--sample auto`, which processes at most N (or the derived count) resolved
files and marks its receipt and cubes as an unverified sample, or
`--drift-verdict`, which spends the sampled evidence the audit wrote. A run
with an explicit `--config` calibration and no registry stays legal and is
recorded as `ungated (no registry)`.

## Where the bench writes

`echelle-calib` derives both of its output folders from the folder argument it
was launched at. Launch it at `D:\NIFS\20250926_calib` and it writes:

```text
D:\NIFS\20250926_calib\calibrations\<snapshot-id>\         the immutable snapshot
D:\NIFS\20250926_calib\calibrations\configs\<snapshot-id>\ campaign/alignment/export settings
```

Inside, because that folder is the calibration: the measured lamp and sphere
frames and the calibration computed from them, side by side in the one folder
that holds it all. That is also the fix for a bench started from a shortcut
scattering its output into whatever directory the shortcut happened to begin
in. Everything the bench generates then lives under one tidy `calibrations\`
folder — the snapshot identities, plus a `configs\` subfolder for the settings
bundles — so a calibration folder gains one generated child, not two. To put
them somewhere else — the trip loop above keeps them one level up, beside the
shots — `--output-root` and `--config-root` override each independently.

The snapshot identity the Save tab prefills is `YYYYMMDD_<detector>` dated by
**acquisition**: the `YYYYMMDD` leading the launch folder's name, else the date
in the loaded SIFs' headers, else today with the Save tab saying it could not
be derived and must be checked. `--snapshot-id`, and your own typing in the ID
field, always win.

Either way the Save tab states both paths **in full**, so you can read where
the files will go before pressing anything, and once a snapshot is saved an
**Open folder** button appears beside the confirmation and shows it in the file
manager. No path has to be retyped, which matters most on a network share.

What each of those files is for, and which of them a cube is actually built
from, is [From calibration to cube](calibration-to-cube.md).

## The campaign page

`echelle web` writes one self-contained `index.html` and always prints its
absolute path. There is no server and nothing is fetched: **double-click the
file** to open it in a browser, or pass `--open` and let the command do that
for you.

The simplest form reads a hand-editable `campaign.toml` beside the campaign —
in the current directory, or wherever `--home DIR-or-campaign.toml` points —
for whichever of `--catalog`/`--output`/`--registry`/`--calibrations`/`--drift`
are missing from the command line; paths inside it resolve against the file's
own folder, and any explicit flag still wins:

```powershell
.\echelle.ps1 web --open
```

Without a `campaign.toml`, name everything explicitly:

```powershell
.\echelle.ps1 web `
  --catalog "$Data\all-years.json" `
  --output "$Data\campaign-page" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --drift "$Data\epoch-drift.json" `
  --open
```

- `--catalog` is required (from the command line or `campaign.toml`): the
  merged all-years index from `echelle catalog merge`, or a single drive's
  catalog.
- `--output` is required (from the command line or `campaign.toml`): the
  folder `index.html` is written into.
- `--registry` and `--calibrations` name the epochs, and `--drift` supplies the
  verdicts; together they are what let the Now tab place the campaign and say
  what to do next. Without them the page still builds, with less to say.
- `--document` appends one extra Markdown file after the packaged reading room.

The page's own composer asks for exactly two things — the data folder and the
calibration epoch — and derives the rest, with the derived values editable
under an "Advanced" fold; see [the campaign page
composer](harbor-candidate.md#the-campaign-page).

The page is a snapshot of the moment it was built. After processing more shots,
rebuild the catalog and rerun `echelle web` over the same `--output` folder to
refresh it. The page never controls a worker and never writes campaign data.

Prefer clicking to pasting? `echelle serve` (from the campaign folder, or
`--home` pointing at one) serves the same page at a local address it prints,
with a Browse… folder picker on the data-folder field; a cold start with no
campaign home lets you pick the campaign folder in the page, which then
writes `campaign.toml` for you. Loopback only; the page composes and checks,
the terminal runs the work. Add `notes = "<folder>"` to `campaign.toml` — a
relative or absolute path — and the served page links that folder's own pages
in its header and serves them at `/notes/<name>`; without the key there are
none, and pointing several campaigns at one folder keeps the pages in a single
place, on your machine and out of this repository.

Learn the page with nothing: `echelle web --practice --open` — an invented campaign
built to a temp folder, labels say PRACTICE.

## Safe inspection versus writing

These commands do not intentionally create campaign artifacts:

```text
echelle --help
echelle --version
echelle status
echelle snapshot show DIR
echelle snapshot validate DIR
echelle process INPUT --dry-run
echelle catalog
echelle drift
echelle-pattern SPHERE BACKGROUND
echelle-align SIGNAL BACKGROUND SPHERE SPHERE_BG
echelle-nist-overlay --list-lamps
```

Bare `echelle snapshot`, `echelle catalog`, and `echelle drift` print their own
help instead of failing, so a subcommand name never has to be remembered.

Writing begins with `snapshot create`, a real `process` run, `echelle-pattern
--output`, `echelle-align --save`, wavelength-QC/overlay output commands, or the
bench's explicit generate/save controls. Existing snapshots and specialist
outputs are refused unless the relevant command explicitly receives its
overwrite option.

## Common command failures

### “Command not found” or “not recognized”

- In a portable kit, bare commands are not placed on `PATH`. Use
  `.\echelle.ps1`/`./echelle` or the explicit `.venv` command paths above.
- In a development checkout, start with `lab activate echelle-spectra`. If an
  older command works but a newly shipped command does not, rerun
  `uv pip install --no-deps -e .`. Editable Python source changes appear
  immediately, but a newly declared console executable is created only when
  the install metadata is refreshed. The Lab-managed environment is uv-native
  and may intentionally have no `pip` module.
- `echelle_spectra` contains an underscore. `echelle-spectra` is not a command.

As a development fallback, every command can be reached as a module. For the
new bench:

```powershell
python -m echelle_spectra.calibration_bench_gui --help
```

### A copied multiline command fails

PowerShell continues a line with a backtick (`` ` ``); macOS shells use a
backslash (`\`). Do not copy the continuation character from the other shell.
Always quote paths containing spaces.

### Status says no snapshots or registry

That is an honest empty state, not an installation failure. Supply explicit
`--calibrations`, `--registry`, and `--runs` paths when the campaign evidence
is not below the current working directory. Status exits nonzero when it finds
an invalid snapshot or a registry whose referenced snapshots, bounds, or
digests do not validate.

### Processing refuses a registry run

```text
ERROR: registry-backed processing requires a sampled epoch audit
(--drift-verdict) or an explicit unverified first sample (--sample N).
  Take the sampled evidence this gate needs:
    echelle process "D:\NIFS\shots" --registry "D:\NIFS\calibration_registry.toml" ...
    echelle drift audit "D:\NIFS\cubes" -o drift-evidence.json
  Then repeat this command with --drift-verdict drift-evidence.json
```

This is the calibration gate, not a broken installation. No drive is processed
under an epoch nothing has checked. The refusal prints the exact sample and
audit commands that produce the missing evidence; run them, then repeat the
processing command with `--drift-verdict`. Every registry-backed form is gated,
including a single file and a folder holding one file.

Two related refusals:

- `registry still selects 20250926_cmos; the accepted correction produced
  20250926_cmos-r1` — the audit's shift was accepted, so only the refinement is
  authorized. Repoint the registry entry at the `-rN` snapshot.
- `bulk processing refused: sampled verdict is insufficient-data` — the sample
  measured too few usable lines. A verdict needs at least three resolved lines
  (blends excluded) from three different Echelle orders, with at least one in
  each half of the audited range. Sample more or different shots; insufficient
  evidence is never read as aligned.

### Processing finds no files

Folder processing defaults to `*.SIF` and, because that is a plain filename
pattern rather than one written with `/` or `**`, walks the whole source tree
under the folder you named — real campaign drives are date-named day folders,
so shots several levels down are still found. It prunes any folder named
`calibrations` and any folder holding `snapshot.toml`, so lamp frames are
never swept up as science shots; junctions and symbolic links are matched
where they sit but never descended into. For lowercase acquisition filenames,
add `--pattern "*.sif"` (also tried automatically as a fallback). Start with
`--dry-run`, which prints a `Would convert:` block naming every file it found
relative to the source root, and inspect that list before allowing writes.

### Fewer files than the drive holds

Read the `Skipped:` line in the batch header. Every run that pruned a folder
names it there — count first, then the folders relative to the source root,
bounded to five plus a count of the rest — and records the whole list under
`pruned_dirs` in the run receipt's `run.toml`. The usual cause is a bench
launched at a day folder that also holds shots: the `snapshot.toml` it saved
there marks the whole day as calibration. Move or rename that file, or name a
path pattern with `--pattern`, and rerun.

### The kit installer refuses to continue

- Wrong-platform and wrong-architecture refusal means the other kit payload is
  required; do not bypass it.
- An incomplete `.runtime` or `.venv` may be removed by the exact name reported
  by the installer, then regenerated.
- A checksum error means the payload must be recopied or rebuilt. Do not remove
  or edit `checksums.sha256`, `runtime`, or `wheelhouse` to force an install.

## Version check

Ask the command itself:

```powershell
.\echelle.ps1 --version
```

The fallbacks below answer the same question without the umbrella command, for
a kit whose install has not finished or an environment where only the library
imports. Inside the active source environment:

```bash
python -c "import echelle_spectra; print(echelle_spectra.__version__)"
```

For an unactivated portable Windows kit:

```powershell
& ".\.venv\Scripts\python.exe" -c "import echelle_spectra; print(echelle_spectra.__version__)"
```

For an unactivated portable macOS kit:

```bash
./.venv/bin/python -c 'import echelle_spectra; print(echelle_spectra.__version__)'
```
