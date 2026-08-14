# Echelle operator cheat sheet

This is the short manual for the installed surfaces in Echelle Spectra 1.6.
It separates the portable NIFS kit from a source checkout because their command
prefixes are different. The portable kit does **not** add commands to the
machine-wide `PATH`.

Lab is optional development convenience on the owner's own machines. It is not
an Echelle dependency and is not required on an institute workstation. For the
trip, use the portable-kit section below. Administrator access may be available
there, but the kit does not need to use it or change the machine-wide Python.

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

The main GUI, live bench, and specialist commands live inside the kit-local
environment:

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
| `echelle_spectra` | Main five-tab analysis GUI for loading, calibrating, viewing, fitting, and exporting spectra | Yes |
| `echelle-calib FOLDER` | Separate live calibration bench; watches stable SIF arrivals, classifies measurements, fits alignment, compares sphere response, and builds configuration/snapshot evidence | Yes, during calibration |
| `echelle status` | Summarizes snapshots, registry presence, and durable processing receipts | Yes; safest first command |
| `echelle snapshot create` | Copies role-named calibration inputs into one immutable, digested snapshot | Yes, usually through the bench |
| `echelle snapshot validate DIR` | Rechecks snapshot schema, paths, sizes, and SHA-256 digests | Yes |
| `echelle snapshot show DIR` | Prints a compact snapshot summary | Yes |
| `echelle process INPUT -o OUTPUT` | Converts one SIF, a folder, or several drives to SpectroCube NetCDF and records resumable receipts | Yes; a registry run needs `--sample N` or `--drift-verdict` |
| `echelle catalog build/merge` | Writes per-drive catalogs and a durable all-years index with volume labels | Candidate; audit/catalog work |
| `echelle txt CUBE OUTPUT` / `echelle-cube2txt` | Writes canonical provenance-complete LHD text from a saved cube | Candidate; no raw SIF needed |
| `echelle recal-cube CUBE --new-snapshot DIR` | Applies safe wavelength/factor snapshot deltas and refuses geometry changes | Candidate; reviewed repair only |
| `echelle drift audit/refine` | Samples Balmer/Fulcher centroids, emits a four-state verdict, and accepts immutable `-rN` refinements | Candidate; required before any registry run |
| `echelle web --catalog INDEX --output DIR` | Builds the static read-only reading room and command composer | Candidate; never controls workers |
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

### Windows PowerShell

```powershell
$Data = "D:\NIFS"

# 1. Confirm the installation and current evidence.
.\echelle.ps1 status `
  --calibrations "$Data\calibrations" `
  --registry "$Data\calibration_registry.toml" `
  --runs "$Data\runs"

# 2. Watch the acquisition folder in the live bench.
& ".\.venv\Scripts\echelle-calib.exe" "$Data\incoming" `
  --output-root "$Data\calibrations" `
  --config-root "$Data\calibration-configs"

# 3. Recheck the snapshot produced by the completed bench procedure.
.\echelle.ps1 snapshot validate "$Data\calibrations\20260814_cmos"

# 4. Preview the first sample without writing files or receipts.
.\echelle.ps1 process "$Data\shots" `
  -o "$Data\cubes" `
  --runs-dir "$Data\runs" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --volume-label NIFS-A `
  --sample 5 `
  --dry-run

# 5. Take the sample for real by repeating step 4 without --dry-run.
#    Nothing else may run yet: these cubes are the audit's input.

# 6. Audit the sample into one immutable verdict file.
.\echelle.ps1 drift audit (Get-ChildItem "$Data\cubes\*.nc").FullName `
  -o "$Data\epoch-drift.json"

# 7. Process the whole drive under that verdict.
.\echelle.ps1 process "$Data\shots" `
  -o "$Data\cubes" `
  --runs-dir "$Data\runs" `
  --registry "$Data\calibration_registry.toml" `
  --calibrations "$Data\calibrations" `
  --volume-label NIFS-A `
  --drift-verdict "$Data\epoch-drift.json"
```

### macOS Terminal

```bash
data="/Volumes/NIFS"

# 1. Confirm the installation and current evidence.
./echelle status \
  --calibrations "$data/calibrations" \
  --registry "$data/calibration_registry.toml" \
  --runs "$data/runs"

# 2. Watch the acquisition folder in the live bench.
./.venv/bin/echelle-calib "$data/incoming" \
  --output-root "$data/calibrations" \
  --config-root "$data/calibration-configs"

# 3. Recheck the completed snapshot.
./echelle snapshot validate "$data/calibrations/20260814_cmos"

# 4. Preview the first sample before writing cubes or receipts.
./echelle process "$data/shots" \
  -o "$data/cubes" \
  --runs-dir "$data/runs" \
  --registry "$data/calibration_registry.toml" \
  --calibrations "$data/calibrations" \
  --volume-label NIFS-A \
  --sample 5 \
  --dry-run

# 5. Take the sample for real by repeating step 4 without --dry-run.
#    Nothing else may run yet: these cubes are the audit's input.

# 6. Audit the sample into one immutable verdict file.
./echelle drift audit "$data"/cubes/*.nc -o "$data/epoch-drift.json"

# 7. Process the whole drive under that verdict.
./echelle process "$data/shots" \
  -o "$data/cubes" \
  --runs-dir "$data/runs" \
  --registry "$data/calibration_registry.toml" \
  --calibrations "$data/calibrations" \
  --volume-label NIFS-A \
  --drift-verdict "$data/epoch-drift.json"
```

Put every accepted snapshot ID into the ordered registry and verify its
inclusive bounds with `echelle status` before processing. Use one
`--volume-label` for each input folder when processing several drives.

Steps 4–7 are the whole gate: a registry-backed run needs either `--sample N`,
which processes at most N resolved files and marks its receipt and cubes as an
unverified sample, or `--drift-verdict`, which spends the sampled evidence the
audit wrote. A run with an explicit `--config` calibration and no registry stays
legal and is recorded as `ungated (no registry)`.

## Safe inspection versus writing

These commands do not intentionally create campaign artifacts:

```text
echelle --help
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
    echelle drift audit "D:\NIFS\cubes\*.nc" -o drift-evidence.json
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
  measured too few usable lines. Sample more or different shots; insufficient
  evidence is never read as aligned.

### Processing finds no files

Folder processing defaults to `*.SIF`. For lowercase acquisition filenames,
add `--pattern "*.sif"`. Start with `--dry-run` and inspect the listed sources
before allowing writes.

### The kit installer refuses to continue

- Wrong-platform and wrong-architecture refusal means the other kit payload is
  required; do not bypass it.
- An incomplete `.runtime` or `.venv` may be removed by the exact name reported
  by the installer, then regenerated.
- A checksum error means the payload must be recopied or rebuilt. Do not remove
  or edit `checksums.sha256`, `runtime`, or `wheelhouse` to force an install.

## Version check

Run this inside the active source environment:

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
