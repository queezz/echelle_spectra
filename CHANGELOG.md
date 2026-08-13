# Changelog

User-visible release history for Echelle Spectra, newest first. Versioning
follows the Fleet convention from 0.3.0 onward: a substantial capability earns
a minor release, a compatible correction earns a patch, and documentation,
tests, or internal refactoring alone do not move the number.

## 1.0.0 — 2026-08-13

**Campaign processing becomes dependable automation.** Every non-dry folder
run writes an atomic `run.toml` summary and append-only `records.jsonl` ledger.
Each source attempt records content identity, output, snapshot ID, volume label,
terminal status, reason, and timing; completed outputs also carry size and
SHA-256 identity.

**Interruptions are safe and resumable.** SpectroCubes are published through a
temporary sibling and atomic replace, Ctrl-C records an interrupted run, and
rerunning the same source/destination resumes automatically. Previously
completed work is skipped only after both source and output still match their
recorded sizes and digests. Ordinary source failures remain isolated and the
batch accounts for later files before exiting nonzero.

**Progress and status come from evidence.** Batch output reports count, measured
rate, and ETA. `echelle status --runs` reads durable receipts to show the latest
run, its state, result counts, and selected calibration snapshot instead of
guessing from output filenames.

The installed CLI now keeps help/progress text compatible with the ordinary
Windows console and reads packaged GUI defaults without trying to create a
configuration file inside the installed package.

Version 1 begins the deliberate campaign-automation epoch. It does not imply
that every later live-bench, registry, recalibration, drift, and catalog feature
is already present; those continue as compatible 1.x releases.

## 0.3.0 — 2026-08-13

**Calibration snapshots become first-class artifacts.** `echelle snapshot
create` assembles role-named calibration inputs under one dated snapshot ID,
records lamp species and source filenames, SHA-256 digests, byte sizes,
validity, alignment/QC summaries, and an optional base-snapshot provenance
chain. Construction is atomic and refuses to replace an existing snapshot.
`validate` rechecks paths, required roles, sizes, and content digests; `show`
provides a compact human summary.

**The campaign gets one front door.** The new `echelle` command explains the
workflow when run bare, reports the honest current snapshot/registry state with
`status`, exposes snapshot commands, and routes `process` to the existing
SpectroCube batch exporter. Existing `echelle-*` commands remain compatible.

**The release contract is now enforced.** Package and manifest versions must
match this newest changelog heading, the complete test suite is the commit
gate, and Ruff now uses the Fleet science-package policy. Pre-0.3.0 numerical
modules retain narrow grandfathered complexity/modernization exceptions while
new campaign code receives the full rule set.

## 0.2.0 — 2026-06-25

- Propagated export and calibration metadata into generated SpectroCubes,
  including calibration file provenance and order coverage.

## 0.1.0 — 2026-06-15

- Established the first versioned release after adding headless order-pattern,
  wavelength alignment/QC, NIST overlay, and SpectroCube export workflows.

Earlier history remains available in Git; 0.3.0 is where the current Fleet
commit, version, and changelog contract begins.
