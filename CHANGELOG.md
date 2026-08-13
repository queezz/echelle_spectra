# Changelog

User-visible release history for Echelle Spectra, newest first. Versioning
follows the Fleet convention from 0.3.0 onward: a substantial capability earns
a minor release, a compatible correction earns a patch, and documentation,
tests, or internal refactoring alone do not move the number.

## 1.4.0 — 2026-08-13

**The live bench now carries the calibration campaign procedure.** Stable SIFs
receive non-binding filename suggestions but only explicit sphere, background,
and lamp-role confirmation completes the self-ticking checklist. Raw-count
exposure guidance names the next safe acquisition action, and the order view
maps the shared packaged ThAr, Ne, Hg, and Fulcher H2 catalogs into a dedicated
line-identification table and overlay.

**Integrating-sphere measurements become an immediate campaign check.** The
bench runs the established absolute-calibration engine for the new sphere pair,
compares its factor curve with the previous packaged campaign, and reports
`insufficient data` as its own honest result when no defensible comparison is
available.

**A completed bench session becomes portable configuration and a validated
snapshot.** Campaign, alignment, and SpectroCube-export TOMLs are generated as
commented ordinary text from explicitly measured inputs. Snapshot save calls
the existing atomic Packet 0 creation API and then its existing validator;
replacement refusal, source immutability, failure detail, and corrected-identity
recovery remain intact.

## 1.3.0 — 2026-08-13

**A separate live calibration bench now turns stable SIF arrivals into an
interactive rigid fit.** `echelle-calib` watches an acquisition folder with an
explicit repeated-size/mtime stability rule, preserves the last good frame
through a failed load, and presents the detector with order traces beside the
selected-order spectrum.

**Calibration interaction is deterministic below Qt.** Click-guided Gaussian
centroids reuse the existing calibration rows and fitting tools; raw detector
windows produce clear/saturated verdicts before an anchor is accepted. Anchor
add, replace, remove, and clear transitions re-solve or invalidate the rigid
translation/rotation explicitly, with RMS and per-anchor residuals shown live.
The file watcher, fit state, error paths, and recovery behavior are unit-tested
without an event loop, while focused off-screen tests pin the separate GUI.

## 1.2.0 — 2026-08-13

**Known spectral lines become shared package knowledge.** One tested API now
serves provenance-carrying Balmer, Fulcher H2, ThAr, Ne, and Hg tables. The
Fulcher Q-branch anchors are bundled from the Fulcher extraction work, while
lamp rows come from the existing package-cached NIST ASD exports. Default line
validation no longer depends on a sibling `fulcheranalyzer` checkout.

**The existing GUI gains optional labeled overlays.** Independent Balmer,
Fulcher, ThAr, Ne, and Hg toggles draw wavelength markers on the main spectrum
plots. Every family starts disabled; labels thin in broad views, expand on zoom,
and disappear immediately when switched off, preserving uncluttered normal use.

## 1.1.0 — 2026-08-13

**Campaign processing now spans several drives safely.** `echelle process`
accepts several source folders and runs one sequential worker per source while
processing independent sources concurrently. Each target has an isolated output
directory, receipt tree, progress prefix, calibration instance, and failure
domain; an ordinary failure on one target does not stop the others.

**Status reconciles the campaign rather than double-counting retries.**
`echelle status` discovers nested per-target receipt trees, selects the newest
receipt for each source/output/pattern target, and reports individual plus
combined accounting. Ctrl-C requests a safe stop across workers and preserves
the established atomic-output and resumable-receipt guarantees.

Bare `echelle`, `echelle process`, `echelle snapshot`, and `echelle status`
surfaces now explain or report their next safe action, while every established
`echelle-*` entry point remains available.

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
