# Changelog

User-visible release history for Echelle Spectra, newest first. Versioning
follows the Fleet convention from 0.3.0 onward: a substantial capability earns
a minor release, a compatible correction earns a patch, and documentation,
tests, or internal refactoring alone do not move the number.

## 1.7.0 — 2026-08-28

This release closes the held campaign-automation era: everything since 1.6.0
shipped under a deliberate version hold while the multi-year LHD conversion
campaign was being built and run, and lands here as one release — the version
that converts the archive's remaining drives.

**Every cube now carries its noise provenance and its maker.** The dark
subtraction in `Spectrum` destroyed the total-counts information a per-pixel
shot-noise model needs; the exporter now writes the subtracted per-pixel dark
level as a wavelength-aligned `background_counts` variable, and
`reconstruct_counts()` recovers net and total detector counts from a cube
alone. Cubes also record `source_package_version`, so a campaign spanning two
code states names each cube's maker. A cube whose `background_frames`
attribute names subtracted frames but which lacks the variable predates this
release.

**Campaigns run from read-only archive drives at scale.** A campaign home can
live beside a read-only NTFS data drive; batch discovery is recursive with
calibration folders pruned and AppleDouble siblings skipped; drive and run
identity cross the input/output seam into receipts and catalogs; resume cost
follows the work remaining; the registry loads snapshots on owned artifacts;
and the field verbs (`echelle inventory`, drift geometry survey) give the
package its own eyes on an unknown drive.

**The drift verdict earned physics discipline.** The auditor reads the plasma
background off the dark frames, decides the verdict on lines strong enough to
carry it, reports weak scatter beside rather than inside the decision, and
raises a loud `era-misassigned-calibration` verdict when an all-Balmer isotope
flip betrays the wrong epoch. Cubes recalibrate onto another era's snapshot
with honestly dropped factor columns counted and reported.

**The served campaign page walks the owner's own procedure.** `echelle serve`
opens real folders with a picker, writes the campaign home from its setup
page, composes every command from two inputs (data folder and calibration
epoch), launches the campaign's own verbs as detached processes with durable
receipts, renders each snapshot as a physicist's record on the Calibration
tab, and shows saved snapshots even without a registry.

**The calibration bench closes the loop at the instrument.** One press
auto-anchors to the CLI's numbers, the bench extracts a pattern from the
folder's own sphere and rebases without a restart, saves campaign-aware with
the destination named before the save, accepts historical folders as they
were shot (signal-only spheres included), and snapshots reference their
source light instead of copying it, dated by the day the light was taken.

**LHD deliverables state their calibration.** `echelle txt` refuses a
non-absolute cube, and the text records the cube's calibration source.

## 1.6.0 — 2026-08-14

**Every registry-backed export now resolves one reviewed immutable calibration
snapshot.** The ordered registry points only to snapshot IDs; the snapshots'
existing validity tables remain the single authority for inclusive shot/date
boundaries. Selection is deterministic, requires all declared identity
dimensions, validates every referenced artifact and digest before processing,
and refuses overlaps, ambiguity, missing identity, no-match sources, or any
fallback to another calibration. Batch receipts record the snapshot selected
for each source and will not resume against a changed epoch decision.

**Flagship SpectroCubes retain the complete recalibration state defined by
SpectroCube 0.2.0.** Raw zero-origin detector columns and integer echelle-order
IDs remain aligned with intensity through partial-order handling, dispersion
reversal, finite masking, crop, and the final stable wavelength sort. The
writer serializes the existing per-order `numpy.polyfit` coefficients in
descending-power order and verifies that they reconstruct every retained
wavelength sample within `5e-10 nm`.

**Absolute cubes state and prove the applied counts-per-second factor.** The
wavelength-aligned factor records `source_units = "counts/s"`, output-mapping
units, and the exact neutral application equation. Export numerically verifies
that count rate multiplied by the factor reproduces stored absolute intensity.
The selected snapshot ID, exact manifest digest and JSON, every calibration
artifact digest, registry digest and ordered position, plus the established
extraction geometry and crop provenance survive NetCDF save/load validation.

**The portable dependency contract now carries SpectroCube 0.2.0.** The public
extra, immutable Git source, uv lock, companion-wheel expectations, manifests,
installers, documentation, and installed-package gates agree on the remotely
reachable neutral contract while retaining the established Python 3.12.13
travel runtime.

## 1.5.0 — 2026-08-13

**The calibration campaign can travel as a pinned, self-contained kit.** One
manifest fixes CPython, uv, immutable companion sources, and exact platform
assets for Windows x86-64, Apple Silicon macOS, and Intel macOS. Each generated
payload includes a cached interpreter, target-specific hash-locked dependency
wheelhouse, verified project and companion wheels, narrow launchers, and a
complete checksum inventory without committing generated payloads.

**Online and offline installation share one reproducible contract.** The
installers verify the payload and native platform before extraction, keep the
runtime, environment, and uv cache inside a fresh kit copy, and refuse corrupt,
incomplete, or wrong-platform inputs. The offline route forbids package-index
access; the three-command guide reaches `echelle status` without an existing
Python installation, administrator rights, Fleet, or lab-cli.

**Release archives are reproducible outside the working checkout.** Independent
clean source exports build normalized wheel and source archives with fixed
ordering, timestamps, ownership, permissions, and gzip/ZIP metadata. Release
validation requires byte-identical artifacts, Twine checks, native Windows and
macOS execution, installed resources and entry points, both GUIs, campaign
memory, and snapshot creation and validation.

**One operator sheet now names every surface and its real invocation.** It
separates portable Windows, portable macOS, and editable-checkout command
prefixes; gives the recommended calibration-to-processing trip loop; marks
inspection versus writing actions; and explains missing entry points, shell
continuations, empty status, SIF filename matching, and installer refusal.
The sheet travels inside every checksummed kit. POSIX launchers are declared
LF-only in Git and normalized again during assembly, preventing a Windows
checkout's line-ending policy from producing an unexecutable macOS shebang.
Installer executable paths are quoted throughout, so fresh kit copies remain
usable from ordinary scratch directories whose names contain spaces.

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
