# Changelog

User-visible release history for Echelle Spectra, newest first. Versioning
follows the Fleet convention from 0.3.0 onward: a substantial capability earns
a minor release, a compatible correction earns a patch, and documentation,
tests, or internal refactoring alone do not move the number.

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
