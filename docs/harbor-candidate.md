# Harbor Packets 9–13 implementation candidate

These surfaces are an **unreleased implementation candidate** for review. The
package version remains 1.6.0. No claim is made that versions 1.7.0–1.11.0 were
released or that the Harbor train shipped.

## Catalog and text

Every completed or partial batch target writes `echelle-catalog.json` beside
its cubes. The JSON records the stable volume label, last run state, relative
cube paths, cube identities, shot/year, snapshot, wavelength coverage, and any
unreadable cube error. Merge drive catalogs into a durable all-years index:

```console
echelle catalog merge /usb-a/cubes/echelle-catalog.json /usb-b/cubes/echelle-catalog.json -o all-years.json
```

The merged index retains the last catalog and volume label when a drive is
disconnected. Availability means the catalog path was reachable when the
index/reading room was read; it never implies that a missing cube is online.

Convert a saved cube without reopening raw data:

```console
echelle txt shot.nc shot.txt
echelle-cube2txt shot.nc shot.txt
```

### The text header is frozen

The LHD-side header is a contract with the machine's data system, so it stays
at its pre-unification byte shape (owner ruling, 2026-08-14). One writer
renders the two legacy dialects, both recovered verbatim into
`resources/header_template.txt` and `resources/header_template_spectrum.txt`
and pinned line for line by golden-file tests:

| Dialect | Written by | Frozen particulars |
| --- | --- | --- |
| `spec_div1` | GUI band save, `echelle txt`, `echelle-cube2txt` | `DimUnit` (singular); the LHD viewing-geometry, PFR reference, CH/D-band and contact lines in `[Comments]` |
| `spectrum` | `Spectrum.save` | `DimUnits` (plural); fixed `Name`/`DimName`; an `exposure` comment |

`ShotNo` is unquoted and `Date` is local `%m/%d/%Y %H:%M` in both. No field is
added outside the templates. Snapshot, digest and schema lines ride as free
text appended inside the existing `[Comments]` block; the bulky JSON payloads
(snapshot manifest, wavelength polynomials, recalibration history) stay in the
cube and its receipts rather than bloating the deliverable.

Cube-derived text states its timing from the cube's own `trigger_delay_s`,
`frame_interval_s` and `exposure_s` attributes. A cube missing any of the
three is refused by name — naming the attribute and the field that would have
supplied it (`[metadata] trigger_delay_s` in the export config; `CycleTime`
and `ExposureTime` recorded by `export_spectrocube`) — because a frozen header
must state its time formula and may never omit it silently. Bench-composed
export configs therefore carry trigger delay, time-axis reference, frame-time
formula and the inherited crop, so bench-born cubes can always produce text.

## Post-hoc recalibration

```console
echelle recal-cube old.nc -o revised.nc --new-snapshot calibrations/20250926_cmos-r1
```

The command evaluates the new snapshot's wavelength solution at the cube's
stored raw `detector_pixel`/`echelle_order`, recovers pre-factor signal by
dividing by the old aligned factor, and applies the new snapshot factor. Use
`--wavelength-only` or `--factor-only` for one representation. The output
retains old/new snapshot provenance and gets an adjacent immutable
`.nc.recalibration.json` manifest with input/output hashes.

The old embedded pattern digest must equal the new snapshot pattern digest.
Any change is refused with an instruction to re-extract from the read-only raw
SIF data. There is deliberately no approximation for changed geometry.

## Drift audit and refinement

```console
echelle drift audit cubes --every 20 --shot 193778 -o epoch-drift.json
echelle drift refine epoch-drift.json --calibrations calibrations --accept-shift 5.0
```

The audit fits baseline-subtracted Balmer and Fulcher centroids over the
plasma-bright frames only, then leaves wavelength for detector pixels. Each
line's wavelength residual is divided by the dispersion the cube's own stored
per-order polynomial has at that line's `detector_pixel`, so one rigid shift is
fitted in pixels and every boundary is a pixel bound. A rigid detector shift
moves a blue order by 0.0066 nm/px and H-alpha by 0.0108 nm/px; judging that in
nanometres condemned repairable shifts and under-corrected accepted ones.

| Verdict | Rule |
| --- | --- |
| `aligned` | every line's pixel residual is within 0.5 px, or within the cube's own `wavelength_accuracy_nm` where that is wider |
| `shifted` | median shift ≤ 25 px and every line lies within 1 px of that median |
| `misaligned-beyond-repair` | quorum reached but the aligned/rigid-shift rules fail |
| `insufficient-data` | fewer than three resolved lines, fewer than three orders, or all lines in one half of the audited range |

Blended catalog rows are skipped, never measured: a duplicated or
sub-resolution line can no longer manufacture a quorum. Thresholds, per-shot
medians, per-order corrections, and the skipped cubes are serialized with every
verdict. Residuals that fall into two separated per-shot groups add an
`interval_warning` naming the boundary to split at, rather than reading an
epoch step as beyond repair. `--from`/`--to` select by acquisition date and
`--shot` matches a whole shot token. A beyond-repair verdict names the drives
holding the affected shots when `--catalog` supplies the merged index.

`insufficient-data` is never treated as aligned. A shifted verdict composes the
real repair sequence — refine, repoint the registry, then `recal-cube` for the
cubes already exported. Exact acknowledgement of the pixel shift creates the
next immutable `-rN` snapshot by sliding the wavelength table's anchors along
the detector, so each order's refit turns that one shift into its own
wavelength correction; the base snapshot/evidence digest is recorded and a
separate accepted verdict is written.

Every registry-backed run — a folder of any size and the single-file path —
requires either `--drift-verdict` or an explicit unverified `--sample N` first
run. Aligned evidence authorizes the snapshots it audited; an accepted shifted
verdict authorizes only the `-rN` refinement it created; insufficient or
misaligned evidence is refused. See [Durable campaign runs](campaign-runs.md)
for the gate and the authorization each receipt records.

## Reading room

```console
echelle web --catalog all-years.json --drift epoch-drift.json --output reading-room
```

The static reading room filters years, epochs, volume labels, and run states;
shows missing-drive cards and the complete four-state drift vocabulary; drills
into per-line evidence; embeds procedure/provenance Markdown; and composes
editable plan TOML plus paste-ready terminal commands. It has no endpoint or
code path that executes a plan, starts a worker, or controls a running batch.

## Historical binders

`echelle historical` validates thin, digested manifests for 2019-05-29,
2024-03-05, and 2025-09-26. They point to the existing packaged filenames and
refuse digest/name drift; they do not rename the historical sources. The 2025
binder states explicitly that its accepted pattern/alignment inherits the
packaged 2024 sphere pair.

## Candidate limitations

The focused fixture connects snapshot provenance, a durable process receipt,
drive/merged catalog, cube text, drift/refinement, optional recalibration
domain logic, and reading-room generation. Release archives, clean/offline
installation, native platform/portable-kit rehearsals, reproducibility builds,
real historical/2026 data, full visual/perimeter review, final documentation
gates, release certification, versions, tags, and pushes remain deferred.
