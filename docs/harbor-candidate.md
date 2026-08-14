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

Both commands, legacy spectrum saves, and GUI intensity saves use the same
`echelle-lhd-text/v1` writer. Its header records cube identity, snapshot and
registry identity, calibration digests, polynomial/factor interpretation, and
recalibration history where available.

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
echelle drift audit cubes/*.nc --every 20 --shot 193778 -o epoch-drift.json
echelle drift refine epoch-drift.json --calibrations calibrations --accept-shift 0.0842
```

The audit fits baseline-subtracted Balmer and Fulcher centroids and retains
per-cube/per-line status, centroid, residual, SNR, source reference, and sample
rule. Thresholds are serialized with every verdict:

| Verdict | Rule |
| --- | --- |
| `aligned` | at least two centroids; median absolute shift ≤ 0.03 nm and every absolute residual ≤ 0.06 nm |
| `shifted` | median absolute shift ≤ 0.25 nm and every measured shift lies within 0.04 nm of the median |
| `misaligned-beyond-repair` | sufficient centroids but the aligned/rigid-shift rules fail |
| `insufficient-data` | fewer than two centroids at SNR ≥ 4 |

`insufficient-data` is never treated as aligned. A shifted verdict prints a
paste-ready refinement command. Exact acknowledgement creates the next
immutable `-rN` snapshot, adjusts only its copied wavelength table, records the
base snapshot/evidence digest, and writes a separate accepted verdict.

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
