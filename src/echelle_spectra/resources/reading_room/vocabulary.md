# Vocabulary

Every word this page prints is defined here. This file is the single source:
it ships inside the installed package as
`echelle_spectra/resources/reading_room/vocabulary.md`, and `echelle web`
renders it from the installed package rather than from a repository checkout,
so a travel kit shows the same text a developer sees.

## Identities

- **Snapshot** — one immutable, digested calibration binder: an order pattern,
  a wavelength table, an integrating-sphere pair and an integral curve, copied
  once and never edited in place. A snapshot id looks like `20250926_cmos`.
- **Refinement** — a snapshot created by accepting a measured shift, named
  `<base>-r1`, `<base>-r2`. It is a new immutable snapshot, never an edit of
  the snapshot it came from.
- **Registry** — an ordered list of snapshot ids. Each referenced snapshot's
  own `[validity]` table carries the shot and date bounds, so an epoch
  boundary has exactly one authority.
- **Epoch** — one registry entry: the snapshot that a shot of a given number
  or date resolves to.
- **Catalog** — the last-known index of the cubes on one drive, written beside
  them. Merging catalogs makes an all-years index that keeps a disconnected
  drive's rows.
- **Run receipt** — the durable record of one batch: its source, output,
  expected file count, per-file outcomes, and the authorization it ran under.

## How a run was authorized

The `gate` word on a run says *why* that run was allowed to touch its
calibration epoch. The four words are not interchangeable.

| Gate | Meaning |
| --- | --- |
| `verdict` | A sampled drift audit produced evidence, and that evidence authorized this run. |
| `sample` | The legal first run of an epoch: an explicit unverified sample, marked as such in the receipt and in every cube it produced. |
| `ungated (no registry)` | The run used a hand-named configuration rather than the registry, so no epoch gate applied to it. |
| `unrecorded (pre-gate receipt)` | The receipt predates the gate. Nothing is claimed about how that run was authorized. |

A sampled drive is not a verified drive, and an ungated run is not a refused
one. The page renders all four differently on purpose.

## Drift verdicts

A drift audit measures Balmer and Fulcher line centroids on plasma-bright
frames, converts each wavelength residual into detector pixels using the
cube's own stored per-order dispersion, and classifies the interval:

- **aligned** — every measured line sits within the pixel tolerance.
- **shifted** — one rigid detector shift explains every line, and that shift
  is inside the repair limit. Repairable; it authorizes nothing until the
  refinement it names has been created and the registry repointed at it.
- **misaligned-beyond-repair** — the quorum was reached and neither rule
  holds. The affected shots need reprocessing from the raw SIF data.
- **insufficient-data** — the measured lines could not carry a verdict at
  all: too few resolved lines, too few distinct orders, or every line in one
  half of the audited range. **This never means aligned.**

Any other value in an evidence file is rendered as `unrecognized: <value>`
and is never dressed as one of the four.

## What this page says about a drive

Four different facts wear four different treatments here, because reading one
as another is how a wrong cube gets made:

- **missing drive** — the drive's catalog was not reachable when this page was
  built. Its rows are shown from the last merge. Nothing is claimed about the
  files themselves.
- **unmeasured** — no run receipt exists for that drive, so this page cannot
  say what was processed there or under what authorization.
- **empty** — a run was recorded and it produced no cubes. That is a measured
  zero, not an absence of measurement.
- **insufficient-data** — a real, judged drift verdict, described above.

Availability means the catalog file answered when the page was built. It never
means a cube is online, and a missing drive never means a deleted one.
