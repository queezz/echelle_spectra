# Durable campaign runs

`echelle process` turns a folder of SIF files into SpectroCubes without making
you keep the run state in your head. Every non-dry batch creates a durable run
receipt under `local/runs/` by default.

## Batch discovery

A plain filename pattern — the default `*.SIF`, or anything else with no `/`
and no `**` — walks the *whole* source tree under the folder you named, not
just its top level: real campaign drives are date-named day folders, and a
shot several levels down is still found. Two kinds of folder are pruned from
that walk so lamp frames are never eaten as science shots: any folder named
`calibrations`, and any folder holding a `snapshot.toml`. A pattern containing
`/` or `**` is used exactly as typed, with no pruning applied. A case-
insensitive `*.sif` fallback still runs when the primary pattern finds
nothing, so lowercase acquisition filenames are not a separate step.

A registry-backed drive takes three commands: sample the epoch, audit the
sample, then process the drive under the verdict the audit wrote.

```powershell
# 1. Sample the epoch. --sample auto derives the count from what the folder
#    holds (max(5, min(30, ceil(files/25))), or every file when fewer than
#    5); no verdict is needed yet, and every produced cube is marked
#    drift_sample.
echelle process D:\NIFS\shots `
  -o D:\NIFS\spectrocubes `
  --registry D:\NIFS\calibration_registry.toml `
  --calibrations D:\NIFS\calibrations `
  --volume-label NIFS-A `
  --sample auto

# 2. Audit the sample into one immutable verdict file. The audit takes the
#    cube folder itself, so the same line works in any shell. -o is optional:
#    left out, the evidence is written beside the audited cubes as the next
#    free drift-evidence-NNN.json, numbered from 001 (evidence is immutable,
#    so a rerun takes the next name rather than overwriting).
echelle drift audit D:\NIFS\spectrocubes `
  --calibrations D:\NIFS\calibrations

# 3. Process the whole drive under that verdict, naming the file the audit
#    printed.
echelle process D:\NIFS\shots `
  -o D:\NIFS\spectrocubes `
  --registry D:\NIFS\calibration_registry.toml `
  --calibrations D:\NIFS\calibrations `
  --volume-label NIFS-A `
  --drift-verdict D:\NIFS\spectrocubes\drift-evidence-001.json
```

Step 3 resumes the receipt written by step 1: the sampled files are verified and
skipped rather than exported twice.

The terminal reports the current count, measured file rate, and estimated time
remaining. An ordinary file failure is recorded and processing continues with
the next source. The command exits nonzero after the run if any source failed.

## The gate

Nothing reaches a registry-selected calibration without saying how it was
authorized. Every run records one of three gates in its receipt:

| Gate | Meaning |
| --- | --- |
| `verdict` | `--drift-verdict` named sampled evidence that covers every snapshot this run selected. The receipt also stores the evidence path, its SHA-256, and the verdict word. |
| `sample` | `--sample N` (or `--sample auto`, which derives N from the folder's file count) processed at most the first N resolved files (sorted by path) with no verdict at all. The receipt stores `sample = true` and the file count; every cube carries `drift_sample = 1`. |
| `ungated (no registry)` | The run used an explicit `--config`/`--snapshot-id`/`--camera` calibration. That is legal, and the receipt says forever that no sampled audit stood behind it. |

`--sample` and `--drift-verdict` cannot be combined: one takes unverified
evidence, the other spends it. `--sample` requires `--registry`.

The gate covers every registry-backed form — a folder of one file, a folder of
ten thousand, and the single-file path. A refusal names the commands that
produce the missing evidence instead of only stating that it is missing.

Aligned evidence authorizes the snapshots it audited. An accepted shifted
verdict authorizes **only** the `-rN` refinement it created: while the registry
still selects the condemned base snapshot, the gate refuses and tells you to
repoint the registry's validity at the refinement.

Every audited shot also carries an `isotope` tag: each Balmer window is judged
against both the hydrogen and the deuterium reference and assigned to whichever
is nearer in detector pixels, so a deuterium shot reads as aligned rather than
shifted by the 0.178 nm isotope offset, drops the H2 Fulcher anchors it has no
D2 table for, and is flagged where the bundled LHD deuterium calendar expects
hydrogen instead.

Receipts written before the gate existed remain readable. They report
`unrecorded (pre-gate receipt)` rather than borrowing an authorization they
never had.

## Several drives at once

Pass one source folder per drive. Echelle runs one sequential reader for each
source and processes those sources concurrently:

```powershell
echelle process D:\NIFS\shots E:\NIFS\shots `
  -o F:\NIFS\spectrocubes `
  --runs-dir F:\NIFS\runs `
  --volume-label NIFS-A `
  --volume-label NIFS-B `
  --registry F:\NIFS\calibration_registry.toml `
  --calibrations F:\NIFS\calibrations `
  --drift-verdict F:\NIFS\epoch-drift.json
```

Each target is gated on its own. `--sample 5` samples each source folder
separately, so one campaign command can produce the audit input for every drive
at once; the verdict that covers all of their snapshots then authorizes the
full multi-drive run.

`-o` is a destination root in this form. Each source receives its own named
child directory, its own sequential calibration/export worker, and its own
receipt tree under `--runs-dir`. Repeat `--volume-label` in the same order as
the source folders, or omit all labels to use the drive/root identity. A file
failure or calibration failure on one target does not stop another target.
`--run-dir` remains the exact-resume control for a single source; rerun the same
multi-source command to resume all matching unfinished target receipts.

Progress lines are prefixed with their drive label. The command exits `0` only
when every target completes, `1` after any ordinary target failure, and `130`
after Ctrl-C. On interrupt, active atomic exports finish or clean up safely and
all workers stop before taking another source.

## What is recorded

Each run directory contains:

- `run.toml` — atomic summary containing run state, roots, volume label,
  calibration authority (`per-source-registry` for an epoch run), expected file
  count, the gate that authorized the run with its evidence path, digest, and
  verdict, and current status counts;
- `records.jsonl` — an append-only, one-record-per-attempt ledger with source
  path, source size and SHA-256, output path, result, reason, timing, volume,
  and the snapshot selected for that source. Successful exports also record
  output size and SHA-256.

These files belong beside local campaign data and must not be committed. A
snapshot ID can be omitted for legacy work; receipts then say `unassigned`
rather than inventing provenance.

## Interrupt and resume

Pressing Ctrl-C marks the run `interrupted` and exits with code 130. Outputs are
first written to a temporary sibling and published atomically, so interruption
does not leave a partial file at the final `.nc` path.

Rerun the same command. Echelle automatically finds the newest unfinished run
with the same source and destination. Before skipping prior work it verifies:

1. the source path, size, and SHA-256 still match;
2. the output path, size, and SHA-256 still match;
3. the prior terminal record says the export completed;
4. for registry runs, the currently selected snapshot ID still matches the
   recorded snapshot.

Only then is the source recorded as skipped because its completed output was
verified. A changed epoch selection, changed source, or changed/missing output
is not accepted as completed. To name the resume explicitly, use the path
printed at interruption:

```powershell
echelle process D:\NIFS\shots -o D:\NIFS\spectrocubes `
  --run-dir local\runs\2026-08-13_12-00-00-shots
```

Use `--new-run` when you intentionally want a separate receipt. Use
`--overwrite` only when existing outputs should be replaced.

Repeat `--registry` and its `--drift-verdict` when resuming a registry-backed
run. The gate is evaluated on every invocation, including a resume, and the
receipt records the authorization of the run in front of it.

## Drive identity

On the first processing of a source, the run writes a small
`echelle-drive-id.toml` at the source root you named:

```toml
schema = "echelle-drive-id/v1"

[drive]
id = "0b6c1e2f9a5b4c7d8e9f0a1b2c3d4e5f"
label = "NIFS-A"
created_at = "2026-08-14T09:12:04.118Z"
```

Later runs read it instead of writing a new one, and catalogs key on the `id`.
A USB disk that returns as `D:` on Windows and under `/Volumes/NIFS` on macOS
therefore stays one drive with one history. The `label` is display only: edit it
freely, and never copy an `id` onto a second drive. `--volume-label` sets the
label a run reports; it never changes an announced id.

A source root that cannot be written — a read-only archive mount — is not an
error. The run continues, keys that drive on its volume label alone, prints a
warning, and records it as `drive_warning` in the receipt and in the catalog's
`run` block, so a reader can see why the drive has no stable id.

## Catalog beside the cubes

Every non-dry batch target writes `echelle-catalog.json` beside its cubes. Each
cube row records its digest, size, snapshot ID, `t_start`, and the `gate` that
authorized the run that produced it, so a sampled drive never reads as a
verified one. Write or refresh one by hand with the same command the run uses:

```powershell
echelle catalog build D:\NIFS\spectrocubes `
  --volume-label NIFS-A `
  --receipt-dir local\runs\2026-08-13_12-00-00-shots
```

`--receipt-dir` is optional and adds the run's state, counts, calibration
authority, and gate to the catalog. `--drive-id` is optional too: by default the
catalog uses the id announced by the nearest `echelle-drive-id.toml` at or above
the cube folder.

Merge per-drive catalogs into a durable all-years index with
`echelle catalog merge CATALOG [CATALOG ...] -o all-years.json`. Two rules make
that merge safe to repeat:

- **Recency, never argument order.** For each drive the row with the newest
  `generated_at` wins, so merging a stale all-years index after a fresh
  per-drive catalog can no longer revert that drive to older rows.
- **Portable locations.** A merged row stores the catalog path *relative* to its
  drive root plus the `drive_root` this machine last saw it at — the catalog
  file's own folder. A drive that comes back at another letter or mount point is
  relocated by merging its catalog again; nothing is marked permanently
  unavailable by an absolute path written on a different operating system.

Merge automatically after every completed target by naming the index on the run:

```powershell
echelle process D:\NIFS\shots -o D:\NIFS\spectrocubes `
  --volume-label NIFS-A `
  --central-index C:\NIFS\all-years.json
```

Only a completed target is folded in; an interrupted or partial run still writes
its per-drive catalog but is not published into the index as though the drive
were finished. Without `--central-index`, the merge stays manual.

`echelle recal-cube` refreshes the per-drive catalog beside its output when one
exists: the revised cube's row gets the new digest and snapshot ID, and the
catalog's `generated_at` is bumped. That bump is what makes the drive outrank
whatever the all-years index still remembers, so the next auto or manual merge
propagates the correction without rebuilding the index.

## Status

```powershell
echelle status --runs local\runs
```

Status reports the newest run's state, accounted sources, result counts,
calibration authority, and the gate that authorized it. Individual records
retain their selected snapshot IDs.
When several source/output targets have receipts, it also reports
their newest states individually and reconciles their combined totals. It reads
receipt files recursively; it does not infer progress from filenames or count
superseded retries twice.

Dry runs remain side-effect free: they create neither outputs nor receipts.
