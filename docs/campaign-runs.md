# Durable campaign runs

`echelle process` turns a folder of SIF files into SpectroCubes without making
you keep the run state in your head. Every non-dry batch creates a durable run
receipt under `local/runs/` by default.

```powershell
echelle process D:\nifs\shots `
  -o D:\nifs\spectrocubes `
  --snapshot-id 20260813_cmos `
  --volume-label NIFS-A
```

The terminal reports the current count, measured file rate, and estimated time
remaining. An ordinary file failure is recorded and processing continues with
the next source. The command exits nonzero after the run if any source failed.

## Several drives at once

Pass one source folder per drive. Echelle runs one sequential reader for each
source and processes those sources concurrently:

```powershell
echelle process D:\nifs-a\shots E:\nifs-b\shots `
  -o F:\nifs-cubes `
  --runs-dir F:\nifs-runs `
  --volume-label NIFS-A `
  --volume-label NIFS-B `
  --snapshot-id 20260813_cmos
```

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
  snapshot ID, expected file count, and current status counts;
- `records.jsonl` — an append-only, one-record-per-attempt ledger with source
  path, source size and SHA-256, output path, result, reason, timing, volume,
  and snapshot. Successful exports also record output size and SHA-256.

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
3. the prior terminal record says the export completed.

Only then is the source recorded as skipped because its completed output was
verified. Changed or missing output is not accepted as completed. To name the
resume explicitly, use the path printed at interruption:

```powershell
echelle process D:\nifs\shots -o D:\nifs\spectrocubes `
  --run-dir local\runs\2026-08-13_12-00-00-shots
```

Use `--new-run` when you intentionally want a separate receipt. Use
`--overwrite` only when existing outputs should be replaced.

## Status

```powershell
echelle status --runs local\runs
```

Status reports the newest run's state, accounted sources, result counts, and
snapshot ID. When several source/output targets have receipts, it also reports
their newest states individually and reconciles their combined totals. It reads
receipt files recursively; it does not infer progress from filenames or count
superseded retries twice.

Dry runs remain side-effect free: they create neither outputs nor receipts.
