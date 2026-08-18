# Procedure

The campaign loop, in the order it is actually walked. This file is the single
source: it ships as `echelle_spectra/resources/reading_room/procedure.md` and
is rendered from the installed package, so the kit on the trip laptop shows
this same text without a repository checkout.

Nothing on this page runs. Every command below is text to read, copy, and
paste into a terminal you are watching.

## 1. Ask what is already known

```
echelle status --calibrations calibrations --registry calibration_registry.toml --runs local/runs
```

This prints the valid snapshots, the registry's epochs and their bounds, the
most recent run receipt, and the gate that receipt records. It is the safest
first command of any session, and it never writes anything.

## 2. Make the first evidence legally

A registry-backed run of any size — a whole folder or a single file — refuses
to start without either an accepted drift verdict or an explicit unverified
sample. The sample is how the first evidence of a new epoch is produced:

```
echelle process shots -o cubes --registry calibration_registry.toml --sample 5
```

Every cube from that run is marked as an unverified sample, and the receipt
records `gate = "sample"`. A sample authorizes itself and nothing else.

## 3. Audit the epoch

```
echelle drift audit cubes --every 20 --catalog all-years.json -o drift-evidence.json
```

The audit fits line centroids over plasma-bright frames only, divides each
residual by the dispersion the cube's own polynomial has at that pixel, and
judges one rigid detector shift in pixels. It writes one immutable evidence
file — a second write to the same path is refused rather than overwritten.

Read the verdict card on this page before continuing. Look at the per-shot
table: if the shots fall into two groups, the evidence says so in an interval
warning, and the honest move is to split the interval and audit each half,
not to accept a shift that fits neither.

## 4. Repair only what is repairable

A `shifted` verdict composes its own repair sequence, and the sequence has
more than one step. Accepting the shift creates the next immutable `-rN`
snapshot; the registry must then be repointed at that refinement, because the
accepted shift authorizes the refinement and never the base snapshot it
corrected; cubes already exported are revised with `echelle recal-cube`.

`misaligned-beyond-repair` is not a repair path. Those shots need the raw SIF
data and a re-measured calibration. The verdict names the drives that hold
them when the merged catalog was supplied.

## 5. Process in bulk

With accepted evidence in hand, the bulk run states which evidence authorized
it, and the receipt keeps that path and its digest:

```
echelle process shots -o cubes --registry calibration_registry.toml --drift-verdict drift-evidence.json
```

The composer in the left rail asks two things — the folder holding the SIF
shots and the calibration to use — and derives everything else, writing this as
a plan file plus a short command, so the long argument list lives in a file you
can read and keep beside the run.

## 6. Catalog and re-read

```
echelle catalog build cubes --volume-label NIFS-A
echelle catalog merge cubes/echelle-catalog.json -o all-years.json
echelle web --catalog all-years.json --drift drift-evidence.json --output reading-room
```

Building the page again is how the loop closes: the merged index carries the
drives that were not connected, and this page renders that absence as absence.
