# Provenance

What each artifact records about where it came from, and what this page is
allowed to say about it. This file is the single source: it ships as
`echelle_spectra/resources/reading_room/provenance.md` and is rendered from
the installed package.

## The chain

1. **Raw SIF data** is read and never written. Every identity below leads back
   to it, and nothing in this toolchain edits it.
2. **A snapshot** copies the calibration inputs once and digests each of them.
   Its manifest carries the role filenames, sizes and SHA-256 digests, so a
   snapshot that was quietly edited afterwards no longer validates.
3. **A cube** stamps the snapshot id, the snapshot manifest digest, its own
   source file, the wavelength medium, the intensity units, and — for cubes
   written by the current exporter — the raw `detector_pixel` and
   `echelle_order` of every retained wavelength column plus the per-order
   polynomial that produced it.
4. **A run receipt** records the source root, output root, expected file
   count, per-file outcome, the snapshot id, the gate, and, when a verdict
   authorized the run, the evidence path and its digest.
5. **A catalog** records each cube's relative path, digest, size, shot, year,
   snapshot id, coverage, and the gate under which it was published.
6. **Drift evidence** records the sampled cubes and their digests, the sample
   rule, the thresholds it judged against, the per-line and per-shot
   measurements, the verdict, and the repair sequence it composed.
7. **A recalibration manifest** sits beside a revised cube and records the old
   and new snapshot manifests plus the input and output digests.

## What "available" means here

The catalog rows on this page come from the merged all-years index. A drive is
shown as available when its own catalog file answered at the moment this page
was built, and as a missing drive when it did not. Neither statement is a
claim about the cubes themselves: a reachable catalog does not prove the files
beside it are readable, and an unreachable one does not mean anything was
deleted. Local availability never rewrites recorded identity.

## What this page cannot tell you

- Whether a cube's bytes still match its recorded digest. Nothing here
  re-hashes anything; the page renders saved records.
- Whether the registry currently selects the snapshot a past run used. The
  receipt says what that run resolved; the registry may have been repointed
  since.
- Whether a drive's absence is temporary. Absence is rendered as absence.

## Reading a repair command

The repair rows are composed from the evidence, not typed by hand, and they
are shown meaning first: the plain-words purpose leads, the literal command
sits behind a toggle, and the copy button carries the whole command whether
that toggle is open or closed. A step with no command is a step you perform
yourself — repointing a registry entry is a decision, not a paste.
