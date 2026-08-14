# Calibration snapshots

A calibration snapshot is one immutable folder containing the files needed to
extract and calibrate spectra. Its folder name is the identity used by cubes,
registries, catalogs, and run receipts:

```text
calibrations/20260901_cmos/
├── snapshot.toml
├── pattern.txt
├── wavelength.txt
├── sphere.sif
├── sphere_bg.sif
├── integral.txt
└── lamps/
    ├── h2.sif
    └── thar.sif
```

The ID format is `YYYYMMDD_<detector>[-rev]`. A refinement derived later from
science lines can use a revision suffix such as `20260901_cmos-r1`.

## Create a snapshot

The command copies source files into role-named locations, digests every file,
writes `snapshot.toml`, validates the finished folder, and refuses to replace an
existing snapshot:

```powershell
echelle snapshot create 20260901_cmos `
  --detector cmos `
  --pattern pattern-source.txt `
  --wavelength wavelength-source.txt `
  --sphere sphere-source.sif `
  --sphere-background sphere-background-source.sif `
  --integral integrating-sphere-source.txt `
  --lamp-used H2 `
  --lamp-used Th-Ar `
  --lamp-file H2=h2-source.sif `
  --lamp-file Th-Ar=thar-source.sif `
  --notes "NIFS/LHD 2026 campaign"
```

Use `--base-snapshot` when the new wavelength solution was aligned from an
earlier snapshot. Optional shot/date validity, alignment-transform, and QC
summary arguments are shown by `echelle snapshot create --help`.

The source files are opened read-only. Construction happens in a temporary
sibling directory and the final folder appears only after every copy, digest,
manifest write, and integrity check succeeds.

## Validate or inspect

```powershell
echelle snapshot validate calibrations/20260901_cmos
echelle snapshot show calibrations/20260901_cmos
echelle status
```

Validation checks:

- the `echelle-snapshot/v1` schema and folder/ID naming contract;
- detector consistency with the snapshot ID;
- required roles (`pattern`, `wavelength`, `sphere`, `sphere_background`, and
  `integral`);
- relative paths that cannot leave the snapshot directory;
- presence, byte size, and SHA-256 digest of every artifact;
- base-snapshot identity and self-reference.

`snapshot.toml` begins with comments explaining that it remains ordinary
hand-editable TOML. If it is edited, run the validator again. A changed artifact
needs a new snapshot ID; an already used snapshot should never be silently
rewritten.

The live calibration bench can assemble the same request after its procedure is
complete. It calls this creation API and then this validator directly, so the GUI
inherits the same atomic construction, digest checks, and replacement refusal.
Its separately generated campaign, alignment, and export TOMLs remain editable
configuration evidence; they do not replace `snapshot.toml` as the binder.

## Historical calibrations

Existing 2019, 2024, and 2025 calibration files are historical evidence and are
not renamed. Thin bundled binders point at them under their original names; list
them with `echelle historical`. New campaign snapshots use the role-named folder
convention from birth.

A binder is evidence, not a snapshot: the registry reads only
`echelle-snapshot/v1` with a `[validity]` table. Convert one with

```powershell
echelle snapshot import-historical 20250926_cmos `
  --calibrations calibrations `
  --artifact-root local\20250926_calib `
  --valid-from 2025-09-26
```

which copies the named artifacts into an immutable, digested snapshot carrying
the epoch you declare and an `[imported]` table naming the binder it came from.
See [the epoch registry guide](calibration-epoch-registry.md) for the artifact
roots and the refusal that names a missing one.
