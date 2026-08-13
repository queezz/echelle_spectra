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

`snapshot.toml` remains ordinary hand-editable TOML. If it is edited, run the
validator again. A changed artifact needs a new snapshot ID; an already used
snapshot should never be silently rewritten.

## Historical calibrations

Existing 2019, 2024, and 2025 calibration files are historical evidence and are
not renamed. They will receive thin manifests that point to their existing role
files. New campaign snapshots use the role-named folder convention from birth.
