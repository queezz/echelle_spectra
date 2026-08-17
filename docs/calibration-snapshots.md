# Calibration snapshots

A calibration snapshot is one immutable folder naming every file needed to
extract and calibrate spectra. Its folder name is the identity used by cubes,
registries, catalogs, and run receipts:

```text
20260901_calib/                        the calibration folder
├── sphere-0.1s-x3.sif                 raw frames, referenced from below
├── sphere-0.1s-x3-bg.sif
├── thar.sif
└── calibrations/20260901_cmos/        the snapshot
    ├── snapshot.toml
    ├── pattern.txt
    ├── wavelength.txt
    └── integral.txt
```

The ID format is `YYYYMMDD_<detector>[-rev]`. A refinement derived later from
science lines can use a revision suffix such as `20260901_cmos-r1`.

## The folder holds the light; the snapshot holds what was computed

A snapshot records two kinds of artifact, and `snapshot.toml` says plainly which
each one is:

| `kind` | Where the bytes are | Which files |
| --- | --- | --- |
| `copied` (the default, written by saying nothing) | inside the snapshot folder, at a relative path that cannot leave it | `pattern.txt`, `wavelength.txt`, `integral.txt` |
| `referenced` | wherever they were measured, named by `path` | the sphere pair and every lamp signal and background SIF |

```toml
[[artifacts]]
role = "wavelength"
path = "wavelength.txt"
sha256 = "1f0c…"
size_bytes = 20416

[[artifacts]]
role = "sphere"
kind = "referenced"
path = "../../sphere-0.1s-x3.sif"
sha256 = "9ab4…"
size_bytes = 398458880
source_name = "sphere-0.1s-x3.sif"
```

The digest is the identity either way. The path only says where the bytes live,
and validation reads them there and hashes them again. A reference is stored
relative to the snapshot folder when the source is inside the same calibration
folder, so the whole folder can be moved or copied to another machine with no
edit; a source from anywhere else is stored as an absolute path, because a
relative path across unrelated trees would only pretend to be portable.

Referencing is why raw frames are not duplicated. A calibration folder already
holds its lamp and sphere SIFs — hundreds of megabytes of them — and a snapshot
two levels below has nothing to gain from a second copy sitting beside the
first. Older snapshots that do hold copies are still perfectly valid and are
read exactly as before: an entry with no `kind` is a copied one.

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

`echelle snapshot create` copies every role, which is what assembling a snapshot
out of scattered files should do. The live bench references the raw frames
instead, because it is standing in the folder that already holds them; see
[the calibration bench](calibration-bench.md).

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
- relative paths that cannot leave the snapshot directory, for copied artifacts;
- presence, byte size, and SHA-256 digest of every artifact, copied or
  referenced — a referenced source is re-read and re-hashed where it is named;
- base-snapshot identity and self-reference.

A referenced source that has gone missing or changed is a plain validation
failure that names the full path it looked in:

```text
Snapshot is invalid:
  - referenced sphere source not found: T:\2025-LHD-BH\Echelle\20260901_calib\sphere-0.1s-x3.sif
```

Relative references are read from the snapshot folder, never from the working
directory, so `echelle snapshot validate` gives the same answer from anywhere.
`echelle snapshot show` prints each referenced frame with the full path it
resolves to, so you can see what a snapshot actually uses without opening the
binder.

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
