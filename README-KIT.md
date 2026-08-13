# Portable NIFS kit

This folder is the source contract for the travel kit. A generated kit contains
the files described here plus a target-specific uv binary, cached CPython,
project wheel, dependency wheelhouse, and `checksums.sha256`. It does not need
Dropbox, Git, an existing Python installation, administrator rights, Fleet, or
lab-cli.

> **Development status (2026-08-13):** the pinned manifest, assembly path,
> installers, offline wheelhouse contract, and reproducible-distribution tooling
> are implemented for the existing 1.4.0 package. Packet 6 and release 1.5.0 are
> not shipped until clean Windows, clean macOS, and network-disabled offline
> rehearsals all pass. A kit built from the current source is a rehearsal kit.

## Windows x86-64: three commands

Open ordinary Windows PowerShell in the generated kit folder:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\install.ps1 offline
.\echelle.ps1 status
```

Use `.\install.ps1 online` for the network-backed route. Both modes use the
same exact, hash-locked dependency versions and the same local project wheel.
Offline mode adds `--offline --no-index` and can succeed only from the packaged
wheelhouse.

## macOS: three commands

Open Terminal in the generated kit folder:

```bash
chmod +x install.sh echelle bin/uv
./install.sh offline
./echelle status
```

Use `./install.sh online` for the network-backed route. Apple Silicon and Intel
Macs have separate payloads; the installer refuses the wrong architecture.

## What the installer proves

Before changing the kit folder, the installer verifies every packaged file
against `checksums.sha256`, refuses the wrong OS/architecture, extracts the
exact cached CPython 3.12.13 runtime, creates a local `.venv` with uv 0.11.31,
installs the hash-locked dependencies, installs the local Echelle wheel without
resolving new dependencies, and runs `echelle status`.

The following states are intentionally distinct:

- **online bootstrap** — dependency artifacts come from the index but must match
  the committed hashes;
- **offline bootstrap** — network access is disabled in uv and every dependency
  must be present in `wheelhouse/`;
- **installed verification** — the installed `echelle status` command succeeds;
- **kit validation** — platform, completeness, file sizes, and SHA-256 digests
  match the generated payload inventory.

If `.runtime/` or `.venv/` exists but is incomplete, the installer refuses to
guess. Remove only the named generated folder, leave `runtime/`, `wheelhouse/`,
and `checksums.sha256` unchanged, and rerun. A checksum failure means the kit
must be rebuilt or recopied; do not install from it.

uv's disposable cache is kept under `.cache/uv` inside the generated kit, so
installation neither depends on nor mutates a user-global uv cache. It may be
removed at any time without changing the checksummed payload.

## Supported payloads

| Kit key | Operating system | Architecture | Python |
| --- | --- | --- | --- |
| `windows-x86_64` | Windows 10/11 | x86-64 (AMD64) | CPython 3.12.13 |
| `macos-aarch64` | macOS 14 or newer | Apple Silicon | CPython 3.12.13 |
| `macos-x86_64` | macOS 13 or newer | Intel x86-64 | CPython 3.12.13 |

Windows ARM64, 32-bit Windows, Linux, and macOS older than the rows above are
not supported by this kit contract.

## Maintainer assembly

Each payload carries its own target-specific, hash-locked
`requirements-runtime.txt`; the three source locks live under
`kit/requirements/`. This is necessary because the supported binary wheels
differ by target. The pure-Python `tzdata` pin is carried in every payload to
keep pip hash checking deterministic even when cross-assembling on Windows.
Generated binaries, Python
archives, wheels, caches, kits, and build outputs
must stay in OS-local scratch outside the repository. Build the project wheel
from an exported clean source copy first, then assemble one target:

```powershell
$Scratch = Join-Path $env:TEMP "echelle-nifs-kit"
$Source = Join-Path $Scratch "source-a"
$Dist = Join-Path $Scratch "dist-a"
$CompanionDist = Join-Path $Scratch "spectrocube-dist"
$Kit = Join-Path $Scratch "kit windows x86_64"
$Cache = Join-Path $Scratch "verified-cache"

git archive --format=zip --output=(Join-Path $Scratch "source-a.zip") HEAD
Expand-Archive -LiteralPath (Join-Path $Scratch "source-a.zip") -DestinationPath $Source
Push-Location $Source
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m scripts.reproducible_build build --source $Source --destination $Dist --manifest (Join-Path $Source "kit\manifest.toml")
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m scripts.companion_wheels --manifest (Join-Path $Source "kit\manifest.toml") --spectrocube-source (Join-Path $Scratch "spectrocube-source") --sif-parser-sdist (Join-Path $Cache "downloads\sif_parser-0.3.6.tar.gz") --destination $CompanionDist
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m scripts.nifs_kit --manifest (Join-Path $Source "kit\manifest.toml") --platform windows-x86_64 --destination $Kit --cache $Cache --project-wheel (Get-ChildItem -LiteralPath $Dist -Filter "echelle_spectra-*.whl" -File).FullName --companion-wheel (Get-ChildItem -LiteralPath $CompanionDist -Filter "spectrocube-*.whl" -File).FullName --companion-wheel (Get-ChildItem -LiteralPath $CompanionDist -Filter "sif_parser-*.whl" -File).FullName
Pop-Location
```

`$CompanionDist` contains two normalized universal wheels. SpectroCube is built
from a clean archive of the manifest's exact commit,
`c46b0c4621ec8d56d7d7209d122b571fb109405e`. `sif_parser` is built from its
manifest-pinned PyPI source archive after its 18,722-byte input and SHA-256 are
verified. SpectroCube 0.1.0 is not published in the configured package
registry, while `sif_parser` 0.3.6 publishes no wheel; carrying both verified
wheels removes the hidden checkout and compiler dependencies from target
machines. Pass both paths as repeated `--companion-wheel` arguments.

Before the example, export the SpectroCube commit as `spectrocube-source`
without a `.git` directory and download the manifest's exact `sif_parser`
source URL to the shown cache path. The companion builder refuses a live Git
checkout, a different source-archive size or digest, a wrong wheel identity,
or a non-empty output directory.

Supply `--wheelhouse-source <verified-folder> --offline-assets` to repeat
assembly without any downloads. An existing complete destination is accepted
unchanged; an existing incomplete or corrupt destination is refused.

## Distribution reproducibility contract

The release contract is byte identity of the **normalized final wheel and
sdist**, not the build backend's temporary archives. The build command pins
`build`, setuptools, wheel, Python, `PYTHONHASHSEED`, timezone, and
`SOURCE_DATE_EPOCH`, then normalizes archive ordering, timestamps, ownership,
permissions, gzip headers, and ZIP metadata. Two independent exported source
copies must produce matching filenames and SHA-256 digests:

```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\python.exe" -m scripts.reproducible_build compare $DistA $DistB
```

Any differing digest fails the release gate. `twine check`, installed-artifact
smokes, both native platform rehearsals, and a demonstrably network-disabled
offline install remain separate required gates.
