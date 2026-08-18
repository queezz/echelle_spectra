
<p align="center">
  <a href="https://doi.org/10.5281/zenodo.21371984">
    <img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21371984-blue.svg" alt="DOI">
  </a>
  <a href="https://www.python.org/downloads/release/python-395">
    <img src="https://img.shields.io/badge/python-3.9+-brightgreen.svg" alt="Python 3.9+">
  </a>
  <a href="https://github.com/queezz/echelle_spectra/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/queezz/echelle_spectra" alt="MIT License">
  </a>
  <a href="https://github.com/queezz/echelle_spectra/releases/latest">
    <img src="https://img.shields.io/github/v/release/queezz/echelle_spectra?include_prereleases&sort=semver" alt="Latest release">
  </a>
</p>

<h2 align="center">
  <br>
  <img src="src/echelle_spectra/resources/graphics/echelle.png" alt="Echelle Spectra Logo" width="60">
  <br>
  Echelle Spectra
  <br>
</h2>

<h4 align="center">Graphical tool for the extraction and analysis of calibrated spectra from 2D Echelle spectrometer images</h4>

<p align="center">
  <a href="https://queezz.github.io/echelle_spectra">Documentation</a> •
  <a href="CHANGELOG.md">Changelog</a> •
  <a href="docs/operator-cheat-sheet.md">Operator Cheat Sheet</a> •
  <a href="#where-to-start">Where to start</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#venv">Venv</a>
</p>

---

## Where to start

| I want to… | Use | Guide |
| --- | --- | --- |
| Make a new calibration at the instrument | `echelle-calib path/to/todays-calibration-folder` — the live bench | [Live calibration bench](docs/calibration-bench.md) |
| Turn a drive full of SIFs into cubes and LHD text | `echelle` — the campaign umbrella command | [Operator cheat sheet](docs/operator-cheat-sheet.md) |
| Look at one SIF | `echelle_spectra` — the single-SIF viewer | [Single-SIF viewer](docs/usage.md) |
| See the campaign as a page I can read and share | `echelle web --open`, reading a `campaign.toml` beside the campaign | [Operator cheat sheet](docs/operator-cheat-sheet.md#the-campaign-page) |

The primary workflow is **snapshot → registry → process**: the bench saves an
immutable calibration snapshot, the epoch registry says which shots it covers,
and `echelle process` converts a drive under it with the calibration's identity
and digests written into every cube. Which calibration file supplies which part
of a cube is [From calibration to cube](docs/calibration-to-cube.md).

---

## Quick Start

```bash
git clone https://github.com/queezz/echelle_spectra.git
cd echelle_spectra
python -m pip install -e .
echelle status           # what calibrations and runs exist
echelle_spectra          # open the single-SIF viewer
```

See the [operator cheat sheet](docs/operator-cheat-sheet.md) for every installed
surface, portable-kit command paths, the recommended NIFS trip loop, and fixes
for common command failures. The full [documentation](https://queezz.github.io/echelle_spectra)
contains the calibration and processing details.

For the pinned, no-admin travel payload, see the
[portable NIFS kit instructions](README-KIT.md). The released 1.6 kit preserves
its Python 3.12.13 runtime and pinned offline installation contract.

> `lab` is queezz's own private development helper for switching between
> environments on his own machines (`lab activate echelle-spectra`). It is not
> a dependency of anything here: the source package, the portable kit, and every
> command run perfectly well without it. Refresh an editable install with
> `uv pip install --no-deps -e .` after pulling a release that adds a new
> command.

---

## Live calibration bench

Use the separate pyqtgraph bench while acquiring calibration lamps:

```bash
echelle-calib path/to/calibration-folder
```

Drag SIFs onto the window (or press **Add SIF files…**; `--watch` adds optional
folder polling). Every file is triaged for exposure the moment it lands —
clustered saturation, counted cosmic-ray anomalies, remaining headroom, and the
counts histogram — before it has any role. Then carry the complete bench
procedure: assign sphere/background/lamp roles by hand with any lamp name,
follow the checklist built from your own data and its exposure advice, identify
lines from the shared packaged catalogs, fit the live rigid alignment, compare
new absolute factors with the previous sphere pair, save the alignment settings
as commented hand-editable files, then save and validate a complete snapshot
through the established snapshot API.

The folder argument is also where the output goes: by default the bench writes
`calibrations/` inside it (settings bundles in its `configs/` subfolder), and the Save tab shows
both paths in full with an **Open folder** button once a snapshot is saved. Use
`--output-root` and `--config-root` to send them elsewhere. See the
[live calibration bench guide](docs/calibration-bench.md).

---

## Calibration campaign commands

The `echelle` umbrella command is the front door for calibration campaigns:

```bash
echelle status
echelle snapshot --help
echelle process --help
echelle web --help
```

Calibration files can be assembled into an immutable, digested snapshot with
`echelle snapshot create` — usually the bench has already done this. Batch
processing writes durable per-source receipts, continues after ordinary
failures, reports rate and ETA, and safely resumes an interrupted run after
verifying completed source and output digests. Several source folders can be
processed together with one sequential worker per drive and isolated outputs,
receipts, progress, and failures. Existing `echelle-*` commands remain available
for compatibility. See the
[calibration snapshot guide](docs/calibration-snapshots.md) and
[campaign run guide](docs/campaign-runs.md). For multi-epoch processing, the
[calibration epoch registry](docs/calibration-epoch-registry.md) resolves and
verifies one immutable snapshot per source and writes complete recalibration
provenance into each cube.

### The campaign page

`echelle web` builds one self-contained static page from a cube catalog — the
Now stepper, the drives, the calibration evidence, and the packaged reading
room. It always prints the absolute path of the `index.html` it writes.

The simple form reads a hand-editable `campaign.toml` sitting beside the
campaign and opens the result in your browser in one step:

```bash
echelle web --open
```

```toml
# campaign.toml — the defaults echelle web reads when --catalog or --output
# is not given on the command line. Paths resolve against this file's own
# folder; any explicit --flag on the command line wins over these.
catalog = "all-years.json"
output = "campaign-page"
registry = "calibration_registry.toml"
calibrations = "calibrations"
drift = ["epoch-drift.json"]
```

Point `--home DIR-or-campaign.toml` at one kept elsewhere. The fully-flagged
form remains the explicit alternative — useful before a `campaign.toml`
exists, or to override every path from the command line:

```bash
echelle web \
  --catalog /data/all-years.json \
  --output /data/campaign-page \
  --registry /data/calibration_registry.toml \
  --calibrations /data/calibrations \
  --drift /data/epoch-drift.json
```

It writes `/data/campaign-page/index.html`. **Double-click that file to open
it**, or pass `--open` above to have the command do that for you — either way
there is no server and nothing is fetched. It is a snapshot of the moment it
was built, so rebuild the catalog and rerun `echelle web` after further
processing to refresh what it shows.

The campaign page's own composer asks for exactly two things — the data
folder and the calibration epoch — and derives the rest (cubes folder, volume
label, evidence file name, sample size), with every derived value editable
under an "Advanced" fold; see the
[campaign page composer](docs/harbor-candidate.md#the-campaign-page).

Learn the page with nothing: `echelle web --practice --open` — an invented
campaign built to a temp folder, labels say PRACTICE.

Prefer clicking to pasting? `echelle serve` serves the same page from this
machine (loopback only) and adds a real folder picker: Browse… beside the
data-folder field opens your drives, and a cold start with no campaign home
lets you pick the campaign folder and writes the `campaign.toml` for you.
The page still never executes commands — the terminal runs the work.

The current checkout also contains an **unreleased Packets 9–13 implementation
candidate** for drive catalogs, cube-derived LHD text, post-hoc cube
recalibration, sampled drift/refinement evidence, a read-only reading room, and
historical binders. See the [candidate contract](docs/harbor-candidate.md) and
[requirement reconciliation](docs/harbor-reconciliation.md). These surfaces do
not constitute releases 1.7.0–1.11.0 or a shipped Harbor train.

---

## CMOS notebooks — manual and tuning reference

The notebooks for the Andor CMOS (2560×2160) setup are a **manual and tuning
reference**, not the routine path. Routine calibration belongs at the bench and
routine conversion belongs to `echelle process`; both record provenance the
notebooks do not.

```
examples/workflows/black_cmos/
├── 01_load_image.ipynb             — sanity check: load image, verify dimensions
├── 02_automated_pattern_extraction.ipynb — run packaged pattern extraction and compare traces
├── 02_pattern_calibration.ipynb    — manual/tuning reference for pattern extraction
├── 03_wavelength_calibration.ipynb — manual line ID, polynomial fit  (rare: new setup only)
├── 04_extract_spectrum.ipynb       — one-file walkthrough: load → calibrate → extract → save
└── 05_calibration_alignment.ipynb  — align existing wavelength table with new Ne lamp data
```

`02` and `03` keep their value for a new optical setup or a pattern that needs
talking into place by hand. `04` is a readable walkthrough of one file, and is
superseded for real work by the bench plus `echelle process`. Historical
notebooks are in `examples/obsolete/`.

The packaged pattern extraction workflow can also be run without a notebook:

```bash
echelle-pattern sphere.sif sphere-bg.sif --prior-pattern pattern_CMOS_20240305.txt
```

The packaged wavelength alignment workflow can also be run headlessly:

```bash
echelle-align \
  local/20250926_calib/Ne-0.02s-x3-bright-lines.sif \
  local/20250926_calib/Ne-0.02s-x3-bright-lines_bg.sif \
  local/20250926_calib/sphere-0.1s-x3.sif \
  local/20250926_calib/sphere-0.1s-x3-bg.sif \
  --pattern pattern_CMOS_20250926.txt
```

---

## The viewer

![UI](images/gui.png)

`echelle_spectra` opens **one SIF at a time**: the detector image, the extracted
orders, and the calibrated spectrum. The Image tab can overlay separately
toggleable Balmer, Fulcher H2, ThAr, Ne, and Hg line markers on it. All overlays
start disabled; labels thin automatically in broad views and reveal local lines
as you zoom. The viewer and `echelle-validate-lines` use the same
provenance-carrying tables. See the
[known-line guide](docs/known-line-overlays.md).

For a calibration, use the bench; for a drive, use `echelle process`.

---

## Venv

### Create virtual environment

Linux / macOS:
```bash
python3 -m venv ~/.venvs/echelle-spectra
```

Windows PowerShell:
```powershell
python -m venv "$env:USERPROFILE\.venvs\echelle-spectra"
```

### Activate virtual environment

Linux / macOS:
```bash
source ~/.venvs/echelle-spectra/bin/activate
```

Windows PowerShell:
```powershell
& "$env:USERPROFILE\.venvs\echelle-spectra\Scripts\Activate.ps1"
```

### Install

```bash
pip install -e ".[dev]"
python -m ipykernel install --user --name echelle-spectra --display-name "echelle-spectra"
```
