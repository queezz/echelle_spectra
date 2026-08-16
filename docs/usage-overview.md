# Where to start

Three different jobs, three different tools. Pick the row that matches what you
are about to do.

| What you are doing | Use | Start here |
| --- | --- | --- |
| **Calibration day at the instrument** — lamps and sphere are on the table, you are producing a new calibration | `echelle-calib` — the live calibration bench | [Live calibration bench](calibration-bench.md) |
| **Batch work on a campaign drive** — turning a drive full of SIFs into cubes and LHD text, and keeping track of what was processed under which calibration | the `echelle` umbrella command | [Operator cheat sheet](operator-cheat-sheet.md) |
| **Looking at one SIF** — checking a single shot, overlaying known lines, reading off a band | `echelle_spectra` — the viewer | [Basic usage](usage.md) |

---

## Calibration day: the bench

```bash
echelle-calib /data/incoming
```

Drag SIFs onto the window. The bench triages each exposure the moment it lands
(saturation, cosmic rays, headroom, histogram), takes the roles you assign by
hand, identifies lines, fits the rigid alignment, compares the new sphere
response against the previous one, and ends with one Save that writes an
immutable, validated snapshot plus a hand-editable configuration bundle.

By default both land inside the folder you launched it at —
`/data/incoming/calibrations/` and `/data/incoming/calibration-configs/`. The
Save tab shows both paths in full, and offers **Open folder** once something has
been saved.

→ [Live calibration bench](calibration-bench.md) ·
[Calibration alignment](calibration-alignment.md) ·
[Calibration snapshots](calibration-snapshots.md)

---

## Campaign drives: the `echelle` umbrella

This is the batch path — many SIFs in, cubes and text out, with a receipt for
every run and full calibration provenance inside every cube.

```mermaid
flowchart LR
    A["echelle status"] --> B["snapshot<br/>create / validate"]
    B --> C["registry<br/>calibration_registry.toml"]
    C --> D["process --sample N<br/>first unverified sample"]
    D --> E["drift audit<br/>one verdict file"]
    E --> F["process --drift-verdict<br/>the whole drive"]
    F --> G["txt · catalog"]
    G --> H["web<br/>campaign page"]
```

| Step | Command | What it is for |
| --- | --- | --- |
| Check | `echelle status` | What snapshots, registry, and run receipts exist. Safest first command. |
| Calibrate | `echelle snapshot create` / `validate` / `show` | Build or recheck an immutable snapshot. Usually the bench already did this. |
| Place it | edit `calibration_registry.toml` | Say which shots and dates each snapshot covers. |
| Sample | `echelle process ... --sample N` | Process the first N files, marked as an unverified sample. |
| Audit | `echelle drift audit` | Turn that sample into one immutable verdict file. |
| Process | `echelle process ... --drift-verdict FILE` | Run the whole drive under that verdict. |
| Text | `echelle txt CUBE OUT` | Write LHD text from a saved cube — no raw SIF needed. |
| Catalog | `echelle catalog build` / `merge` | Per-drive catalogs, and one all-years index. |
| Publish | `echelle web` | Build the static campaign page you open in a browser. |

The sample → audit → process gate is not optional for a registry-backed run: no
drive is processed under a calibration epoch nothing has checked. The full
copy-paste loop, in both shells, is in the
[operator cheat sheet](operator-cheat-sheet.md).

→ [From calibration to cube](calibration-to-cube.md) ·
[Calibration epoch registry](calibration-epoch-registry.md) ·
[Durable campaign runs](campaign-runs.md) ·
[Batch SpectroCube workflow](spectrocube-batch-workflow.md)

---

## One SIF: the viewer

```bash
echelle_spectra
```

The `echelle_spectra` window is a **single-SIF viewer**. Open one file, see the
detector image and the extracted, wavelength-calibrated orders, toggle Balmer,
Fulcher H₂, ThAr, Ne, and Hg overlays on top of it, and look at the defined
[band windows](band-data.md).

It is the right tool for inspecting a shot. It is not the tool for producing a
calibration (that is the bench) or for converting a drive (that is
`echelle process`).

→ [Basic usage](usage.md) · [Known-line overlays](known-line-overlays.md)

---

## Reading the results back

- **`echelle txt`** writes the frozen-header LHD text from a saved cube.
- **`echelle catalog build` / `merge`** index what is on each drive and across
  all years.
- **`echelle web`** turns a merged catalog into one self-contained `index.html`
  you double-click to open — no server, nothing fetched. Rebuild it after a run
  to refresh what it shows.

```bash
echelle web \
  --catalog /data/all-years.json \
  --output /data/campaign-page \
  --registry /data/calibration_registry.toml \
  --calibrations /data/calibrations \
  --drift /data/epoch-drift.json
```

Then open `/data/campaign-page/index.html`.

---

## Specialist commands

For recalibration work outside the bench, and for QC:
`echelle-pattern`, `echelle-align`, `echelle-validate-lines`,
`echelle-wavelength-qc`, `echelle-nist-overlay`. Every one of them is listed
with its purpose in the [operator cheat sheet](operator-cheat-sheet.md).
