# Known-Line Tables and GUI Overlays

Echelle Spectra carries one shared API for the known lines used by the existing
GUI and by future calibration-bench tools. The public families are `balmer`,
`fulcher`, `thar`, `ne`, and `hg`.

```python
from echelle_spectra import load_line_table

for line in load_line_table("fulcher"):
    print(line.label, line.wavelength_nm, line.source_reference)
```

Every `SpectralLine` records its family, label, air wavelength, species, source
name and reference, bundled resource, optional relative intensity, and notes.
Callers can use `filter_line_table()` to select a wavelength or intensity range
without losing those provenance fields.

## Sources and conventions

| family | bundled knowledge | provenance and convention |
| --- | --- | --- |
| Balmer | H-alpha through H-delta | Established Echelle air-wavelength validation values, shared with `emissiondata` and `echelle-validate-lines`. |
| Fulcher H2 | 42 Q-branch lines in 0-0 through 3-3 | Numerically copied from the Fulcher Extractor resource archived at [Zenodo](https://doi.org/10.5281/zenodo.21372039). That resource traces the values to its early manual-extraction table but does not identify a primary publication. Echelle interprets them as air wavelengths to preserve its accepted validation convention. |
| ThAr, Ne, Hg | Atomic and ionic lamp lines | Package-cached exports from the [NIST Atomic Spectra Database](https://physics.nist.gov/asd), covering 380–810 nm. Observed air wavelengths are preferred, with Ritz air values used when no observation is present. Relative intensity is normalized within each spectrum and then scaled by ionization stage, so a lamp's species can be ranked against one another — see below. |

### Why relative intensity is not NIST's number verbatim

NIST prints each spectrum's `Rel. Int.` on the scale its own source reference
used, and those scales do not line up between ionization stages. Ne I and Ne II
come from different line references and top out at 100000 and 400; Hg is worse
and inverted, with Hg II printed against a 25,000,000 maximum while Hg I's
strongest line in the same range is 12,000.

So neither reading works on its own. Normalizing per spectrum makes every
stage's brightest line 1.0, which is how Ne II rows a neon lamp barely excites
came to outrank Ne I lines and put marks on dark detector. Pooling the raw
numbers across stages would hand a mercury lamp over to Hg II.

`relative_intensity` is therefore the per-spectrum fraction scaled by
`LAMP_STAGE_WEIGHTS` — 1.0 for the neutral stage, 0.1 for singly ionized — a
stated lamp-physics prior rather than a claim about the database. A
low-pressure discharge radiates overwhelmingly in its neutral spectrum, and its
ions appear roughly a decade or two down: weakly, but really, which is why the
curated 20240305 table anchors on Hg II 794.4555 nm and marks it OK. Read the
field as a selection aid; do not take a line ratio off two of them.

Regenerate the caches with `python -m echelle_spectra.tools.nist_cache_refresh`, which writes
down the exact ASD query that produced them.

The Fulcher table is deliberately a position/identification aid. Molecular
deblending, intensity extraction, and population analysis remain in the Fulcher
packages.

## GUI use

The Image tab's Controls dock has five independent checkboxes under
**Known-line overlays**. All start unchecked, so ordinary viewing is unchanged.
Enabling a family draws color-coded vertical markers on both main spectrum
plots and compact labels on the calibrated plot.

The renderer follows the visible wavelength range:

- broad views retain only labels that have enough horizontal room;
- zooming reveals the local molecular structure and progressively weaker NIST
  rows;
- disabling a family removes every marker immediately;
- pan and zoom never alter the selected families.

The screenshots below use a synthetic spectrum solely to make the UI states
and marker behavior reproducible; they are not physics-validation evidence.

Disabled normal use stays visually identical to the established spectrum view:

![Image-tab spectrum plots with every known-line toggle disabled](assets/known-line-overlays/disabled-normal.png)

A broad Balmer/Fulcher view keeps only labels that have enough room:

![Sparse Balmer and Fulcher labels across the full wavelength range](assets/known-line-overlays/balmer-fulcher-broad.png)

Zooming exposes the local line structure:

![Fulcher Q-branch labels after zooming to 600–616 nm](assets/known-line-overlays/fulcher-zoom.png)

![NIST Ne I and Ne II labels after zooming to 580–610 nm](assets/known-line-overlays/ne-zoom.png)

The NIST caches cover `380–810 nm`, which brackets the ~401–802 nm the CMOS
wavelength solution reaches, so a lamp overlay has something to say anywhere on
the detector. They covered only `578–640 nm` — the Fulcher window — through
v0.10, which left every lamp overlay empty outside that interval and, less
visibly, left strong lines just past the edge with no strength annotation and
therefore filtered out of both views. Balmer and Fulcher tables cover their own
listed wavelengths.

## Relationship to line validation

`echelle-validate-lines` now reads the same bundled Balmer and Fulcher records
as the GUI. An external legacy Fulcher matrix can still be supplied with
`--fulcher-table`, but no sibling checkout is required for the default command.
