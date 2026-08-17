# Developer Notes

## Local docs setup

To build and preview the documentation locally, set up a virtual environment with the required MkDocs dependencies.

=== "Unix / macOS"

    ```bash
    python -m venv .venv
    source .venv/bin/activate
    pip install mkdocs mkdocs-material pymdown-extensions mdx_truly_sane_lists \
        "mkdocstrings[python]>=0.25" mkdocs-glightbox
    ```

=== "Windows PowerShell"

    ```powershell
    python -m venv .venv
    .venv\Scripts\activate
    pip install mkdocs mkdocs-material pymdown-extensions mdx_truly_sane_lists `
        "mkdocstrings[python]>=0.25" mkdocs-glightbox
    ```

### Serve locally

```bash
mkdocs serve
```

Opens a live-reloading preview at [http://127.0.0.1:8000](http://127.0.0.1:8000).

### Build

```bash
mkdocs build --strict
```

Output is written to `./site/`. The `--strict` flag treats warnings as errors — keep the build clean before pushing.

---

## Deployment

GitHub Actions deploys the site automatically on every push to `main` or `master`. The workflow is defined in `.github/workflows/gh-pages.yml` and uses `peaceiris/actions-gh-pages@v3` to push the built site to the `gh-pages` branch.

The live site is served at: <https://queezz.github.io/echelle_spectra>

---

## Portable kit and release builds

The pinned NIFS kit manifest, external-only assembly commands, three-command
bare-machine route, recovery behavior, and byte-identical wheel/sdist contract
are documented in [Portable NIFS kit](portable-kit.md) and the root
[`README-KIT.md`](https://github.com/queezz/echelle_spectra/blob/master/README-KIT.md).

Do not build kit payloads or release archives in this checkout. Export a clean
source copy into OS-local scratch and use `scripts.reproducible_build` there;
assemble kits into a separate external destination with `scripts.nifs_kit`.

---

## Refreshing the packaged NIST caches

The lamp overlays and the bench's expected-line annotations read cached NIST ASD
exports from `src/echelle_spectra/resources/nist_asd_cache/`. A calibration run
never queries NIST; refreshing the cache is a separate, deliberate act:

```bash
python -m echelle_spectra.tools.nist_cache_refresh                 # every packaged lamp
python -m echelle_spectra.tools.nist_cache_refresh --species NeI NeII --dry-run
```

Run it **before travel**, so an offline instrument carries current line data,
and **when the instrument's range moves** — the packaged spans are 380–810 nm
against a detector reaching ~401–802 nm, and a re-alignment that walks curated
rows past a cache edge leaves those lines with no strength annotation, which
silently removes them from both overlay views.

The module carries the exact ASD query as its parameter block, recovered by
replaying the web form until the response reproduced the packaged files byte for
byte. It refuses to write anything that is not a tab-delimited export, so a
malformed query cannot cache an HTML error page over packaged data. Widening a
span writes new `nist_<element>_<ion>_<low>_<high>.csv` files and removes the
narrower ones they supersede. It has no `echelle-*` console script on purpose:
it reaches the network and rewrites packaged resources, which is maintenance,
not an operator command.

---

## Historical HTML docs

The `html-docs` branch contains the original static HTML documentation pages (`index.html` and `band_data.html`). That branch is **historical** — do not delete it during or after this migration. The new MkDocs site supersedes it but the old pages remain available for reference.
