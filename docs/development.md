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

## Historical HTML docs

The `html-docs` branch contains the original static HTML documentation pages (`index.html` and `band_data.html`). That branch is **historical** — do not delete it during or after this migration. The new MkDocs site supersedes it but the old pages remain available for reference.
