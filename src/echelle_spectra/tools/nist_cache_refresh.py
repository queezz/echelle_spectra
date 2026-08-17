"""Re-extract the packaged NIST ASD line caches from physics.nist.gov.

The caches under ``src/echelle_spectra/resources/nist_asd_cache`` were first
produced by hand from the NIST ASD "Lines" web form, and nothing in the tree
recorded which query produced them.  That is how they came to span only
580--638 nm -- the Fulcher window -- on an instrument that reaches 400--800 nm,
with no way to widen them short of repeating a forgotten form.

This script *is* that query, written down.  Its parameter block was recovered
by replaying the ASD form until the response reproduced the packaged files
byte for byte, so re-running it over the original ``580--640`` window
regenerates the shipped caches exactly and re-running it over the instrument's
own range widens them without changing anything else about their format.

This is the one module in the package that reaches the network, and it is the
only one that writes to ``resources/``.  Everything else reads the cache, which
is what keeps a calibration run reproducible; refreshing it is a deliberate,
dated act somebody performs — before travel, or when the instrument's range
moves and the packaged spans no longer cover it.

It deliberately has no ``echelle-*`` console script.  The others are operator
commands; this one fetches from the internet and rewrites packaged data, and it
does not belong one tab-completion away from ``echelle-align`` on a lab
machine.  Run it as a module::

    python -m echelle_spectra.tools.nist_cache_refresh                 # every lamp
    python -m echelle_spectra.tools.nist_cache_refresh --species NeI NeII
    python -m echelle_spectra.tools.nist_cache_refresh --low-nm 580 --high-nm 640 --dry-run

Every written file is named ``nist_<element>_<ion>_<low>_<high>.csv`` so the
span a file covers is legible without opening it, and so a stale narrower file
is visibly a different file rather than a silent overwrite.
"""

from __future__ import annotations

import argparse
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from .nist_lamp_calibration import (
    COMMON_NIST_SPECIES,
    default_nist_cache_dir,
    normalize_species_key,
)

ASD_LINES_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"

#: The instrument reaches roughly 401.5--801.9 nm on the 20250926 CMOS
#: calibration, and the curated wavelength tables carry rows out to 800.7 nm.
#: These bounds cover both with margin so a re-alignment that shifts the
#: solution by a nanometre does not walk a row off the end of the cache.
DEFAULT_LOW_NM = 380.0
DEFAULT_HIGH_NM = 810.0

#: Species the package ships a cache for, in the order the files are written.
PACKAGED_SPECIES = (
    "NeI",
    "NeII",
    "HgI",
    "HgII",
    "ArI",
    "ArII",
    "ThI",
    "ThII",
    "XeI",
    "XeII",
)

#: The exact ASD form settings that produce the packaged six-column export:
#: ``obs_wl_air(nm) ritz_wl_air(nm) intens Aki(s^-1) Acc Type``.
#:
#: ASD reads its column checkboxes as *present or absent* -- passing
#: ``conf_out=0`` is not "off", it is an "Invalid Column Setting" error -- so
#: the columns we do not want are simply missing from this list.
_QUERY_TEMPLATE: tuple[tuple[str, str], ...] = (
    ("limits_type", "0"),
    ("unit", "1"),  # nanometres
    ("submit", "Retrieve Data"),
    ("de", "0"),
    ("format", "3"),  # tab-delimited
    ("line_out", "0"),  # every line, not only those with transition data
    ("en_unit", "0"),
    ("output", "0"),  # one page, not paginated
    ("page_size", "15"),
    ("show_obs_wl", "1"),
    ("show_calc_wl", "1"),  # Ritz, used when no observed wavelength exists
    ("order_out", "0"),
    ("max_low_enrg", ""),
    ("show_av", "2"),  # air wavelengths
    ("max_upp_enrg", ""),
    ("tsb_value", "0"),
    ("min_str", ""),
    ("A_out", "0"),
    ("intens_out", "on"),  # the relative-intensity column this cache exists for
    ("max_str", ""),
    ("allowed_out", "1"),
    ("forbid_out", "1"),
    ("min_accur", ""),
    ("min_intens", ""),
)

USER_AGENT = "echelle_spectra NIST ASD cache refresh (https://github.com/queezz)"


def _span_token(value: float) -> str:
    """Render a bound the way the packaged filenames do: a bare integer."""

    return str(int(round(value)))


def cache_filename(species_key: str, low_nm: float, high_nm: float) -> str:
    """Return the cache filename for one species over one span."""

    species = COMMON_NIST_SPECIES[species_key]
    return (
        f"nist_{species.element.lower()}_{species.ion.lower()}"
        f"_{_span_token(low_nm)}_{_span_token(high_nm)}.csv"
    )


def build_query(species_key: str, low_nm: float, high_nm: float) -> str:
    """Return the full ASD URL for one species over one wavelength span."""

    species = COMMON_NIST_SPECIES[species_key]
    params = [
        ("spectra", species.nist_name),
        ("low_w", _span_token(low_nm)),
        ("upp_w", _span_token(high_nm)),
        *_QUERY_TEMPLATE,
    ]
    return f"{ASD_LINES_URL}?{urllib.parse.urlencode(params)}"


def fetch_species(species_key: str, low_nm: float, high_nm: float, *, timeout: float) -> str:
    """Fetch one ASD export, refusing anything that is not the tabular answer."""

    url = build_query(species_key, low_nm, high_nm)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read().decode("utf-8", errors="replace")
    if not payload.startswith("obs_wl_air(nm)"):
        # ASD answers a malformed query with an HTML error page, and writing
        # that into the cache would be a silent corruption of packaged data.
        head = " ".join(payload.split())[:200]
        raise RuntimeError(
            f"NIST ASD did not return a tab-delimited export for "
            f"{COMMON_NIST_SPECIES[species_key].nist_name}: {head}"
        )
    return payload


def refresh(
    species_keys: list[str],
    *,
    low_nm: float,
    high_nm: float,
    cache_dir: Path,
    dry_run: bool = False,
    prune: bool = True,
    timeout: float = 180.0,
    pause_s: float = 1.0,
) -> list[Path]:
    """Fetch and write one cache file per species; return the paths written."""

    cache_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for index, key in enumerate(species_keys):
        name = cache_filename(key, low_nm, high_nm)
        target = cache_dir / name
        if index:
            time.sleep(pause_s)
        payload = fetch_species(key, low_nm, high_nm, timeout=timeout)
        rows = max(0, len(payload.splitlines()) - 1)
        print(f"{COMMON_NIST_SPECIES[key].nist_name:6s} {rows:6d} rows -> {name}")
        if dry_run:
            continue
        # Newline="" keeps the LF endings NIST sends; the packaged caches are
        # LF files and a CRLF rewrite would be a whole-file diff saying nothing.
        target.write_text(payload, encoding="utf-8", newline="")
        written.append(target)
        if prune:
            _prune_superseded(cache_dir, key, target)
    return written


def _prune_superseded(cache_dir: Path, species_key: str, keep: Path) -> None:
    """Delete narrower caches for the same species that this one replaces.

    Leaving them behind is not harmless: cache discovery matches on the species
    tokens in the filename and would then have two files to choose between,
    with the span deciding nothing.
    """

    species = COMMON_NIST_SPECIES[species_key]
    prefix = f"nist_{species.element.lower()}_{species.ion.lower()}_"
    for path in sorted(cache_dir.glob(f"{prefix}*.csv")):
        if path.resolve() != keep.resolve():
            print(f"  superseded: removing {path.name}")
            path.unlink()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m echelle_spectra.tools.nist_cache_refresh",
        description="Re-extract packaged NIST ASD lamp-line caches over a wavelength span.",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=list(PACKAGED_SPECIES),
        help="Species to refresh, e.g. NeI NeII. Default: every packaged species.",
    )
    parser.add_argument("--low-nm", type=float, default=DEFAULT_LOW_NM)
    parser.add_argument("--high-nm", type=float, default=DEFAULT_HIGH_NM)
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Destination directory. Default: the package-bundled cache.",
    )
    parser.add_argument(
        "--keep-superseded",
        action="store_true",
        help="Leave narrower caches for the same species in place.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch and report row counts without writing anything.",
    )
    parser.add_argument("--timeout", type=float, default=180.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    keys = [normalize_species_key(item) for item in args.species]
    cache_dir = Path(args.cache_dir) if args.cache_dir else default_nist_cache_dir()
    if args.low_nm >= args.high_nm:
        raise SystemExit("--low-nm must be below --high-nm")
    print(f"NIST ASD {args.low_nm:g}-{args.high_nm:g} nm -> {cache_dir}")
    try:
        refresh(
            keys,
            low_nm=args.low_nm,
            high_nm=args.high_nm,
            cache_dir=cache_dir,
            dry_run=args.dry_run,
            prune=not args.keep_superseded,
            timeout=args.timeout,
        )
    except (urllib.error.URLError, RuntimeError) as exc:
        raise SystemExit(f"NIST ASD refresh failed: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
