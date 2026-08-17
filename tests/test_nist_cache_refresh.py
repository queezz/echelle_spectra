"""The packaged NIST cache refresher — offline behaviour only.

Nothing here reaches the network. What is worth pinning is the shape of the
query, the naming that makes a cache's span legible, and above all the refusal
to write a response that is not a tab-delimited export: this is the one module
that overwrites packaged resources, and ASD answers a malformed query with an
HTML error page that would otherwise be cached over real data.
"""

from __future__ import annotations

import urllib.parse
from pathlib import Path

import pytest

from echelle_spectra.tools import nist_cache_refresh as refresh


def test_the_filename_says_which_span_the_file_covers():
    assert refresh.cache_filename("NeI", 380.0, 810.0) == "nist_ne_i_380_810.csv"
    assert refresh.cache_filename("ThII", 380.0, 810.0) == "nist_th_ii_380_810.csv"


def test_the_query_asks_for_the_columns_the_loader_reads():
    query = urllib.parse.parse_qs(
        urllib.parse.urlparse(refresh.build_query("NeII", 380.0, 810.0)).query
    )
    assert query["spectra"] == ["Ne II"]
    assert query["low_w"] == ["380"]
    assert query["upp_w"] == ["810"]
    assert query["unit"] == ["1"]  # nanometres
    assert query["format"] == ["3"]  # tab-delimited
    assert query["show_av"] == ["2"]  # air wavelengths
    assert query["intens_out"] == ["on"]  # the column this cache exists for
    assert query["show_obs_wl"] == ["1"]
    assert query["show_calc_wl"] == ["1"]
    # ASD reads its column checkboxes as present-or-absent: sending "0" is not
    # "off", it is an Invalid Column Setting error. So the columns the packaged
    # format does not carry must be missing rather than disabled.
    for absent in ("conf_out", "term_out", "enrg_out", "J_out", "unc_out", "bibrefs"):
        assert absent not in query


def test_an_html_error_page_is_never_written_over_packaged_data(monkeypatch):
    class _Response:
        def __init__(self, payload: str):
            self._payload = payload

        def read(self) -> bytes:
            return self._payload.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *_exc):
            return False

    monkeypatch.setattr(
        refresh.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: _Response("<html><title>NIST ASD : Input Error</title>"),
    )
    with pytest.raises(RuntimeError, match="did not return a tab-delimited export"):
        refresh.fetch_species("NeI", 380.0, 810.0, timeout=1.0)


def test_a_widened_cache_removes_the_narrower_one_it_supersedes(tmp_path: Path):
    """Two files for one species would leave discovery picking by filename luck."""

    stale = tmp_path / "nist_ne_i_578_640.csv"
    stale.write_text("obs_wl_air(nm)\tintens\n", encoding="utf-8")
    fresh = tmp_path / "nist_ne_i_380_810.csv"
    fresh.write_text("obs_wl_air(nm)\tintens\n", encoding="utf-8")
    other = tmp_path / "nist_ne_ii_578_640.csv"
    other.write_text("obs_wl_air(nm)\tintens\n", encoding="utf-8")

    refresh._prune_superseded(tmp_path, "NeI", fresh)

    assert fresh.exists()
    assert not stale.exists()
    assert other.exists(), "another species' cache is not this one's to remove"


def test_every_packaged_species_is_one_the_loader_knows():
    for key in refresh.PACKAGED_SPECIES:
        assert refresh.normalize_species_key(key) == key
    assert refresh.DEFAULT_LOW_NM < 401.0
    assert refresh.DEFAULT_HIGH_NM > 802.0
