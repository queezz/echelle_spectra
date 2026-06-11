from pathlib import Path

import numpy as np

from echelle_spectra.tools.nist_lamp_calibration import (
    default_nist_cache_dir,
    lamp_species,
    load_nist_asd_exports,
    normalize_species_key,
    resolve_cached_line_lists,
)


def test_normalize_species_key_accepts_common_spellings():
    assert normalize_species_key("Th I") == "ThI"
    assert normalize_species_key("th_i") == "ThI"
    assert normalize_species_key("NeII") == "NeII"
    assert normalize_species_key("HII") == "HII"


def test_lamp_species_expands_presets_without_duplicates():
    assert lamp_species(["ThAr"]) == ("ThI", "ThII", "ArI", "ArII")
    assert lamp_species(["Ne", "ne"]) == ("NeI", "NeII")


def test_resolve_cached_line_lists_from_lamp_preset(tmp_path: Path):
    th_i = tmp_path / "nist_th_i_578_620.csv"
    th_ii = tmp_path / "nist_th_ii_578_620.csv"
    ar_i = tmp_path / "nist_ar_i_578_620_normalized.csv"
    ar_ii = tmp_path / "nist_arii_578_620.csv"
    for path in (th_i, th_ii, ar_i, ar_ii):
        path.write_text("obs_wl_air(nm),intens\n600.0,1\n", encoding="utf-8")

    resolved = resolve_cached_line_lists(lamps=["thar"], cache_dir=tmp_path)

    assert resolved == {
        "ThI": th_i.resolve(),
        "ThII": th_ii.resolve(),
        "ArI": ar_i.resolve(),
        "ArII": ar_ii.resolve(),
    }


def test_resolve_cached_line_lists_explicit_overrides_cache(tmp_path: Path):
    cache = tmp_path / "cache"
    cache.mkdir()
    auto = cache / "nist_ne_i_580_620.csv"
    explicit = tmp_path / "custom_ne_i.csv"
    ne_ii = cache / "nist_ne_ii_580_620.csv"
    for path in (auto, explicit, ne_ii):
        path.write_text("obs_wl_air(nm),intens\n600.0,1\n", encoding="utf-8")

    resolved = resolve_cached_line_lists(
        lamps=["ne"],
        cache_dir=cache,
        explicit={"Ne I": explicit},
    )

    assert resolved["NeI"] == explicit.resolve()
    assert resolved["NeII"] == ne_ii.resolve()


def test_default_nist_cache_contains_curated_th_ar_exports():
    resolved = resolve_cached_line_lists(lamps=["thar"], cache_dir=default_nist_cache_dir())

    assert set(resolved) == {"ThI", "ThII", "ArI", "ArII"}
    assert {path.name for path in resolved.values()} == {
        "nist_ar_i_578_640.csv",
        "nist_ar_ii_578_640.csv",
        "nist_th_i_578_640.csv",
        "nist_th_ii_578_640.csv",
    }


def test_default_nist_cache_contains_curated_hg_exports():
    resolved = resolve_cached_line_lists(lamps=["hg"], cache_dir=default_nist_cache_dir())

    assert set(resolved) == {"HgI", "HgII"}
    assert {path.name for path in resolved.values()} == {
        "nist_hg_i_578_640.csv",
        "nist_hg_ii_578_640.csv",
    }


def test_default_nist_cache_contains_curated_ne_exports():
    resolved = resolve_cached_line_lists(lamps=["ne"], cache_dir=default_nist_cache_dir())

    assert set(resolved) == {"NeI", "NeII"}
    assert {path.name for path in resolved.values()} == {
        "nist_ne_i_578_640.csv",
        "nist_ne_ii_578_640.csv",
    }


def test_load_nist_asd_exports_reads_csv_and_tsv(tmp_path: Path):
    csv_path = tmp_path / "nist_hg_i.csv"
    csv_path.write_text(
        "obs_wl_air(nm),ritz_wl_air(nm),intens,Aki(s^-1)\n"
        '="579.0663",="579.0663",100,\n'
        ',"580.0",,1.0e+05\n',
        encoding="utf-8",
    )
    tsv_path = tmp_path / "nist_ne_ii.tsv.csv"
    tsv_path.write_text(
        "obs_wl_air(nm)\tritz_wl_air(nm)\tintens\tAki(s^-1)\n" '"588.189"\t"588.189"\t"50"\t""\n',
        encoding="utf-8",
    )

    lines = load_nist_asd_exports(
        {"HgI": csv_path, "Ne II": tsv_path},
        min_wavelength_nm=578.0,
        max_wavelength_nm=590.0,
    )

    assert set(lines["species"]) == {"HgI", "NeII"}
    assert np.isclose(lines.loc[lines["species"] == "HgI", "weight"].max(), 1.0)
    assert np.isclose(lines.loc[lines["species"] == "NeII", "wavelength_nm"].iloc[0], 588.189)
