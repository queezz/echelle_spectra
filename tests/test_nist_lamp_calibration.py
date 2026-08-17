from pathlib import Path

import numpy as np

from echelle_spectra.tools.nist_lamp_calibration import (
    LAMP_STAGE_WEIGHTS,
    default_nist_cache_dir,
    lamp_species,
    load_nist_asd_exports,
    normalize_species_key,
    resolve_cached_line_lists,
    stage_weight_for_species,
)

#: The span the 20250926 CMOS wavelength solution reaches, rounded inward so a
#: re-alignment of a few pixels cannot make this test lie either way.
INSTRUMENT_LOW_NM = 401.0
INSTRUMENT_HIGH_NM = 802.0


def _species_from_filename(path: Path) -> str:
    """``nist_ne_ii_380_810.csv`` -> ``NeII``."""

    element, ion = path.stem.split("_")[1:3]
    return f"{element.capitalize()}{ion.upper()}"


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
        "nist_ar_i_380_810.csv",
        "nist_ar_ii_380_810.csv",
        "nist_th_i_380_810.csv",
        "nist_th_ii_380_810.csv",
    }


def test_default_nist_cache_contains_curated_hg_exports():
    resolved = resolve_cached_line_lists(lamps=["hg"], cache_dir=default_nist_cache_dir())

    assert set(resolved) == {"HgI", "HgII"}
    assert {path.name for path in resolved.values()} == {
        "nist_hg_i_380_810.csv",
        "nist_hg_ii_380_810.csv",
    }


def test_default_nist_cache_contains_curated_ne_exports():
    resolved = resolve_cached_line_lists(lamps=["ne"], cache_dir=default_nist_cache_dir())

    assert set(resolved) == {"NeI", "NeII"}
    assert {path.name for path in resolved.values()} == {
        "nist_ne_i_380_810.csv",
        "nist_ne_ii_380_810.csv",
    }


def test_every_packaged_cache_covers_the_whole_instrument_range():
    """No packaged cache may be narrower than the light the instrument sees.

    The caches shipped covering 580.4--638.3 nm, the Fulcher window, on an
    instrument whose 20250926 CMOS solution runs 401.5--801.9 nm and whose
    curated wavelength tables carry rows out to 800.7 nm.  Every consequence of
    that followed silently: a lamp overlay simply had nothing to draw outside
    the window, and the brightest line on the owner's Ne frame — Ne I 640.2248,
    two nanometres past the edge — carried no strength annotation and was
    filtered out of both views for lack of one.
    """

    cache_dir = default_nist_cache_dir()
    files = sorted(cache_dir.glob("*.csv"))
    assert files, f"no packaged NIST caches in {cache_dir}"
    for path in files:
        low, high = (float(token) for token in path.stem.split("_")[3:5])
        # The span is asserted on the *extraction*, not on where lines happen
        # to fall: Hg II has nothing above 794.4555 nm and that is mercury's
        # answer, not a clipped file.
        assert low <= INSTRUMENT_LOW_NM, path.name
        assert high >= INSTRUMENT_HIGH_NM, path.name

        table = load_nist_asd_exports(
            {_species_from_filename(path): path},
            min_wavelength_nm=0.0,
            max_wavelength_nm=10_000.0,
        )
        assert not table.empty, path.name
        assert table["wavelength_nm"].min() >= low, path.name
        assert table["wavelength_nm"].max() <= high, path.name
        # And every file really does reach past the old Fulcher-window edges,
        # which is the clipping this exists to catch.
        assert table["wavelength_nm"].min() < 578.0, path.name
        assert table["wavelength_nm"].max() > 640.0, path.name


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
    assert np.isclose(lines.loc[lines["species"] == "HgI", "spectrum_weight"].max(), 1.0)
    assert np.isclose(lines.loc[lines["species"] == "HgI", "weight"].max(), 1.0)
    assert np.isclose(lines.loc[lines["species"] == "NeII", "wavelength_nm"].iloc[0], 588.189)


def test_the_strongest_line_of_every_spectrum_is_no_longer_worth_the_same():
    """Each stage's brightest line used to read 1.0, whatever the lamp does.

    ``spectrum_weight`` is still that per-spectrum reading, because that is all
    NIST's own number supports.  ``weight`` is the one selections rank by, and
    it carries the ionization stage, so a Ne II line cannot present itself as
    the equal of the brightest Ne I line in the frame.
    """

    cache = default_nist_cache_dir()
    lines = load_nist_asd_exports(
        {
            "NeI": cache / "nist_ne_i_380_810.csv",
            "NeII": cache / "nist_ne_ii_380_810.csv",
        },
        min_wavelength_nm=401.0,
        max_wavelength_nm=802.0,
    )
    peaks = lines.groupby("species")[["spectrum_weight", "weight"]].max()
    assert np.isclose(peaks.loc["NeI", "spectrum_weight"], 1.0)
    assert np.isclose(peaks.loc["NeII", "spectrum_weight"], 1.0)
    assert np.isclose(peaks.loc["NeI", "weight"], 1.0)
    assert np.isclose(peaks.loc["NeII", "weight"], LAMP_STAGE_WEIGHTS[2])

    # And the neutral stage now wins on the real number a selection sorts by:
    # Ne I 640.2248, the brightest blob on the owner's frame, outranks every
    # cached Ne II row rather than tying with the strongest of them.
    brightest_ion = lines.loc[lines["species"] == "NeII", "weight"].max()
    line_640 = lines.loc[(lines["wavelength_nm"] - 640.2248).abs() < 0.01]
    assert len(line_640) == 1
    assert float(line_640["weight"].iloc[0]) > brightest_ion


def test_stage_weight_places_ions_below_neutrals_without_hiding_them():
    """The curated 20240305 table anchors on Hg II 794.4555, so 'ion' is not 'gone'."""

    assert stage_weight_for_species("NeI") == 1.0
    assert stage_weight_for_species("Hg II") == LAMP_STAGE_WEIGHTS[2]
    assert 0.0 < LAMP_STAGE_WEIGHTS[2] < LAMP_STAGE_WEIGHTS[1]
    assert LAMP_STAGE_WEIGHTS[3] < LAMP_STAGE_WEIGHTS[2]
