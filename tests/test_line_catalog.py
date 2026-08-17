"""Tests for the shared provenance-carrying line catalog."""

from __future__ import annotations

import pytest

from echelle_spectra import load_line_table
from echelle_spectra.tools.line_catalog import (
    CATALOG_LINE_FAMILIES,
    CURATED_LINE_TABLE,
    LINE_FAMILIES,
    filter_line_table,
)


def test_all_shared_line_families_load_in_wavelength_order():
    assert LINE_FAMILIES == ("balmer", "fulcher", "thar", "ne", "hg")
    for family in LINE_FAMILIES:
        lines = load_line_table(family)
        assert lines
        assert {line.family for line in lines} == {family}
        assert [line.wavelength_nm for line in lines] == sorted(
            line.wavelength_nm for line in lines
        )
        assert all(line.source_name and line.source_reference for line in lines)


def test_balmer_table_preserves_established_air_values():
    lines = {line.label: line for line in load_line_table("balmer")}
    assert lines["H-alpha"].wavelength_nm == pytest.approx(656.279)
    assert lines["H-beta"].wavelength_nm == pytest.approx(486.135)
    assert lines["H-gamma"].wavelength_nm == pytest.approx(434.047)
    assert lines["H-delta"].wavelength_nm == pytest.approx(410.1734)
    assert {line.wavelength_medium for line in lines.values()} == {"air"}


def test_fulcher_table_matches_authoritative_extractor_resource_with_provenance():
    lines = load_line_table("fulcher")
    assert len(lines) == 42
    by_label = {line.label: line for line in lines}
    assert by_label["Q1(0-0)"].wavelength_nm == pytest.approx(601.8299)
    assert by_label["Q11(2-2)"].wavelength_nm == pytest.approx(639.1882)
    assert by_label["Q9(3-3)"].wavelength_nm == pytest.approx(644.1498)
    assert {line.source_resource for line in lines} == {
        "fulcher_extractor/resources/fulcher_alpha_lines.toml"
    }
    assert all("upstream table does not identify" in line.notes for line in lines)


@pytest.mark.parametrize(
    ("family", "species"),
    [
        ("thar", {"Th I", "Th II", "Ar I", "Ar II"}),
        ("ne", {"Ne I", "Ne II"}),
        ("hg", {"Hg I", "Hg II"}),
    ],
)
def test_atomic_tables_map_cached_nist_rows_without_duplicates(family, species):
    lines = load_line_table(family)
    assert {line.species for line in lines} == species
    assert len({(line.species, line.wavelength_nm) for line in lines}) == len(lines)
    # Every row is a cached NIST export except the handful the curated table
    # names and the cache has no counterpart for; those carry the curated
    # table's own provenance and no invented strength.
    cached = [line for line in lines if line.relative_intensity is not None]
    curated_only = [line for line in lines if line.relative_intensity is None]
    assert all(line.source_name == "NIST Atomic Spectra Database (ASD)" for line in cached)
    assert all(line.source_resource.endswith(".csv") for line in cached)
    assert all(line.curated for line in curated_only)
    assert all(CURATED_LINE_TABLE in line.source_reference for line in curated_only)


def test_filter_line_table_keeps_full_records():
    selected = filter_line_table(
        load_line_table("ne"),
        minimum_nm=600.0,
        maximum_nm=610.0,
        minimum_relative_intensity=0.02,
    )
    assert selected
    assert all(600.0 <= line.wavelength_nm <= 610.0 for line in selected)
    assert all(line.relative_intensity >= 0.02 for line in selected)


def test_every_lamp_row_across_the_instrument_carries_a_strength():
    """The annotation has to reach every row, not only the ones in one window.

    Selection filters on ``relative_intensity``, so a cached row without one is
    a row that cannot be drawn — which is how Ne I 640.2248, the brightest line
    on the owner's frame, went unmarked while sitting two nanometres outside
    the packaged cache.  A curated row is the one exception, and it is not the
    same defect: it is drawn on its vetting rather than on a strength, so the
    field being empty is the honest answer rather than a hole.
    """

    for family in ("ne", "hg", "thar"):
        lines = load_line_table(family)
        for line in lines:
            if line.relative_intensity is None:
                assert line.curated, (family, line.label)
                continue
            assert 0.0 < line.relative_intensity <= 1.0, (family, line.label)
        assert lines[0].wavelength_nm <= 401.0, family
        assert lines[-1].wavelength_nm >= 802.0, family


def test_a_lamp_prefers_its_neutral_lines_to_its_ions():
    """One footing across a lamp's species, so strength ranks the way light does.

    Cached Ne II rows were normalized on their own scale and so reported the
    same 1.0 as the brightest Ne I line, and a selection sorted on strength
    duly preferred ions a neon lamp barely excites.  On the owner's frame that
    put 29 of the Ne boxes on dark detector.
    """

    for family, neutral, ion in (("ne", "Ne I", "Ne II"), ("hg", "Hg I", "Hg II")):
        lines = load_line_table(family)
        strongest = {
            species: max(
                line.relative_intensity
                for line in lines
                if line.species == species and line.relative_intensity is not None
            )
            for species in (neutral, ion)
        }
        assert strongest[neutral] > strongest[ion], (family, strongest)

    # An ion is ranked below the neutral stage, not deleted: the curated
    # 20240305 table anchors on Hg II 794.4555 and marks it OK.
    hg_ii = [
        line
        for line in load_line_table("hg")
        if line.species == "Hg II" and abs(line.wavelength_nm - 794.4555) < 0.01
    ]
    assert len(hg_ii) == 1
    assert hg_ii[0].relative_intensity > 0.0


def test_unknown_line_family_is_rejected():
    with pytest.raises(ValueError, match="known families"):
        load_line_table("xenon")


def test_xenon_is_a_served_catalog_family_without_an_overlay_toggle():
    """The bench serves Xe lines; the viewer's five-hue palette is untouched.

    Which lamps the catalog can answer for and how many overlay colours fit
    between the two spectrum curves are separate questions, and the second one
    is not allowed to decide the first.
    """

    assert "xe" not in LINE_FAMILIES
    assert "xe" in CATALOG_LINE_FAMILIES
    assert set(LINE_FAMILIES) < set(CATALOG_LINE_FAMILIES)
    assert load_line_table("xe")


def test_the_packaged_xenon_cache_spans_the_instrument_with_a_strength_on_every_row():
    """Xe I and Xe II, 380-810 nm, every row rankable — the wall's other half.

    A cache that stopped short of the detector's range, or a row with no
    strength, would put the xenon lamp back where it was: assignable in name
    and useless in the panel that has to choose which lines to show.
    """

    lines = load_line_table("xe")
    species = {line.species for line in lines}
    assert species == {"Xe I", "Xe II"}
    assert all(line.relative_intensity is not None for line in lines)
    assert all(line.relative_intensity > 0.0 for line in lines)
    wavelengths = [line.wavelength_nm for line in lines]
    assert min(wavelengths) < 385.0
    assert max(wavelengths) > 800.0
    # The lamp-stage prior applies to xenon exactly as to neon and mercury:
    # the ion sits a decade below the neutral rather than beside it.
    strongest = {
        name: max(
            line.relative_intensity for line in lines if line.species == name
        )
        for name in species
    }
    assert strongest["Xe I"] > strongest["Xe II"]
    assert strongest["Xe II"] == pytest.approx(0.1 * strongest["Xe I"])
    # Cached from the packaged files, named so their span is legible.
    assert {
        line.source_resource.rsplit("/", 1)[-1] for line in lines
    } == {"nist_xe_i_380_810.csv", "nist_xe_ii_380_810.csv"}


def test_no_curated_row_vouches_for_a_xenon_line():
    """Nobody has vetted a xenon line on this instrument, and the table says so."""

    assert not any(line.curated for line in load_line_table("xe"))
