"""Tests for the shared provenance-carrying line catalog."""

from __future__ import annotations

import pytest

from echelle_spectra import load_line_table
from echelle_spectra.tools.line_catalog import LINE_FAMILIES, filter_line_table


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
    assert all(line.source_name == "NIST Atomic Spectra Database (ASD)" for line in lines)
    assert all(line.source_resource.endswith(".csv") for line in lines)
    assert all(line.relative_intensity is not None for line in lines)


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

    Selection filters on ``relative_intensity``, so a row without one is a row
    that cannot be drawn — which is how Ne I 640.2248, the brightest line on
    the owner's frame, went unmarked while sitting two nanometres outside the
    packaged cache.
    """

    for family in ("ne", "hg", "thar"):
        lines = load_line_table(family)
        assert all(line.relative_intensity is not None for line in lines), family
        assert all(0.0 < line.relative_intensity <= 1.0 for line in lines), family
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
                if line.species == species
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
