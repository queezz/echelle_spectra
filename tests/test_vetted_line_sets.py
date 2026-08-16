"""Packet F19, first rider — the trusted set is the BH-paper set.

An ``OK`` beside a row in a wavelength table is two letters; what makes it
mean something is the work behind it — the alignment vetted during the BH
paper's own calibration, tried and tested for its Balmer and Fulcher
analysis.  So the auto-anchor may not simply trust whatever OK marks a
loaded table happens to carry: it states whose vetting they are, or states
that they carry none, and the snapshot it anchors records which set stood
behind it.  These tests pin that lineage, and pin that it is read from the
tables rather than assumed from their names.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.calibration_bench import BenchFrame, CalibrationBenchSession
from echelle_spectra.calibration_bench_gui import CalibrationBenchWindow
from echelle_spectra.calibration_campaign import lamp_reference_set
from echelle_spectra.tools.calibration_alignment import (
    BH_PAPER_WAVELENGTH_TABLE,
    CalibrationTableLine,
    read_table_metadata,
    table_vetting,
)

_CALIBRATION_FILES = (
    Path(__file__).parents[1]
    / "src"
    / "echelle_spectra"
    / "resources"
    / "calibration_files"
)
_BASE_TABLE = _CALIBRATION_FILES / BH_PAPER_WAVELENGTH_TABLE
_ALIGNED_TABLE = (
    _CALIBRATION_FILES
    / "alignments"
    / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
)


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


class TestTheLineageIsReadFromTheTables:
    def test_the_curated_table_is_the_vetted_set_itself(self):
        vetting = table_vetting(_BASE_TABLE)

        assert vetting.is_vetted
        assert vetting.vetted_set == "BH paper"
        assert vetting.lineage == (BH_PAPER_WAVELENGTH_TABLE,)
        assert "Balmer and Fulcher" in vetting.description

    def test_the_aligned_table_inherits_the_vetting_it_was_derived_from(self):
        # This is the rider's own wording: the 20240305 table's OK marks "as
        # carried through the 20250926 alignment".  The bench opens on this
        # table by default, so the inheritance is the common case, not a corner.
        vetting = table_vetting(_ALIGNED_TABLE)

        assert vetting.is_vetted
        assert vetting.vetted_set == "BH paper"
        assert vetting.vetted_table == BH_PAPER_WAVELENGTH_TABLE
        assert vetting.lineage == (
            _ALIGNED_TABLE.name,
            BH_PAPER_WAVELENGTH_TABLE,
        )
        assert f"inherited from {BH_PAPER_WAVELENGTH_TABLE}" in vetting.message

    def test_the_lineage_comes_from_the_header_not_the_filename(self, tmp_path):
        # Renaming a file must not launder its pedigree in either direction.
        renamed = tmp_path / "something_entirely_else.txt"
        renamed.write_text(_ALIGNED_TABLE.read_text(), encoding="utf-8")

        vetting = table_vetting(renamed, [_CALIBRATION_FILES])

        assert vetting.is_vetted
        assert vetting.table == "something_entirely_else.txt"
        assert vetting.vetted_table == BH_PAPER_WAVELENGTH_TABLE

    def test_a_table_named_like_the_vetted_one_still_needs_the_lineage(
        self, tmp_path
    ):
        impostor = tmp_path / "Th_wavelength_CMOS_20240305_aligned_to_someday.txt"
        impostor.write_text(
            "# Adjusted wavelength calibration lookup table\n"
            "# Base wavelength file: invented_table.txt\n"
            "0 10 20 15.0 650.0 NeI  # ok\n",
            encoding="utf-8",
        )

        vetting = table_vetting(impostor)

        assert not vetting.is_vetted
        assert vetting.lineage == (impostor.name, "invented_table.txt")

    def test_an_unrelated_table_borrows_nobodys_authority(self):
        vetting = table_vetting(_CALIBRATION_FILES / "wavelength_fujii.txt")

        assert not vetting.is_vetted
        assert vetting.vetted_set == ""
        assert "carry no recorded vetting" in vetting.message
        assert "not the BH paper's" in vetting.message

    def test_a_lineage_that_loops_terminates(self, tmp_path):
        looping = tmp_path / "loop.txt"
        looping.write_text(
            "# Base wavelength file: loop.txt\n0 10 20 15.0 650.0 NeI  # ok\n",
            encoding="utf-8",
        )

        vetting = table_vetting(looping)

        assert not vetting.is_vetted
        assert vetting.lineage == ("loop.txt",)


class TestOnlyTheLeadingHeaderIsMetadata:
    def test_commented_out_data_rows_are_not_read_as_metadata(self):
        # The curated table disables rows by commenting them out, and some of
        # those carry colons in their notes.  A row is a row however much it
        # looks like a header once a "#" is in front of it.
        metadata = read_table_metadata(_BASE_TABLE)

        assert "lamps used" in metadata
        assert not any(key.strip().isdigit() for key in metadata)

    def test_the_aligned_tables_provenance_survives_the_round_trip(self):
        metadata = read_table_metadata(_ALIGNED_TABLE)

        assert metadata["base wavelength file"] == BH_PAPER_WAVELENGTH_TABLE
        assert metadata["alignment dataset"] == "20250926"


_COLUMNS = 120


def _pattern() -> np.ndarray:
    pattern = np.empty((_COLUMNS, 2), dtype=float)
    pattern[:, 0] = 12
    pattern[:, 1] = 30
    return pattern


def _rows() -> tuple[CalibrationTableLine, ...]:
    return (
        CalibrationTableLine(0, 25, 35, 30.0, 585.249, "NeI", "ok"),
        CalibrationTableLine(1, 65, 75, 70.0, 600.062, "NeI", "ok"),
    )


def _frame(path: Path) -> BenchFrame:
    x = np.arange(_COLUMNS, dtype=float)
    spectra = tuple(
        10.0 + 500.0 * np.exp(-0.5 * ((x - center) / 1.6) ** 2)
        for center in (30.0, 70.0)
    )
    images = np.zeros((1, 44, _COLUMNS), dtype=float)
    pattern = _pattern()
    for order, spectrum in enumerate(spectra):
        for column, value in enumerate(spectrum):
            images[0, int(pattern[column, order]), column] = value
    return BenchFrame(path, images, images[0], spectra, {"ExposureTime": 0.1})


def _window(tmp_path: Path, table: Path | None):
    session = CalibrationBenchSession(
        _pattern(),
        _rows(),
        minimum_snr=3.0,
        vetting=None if table is None else table_vetting(table, [_CALIBRATION_FILES]),
    )
    session.accept_frame(_frame(tmp_path / "Ne.sif"))
    session.use_lamp_reference(lamp_reference_set("Ne", session.lines))
    window = CalibrationBenchWindow(session, start_timer=False)
    window.refresh()
    return window


class TestTheBenchStatesWhoseVettingItIsUsing:
    def test_a_vetted_table_is_named_before_the_pass_is_pressed(
        self, qt_app, tmp_path
    ):
        window = _window(tmp_path, _ALIGNED_TABLE)

        assert "carrying the BH paper vetting" in window.auto_anchor_value.text()
        window.close()

    def test_an_unvetted_table_says_so_instead_of_staying_quiet(
        self, qt_app, tmp_path
    ):
        window = _window(tmp_path, _CALIBRATION_FILES / "wavelength_fujii.txt")

        caption = window.auto_anchor_value.text()
        assert "carry no recorded vetting" in caption
        assert "not the BH paper's" in caption
        window.close()

    def test_the_pass_reports_which_set_it_anchored_against(self, qt_app, tmp_path):
        window = _window(tmp_path, _ALIGNED_TABLE)

        window.auto_anchor_button.click()
        qt_app.processEvents()

        assert "against the BH paper vetted set" in window.message_value.text()
        window.close()

    def test_an_unvetted_pass_does_not_borrow_the_papers_authority(
        self, qt_app, tmp_path
    ):
        window = _window(tmp_path, _CALIBRATION_FILES / "wavelength_fujii.txt")

        window.auto_anchor_button.click()
        qt_app.processEvents()

        message = window.message_value.text()
        assert "against OK marks carrying no recorded vetting" in message
        assert "BH paper" not in message
        window.close()

    def test_a_bench_that_was_never_told_claims_nothing_either_way(
        self, qt_app, tmp_path
    ):
        window = _window(tmp_path, None)

        caption = window.auto_anchor_value.text()
        assert "vetting" not in caption
        assert "ready to be measured in one pass" in caption
        window.close()
