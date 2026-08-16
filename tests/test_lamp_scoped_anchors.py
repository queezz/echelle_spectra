"""Packet F13 — the bench fits the lamp it was given.

Click-to-fit used to measure every accepted anchor against the whole curated
wavelength table, which is a ThAr table for the packaged defaults.  On the real
Ne-only 2025 data that produced a green but scientifically poor solution.  These
tests pin the lamp boundary: the nearest-row lookup never crosses lamps, a lamp
with no catalog says so, a mismatched line-help catalog warns, and the procedure
row names the catalog that actually anchored the fit.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.calibration_bench import (
    AlignmentState,
    BenchFrame,
    CalibrationBenchSession,
)
from echelle_spectra.calibration_bench_gui import CalibrationBenchWindow
from echelle_spectra.calibration_campaign import (
    CalibrationCampaignSession,
    ChecklistState,
    MeasurementRole,
    ReferenceState,
    catalog_mismatch_warning,
    lamp_reference_set,
)
from echelle_spectra.tools.calibration_alignment import (
    CalibrationTableLine,
    load_wavelength_table,
    select_candidate_lines,
)

_COLUMNS = 120

#: One order's worth of interleaved lamps, the situation the curated ThAr table
#: presents to a neon frame: every real Ne line has a thorium or argon row a few
#: pixels away, close enough for a click to land on the wrong one.
_INTERLEAVED = (
    CalibrationTableLine(0, 25, 35, 30.0, 585.249, "NeI", "ok"),
    CalibrationTableLine(0, 29, 39, 34.0, 586.104, "ThI", "ok"),
    CalibrationTableLine(1, 65, 75, 70.0, 600.062, "NeI", "ok"),
    CalibrationTableLine(1, 69, 79, 74.0, 601.113, "ArI", "ok"),
)

#: Where the neon frame really emits: on the two NeI rows and nowhere else.
_NE_PEAKS = (30.0, 70.0)


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


def _pattern(columns: int = _COLUMNS) -> np.ndarray:
    pattern = np.empty((columns, 2), dtype=float)
    pattern[:, 0] = 12
    pattern[:, 1] = 30
    return pattern


def _neon_frame(path: Path, shift_px: float = 0.0) -> BenchFrame:
    """A frame whose only lines sit on the NeI rows, moved by *shift_px*."""

    x = np.arange(_COLUMNS, dtype=float)
    spectra = tuple(
        10.0 + 500.0 * np.exp(-0.5 * ((x - (center + shift_px)) / 1.6) ** 2)
        for center in _NE_PEAKS
    )
    images = np.zeros((1, 44, _COLUMNS), dtype=float)
    pattern = _pattern()
    for order, spectrum in enumerate(spectra):
        for column, value in enumerate(spectrum):
            images[0, int(pattern[column, order]), column] = value
    return BenchFrame(path, images, images[0], spectra, {"ExposureTime": 0.1})


def _session(tmp_path: Path, lamp: str | None = None) -> CalibrationBenchSession:
    session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
    session.accept_frame(_neon_frame(tmp_path / "Ne-0.02s-x3-bright-lines.sif"))
    if lamp is not None:
        session.use_lamp_reference(lamp_reference_set(lamp, session.lines))
    return session


class TestTheLookupNeverCrossesLamps:
    def test_a_click_between_two_lamps_snaps_to_the_assigned_one(self, tmp_path):
        # Pixel 33 is one pixel from the ThI row and three from the NeI row, so
        # a lamp-blind lookup takes thorium and a Ne-scoped one cannot.
        scoped = _session(tmp_path, "Ne")
        blind = _session(tmp_path)

        chosen = scoped.fit_anchor_at(0, 33.0)
        crossed = blind.fit_anchor_at(0, 33.0)

        assert chosen.accepted and crossed.accepted
        assert chosen.anchor.line.species == "NeI"
        assert chosen.anchor.line.wavelength_nm == pytest.approx(585.249)
        assert crossed.anchor.line.species == "ThI"

    def test_the_other_lamp_selects_the_other_rows_from_the_same_table(self, tmp_path):
        session = _session(tmp_path, "ThAr")

        result = session.fit_anchor_at(0, 33.0)

        assert session.reference.state is ReferenceState.MATCHED
        assert {line.species for line in session.reference_lines()} == {"ThI", "ArI"}
        assert result.accepted
        assert result.anchor.line.species == "ThI"

    def test_only_the_assigned_lamps_rows_are_offered_for_clicking(self, tmp_path):
        session = _session(tmp_path, "Ne")

        assert [line.species for line in session.reference_lines()] == ["NeI", "NeI"]
        assert [line.species for line in session.lines_for_order(0)] == ["NeI"]
        assert [line.species for line in session.lines_for_order(1)] == ["NeI"]
        # The full table is still there; only what a click may snap to narrowed.
        assert len(session.lines) == 4

    def test_the_solved_transform_is_the_lamps_own(self, tmp_path):
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(
            _neon_frame(tmp_path / "Ne-0.02s.sif", shift_px=1.0)
        )
        session.use_lamp_reference(lamp_reference_set("Ne", session.lines))

        assert session.fit_anchor_at(0, 31.0).accepted
        assert session.fit_anchor_at(1, 71.0).accepted

        assert session.alignment_state is AlignmentState.ALIGNED
        assert session.transform.dx_px == pytest.approx(1.0, abs=0.05)
        assert session.rms_px < 0.05

    def test_changing_lamp_drops_anchors_the_old_lamp_referenced(self, tmp_path):
        session = _session(tmp_path)
        assert session.fit_anchor_at(0, 33.0).accepted
        assert session.fit_anchor_at(1, 73.0).accepted
        assert {anchor.line.species for anchor in session.anchor_rows()} == {"ThI", "ArI"}

        session.use_lamp_reference(lamp_reference_set("Ne", session.lines))

        assert not session.anchors
        assert session.alignment_state is AlignmentState.EMPTY

    def test_an_order_without_the_lamps_rows_says_which_catalog_is_missing(
        self, tmp_path
    ):
        rows = (_INTERLEAVED[0], _INTERLEAVED[3])
        session = CalibrationBenchSession(_pattern(), rows, minimum_snr=3.0)
        session.accept_frame(_neon_frame(tmp_path / "Ne.sif"))
        session.use_lamp_reference(lamp_reference_set("Ne", rows))

        result = session.fit_anchor_at(1, 70.0)

        assert not result.accepted
        assert result.reason == "no Ne rows for this order"

    def test_the_packaged_thar_table_yields_exactly_its_curated_neon_rows(self):
        resources = Path(__file__).parents[1] / "src" / "echelle_spectra" / "resources"
        table = (
            resources
            / "calibration_files"
            / "alignments"
            / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
        )
        rows = load_wavelength_table(table)

        reference = lamp_reference_set("Ne", rows)

        assert reference.state is ReferenceState.MATCHED
        assert reference.catalog_label == "Ne"
        assert {line.species for line in reference.lines} == {"NeI"}
        # The same rows `echelle-align` reaches 0.67 px with, before its extra
        # OK/width quality filter narrows them further.
        curated = select_candidate_lines(rows, species=("NeI",))
        assert set(curated).issubset(set(reference.lines))
        assert len(reference.lines) < len(rows)


class TestAbsenceExplainsAndMismatchWarns:
    def test_a_free_text_lamp_with_no_catalog_says_so(self, tmp_path):
        session = _session(tmp_path, "Kr")
        reference = session.reference

        assert reference.state is ReferenceState.NO_CATALOG
        assert not reference.is_referenceable
        assert reference.message == (
            "no line catalog for Kr — anchors cannot be auto-referenced"
        )
        assert session.reference_lines() == ()
        assert session.lines_for_order(0) == ()

    def test_a_click_under_a_catalogless_lamp_is_refused_with_that_reason(
        self, tmp_path
    ):
        session = _session(tmp_path, "Kr")

        result = session.fit_anchor_at(0, 30.0)

        assert not result.accepted
        assert "no line catalog for Kr" in result.reason
        assert not session.anchors

    def test_a_catalog_this_table_never_carries_explains_itself(self, tmp_path):
        reference = lamp_reference_set("Hg", _INTERLEAVED)

        assert reference.state is ReferenceState.NO_ROWS
        assert not reference.is_referenceable
        assert "packaged Hg catalog" in reference.message
        assert "HgI, HgII" in reference.message

    def test_no_lamp_at_all_keeps_every_row_and_asks_for_one(self):
        reference = lamp_reference_set("", _INTERLEAVED)

        assert reference.state is ReferenceState.UNSCOPED
        assert reference.is_referenceable
        assert reference.lines == _INTERLEAVED
        assert "assign a lamp role" in reference.message

    @pytest.mark.parametrize(
        ("active", "assigned"),
        [("ThAr", "Ne"), ("Ne", "ThAr"), ("Hg", "Kr"), ("H2", "Ne")],
    )
    def test_a_different_line_help_catalog_warns(self, active, assigned):
        warning = catalog_mismatch_warning(active, assigned)

        assert active in warning and assigned in warning
        assert "referenced against" in warning

    @pytest.mark.parametrize(
        ("active", "assigned"),
        [("Ne", "Ne"), ("ne", "Ne"), ("ThAr", "th-ar"), ("", "Ne"), ("Ne", "")],
    )
    def test_the_lamps_own_catalog_never_warns(self, active, assigned):
        assert catalog_mismatch_warning(active, assigned) == ""


def _campaign(tmp_path: Path) -> CalibrationCampaignSession:
    references = {}
    for name in ("pattern.txt", "wavelength.txt", "integral.txt"):
        path = tmp_path / name
        path.write_text(name, encoding="utf-8")
        references[name] = path
    return CalibrationCampaignSession(
        pattern_source=references["pattern.txt"],
        wavelength_source=references["wavelength.txt"],
        integral_source=references["integral.txt"],
    )


def _assign_lamp(
    campaign: CalibrationCampaignSession, tmp_path: Path, lamp: str
) -> Path:
    source = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    source.write_bytes(b"sif\n")
    campaign.classify_file(source, MeasurementRole.LAMP, lamp_family=lamp)
    return source


class TestTheChecklistNamesItsCatalog:
    def test_the_alignment_row_reports_anchors_catalog_and_rms(self, tmp_path):
        campaign = _campaign(tmp_path)
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(_neon_frame(_assign_lamp(campaign, tmp_path, "Ne"), 1.0))
        campaign.scope_alignment_to_lamp(session)
        assert session.fit_anchor_at(0, 31.0).accepted
        assert session.fit_anchor_at(1, 71.0).accepted

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        assert row.state is ChecklistState.DONE
        assert row.detail.startswith("2 anchors vs Ne catalog, RMS ")
        # F19's second rider deliberately supersedes the old "ends with px":
        # a solved fit now has to say whether it was ever held against the
        # science lines, and a neon frame honestly says it cannot be.
        assert "emits no Balmer or Fulcher light" in row.detail
        assert "validate against Fulcher" in row.unblocked_by

    def test_one_anchor_reads_as_one_anchor(self, tmp_path):
        campaign = _campaign(tmp_path)
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(_neon_frame(_assign_lamp(campaign, tmp_path, "Ne")))
        campaign.scope_alignment_to_lamp(session)
        assert session.fit_anchor_at(0, 30.0).accepted

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        assert row.detail == "1 anchor vs Ne catalog"
        assert row.state is ChecklistState.WAITING

    def test_a_catalogless_lamp_turns_the_row_into_an_attention(self, tmp_path):
        campaign = _campaign(tmp_path)
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(_neon_frame(_assign_lamp(campaign, tmp_path, "Kr")))
        campaign.scope_alignment_to_lamp(session)

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        assert row.state is ChecklistState.ATTENTION
        assert "no line catalog for Kr" in row.detail
        assert "ThAr, Ne, Hg, H2" in row.unblocked_by

    def test_an_unscoped_bench_admits_it_is_unscoped(self, tmp_path):
        campaign = _campaign(tmp_path)
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(_neon_frame(tmp_path / "unassigned.sif"))

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        assert "no lamp catalog is scoping the fit yet" in row.detail

    def test_the_open_frames_own_lamp_wins_over_the_other_assigned_one(self, tmp_path):
        campaign = _campaign(tmp_path)
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        thar = tmp_path / "ThAr-0.3s-x3.sif"
        thar.write_bytes(b"sif\n")
        _assign_lamp(campaign, tmp_path, "Ne")
        campaign.classify_file(thar, MeasurementRole.LAMP, lamp_family="ThAr")
        session.accept_frame(_neon_frame(thar))

        reference = campaign.scope_alignment_to_lamp(session)

        assert reference.lamp == "ThAr"
        assert {line.species for line in reference.lines} == {"ThI", "ArI"}


class TestTheFitPanelSaysWhichCatalogAnchorsIt:
    def _window(self, tmp_path: Path) -> CalibrationBenchWindow:
        session = CalibrationBenchSession(_pattern(), _INTERLEAVED, minimum_snr=3.0)
        session.accept_frame(_neon_frame(tmp_path / "bench-fixture.sif"))
        window = CalibrationBenchWindow(session, start_timer=False)
        window.campaign = _campaign(tmp_path)
        return window

    def test_an_unscoped_bench_asks_for_a_lamp(self, qt_app, tmp_path):
        window = self._window(tmp_path)
        window.refresh()
        qt_app.processEvents()

        assert "No lamp catalog is scoping the fit yet" in window.reference_value.text()
        window.close()

    def test_the_assigned_lamp_is_named_and_a_foreign_catalog_warns(
        self, qt_app, tmp_path
    ):
        window = self._window(tmp_path)
        window.campaign.scope_alignment_to_lamp(window.session)
        source = _assign_lamp(window.campaign, tmp_path, "Ne")
        window.session.accept_frame(_neon_frame(source))
        window.campaign.scope_alignment_to_lamp(window.session)
        window.line_family_combo.setCurrentText("Ne")
        window.refresh()
        qt_app.processEvents()
        assert "anchors reference Ne only" in window.reference_value.text()
        assert "WARNING" not in window.reference_value.text()

        window.line_family_combo.setCurrentText("ThAr")
        qt_app.processEvents()

        text = window.reference_value.text()
        assert "WARNING" in text
        assert "the ThAr line help on screen is not Ne's own catalog" in text
        window.close()

    def test_a_catalogless_lamp_explains_itself_in_the_panel(self, qt_app, tmp_path):
        window = self._window(tmp_path)
        source = _assign_lamp(window.campaign, tmp_path, "Kr")
        window.session.accept_frame(_neon_frame(source))
        window.campaign.scope_alignment_to_lamp(window.session)
        window.refresh()
        qt_app.processEvents()

        assert (
            "no line catalog for Kr — anchors cannot be auto-referenced"
            in window.reference_value.text()
        )
        window.close()
