"""Packet F19 — the bench anchors the lines we already trust.

``echelle-align`` has always selected its own anchors: the curated table's
``OK``-marked rows, ranked, saturation-guarded, centroid-fitted.  That path
reached 0.67 px over 37 lines with nobody clicking anything.  The bench asked
the operator to click every line by hand instead.  These tests pin the pass
that closes the gap: it trusts the same rows the CLI trusts, it reaches the
same numbers over the same data, it declines out loud, and it composes with
the manual clicking rather than replacing it.
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
from echelle_spectra.calibration_campaign import lamp_reference_set
from echelle_spectra.tools.calibration_alignment import (
    CalibrationTableLine,
    detector_points_from_lines,
    fit_rigid_transform,
    measure_detector_window_saturation,
    measure_line_centroids,
    select_candidate_lines,
)

_COLUMNS = 200
_ORDERS = 4
_ROWS = 60
_SATURATION = 60000.0

#: The detector shift the fixture frame carries, which both paths must find.
_SHIFT_PX = 1.4


def _rows() -> tuple[CalibrationTableLine, ...]:
    """A curated table in miniature: trusted rows, and the traps beside them.

    Four vetted neon rows spread over four orders are what an auto pass should
    take.  Everything else here is something it must leave: an unvetted neon
    row, another lamp's vetted row, a vetted row whose line this frame never
    emitted, a vetted row sitting on saturated pixels, and a vetted row whose
    hand-picked interval is too wide to be one line.
    """

    def line(order, center, species="NeI", comment="ok", half_width=5.0):
        return CalibrationTableLine(
            order, center - half_width, center + half_width, center, 600.0 + center,
            species, comment,
        )

    return (
        line(0, 40.0),
        line(1, 80.0),
        line(2, 120.0),
        line(3, 160.0),
        line(0, 70.0, comment="?"),          # real light, never vetted
        line(1, 100.0, species="ThI"),       # vetted, but another lamp's
        line(2, 30.0),                       # vetted, but dark in this frame
        line(3, 60.0),                       # vetted, but saturated here
        line(0, 150.0, half_width=25.0),     # vetted, but 50 px is not one line
    )


#: Where the fixture frame actually emits, before the shift is applied.
_EMITTING = ((0, 40.0), (1, 80.0), (2, 120.0), (3, 160.0), (0, 70.0), (1, 100.0),
             (3, 60.0), (0, 150.0))

#: The rows a correct pass anchors: vetted, neon, present, clean, one line wide.
_EXPECTED_ANCHORS = ((0, 40.0), (1, 80.0), (2, 120.0), (3, 160.0))


def _pattern() -> np.ndarray:
    pattern = np.empty((_COLUMNS, _ORDERS), dtype=float)
    for order in range(_ORDERS):
        pattern[:, order] = 8 + 12 * order
    return pattern


def _frame(path: Path, shift_px: float = _SHIFT_PX) -> BenchFrame:
    """A neon acquisition carrying a rigid shift and one saturated line."""

    x = np.arange(_COLUMNS, dtype=float)
    spectra = []
    for order in range(_ORDERS):
        trace = np.full(_COLUMNS, 10.0)
        for line_order, center in _EMITTING:
            if line_order != order:
                continue
            height = 40000.0 if (order, center) == (3, 60.0) else 900.0
            trace = trace + height * np.exp(
                -0.5 * ((x - (center + shift_px)) / 1.7) ** 2
            )
        spectra.append(trace)

    pattern = _pattern()
    images = np.zeros((1, _ROWS, _COLUMNS), dtype=float)
    for order, spectrum in enumerate(spectra):
        for column, value in enumerate(spectrum):
            row = int(pattern[column, order])
            images[0, row, column] = value
    # The one line the raw detector clips, which the extracted trace alone
    # would never reveal: an integrated trace can exceed full scale honestly.
    saturated_column = int(round(60.0 + shift_px))
    saturated_row = int(pattern[saturated_column, 3])
    images[0, saturated_row - 1 : saturated_row + 2,
           saturated_column - 2 : saturated_column + 3] = _SATURATION + 100.0

    return BenchFrame(path, images, images[0], tuple(spectra), {"ExposureTime": 0.02})


def _session(tmp_path: Path, lamp: str | None = "Ne") -> CalibrationBenchSession:
    session = CalibrationBenchSession(
        _pattern(), _rows(), saturation_level=_SATURATION, minimum_snr=5.0
    )
    session.accept_frame(_frame(tmp_path / "Ne-0.02s-x3-bright-lines.sif"))
    if lamp is not None:
        session.use_lamp_reference(lamp_reference_set(lamp, session.lines))
    return session


def _cli_alignment(session: CalibrationBenchSession):
    """Run the CLI's own primitives over the bench's data.

    ``run_calibration_alignment`` is this sequence with SIF loading in front of
    it.  Comparing against the sequence rather than the wrapper keeps the pin
    on the arithmetic both paths must share, with no raw data in the tests.
    """

    spectra = session.active_order_spectra()
    candidates = select_candidate_lines(
        [line for line in session.lines if line.species == "NeI"],
        species=("NeI",),
    )
    saturation = measure_detector_window_saturation(
        session.active_images(),
        session.pattern,
        candidates,
        x_radius_px=session.fit_window_radius_px,
        saturation_level=session.saturation_level,
    )
    fits = measure_line_centroids(
        spectra,
        candidates,
        window_radius_px=session.fit_window_radius_px,
        saturation_level=None,
        min_snr=session.minimum_snr,
    )
    good = [
        fit
        for fit, qc in zip(fits, saturation)
        if fit.success and not qc.is_saturated
    ]
    lines = [fit.line for fit in good]
    expected = detector_points_from_lines(lines, session.pattern)
    measured = detector_points_from_lines(
        lines, session.pattern, [fit.center_pixel for fit in good]
    )
    transform, rms_px = fit_rigid_transform(expected, measured)
    return len(good), transform, rms_px


class TestTheAutoPassReachesTheClisNumbers:
    def test_one_press_anchors_every_trusted_line_in_every_order(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        assert result.ran
        anchored = {
            (anchor.line.order_idx, anchor.line.center_pixel)
            for anchor in result.accepted
        }
        assert anchored == set(_EXPECTED_ANCHORS)
        assert session.alignment_state is AlignmentState.ALIGNED

    def test_the_anchor_count_and_rms_match_the_cli_path(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()
        cli_count, cli_transform, cli_rms = _cli_alignment(session)

        assert len(result.accepted) == cli_count
        assert session.rms_px == pytest.approx(cli_rms, abs=1e-9)
        assert session.transform.dx_px == pytest.approx(cli_transform.dx_px, abs=1e-9)
        assert session.transform.dy_px == pytest.approx(cli_transform.dy_px, abs=1e-9)

    def test_the_pass_finds_the_shift_the_frame_carries(self, tmp_path):
        session = _session(tmp_path)

        session.auto_anchor()

        assert session.transform.dx_px == pytest.approx(_SHIFT_PX, abs=0.05)
        assert session.rms_px < 0.1

    def test_the_solve_runs_once_for_the_whole_pass(self, tmp_path, monkeypatch):
        session = _session(tmp_path)
        calls = []
        original = CalibrationBenchSession._recompute_alignment
        monkeypatch.setattr(
            CalibrationBenchSession,
            "_recompute_alignment",
            lambda self: (calls.append(1), original(self))[1],
        )

        session.auto_anchor()

        assert len(calls) == 1


class TestOnlyVettedRowsAreAnchored:
    def test_a_row_without_an_ok_mark_is_never_taken(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        # Pixel 70 of order 0 emits as brightly as any anchor and would fit.
        assert (0, 70.0) not in {
            (anchor.line.order_idx, anchor.line.center_pixel)
            for anchor in result.accepted
        }
        assert all(line.comment.lower() == "ok" for line in
                   [anchor.line for anchor in result.accepted])

    def test_the_pass_stays_inside_the_assigned_lamps_catalog(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        assert {anchor.line.species for anchor in result.accepted} == {"NeI"}

    def test_the_other_lamp_anchors_its_own_row_from_the_same_table(self, tmp_path):
        session = _session(tmp_path, lamp="ThAr")

        result = session.auto_anchor()

        assert [anchor.line.species for anchor in result.accepted] == ["ThI"]

    def test_an_interval_too_wide_to_be_one_line_is_skipped(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        measured = {line.center_pixel for line in
                    [anchor.line for anchor in result.accepted]}
        measured |= {rejection.line.center_pixel for rejection in result.rejected}
        assert 150.0 not in measured


class TestEveryDeclineSaysWhy:
    def test_a_saturated_line_is_declined_on_its_raw_pixels(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        declined = {
            rejection.line.center_pixel: rejection.reason
            for rejection in result.rejected
        }
        assert "saturated on the raw detector" in declined[60.0]
        assert "full-scale pixel" in declined[60.0]

    def test_a_line_this_frame_never_emitted_is_declined_with_the_fit_reason(
        self, tmp_path
    ):
        session = _session(tmp_path)

        result = session.auto_anchor()

        declined = {
            rejection.line.center_pixel: rejection.reason
            for rejection in result.rejected
        }
        assert "snr" in declined[30.0].lower()

    def test_the_pass_reports_everything_it_measured(self, tmp_path):
        session = _session(tmp_path)

        result = session.auto_anchor()

        assert result.considered == len(result.accepted) + len(result.rejected)
        assert result.considered == 6  # four vetted neon lines, plus the two traps


class TestThePassIsHonestWhenItCannotRun:
    def test_without_a_frame_it_says_so_and_anchors_nothing(self, tmp_path):
        session = CalibrationBenchSession(_pattern(), _rows())

        result = session.auto_anchor()

        assert not result.ran
        assert result.reason == "no SIF frame loaded"
        assert not session.anchors

    def test_a_lamp_with_no_catalog_repeats_its_own_reason(self, tmp_path):
        session = _session(tmp_path, lamp="Kr")

        result = session.auto_anchor()

        assert not result.ran
        assert "no line catalog for Kr" in result.reason
        assert not session.anchors

    def test_a_table_carrying_no_vetted_rows_offers_the_manual_path(self, tmp_path):
        unvetted = tuple(
            CalibrationTableLine(0, 35.0, 45.0, 40.0, 640.0, "NeI", "?")
            for _ in range(1)
        )
        session = CalibrationBenchSession(_pattern(), unvetted)
        session.accept_frame(_frame(tmp_path / "Ne.sif"))
        session.use_lamp_reference(lamp_reference_set("Ne", session.lines))

        result = session.auto_anchor()

        assert not result.ran
        assert "carry the OK marks the auto-anchor trusts" in result.reason
        assert "click the lines you trust instead" in result.reason


class TestAutoAndManualCompose:
    def test_a_hand_placed_anchor_survives_an_auto_pass(self, tmp_path):
        session = _session(tmp_path)
        assert session.fit_anchor_at(0, 71.4).accepted
        manual_key = next(iter(session.anchors))

        session.auto_anchor()

        assert manual_key in session.anchors
        assert len(session.anchors) == len(_EXPECTED_ANCHORS) + 1

    def test_an_auto_anchor_comes_off_the_way_it_went_on(self, tmp_path):
        session = _session(tmp_path)
        result = session.auto_anchor()

        removed = session.remove_anchor(result.accepted[0].key)

        assert removed
        assert len(session.anchors) == len(_EXPECTED_ANCHORS) - 1
        assert session.alignment_state is AlignmentState.ALIGNED

    def test_a_second_pass_re_measures_rather_than_duplicates(self, tmp_path):
        session = _session(tmp_path)
        session.auto_anchor()

        second = session.auto_anchor()

        assert len(session.anchors) == len(_EXPECTED_ANCHORS)
        assert len(second.accepted) == len(_EXPECTED_ANCHORS)

    def test_changing_lamp_still_drops_what_the_old_lamp_referenced(self, tmp_path):
        session = _session(tmp_path, lamp="ThAr")
        assert session.auto_anchor().accepted

        session.use_lamp_reference(lamp_reference_set("Ne", session.lines))

        assert not session.anchors


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


class TestTheBenchOffersThePassAndSaysWhatItDid:
    def _window(self, tmp_path: Path, lamp: str | None = "Ne"):
        window = CalibrationBenchWindow(_session(tmp_path, lamp), start_timer=False)
        window.refresh()
        return window

    def test_the_caption_says_what_a_press_would_measure(self, qt_app, tmp_path):
        window = self._window(tmp_path)

        caption = window.auto_anchor_value.text()

        # Six is what the pass will MEASURE, which is all it can honestly
        # promise: whether a line fits is not knowable until it is fitted.
        assert "6 Ne line(s)" in caption
        assert "4 order(s)" in caption
        window.close()

    def test_the_pass_is_offered_before_a_single_anchor_exists(self, qt_app, tmp_path):
        window = self._window(tmp_path)

        # Remove and Clear act on a row, so an empty table disables them; the
        # pass that FILLS the table must not be gated the same way.
        assert not window.session.anchors
        assert window.auto_anchor_button.isEnabled()
        assert not window.remove_button.isEnabled()
        assert not window.clear_button.isEnabled()
        window.close()

    def test_pressing_it_fills_the_table_and_reports_the_alignment(
        self, qt_app, tmp_path
    ):
        window = self._window(tmp_path)

        window.auto_anchor_button.click()
        qt_app.processEvents()

        assert window.anchor_table.rowCount() == len(_EXPECTED_ANCHORS)
        message = window.message_value.text()
        assert "anchored 4" in message
        assert "RMS" in message
        assert "right-click any line on the spectrum to drop it" in message
        window.close()

    def test_the_message_groups_the_declines_by_reason(self, qt_app, tmp_path):
        window = self._window(tmp_path)

        window.auto_anchor_button.click()
        qt_app.processEvents()

        message = window.message_value.text()
        assert "Declined 2" in message
        assert "saturated on the raw detector" in message
        window.close()

    def test_a_hand_placed_anchor_is_reported_as_kept(self, qt_app, tmp_path):
        window = self._window(tmp_path)
        assert window.session.fit_anchor_at(0, 71.4).accepted

        window.auto_anchor_button.click()
        qt_app.processEvents()

        assert "Your 1 hand-placed anchor(s) were kept." in window.message_value.text()
        window.close()

    def test_a_lamp_with_no_catalog_disables_the_pass_and_says_why(
        self, qt_app, tmp_path
    ):
        window = self._window(tmp_path, lamp="Kr")

        assert not window.auto_anchor_button.isEnabled()
        assert "no line catalog for Kr" in window.auto_anchor_value.text()
        window.close()

    def test_without_a_frame_the_pass_is_not_offered(self, qt_app, tmp_path):
        session = CalibrationBenchSession(_pattern(), _rows())
        window = CalibrationBenchWindow(session, start_timer=False)
        window.refresh()

        assert not window.auto_anchor_button.isEnabled()
        assert window.auto_anchor_value.text() == "no acquisition open"
        window.close()
