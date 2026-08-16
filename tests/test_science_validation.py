"""Packet F19, second rider — good alignment is not good residuals.

The owner: "good alignment is not good residuals — it is agreement with the
science lines."  Anchor RMS is self-consistency, measured in pixels among
the anchors themselves; the BH paper's standard was agreement with
Fulcher-alpha, measured in nanometres against physics.  These tests pin that
the bench closes its alignment step on the second question: it answers it
where a hydrogen frame allows, it states plainly why it cannot where a neon
frame does not, and neither the checklist nor a saved snapshot is ever
allowed to look validated on the strength of RMS alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from echelle_spectra.calibration_bench import (
    BenchFrame,
    CalibrationBenchSession,
    ScienceValidationState,
)
from echelle_spectra.calibration_campaign import lamp_reference_set
from echelle_spectra.tools.calibration_alignment import CalibrationTableLine
from echelle_spectra.tools.line_validation import LineValidationTarget

_COLUMNS = 400
_ROWS = 40

#: A pair of hydrogen rows per order, far enough apart to fix a straight
#: wavelength solution across the whole detector width.
_H_ALPHA_NM = 656.279
_H_BETA_NM = 486.135


def _table() -> tuple[CalibrationTableLine, ...]:
    """Two orders carrying both lamps' rows, as a real curated table does.

    The neon rows sit exactly on the same straight solution as the hydrogen
    ones, so which lamp is assigned changes which rows anchor the fit without
    changing the wavelength axis the fit produces.  That is what lets one
    fixture ask both halves of the question: a hydrogen frame that can be
    validated, and a neon frame that honestly cannot.
    """

    def row(order, center, wavelength, species="H2"):
        return CalibrationTableLine(
            order, center - 5, center + 5, center, wavelength, species, "ok"
        )

    return (
        # Order 0 runs 650 → 660 nm, so H-alpha at 656.279 sits inside it.
        row(0, 40.0, 651.0),
        row(0, 360.0, 659.0),
        row(0, 120.0, 653.0, "NeI"),
        row(0, 280.0, 657.0, "NeI"),
        # Order 1 runs 482 → 490 nm, which is where H-beta lives.
        row(1, 40.0, 482.5),
        row(1, 360.0, 489.5),
        row(1, 120.0, 484.25, "NeI"),
        row(1, 280.0, 487.75, "NeI"),
    )


def _pattern() -> np.ndarray:
    pattern = np.empty((_COLUMNS, 2), dtype=float)
    pattern[:, 0] = 10
    pattern[:, 1] = 26
    return pattern


def _wavelength_of(order: int, column: np.ndarray) -> np.ndarray:
    """The same straight solution the table's two rows imply, for the fixture."""

    rows = [line for line in _table() if line.order_idx == order]
    pixels = np.asarray([line.center_pixel for line in rows], dtype=float)
    waves = np.asarray([line.wavelength_nm for line in rows], dtype=float)
    # Every row of an order is collinear by construction, so one straight fit
    # describes the whole order however many rows it holds.
    return np.poly1d(np.polyfit(pixels, waves, 1))(column)


def _hydrogen_frame(path: Path) -> BenchFrame:
    """A frame emitting the table's own rows plus real Balmer light."""

    column = np.arange(_COLUMNS, dtype=float)
    spectra = []
    for order in range(2):
        axis = _wavelength_of(order, column)
        trace = np.full(_COLUMNS, 20.0)
        # The curated rows themselves, so anchors can be placed.
        for line in _table():
            if line.order_idx == order:
                trace = trace + 800.0 * np.exp(
                    -0.5 * ((column - line.center_pixel) / 2.0) ** 2
                )
        # The science line this order carries, placed by wavelength.
        science = _H_ALPHA_NM if order == 0 else _H_BETA_NM
        if axis.min() <= science <= axis.max():
            centre = float(np.interp(science, axis, column))
            trace = trace + 1500.0 * np.exp(-0.5 * ((column - centre) / 2.2) ** 2)
        spectra.append(trace)

    images = np.zeros((1, _ROWS, _COLUMNS), dtype=float)
    pattern = _pattern()
    for order, spectrum in enumerate(spectra):
        for index, value in enumerate(spectrum):
            images[0, int(pattern[index, order]), index] = value
    return BenchFrame(path, images, images[0], tuple(spectra), {"ExposureTime": 0.5})


def _session(tmp_path: Path, lamp: str) -> CalibrationBenchSession:
    session = CalibrationBenchSession(_pattern(), _table(), minimum_snr=3.0)
    session.accept_frame(_hydrogen_frame(tmp_path / f"{lamp}.sif"))
    session.use_lamp_reference(lamp_reference_set(lamp, session.lines))
    return session


def _solved(tmp_path: Path, lamp: str = "H2") -> CalibrationBenchSession:
    session = _session(tmp_path, lamp)
    assert session.auto_anchor().accepted
    assert session.transform is not None
    return session


class TestTheBenchAnswersTheQuestionItCan:
    def test_a_hydrogen_frame_is_measured_against_its_science_lines(self, tmp_path):
        session = _solved(tmp_path)

        validation = session.validate_science_lines()

        assert validation.state is ScienceValidationState.MEASURED
        assert validation.measured
        assert validation.line_count >= 1
        # The fixture places the science lines exactly where the solution
        # predicts, so agreement is limited by the fit, not by the placement.
        assert validation.rms_residual_nm == pytest.approx(0.0, abs=0.02)
        assert "nm RMS" in validation.message

    def test_a_real_disagreement_is_reported_rather_than_smoothed(self, tmp_path):
        session = _solved(tmp_path)

        # Ask for a target that is genuinely NOT where the solution says: a
        # line half a nanometre off must show up as half a nanometre off.
        moved = [LineValidationTarget("H-alpha (moved)", _H_ALPHA_NM + 0.25)]
        validation = session.validate_science_lines(moved)

        assert validation.measured
        assert abs(validation.rms_residual_nm) == pytest.approx(0.25, abs=0.05)

    def test_the_wavelength_axis_is_built_per_order_from_the_solved_rows(
        self, tmp_path
    ):
        session = _solved(tmp_path)

        axes = session.order_wavelengths()

        assert set(axes) == {0, 1}
        for order, axis in axes.items():
            assert axis.size == _COLUMNS
            expected = _wavelength_of(order, np.arange(_COLUMNS, dtype=float))
            # The frame carries no shift, so the solved axis is the table's.
            assert np.allclose(axis, expected, atol=0.05)


class TestTheBenchAdmitsTheQuestionItCannot:
    def test_a_neon_frame_says_why_it_cannot_validate(self, tmp_path):
        session = _solved(tmp_path, lamp="Ne")

        validation = session.validate_science_lines()

        assert validation.state is ScienceValidationState.NO_FRAME
        assert not validation.measured
        assert "emits no Balmer or Fulcher light" in validation.message
        assert "validate against Fulcher when the first plasma data exists" in (
            validation.message
        )

    def test_an_unsolved_bench_has_nothing_to_validate_yet(self, tmp_path):
        session = _session(tmp_path, "H2")

        validation = session.validate_science_lines()

        assert validation.state is ScienceValidationState.NO_ALIGNMENT
        assert "nothing to validate" in validation.message

    def test_a_hydrogen_frame_with_no_fittable_line_says_that_instead(self, tmp_path):
        session = _solved(tmp_path)

        # Targets nowhere near this detector's range: nothing can be fitted,
        # which is a different answer from "no frame" and says so.
        validation = session.validate_science_lines(
            [LineValidationTarget("invented", 1200.0)]
        )

        assert validation.state is ScienceValidationState.NO_LINES
        assert "could be fitted" in validation.message


class TestNeitherSurfaceLooksValidatedOnRmsAlone:
    def test_the_checklist_row_carries_the_validation_it_has(self, tmp_path):
        from echelle_spectra.calibration_campaign import (
            CalibrationCampaignSession,
            ChecklistState,
            MeasurementRole,
        )

        references = {}
        for name in ("pattern.txt", "wavelength.txt", "integral.txt"):
            path = tmp_path / name
            path.write_text(name, encoding="utf-8")
            references[name] = path
        campaign = CalibrationCampaignSession(
            pattern_source=references["pattern.txt"],
            wavelength_source=references["wavelength.txt"],
            integral_source=references["integral.txt"],
        )
        source = tmp_path / "H2.sif"
        source.write_bytes(b"sif\n")
        campaign.classify_file(source, MeasurementRole.LAMP, lamp_family="H2")
        session = _solved(tmp_path)
        campaign.scope_alignment_to_lamp(session)

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        assert row.state is ChecklistState.DONE
        assert "nm RMS" in row.detail
        # Validated here, so nothing is left to carry to first plasma.
        assert row.unblocked_by == ""

    def test_an_unvalidated_solve_states_what_is_still_owed(self, tmp_path):
        from echelle_spectra.calibration_campaign import (
            CalibrationCampaignSession,
            ChecklistState,
            MeasurementRole,
        )

        references = {}
        for name in ("pattern.txt", "wavelength.txt", "integral.txt"):
            path = tmp_path / name
            path.write_text(name, encoding="utf-8")
            references[name] = path
        campaign = CalibrationCampaignSession(
            pattern_source=references["pattern.txt"],
            wavelength_source=references["wavelength.txt"],
            integral_source=references["integral.txt"],
        )
        source = tmp_path / "Ne.sif"
        source.write_bytes(b"sif\n")
        campaign.classify_file(source, MeasurementRole.LAMP, lamp_family="Ne")
        session = _solved(tmp_path, lamp="Ne")
        campaign.scope_alignment_to_lamp(session)

        row = {item.key: item for item in campaign.checklist(session)}["alignment"]

        # Still DONE — the fit really is solved and reviewed — but it does not
        # pretend the science question was asked, and it says who must ask it.
        assert row.state is ChecklistState.DONE
        assert "emits no Balmer or Fulcher light" in row.detail
        assert "validate against Fulcher there" in row.unblocked_by
