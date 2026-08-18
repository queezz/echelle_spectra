"""The 2026-08-18 rehearsal on real light, encoded as the judge's ground truth.

Three cases came off the first real-data rehearsal, and the verdict logic of
the day failed two of them.  Each one is written here as the residual set the
audit actually measured, rather than as a synthetic cube, because what is under
test is the arithmetic that turns residuals into a word -- and because the
numbers themselves are the evidence.  Every figure the rehearsal recorded is
marked as recorded; the two places where a number had to be constructed to
reproduce a recorded *statistic* say so at the point of construction.

  1. Shot 193778, 2025 cube through the matched 2025 calibration.  A flawless
     era match.  It read ``misaligned-beyond-repair``.
  2. The same, for 2019 through the matched 2019 calibration.  Also
     ``misaligned-beyond-repair``.
  3. The 2019 cube recalibrated onto the *wrong* era, whose hydrogen Balmer
     lines were re-identified as deuterium -- and whose flipped fit then scored
     better than the correct one.

Cases 1 and 2 were condemned by weak lines: strong Balmer lines agreeing to
about a pixel, drowned by faint Fulcher centroids scattering over tens of
pixels inside a fit window wide enough to reach their neighbours.  Case 3 is
the science hazard: a silent isotope flip that made a wrong-era calibration
look like the better answer.
"""

from __future__ import annotations

from datetime import date
from typing import Any

import pytest

from echelle_spectra.drift import (
    ERA_MISASSIGNED_VERDICT,
    GEOMETRIC_DAMAGE_PX,
    SHIFT_CONSISTENCY_PX,
    VERDICT_SNR_FLOOR,
    deuterium_prior,
    isotope_ambiguity,
    verdict_from_evidence,
)

# The audited orders and the dispersion each carries, from the shipped
# instrument geometry the detector-space judge already works in.
DISPERSION_NM_PER_PX = {6: 0.0108, 9: 0.0102, 12: 0.0080, 16: 0.0071, 18: 0.0068}

# H-gamma to D-gamma is 0.1187 nm, which is 16.5 px at order 16's dispersion --
# and the same 16.5 px in every audited order, because the gap and the
# dispersion grow together.  That is the degeneracy this file's third case is
# about.
ISOTOPE_SEPARATION_PX = 16.5


def _line(
    *,
    order: int,
    expected_nm: float,
    residual_px: float,
    snr: float,
    label: str,
    shot: str,
    family: str = "balmer",
    isotope: str = "H",
    flipped: bool = False,
) -> dict[str, Any]:
    """One measured line record, shaped exactly as ``audit_cubes`` writes it."""

    dispersion = DISPERSION_NM_PER_PX[order]
    record: dict[str, Any] = {
        "status": "measured",
        "family": family,
        "line": label,
        "cube": f"{shot}.nc",
        "shot_number": shot,
        "date": "",
        "blended": False,
        "isotope": isotope,
        "expected_nm": expected_nm,
        "echelle_order": order,
        "detector_pixel": 1152.0,
        "dispersion_nm_per_px": dispersion,
        "residual_nm": residual_px * dispersion,
        "pixel_residual_px": residual_px,
        "snr": snr,
    }
    if flipped:
        # Both references the window was judged against, as the audit keeps
        # them: deuterium sits ISOTOPE_SEPARATION_PX blueward of hydrogen, so
        # the hydrogen reading of this same centroid is that much lower.
        record["isotope_candidates"] = [
            {
                "isotope": "D",
                "line": label,
                "expected_nm": expected_nm,
                "pixel_residual_px": residual_px,
            },
            {
                "isotope": "H",
                "line": label.replace("D-", "H-"),
                "expected_nm": expected_nm + ISOTOPE_SEPARATION_PX * dispersion,
                "pixel_residual_px": residual_px - ISOTOPE_SEPARATION_PX,
            },
        ]
    return record


def _as_if_all_lines_were_strong(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return the same residuals with every line above the verdict's SNR floor.

    This reconstructs the pre-rehearsal rule exactly -- every resolved line
    entering the arithmetic on equal terms -- so each case below can show what
    it used to return without keeping a second copy of the verdict code.
    """

    return [{**item, "snr": 10.0 * VERDICT_SNR_FLOOR} for item in lines]


# ---------------------------------------------------------------------------
# Case 1 — shot 193778, 2025 cube through the matched 2025 calibration
# ---------------------------------------------------------------------------

# Recorded that night: four Balmer residuals in pixels, and a peak SNR of 440.
SHOT_193778_BALMER_PX = ((6, 656.2790, 1.21), (12, 486.1350, 0.66), (16, 434.0470, 0.65), (18, 410.1734, 0.38))  # fmt: skip

#: Recorded only as "SNR up to 440".  Every Balmer line that night was far
#: above any plausible metrology floor, so the fixture gives them the one
#: figure that was written down rather than inventing four.  What is under test
#: is which lines the verdict rests on, not exactly how bright each one was.
SHOT_193778_BALMER_SNR = 440.0

#: Recorded as the extremes of the weak Fulcher scatter: -10.1 px to +16.9 px,
#: centroided in +/-0.4 nm windows wide enough to catch a neighbour.  The two
#: endpoints are what the rehearsal wrote down, so they are what is encoded.
SHOT_193778_FULCHER_PX = (-10.1, 16.9)
SHOT_193778_FULCHER_SNR = 6.0


def _shot_193778() -> list[dict[str, Any]]:
    lines = [
        _line(
            order=order,
            expected_nm=nm,
            residual_px=residual,
            snr=SHOT_193778_BALMER_SNR,
            label=f"H-{order}",
            shot="193778",
        )
        for order, nm, residual in SHOT_193778_BALMER_PX
    ]
    lines += [
        _line(
            order=9,
            expected_nm=620.0 + index,
            residual_px=residual,
            snr=SHOT_193778_FULCHER_SNR,
            label=f"Q-{index}",
            shot="193778",
            family="fulcher",
        )
        for index, residual in enumerate(SHOT_193778_FULCHER_PX)
    ]
    return lines


def test_a_flawless_2025_era_match_is_no_longer_condemned_by_its_weak_lines() -> None:
    lines = _shot_193778()

    verdict, summary = verdict_from_evidence(lines)

    # The acceptance the rehearsal set: a matched era must not read as damage.
    assert verdict == "shifted"
    # The Balmer median, recorded that night as 0.645 px; these four residuals
    # give 0.655, the same number to the precision they were written down in.
    assert summary["median_shift_px"] == pytest.approx(0.65, abs=0.02)
    assert summary["verdict_lines"] == len(SHOT_193778_BALMER_PX)
    assert summary["maximum_pixel_deviation_px"] < SHIFT_CONSISTENCY_PX
    assert summary["strong_lines_inconsistent"] is False


def test_the_old_all_lines_rule_condemned_that_same_flawless_match() -> None:
    """Non-vacuity: without the floor these residuals still read beyond repair.

    If this ever returned ``shifted`` the test above would be proving nothing,
    and the rehearsal's first failure would have been a story about a threshold
    rather than about which lines are allowed to decide.
    """

    verdict, summary = verdict_from_evidence(_as_if_all_lines_were_strong(_shot_193778()))

    assert verdict == "misaligned-beyond-repair"
    assert summary["maximum_pixel_deviation_px"] > SHIFT_CONSISTENCY_PX
    # And the advice that verdict carries -- reprocess the raw SIF -- could not
    # have repaired any of it, because none of it was detector motion.


def test_the_weak_lines_are_reported_beside_the_verdict_in_their_own_words() -> None:
    _verdict, summary = verdict_from_evidence(_shot_193778())

    scatter = summary["weak_line_scatter"]
    assert scatter["lines"] == len(SHOT_193778_FULCHER_PX)
    assert scatter["snr_floor"] == VERDICT_SNR_FLOOR
    assert scatter["pixel_residual_range_px"] == [
        min(SHOT_193778_FULCHER_PX),
        max(SHOT_193778_FULCHER_PX),
    ]
    # It says what it is: a figure for how crowded and how clean the spectrum
    # is, and explicitly not the thing that decided the verdict.
    assert "did not decide the verdict" in scatter["note"]
    assert "crowded" in scatter["note"]
    assert summary["verdict_basis"].startswith("decided on the 4 resolved line(s)")


# ---------------------------------------------------------------------------
# Case 2 — the 2019 cube through the matched 2019 calibration
# ---------------------------------------------------------------------------

# Recorded that night: H-alpha, H-beta and H-gamma at these residuals and these
# signal-to-noise ratios, on a calibration that was the correct era for the
# shot.
MATCHED_2019_BALMER = (
    (6, 656.2790, 1.82, 83.0),
    (12, 486.1350, 3.20, 105.0),
    (16, 434.0470, 2.22, 68.0),
)

#: Recorded for the whole line set, weak lines included: median 0.397 px and a
#: spread of 13.5 px.  The individual weak residuals were not written down, so
#: these four are constructed to reproduce exactly those two recorded
#: statistics against the three recorded Balmer residuals -- sorted, the seven
#: put 0.397 in the middle and put the farthest line 13.5 px from it.
MATCHED_2019_WEAK_PX = (-13.10, -4.2, -1.6, 0.397)
MATCHED_2019_WEAK_SNR = 7.0


def _matched_2019() -> list[dict[str, Any]]:
    lines = [
        _line(order=order, expected_nm=nm, residual_px=residual, snr=snr, label=f"H-{order}", shot="2019")
        for order, nm, residual, snr in MATCHED_2019_BALMER
    ]
    lines += [
        _line(
            order=9,
            expected_nm=620.0 + index,
            residual_px=residual,
            snr=MATCHED_2019_WEAK_SNR,
            label=f"Q-{index}",
            shot="2019",
            family="fulcher",
        )
        for index, residual in enumerate(MATCHED_2019_WEAK_PX)
    ]
    return lines


def test_the_constructed_weak_set_reproduces_the_recorded_whole_set_figures() -> None:
    """The one place a number was constructed, checked against what it stands for."""

    import numpy as np

    everything = np.array(
        [residual for _order, _nm, residual, _snr in MATCHED_2019_BALMER]
        + list(MATCHED_2019_WEAK_PX)
    )
    median = float(np.median(everything))

    assert median == pytest.approx(0.397, abs=0.001)
    assert float(np.max(np.abs(everything - median))) == pytest.approx(13.5, abs=0.05)


def test_the_matched_2019_era_reads_as_a_small_stated_correction() -> None:
    lines = _matched_2019()

    verdict, summary = verdict_from_evidence(lines)

    assert verdict == "shifted"
    # The strong lines agree on ~2.2 px, and they agree to 0.98 px -- inside the
    # one pixel a single translation allows, but only just.  That narrowness is
    # the measurement's, not the rule's: it is what three real Balmer lines at
    # SNR 68..105 looked like on a correctly matched era.
    assert summary["median_shift_px"] == pytest.approx(2.22, abs=0.01)
    assert summary["maximum_pixel_deviation_px"] == pytest.approx(0.98, abs=0.01)
    assert summary["verdict_lines"] == len(MATCHED_2019_BALMER)
    assert summary["verdict_snr_range"] == [68.0, 105.0]
    # The misleading number -- the whole set's 0.397 px median -- is nowhere
    # near the verdict, and the scatter that produced it is reported as quality.
    assert summary["weak_line_scatter"]["lines"] == len(MATCHED_2019_WEAK_PX)


def test_the_old_all_lines_rule_condemned_the_matched_2019_era_too() -> None:
    verdict, summary = verdict_from_evidence(_as_if_all_lines_were_strong(_matched_2019()))

    assert verdict == "misaligned-beyond-repair"
    assert summary["median_shift_px"] == pytest.approx(0.397, abs=0.001)
    assert summary["maximum_pixel_deviation_px"] == pytest.approx(13.5, abs=0.05)


def test_no_real_drive_could_pass_the_epoch_gate_under_the_old_rule() -> None:
    """Both matched eras, judged the old way, failed. That was the whole defect."""

    for lines in (_shot_193778(), _matched_2019()):
        old, _summary = verdict_from_evidence(_as_if_all_lines_were_strong(lines))
        new, _summary = verdict_from_evidence(lines)

        assert old == "misaligned-beyond-repair"
        assert new == "shifted"


# ---------------------------------------------------------------------------
# Case 3 — the silent isotope flip
# ---------------------------------------------------------------------------

# Recorded that night: the 2019 cube recalibrated onto the wrong era had its
# hydrogen Balmer lines re-identified as deuterium at these residuals and these
# signal-to-noise ratios.  Their median, 2.59 px, looked *better* than the
# correct baseline's 2.22 px -- which is exactly how a wrong-era calibration
# hides.
FLIPPED_2019_BALMER = (
    (16, 433.9283, 2.10, 46.0, "D-gamma"),
    (12, 486.0001, 2.60, 35.0, "D-beta"),
    (6, 656.1012, 2.71, 5.0, "D-alpha"),
)

#: The rehearsal's 2019 shot: inside the LHD deuterium phase, so the calendar
#: allows deuterium and cannot settle the question by itself.
FLIPPED_SHOT_DATE = date(2019, 11, 21)

#: The same residuals dated after the deuterium phase closed on 2022-12-02.
#: Nothing measured is changed; only the calendar's answer is, which is the
#: one thing that turns the finding from a warning into a refusal.
AFTER_THE_DEUTERIUM_PHASE = date(2025, 6, 1)


def _flipped_2019() -> list[dict[str, Any]]:
    return [
        _line(
            order=order,
            expected_nm=nm,
            residual_px=residual,
            snr=snr,
            label=label,
            shot="2019-wrong-era",
            isotope="D",
            flipped=True,
        )
        for order, nm, residual, snr, label in FLIPPED_2019_BALMER
    ]


def test_the_flipped_fit_looks_better_than_the_correct_one() -> None:
    """Why the flip is a hazard and not an inconvenience, stated as numbers.

    "Better" is not a smaller median -- it is a tighter one.  Re-identified as
    deuterium, the three lines describe one clean rigid shift of 2.59 px that
    they agree on to 0.3 px, while the correct era's own reading that night
    scattered by 13.5 px and was condemned outright.  Judged on how well the
    residuals hang together, the wrong-era calibration wins, and an audit that
    ranked fits without asking which references they were fitted against would
    have handed the operator the wrong epoch with a confident word attached.
    """

    import numpy as np

    flipped = np.array([residual for _o, _nm, residual, _s, _l in FLIPPED_2019_BALMER])
    flipped_median = float(np.median(flipped))

    assert flipped_median == pytest.approx(2.59, abs=0.02)
    assert float(np.max(np.abs(flipped - flipped_median))) == pytest.approx(0.5, abs=0.2)

    # Both read through the pre-rehearsal rule: the wrong era earns a confident,
    # repairable shift and the right era is condemned.
    wrong_era, _summary = verdict_from_evidence(_as_if_all_lines_were_strong(_flipped_2019()))
    right_era, _summary = verdict_from_evidence(_as_if_all_lines_were_strong(_matched_2019()))

    assert wrong_era == "shifted"
    assert right_era == "misaligned-beyond-repair"


def test_the_flip_raises_the_degeneracy_finding_naming_both_readings() -> None:
    prior = deuterium_prior(FLIPPED_SHOT_DATE, "2019")
    assert prior["expectation"] == "deuterium possible"

    ambiguity = isotope_ambiguity(_flipped_2019(), prior)

    assert ambiguity is not None
    assert ambiguity["finding"] == (
        "isotope flip or era shift, degenerate; calendar says D-possible"
    )
    assert ambiguity["flipped_lines"] == len(FLIPPED_2019_BALMER)
    # Read as hydrogen, the same three centroids are one common ~-13.9 px shift
    # against a 16.5 px H/D separation: the two explanations differ by the 2.6 px
    # the deuterium fit is itself out by, well inside the degeneracy window.
    assert ambiguity["implied_common_shift_px"] == pytest.approx(-13.9, abs=0.05)
    assert ambiguity["isotope_separation_px"] == pytest.approx(ISOTOPE_SEPARATION_PX, abs=0.01)
    assert ambiguity["flipped_fit_residual_px"] == pytest.approx(2.6, abs=0.05)
    assert ambiguity["excludes_deuterium"] is False
    assert "confirm the epoch" in ambiguity["detail"]


def test_the_better_looking_flipped_fit_is_never_reported_as_a_verdict() -> None:
    """Inside the deuterium window the calendar cannot decide -- and neither does this.

    Only two of the three flipped lines reach the SNR floor (D-alpha came in at
    5), so the strong lines cannot carry a verdict at all, and the audit says
    insufficient-data rather than publishing the flattering 2.59 px median.
    """

    prior = deuterium_prior(FLIPPED_SHOT_DATE, "2019")
    lines = _flipped_2019()

    verdict, summary = verdict_from_evidence(
        lines, isotope_ambiguity=isotope_ambiguity(lines, prior)
    )

    assert verdict == "insufficient-data"
    assert "median_shift_px" not in summary
    quorum = summary["quorum"]
    assert quorum["satisfied"] is False
    assert quorum["below_snr_floor_lines"] == 1
    assert f"SNR {VERDICT_SNR_FLOOR:g} floor" in quorum["reason"]


def test_a_calendar_that_excludes_deuterium_names_the_misassigned_era() -> None:
    """The forcing clause: same measurement, a date the calendar can rule on.

    Nothing about the spectroscopy changes -- the lines are still assigned to
    deuterium and their residuals are untouched.  What changes is that the
    calendar says the plasma cannot have been deuterium, which makes the
    hydrogen reading of the same centroids the likelier one, and that reading
    is a calibration from the wrong epoch.
    """

    prior = deuterium_prior(AFTER_THE_DEUTERIUM_PHASE, "2019")
    assert prior["expectation"] == "hydrogen expected"

    # Every line strong, so quorum is not what withholds the verdict here.
    lines = _as_if_all_lines_were_strong(_flipped_2019())
    ambiguity = isotope_ambiguity(lines, prior)
    assert ambiguity is not None
    assert ambiguity["excludes_deuterium"] is True
    assert ambiguity["calendar"] == "H-only"

    plain, _summary = verdict_from_evidence(lines)
    verdict, summary = verdict_from_evidence(lines, isotope_ambiguity=ambiguity)

    # Judged on the residuals alone this reads as an ordinary small shift, which
    # is precisely the silence the rehearsal caught.
    assert plain == "shifted"
    assert verdict == ERA_MISASSIGNED_VERDICT
    assert summary["verdict_reason"] == ambiguity["finding"]
    assert summary["isotope_ambiguity"]["calendar_expectation"] == "hydrogen expected"


def test_a_mixed_isotope_shot_is_not_an_era_shift_and_raises_nothing() -> None:
    """One hydrogen assignment refutes the rigid-shift hypothesis outright.

    A translation of the detector moves every line by the same pixels; it
    cannot move two lines onto deuterium and leave a third on hydrogen.  So a
    mixed shot is mixed, not ambiguous, and gets no finding.
    """

    lines = _flipped_2019()
    lines[-1] = {**lines[-1], "isotope": "H"}

    assert isotope_ambiguity(lines, deuterium_prior(FLIPPED_SHOT_DATE, "2019")) is None


def test_a_flip_whose_own_fit_is_poor_is_a_shift_question_not_an_isotope_one() -> None:
    """Outside the degeneracy window neither reading fits, so nothing is degenerate."""

    lines = [
        {**item, "pixel_residual_px": item["pixel_residual_px"] + 9.0}
        for item in _flipped_2019()
    ]
    for item in lines:
        item["isotope_candidates"] = [
            {**entry, "pixel_residual_px": entry["pixel_residual_px"] + 9.0}
            for entry in item["isotope_candidates"]
        ]

    assert isotope_ambiguity(lines, deuterium_prior(FLIPPED_SHOT_DATE, "2019")) is None


# ---------------------------------------------------------------------------
# What the words now mean
# ---------------------------------------------------------------------------


def test_beyond_repair_needs_the_disagreement_to_be_large_as_well_as_real() -> None:
    """Damage is inconsistency plus size; inconsistency alone is quality.

    Strong lines that disagree by more than one pixel but sit within
    GEOMETRIC_DAMAGE_PX of zero are still worth a correction: sliding the table
    lands every one of them inside about a pixel.  Condemning them would send
    an operator back to the raw SIF for something a table shift repairs.
    """

    near = [
        _line(order=order, expected_nm=nm, residual_px=residual, snr=200.0, label="H", shot="s")
        for order, nm, residual in ((6, 656.2790, 1.4), (12, 486.1350, -0.2), (16, 434.0470, 0.9))
    ]
    far = [
        _line(order=order, expected_nm=nm, residual_px=residual, snr=200.0, label="H", shot="s")
        for order, nm, residual in ((6, 656.2790, 2.4), (12, 486.1350, -2.0), (16, 434.0470, 1.9))
    ]

    near_verdict, near_summary = verdict_from_evidence(near)
    far_verdict, far_summary = verdict_from_evidence(far)

    assert near_summary["maximum_pixel_deviation_px"] > SHIFT_CONSISTENCY_PX
    assert near_summary["maximum_absolute_pixel_residual_px"] < GEOMETRIC_DAMAGE_PX
    assert near_verdict == "shifted"
    assert "inside the" in near_summary["verdict_reason"]

    assert far_summary["maximum_absolute_pixel_residual_px"] > GEOMETRIC_DAMAGE_PX
    assert far_verdict == "misaligned-beyond-repair"
    assert "the detector geometry" in far_summary["verdict_reason"]


def test_a_consistent_shift_past_the_repair_limit_says_why_it_cannot_be_slid() -> None:
    lines = [
        _line(order=order, expected_nm=nm, residual_px=31.0, snr=200.0, label="H", shot="s")
        for order, nm in ((6, 656.2790), (12, 486.1350), (16, 434.0470))
    ]

    verdict, summary = verdict_from_evidence(lines)

    assert verdict == "misaligned-beyond-repair"
    assert summary["strong_lines_inconsistent"] is False
    assert "re-identified against lamp data" in summary["verdict_reason"]
