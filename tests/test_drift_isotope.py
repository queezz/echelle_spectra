"""Packet F10 — the drift judge knows which hydrogen it is looking at.

LHD ran deuterium from March 2017 to December 2022, and those are the shots
queezz cares most about.  D-alpha sits 0.178 nm blueward of H-alpha, inside the
+/-0.4 nm window this audit fits, so a perfectly calibrated deuterium cube
judged against hydrogen alone reads as a confident ~16.5 px shift and invites a
refinement wrong by exactly one isotope.

Every cube here is built by putting Gaussians at the wavelengths a real plasma
would put them: the deuterium fixtures are lit at the *published D I* air
wavelengths on a calibration that is exactly correct.  Nothing is shifted.  A
fixture that moved the detector instead would test Packet F4 again.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from echelle_spectra.drift import (
    ALIGNMENT_MAX_PIXEL_RESIDUAL,
    DEUTERIUM_CALENDAR_RESOURCE,
    ERA_MISASSIGNED_VERDICT,
    REPAIR_LIMIT_PX,
    DriftError,
    audit_cubes,
    deuterium_prior,
    load_deuterium_calendar,
    require_sampled_verdict,
    verdict_from_evidence,
    write_drift_evidence,
)
from echelle_spectra.tools.line_catalog import (
    LINE_FAMILY_ISOTOPES,
    load_line_table,
)

# The four audited Balmer orders, plus one order over the Fulcher Q branch so a
# shot can be judged with and without the H2 anchors.
BALMER_ORDERS: tuple[tuple[int, float], ...] = (
    (6, 656.2790),
    (12, 486.1350),
    (16, 434.0470),
    (18, 410.1734),
)
FULCHER_ORDER: tuple[int, float] = (9, 626.0)
DETECTOR_WIDTH_PX = 2304
CENTRE_PX = (DETECTOR_WIDTH_PX - 1) / 2.0
ALPHA_DISPERSION_NM_PER_PX = 0.0108
CURVATURE_NM_PER_PX2 = 5e-8
LINE_SIGMA_NM = 0.035

# The catalog's own numbers, restated here so a silent edit to either table is
# a test failure rather than a quiet change of what "deuterium" means.
HYDROGEN_NM = {"3-2": 656.2790, "4-2": 486.1350, "5-2": 434.0470, "6-2": 410.1734}
DEUTERIUM_NM = {"3-2": 656.1012, "4-2": 486.0001, "5-2": 433.9283, "6-2": 410.0619}
TRANSITION_ORDER = {"3-2": 6, "4-2": 12, "5-2": 16, "6-2": 18}

# Resolved Fulcher rows inside the 626 nm order, each clear of its neighbours by
# more than the +/-0.4 nm fit window.
FULCHER_NM = (618.6530, 621.4507, 628.9709)

INSIDE_THE_DEUTERIUM_PHASE = "2018-05-15T09:00:00+09:00"
AFTER_THE_DEUTERIUM_PHASE = "2025-06-01T09:00:00+09:00"


def _dispersion(centre_nm: float) -> float:
    return ALPHA_DISPERSION_NM_PER_PX * centre_nm / 656.2790


def _polynomial(centre_nm: float) -> np.ndarray:
    slope = _dispersion(centre_nm)
    return np.array(
        [
            CURVATURE_NM_PER_PX2,
            slope - 2.0 * CURVATURE_NM_PER_PX2 * CENTRE_PX,
            centre_nm - slope * CENTRE_PX + CURVATURE_NM_PER_PX2 * CENTRE_PX**2,
        ]
    )


def _polynomials(orders: Sequence[tuple[int, float]]) -> dict[int, np.ndarray]:
    return {order: _polynomial(centre) for order, centre in orders}


def _axis(polynomials: dict[int, np.ndarray]):
    pixels = np.arange(DETECTOR_WIDTH_PX, dtype=float)
    wavelength = np.concatenate([np.polyval(polynomials[o], pixels) for o in polynomials])
    detector_pixel = np.tile(pixels, len(polynomials))
    echelle_order = np.repeat(list(polynomials), DETECTOR_WIDTH_PX)
    permutation = np.argsort(wavelength, kind="stable")
    return wavelength[permutation], detector_pixel[permutation], echelle_order[permutation]


def _polynomial_payload(polynomials: dict[int, np.ndarray]) -> str:
    return json.dumps(
        {
            "schema": "spectrocube.wavelength-polynomials/v1",
            "coefficient_order": "descending_power",
            "input": "detector_pixel",
            "input_units": "pixel",
            "output": "wavelength",
            "output_units": "nm",
            "orders": [
                {"order": order, "coefficients": coefficients.tolist()}
                for order, coefficients in sorted(polynomials.items())
            ],
        },
        sort_keys=True,
    )


def _cube(
    path: Path,
    *,
    lit_nm: Sequence[float],
    orders: Sequence[tuple[int, float]] = BALMER_ORDERS,
    shot_number: str = "133500",
    t_start: str = INSIDE_THE_DEUTERIUM_PHASE,
) -> Path:
    """Write a cube whose calibration is correct and whose lines sit at ``lit_nm``."""

    polynomials = _polynomials(orders)
    wavelength, detector_pixel, echelle_order = _axis(polynomials)
    rng = np.random.default_rng(7)
    values = rng.normal(0.0, 0.02, size=wavelength.size)
    for centre in lit_nm:
        values += 50.0 * np.exp(-0.5 * ((wavelength - float(centre)) / LINE_SIGMA_NM) ** 2)
    dataset = xr.Dataset(
        {"intensity": (("frame", "wavelength"), values[np.newaxis, :])},
        coords={
            "frame": np.zeros(1, dtype=float),
            "wavelength": wavelength,
            "detector_pixel": ("wavelength", detector_pixel),
            "echelle_order": ("wavelength", echelle_order.astype(np.int64)),
        },
        attrs={
            "spectrocube_version": "0.2.0",
            "instrument_id": "echelle",
            "calibration_type": "counts",
            "intensity_units": "counts",
            "wavelength_medium": "air",
            "created_at": "2025-01-01T00:00:00+00:00",
            "snapshot_id": "20170301_ccd",
            "shot_number": shot_number,
            "t_start": t_start,
            "wavelength_polynomials_json": _polynomial_payload(polynomials),
        },
    )
    dataset.coords["wavelength"].attrs.update(units="nm", medium="air")
    dataset.coords["detector_pixel"].attrs.update(
        units="pixel", detector_axis="column", reference_frame="raw_detector", index_origin=0
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path)
    return path


def _balmer(payload: dict) -> list[dict]:
    return [
        item
        for item in payload["lines"]
        if item.get("family") == "balmer" and item.get("status") == "measured"
    ]


def _fulcher(payload: dict) -> list[dict]:
    return [item for item in payload["lines"] if item.get("family") == "fulcher"]


# ---------------------------------------------------------------------------
# The trap: a pure deuterium shot on a correct calibration
# ---------------------------------------------------------------------------


def test_pure_deuterium_cube_on_a_correct_calibration_reads_aligned(tmp_path: Path) -> None:
    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133500.nc", lit_nm=tuple(DEUTERIUM_NM.values()))],
        families=("balmer",),
    )

    assert payload["verdict"] == "aligned"
    assert payload["sampled_cubes"][0]["isotope"] == "D"
    assert [item["isotope"] for item in payload["per_shot"]] == ["D"]
    assert payload["summary"]["isotope_tags"] == {"D": 1}

    measured = _balmer(payload)
    assert len(measured) == len(DEUTERIUM_NM)
    for item in measured:
        assert item["isotope"] == "D"
        assert item["line"].startswith("D-")
        assert item["expected_nm"] == pytest.approx(DEUTERIUM_NM[_transition_of(item)])
        assert abs(item["pixel_residual_px"]) < ALIGNMENT_MAX_PIXEL_RESIDUAL


def _transition_of(item: dict) -> str:
    """Recover a measured line's transition from the order it was found in."""

    return next(
        transition
        for transition, order in TRANSITION_ORDER.items()
        if order == item["echelle_order"]
    )


def test_the_hydrogen_only_reading_of_that_deuterium_cube_is_a_false_shifted(
    tmp_path: Path,
) -> None:
    """Non-vacuity: the same evidence, judged as before F10, condemns a good cube.

    Every measured record keeps both candidates, so the pre-F10 reading can be
    reconstructed exactly: take each window's hydrogen candidate and run the
    unchanged verdict arithmetic over it.  It returns ``shifted`` at ~-16.5 px,
    inside the repair limit — a confident, wrong verdict that would have
    produced a refinement snapshot off by one isotope.  If this test ever
    reported ``aligned``, the trap above would be proving nothing.
    """

    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133501.nc", lit_nm=tuple(DEUTERIUM_NM.values()))],
        families=("balmer",),
    )

    hydrogen_only = []
    for item in _balmer(payload):
        candidates = {entry["isotope"]: entry for entry in item["isotope_candidates"]}
        assert set(candidates) == {"H", "D"}
        hydrogen_only.append({**item, **candidates["H"]})

    verdict, summary = verdict_from_evidence(hydrogen_only)

    assert verdict == "shifted"
    assert summary["median_shift_px"] == pytest.approx(-16.5, abs=0.5)
    assert abs(summary["median_shift_px"]) < REPAIR_LIMIT_PX
    assert summary["maximum_pixel_deviation_px"] < 1.0


def test_the_isotope_offset_is_one_rigid_pixel_offset_in_every_order() -> None:
    """Why nearest-assignment cannot be finer than it is, stated as a number.

    The D/H wavelength gap and the dispersion both grow with wavelength, so the
    gap is the same ~16.5 px in the blue as at H-alpha.  That is what makes the
    isotope ambiguity exactly degenerate with a rigid detector shift, and why
    the calendar prior flags rather than resolves.
    """

    offsets = [
        (HYDROGEN_NM[key] - DEUTERIUM_NM[key]) / _dispersion(HYDROGEN_NM[key])
        for key in HYDROGEN_NM
    ]

    assert min(offsets) == pytest.approx(16.5, abs=0.5)
    assert max(offsets) - min(offsets) < 0.5


# ---------------------------------------------------------------------------
# Per-line assignment
# ---------------------------------------------------------------------------


def test_a_mixed_shot_assigns_each_line_to_its_nearest_isotope(tmp_path: Path) -> None:
    lit = [DEUTERIUM_NM["3-2"], DEUTERIUM_NM["4-2"], DEUTERIUM_NM["5-2"], HYDROGEN_NM["6-2"]]

    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133600.nc", lit_nm=lit)], families=("balmer",)
    )

    assigned = {_transition_of(item): item["isotope"] for item in _balmer(payload)}
    assert assigned == {"3-2": "D", "4-2": "D", "5-2": "D", "6-2": "H"}
    assert payload["sampled_cubes"][0]["isotope"] == "D"
    assert payload["verdict"] == "aligned"

    # Each line is judged against its own reference, so all four residuals are
    # zero at once — which a single-isotope audit can never report.
    for item in _balmer(payload):
        assert abs(item["pixel_residual_px"]) < ALIGNMENT_MAX_PIXEL_RESIDUAL


def test_an_even_split_is_reported_as_mixed_rather_than_broken(tmp_path: Path) -> None:
    lit = [DEUTERIUM_NM["3-2"], DEUTERIUM_NM["4-2"], HYDROGEN_NM["5-2"], HYDROGEN_NM["6-2"]]

    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133601.nc", lit_nm=lit)], families=("balmer",)
    )

    assert payload["sampled_cubes"][0]["isotope"] == "mixed"
    assert payload["summary"]["isotope_tags"] == {"mixed": 1}


# ---------------------------------------------------------------------------
# The H2 anchors on a deuterium shot
# ---------------------------------------------------------------------------


def test_a_hydrogen_shot_keeps_the_fulcher_anchors_and_is_tagged_hydrogen(
    tmp_path: Path,
) -> None:
    payload = audit_cubes(
        [
            _cube(
                tmp_path / "cubes" / "133700.nc",
                lit_nm=tuple(HYDROGEN_NM.values()) + FULCHER_NM,
                orders=(*BALMER_ORDERS, FULCHER_ORDER),
            )
        ]
    )

    assert payload["sampled_cubes"][0]["isotope"] == "H"
    assert all(item["isotope"] == "H" for item in _balmer(payload))
    measured = [item for item in _fulcher(payload) if item.get("status") == "measured"]
    assert {round(item["expected_nm"], 4) for item in measured} == set(FULCHER_NM)
    assert all(item["isotope"] == "H" for item in measured)
    assert payload["summary"]["isotope_excluded_lines"] == 0


def test_a_deuterium_shot_drops_the_fulcher_anchors_from_the_quorum(
    tmp_path: Path,
) -> None:
    payload = audit_cubes(
        [
            _cube(
                tmp_path / "cubes" / "133701.nc",
                lit_nm=tuple(DEUTERIUM_NM.values()) + FULCHER_NM,
                orders=(*BALMER_ORDERS, FULCHER_ORDER),
            )
        ]
    )

    assert payload["sampled_cubes"][0]["isotope"] == "D"
    fulcher = _fulcher(payload)
    assert fulcher, "the fixture must cover the Fulcher band for this to mean anything"
    assert all(item["status"] == "skipped" for item in fulcher)
    assert all(item["isotope_excluded"] for item in fulcher)
    assert all("no D2 table exists here" in item["reason"] for item in fulcher)
    assert payload["summary"]["isotope_excluded_lines"] == len(fulcher)

    # The dropped lines reach neither the quorum nor the fitted shift.
    assert payload["summary"]["quorum"]["orders"] == [
        order for order, _centre in BALMER_ORDERS
    ]
    assert FULCHER_ORDER[0] not in payload["summary"]["quorum"]["orders"]
    assert payload["verdict"] == "aligned"


# ---------------------------------------------------------------------------
# The calendar, as a prior
# ---------------------------------------------------------------------------


def test_a_deuterium_reading_outside_the_calendar_window_flags_the_shot(
    tmp_path: Path,
) -> None:
    payload = audit_cubes(
        [
            _cube(
                tmp_path / "cubes" / "190000.nc",
                lit_nm=tuple(DEUTERIUM_NM.values()),
                shot_number="190000",
                t_start=AFTER_THE_DEUTERIUM_PHASE,
            )
        ],
        families=("balmer",),
    )

    shot = payload["sampled_cubes"][0]
    assert shot["isotope_prior"] == "hydrogen expected"
    assert shot["isotope"] == "D"
    assert "outside every deuterium window" in shot["isotope_prior_basis"]
    assert "expects hydrogen" in shot["isotope_flag"]
    assert payload["summary"]["isotope_flagged_shots"] == ["190000"]
    # CHANGED after the 2026-08-18 real-light rehearsal.  This used to assert
    # ``aligned``: the flag was advisory and the verdict was "the one the
    # residuals earned".  That night showed what earning it costs — a 2019 cube
    # recalibrated onto the wrong era had its hydrogen lines re-read as
    # deuterium and the flipped fit scored *better* than the correct one, so an
    # aligned verdict here would authorize a bulk run over cubes processed
    # against the wrong epoch.  Every measurement below is untouched: the
    # assignment is still D, the residuals are still sub-pixel against D.  Only
    # the interval's authorization changed, and it changed to a word that names
    # the cause.
    assert payload["verdict"] == ERA_MISASSIGNED_VERDICT
    assert all(item["isotope"] == "D" for item in _balmer(payload))
    for item in _balmer(payload):
        assert abs(item["pixel_residual_px"]) < ALIGNMENT_MAX_PIXEL_RESIDUAL

    ambiguity = payload["summary"]["isotope_ambiguity"]
    assert ambiguity["finding"] == "isotope flip or era shift, degenerate; calendar says H-only"
    assert ambiguity["excludes_deuterium"] is True
    assert ambiguity["implied_common_shift_px"] == pytest.approx(-16.5, abs=0.5)
    assert ambiguity["isotope_separation_px"] == pytest.approx(16.5, abs=0.5)
    assert payload["summary"]["isotope_ambiguous_shots"] == ["190000"]
    assert "wrong epoch" in ambiguity["detail"]
    # The repair is a cross-era recalibration, so the advice names the full
    # 'recal-cube' and warns off the --wavelength-only form the shifted verdict
    # composes: that flag reuses the base snapshot's sphere pair, and another
    # era's snapshot does not share it.
    assert "190000" in payload["verdict_advice"]
    assert "NOT the --wavelength-only form" in payload["verdict_advice"]

    # And it is refused at the gate rather than merely displayed.
    evidence = write_drift_evidence(tmp_path / "flipped.json", payload)
    with pytest.raises(DriftError, match="recalibrate these cubes"):
        require_sampled_verdict(evidence, {"20170301_ccd"})


def test_a_deuterium_reading_inside_the_window_is_expected_and_unflagged(
    tmp_path: Path,
) -> None:
    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133800.nc", lit_nm=tuple(DEUTERIUM_NM.values()))],
        families=("balmer",),
    )

    shot = payload["sampled_cubes"][0]
    assert shot["isotope_prior"] == "deuterium possible"
    assert "LHD deuterium phase" in shot["isotope_prior_basis"]
    assert "isotope_flag" not in shot
    assert payload["summary"]["isotope_flagged_shots"] == []


def test_hydrogen_inside_a_deuterium_window_is_not_a_disagreement(
    tmp_path: Path,
) -> None:
    """LHD ran hydrogen shots inside its deuterium cycles; that is not news."""

    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133801.nc", lit_nm=tuple(HYDROGEN_NM.values()))],
        families=("balmer",),
    )

    shot = payload["sampled_cubes"][0]
    assert shot["isotope_prior"] == "deuterium possible"
    assert shot["isotope"] == "H"
    assert "isotope_flag" not in shot


def test_the_prior_places_shots_by_date_and_refuses_to_invent_a_ceiling() -> None:
    from datetime import date

    assert deuterium_prior(date(2017, 3, 7))["expectation"] == "deuterium possible"
    assert deuterium_prior(date(2022, 12, 2))["expectation"] == "deuterium possible"
    assert deuterium_prior(date(2017, 3, 6))["expectation"] == "hydrogen expected"
    assert deuterium_prior(date(2022, 12, 3))["expectation"] == "hydrogen expected"

    # Without a date, a shot below the first deuterium shot is still safely
    # hydrogen; one above it cannot be placed, because the calendar publishes no
    # closing shot number and guessing one would be a fabricated prior.
    assert deuterium_prior(None, "133269")["expectation"] == "hydrogen expected"
    assert deuterium_prior(None, "133270")["expectation"] == "unknown"
    assert deuterium_prior(None, "")["expectation"] == "unknown"


def test_the_bundled_calendar_parses_and_carries_its_citations() -> None:
    calendar = load_deuterium_calendar()

    assert calendar["schema"] == "echelle-deuterium-calendar/v1"
    assert calendar["facility"] == "LHD"
    window = calendar["windows"][0]
    assert window["date_from"] == "2017-03-07"
    assert window["date_to"] == "2022-12-02"
    assert window["cycles"] == 6
    assert window["shot_from"] == 133270
    assert calendar["campaigns"][0]["shot_from"] == 133270

    from importlib.resources import files

    text = (
        files("echelle_spectra.resources")
        .joinpath("lhd_deuterium_campaign.toml")
        .read_text(encoding="utf-8")
    )
    for citation in (
        "10.1016/j.fusengdes.2019.03.062",
        "10.1088/1741-4326/ac3cda",
        "10.1007/s10894-025-00532-0",
    ):
        assert citation in text


def test_the_evidence_names_the_calendar_it_used(tmp_path: Path) -> None:
    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "133802.nc", lit_nm=tuple(HYDROGEN_NM.values()))],
        families=("balmer",),
    )

    prior = payload["isotope_prior"]
    assert prior["calendar"] == DEUTERIUM_CALENDAR_RESOURCE
    assert prior["facility"] == "LHD"
    assert "prior only" in prior["role"]
    assert prior["windows"][0]["date_from"] == "2017-03-07"


# ---------------------------------------------------------------------------
# The catalog facet a future D2 table slots into
# ---------------------------------------------------------------------------


def test_the_balmer_table_carries_both_isotopologues_with_provenance() -> None:
    hydrogen = load_line_table("balmer", isotope="H")
    deuterium = load_line_table("balmer", isotope="D")

    assert load_line_table("balmer") == hydrogen, "the default table must not move"
    assert {line.transition for line in deuterium} == set(DEUTERIUM_NM)
    assert {line.species for line in deuterium} == {"D I"}
    for line in deuterium:
        assert line.isotope == "D"
        assert line.wavelength_medium == "air"
        assert line.wavelength_nm == pytest.approx(DEUTERIUM_NM[line.transition])
        assert "NIST" in line.source_name
        assert "physics.nist.gov" in line.source_reference
    assert {line.transition for line in hydrogen} == {line.transition for line in deuterium}


def test_the_fulcher_table_answers_deuterium_with_nothing_rather_than_h2() -> None:
    assert LINE_FAMILY_ISOTOPES["fulcher"] == ("H",)
    assert load_line_table("fulcher", isotope="D") == ()
    assert load_line_table("fulcher", isotope="H") == load_line_table("fulcher")
    assert all(line.isotope == "H" for line in load_line_table("fulcher"))


def test_a_lamp_family_refuses_an_isotope_it_cannot_have() -> None:
    with pytest.raises(ValueError, match="no hydrogen isotopologue"):
        load_line_table("ne", isotope="D")
