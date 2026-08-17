"""Packet F4 — the drift judge measures detector pixels, not nanometres.

Every cube here is built the way the instrument builds one: four Echelle orders
with their own second-order wavelength polynomials, a dispersion running from
0.0068 nm/px in the blue to 0.0108 nm/px at H-alpha (a 60% spread), and lines
placed by moving the *detector* rather than by adding one wavelength offset to
everything.  A fixture that shifted all lines by a constant wavelength would
only reproduce the model this packet replaced.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from echelle_spectra.catalog import build_drive_catalog, merge_catalogs
from echelle_spectra.drift import (
    ALIGNMENT_MAX_PIXEL_RESIDUAL,
    DRIFT_SCHEMA,
    DriftError,
    audit_cubes,
    create_refinement_snapshot,
    select_sample_paths,
    write_drift_evidence,
)
from echelle_spectra.recalibration import recalibrate_cube
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot

# One Balmer line per order, with the shipped instrument's dispersion scale.
ORDERS: tuple[tuple[int, float], ...] = (
    (6, 656.2790),  # H-alpha
    (12, 486.1350),  # H-beta
    (16, 434.0470),  # H-gamma
    (18, 410.1734),  # H-delta
)
DETECTOR_WIDTH_PX = 2304
CENTRE_PX = (DETECTOR_WIDTH_PX - 1) / 2.0
ALPHA_DISPERSION_NM_PER_PX = 0.0108
CURVATURE_NM_PER_PX2 = 5e-8
LINE_SIGMA_NM = 0.035
ANCHOR_PIXELS = (100.0, 600.0, 1152.0, 1700.0, 2200.0)
SNAPSHOT_ID = "20250101_cmos"


def _dispersion(centre_nm: float) -> float:
    """Return the order's mid-detector dispersion, scaled with wavelength."""

    return ALPHA_DISPERSION_NM_PER_PX * centre_nm / 656.2790


def _polynomial(centre_nm: float) -> np.ndarray:
    """Return descending-power coefficients of one order's lambda(pixel)."""

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


POLYNOMIALS = _polynomials(ORDERS)


def _pixel_at(polynomials: dict[int, np.ndarray], order: int, wavelength_nm: float) -> float:
    coefficients = polynomials[order].copy()
    coefficients[-1] -= wavelength_nm
    roots = [float(value.real) for value in np.roots(coefficients) if abs(value.imag) < 1e-9]
    inside = [value for value in roots if -50.0 <= value <= DETECTOR_WIDTH_PX + 50.0]
    assert inside, f"order {order} does not cover {wavelength_nm} nm"
    return inside[0]


def _axis(
    polynomials: dict[int, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the stitched, ascending wavelength axis and its raw geometry."""

    pixels = np.arange(DETECTOR_WIDTH_PX, dtype=float)
    wavelength = np.concatenate(
        [np.polyval(polynomials[order], pixels) for order in polynomials]
    )
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


def _frame(
    wavelength: np.ndarray,
    polynomials: dict[int, np.ndarray],
    lines: Sequence[tuple[int, float]],
    shifts_px: Sequence[float],
    amplitude: float,
    *,
    seed: int,
) -> np.ndarray:
    """Emit one frame whose lines sit where a shifted detector would put them."""

    rng = np.random.default_rng(seed)
    values = rng.normal(0.0, 0.02, size=wavelength.size)
    if amplitude <= 0:
        return values
    for (order, centre), shift in zip(lines, shifts_px):
        observed = float(
            np.polyval(polynomials[order], _pixel_at(polynomials, order, centre) + shift)
        )
        values += amplitude * np.exp(-0.5 * ((wavelength - observed) / LINE_SIGMA_NM) ** 2)
    return values


def _cube(
    path: Path,
    *,
    shifts_px: float | Sequence[float] = 0.0,
    shot_number: str = "100100",
    snapshot_id: str = SNAPSHOT_ID,
    provenance: dict[str, str] | None = None,
    frame_specs: Sequence[tuple[Sequence[float] | float, float]] | None = None,
    orders: Sequence[tuple[int, float]] = ORDERS,
    lines: Sequence[tuple[int, float]] | None = None,
    t_start: str = "",
    accuracy_nm: float | None = None,
) -> Path:
    """Write one SpectroCube-shaped fixture cube with real detector geometry."""

    polynomials = _polynomials(orders)
    lines = list(lines if lines is not None else orders)
    wavelength, detector_pixel, echelle_order = _axis(polynomials)
    if not isinstance(shifts_px, Sequence):
        shifts_px = [float(shifts_px)] * len(lines)
    specs = list(frame_specs) if frame_specs is not None else [(shifts_px, 50.0)]
    frames = np.vstack(
        [
            _frame(
                wavelength,
                polynomials,
                lines,
                [float(shift)] * len(lines) if not isinstance(shift, Sequence) else shift,
                amplitude,
                seed=index,
            )
            for index, (shift, amplitude) in enumerate(specs)
        ]
    )
    attrs = {
        "spectrocube_version": "0.2.0",
        "instrument_id": "echelle",
        "calibration_type": "counts",
        "intensity_units": "counts",
        "wavelength_medium": "air",
        "created_at": "2025-01-01T00:00:00+00:00",
        "snapshot_id": snapshot_id,
        "shot_number": shot_number,
        "wavelength_polynomials_json": _polynomial_payload(polynomials),
        **(provenance or {}),
    }
    if t_start:
        attrs["t_start"] = t_start
    if accuracy_nm is not None:
        attrs["wavelength_accuracy_nm"] = float(accuracy_nm)
    dataset = xr.Dataset(
        {"intensity": (("frame", "wavelength"), frames)},
        coords={
            "frame": np.arange(frames.shape[0], dtype=float),
            "wavelength": wavelength,
            "detector_pixel": ("wavelength", detector_pixel),
            "echelle_order": ("wavelength", echelle_order.astype(np.int64)),
        },
        attrs=attrs,
    )
    dataset.coords["wavelength"].attrs.update(units="nm", medium="air")
    dataset.coords["detector_pixel"].attrs.update(
        units="pixel", detector_axis="column", reference_frame="raw_detector", index_origin=0
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path)
    return path


def _snapshot(root: Path, snapshot_id: str = SNAPSHOT_ID, **manifest):
    """Create a snapshot whose wavelength anchors sit on the fixture polynomials."""

    sources = root / f"sources-{snapshot_id}"
    sources.mkdir(parents=True)
    files = {}
    for role in ROLE_FILENAMES:
        source = sources / f"{role}.dat"
        if role == "wavelength":
            rows = [
                f"{order}\t{pixel - 8:09.3f}\t{pixel + 8:09.3f}\t{pixel:010.4f}\t"
                f"{float(np.polyval(POLYNOMIALS[order], pixel)):014.9f}\tHI"
                for order, _ in ORDERS
                for pixel in ANCHOR_PIXELS
            ]
            source.write_text(
                "# order from to center wavelength species\n" + "\n".join(rows) + "\n",
                encoding="utf-8",
            )
        else:
            source.write_text(f"{role} fixture\n", encoding="utf-8")
        files[role] = source
    return create_snapshot(
        root / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity={"shot_from": 1, "shot_to": 999999},
        **manifest,
    )


def _measured(payload: dict) -> list[dict]:
    return [item for item in payload["lines"] if item.get("status") == "measured"]


# ---------------------------------------------------------------------------
# Detector-space classification
# ---------------------------------------------------------------------------


def test_rigid_five_pixel_shift_is_shifted_and_recovered_in_pixels(tmp_path: Path) -> None:
    payload = audit_cubes([_cube(tmp_path / "cubes" / "100100.nc", shifts_px=5.0)])

    assert payload["schema"] == DRIFT_SCHEMA
    assert payload["verdict"] == "shifted"
    assert payload["summary"]["median_shift_px"] == pytest.approx(5.0, abs=0.2)
    assert payload["summary"]["maximum_pixel_deviation_px"] < 1.0

    measured = _measured(payload)
    assert {item["echelle_order"] for item in measured} == {order for order, _ in ORDERS}
    # The wavelength residuals are NOT constant: that is the defect being fixed.
    residuals_nm = sorted(abs(item["residual_nm"]) for item in measured)
    assert residuals_nm[-1] / residuals_nm[0] == pytest.approx(
        _dispersion(656.2790) / _dispersion(410.1734), rel=0.05
    )
    for item in measured:
        assert item["pixel_residual_px"] == pytest.approx(5.0, abs=0.2)

    # The emitted correction is per-order, through each order's own dispersion.
    corrections = {item["order"]: item for item in payload["order_corrections"]}
    assert set(corrections) == {order for order, _ in ORDERS}
    for order, centre in ORDERS:
        assert corrections[order]["predicted_shift_nm"] == pytest.approx(
            -5.0 * _dispersion(centre), rel=0.02
        )


def test_eighteen_pixel_shift_stays_repairable(tmp_path: Path) -> None:
    """The band the wavelength-space judge misread as beyond repair."""

    payload = audit_cubes([_cube(tmp_path / "cubes" / "100101.nc", shifts_px=18.0)])

    assert payload["verdict"] == "shifted"
    assert payload["summary"]["median_shift_px"] == pytest.approx(18.0, abs=0.3)
    # In wavelength the same residuals scatter by more than the retired 0.04 nm
    # consistency bound, which is exactly why they used to be condemned.
    residuals_nm = np.array([item["residual_nm"] for item in _measured(payload)])
    assert np.max(np.abs(residuals_nm - np.median(residuals_nm))) > 0.04


def test_inconsistent_residuals_are_beyond_repair(tmp_path: Path) -> None:
    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "100102.nc", shifts_px=[2.0, -2.0, 1.6, -1.9])]
    )

    assert payload["verdict"] == "misaligned-beyond-repair"
    assert payload["summary"]["maximum_pixel_deviation_px"] > 1.0
    assert not payload["repair_commands"]


def test_aligned_cube_reports_no_repair(tmp_path: Path) -> None:
    payload = audit_cubes([_cube(tmp_path / "cubes" / "100103.nc", shifts_px=0.0)])

    assert payload["verdict"] == "aligned"
    assert payload["summary"]["maximum_absolute_pixel_residual_px"] < ALIGNMENT_MAX_PIXEL_RESIDUAL
    assert payload["repair_command"] == ""


# ---------------------------------------------------------------------------
# The round trip that kills the under-correction
# ---------------------------------------------------------------------------


def test_accepted_shift_recalibrates_every_order_and_re_audits_aligned(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path)
    cube = _cube(
        tmp_path / "cubes" / "100200.nc",
        shifts_px=5.0,
        provenance=snapshot.provenance_attrs(),
    )
    evidence_path = tmp_path / "drift.json"
    payload = audit_cubes(
        [cube], evidence_path=evidence_path, calibrations_root=tmp_path / "calibrations"
    )
    assert payload["verdict"] == "shifted"
    assert payload["refined_snapshot_id"] == f"{SNAPSHOT_ID}-r1"
    assert str(evidence_path) in payload["repair_command"]
    assert "DRIFT_EVIDENCE.json" not in json.dumps(payload)
    shells = {step["shell"] for step in payload["repair_commands"]}
    assert {"any", "powershell", "posix"} <= shells
    recal = [step for step in payload["repair_commands"] if "recal-cube" in step["command"]]
    assert len(recal) == 2
    # The refinement changed only the wavelength role, so that is all the
    # composed revision claims — and it is what a counts cube can even accept.
    assert all("--wavelength-only" in step["command"] for step in recal)
    assert all(f"{SNAPSHOT_ID}-r1" in step["command"] for step in recal)
    written = write_drift_evidence(evidence_path, payload)

    # Paste the composed command's own shift back: a rounded one would be
    # refused by the exact acknowledgement it is composed for.
    tokens = payload["repair_command"].split()
    pasted_shift = float(tokens[tokens.index("--accept-shift") + 1])
    refined, accepted = create_refinement_snapshot(
        written,
        calibrations_root=tmp_path / "calibrations",
        accepted_shift_px=pasted_shift,
    )
    assert refined.snapshot_id == f"{SNAPSHOT_ID}-r1"
    assert refined.manifest["alignment"]["applied_pixel_shift_px"] == pytest.approx(
        payload["summary"]["median_shift_px"]
    )
    assert accepted.is_file()
    # A refinement re-solves the wavelength table and nothing else.  Its shift
    # was measured through the base snapshot's own order pattern, so carrying
    # that pattern across untouched is what keeps the measurement and the file
    # it was measured in agreement — unlike a bench save, which solves a rigid
    # transform against a reference and must move both halves together.
    assert (
        refined.source_path("pattern").read_bytes()
        == snapshot.source_path("pattern").read_bytes()
    )

    revised, _manifest = recalibrate_cube(
        cube,
        tmp_path / "cubes" / "100200-r1.nc",
        new_snapshot_path=refined.root,
        update_factor=False,
    )

    # What the audit predicted per order is what recalibration actually applied.
    with xr.open_dataset(cube) as before, xr.open_dataset(revised) as after:
        for correction in payload["order_corrections"]:
            order = correction["order"]
            selected = before["echelle_order"].values == order
            index = int(
                np.argmin(np.abs(before["detector_pixel"].values[selected] - correction["reference_pixel"]))
            )
            applied = (
                after["wavelength"].values[selected][index]
                - before["wavelength"].values[selected][index]
            )
            # 1e-5 nm is 1.5e-3 px, the refit's own conditioning: four orders
            # below the 0.01 nm spread between these per-order corrections.
            assert applied == pytest.approx(correction["predicted_shift_nm"], abs=1e-5)

    re_audit = audit_cubes([revised])
    assert re_audit["verdict"] == "aligned"
    for item in _measured(re_audit):
        assert abs(item["pixel_residual_px"]) < ALIGNMENT_MAX_PIXEL_RESIDUAL


def test_v1_evidence_cannot_be_accepted_as_a_pixel_shift(tmp_path: Path) -> None:
    _snapshot(tmp_path)
    legacy = write_drift_evidence(
        tmp_path / "legacy.json",
        {
            "schema": "echelle-drift-evidence/v1",
            "verdict": "shifted",
            "snapshot_ids": [SNAPSHOT_ID],
            "summary": {"median_shift_nm": 0.05},
        },
    )
    with pytest.raises(DriftError, match="measure it in pixels"):
        create_refinement_snapshot(
            legacy, calibrations_root=tmp_path / "calibrations", accepted_shift_px=0.05
        )


# ---------------------------------------------------------------------------
# Blends, duplicates, and quorum
# ---------------------------------------------------------------------------

# One order over the Fulcher blend region, lit only where the catalog's
# coincident rows are: the exact 626.2495 pair and one sub-resolution pair.
FULCHER_ORDER: tuple[tuple[int, float], ...] = ((9, 626.0),)
FULCHER_BLEND_LINES: tuple[tuple[int, float], ...] = (
    (9, 626.2495),
    (9, 623.7457),
    (9, 623.8391),
)


def test_blended_and_duplicated_fulcher_rows_cannot_carry_a_verdict(tmp_path: Path) -> None:
    payload = audit_cubes(
        [
            _cube(
                tmp_path / "cubes" / "100300.nc",
                shifts_px=5.0,
                orders=FULCHER_ORDER,
                lines=FULCHER_BLEND_LINES,
            )
        ]
    )

    assert payload["verdict"] == "insufficient-data"
    assert "resolved line" in payload["summary"]["quorum"]["reason"]
    duplicated = [item for item in payload["lines"] if item.get("expected_nm") == 626.2495]
    # The catalog keeps both coincident rows; neither may be measured.
    assert len(duplicated) == 2
    assert {item["status"] for item in duplicated} == {"skipped"}
    for wavelength in (623.7457, 623.8391):
        pair = [item for item in payload["lines"] if item.get("expected_nm") == wavelength]
        assert pair and all(item["status"] == "skipped" for item in pair)
    assert payload["summary"]["blended_lines_skipped"] >= 4


def test_three_resolved_order_diverse_lines_reach_quorum(tmp_path: Path) -> None:
    payload = audit_cubes([_cube(tmp_path / "cubes" / "100301.nc", shifts_px=0.0)])

    quorum = payload["summary"]["quorum"]
    assert quorum["satisfied"] is True
    assert quorum["resolved_lines"] >= 3
    assert quorum["distinct_orders"] >= 3
    assert quorum["bluer_half"] and quorum["redder_half"]


# ---------------------------------------------------------------------------
# Per-shot structure
# ---------------------------------------------------------------------------


def test_two_shot_groups_are_flagged_as_a_split_interval(tmp_path: Path) -> None:
    cubes = [
        _cube(tmp_path / "cubes" / f"{shot}.nc", shifts_px=shift, shot_number=str(shot))
        for shot, shift in ((100400, 2.0), (100401, 2.0), (100500, 8.0), (100501, 8.0))
    ]

    payload = audit_cubes(cubes)

    assert "consider splitting the interval" in payload["interval_warning"]
    assert "100400..100401" in payload["interval_warning"]
    assert "100500..100501" in payload["interval_warning"]
    assert payload["verdict"] in {
        "aligned",
        "shifted",
        "misaligned-beyond-repair",
        "insufficient-data",
    }
    groups = {item["shot_number"]: item["group"] for item in payload["per_shot"]}
    assert groups == {"100400": 1, "100401": 1, "100500": 2, "100501": 2}
    per_shot = {item["shot_number"]: item for item in payload["per_shot"]}
    assert per_shot["100400"]["median_shift_px"] == pytest.approx(2.0, abs=0.2)
    assert per_shot["100500"]["median_shift_px"] == pytest.approx(8.0, abs=0.3)
    assert all(item["lines"] == len(ORDERS) for item in payload["per_shot"])


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def test_directory_argument_audits_every_cube_below_it(tmp_path: Path) -> None:
    root = tmp_path / "cubes"
    for shot in (100600, 100601):
        _cube(root / f"{shot}.nc", shifts_px=0.0, shot_number=str(shot))
    _cube(root / "nested" / "100602.nc", shifts_px=0.0, shot_number="100602")

    payload = audit_cubes([root])

    assert [item["cube"] for item in payload["sampled_cubes"]] == [
        "100600.nc",
        "100601.nc",
        "100602.nc",
    ]


def test_date_selection_uses_the_documented_attribute(tmp_path: Path) -> None:
    root = tmp_path / "cubes"
    for shot, day in ((100700, "2026-08-10"), (100701, "2026-08-14")):
        _cube(
            root / f"{shot}.nc",
            shifts_px=0.0,
            shot_number=str(shot),
            t_start=f"{day}T09:00:00+09:00",
        )

    payload = audit_cubes([root], date_from="2026-08-12", date_to="2026-08-20")

    assert [item["shot_number"] for item in payload["sampled_cubes"]] == ["100701"]
    assert payload["sampled_cubes"][0]["date_attribute"] == "t_start"
    assert payload["sample_rule"]["date_from"] == "2026-08-12"

    with pytest.raises(DriftError, match="no sampled cube falls between"):
        audit_cubes([root], date_from="2020-01-01", date_to="2020-12-31")


def test_shot_selection_matches_whole_tokens_only(tmp_path: Path) -> None:
    paths = [tmp_path / name for name in ("shot_042.nc", "shot_142.nc", "shot_242.nc")]
    for path in paths:
        path.write_bytes(b"")

    selected = select_sample_paths(paths, every=3, shots={"42"})

    assert [path.name for path in selected] == ["shot_042.nc"]
    assert [path.name for path in select_sample_paths(paths, every=3, shots={"142"})] == [
        "shot_042.nc",
        "shot_142.nc",
    ]


# ---------------------------------------------------------------------------
# Beyond repair names the drives
# ---------------------------------------------------------------------------


def test_beyond_repair_names_the_drives_from_a_merged_catalog(tmp_path: Path) -> None:
    root = tmp_path / "cubes"
    cube = _cube(root / "100800.nc", shifts_px=[2.0, -2.0, 1.6, -1.9], shot_number="100800")
    merged = merge_catalogs(
        [build_drive_catalog(root, volume_label="NIFS-A")], tmp_path / "all-years.json"
    )

    named = audit_cubes([cube], catalog=merged)
    anonymous = audit_cubes([cube])

    assert named["verdict"] == "misaligned-beyond-repair"
    assert "NIFS-A" in named["data_requirement"]
    assert "--catalog" in anonymous["data_requirement"]
    assert "raw SIF data" in anonymous["data_requirement"]


# ---------------------------------------------------------------------------
# Plasma-presence frame selection
# ---------------------------------------------------------------------------


def test_only_plasma_bright_frames_reach_the_centroid(tmp_path: Path) -> None:
    # Eighteen dark frames carry a weak decoy line at zero shift; the two plasma
    # frames carry the real lines five pixels away.  A median over all frames
    # would report the decoy and call the epoch aligned.
    specs = [(0.0, 1.5)] * 18 + [(5.0, 400.0)] * 2

    payload = audit_cubes(
        [_cube(tmp_path / "cubes" / "100900.nc", frame_specs=specs, shot_number="100900")]
    )

    assert payload["sampled_cubes"][0]["frame_selection"]["selected_frames"] == 2
    assert payload["verdict"] == "shifted"
    assert payload["summary"]["median_shift_px"] == pytest.approx(5.0, abs=0.2)


def test_a_cube_without_plasma_frames_is_insufficient_and_counted(tmp_path: Path) -> None:
    payload = audit_cubes(
        [
            _cube(
                tmp_path / "cubes" / "101000.nc",
                frame_specs=[(0.0, 0.0)] * 20,
                shot_number="101000",
            )
        ]
    )

    assert payload["verdict"] == "insufficient-data"
    assert payload["summary"]["skipped_cubes"] == 1
    assert payload["skipped_cubes"][0]["reason"] == "no plasma-bright frames"
    assert any(item.get("reason") == "no plasma-bright frames" for item in payload["lines"])


# ---------------------------------------------------------------------------
# Thresholds derived from the calibration's own accuracy
# ---------------------------------------------------------------------------


def test_export_writes_wavelength_accuracy_from_the_alignment_rms(tmp_path: Path) -> None:
    pytest.importorskip("spectrocube")
    from echelle_spectra.tools.spectrocube_export import to_spectrocube

    snapshot = _snapshot(tmp_path, alignment={"dx_px": 1.0, "rms_px": 0.67})

    class FakeSpectrum:
        pass

    spectrum = FakeSpectrum()
    wavelength, detector_pixel, echelle_order = _axis(POLYNOMIALS)
    spectrum.wavelength = wavelength
    spectrum.detector_pixel = detector_pixel
    spectrum.echelle_order = echelle_order.astype(np.int64)
    spectrum.wavelength_polynomials = [
        {"order": order, "coefficients": POLYNOMIALS[order].tolist()} for order, _ in ORDERS
    ]
    spectrum.counts = np.ones((1, wavelength.size))
    spectrum.info = {"ExposureTime": 1.0}
    spectrum.fpth = Path("100100_Echelle.SIF")
    spectrum.shotnumber = "100100"

    cube = to_spectrocube(spectrum, snapshot=snapshot, units="counts")

    dispersions = [abs(_dispersion(centre)) for _order, centre in ORDERS]
    assert cube.ds.attrs["wavelength_accuracy_nm"] == pytest.approx(
        0.67 * float(np.median(dispersions)), rel=1e-3
    )
    assert cube.ds.attrs["wavelength_accuracy_source"] == "snapshot alignment rms_px"


def test_audit_relaxes_alignment_to_the_stated_accuracy(tmp_path: Path) -> None:
    strict = audit_cubes([_cube(tmp_path / "strict" / "101100.nc", shifts_px=0.8)])
    relaxed = audit_cubes(
        [_cube(tmp_path / "relaxed" / "101100.nc", shifts_px=0.8, accuracy_nm=0.012)]
    )

    assert strict["verdict"] == "shifted"
    assert relaxed["verdict"] == "aligned"
    tolerances = {item["alignment_tolerance_px"] for item in _measured(relaxed)}
    assert min(tolerances) > ALIGNMENT_MAX_PIXEL_RESIDUAL
