"""Focused Packets 9--13 implementation-candidate acceptance."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from echelle_spectra.campaign_run import RunReceipt, utc_now
from echelle_spectra.catalog import build_drive_catalog, load_catalog, merge_catalogs
from echelle_spectra.drift import (
    DriftError,
    audit_cubes,
    create_refinement_snapshot,
    require_sampled_verdict,
    verdict_from_evidence,
    write_drift_evidence,
)
from echelle_spectra.historical import bundled_historical_manifests
from echelle_spectra.lhd_text import TEXT_SCHEMA, write_cube_text
from echelle_spectra.reading_room import build_reading_room
from echelle_spectra.recalibration import RecalibrationError, recalibrate_cube, recalibrate_dataset
from echelle_spectra.snapshot import create_snapshot


def _snapshot(root: Path, snapshot_id: str, *, shift: float = 0.0, pattern: str = "same"):
    sources = root / f"sources-{snapshot_id}"
    sources.mkdir()
    files = {}
    for role, name in {
        "pattern": "legacy-pattern.txt",
        "sphere": "legacy-sphere.sif",
        "sphere_background": "legacy-sphere-background.sif",
        "integral": "legacy-integral.txt",
    }.items():
        path = sources / name
        path.write_text(pattern if role == "pattern" else "fixture", encoding="utf-8")
        files[role] = path
    wavelength = sources / "legacy-wavelength.txt"
    wavelength.write_text(
        "# fixture\n# order from to center wavelength\n"
        + "\n".join(
            f"0 {index} {index + 1} {index} {410.0 + index + shift} line"
            for index in range(5)
        )
        + "\n",
        encoding="utf-8",
    )
    files["wavelength"] = wavelength
    return create_snapshot(
        root / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=["ThAr"],
        validity={"shot_from": 1, "shot_to": 999999},
    )


def _packet8_dataset(snapshot, *, wavelengths: np.ndarray | None = None) -> xr.Dataset:
    wavelength = np.asarray(wavelengths if wavelengths is not None else np.arange(410.0, 415.0))
    factor = np.full(wavelength.size, 2.0)
    intensity = np.arange(1.0, wavelength.size + 1.0) * factor
    attrs = {
        "spectrocube_version": "0.2.0",
        "instrument_id": "echelle",
        "calibration_type": "absolute",
        "intensity_units": "W/m2/nm/sr",
        "wavelength_medium": "air",
        "shot_number": "12345",
        "created_at": "2025-01-01T00:00:00+00:00",
        "trigger_delay_s": 2.5,
        "frame_interval_s": 0.5,
        "exposure_s": 0.25,
        **snapshot.provenance_attrs(),
        "wavelength_polynomials_json": json.dumps(
            {
                "schema": "spectrocube.wavelength-polynomials/v1",
                "coefficient_order": "descending_power",
                "input": "detector_pixel",
                "input_units": "pixel",
                "output": "wavelength",
                "output_units": "nm",
                "orders": [{"order": 0, "coefficients": [1.0, 410.0]}],
            }
        ),
    }
    ds = xr.Dataset(
        {"intensity": (("frame", "wavelength"), intensity[None, :])},
        coords={
            "frame": [0],
            "wavelength": wavelength,
            "detector_pixel": ("wavelength", np.arange(wavelength.size, dtype=float)),
            "echelle_order": ("wavelength", np.zeros(wavelength.size, dtype=np.int64)),
        },
        attrs=attrs,
    )
    ds.coords["wavelength"].attrs.update(units="nm", medium="air")
    ds.coords["detector_pixel"].attrs.update(
        units="pixel", detector_axis="column", reference_frame="raw_detector", index_origin=0
    )
    ds["applied_absolute_calibration_factor"] = (
        ("wavelength",),
        factor,
        {
            "units": "W/m2/nm/sr per (counts/s)",
            "source_units": "counts/s",
            "application": "stored_intensity = source_signal * applied_absolute_calibration_factor",
        },
    )
    return ds


#: One Balmer line per Echelle order, with the instrument's dispersion scale.
_DRIFT_ORDERS = ((6, 656.2790), (12, 486.1350), (16, 434.0470), (18, 410.1734))
_DRIFT_WIDTH_PX = 1024


def _drift_dispersion(centre_nm: float) -> float:
    return 0.0108 * centre_nm / 656.2790


def _line_cube(path: Path, shift_px: float | None, snapshot_id: str = "20250101_cmos") -> None:
    """Write a cube whose lines sit where a shifted detector would put them."""

    pixels = np.arange(_DRIFT_WIDTH_PX, dtype=float)
    centre_px = (_DRIFT_WIDTH_PX - 1) / 2.0
    segments, signal, detector_pixel, echelle_order = [], [], [], []
    for order, centre in _DRIFT_ORDERS:
        dispersion = _drift_dispersion(centre)
        wavelength = centre + dispersion * (pixels - centre_px)
        values = 0.01 * np.sin(wavelength)
        if shift_px is not None:
            observed = centre + dispersion * shift_px
            values = values + 50.0 * np.exp(-0.5 * ((wavelength - observed) / 0.035) ** 2)
        segments.append(wavelength)
        signal.append(values)
        detector_pixel.append(pixels)
        echelle_order.append(np.full(pixels.size, order, dtype=np.int64))
    wavelength = np.concatenate(segments)
    permutation = np.argsort(wavelength, kind="stable")
    xr.Dataset(
        {"intensity": ("wavelength", np.concatenate(signal)[permutation])},
        coords={
            "wavelength": wavelength[permutation],
            "detector_pixel": ("wavelength", np.concatenate(detector_pixel)[permutation]),
            "echelle_order": ("wavelength", np.concatenate(echelle_order)[permutation]),
        },
        attrs={
            "spectrocube_version": "0.2.0",
            "instrument_id": "echelle",
            "calibration_type": "counts",
            "intensity_units": "counts",
            "wavelength_medium": "air",
            "snapshot_id": snapshot_id,
            "shot_number": path.stem,
            "trigger_delay_s": 2.5,
            "frame_interval_s": 0.5,
            "exposure_s": 0.25,
            "wavelength_polynomials_json": json.dumps(
                {
                    "schema": "spectrocube.wavelength-polynomials/v1",
                    "coefficient_order": "descending_power",
                    "input": "detector_pixel",
                    "input_units": "pixel",
                    "output": "wavelength",
                    "output_units": "nm",
                    "orders": [
                        {
                            "order": order,
                            "coefficients": [
                                _drift_dispersion(centre),
                                centre - _drift_dispersion(centre) * centre_px,
                            ],
                        }
                        for order, centre in _DRIFT_ORDERS
                    ],
                }
            ),
        },
    ).to_netcdf(path)


def _pixel_residual(order: int, centre_nm: float, shift_px: float) -> dict:
    """One resolved detector-space residual, shaped as the audit records it."""

    dispersion = _drift_dispersion(centre_nm)
    return {
        "cube": "fixture.nc",
        "shot_number": "100100",
        "status": "measured",
        "echelle_order": order,
        "expected_nm": centre_nm,
        "dispersion_nm_per_px": dispersion,
        "residual_nm": shift_px * dispersion,
        "pixel_residual_px": shift_px,
        "alignment_tolerance_px": 0.5,
    }


def test_catalog_cube_text_and_missing_drive_reading_room(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path, "20250101_cmos")
    drive_a = tmp_path / "drive-a"
    drive_b = tmp_path / "drive-b"
    drive_a.mkdir()
    drive_b.mkdir()
    cube = drive_a / "12345.nc"
    _packet8_dataset(snapshot).to_netcdf(cube)
    cube_b = drive_b / "12346.nc"
    _packet8_dataset(snapshot).to_netcdf(cube_b)
    catalog_a = build_drive_catalog(drive_a, volume_label="USB-A")
    catalog_b = build_drive_catalog(drive_b, volume_label="USB-B")
    merged = merge_catalogs([catalog_a, catalog_b], tmp_path / "all-years.json")
    catalog_b.unlink()

    text = write_cube_text(cube, tmp_path / "12345.txt")
    header = text.read_text(encoding="utf-8")
    # The header is frozen: provenance rides inside [Comments], never as fields.
    assert "# DimUnit = 'nm'" in header
    assert "# ShotNo = 12345" in header
    assert "# time = 2.5 + frameNo*0.5 (s)" in header
    assert f"# format_schema = {TEXT_SCHEMA}" in header
    assert f"# snapshot_id = {snapshot.snapshot_id}" in header
    assert "# snapshot_manifest_sha256 = " in header
    assert "[Provenance]" not in header

    page = build_reading_room(merged, tmp_path / "web", document_paths=[])
    rendered = page.read_text(encoding="utf-8")
    assert "never executes commands" in rendered
    assert "USB-B" in rendered
    data = json.loads((page.parent / "reading-room.json").read_text(encoding="utf-8"))
    missing = {source["volume_label"]: source["available"] for source in data["catalog"]["sources"]}
    assert missing == {"USB-A": True, "USB-B": False}


def test_snapshot_delta_recalibration_and_geometry_refusal(tmp_path: Path) -> None:
    old = _snapshot(tmp_path, "20240101_cmos")
    new = _snapshot(tmp_path, "20250101_cmos", shift=0.1)
    original = _packet8_dataset(old)
    revised, event = recalibrate_dataset(original, new, new_factor=np.full(5, 4.0))
    assert np.allclose(revised.wavelength, np.arange(410.1, 415.1))
    assert np.allclose(revised.intensity, original.intensity * 2.0)
    assert revised.attrs["snapshot_id"] == new.snapshot_id
    assert event["old_snapshot_id"] == old.snapshot_id
    assert json.loads(revised.attrs["recalibration_history_json"])[0]["changes"] == [
        "wavelength",
        "absolute-factor",
    ]

    source_cube = tmp_path / "old.nc"
    output_cube = tmp_path / "revised.nc"
    original.to_netcdf(source_cube)
    output, manifest = recalibrate_cube(
        source_cube,
        output_cube,
        new_snapshot_path=new.root,
        update_factor=False,
    )
    retroactive = json.loads(manifest.read_text(encoding="utf-8"))
    assert output.is_file()
    assert retroactive["old_snapshot_manifest_json"] == old.provenance_attrs()[
        "snapshot_manifest_json"
    ]
    assert retroactive["new_snapshot_manifest_json"] == new.provenance_attrs()[
        "snapshot_manifest_json"
    ]
    assert len(retroactive["input_sha256"]) == len(retroactive["output_sha256"]) == 64

    changed = _snapshot(tmp_path, "20260101_cmos", pattern="changed geometry")
    with pytest.raises(RecalibrationError, match="raw SIF data"):
        recalibrate_dataset(original, changed, update_factor=False)


def test_all_drift_verdicts_and_bulk_gate(tmp_path: Path) -> None:
    shifted_cube = tmp_path / "shot_42.nc"
    _line_cube(shifted_cube, 8.0)
    shifted = audit_cubes([shifted_cube])
    assert shifted["verdict"] == "shifted"
    assert "echelle drift refine" in shifted["repair_command"]
    coverage = (410.0, 660.0)
    aligned = [
        _pixel_residual(order, centre, shift)
        for (order, centre), shift in zip(_DRIFT_ORDERS, (0.0, 0.1, -0.1, 0.2))
    ]
    assert verdict_from_evidence(aligned, coverage_nm=coverage)[0] == "aligned"
    inconsistent = [
        _pixel_residual(order, centre, shift)
        for (order, centre), shift in zip(_DRIFT_ORDERS, (4.0, -4.0, 3.5, -3.2))
    ]
    assert (
        verdict_from_evidence(inconsistent, coverage_nm=coverage)[0]
        == "misaligned-beyond-repair"
    )
    assert verdict_from_evidence(aligned[:2], coverage_nm=coverage)[0] == "insufficient-data"
    assert verdict_from_evidence([{"status": "insufficient-data"}])[0] == "insufficient-data"

    evidence = write_drift_evidence(tmp_path / "shifted.json", shifted)
    with pytest.raises(DriftError, match="bulk processing refused"):
        require_sampled_verdict(evidence, {"20250101_cmos"})
    aligned = {**shifted, "verdict": "aligned"}
    aligned_path = write_drift_evidence(tmp_path / "aligned.json", aligned)
    assert require_sampled_verdict(aligned_path, {"20250101_cmos"})["verdict"] == "aligned"


def test_refinement_historical_and_connected_fixture_path(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path, "20250101_cmos")
    cube_root = tmp_path / "cubes"
    cube_root.mkdir()
    cube = cube_root / "shot_42.nc"
    _line_cube(cube, 8.0, snapshot.snapshot_id)
    evidence_payload = audit_cubes([cube], evidence_path=tmp_path / "drift.json")
    evidence = write_drift_evidence(tmp_path / "drift.json", evidence_payload)
    refined, accepted = create_refinement_snapshot(
        evidence,
        calibrations_root=tmp_path / "calibrations",
        accepted_shift_px=evidence_payload["summary"]["median_shift_px"],
    )
    assert refined.snapshot_id == "20250101_cmos-r1"
    accepted_payload = json.loads(accepted.read_text(encoding="utf-8"))
    assert accepted_payload["accepted"] is True
    assert require_sampled_verdict(accepted, {refined.snapshot_id})["accepted"] is True

    source = tmp_path / "shot_42.SIF"
    source.write_bytes(b"read-only-fixture-identity")
    receipt = RunReceipt.create(
        tmp_path / "runs" / "fixture-run",
        source_root=tmp_path,
        output_root=cube_root,
        pattern="*.SIF",
        volume_label="FIXTURE-USB",
        snapshot_id=snapshot.snapshot_id,
        expected_files=1,
    )
    identity = receipt.identity_for(source)
    receipt.append(
        identity,
        cube,
        status="exported",
        started_at=utc_now(),
        finished_at=utc_now(),
        duration_s=0.01,
        snapshot_id=snapshot.snapshot_id,
    )
    receipt.finish("completed")
    drive_catalog = build_drive_catalog(
        cube_root,
        volume_label=receipt.volume_label,
        receipt_dir=receipt.directory,
    )
    merged = merge_catalogs([drive_catalog], tmp_path / "all-years.json")
    write_cube_text(cube, tmp_path / "fixture.txt")
    page = build_reading_room(merged, tmp_path / "reading-room", drift_paths=[evidence])
    assert page.is_file()
    assert load_catalog(drive_catalog)["run"]["state"] == "completed"
    assert {item.snapshot_id for item in bundled_historical_manifests()} == {
        "20190529_cmos",
        "20240305_cmos",
        "20250926_cmos",
    }
