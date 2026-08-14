from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from spectrocube import APPLIED_FACTOR_APPLICATION, SpectroCube

from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot
from echelle_spectra.tools.spectrocube_export import export_spectrocube, to_spectrocube


def _snapshot(tmp_path: Path):
    sources = tmp_path / "sources"
    sources.mkdir()
    files = {}
    for role in ROLE_FILENAMES:
        path = sources / f"{role}.input"
        path.write_text(f"fixture {role}\n", encoding="utf-8")
        files[role] = path
    return create_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20260814_cmos",
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity={"shot_from": 200000, "shot_to": 299999},
        alignment={"dx_px": 0.25, "dy_px": -0.5, "rms_px": 0.1},
        qc={"accepted": True, "line_count": 8},
    )


def _spectrum(*, decreasing_dispersion: bool) -> object:
    class FakeSpectrum:
        pass

    spectrum = FakeSpectrum()
    if decreasing_dispersion:
        # This is the state after Spectrum's whole-axis negative-dispersion
        # reversal. Raw detector columns remain their original numeric values.
        spectrum.wavelength = np.array([500.0, 500.5, 501.0, 400.0, 400.5, 401.0])
        spectrum.detector_pixel = np.array([2, 1, 0, 2, 1, 0], dtype=float)
        spectrum.echelle_order = np.array([1, 1, 1, 0, 0, 0], dtype=np.int16)
        spectrum.wavelength_polynomials = (
            {"order": 0, "coefficients": [-0.5, 401.0]},
            {"order": 1, "coefficients": [-0.5, 501.0]},
        )
    else:
        spectrum.wavelength = np.array([400.0, 400.5, 401.0, 500.0, 500.5, 501.0])
        spectrum.detector_pixel = np.array([0, 1, 2, 0, 1, 2], dtype=float)
        spectrum.echelle_order = np.array([0, 0, 0, 1, 1, 1], dtype=np.int16)
        spectrum.wavelength_polynomials = (
            {"order": 0, "coefficients": [0.5, 400.0]},
            {"order": 1, "coefficients": [0.5, 500.0]},
        )

    spectrum.counts = np.vstack(
        [
            100.0 + spectrum.echelle_order * 10.0 + spectrum.detector_pixel,
            200.0 + spectrum.echelle_order * 10.0 + spectrum.detector_pixel,
        ]
    )
    factor = 0.01 + spectrum.wavelength * 1e-5
    spectrum.absolute = {
        "wm": factor * 4.0 * np.pi,
        "wmsr": factor,
        "phmsr": factor * spectrum.wavelength * 1e-9 / (299792458.0 * 6.62607015e-34),
    }
    spectrum.info = {
        "NumberOfFrames": 2,
        "ExposureTime": 2.0,
        "CycleTime": 3.0,
        "BackgroundFrames": [],
    }
    spectrum.wm = spectrum.counts / 2.0 * spectrum.absolute["wm"]
    spectrum.wmsr = spectrum.counts / 2.0 * spectrum.absolute["wmsr"]
    spectrum.phmsr = spectrum.counts / 2.0 * spectrum.absolute["phmsr"]
    spectrum.fpth = Path("193778_Echelle.SIF")
    spectrum.shotnumber = "193778"
    spectrum.calibration_order_count = 2
    spectrum.calibration_detector_width_px = 3
    spectrum.calibration_order_half_width_px = 1
    spectrum.order_border_pixel_ranges = [
        {"order": 0, "start_px": 0, "end_px": 2, "n_px": 3},
        {"order": 1, "start_px": 0, "end_px": 2, "n_px": 3},
    ]
    spectrum.order_wavelength_ranges_nm = [
        {"order": 0, "min_nm": 400.0, "max_nm": 401.0, "n_px": 3},
        {"order": 1, "min_nm": 500.0, "max_nm": 501.0, "n_px": 3},
    ]
    return spectrum


def test_decreasing_dispersion_raw_pixels_share_crop_mask_and_sort(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path)
    spectrum = _spectrum(decreasing_dispersion=True)
    cube = to_spectrocube(
        spectrum,
        snapshot=snapshot,
        units="wmsr",
        squeeze_single_frame=False,
        wavelength_min_nm=400.5,
        wavelength_max_nm=500.5,
    )

    np.testing.assert_array_equal(cube.wavelength, [400.5, 401.0, 500.0, 500.5])
    np.testing.assert_array_equal(cube.ds["detector_pixel"], [1.0, 0.0, 2.0, 1.0])
    np.testing.assert_array_equal(cube.ds["echelle_order"], [0, 0, 1, 1])
    assert np.all(np.diff(cube.wavelength) > 0)

    payload = json.loads(cube.ds.attrs["wavelength_polynomials_json"])
    assert {record["order"] for record in payload["orders"]} == set(
        np.unique(cube.ds["echelle_order"]).tolist()
    )
    for record in payload["orders"]:
        selected = cube.ds["echelle_order"].values == record["order"]
        reconstructed = np.polyval(
            record["coefficients"], cube.ds["detector_pixel"].values[selected]
        )
        np.testing.assert_allclose(
            reconstructed,
            cube.wavelength[selected],
            rtol=0.0,
            atol=5e-10,
        )

    factor = cube.ds["applied_absolute_calibration_factor"]
    assert factor.attrs == {
        "units": "W/m2/nm/sr per (counts/s)",
        "source_units": "counts/s",
        "absolute_kind": "wmsr",
        "application": APPLIED_FACTOR_APPLICATION,
    }
    assert np.all(np.isfinite(factor)) and np.all(factor.values > 0)
    source_count_rate = cube.intensity / factor.values
    expected_counts = np.vstack(
        [
            [101.0, 100.0, 112.0, 111.0],
            [201.0, 200.0, 212.0, 211.0],
        ]
    )
    np.testing.assert_allclose(source_count_rate, expected_counts / 2.0, rtol=2e-12)
    assert cube.validate().ok, str(cube.validate())


def test_increasing_dispersion_finite_mask_is_identical_for_every_field(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path)
    spectrum = _spectrum(decreasing_dispersion=False)
    spectrum.wmsr[:, 1] = np.nan
    spectrum.absolute["wmsr"][1] = np.nan
    cube = to_spectrocube(
        spectrum,
        snapshot=snapshot,
        units="wmsr",
        squeeze_single_frame=False,
    )
    np.testing.assert_array_equal(cube.wavelength, [400.0, 401.0, 500.0, 500.5, 501.0])
    np.testing.assert_array_equal(cube.ds["detector_pixel"], [0, 2, 0, 1, 2])
    np.testing.assert_array_equal(cube.ds["echelle_order"], [0, 0, 1, 1, 1])
    assert cube.intensity.shape[-1] == cube.ds["applied_absolute_calibration_factor"].size
    assert cube.ds.attrs["dropped_nonfinite_wavelength_columns"] == 1
    assert cube.validate().ok, str(cube.validate())


def test_snapshot_and_all_calibration_digests_survive_netcdf(tmp_path: Path) -> None:
    snapshot = _snapshot(tmp_path)
    spectrum = _spectrum(decreasing_dispersion=True)
    output = tmp_path / "complete.nc"
    export_spectrocube(
        spectrum,
        str(output),
        snapshot=snapshot,
        units="wmsr",
        squeeze_single_frame=False,
    )
    loaded = SpectroCube.load(str(output))
    assert loaded.validate().ok, str(loaded.validate())
    assert loaded.ds.attrs["snapshot_id"] == snapshot.snapshot_id
    assert len(loaded.ds.attrs["snapshot_manifest_sha256"]) == 64
    manifest = json.loads(loaded.ds.attrs["snapshot_manifest_json"])
    assert manifest["id"] == snapshot.snapshot_id
    assert manifest["alignment"]["rms_px"] == 0.1
    digests = json.loads(loaded.ds.attrs["calibration_file_digests_json"])
    assert {record["role"] for record in digests} == set(ROLE_FILENAMES)
    assert all(len(record["sha256"]) == 64 for record in digests)
    assert "calibration_folder" not in loaded.ds.attrs


def test_snapshot_backed_counts_export_keeps_independent_optional_fields(tmp_path: Path) -> None:
    cube = to_spectrocube(
        _spectrum(decreasing_dispersion=False),
        snapshot=_snapshot(tmp_path),
        units="counts",
    )
    assert "detector_pixel" in cube.ds.coords
    assert "echelle_order" in cube.ds.coords
    assert "wavelength_polynomials_json" in cube.ds.attrs
    assert "applied_absolute_calibration_factor" not in cube.ds
    assert cube.validate().ok, str(cube.validate())
