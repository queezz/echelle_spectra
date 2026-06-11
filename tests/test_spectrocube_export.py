"""Tests for echelle_spectra.tools.spectrocube_export.

All tests use a small synthetic Spectrum-like object so no real SIF files or
calibration resources are needed.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

pytest.importorskip("spectrocube", reason="spectrocube package not installed")

from echelle_spectra.tools.spectrocube_export import export_spectrocube, to_spectrocube

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_synthetic_spectrum(
    n_frames: int = 3,
    n_wavelengths: int = 50,
    wl_start: float = 400.0,
    wl_end: float = 700.0,
    fpth: str = "/data/shot_042.sif",
) -> object:
    """Return a minimal duck-typed object that matches the Spectrum interface."""

    class _FakeSpectrum:
        pass

    sp = _FakeSpectrum()
    sp.wavelength = np.linspace(wl_start, wl_end, n_wavelengths)

    rng = np.random.default_rng(0)
    sp.counts = rng.integers(100, 10_000, size=(n_frames, n_wavelengths)).astype(float)

    # Physically calibrated versions (crude scaling so units make sense)
    scale_wm = 1e-6
    scale_wmsr = scale_wm / (4 * np.pi)
    scale_phmsr = scale_wmsr * sp.wavelength[None, :] * 1e-9 / (3e8 * 6.626e-34)

    sp.wm = sp.counts * scale_wm
    sp.wmsr = sp.counts * scale_wmsr
    sp.phmsr = sp.counts * scale_phmsr

    sp.info = {
        "NumberOfFrames": n_frames,
        "ExposureTime": 0.5,
        "CycleTime": 1.0,
        "BackgroundFrames": [0],
    }
    sp.fpth = fpth
    sp.shotnumber = "042"
    sp.absolute = {}  # present but empty for this synthetic case
    sp.calibration_folder = "/cal"
    sp.calibration_files = {
        "orders": "pattern_CMOS_20250926.txt",
        "wavelength": "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
    }
    sp.calibration_order_count = 2
    sp.calibration_detector_width_px = 2304
    sp.calibration_order_half_width_px = 8
    sp.order_border_pixel_ranges = [
        {"order": 0, "start_px": 0, "end_px": 1000, "n_px": 1001},
        {"order": 1, "start_px": 1001, "end_px": 2303, "n_px": 1303},
    ]
    sp.order_wavelength_ranges_nm = [
        {"order": 0, "min_nm": 600.0, "max_nm": 650.0, "n_px": 2304},
        {"order": 1, "min_nm": 550.0, "max_nm": 610.0, "n_px": 2304},
    ]

    return sp


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestToSpectrocube:
    def test_returns_spectrocube_type(self):
        from spectrocube import SpectroCube

        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert isinstance(sc, SpectroCube)

    def test_wavelength_preserved(self):
        sp = _make_synthetic_spectrum(n_wavelengths=80)
        sc = to_spectrocube(sp)
        np.testing.assert_allclose(sc.wavelength, sp.wavelength)

    def test_single_frame_squeezed_to_1d(self):
        sp = _make_synthetic_spectrum(n_frames=1)
        sc = to_spectrocube(sp, squeeze_single_frame=True)
        assert sc.intensity.ndim == 1
        assert sc.dims == ("wavelength",)

    def test_single_frame_not_squeezed_when_disabled(self):
        sp = _make_synthetic_spectrum(n_frames=1)
        sc = to_spectrocube(sp, squeeze_single_frame=False)
        assert sc.intensity.ndim == 2
        assert sc.dims == ("frame", "wavelength")

    def test_multiframe_is_2d(self):
        sp = _make_synthetic_spectrum(n_frames=5)
        sc = to_spectrocube(sp)
        assert sc.intensity.ndim == 2
        assert sc.dims == ("frame", "wavelength")
        assert sc.intensity.shape == (5, 50)

    def test_intensity_counts_preserved(self):
        sp = _make_synthetic_spectrum(n_frames=2)
        sc = to_spectrocube(sp, units="counts", squeeze_single_frame=False)
        np.testing.assert_allclose(sc.intensity, sp.counts)

    def test_intensity_wm_preserved(self):
        sp = _make_synthetic_spectrum(n_frames=2)
        sc = to_spectrocube(sp, units="wm", squeeze_single_frame=False)
        np.testing.assert_allclose(sc.intensity, sp.wm)

    def test_intensity_wmsr_preserved(self):
        sp = _make_synthetic_spectrum(n_frames=2)
        sc = to_spectrocube(sp, units="wmsr", squeeze_single_frame=False)
        np.testing.assert_allclose(sc.intensity, sp.wmsr)

    def test_intensity_phmsr_preserved(self):
        sp = _make_synthetic_spectrum(n_frames=2)
        sc = to_spectrocube(sp, units="phmsr", squeeze_single_frame=False)
        np.testing.assert_allclose(sc.intensity, sp.phmsr)

    def test_calibration_type_counts(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, units="counts")
        assert sc.calibration_type == "counts"

    def test_calibration_type_absolute_for_wm(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, units="wm")
        assert sc.calibration_type == "absolute"

    def test_instrument_id_default(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["instrument_id"] == "echelle"

    def test_instrument_id_custom(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, instrument_id="BlackEchelle")
        assert sc.ds.attrs["instrument_id"] == "BlackEchelle"

    def test_metadata_exposure_preserved(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["exposure_s"] == pytest.approx(0.5)

    def test_metadata_cycle_time_preserved(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["frame_interval_s"] == pytest.approx(1.0)

    def test_metadata_source_package(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["source_package"] == "echelle_spectra"

    def test_metadata_source_file(self):
        sp = _make_synthetic_spectrum(fpth="/data/shot_042.sif")
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["source_file"] == "/data/shot_042.sif"

    def test_metadata_shot_number(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["shot_number"] == "042"

    def test_extra_attrs_forwarded(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, notes="test run", grating="316 l/mm")
        assert sc.ds.attrs["notes"] == "test run"
        assert sc.ds.attrs["grating"] == "316 l/mm"

    def test_calibration_file_metadata_preserved(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["calibration_folder"] == "/cal"
        assert sc.ds.attrs["calibration_order_pattern_file"] == "pattern_CMOS_20250926.txt"
        assert (
            sc.ds.attrs["wavelength_calibration_file"]
            == "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
        )

    def test_order_border_metadata_preserved(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp)
        assert sc.ds.attrs["calibration_order_count"] == 2
        assert sc.ds.attrs["calibration_detector_width_px"] == 2304
        assert sc.ds.attrs["calibration_order_half_width_px"] == 8
        assert json.loads(sc.ds.attrs["order_border_pixel_ranges_json"]) == [
            {"end_px": 1000, "n_px": 1001, "order": 0, "start_px": 0},
            {"end_px": 2303, "n_px": 1303, "order": 1, "start_px": 1001},
        ]
        assert json.loads(sc.ds.attrs["order_wavelength_ranges_nm_json"]) == [
            {"max_nm": 650.0, "min_nm": 600.0, "n_px": 2304, "order": 0},
            {"max_nm": 610.0, "min_nm": 550.0, "n_px": 2304, "order": 1},
        ]

    def test_invalid_units_raises(self):
        sp = _make_synthetic_spectrum()
        with pytest.raises(ValueError, match="Unknown units"):
            to_spectrocube(sp, units="bad_unit")

    def test_validation_passes_for_counts(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, units="counts")
        report = sc.validate()
        assert report.ok, str(report)

    def test_validation_passes_for_absolute(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, units="wm")
        report = sc.validate()
        assert report.ok, str(report)

    def test_nonfinite_intensity_column_is_dropped(self):
        sp = _make_synthetic_spectrum(n_frames=2, n_wavelengths=5)
        sp.wmsr[:, 2] = np.inf
        sc = to_spectrocube(sp, units="wmsr", squeeze_single_frame=False)
        assert sc.intensity.shape == (2, 4)
        assert np.all(np.isfinite(sc.intensity))
        assert sc.ds.attrs["dropped_nonfinite_wavelength_columns"] == 1
        assert 550.0 not in sc.wavelength

    def test_wavelength_min_crop_updates_metadata(self):
        sp = _make_synthetic_spectrum(n_frames=2, n_wavelengths=5)
        sc = to_spectrocube(
            sp,
            units="counts",
            squeeze_single_frame=False,
            wavelength_min_nm=475.0,
        )
        np.testing.assert_allclose(sc.wavelength, [475.0, 550.0, 625.0, 700.0])
        assert sc.intensity.shape == (2, 4)
        assert sc.ds.attrs["wavelength_crop_min_nm"] == pytest.approx(475.0)
        assert sc.ds.attrs["original_wavelength_min_nm"] == pytest.approx(400.0)
        assert sc.ds.attrs["original_wavelength_max_nm"] == pytest.approx(700.0)
        assert sc.ds.attrs["original_wavelength_points"] == 5
        assert sc.ds.attrs["dropped_wavelength_crop_columns"] == 1

    def test_wavelength_crop_rejects_empty_result(self):
        sp = _make_synthetic_spectrum(n_frames=2, n_wavelengths=5)
        with pytest.raises(ValueError, match="removed all columns"):
            to_spectrocube(sp, wavelength_min_nm=900.0)

    def test_wavelength_medium_stored(self):
        sp = _make_synthetic_spectrum()
        sc = to_spectrocube(sp, wavelength_medium="vacuum")
        assert sc.ds.attrs["wavelength_medium"] == "vacuum"


class TestExportSpectrocube:
    def test_saves_and_reloads(self, tmp_path):
        from spectrocube import SpectroCube

        sp = _make_synthetic_spectrum(n_frames=2)
        out = tmp_path / "spectrum.nc"
        sc = export_spectrocube(sp, str(out), units="counts", squeeze_single_frame=False)

        assert out.exists()
        sc2 = SpectroCube.load(str(out))
        np.testing.assert_allclose(sc2.wavelength, sc.wavelength)
        np.testing.assert_allclose(sc2.intensity, sc.intensity)

    def test_roundtrip_metadata(self, tmp_path):
        sp = _make_synthetic_spectrum()
        out = tmp_path / "spectrum.nc"
        export_spectrocube(sp, str(out))

        from spectrocube import SpectroCube

        sc2 = SpectroCube.load(str(out))
        assert sc2.ds.attrs["source_package"] == "echelle_spectra"
        assert sc2.ds.attrs["shot_number"] == "042"
        assert sc2.ds.attrs["exposure_s"] == pytest.approx(0.5)

    def test_roundtrip_absolute_units(self, tmp_path):
        from spectrocube import SpectroCube

        sp = _make_synthetic_spectrum(n_frames=3)
        out = tmp_path / "spectrum_wm.nc"
        export_spectrocube(sp, str(out), units="wm", squeeze_single_frame=False)

        sc2 = SpectroCube.load(str(out))
        assert sc2.calibration_type == "absolute"
        np.testing.assert_allclose(sc2.intensity, sp.wm, rtol=1e-6)

    def test_returns_spectrocube_object(self, tmp_path):
        from spectrocube import SpectroCube

        sp = _make_synthetic_spectrum()
        out = tmp_path / "spectrum.nc"
        result = export_spectrocube(sp, str(out))
        assert isinstance(result, SpectroCube)


class TestWavelengthAxisNormalization:
    """The exporter must handle descending wavelength axes and reject non-monotonic ones."""

    def _make_descending_spectrum(self, n_frames: int = 2, n_wavelengths: int = 5):
        """Spectrum whose wavelength runs high-to-low (descending), as produced
        by some echelle calibrations before the Spectrum flip is applied."""

        class _FakeSpectrum:
            pass

        sp = _FakeSpectrum()
        # descending: 700 → 400 nm
        sp.wavelength = np.linspace(700.0, 400.0, n_wavelengths)
        rng = np.random.default_rng(7)
        sp.counts = rng.integers(1, 100, size=(n_frames, n_wavelengths)).astype(float)
        sp.wm = sp.counts * 1e-6
        sp.wmsr = sp.counts * 1e-7
        sp.phmsr = sp.counts * 1e-12
        sp.info = {"NumberOfFrames": n_frames, "ExposureTime": 0.5, "CycleTime": 1.0,
                   "BackgroundFrames": []}
        sp.fpth = "/data/shot_desc.sif"
        sp.shotnumber = None
        return sp

    # -- descending wavelength ------------------------------------------------

    def test_descending_wavelength_exports_without_error(self):
        sp = self._make_descending_spectrum()
        sc = to_spectrocube(sp)  # must not raise
        report = sc.validate()
        assert report.ok, str(report)

    def test_descending_wavelength_is_reversed_to_ascending(self):
        sp = self._make_descending_spectrum(n_wavelengths=5)
        sc = to_spectrocube(sp, squeeze_single_frame=False)
        wl = sc.wavelength
        assert np.all(np.diff(wl) > 0), f"Expected ascending, got {wl}"

    def test_descending_wavelength_values_preserved(self):
        """After reversal, the wavelength values must be the same set."""
        sp = self._make_descending_spectrum(n_wavelengths=5)
        sc = to_spectrocube(sp, squeeze_single_frame=False)
        np.testing.assert_array_equal(
            np.sort(sc.wavelength), np.sort(sp.wavelength)
        )

    def test_descending_1d_intensity_reversed_with_wavelength(self):
        """Single frame (1D): reversed wavelength → intensity also reversed."""

        class _Sp:
            pass

        sp = _Sp()
        sp.wavelength = np.array([500.0, 400.0, 300.0])
        sp.counts = np.array([[1.0, 2.0, 3.0]])  # 1 frame × 3 wavelengths
        sp.wm = sp.counts * 1e-6
        sp.wmsr = sp.counts * 1e-7
        sp.phmsr = sp.counts * 1e-12
        sp.info = {"ExposureTime": 1.0, "CycleTime": 1.0, "BackgroundFrames": []}
        sp.fpth = None
        sp.shotnumber = None

        sc = to_spectrocube(sp, units="counts", squeeze_single_frame=True)
        np.testing.assert_array_equal(sc.wavelength, [300.0, 400.0, 500.0])
        np.testing.assert_array_equal(sc.intensity, [3.0, 2.0, 1.0])

    def test_descending_2d_intensity_reversed_along_wavelength_axis(self):
        """Multi-frame (2D): each frame's intensity is reversed along wavelength."""

        class _Sp:
            pass

        sp = _Sp()
        sp.wavelength = np.array([500.0, 400.0, 300.0])
        sp.counts = np.array([
            [1.0, 2.0, 3.0],   # frame 0
            [4.0, 5.0, 6.0],   # frame 1
        ])
        sp.wm = sp.counts * 1e-6
        sp.wmsr = sp.counts * 1e-7
        sp.phmsr = sp.counts * 1e-12
        sp.info = {"ExposureTime": 1.0, "CycleTime": 1.0, "BackgroundFrames": []}
        sp.fpth = None
        sp.shotnumber = None

        sc = to_spectrocube(sp, units="counts", squeeze_single_frame=False)
        np.testing.assert_array_equal(sc.wavelength, [300.0, 400.0, 500.0])
        # Each row must be reversed
        np.testing.assert_array_equal(sc.intensity[0], [3.0, 2.0, 1.0])
        np.testing.assert_array_equal(sc.intensity[1], [6.0, 5.0, 4.0])

    def test_ascending_wavelength_intensity_unchanged(self):
        """Ascending input must pass through unmodified."""

        class _Sp:
            pass

        sp = _Sp()
        sp.wavelength = np.array([300.0, 400.0, 500.0])
        sp.counts = np.array([[1.0, 2.0, 3.0]])
        sp.wm = sp.counts * 1e-6
        sp.wmsr = sp.counts * 1e-7
        sp.phmsr = sp.counts * 1e-12
        sp.info = {"ExposureTime": 1.0, "CycleTime": 1.0, "BackgroundFrames": []}
        sp.fpth = None
        sp.shotnumber = None

        sc = to_spectrocube(sp, units="counts", squeeze_single_frame=True)
        np.testing.assert_array_equal(sc.wavelength, [300.0, 400.0, 500.0])
        np.testing.assert_array_equal(sc.intensity, [1.0, 2.0, 3.0])

    # -- echelle seam stitching (mostly increasing with tiny reversals) -------

    def test_small_seam_reversal_is_cosorted(self):
        """Realistic echelle case: globally increasing with a single tiny seam reversal."""

        class _Sp:
            pass

        sp = _Sp()
        # ascending 400 → 800 nm with one out-of-order point at index 50
        wl = np.linspace(400.0, 800.0, 200)
        wl[50] = wl[51] - 0.001  # tiny dip; still inside (wl[49], wl[51])
        wl[51], wl[50] = wl[50], wl[51]  # swap so dw[50] < 0
        sp.wavelength = wl
        sp.counts = np.arange(200, dtype=float)[None, :]  # 1 frame
        sp.wm = sp.counts * 1e-6
        sp.wmsr = sp.counts * 1e-7
        sp.phmsr = sp.counts * 1e-12
        sp.info = {"ExposureTime": 1.0, "CycleTime": 1.0, "BackgroundFrames": []}
        sp.fpth = None
        sp.shotnumber = None

        sc = to_spectrocube(sp, units="counts", squeeze_single_frame=True)
        assert np.all(np.diff(sc.wavelength) > 0)
        # validation must pass
        assert sc.validate().ok, str(sc.validate())

    def test_seam_cosort_preserves_pairing(self):
        """After co-sort, each wavelength must still be paired with its intensity."""
        from echelle_spectra.tools.spectrocube_export import _normalize_wavelength_axis

        # Build paired data where intensity == wavelength * 1000 so we can check
        # the pairing was preserved across the sort.
        wl = np.array([400.0, 401.0, 402.0, 403.5, 403.0, 404.0, 405.0])
        intensity = (wl * 1000.0)[None, :]  # shape (1, 7)

        wl2, intensity2 = _normalize_wavelength_axis(wl.copy(), intensity.copy())

        assert np.all(np.diff(wl2) > 0)
        np.testing.assert_allclose(intensity2[0], wl2 * 1000.0)

    def test_seam_cosort_2d_preserves_pairing(self):
        """Same as above for a multi-frame array."""
        from echelle_spectra.tools.spectrocube_export import _normalize_wavelength_axis

        wl = np.array([400.0, 401.0, 402.0, 403.5, 403.0, 404.0])
        intensity = np.array([
            wl * 10.0,
            wl * 100.0,
        ])  # shape (2, 6)

        wl2, intensity2 = _normalize_wavelength_axis(wl.copy(), intensity.copy())

        assert np.all(np.diff(wl2) > 0)
        np.testing.assert_allclose(intensity2[0], wl2 * 10.0)
        np.testing.assert_allclose(intensity2[1], wl2 * 100.0)

    def test_duplicate_wavelength_after_cosort_raises(self):
        from echelle_spectra.tools.spectrocube_export import _normalize_wavelength_axis

        wl = np.array([400.0, 401.0, 402.0, 401.5, 401.5, 403.0])
        intensity = np.ones((1, 6))
        with pytest.raises(ValueError, match="duplicate"):
            _normalize_wavelength_axis(wl, intensity)

    # -- non-monotonic raises -------------------------------------------------

    def test_non_monotonic_wavelength_raises_value_error(self):
        """Genuinely scrambled axis (random permutation) must raise."""
        sp = _make_synthetic_spectrum(n_wavelengths=500)
        rng = np.random.default_rng(0)
        sp.wavelength = rng.permutation(sp.wavelength.copy())
        with pytest.raises(ValueError, match="reversals"):
            to_spectrocube(sp)

    def test_non_monotonic_error_message_is_clear(self):
        sp = _make_synthetic_spectrum(n_wavelengths=500)
        rng = np.random.default_rng(1)
        sp.wavelength = rng.permutation(sp.wavelength.copy())
        with pytest.raises(ValueError, match="reversals"):
            to_spectrocube(sp)

    # -- round-trip with descending input -------------------------------------

    def test_descending_export_roundtrip(self, tmp_path):
        from spectrocube import SpectroCube

        sp = self._make_descending_spectrum(n_frames=2, n_wavelengths=20)
        out = tmp_path / "desc.nc"
        export_spectrocube(sp, str(out), squeeze_single_frame=False)

        sc2 = SpectroCube.load(str(out))
        assert np.all(np.diff(sc2.wavelength) > 0), "Reloaded wavelength must be ascending"
        assert sc2.intensity.shape == (2, 20)


class TestEdgeCases:
    def test_spectrum_without_shotnumber(self):
        sp = _make_synthetic_spectrum()
        sp.shotnumber = None
        sc = to_spectrocube(sp)
        assert "shot_number" not in sc.ds.attrs

    def test_spectrum_without_fpth(self):
        sp = _make_synthetic_spectrum()
        sp.fpth = None
        sc = to_spectrocube(sp)
        assert "source_file" not in sc.ds.attrs

    def test_spectrum_without_background_frames(self):
        sp = _make_synthetic_spectrum()
        sp.info["BackgroundFrames"] = []
        sc = to_spectrocube(sp)
        assert "background_frames" not in sc.ds.attrs

    def test_wavelength_dims_correct_size(self):
        sp = _make_synthetic_spectrum(n_wavelengths=120)
        sc = to_spectrocube(sp, squeeze_single_frame=False)
        assert sc.sizes["wavelength"] == 120

    def test_frame_coord_correct_size(self):
        sp = _make_synthetic_spectrum(n_frames=7)
        sc = to_spectrocube(sp)
        assert sc.sizes["frame"] == 7
