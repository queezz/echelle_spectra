"""Tests for echelle_spectra.tools.spectrocube_export.

All tests use a small synthetic Spectrum-like object so no real SIF files or
calibration resources are needed.
"""

from __future__ import annotations

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
