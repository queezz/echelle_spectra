"""Tests for Spectrum dark/background subtraction edge cases.

No real SIF files or calibration resources are needed – a small duck-typed
image object is used to drive Spectrum.__init__ directly.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from echelle_spectra.tools.echelle import Spectrum


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_fake_image(spectra: np.ndarray) -> object:
    """Return a minimal duck-typed EchelleImage-like object."""

    n_wl = spectra.shape[1]

    class _FakeClbr:
        absolute = {
            "wm": np.ones(n_wl),
            "wmsr": np.ones(n_wl),
            "phmsr": np.ones(n_wl),
        }
        direction = 1  # positive: no flip

    class _FakeImage:
        pass

    img = _FakeImage()
    img.spectra = spectra
    img.wavelength = np.linspace(400.0, 700.0, n_wl)
    img.fpth = "/fake/path.sif"
    img.info = {
        "NumberOfFrames": spectra.shape[0],
        "ExposureTime": 1.0,
        "CycleTime": 0.5,
        "DetectorTemperature": -70,
        "xbin": 1,
        "ybin": 1,
        "BackgroundFrames": [],
    }
    img.clbr = _FakeClbr()
    return img


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestNoDarkFrame:
    """Spectrum built from data where every frame has strong signal (no dark)."""

    def _make_all_bright(self, n_frames: int = 5, n_wl: int = 40) -> np.ndarray:
        """All frames have a prominent peak so none qualify as background."""
        rng = np.random.default_rng(1)
        base = rng.integers(1000, 5000, size=(n_frames, n_wl)).astype(float)
        # inject a tall spike into every frame so no frame passes the 5-sigma test
        base[:, n_wl // 2] += 1e6
        return base

    def test_counts_has_no_nan(self):
        spectra = self._make_all_bright()
        sp = Spectrum(_make_fake_image(spectra))
        assert np.all(np.isfinite(sp.counts)), "counts must be finite when no dark frame"

    def test_background_frames_is_empty(self):
        spectra = self._make_all_bright()
        sp = Spectrum(_make_fake_image(spectra))
        assert sp.info["BackgroundFrames"] == []

    def test_counts_equal_to_raw_spectra(self):
        """With no subtraction the counts must equal the original spectra."""
        spectra = self._make_all_bright()
        sp = Spectrum(_make_fake_image(spectra))
        np.testing.assert_array_equal(sp.counts, spectra)

    def test_no_runtime_warning(self):
        spectra = self._make_all_bright()
        img = _make_fake_image(spectra)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            # Should not raise – previously emitted "invalid value in divide"
            sp = Spectrum(img)
        assert sp is not None

    def test_calibrated_arrays_finite(self):
        spectra = self._make_all_bright()
        sp = Spectrum(_make_fake_image(spectra))
        for attr in ("wm", "wmsr", "phmsr"):
            arr = getattr(sp, attr)
            assert np.all(np.isfinite(arr)), f"{attr} must be finite when no dark frame"


class TestWithDarkFrame:
    """Sanity-check that normal dark subtraction still works correctly."""

    def _make_mixed(self, n_frames: int = 6, n_wl: int = 40) -> np.ndarray:
        """First 2 frames are dark (low, flat), rest have bright peaks."""
        rng = np.random.default_rng(2)
        data = np.zeros((n_frames, n_wl))
        # dark frames: low uniform values
        data[:2] = rng.integers(10, 50, size=(2, n_wl)).astype(float)
        # bright frames: dark level + a large spike
        data[2:] = rng.integers(100, 500, size=(n_frames - 2, n_wl)).astype(float)
        data[2:, n_wl // 2] += 1e5
        return data

    def test_background_frames_detected(self):
        spectra = self._make_mixed()
        sp = Spectrum(_make_fake_image(spectra))
        assert len(sp.info["BackgroundFrames"]) > 0

    def test_counts_finite(self):
        spectra = self._make_mixed()
        sp = Spectrum(_make_fake_image(spectra))
        assert np.all(np.isfinite(sp.counts))

    def test_subtraction_reduces_dark_frames(self):
        """After subtraction, dark frames should be close to zero."""
        spectra = self._make_mixed()
        sp = Spectrum(_make_fake_image(spectra))
        bg_idx = sp.info["BackgroundFrames"]
        # The mean of dark-subtracted background frames should be near zero
        assert np.abs(sp.counts[bg_idx].mean()) < 50


class TestSingleFrame:
    """Single-frame SIF files skip background detection entirely."""

    def test_single_frame_no_subtraction(self):
        rng = np.random.default_rng(3)
        spectra = rng.integers(100, 5000, size=(1, 30)).astype(float)
        sp = Spectrum(_make_fake_image(spectra))
        assert sp.info["BackgroundFrames"] == []
        np.testing.assert_array_equal(sp.counts, spectra)

    def test_single_frame_counts_finite(self):
        rng = np.random.default_rng(4)
        spectra = rng.integers(100, 5000, size=(1, 30)).astype(float)
        sp = Spectrum(_make_fake_image(spectra))
        assert np.all(np.isfinite(sp.counts))
