"""Tests for the echelle-pattern CLI."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from echelle_spectra.pattern_cli import _build_parser, _parse_csv_floats, _parse_csv_ints, main


class _FakeDetection:
    def __init__(self, count=3):
        self.n_peaks = count


class _FakeResult:
    def __init__(self):
        self.pattern = np.arange(18, dtype=int).reshape(6, 3)
        self.detections = [_FakeDetection(3), _FakeDetection(3)]


class _FakeTrial:
    success = True
    threshold = 0.12
    columns_px = np.array([1, 4])
    result = _FakeResult()

    @property
    def peak_counts(self):
        return np.array([3, 3])


def test_parse_csv_values():
    assert _parse_csv_ints("530, 680,750") == [530, 680, 750]
    assert _parse_csv_floats("0.10, 0.12") == [0.10, 0.12]


def test_parser_defaults():
    args = _build_parser().parse_args(["sphere.sif", "bg.sif"])
    assert args.expected_orders == 29
    assert args.output is None
    assert args.prior_pattern is None
    assert not args.no_prior_fit


def test_missing_input_exits(tmp_path):
    bg = tmp_path / "bg.sif"
    bg.touch()
    with pytest.raises(SystemExit) as exc:
        main([str(tmp_path / "missing.sif"), str(bg)])
    assert exc.value.code == 1


def test_preview_mode_does_not_write(tmp_path, capsys):
    sphere = tmp_path / "sphere.sif"
    bg = tmp_path / "bg.sif"
    sphere.touch()
    bg.touch()

    with patch("echelle_spectra.pattern_cli._load_image_pair", return_value=np.ones((5, 6))), \
        patch(
            "echelle_spectra.pattern_cli.trial_order_pattern_extraction",
            return_value=[_FakeTrial()],
        ), \
        pytest.raises(SystemExit) as exc:
        main([str(sphere), str(bg), "--expected-orders", "3"])

    assert exc.value.code == 0
    assert "Preview only" in capsys.readouterr().out


def test_output_writes_pattern(tmp_path):
    sphere = tmp_path / "sphere.sif"
    bg = tmp_path / "bg.sif"
    out = tmp_path / "pattern.txt"
    sphere.touch()
    bg.touch()

    with patch("echelle_spectra.pattern_cli._load_image_pair", return_value=np.ones((5, 6))), \
        patch(
            "echelle_spectra.pattern_cli.trial_order_pattern_extraction",
            return_value=[_FakeTrial()],
        ), \
        pytest.raises(SystemExit) as exc:
        main([str(sphere), str(bg), "--expected-orders", "3", "--output", str(out)])

    assert exc.value.code == 0
    loaded = np.loadtxt(out, dtype=int)
    np.testing.assert_array_equal(loaded, _FakeTrial.result.pattern)


def test_prior_pattern_runs_prior_guided_fit(tmp_path):
    sphere = tmp_path / "sphere.sif"
    bg = tmp_path / "bg.sif"
    prior = tmp_path / "prior.txt"
    sphere.touch()
    bg.touch()
    np.savetxt(prior, _FakeTrial.result.pattern, fmt="%d")

    with patch("echelle_spectra.pattern_cli._load_image_pair", return_value=np.ones((5, 6))), \
        patch(
            "echelle_spectra.pattern_cli.trial_order_pattern_extraction",
            return_value=[_FakeTrial()],
        ), \
        patch(
            "echelle_spectra.pattern_cli.extract_order_pattern_near_prior",
            return_value=_FakeResult(),
        ) as prior_fit, \
        pytest.raises(SystemExit) as exc:
        main([str(sphere), str(bg), "--expected-orders", "3", "--prior-pattern", str(prior)])

    assert exc.value.code == 0
    prior_fit.assert_called_once()
