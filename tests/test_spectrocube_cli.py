"""Tests for the echelle-spectrocube CLI.

All tests use mocks / monkeypatching so no real SIF files or calibration
resources are required.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from echelle_spectra.spectrocube_cli import _build_parser, _output_path_for, main

# ---------------------------------------------------------------------------
# Parser tests — pure argparse, no side effects
# ---------------------------------------------------------------------------


class TestParser:
    def test_single_file_defaults(self):
        args = _build_parser().parse_args(["shot.sif"])
        assert args.input == "shot.sif"
        assert args.units == "counts"
        assert args.camera == "CMOS"
        assert args.output is None
        assert args.overwrite is False
        assert args.dry_run is False
        assert args.verbose is False
        assert args.instrument_id == "echelle"
        assert args.pattern == "*.SIF"

    @pytest.mark.parametrize("units", ["counts", "wm", "wmsr", "phmsr"])
    def test_units_choices(self, units):
        args = _build_parser().parse_args(["shot.sif", "--units", units])
        assert args.units == units

    def test_invalid_units_exits(self):
        with pytest.raises(SystemExit):
            _build_parser().parse_args(["shot.sif", "--units", "bad"])

    def test_output_short_flag(self):
        args = _build_parser().parse_args(["shot.sif", "-o", "out.nc"])
        assert args.output == "out.nc"

    def test_output_long_flag(self):
        args = _build_parser().parse_args(["shot.sif", "--output", "out.nc"])
        assert args.output == "out.nc"

    def test_camera_ccd(self):
        args = _build_parser().parse_args(["shot.sif", "--camera", "CCD"])
        assert args.camera == "CCD"

    def test_invalid_camera_exits(self):
        with pytest.raises(SystemExit):
            _build_parser().parse_args(["shot.sif", "--camera", "BAD"])

    def test_overwrite_flag(self):
        args = _build_parser().parse_args(["shot.sif", "--overwrite"])
        assert args.overwrite is True

    def test_dry_run_flag(self):
        args = _build_parser().parse_args(["shot.sif", "--dry-run"])
        assert args.dry_run is True

    def test_verbose_flag(self):
        args = _build_parser().parse_args(["shot.sif", "--verbose"])
        assert args.verbose is True

    def test_instrument_id(self):
        args = _build_parser().parse_args(["shot.sif", "--instrument-id", "BlackEchelle"])
        assert args.instrument_id == "BlackEchelle"

    def test_custom_pattern(self):
        args = _build_parser().parse_args(["/data/", "--pattern", "*.sif"])
        assert args.pattern == "*.sif"

    def test_calibration_dir(self):
        args = _build_parser().parse_args(["shot.sif", "--calibration-dir", "/cal"])
        assert args.calibration_dir == "/cal"

    def test_frame_flag_accepted(self):
        args = _build_parser().parse_args(["shot.sif", "--frame", "3"])
        assert args.frame == 3


# ---------------------------------------------------------------------------
# _output_path_for
# ---------------------------------------------------------------------------


class TestOutputPathFor:
    def test_stem_preserved(self):
        sif = Path("/data/shot_042_Echelle.SIF")
        result = _output_path_for(sif, Path("/out"))
        assert result == Path("/out/shot_042_Echelle_spectrocube.nc")

    def test_suffix_is_nc(self):
        result = _output_path_for(Path("foo.sif"), Path("/out"))
        assert result.suffix == ".nc"

    def test_output_in_given_dir(self):
        result = _output_path_for(Path("/data/abc.sif"), Path("/tmp/out"))
        assert result.parent == Path("/tmp/out")


# ---------------------------------------------------------------------------
# Dry-run integration tests (no real SIF loading, spectrocube mocked)
# ---------------------------------------------------------------------------


@pytest.fixture()
def mock_spectrocube():
    """Inject a fake spectrocube module into sys.modules."""
    fake = MagicMock()
    with patch.dict(sys.modules, {"spectrocube": fake}):
        yield fake


class TestDryRunSingleFile:
    def test_dry_run_prints_dry_line(self, tmp_path, capsys, mock_spectrocube):
        sif = tmp_path / "shot_042.SIF"
        sif.touch()
        out = tmp_path / "shot_042.nc"

        with pytest.raises(SystemExit) as exc:
            main([str(sif), "--dry-run", "-o", str(out)])

        assert exc.value.code == 0
        assert "DRY" in capsys.readouterr().out

    def test_dry_run_does_not_create_file(self, tmp_path, mock_spectrocube):
        sif = tmp_path / "shot_042.SIF"
        sif.touch()
        out = tmp_path / "shot_042.nc"

        with pytest.raises(SystemExit):
            main([str(sif), "--dry-run", "-o", str(out)])

        assert not out.exists()


class TestDryRunBatch:
    def test_batch_dry_run_two_files(self, tmp_path, capsys, mock_spectrocube):
        (tmp_path / "a.SIF").touch()
        (tmp_path / "b.SIF").touch()
        out_dir = tmp_path / "out"

        with pytest.raises(SystemExit) as exc:
            main([str(tmp_path), "--dry-run", "-o", str(out_dir)])

        assert exc.value.code == 0
        stdout = capsys.readouterr().out
        assert "DRY" in stdout
        assert not out_dir.exists()

    def test_batch_lowercase_sif_fallback(self, tmp_path, capsys, mock_spectrocube):
        (tmp_path / "a.sif").touch()
        (tmp_path / "b.sif").touch()

        with pytest.raises(SystemExit) as exc:
            main([str(tmp_path), "--dry-run"])

        assert exc.value.code == 0
        stdout = capsys.readouterr().out
        assert "DRY" in stdout

    def test_batch_no_files_exits_1(self, tmp_path, mock_spectrocube):
        with pytest.raises(SystemExit) as exc:
            main([str(tmp_path), "--dry-run"])
        assert exc.value.code == 1

    def test_batch_custom_pattern_matched(self, tmp_path, capsys, mock_spectrocube):
        (tmp_path / "shot.tif").touch()
        with pytest.raises(SystemExit) as exc:
            main([str(tmp_path), "--pattern", "*.tif", "--dry-run"])
        assert exc.value.code == 0
        assert "DRY" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Missing spectrocube
# ---------------------------------------------------------------------------


class TestMissingSpectrocube:
    def test_exits_1_when_spectrocube_missing(self, tmp_path, capsys):
        sif = tmp_path / "shot.sif"
        sif.touch()

        # Remove spectrocube from sys.modules so the import inside main fails
        with patch.dict(sys.modules, {"spectrocube": None}):
            with pytest.raises(SystemExit) as exc:
                main([str(sif)])

        assert exc.value.code == 1
        assert "spectrocube" in capsys.readouterr().err.lower()


# ---------------------------------------------------------------------------
# Missing input
# ---------------------------------------------------------------------------


class TestMissingInput:
    def test_nonexistent_path_exits_1(self, tmp_path, mock_spectrocube):
        fake = str(tmp_path / "does_not_exist.sif")
        with pytest.raises(SystemExit) as exc:
            main([fake])
        assert exc.value.code == 1


# ---------------------------------------------------------------------------
# Skip-existing behaviour in dry-run
# ---------------------------------------------------------------------------


class TestSkipExisting:
    def test_existing_output_skipped_without_overwrite(
        self, tmp_path, capsys, mock_spectrocube
    ):
        sif = tmp_path / "shot.SIF"
        sif.touch()
        out = tmp_path / "shot_spectrocube.nc"
        out.touch()  # already exists

        with pytest.raises(SystemExit) as exc:
            main([str(sif), "--dry-run", "-o", str(out), "--verbose"])

        assert exc.value.code == 0
        # Should print SKIP, not DRY (file exists, no overwrite)
        assert "SKIP" in capsys.readouterr().out

    def test_existing_output_overwritten_with_flag(
        self, tmp_path, capsys, mock_spectrocube
    ):
        sif = tmp_path / "shot.SIF"
        sif.touch()
        out = tmp_path / "shot_spectrocube.nc"
        out.touch()

        with pytest.raises(SystemExit) as exc:
            main([str(sif), "--dry-run", "--overwrite", "-o", str(out)])

        assert exc.value.code == 0
        assert "DRY" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# export_spectrocube backend invocation (monkeypatched, no real SIF)
# ---------------------------------------------------------------------------


def _make_synthetic_spectrum():
    """Return a minimal duck-typed Spectrum-like object."""

    class _FakeSpectrum:
        pass

    sp = _FakeSpectrum()
    n = 50
    sp.wavelength = np.linspace(400.0, 700.0, n)
    sp.counts = np.random.default_rng(1).integers(100, 10000, (2, n)).astype(float)
    sp.wm = sp.counts * 1e-6
    sp.wmsr = sp.counts * 1e-7
    sp.phmsr = sp.counts * 1e-12
    sp.info = {"NumberOfFrames": 2, "ExposureTime": 0.5, "CycleTime": 1.0,
               "BackgroundFrames": []}
    sp.fpth = Path("/fake/shot.sif")
    sp.shotnumber = None
    return sp


class TestExportBackendCalled:
    def test_export_spectrocube_called_on_success(self, tmp_path, mock_spectrocube):
        """When load_spectrum succeeds, export_spectrocube must be called."""
        sif = tmp_path / "shot.SIF"
        sif.touch()
        out = tmp_path / "shot_spectrocube.nc"

        with patch(
            "echelle_spectra.spectrocube_cli._export_one",
            return_value=True,
        ) as mock_export:
            with pytest.raises(SystemExit) as exc:
                main([str(sif), "-o", str(out)])

        assert exc.value.code == 0
        mock_export.assert_called_once()

    def test_units_forwarded_to_export_one(self, tmp_path, mock_spectrocube):
        sif = tmp_path / "shot.SIF"
        sif.touch()

        with patch(
            "echelle_spectra.spectrocube_cli._export_one",
            return_value=True,
        ) as mock_export:
            with pytest.raises(SystemExit):
                main([str(sif), "--units", "wmsr"])

        call_kwargs = mock_export.call_args.kwargs
        assert call_kwargs["units"] == "wmsr"

    def test_instrument_id_forwarded(self, tmp_path, mock_spectrocube):
        sif = tmp_path / "shot.SIF"
        sif.touch()

        with patch(
            "echelle_spectra.spectrocube_cli._export_one",
            return_value=True,
        ) as mock_export:
            with pytest.raises(SystemExit):
                main([str(sif), "--instrument-id", "BlackEchelle"])

        call_kwargs = mock_export.call_args.kwargs
        assert call_kwargs["instrument_id"] == "BlackEchelle"
