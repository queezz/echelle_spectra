"""Tests for the echelle-align CLI."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from echelle_spectra.alignment_cli import _build_parser, main
from echelle_spectra.tools.calibration_alignment import (
    AlignmentRunResult,
    AlignmentSettings,
    CalibrationTableLine,
    DetectorWindowSaturation,
    LineCentroidFit,
    LineWindowStats,
    RigidTransform,
)


def _touch_inputs(tmp_path):
    paths = {
        "signal": tmp_path / "Ne.sif",
        "background": tmp_path / "Ne_bg.sif",
        "sphere": tmp_path / "sphere.sif",
        "sphere_background": tmp_path / "sphere_bg.sif",
        "wavelength": tmp_path / "Th_wavelength_CMOS_20240305.txt",
        "pattern": tmp_path / "pattern_CMOS_20250926.txt",
        "integral": tmp_path / "integrating_sphere.txt",
    }
    for path in paths.values():
        path.touch()
    return paths


def _fake_result():
    line = CalibrationTableLine(0, 10, 20, 15, 650.0, "NeI", "ok")
    settings = AlignmentSettings(
        instrument_id="lhd_cmos",
        created_at="2026-06-04",
        alignment_dataset_id="20250926",
        alignment_source_dir="local/20250926_calib",
        alignment_lamp="Ne",
        signal_file="Ne.sif",
        background_file="Ne_bg.sif",
        base_wavelength_file="Th_wavelength_CMOS_20240305.txt",
        base_pattern_file="pattern_CMOS_20250926.txt",
        sphere_file="sphere.sif",
        sphere_background_file="sphere_bg.sif",
        output_wavelength_file="adjusted.txt",
        output_pattern_file="adjusted_pattern.txt",
        transform=RigidTransform(-14.0, 0.3, 0.001),
        n_lines=2,
        rms_px=0.5,
        notes="unit test",
    )
    fit = LineCentroidFit(line, 14.0, 1.0, 10.0, 0.0, 20.0, True)
    window = LineWindowStats(line, 14.0, 10.0, 0.0, 1.0, 10.0, 20.0, 10, 0, 0.0, True)
    saturation = DetectorWindowSaturation(line, 100.0, 10, 0, 0.0)
    return AlignmentRunResult(
        settings=settings,
        adjusted_rows=[line],
        rows=[line],
        candidates=[line],
        fits=[fit],
        good_fits=[fit],
        residual_xy=np.array([[0.3, 0.4]]),
        residual_px=np.array([0.5]),
        expected_xy=np.array([[15.0, 20.0]]),
        measured_xy=np.array([[14.0, 20.0]]),
        predicted_xy=np.array([[14.3, 20.4]]),
        detector_saturation=[saturation],
        ranked_stats=[window],
    )


def test_parser_defaults():
    args = _build_parser().parse_args(["Ne.sif", "Ne_bg.sif", "sphere.sif", "sphere_bg.sif"])
    assert args.pattern == "pattern_CMOS_20250926.txt"
    assert args.wavelength == "Th_wavelength_CMOS_20240305.txt"
    assert not args.save


def test_missing_input_exits(tmp_path):
    paths = _touch_inputs(tmp_path)
    paths["signal"].unlink()
    with pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
            ]
        )
    assert exc.value.code == 1


def test_preview_mode_does_not_write(tmp_path, capsys):
    paths = _touch_inputs(tmp_path)
    with patch(
        "echelle_spectra.alignment_cli.run_calibration_alignment",
        return_value=_fake_result(),
    ), \
        pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
            ]
        )

    assert exc.value.code == 0
    assert "Preview only" in capsys.readouterr().out
    assert not (tmp_path / "alignments").exists()


def test_save_writes_settings_and_table(tmp_path):
    paths = _touch_inputs(tmp_path)
    settings = tmp_path / "settings.toml"
    table = tmp_path / "adjusted.txt"

    with patch(
        "echelle_spectra.alignment_cli.run_calibration_alignment",
        return_value=_fake_result(),
    ), \
        pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
                "--settings-out",
                str(settings),
                "--table-out",
                str(table),
                "--save",
            ]
        )

    assert exc.value.code == 0
    assert "base_pattern_file = \"pattern_CMOS_20250926.txt\"" in settings.read_text()
    assert "# Settings file: settings.toml" in table.read_text()


def test_save_refuses_existing_outputs_without_overwrite(tmp_path):
    paths = _touch_inputs(tmp_path)
    settings = tmp_path / "settings.toml"
    table = tmp_path / "adjusted.txt"
    settings.touch()
    table.touch()

    with patch("echelle_spectra.alignment_cli.run_calibration_alignment") as run_alignment, \
        pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
                "--settings-out",
                str(settings),
                "--table-out",
                str(table),
                "--save",
            ]
        )

    assert exc.value.code == 1
    run_alignment.assert_not_called()


def _write_pattern(path, columns=40, orders=2):
    rows = np.arange(columns, dtype=float)
    pattern = np.rint(
        np.column_stack([20.0 + 0.2 * rows, 60.0 + 0.05 * rows])[:, :orders]
    ).astype(int)
    np.savetxt(path, pattern, fmt="%d")
    return pattern


def test_save_writes_a_corrected_pattern_beside_the_corrected_table(tmp_path):
    """The CLI and the bench must not disagree about which detector they mean."""

    paths = _touch_inputs(tmp_path)
    base_pattern = _write_pattern(paths["pattern"])
    settings = tmp_path / "settings.toml"
    table = tmp_path / "adjusted.txt"
    pattern_out = tmp_path / "adjusted_pattern.txt"

    with patch(
        "echelle_spectra.alignment_cli.run_calibration_alignment",
        return_value=_fake_result(),
    ), \
        pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
                "--settings-out",
                str(settings),
                "--table-out",
                str(table),
                "--pattern-out",
                str(pattern_out),
                "--save",
            ]
        )

    assert exc.value.code == 0
    saved = np.loadtxt(pattern_out, dtype=int)
    assert saved.shape == base_pattern.shape
    # dx -14 px against a sloping trace, so the pattern really moved.
    assert not np.array_equal(saved, base_pattern)
    assert "# Base pattern file: pattern_CMOS_20250926.txt" in pattern_out.read_text()
    assert 'output_pattern_file = "adjusted_pattern.txt"' in settings.read_text()


def test_save_names_the_pattern_after_the_dataset_by_default(tmp_path):
    paths = _touch_inputs(tmp_path)
    _write_pattern(paths["pattern"])

    with patch(
        "echelle_spectra.alignment_cli.run_calibration_alignment",
        return_value=_fake_result(),
    ), \
        pytest.raises(SystemExit):
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
                "--save",
            ]
        )

    expected = tmp_path / "alignments" / "pattern_CMOS_20250926_aligned_to_20250926.txt"
    assert expected.is_file()


def test_save_refuses_an_existing_pattern_without_overwrite(tmp_path):
    paths = _touch_inputs(tmp_path)
    _write_pattern(paths["pattern"])
    pattern_out = tmp_path / "adjusted_pattern.txt"
    pattern_out.touch()

    with patch("echelle_spectra.alignment_cli.run_calibration_alignment") as run_alignment, \
        pytest.raises(SystemExit) as exc:
        main(
            [
                str(paths["signal"]),
                str(paths["background"]),
                str(paths["sphere"]),
                str(paths["sphere_background"]),
                "--calibration-dir",
                str(tmp_path),
                "--pattern-out",
                str(pattern_out),
                "--save",
            ]
        )

    assert exc.value.code == 1
    run_alignment.assert_not_called()
