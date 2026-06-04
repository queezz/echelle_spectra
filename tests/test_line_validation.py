"""Tests for known-line wavelength validation helpers."""

from __future__ import annotations

import numpy as np

from echelle_spectra.tools.line_validation import (
    LineValidationTarget,
    balmer_air_targets,
    fit_validation_line,
    load_fulcher_h2_q_branch_targets,
    summarize_validation,
    validate_lines,
)


def test_balmer_targets_are_air_wavelengths():
    targets = balmer_air_targets()
    assert [target.label for target in targets] == ["H-alpha", "H-beta", "H-gamma"]
    assert targets[0].wavelength_nm == 656.279
    assert {target.wavelength_medium for target in targets} == {"air"}


def test_fulcher_table_loader_labels_q_branch_rows(tmp_path):
    table = tmp_path / "fulcher.txt"
    table.write_text(
        "# Q rows, band columns\n"
        "601.0 612.0 0 632.0\n"
        "602.0 0 623.0 633.0\n",
        encoding="utf-8",
    )
    targets = load_fulcher_h2_q_branch_targets(table)
    assert [target.label for target in targets] == [
        "Fulcher H2 Q1(0-0)",
        "Fulcher H2 Q2(0-0)",
        "Fulcher H2 Q1(1-1)",
        "Fulcher H2 Q2(2-2)",
        "Fulcher H2 Q1(3-3)",
        "Fulcher H2 Q2(3-3)",
    ]


def test_fit_validation_line_recovers_wavelength_centroid():
    wavelength = np.linspace(655.8, 656.6, 300)
    intensity = 3.0 + 80.0 * np.exp(-0.5 * ((wavelength - 656.2812) / 0.025) ** 2)
    result = fit_validation_line(
        wavelength,
        intensity,
        LineValidationTarget("H-alpha", 656.279),
        window_nm=0.2,
        min_snr=3.0,
    )
    assert result is not None
    assert abs(result.measured_nm - 656.2812) < 0.002
    assert result.peak_snr > 3


def test_validate_lines_reports_multiframe_fit_count_and_summary():
    wavelength = np.linspace(485.9, 486.3, 240)
    centers = [486.134, 486.135, 486.136, 486.137]
    spectra = np.array(
        [
            2.0 + 50.0 * np.exp(-0.5 * ((wavelength - center) / 0.018) ** 2)
            for center in centers
        ]
    )
    results = validate_lines(
        wavelength,
        spectra,
        [LineValidationTarget("H-beta", 486.135)],
        signal_frames=[0, 1, 2],
        window_nm=0.15,
        min_snr=3.0,
    )
    assert len(results) == 1
    assert results[0].frames_fit == 3
    assert results[0].frames_total == 3
    assert results[0].frame_centroid_std_nm is not None
    summary = summarize_validation(results)
    assert summary.n_lines == 1
    assert summary.rms_residual_nm < 0.002
