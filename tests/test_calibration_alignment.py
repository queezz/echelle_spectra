"""Tests for headless calibration-alignment helpers."""

from __future__ import annotations

import numpy as np

from echelle_spectra.tools.calibration_alignment import (
    AlignmentSettings,
    CalibrationTableLine,
    RigidTransform,
    apply_rigid_correction_to_lines,
    detector_points_from_lines,
    fit_rigid_transform,
    fit_single_gaussian_centroid,
    load_alignment_settings,
    load_wavelength_table,
    measure_detector_window_saturation,
    measure_line_window_stats,
    measure_line_centroids,
    rank_candidate_lines,
    save_alignment_settings,
    select_candidate_lines,
    write_wavelength_table,
)


def _line(order=0, center=50.0, species="NeI", comment="ok"):
    return CalibrationTableLine(
        order_idx=order,
        pixel_from=center - 5,
        pixel_to=center + 5,
        center_pixel=center,
        wavelength_nm=650.0,
        species=species,
        comment=comment,
    )


def test_wavelength_table_preserves_inline_comment(tmp_path):
    p = tmp_path / "wave.txt"
    p.write_text(
        "# header\n"
        "0 10 20 15.5 650.123 NeI # ok single\n"
        "#1 10 20 15.5 651.0 ArI # disabled\n",
        encoding="utf-8",
    )
    rows = load_wavelength_table(p)
    assert len(rows) == 1
    assert rows[0].species == "NeI"
    assert rows[0].comment == "ok single"


def test_candidate_selection_defaults_to_curated_ne_lines():
    rows = [
        _line(species="NeI", comment="ok"),
        _line(species="NeI", comment="?"),
        _line(species="ArI", comment="ok"),
        _line(species="NeI", center=50.0, comment="ok"),
    ]
    rows[-1] = CalibrationTableLine(0, 49, 50, 49.5, 650.0, "NeI", "ok")
    selected = select_candidate_lines(rows)
    assert selected == [rows[0]]


def test_single_gaussian_centroid_recovers_center():
    x = np.arange(101, dtype=float)
    y = 3.0 + 80.0 * np.exp(-0.5 * ((x - 42.35) / 1.8) ** 2)
    ok, center, sigma, amp, baseline, snr, reason = fit_single_gaussian_centroid(
        y,
        expected_center_px=42.0,
        min_snr=3.0,
    )
    assert ok, reason
    assert abs(center - 42.35) < 0.05
    assert 1.7 < sigma < 1.9
    assert amp > 70
    assert baseline > 0
    assert snr > 3


def test_single_gaussian_rejects_saturated_peak():
    y = np.zeros(100)
    y[40] = 65535
    ok, *_rest, reason = fit_single_gaussian_centroid(
        y,
        expected_center_px=40,
        saturation_level=60000,
    )
    assert not ok
    assert reason == "saturated"


def test_line_window_stats_report_saturation_fraction():
    line = _line(center=40)
    y = np.zeros(100)
    y[38:42] = 65535
    stats = measure_line_window_stats(
        [y],
        [line],
        window_radius_px=10,
        saturation_level=60000,
    )[0]
    assert stats.is_saturated
    assert stats.saturated_pixels == 4
    assert stats.saturated_fraction > 0
    assert not stats.fit_candidate
    assert stats.reason == "saturated"


def test_detector_window_saturation_uses_raw_2d_pixels():
    line = _line(order=0, center=10)
    images = np.zeros((2, 20, 30))
    pattern = np.full((30, 1), 8)
    images[1, 8, 10] = 65535
    stats = measure_detector_window_saturation(
        images,
        pattern,
        [line],
        x_radius_px=2,
        y_radius_px=1,
        saturation_level=60000,
    )[0]
    assert stats.is_saturated
    assert stats.saturated_pixels == 1
    assert stats.finite_pixels == 30


def test_rank_candidate_lines_prefers_high_snr_unsaturated_windows():
    strong = _line(center=40)
    weak = _line(center=70)
    x = np.arange(120, dtype=float)
    y = (
        2.0
        + 100.0 * np.exp(-0.5 * ((x - 40) / 1.5) ** 2)
        + 8.0 * np.exp(-0.5 * ((x - 70) / 1.5) ** 2)
    )
    ranked = rank_candidate_lines([y], [weak, strong], window_radius_px=10, min_snr=3.0)
    assert ranked[0].line == strong
    assert ranked[0].fit_candidate


def test_measure_line_centroids_uses_order_index():
    line = _line(order=1, center=30)
    spectra = [np.zeros(80), np.zeros(80)]
    x = np.arange(80)
    spectra[1] = 5 + 50 * np.exp(-0.5 * ((x - 31.2) / 1.5) ** 2)
    result = measure_line_centroids(spectra, [line], min_snr=3.0)
    assert result[0].success
    assert abs(result[0].center_pixel - 31.2) < 0.1


def test_fit_rigid_transform_recovers_translation_and_rotation():
    expected = np.array([[0, 0], [10, 0], [0, 20], [30, 40]], dtype=float)
    transform = RigidTransform(dx_px=3.0, dy_px=-2.0, theta_rad=np.deg2rad(1.5))
    measured = transform.apply(expected)
    fitted, rms = fit_rigid_transform(expected, measured)
    assert rms < 1e-10
    assert abs(fitted.dx_px - 3.0) < 1e-10
    assert abs(fitted.dy_px + 2.0) < 1e-10
    assert abs(fitted.theta_rad - transform.theta_rad) < 1e-10


def test_apply_correction_moves_lookup_center_pixels():
    pattern = np.zeros((100, 2))
    pattern[:, 0] = 20
    pattern[:, 1] = 40
    rows = [_line(order=0, center=10), _line(order=1, center=50)]
    adjusted = apply_rigid_correction_to_lines(
        rows,
        pattern,
        RigidTransform(dx_px=2.5, dy_px=0.0, theta_rad=0.0),
    )
    assert [r.center_pixel for r in adjusted] == [12.5, 52.5]
    assert adjusted[0].wavelength_nm == rows[0].wavelength_nm


def test_detector_points_from_pattern():
    pattern = np.zeros((100, 2))
    pattern[:, 1] = np.arange(100)
    pts = detector_points_from_lines([_line(order=1, center=25.2)], pattern)
    np.testing.assert_allclose(pts, [[25.2, 25.0]])


def test_settings_round_trip(tmp_path):
    p = tmp_path / "lhd_cmos_20240305.settings.toml"
    settings = AlignmentSettings(
        instrument_id="lhd_cmos",
        base_wavelength_file="Th_wavelength_CMOS_20240305.txt",
        n_lines=12,
        rms_px=0.3,
        notes="unit test",
        transform=RigidTransform(1.0, 2.0, 0.01),
    )
    save_alignment_settings(settings, p)
    loaded = load_alignment_settings(p)
    assert loaded == settings


def test_write_wavelength_table_round_trip(tmp_path):
    p = tmp_path / "adjusted.txt"
    write_wavelength_table([_line(order=2, center=123.456)], p)
    rows = load_wavelength_table(p)
    assert len(rows) == 1
    assert rows[0].order_idx == 2
    assert abs(rows[0].center_pixel - 123.456) < 1e-4
