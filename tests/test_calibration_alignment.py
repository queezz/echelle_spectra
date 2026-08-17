"""Tests for headless calibration-alignment helpers."""

from __future__ import annotations

import numpy as np

from echelle_spectra.tools.calibration_alignment import (
    AlignmentSettings,
    CalibrationTableLine,
    RigidTransform,
    apply_rigid_correction_to_lines,
    apply_rigid_correction_to_pattern,
    detector_points_from_lines,
    fit_rigid_transform,
    fit_single_gaussian_centroid,
    load_alignment_settings,
    load_wavelength_table,
    measure_detector_window_saturation,
    measure_line_centroids,
    measure_line_window_stats,
    rank_candidate_lines,
    save_alignment_settings,
    select_candidate_lines,
    write_corrected_pattern_table,
    write_pattern_table,
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


def test_pattern_correction_moves_every_trace_by_dy():
    pattern = np.zeros((100, 2), dtype=int)
    pattern[:, 0] = 20
    pattern[:, 1] = 40
    moved = apply_rigid_correction_to_pattern(
        pattern, RigidTransform(dx_px=0.0, dy_px=3.0, theta_rad=0.0)
    )
    assert moved.dtype == pattern.dtype
    np.testing.assert_array_equal(moved[:, 0], np.full(100, 23))
    np.testing.assert_array_equal(moved[:, 1], np.full(100, 43))


def test_pattern_correction_of_a_flat_trace_ignores_the_dispersion_shift():
    """A flat trace slid along itself is the same trace, dy apart."""

    pattern = np.full((100, 1), 30, dtype=int)
    moved = apply_rigid_correction_to_pattern(
        pattern, RigidTransform(dx_px=7.0, dy_px=-2.0, theta_rad=0.0)
    )
    np.testing.assert_array_equal(moved[:, 0], np.full(100, 28))


def test_pattern_correction_reads_a_sloped_trace_at_its_moved_column():
    """The corrected row at column x is the trace's row, not the column's."""

    columns = np.arange(50, dtype=float)
    pattern = np.column_stack([10.0 + 0.5 * columns])
    moved = apply_rigid_correction_to_pattern(
        pattern, RigidTransform(dx_px=4.0, dy_px=0.0, theta_rad=0.0)
    )
    # Away from the edges the trace slid four columns to the right, which for a
    # 0.5 px/column slope is two rows down where it is read.
    np.testing.assert_allclose(moved[10:, 0], pattern[10:, 0] - 2.0)
    # The four columns the transform slid in from beyond the left edge hold the
    # first corrected row rather than extrapolating the slope off the detector.
    np.testing.assert_allclose(moved[:4, 0], moved[4, 0])


def test_pattern_correction_keeps_rows_that_leave_the_detector():
    """A trace pushed past the top edge is still described where it went."""

    pattern = np.full((20, 1), 5, dtype=int)
    moved = apply_rigid_correction_to_pattern(
        pattern, RigidTransform(dx_px=0.0, dy_px=-9.0, theta_rad=0.0)
    )
    np.testing.assert_array_equal(moved[:, 0], np.full(20, -4))


def test_zero_transform_copies_the_base_pattern_byte_for_byte(tmp_path):
    base = tmp_path / "pattern.txt"
    base.write_text("10 40\r\n11 41\r\n", encoding="utf-8")
    target = tmp_path / "corrected.txt"

    correction = write_corrected_pattern_table(
        base, target, transform=RigidTransform(0.0, 0.0, 0.0)
    )

    assert not correction.applied
    assert "moves no trace" in correction.reason
    assert target.read_bytes() == base.read_bytes()


def test_unreadable_pattern_is_copied_and_says_so(tmp_path):
    base = tmp_path / "pattern.txt"
    base.write_text("pattern\n", encoding="utf-8")
    target = tmp_path / "corrected.txt"

    correction = write_corrected_pattern_table(
        base, target, transform=RigidTransform(0.0, 3.0, 0.0)
    )

    assert not correction.applied
    assert "no order columns" in correction.reason
    assert target.read_bytes() == base.read_bytes()


def test_corrected_pattern_is_written_as_integers_readers_can_load(tmp_path):
    base = tmp_path / "pattern.txt"
    write_pattern_table(np.column_stack([np.full(6, 10), np.full(6, 40)]), base)
    target = tmp_path / "corrected.txt"

    correction = write_corrected_pattern_table(
        base,
        target,
        transform=RigidTransform(0.0, 2.4, 0.0),
        metadata=[("Alignment dataset", "20250926")],
    )

    assert correction.applied
    assert correction.order_count == 2
    assert correction.max_shift_px == 2.0
    assert "# Alignment dataset: 20250926" in target.read_text(encoding="utf-8")
    loaded = np.loadtxt(target, dtype=int)
    np.testing.assert_array_equal(loaded[:, 0], np.full(6, 12))
    np.testing.assert_array_equal(loaded[:, 1], np.full(6, 42))


def test_detector_points_from_pattern():
    pattern = np.zeros((100, 2))
    pattern[:, 1] = np.arange(100)
    pts = detector_points_from_lines([_line(order=1, center=25.2)], pattern)
    np.testing.assert_allclose(pts, [[25.2, 25.0]])


def test_settings_round_trip(tmp_path):
    p = tmp_path / "lhd_cmos_alignment_20250926.settings.toml"
    settings = AlignmentSettings(
        instrument_id="lhd_cmos",
        created_at="2026-06-04",
        alignment_dataset_id="20250926",
        alignment_source_dir="local/20250926_calib",
        alignment_lamp="Ne",
        signal_file="Ne-0.02s-x3-bright-lines.sif",
        background_file="Ne-0.02s-x3-bright-lines_bg.sif",
        base_wavelength_file="Th_wavelength_CMOS_20240305.txt",
        base_pattern_file="pattern_CMOS_20240305.txt",
        sphere_file="sphere-0.1s-x3.sif",
        sphere_background_file="sphere-0.1s-x3-bg.sif",
        output_wavelength_file="Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        output_pattern_file="pattern_CMOS_20240305_aligned_to_20250926.txt",
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
    write_wavelength_table(
        [_line(order=2, center=123.456)],
        p,
        metadata=[("Alignment dataset", "20250926")],
    )
    assert "# Alignment dataset: 20250926" in p.read_text(encoding="utf-8")
    rows = load_wavelength_table(p)
    assert len(rows) == 1
    assert rows[0].order_idx == 2
    assert abs(rows[0].center_pixel - 123.456) < 1e-4
