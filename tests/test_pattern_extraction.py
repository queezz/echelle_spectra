"""Tests for headless order-pattern extraction helpers."""

from __future__ import annotations

import numpy as np

from echelle_spectra.tools.pattern_extraction import (
    PatternExtractionConfig,
    average_detector_frames,
    detect_order_peaks_at_column,
    detect_order_peaks_near_prior_at_column,
    extract_order_pattern,
    extract_order_pattern_near_prior,
    sample_columns,
    subtract_background,
    trial_order_pattern_extraction,
)


def _synthetic_order_image(nrows=180, ncols=240, n_orders=5):
    rows = np.arange(nrows, dtype=float)
    cols = np.arange(ncols, dtype=float)
    image = np.full((nrows, ncols), 10.0)
    truth = np.zeros((ncols, n_orders), dtype=float)

    for order_idx in range(n_orders):
        base = 24.0 + order_idx * 27.0
        slope = 0.015 * (order_idx - 2)
        curve = 0.00018 * (cols - ncols / 2.0) ** 2
        centers = base + slope * cols + curve
        truth[:, order_idx] = centers
        for col, center in enumerate(centers):
            image[:, col] += 300.0 * np.exp(-0.5 * ((rows - center) / 1.8) ** 2)

    return image, truth


def _config(expected_orders=5):
    return PatternExtractionConfig(
        expected_orders=expected_orders,
        sample_step_px=35,
        sample_count=5,
        smooth_window_px=9,
        smooth_polyorder=1,
        amplification_rate=0.0,
        peak_threshold=0.2,
        peak_min_dist_px=15,
        baseline_poly_deg=2,
        trace_poly_degree=2,
    )


def test_average_detector_frames_accepts_2d_or_3d():
    image = np.arange(6).reshape(2, 3)
    stack = np.stack([image, image + 2])

    np.testing.assert_allclose(average_detector_frames(image), image)
    np.testing.assert_allclose(average_detector_frames(stack), image + 1)


def test_subtract_background_averages_frame_stacks():
    signal = np.ones((3, 4, 5)) * 12
    background = np.ones((2, 4, 5)) * 2

    result = subtract_background(signal, background)

    assert result.shape == (4, 5)
    np.testing.assert_allclose(result, 10)


def test_sample_columns_matches_notebook_centered_spacing():
    columns = sample_columns(2560, step_size=150, num_steps=10)

    np.testing.assert_array_equal(
        columns,
        np.array([530, 680, 830, 980, 1130, 1280, 1430, 1580, 1730, 1880]),
    )


def test_detect_order_peaks_at_column_finds_expected_count():
    image, truth = _synthetic_order_image()
    col = 120

    detection = detect_order_peaks_at_column(image, col, config=_config())

    assert detection.column_px == col
    assert detection.n_peaks == 5
    np.testing.assert_allclose(detection.row_peaks_px, np.rint(truth[col]), atol=1)


def test_extract_order_pattern_recovers_synthetic_traces():
    image, truth = _synthetic_order_image()

    result = extract_order_pattern(image, config=_config())

    assert result.pattern.shape == truth.shape
    assert result.n_orders == 5
    np.testing.assert_allclose(result.pattern, np.rint(truth), atol=2)


def test_extract_order_pattern_rejects_wrong_order_count():
    image, _truth = _synthetic_order_image()

    try:
        extract_order_pattern(image, config=_config(expected_orders=6))
    except ValueError as exc:
        assert "expected 6 orders" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_trial_order_pattern_extraction_ranks_successful_fit_first():
    image, _truth = _synthetic_order_image()
    config = _config()

    trials = trial_order_pattern_extraction(
        image,
        config=config,
        threshold_values=[0.99, config.peak_threshold],
        column_start_values=[50],
        column_step_px=35,
        sample_count=5,
    )

    assert trials[0].success
    assert trials[0].threshold == config.peak_threshold
    assert trials[0].result is not None
    assert not trials[-1].success
    assert trials[-1].error
    assert trials[-1].peak_counts.size == 5


def test_detect_order_peaks_near_prior_recovers_shifted_traces():
    image, truth = _synthetic_order_image()
    prior = np.rint(truth - 3).astype(int)
    col = 120

    detection = detect_order_peaks_near_prior_at_column(
        image,
        prior,
        col,
        config=_config(),
        search_radius_px=8,
    )

    assert detection.n_peaks == 5
    np.testing.assert_allclose(detection.row_peaks_px, np.rint(truth[col]), atol=1)


def test_extract_order_pattern_near_prior_recovers_synthetic_traces():
    image, truth = _synthetic_order_image()
    prior = np.rint(truth - 3).astype(int)

    result = extract_order_pattern_near_prior(
        image,
        prior,
        config=_config(),
        search_radius_px=8,
    )

    assert result.pattern.shape == truth.shape
    assert result.n_orders == 5
    np.testing.assert_allclose(result.pattern, np.rint(truth), atol=2)
