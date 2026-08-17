"""Owner rider, 2026-08-17: "can we do anything about auto histogram?"

The viewer levelled the detector image between the 1st and 99th percentile of
the frame.  On an echelle calibration frame both of those percentiles are
*inside the background* — measured on the owner's own pair, 96 and 149 counts
on the bright frame — because fewer than one pixel in a thousand carries line
light at all.  The result is the grey wash he photographed: read noise stretched
across the whole ramp and every line flat white.

These pin the replacement rule and the one thing it must not do, which is fight
the operator's hand on the histogram.
"""

from __future__ import annotations

import numpy as np
import pytest

from echelle_spectra.tools.detector_display import (
    AUTO_LEVEL_HIGH_PERCENTILE,
    AUTO_LEVEL_LOW_PERCENTILE,
    auto_display_levels,
)

BACKGROUND = 120.0
NOISE = 6.0
FAINTEST_SPIKE = 400.0
SPIKE = 20000.0
SHAPE = (2160, 2560)


def _frame(spike_fraction: float = 0.005) -> np.ndarray:
    """A flat background with read noise and a sparse spread of bright spikes.

    The spikes run over two decades, as a lamp's lines do, so the white point
    has somewhere to land *between* the background and the brightest of them.
    """

    generator = np.random.default_rng(20250926)
    frame = generator.normal(BACKGROUND, NOISE, SHAPE)
    count = int(frame.size * spike_fraction)
    rows = generator.integers(0, SHAPE[0], count)
    columns = generator.integers(0, SHAPE[1], count)
    heights = np.exp(
        generator.uniform(np.log(FAINTEST_SPIKE), np.log(SPIKE), count)
    )
    frame[rows, columns] = heights
    return frame


def test_the_levels_land_between_the_background_and_the_spikes():
    frame = _frame()

    levels = auto_display_levels(frame)

    assert levels is not None
    low, high = levels
    # Black at the background floor: the median, not the minimum and not a
    # percentile still deep in the noise.
    assert low == pytest.approx(BACKGROUND, abs=0.5)
    assert low == pytest.approx(np.percentile(frame, AUTO_LEVEL_LOW_PERCENTILE), abs=0.5)
    # White above the noise and below the brightest spikes, so the strong ones
    # clip and everything between the two stays visible.
    assert BACKGROUND + 5 * NOISE < high < SPIKE
    assert high == pytest.approx(
        np.percentile(frame, AUTO_LEVEL_HIGH_PERCENTILE), rel=0.2
    )
    # What the percentile is actually promising: about a thousandth of the
    # frame clips, and it is spikes that clip, never background.
    clipped = float(np.mean(frame > high))
    assert 0.0002 < clipped < 0.004
    assert frame[frame > high].min() > FAINTEST_SPIKE


def test_the_white_point_follows_the_light_rather_than_the_frame_size():
    """A brighter exposure of the same scene levels the same way, scaled up.

    The owner's pair is one lamp at two exposures; the rule has to render them
    alike rather than hand the longer one a washed-out ramp.
    """

    dim = _frame()
    bright = (dim - BACKGROUND) * 5.0 + BACKGROUND

    dim_low, dim_high = auto_display_levels(dim)
    bright_low, bright_high = auto_display_levels(bright)

    assert bright_low == pytest.approx(dim_low, abs=1.0)
    assert (bright_high - bright_low) == pytest.approx(5.0 * (dim_high - dim_low), rel=0.05)


def test_a_subsample_answers_the_same_question_as_the_whole_frame():
    """Levels come off a strided slice, so they must not depend on the stride.

    This is what makes it affordable to level on every frame change.  The floor
    is a median and lands on the same count either way; the white point is a
    far-tail percentile of a two-decade spread and carries real sampling
    scatter, so what has to hold of it is what it is *for* — the same sliver of
    the frame clips whichever stride read it.
    """

    frame = _frame()

    sampled = auto_display_levels(frame, sample_target=50_000)
    full = auto_display_levels(frame, sample_target=frame.size)

    assert sampled is not None and full is not None
    assert sampled[0] == pytest.approx(full[0], abs=0.5)
    for _low, high in (sampled, full):
        assert 0.0002 < float(np.mean(frame > high)) < 0.004


def test_a_frame_with_nothing_to_say_leaves_the_histogram_alone():
    assert auto_display_levels(np.full(SHAPE, 3.0)) is None
    assert auto_display_levels(np.empty((0, 0))) is None
    assert auto_display_levels(np.full(SHAPE, np.nan)) is None
    assert auto_display_levels(np.arange(10.0)) is None
