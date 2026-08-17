"""Display levels for the raw 2-D detector frame.

An echelle calibration frame is not a picture, and the level rule that suits a
picture ruins it.  Measured on the owner's 20250926 pair
(``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif`` and
``Ne-0.1s-x3-dimm-lines.sif``, both 2560x2160): the frame is one enormous
near-background mass — the median sits at 119 and 122 counts, and the 99th
percentile is still only 149 and 227 — with a few dozen lamp blobs two to three
decades above it, peaking at 51 000 and 65 535.  Fewer than one pixel in a
thousand carries line light at all.

So a symmetric quantile window such as the 1st--99th percentile the viewer used
to apply lands *entirely inside the background*: 96--149 counts on the bright
frame.  Every line saturates to white, every blob loses its shape, and the
sensor's read noise is stretched across the whole colour ramp — the grey wash
the owner reported.

The rule below is the shape of the data instead:

* the low level is the **background floor**, the median of a subsample.  It is
  robust — the minimum is one cold pixel and the 1st percentile is still noise
  — and putting it at the median renders the background as black rather than as
  a texture competing with the lines.
* the high level is a **high percentile**, so the top ~0.1% of pixels clip to
  white.  On this data clipping is correct: the clipped pixels are the cores of
  the strong lamp lines, which are the point of the frame and are saturated in
  the sensor as often as not, and refusing to clip them costs every fainter
  line its visibility.

Both are read from a strided subsample rather than the full 5.5-megapixel
frame, because these levels are recomputed on every frame change and a full
sort is not worth a hundredth of a second of the operator's time.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "AUTO_LEVEL_HIGH_PERCENTILE",
    "AUTO_LEVEL_LOW_PERCENTILE",
    "AUTO_LEVEL_SAMPLE_TARGET",
    "auto_display_levels",
]

#: The background floor: half the pixels of an echelle frame are background and
#: nothing else, so their median *is* the floor.
AUTO_LEVEL_LOW_PERCENTILE = 50.0

#: The white point.  Tuned on the owner's real pair: it puts the bright frame at
#: 119--2912 counts and the dim frame at 122--14499.  Of the blobs those frames
#: actually carry (45 and 56 connected regions above three times the floor),
#: every one of the bright frame's and 50 of the dim frame's render above a
#: twentieth of full scale, 41 and 39 above a fifth of it, and the strongest 29
#: and 27 clip to white.  99.5 would leave the white point at 188 counts on the
#: bright frame — inside the noise, the old defect again — and 99.99 pushes it
#: into the saturated cores and takes the faint lines down with it.
AUTO_LEVEL_HIGH_PERCENTILE = 99.9

#: About how many pixels the percentiles are read from.  Two hundred thousand
#: is far more than either percentile needs to be stable and costs a few
#: milliseconds against roughly a hundred for the whole frame.
AUTO_LEVEL_SAMPLE_TARGET = 200_000


def _subsample(image: np.ndarray, target: int) -> np.ndarray:
    """A strided slice of about *target* finite samples, never a copy-sort."""

    flat_size = int(image.size)
    if flat_size <= target:
        sample = np.asarray(image).ravel()
    else:
        stride = max(1, int(np.sqrt(flat_size / float(target))))
        sample = np.asarray(image)[::stride, ::stride].ravel()
    finite = sample[np.isfinite(sample)]
    return finite


def auto_display_levels(
    image,
    *,
    low_percentile: float = AUTO_LEVEL_LOW_PERCENTILE,
    high_percentile: float = AUTO_LEVEL_HIGH_PERCENTILE,
    sample_target: int = AUTO_LEVEL_SAMPLE_TARGET,
):
    """Return ``(low, high)`` display levels for one detector frame.

    ``None`` when the frame carries nothing to level — empty, all-NaN, or
    perfectly flat — so the caller can leave the histogram exactly as the
    operator left it rather than collapse it onto a single value.
    """

    frame = np.asarray(image)
    if frame.ndim != 2 or frame.size == 0:
        return None
    sample = _subsample(frame, int(sample_target))
    if sample.size == 0:
        return None
    low = float(np.percentile(sample, float(low_percentile)))
    high = float(np.percentile(sample, float(high_percentile)))
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return None
    return low, high
