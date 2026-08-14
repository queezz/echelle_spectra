"""Science-line drift sampling, verdicts, and immutable refinements.

The judge works in detector space.  A rigid detector shift moves every line by
the same number of *pixels*, never by the same number of nanometres: across the
audited Balmer/Fulcher lines the dispersion runs from about 0.0066 nm/px in the
blue orders to about 0.0108 nm/px at H-alpha, a ~60% spread.  Each measured
wavelength residual is therefore divided by the dispersion that the cube's own
stored per-order polynomial has at that line's detector pixel, one rigid shift
is fitted in pixels, and every verdict boundary is a pixel bound.

The emitted correction is per-order for the same reason.  An accepted shift
slides the wavelength table's anchors along the detector; each order's refitted
polynomial then turns those pixels into that order's own wavelength correction,
which is what ``echelle recal-cube`` later applies.  Correcting every order by
one scalar wavelength instead would leave the low-dispersion orders several
pixels wrong while the next audit called them aligned.
"""

from __future__ import annotations

import json
import re
import tempfile
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .campaign_run import sha256_file
from .snapshot import REQUIRED_ROLES, Snapshot, create_snapshot, load_snapshot
from .tools.calibration_alignment import (
    load_wavelength_table,
    shift_lines_in_pixels,
    write_wavelength_table,
)
from .tools.line_catalog import load_line_table

DRIFT_SCHEMA = "echelle-drift-evidence/v2"
DRIFT_SCHEMA_V1 = "echelle-drift-evidence/v1"
#: Schemas the bulk-processing gate still reads.  Only v2 is ever written: it
#: carries the detector-space residuals and per-order corrections that v1 could
#: not express.  v1 records stay readable so evidence taken before this packet
#: keeps authorizing the runs it already authorized.
READABLE_DRIFT_SCHEMAS = (DRIFT_SCHEMA, DRIFT_SCHEMA_V1)

# --- Verdict boundaries, in detector pixels --------------------------------
#
# ALIGNMENT_MAX_PIXEL_RESIDUAL — the finest disagreement a sampled audit may
# call a misalignment.  The bench solves its rigid alignment to about 0.67 px
# RMS on real Neon data (the shipped 20250926 alignment settings), so half a
# pixel is already below what the calibration itself is known to.  It is also
# 0.003-0.005 nm across the audited dispersion range, an order inside the
# ~0.03 nm instrumental line width.  Because the median can never exceed the
# maximum, this one bound states alignment completely.  A cube carrying
# ``wavelength_accuracy_nm`` relaxes the bound per line to whatever that
# accuracy is worth in pixels there: the judge never claims a misalignment
# finer than the calibration's own stated accuracy.
ALIGNMENT_MAX_PIXEL_RESIDUAL = 0.5

# SHIFT_CONSISTENCY_PX — how far one line's pixel residual may sit from the
# fitted rigid shift and still belong to the same rigid shift.  One pixel is
# ~1.5x the bench's alignment RMS and covers ordinary centroid noise; beyond it
# the residuals are not one translation of the detector, so no single shift
# repairs them and the interval is beyond repair (or wants splitting).
SHIFT_CONSISTENCY_PX = 1.0

# REPAIR_LIMIT_PX — the largest rigid shift a table correction may still claim.
# Two independent measurement windows bound it: this audit fits centroids in a
# +/-0.4 nm window, which is +/-37 px at H-alpha's 0.0108 nm/px (+/-60 px in
# the blue), and the bench fits its lamp centroids in a +/-18 px window.  At
# 25 px every audited line is still inside the window that measured it, with
# margin, and the corrected anchors still fall inside the pixel intervals the
# table records for them.  Beyond that the anchors must be re-identified
# against lamp data rather than slid.  25 px is also the honest pixel reading
# of the 0.25 nm limit this rule replaces: that limit was 23 px at H-alpha and
# 38 px in the blue, differing only by which order a line happened to fall in.
REPAIR_LIMIT_PX = 25.0

# --- Evidence sufficiency ---------------------------------------------------
#
# MINIMUM_CENTROIDS counts *resolved* centroids: lines the catalog does not
# flag as sub-resolution blends, measured at or above MINIMUM_SNR.  Three of
# them must come from three distinct Echelle orders and must not all sit in one
# half of the audited wavelength coverage.  A rigid shift is a one-parameter
# claim about the whole detector, and one order in one colour cannot see the
# dispersion spread that separates a real shift from a dispersion change; a
# duplicated or blended catalog row must never manufacture that quorum.
MINIMUM_CENTROIDS = 3
MINIMUM_DISTINCT_ORDERS = 3
MINIMUM_SNR = 4.0

# --- Plasma-presence frame selection ---------------------------------------
#
# Cube intensities are already background-subtracted, so a frame without plasma
# hovers about zero and the median of the per-frame totals is the background
# level.  A frame counts as plasma-bright when its total stands PLASMA_FRAME_-
# SIGMA robust deviations above that level, which is the same 5-sigma spirit
# ``Spectrum`` uses to recognise its dark frames.  A median over all frames
# instead would dilute every line with the dark frames beside it.
PLASMA_FRAME_SIGMA = 5.0

#: The cube attributes the audit reads for an acquisition date, in the order it
#: tries them: SpectroCube's ISO acquisition start, then a calendar date inside
#: the recorded source filename, then the export timestamp as a last resort.
DATE_ATTRIBUTES = ("t_start", "source_file", "created_at")

_DATE_PATTERN = re.compile(
    r"(?<!\d)(?P<year>20\d{2})[-_]?(?P<month>0[1-9]|1[0-2])"
    r"[-_]?(?P<day>0[1-9]|[12]\d|3[01])(?!\d)"
)


class DriftError(ValueError):
    """Drift evidence or refinement is unsafe or incomplete."""


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def resolve_cube_paths(inputs: list[str | Path]) -> list[Path]:
    """Expand directories to the cubes under them, shell-neutrally.

    A directory contributes every ``*.nc`` file beneath it in sorted order, so
    one audit command reads the same set on PowerShell and on a POSIX shell
    without either one's glob syntax.
    """

    resolved: list[Path] = []
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            found = sorted(path.rglob("*.nc"))
            if not found:
                raise DriftError(f"no .nc cubes found under {path}")
            resolved.extend(found)
        elif path.is_file():
            resolved.append(path)
        else:
            raise DriftError(f"cube path does not exist: {path}")
    if not resolved:
        raise DriftError("no cubes were selected")
    return sorted(dict.fromkeys(resolved))


def _shot_identity(path: Path) -> tuple[int | None, str]:
    """Return the parsed shot number and the raw stem token for one cube."""

    from .calibration_registry import CalibrationRegistryError, source_identity_from_path

    try:
        identity = source_identity_from_path(path)
    except CalibrationRegistryError:
        return None, path.stem
    return identity.shot_number, path.stem


def _matches_shot(path: Path, shots: set[str]) -> bool:
    """Match a whole shot token, never a substring: 42 is not part of 142."""

    number, stem = _shot_identity(path)
    for shot in shots:
        token = shot.strip()
        if not token:
            continue
        if number is not None and token.isdigit() and int(token) == number:
            return True
        if token == stem:
            return True
    return False


def select_sample_paths(
    paths: list[str | Path], *, every: int = 1, shots: set[str] | None = None
) -> list[Path]:
    """Select interval samples plus explicitly named shots."""

    if every < 1:
        raise DriftError("sample interval must be at least 1")
    ordered = sorted(Path(path) for path in paths)
    chosen = {path for index, path in enumerate(ordered) if index % every == 0}
    if shots:
        chosen.update(path for path in ordered if _matches_shot(path, shots))
    return sorted(chosen)


def _iso_date(value: str) -> date | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        return date.fromisoformat(text[:10])
    except ValueError:
        pass
    match = _DATE_PATTERN.search(text.replace("\\", "/"))
    if match is None:
        return None
    try:
        return date(int(match["year"]), int(match["month"]), int(match["day"]))
    except ValueError:
        return None


def cube_date(attrs: dict[str, Any]) -> tuple[date | None, str]:
    """Return one cube's acquisition date and the attribute it came from."""

    for key in DATE_ATTRIBUTES:
        parsed = _iso_date(str(attrs.get(key, "")))
        if parsed is not None:
            return parsed, key
    return None, ""


def _filter_by_date(
    paths: list[Path], *, date_from: str | None, date_to: str | None
) -> list[Path]:
    """Keep cubes whose acquisition date lies inside an inclusive interval."""

    import xarray as xr

    try:
        low = date.fromisoformat(date_from) if date_from else None
        high = date.fromisoformat(date_to) if date_to else None
    except ValueError as exc:
        raise DriftError(f"--from/--to take ISO YYYY-MM-DD dates: {exc}") from exc
    if low is not None and high is not None and low > high:
        raise DriftError("--from must not be after --to")
    selected = []
    for path in paths:
        with xr.open_dataset(path) as opened:
            attrs = dict(opened.attrs)
        when, _source = cube_date(attrs)
        if when is None:
            raise DriftError(
                f"{path.name} carries no acquisition date in "
                + ", ".join(DATE_ATTRIBUTES)
                + "; date selection would silently drop it, so audit it by shot instead"
            )
        if (low is None or when >= low) and (high is None or when <= high):
            selected.append(path)
    if not selected:
        raise DriftError(
            "no sampled cube falls between "
            f"{date_from or '-inf'} and {date_to or '+inf'}"
        )
    return selected


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


def select_plasma_frames(intensity_2d: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    """Return the plasma-bright frame indices and the evidence for the choice."""

    totals = np.nansum(intensity_2d, axis=1)
    if totals.size == 1:
        return np.array([0]), {
            "rule": "single frame",
            "frames": 1,
            "selected_frames": 1,
        }
    background = float(np.median(totals))
    deviation = 1.4826 * float(np.median(np.abs(totals - background)))
    if not np.isfinite(deviation) or deviation <= 0.0:
        deviation = float(np.std(totals))
    threshold = background + PLASMA_FRAME_SIGMA * deviation
    selected = np.flatnonzero(np.isfinite(totals) & (totals > threshold))
    return selected, {
        "rule": f"frame total > background + {PLASMA_FRAME_SIGMA:g} sigma",
        "frames": int(totals.size),
        "selected_frames": int(selected.size),
        "background_total": background,
        "sigma_total": float(deviation),
        "threshold_total": float(threshold),
    }


def _spectrum(ds) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Return the wavelength axis and the median over plasma-bright frames."""

    wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
    intensity = ds["intensity"]
    reduce_dims = [dim for dim in intensity.dims if dim != "wavelength"]
    if not reduce_dims:
        return (
            wavelength,
            np.asarray(intensity.values, dtype=float),
            {"rule": "single frame", "frames": 1, "selected_frames": 1},
        )
    stacked = np.asarray(
        intensity.transpose(*reduce_dims, "wavelength").values, dtype=float
    ).reshape(-1, wavelength.size)
    frames, selection = select_plasma_frames(stacked)
    if frames.size == 0:
        return wavelength, np.full(wavelength.size, np.nan), selection
    return (
        wavelength,
        np.nanmedian(stacked[frames], axis=0),
        selection,
    )


def centroid_evidence(
    wavelength: np.ndarray,
    intensity: np.ndarray,
    *,
    expected_nm: float,
    window_nm: float = 0.4,
) -> dict[str, Any]:
    """Measure a baseline-subtracted centroid with explicit sufficiency evidence."""

    selected = np.abs(wavelength - expected_nm) <= window_nm
    x = wavelength[selected]
    y = intensity[selected]
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 5:
        return {"status": "insufficient-data", "reason": "fewer than five samples in window"}
    edge_count = max(1, x.size // 5)
    baseline = float(np.median(np.r_[y[:edge_count], y[-edge_count:]]))
    signal = np.clip(y - baseline, 0.0, None)
    noise = float(np.std(np.r_[y[:edge_count], y[-edge_count:]] - baseline))
    peak = float(np.max(signal))
    snr = peak / max(noise, np.finfo(float).eps)
    if not np.isfinite(snr) or snr < MINIMUM_SNR or signal.sum() <= 0:
        return {
            "status": "insufficient-data",
            "reason": f"peak SNR {snr:.3g} is below {MINIMUM_SNR:g}",
            "snr": snr,
        }
    centroid = float(np.sum(x * signal) / np.sum(signal))
    return {
        "status": "measured",
        "centroid_nm": centroid,
        "residual_nm": centroid - expected_nm,
        "snr": snr,
        "window_nm": window_nm,
    }


class CubeGeometry:
    """One cube's own detector geometry: per-order polynomials and coverage.

    Packet 8 stores the wavelength polynomials, the pre-flip ``detector_pixel``
    and the ``echelle_order`` of every sample precisely so a later audit can ask
    what a wavelength residual is worth in pixels without reopening raw data.
    """

    def __init__(self, ds):
        raw = ds.attrs.get("wavelength_polynomials_json")
        if not raw:
            raise DriftError(
                "cube carries no wavelength_polynomials_json; the drift judge measures "
                "detector pixels, so re-export it through a calibration snapshot"
            )
        try:
            payload = json.loads(str(raw))
            self.coefficients = {
                int(record["order"]): [float(value) for value in record["coefficients"]]
                for record in payload["orders"]
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise DriftError(f"cube has malformed wavelength_polynomials_json: {exc}") from exc
        missing = [name for name in ("detector_pixel", "echelle_order") if name not in ds.coords]
        if missing:
            raise DriftError(
                "cube lacks the aligned " + ", ".join(missing) + " coordinate(s) the "
                "detector-space judge needs"
            )
        self.wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
        self.pixel = np.asarray(ds.coords["detector_pixel"].values, dtype=float)
        self.order = np.asarray(ds.coords["echelle_order"].values).astype(int)
        self.pixel_range = {
            int(order): (
                float(np.min(self.pixel[self.order == order])),
                float(np.max(self.pixel[self.order == order])),
            )
            for order in np.unique(self.order)
        }

    @property
    def orders(self) -> tuple[int, ...]:
        return tuple(sorted(self.pixel_range))

    def sample_at(self, wavelength_nm: float) -> tuple[int, float]:
        """Return the ``(order, detector_pixel)`` of the nearest stored sample."""

        index = int(np.argmin(np.abs(self.wavelength - float(wavelength_nm))))
        return int(self.order[index]), float(self.pixel[index])

    def _clamped(self, order: int, pixel: float) -> tuple[float, bool]:
        low, high = self.pixel_range[int(order)]
        clamped = not (low <= pixel <= high)
        return float(min(max(pixel, low), high)), clamped

    def dispersion(self, order: int, pixel: float) -> tuple[float, bool]:
        """Return d(lambda)/dx in nm/px, clamped to the order's covered pixels."""

        if int(order) not in self.coefficients:
            raise DriftError(f"cube has no wavelength polynomial for order {order}")
        clamped_pixel, clamped = self._clamped(order, pixel)
        derivative = np.polyder(np.asarray(self.coefficients[int(order)], dtype=float))
        return float(np.polyval(derivative, clamped_pixel)), clamped

    def wavelength_at(self, order: int, pixel: float) -> tuple[float, bool]:
        """Return the order's wavelength at one pixel, clamped to its coverage."""

        clamped_pixel, clamped = self._clamped(order, pixel)
        return (
            float(np.polyval(np.asarray(self.coefficients[int(order)], dtype=float), clamped_pixel)),
            clamped,
        )


def _alignment_tolerance_px(accuracy_nm: float | None, dispersion_nm_per_px: float) -> float:
    """Return the alignment bound for one line, in detector pixels."""

    tolerance = ALIGNMENT_MAX_PIXEL_RESIDUAL
    if accuracy_nm is not None and np.isfinite(accuracy_nm) and accuracy_nm > 0:
        derived = float(accuracy_nm) / abs(dispersion_nm_per_px)
        tolerance = max(tolerance, derived)
    return float(tolerance)


def _in_detector_space(
    geometry: CubeGeometry, measured: dict[str, Any], *, accuracy_nm: float | None
) -> dict[str, Any]:
    """Turn one measured wavelength residual into a detector-pixel residual."""

    order, pixel = geometry.sample_at(measured["centroid_nm"])
    try:
        dispersion, clamped = geometry.dispersion(order, pixel)
    except DriftError as exc:
        return {"status": "insufficient-data", "reason": str(exc)}
    if not np.isfinite(dispersion) or dispersion == 0.0:
        return {
            "status": "insufficient-data",
            "reason": f"order {order} has no usable dispersion at pixel {pixel:.1f}",
        }
    tolerance_px = _alignment_tolerance_px(accuracy_nm, dispersion)
    result = {
        **measured,
        "echelle_order": int(order),
        "detector_pixel": float(pixel),
        "dispersion_nm_per_px": float(dispersion),
        "pixel_residual_px": float(measured["residual_nm"] / dispersion),
        "alignment_tolerance_px": tolerance_px,
        "alignment_tolerance_nm": float(tolerance_px * abs(dispersion)),
    }
    if clamped:
        low, high = geometry.pixel_range[int(order)]
        result["extrapolated_dispersion"] = True
        result["extrapolation_note"] = (
            f"dispersion evaluated at order {order}'s covered pixel range "
            f"[{low:.1f}, {high:.1f}] instead of extrapolating the polynomial"
        )
    return result


def _cube_accuracy_nm(attrs: dict[str, Any]) -> float | None:
    raw = attrs.get("wavelength_accuracy_nm")
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


def _resolved(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        item
        for item in lines
        if item.get("status") == "measured"
        and not item.get("blended")
        and np.isfinite(float(item.get("pixel_residual_px", np.nan)))
    ]


def _quorum(
    lines: list[dict[str, Any]], coverage_nm: tuple[float, float] | None
) -> dict[str, Any]:
    """Report whether the measured lines can carry a verdict at all."""

    measured = _resolved(lines)
    orders = sorted({int(item["echelle_order"]) for item in measured})
    if coverage_nm is None:
        expected = [float(item["expected_nm"]) for item in lines if "expected_nm" in item]
        coverage_nm = (min(expected), max(expected)) if expected else (0.0, 0.0)
    midpoint = 0.5 * (float(coverage_nm[0]) + float(coverage_nm[1]))
    blue = any(float(item["expected_nm"]) < midpoint for item in measured)
    red = any(float(item["expected_nm"]) > midpoint for item in measured)
    reasons = []
    if len(measured) < MINIMUM_CENTROIDS:
        reasons.append(f"{len(measured)} resolved line(s) below the required {MINIMUM_CENTROIDS}")
    if len(orders) < MINIMUM_DISTINCT_ORDERS:
        reasons.append(
            f"{len(orders)} distinct order(s) below the required {MINIMUM_DISTINCT_ORDERS}"
        )
    if not (blue and red):
        reasons.append(
            "measured lines do not span both halves of the audited coverage "
            f"({coverage_nm[0]:.3f}-{coverage_nm[1]:.3f} nm)"
        )
    return {
        "resolved_lines": len(measured),
        "distinct_orders": len(orders),
        "orders": orders,
        "coverage_nm": [float(coverage_nm[0]), float(coverage_nm[1])],
        "bluer_half": bool(blue),
        "redder_half": bool(red),
        "satisfied": not reasons,
        "reason": "; ".join(reasons),
    }


def verdict_from_evidence(
    evidence: list[dict[str, Any]],
    *,
    coverage_nm: tuple[float, float] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Classify one interval in detector space and report how it was classified."""

    quorum = _quorum(evidence, coverage_nm)
    if not quorum["satisfied"]:
        return "insufficient-data", {"quorum": quorum}
    measured = _resolved(evidence)
    shifts = np.asarray([float(item["pixel_residual_px"]) for item in measured], dtype=float)
    dispersions = np.asarray(
        [abs(float(item["dispersion_nm_per_px"])) for item in measured], dtype=float
    )
    median = float(np.median(shifts))
    maximum = float(np.max(np.abs(shifts)))
    spread = float(np.max(np.abs(shifts - median)))
    tolerances = np.asarray(
        [
            float(item.get("alignment_tolerance_px", ALIGNMENT_MAX_PIXEL_RESIDUAL))
            for item in measured
        ],
        dtype=float,
    )
    summary = {
        "median_shift_px": median,
        "maximum_absolute_pixel_residual_px": maximum,
        "maximum_pixel_deviation_px": spread,
        "alignment_tolerance_px": float(np.max(tolerances)),
        "shift_consistency_px": SHIFT_CONSISTENCY_PX,
        # Reporting only: one rigid pixel shift is worth different wavelengths
        # in different orders, which is the whole reason the judge left nm.
        "median_shift_nm_equivalents": [
            float(median * value) for value in (float(np.min(dispersions)), float(np.max(dispersions)))
        ],
        "quorum": quorum,
    }
    if np.all(np.abs(shifts) <= tolerances):
        return "aligned", summary
    if abs(median) <= REPAIR_LIMIT_PX and spread <= SHIFT_CONSISTENCY_PX:
        return "shifted", summary
    return "misaligned-beyond-repair", summary


def _per_shot_summary(lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group the resolved residuals per shot, ordered by date then shot."""

    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for item in _resolved(lines):
        key = (str(item.get("date", "")), str(item.get("shot_number", "")), str(item["cube"]))
        groups.setdefault(key, []).append(item)
    summary = []
    for (when, shot, cube), items in sorted(groups.items()):
        shifts = np.asarray([float(item["pixel_residual_px"]) for item in items], dtype=float)
        median = float(np.median(shifts))
        summary.append(
            {
                "shot_number": shot,
                "cube": cube,
                "date": when,
                "lines": len(items),
                "median_shift_px": median,
                "pixel_spread_px": float(np.max(np.abs(shifts - median))),
            }
        )
    return summary


def _split_boundary(per_shot: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], int] | None:
    """Return the shift-ordered shots and the widest gap when it is a real split."""

    if len(per_shot) < 2:
        return None
    ordered = sorted(per_shot, key=lambda item: item["median_shift_px"])
    medians = [item["median_shift_px"] for item in ordered]
    gaps = [medians[index + 1] - medians[index] for index in range(len(medians) - 1)]
    widest = int(np.argmax(gaps))
    if gaps[widest] <= SHIFT_CONSISTENCY_PX:
        return None
    return ordered, widest


def _shot_sort_key(shot: str) -> tuple[int, float, str]:
    return (0, float(shot), shot) if shot.isdigit() else (1, 0.0, shot)


def detect_interval_split(per_shot: list[dict[str, Any]]) -> str:
    """Flag two separated groups of per-shot shifts instead of condemning them.

    An interval that straddles a calibration change holds two consistent rigid
    shifts, not one inconsistent set of residuals.  When the per-shot medians
    fall into two clusters further apart than the consistency bound, the
    interval wants splitting at that boundary; the four verdict words still
    classify the interval, and this warning travels alongside them.
    """

    split = _split_boundary(per_shot)
    if split is None:
        return ""
    ordered, widest = split
    low = sorted((item["shot_number"] for item in ordered[: widest + 1]), key=_shot_sort_key)
    high = sorted((item["shot_number"] for item in ordered[widest + 1 :]), key=_shot_sort_key)
    return (
        f"residuals form two groups (shots {low[0]}..{low[-1]} vs "
        f"{high[0]}..{high[-1]}); consider splitting the interval at the boundary"
    )


def _label_shot_groups(per_shot: list[dict[str, Any]]) -> None:
    """Mark which cluster each shot belongs to when the interval looks split."""

    split = _split_boundary(per_shot)
    if split is None:
        return
    ordered, widest = split
    for position, item in enumerate(ordered):
        item["group"] = 1 if position <= widest else 2


# ---------------------------------------------------------------------------
# Repair composition
# ---------------------------------------------------------------------------


def _next_refinement_id(base_snapshot_id: str, root: Path) -> str:
    stem = re.sub(r"-r\d+$", "", str(base_snapshot_id))
    revision = 1
    while (root / f"{stem}-r{revision}").exists():
        revision += 1
    return f"{stem}-r{revision}"


def _drive_requirement(
    shots: list[str], catalog_path: str | Path | None
) -> str:
    """Name the drives holding the affected shots, or the flag that would."""

    if catalog_path is None:
        return (
            "recalibrating beyond repair needs the raw SIF data: pass "
            "--catalog <merged-catalog.json> to name the drives holding these shots"
        )
    from .catalog import load_catalog

    payload = load_catalog(catalog_path)
    sources = payload.get("sources") or [
        {"volume_label": payload.get("volume_label", "unknown"), "cubes": payload.get("cubes", [])}
    ]
    wanted = {str(shot) for shot in shots}
    drives = sorted(
        {
            str(source.get("volume_label", "unknown"))
            for source in sources
            if any(str(cube.get("shot_number", "")) in wanted for cube in source.get("cubes", []))
        }
    )
    if not drives:
        return (
            f"no catalog row matches the affected shot(s) {', '.join(sorted(wanted))}; "
            "rebuild the catalog or name another one with --catalog"
        )
    return (
        "recalibrating beyond repair needs the raw SIF data on drive(s) "
        + ", ".join(drives)
        + "; reprocess those shots against a re-measured calibration snapshot"
    )


def _repair_steps(
    *,
    evidence_path: Path,
    calibrations_root: Path,
    refined_id: str,
    refined_is_preview: bool,
    shift_px: float,
    cube_directories: list[Path],
) -> list[dict[str, str]]:
    """Compose the real, paste-ready repair sequence for a shifted interval."""

    snapshot_dir = calibrations_root / refined_id
    steps = [
        {
            "shell": "any",
            "purpose": "accept the sampled shift and emit the immutable -rN snapshot",
            # Paths are quoted so a campaign folder with a space survives both
            # shells, and the shift carries 12 significant digits: the
            # acknowledgement is exact, so a rounded command would refuse itself.
            "command": (
                f'echelle drift refine "{evidence_path}" '
                f'--calibrations "{calibrations_root}" --accept-shift {shift_px:.12g}'
            ),
        },
        {
            "shell": "any",
            "purpose": (
                f"repoint the registry [validity] entry at {refined_id}"
                + (
                    " (the id 'echelle drift refine' prints)"
                    if refined_is_preview
                    else ""
                )
                + " before any bulk run; the accepted shift authorizes only the refinement"
            ),
            "command": "",
        },
    ]
    # ``recal-cube`` revises one cube at a time, so the loop is written out in
    # both shells.  The refinement changed only the wavelength role, so
    # --wavelength-only is what it authorizes: the sphere pair and the integral
    # curve are the base snapshot's own bytes.
    for directory in cube_directories:
        purpose = f"revise every cube already exported from {directory}"
        steps.append(
            {
                "shell": "powershell",
                "purpose": purpose,
                "command": (
                    f'Get-ChildItem "{directory}\\*.nc" | ForEach-Object {{ '
                    f'echelle recal-cube $_.FullName -o ($_.FullName -replace "\\.nc$", '
                    f'"-{refined_id}.nc") --new-snapshot "{snapshot_dir}" --wavelength-only }}'
                ),
            }
        )
        steps.append(
            {
                "shell": "posix",
                "purpose": purpose,
                "command": (
                    f'for cube in "{directory.as_posix()}"/*.nc; do echelle recal-cube '
                    f'"$cube" -o "${{cube%.nc}}-{refined_id}.nc" '
                    f'--new-snapshot "{snapshot_dir.as_posix()}" --wavelength-only; done'
                ),
            }
        )
    return steps


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------


def audit_cubes(  # noqa: C901 - one readable pass over cubes, lines, and shots
    cube_paths: list[str | Path],
    *,
    every: int = 1,
    shots: set[str] | None = None,
    families: tuple[str, ...] = ("balmer", "fulcher"),
    date_from: str | None = None,
    date_to: str | None = None,
    catalog: str | Path | None = None,
    evidence_path: str | Path | None = None,
    calibrations_root: str | Path = "calibrations",
) -> dict[str, Any]:
    """Audit sampled saved cubes and return an immutable JSON-ready verdict."""

    import xarray as xr

    resolved = resolve_cube_paths(cube_paths)
    if date_from or date_to:
        resolved = _filter_by_date(resolved, date_from=date_from, date_to=date_to)
    selected = select_sample_paths(resolved, every=every, shots=shots)

    per_line: list[dict[str, Any]] = []
    snapshot_ids: set[str] = set()
    sampled: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    coverage: list[float] = []
    geometries: dict[str, CubeGeometry] = {}
    for path in selected:
        with xr.open_dataset(path) as opened:
            ds = opened.load()
        attrs = dict(ds.attrs)
        snapshot_ids.add(str(attrs.get("snapshot_id", "unassigned")))
        shot_number = str(attrs.get("shot_number", path.stem))
        when, date_source = cube_date(attrs)
        wavelength, intensity, frames = _spectrum(ds)
        accuracy_nm = _cube_accuracy_nm(attrs)
        record = {
            "cube": path.name,
            "sha256": sha256_file(path),
            "shot_number": shot_number,
            "date": when.isoformat() if when else "",
            "date_attribute": date_source,
            "frame_selection": frames,
            "wavelength_accuracy_nm": accuracy_nm,
        }
        sampled.append(record)
        common = {
            "cube": path.name,
            "shot_number": shot_number,
            "date": record["date"],
        }
        if frames["selected_frames"] == 0:
            reason = "no plasma-bright frames"
            skipped.append({**common, "reason": reason})
            per_line.append({**common, "status": "insufficient-data", "reason": reason})
            continue
        try:
            geometry = CubeGeometry(ds)
        except DriftError as exc:
            skipped.append({**common, "reason": str(exc)})
            per_line.append({**common, "status": "insufficient-data", "reason": str(exc)})
            continue
        geometries[path.name] = geometry
        minimum, maximum = float(np.min(wavelength)), float(np.max(wavelength))
        coverage.extend((minimum, maximum))
        for family in families:
            for line in load_line_table(family):
                if not minimum <= line.wavelength_nm <= maximum:
                    continue
                entry = {
                    **common,
                    "family": family,
                    "line": line.label,
                    "expected_nm": line.wavelength_nm,
                    "source_reference": line.source_reference,
                    "blended": bool(line.blended),
                }
                if line.blended:
                    per_line.append(
                        {
                            **entry,
                            "status": "skipped",
                            "reason": (
                                "sub-resolution blend: the measured centroid belongs to the "
                                "blend, not to this transition"
                            ),
                        }
                    )
                    continue
                result = centroid_evidence(wavelength, intensity, expected_nm=line.wavelength_nm)
                if result.get("status") == "measured":
                    result = _in_detector_space(geometry, result, accuracy_nm=accuracy_nm)
                per_line.append({**entry, **result})

    coverage_nm = (min(coverage), max(coverage)) if coverage else None
    verdict, summary = verdict_from_evidence(per_line, coverage_nm=coverage_nm)
    per_shot = _per_shot_summary(per_line)
    interval_warning = detect_interval_split(per_shot)
    _label_shot_groups(per_shot)
    summary["sampled_cubes"] = len(sampled)
    summary["skipped_cubes"] = len(skipped)
    summary["blended_lines_skipped"] = sum(
        1 for item in per_line if item.get("status") == "skipped" and item.get("blended")
    )

    order_corrections: list[dict[str, Any]] = []
    repair_steps: list[dict[str, str]] = []
    repair = ""
    refined_preview = ""
    if verdict == "shifted":
        shift_px = float(summary["median_shift_px"])
        order_corrections = _order_corrections(geometries, shift_px)
        root = Path(calibrations_root)
        base_ids = sorted(snapshot_ids)
        refined_is_preview = True
        if len(base_ids) == 1 and (root / base_ids[0]).is_dir():
            refined_preview = _next_refinement_id(base_ids[0], root)
            refined_is_preview = False
        else:
            refined_preview = f"{base_ids[0] if base_ids else 'snapshot'}-rN"
        evidence_target = Path(evidence_path) if evidence_path else Path("drift-evidence.json")
        repair_steps = _repair_steps(
            evidence_path=evidence_target,
            calibrations_root=root,
            refined_id=refined_preview,
            refined_is_preview=refined_is_preview,
            shift_px=shift_px,
            cube_directories=sorted({path.parent for path in selected}),
        )
        repair = repair_steps[0]["command"]

    payload: dict[str, Any] = {
        "schema": DRIFT_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample_rule": {
            "every": every,
            "shots": sorted(shots or set()),
            "date_from": date_from or "",
            "date_to": date_to or "",
            "date_attributes": list(DATE_ATTRIBUTES),
        },
        "snapshot_ids": sorted(snapshot_ids),
        "sampled_cubes": sampled,
        "skipped_cubes": skipped,
        "thresholds_px": {
            "alignment_maximum_residual": ALIGNMENT_MAX_PIXEL_RESIDUAL,
            "shift_consistency": SHIFT_CONSISTENCY_PX,
            "repair_limit": REPAIR_LIMIT_PX,
            "minimum_centroids": MINIMUM_CENTROIDS,
            "minimum_distinct_orders": MINIMUM_DISTINCT_ORDERS,
            "minimum_snr": MINIMUM_SNR,
            "plasma_frame_sigma": PLASMA_FRAME_SIGMA,
        },
        "verdict": verdict,
        "summary": summary,
        "per_shot": per_shot,
        "lines": per_line,
        "order_corrections": order_corrections,
        "repair_command": repair,
        "repair_commands": repair_steps,
    }
    if refined_preview:
        payload["refined_snapshot_id"] = refined_preview
    if interval_warning:
        payload["interval_warning"] = interval_warning
    if verdict == "misaligned-beyond-repair":
        payload["data_requirement"] = _drive_requirement(
            sorted({str(item["shot_number"]) for item in per_shot}), catalog
        )
    return payload


def _order_corrections(
    geometries: dict[str, CubeGeometry], shift_px: float
) -> list[dict[str, Any]]:
    """Predict the wavelength correction the refinement will apply per order.

    Moving the table's anchors by ``shift_px`` makes the refitted solution
    ``p_B(x) = p_A(x - shift_px)``, so this is exactly the change
    ``recal-cube`` re-evaluates from the refinement's own polynomials.
    """

    merged: dict[int, tuple[CubeGeometry, tuple[float, float]]] = {}
    for geometry in geometries.values():
        for order, bounds in geometry.pixel_range.items():
            if order not in merged:
                merged[order] = (geometry, bounds)
    corrections = []
    for order in sorted(merged):
        geometry, (low, high) = merged[order]
        reference = 0.5 * (low + high)
        before, _ = geometry.wavelength_at(order, reference)
        after, clamped = geometry.wavelength_at(order, reference - shift_px)
        dispersion, _ = geometry.dispersion(order, reference)
        corrections.append(
            {
                "order": int(order),
                "reference_pixel": float(reference),
                "pixel_range": [float(low), float(high)],
                "dispersion_nm_per_px": float(dispersion),
                "predicted_shift_nm": float(after - before),
                **({"extrapolated": True} if clamped else {}),
            }
        )
    return corrections


def write_drift_evidence(path: str | Path, payload: dict[str, Any]) -> Path:
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(f"drift evidence is immutable and already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return destination


# ---------------------------------------------------------------------------
# Refinement
# ---------------------------------------------------------------------------


def _shift_wavelength_table(
    source: Path,
    destination: Path,
    *,
    shift_px: float,
    metadata: list[tuple[str, str]],
) -> int:
    """Slide every anchor along the detector, keeping its wavelength.

    The anchor that produced wavelength ``lambda`` at pixel ``c`` is now seen at
    ``c + shift_px``, so the refit solution assigns ``lambda`` there and each
    order's own dispersion turns the one pixel shift into that order's
    wavelength correction.
    """

    rows = load_wavelength_table(source)
    if not rows:
        raise DriftError(f"wavelength table {source.name} holds no correctable rows")
    write_wavelength_table(shift_lines_in_pixels(rows, shift_px), destination, metadata=metadata)
    return len(rows)


def create_refinement_snapshot(
    evidence_path: str | Path,
    *,
    calibrations_root: str | Path,
    accepted_shift_px: float,
) -> tuple[Snapshot, Path]:
    """Accept a shifted verdict and emit a new immutable ``-rN`` snapshot."""

    evidence_file = Path(evidence_path)
    accepted_path = evidence_file.with_suffix(".accepted.json")
    if accepted_path.exists():
        raise FileExistsError(f"accepted drift record already exists: {accepted_path}")
    evidence = json.loads(evidence_file.read_text(encoding="utf-8"))
    if evidence.get("schema") == DRIFT_SCHEMA_V1:
        raise DriftError(
            "this evidence states a constant wavelength offset, which no longer describes a "
            "rigid detector shift; re-run 'echelle drift audit' to measure it in pixels"
        )
    if evidence.get("schema") != DRIFT_SCHEMA or evidence.get("verdict") != "shifted":
        raise DriftError("only a shifted sampled verdict can create a refinement snapshot")
    measured = float(evidence.get("summary", {}).get("median_shift_px", np.nan))
    # Exact acknowledgement, with room only for the decimal text the composed
    # command carries: 1e-9 relative on a pixel is a billionth of a pixel.
    if not np.isfinite(measured) or not np.isclose(
        accepted_shift_px, measured, rtol=1e-9, atol=1e-9
    ):
        raise DriftError("--accept-shift must exactly acknowledge the sampled median shift")
    snapshot_ids = list(evidence.get("snapshot_ids", []))
    if len(snapshot_ids) != 1:
        raise DriftError("refinement evidence must refer to exactly one base snapshot")
    root = Path(calibrations_root)
    base = load_snapshot(root / snapshot_ids[0])
    refined_id = _next_refinement_id(base.snapshot_id, root)
    with tempfile.TemporaryDirectory(prefix="echelle-refinement-") as temporary:
        adjusted = Path(temporary) / "wavelength.txt"
        moved = _shift_wavelength_table(
            base.root / base.artifact_for_role("wavelength").path,
            adjusted,
            shift_px=measured,
            metadata=[
                ("base snapshot", base.snapshot_id),
                ("evidence", evidence_file.name),
                ("evidence sha256", sha256_file(evidence_file)),
                ("model", "rigid detector translation; anchors moved, wavelengths unchanged"),
                ("applied pixel shift", f"{measured:+.4f} px"),
            ],
        )
        files = {
            role: (
                adjusted
                if role == "wavelength"
                else base.root / base.artifact_for_role(role).path
            )
            for role in REQUIRED_ROLES
        }
        lamps = [
            (artifact.label, base.root / artifact.path)
            for artifact in base.artifacts
            if artifact.role == "lamp"
        ]
        snapshot = create_snapshot(
            root,
            snapshot_id=refined_id,
            detector=base.detector,
            files=files,
            lamps=base.lamps,
            lamp_files=lamps,
            notes=f"Accepted science-line refinement from {evidence_file.name}.",
            base_snapshot=base.snapshot_id,
            validity=base.manifest.get("validity"),
            alignment={
                "method": "sampled Balmer/Fulcher rigid detector shift",
                "observed_shift_px": measured,
                "applied_pixel_shift_px": measured,
                "corrected_rows": moved,
                "evidence_sha256": sha256_file(evidence_file),
            },
            qc={"sampled_verdict": "shifted", "accepted": True},
        )
    accepted = {
        **evidence,
        "accepted": True,
        "accepted_snapshot_id": snapshot.snapshot_id,
        # Accepting a shift condemns the audited snapshot: only the refinement
        # carries the corrected wavelength table, so only the refinement is
        # authorized for bulk processing.
        "snapshot_ids": [snapshot.snapshot_id],
        "base_snapshot_ids": sorted(snapshot_ids),
        "base_evidence_sha256": sha256_file(evidence_file),
    }
    write_drift_evidence(accepted_path, accepted)
    return snapshot, accepted_path


# ---------------------------------------------------------------------------
# Gate
# ---------------------------------------------------------------------------


def _require_authorized_snapshots(selected: set[str], authorized: set[str]) -> None:
    if not selected <= authorized:
        missing = ", ".join(sorted(selected - authorized))
        raise DriftError(f"drift evidence does not cover selected snapshot(s): {missing}")


def require_sampled_verdict(path: str | Path, snapshot_ids: set[str]) -> dict[str, Any]:
    """Validate the bulk-processing prerequisite without treating insufficiency as aligned."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") not in READABLE_DRIFT_SCHEMAS:
        raise DriftError("unsupported drift evidence schema")
    verdict = payload.get("verdict")
    if verdict == "aligned":
        _require_authorized_snapshots(
            snapshot_ids, set(map(str, payload.get("snapshot_ids", [])))
        )
        return payload
    if verdict == "shifted" and payload.get("accepted"):
        refined = str(payload.get("accepted_snapshot_id") or "")
        recorded = set(map(str, payload.get("snapshot_ids", [])))
        # Records written before the refinement narrowed its own authorization
        # list both ids; everything that is not the refinement is superseded.
        superseded = set(map(str, payload.get("base_snapshot_ids") or [])) or (
            recorded - {refined}
        )
        condemned = sorted(snapshot_ids & superseded)
        if condemned:
            raise DriftError(
                f"registry still selects {', '.join(condemned)}; the accepted correction "
                f"produced {refined} — update the registry validity to point at the "
                "refined snapshot"
            )
        _require_authorized_snapshots(snapshot_ids, {refined} if refined else set())
        return payload
    if verdict == "insufficient-data":
        raise DriftError("bulk processing refused: sampled verdict is insufficient-data")
    raise DriftError(f"bulk processing refused: sampled verdict is {verdict}")
