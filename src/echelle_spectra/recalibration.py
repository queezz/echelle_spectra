"""Post-hoc wavelength and absolute-factor recalibration of saved cubes."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .campaign_run import sha256_file
from .snapshot import Snapshot, load_snapshot

RECALIBRATION_SCHEMA = "echelle-cube-recalibration/v1"
POLYNOMIAL_SCHEMA = "spectrocube.wavelength-polynomials/v1"


class RecalibrationError(ValueError):
    """A saved cube cannot be revised safely under the requested snapshots."""


def _artifact_digest(records: object, role: str) -> str:
    if isinstance(records, str):
        records = json.loads(records)
    if isinstance(records, dict):
        item = records.get(role, {})
        return str(item.get("sha256", "")) if isinstance(item, dict) else ""
    if isinstance(records, list):
        for item in records:
            if isinstance(item, dict) and item.get("role") == role:
                return str(item.get("sha256", ""))
    return ""


def _new_snapshot_digest(snapshot: Snapshot, role: str) -> str:
    return snapshot.artifact_for_role(role).sha256


def _validated_old_records(ds) -> object:
    try:
        manifest = json.loads(str(ds.attrs["snapshot_manifest_json"]))
        records = json.loads(str(ds.attrs["calibration_file_digests_json"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise RecalibrationError("cube has malformed Packet 8 snapshot provenance") from exc
    if str(manifest.get("id", "")) != str(ds.attrs.get("snapshot_id", "")):
        raise RecalibrationError("cube snapshot_id disagrees with its embedded snapshot manifest")
    manifest_pattern = _artifact_digest(manifest.get("artifacts", []), "pattern")
    records_pattern = _artifact_digest(records, "pattern")
    if not manifest_pattern or manifest_pattern != records_pattern:
        raise RecalibrationError("cube pattern digest disagrees across embedded provenance records")
    return records


def _require_unchanged_unapplied_roles(
    old_records: object,
    new_snapshot: Snapshot,
    *,
    update_wavelength: bool,
    update_factor: bool,
) -> None:
    unchanged_roles = []
    if not update_wavelength:
        unchanged_roles.append("wavelength")
    if not update_factor:
        unchanged_roles.extend(("sphere", "sphere_background", "integral"))
    changed = [
        role
        for role in unchanged_roles
        if _artifact_digest(old_records, role) != _new_snapshot_digest(new_snapshot, role)
    ]
    if changed:
        raise RecalibrationError(
            "new snapshot changes unselected calibration role(s) "
            + ", ".join(changed)
            + "; a full recalibration (neither --wavelength-only nor --factor-only) now "
            "applies a different era's wavelength solution and sphere response together "
            "on the cube's own detector grid, so run that instead -- or choose a snapshot "
            "whose other roles are unchanged"
        )


def _wavelength_solution(snapshot: Snapshot) -> dict[int, list[float]]:
    table = snapshot.source_path("wavelength")
    try:
        rows = np.loadtxt(table, comments="#", usecols=(0, 3, 4), ndmin=2)
    except (OSError, ValueError) as exc:
        raise RecalibrationError(f"cannot read wavelength artifact {table}: {exc}") from exc
    solutions: dict[int, list[float]] = {}
    for raw_order in np.unique(rows[:, 0]).astype(int):
        selected = rows[rows[:, 0].astype(int) == raw_order]
        degree = 1 if len(selected) < 3 else 2
        solutions[int(raw_order)] = np.polyfit(selected[:, 1], selected[:, 2], degree).tolist()
    return solutions


def _polynomial_payload(snapshot: Snapshot, solutions: dict[int, list[float]]) -> dict[str, Any]:
    return {
        "schema": POLYNOMIAL_SCHEMA,
        "coefficient_order": "descending_power",
        "input": "detector_pixel",
        "input_units": "pixel",
        "output": "wavelength",
        "output_units": "nm",
        "orders": [
            {"order": order, "coefficients": solutions[order]} for order in sorted(solutions)
        ],
        "writer": "echelle_spectra",
        "snapshot_id": snapshot.snapshot_id,
        "wavelength_artifact_sha256": _new_snapshot_digest(snapshot, "wavelength"),
    }


def _evaluate_wavelength(
    ds, snapshot: Snapshot
) -> tuple[np.ndarray, dict[str, Any], dict[int, list[float]]]:
    if "detector_pixel" not in ds.coords or "echelle_order" not in ds.coords:
        raise RecalibrationError(
            "wavelength recalibration needs aligned detector_pixel and echelle_order coordinates"
        )
    solutions = _wavelength_solution(snapshot)
    pixels = np.asarray(ds.coords["detector_pixel"].values, dtype=float)
    orders = np.asarray(ds.coords["echelle_order"].values, dtype=int)
    missing = sorted(set(orders.tolist()) - set(solutions))
    if missing:
        raise RecalibrationError(
            "new snapshot has no wavelength polynomial for represented order(s): "
            + ", ".join(map(str, missing))
        )
    represented_solutions = {order: solutions[order] for order in sorted(set(orders.tolist()))}
    wavelength = np.empty_like(pixels, dtype=float)
    for order in np.unique(orders):
        selected = orders == order
        wavelength[selected] = np.polyval(solutions[int(order)], pixels[selected])
    return wavelength, _polynomial_payload(snapshot, represented_solutions), solutions


#: Absolute-factor kinds and the intensity units each one produces.  Ordered
#: longest unit text first so ``wmsr`` is never read as its ``wm`` prefix.
ABSOLUTE_KIND_UNITS = (
    ("phmsr", "ph/s/nm/sr"),
    ("wmsr", "W/m2/nm/sr"),
    ("wm", "W/m2/nm"),
)


def _factor_kind_from_units(units: str) -> str | None:
    """Recover the applied kind from a pre-F2 cube's factor units text."""

    lowered = units.casefold()
    for kind, unit_text in ABSOLUTE_KIND_UNITS:
        if lowered.startswith(unit_text.casefold()):
            return kind
    return None


def _factor_kind(factor) -> str:
    """Return which sphere curve produced a cube's applied absolute factor.

    Cubes written from this version onwards state the kind outright in the
    ``absolute_kind`` attribute.  The units text remains the documented fallback
    for cubes written before it existed, and an unrecognisable units string is
    refused rather than guessed: rescaling by the wrong sphere curve is silent.
    """

    known = {kind for kind, _ in ABSOLUTE_KIND_UNITS}
    declared = str(factor.attrs.get("absolute_kind", "")).strip().casefold()
    if declared:
        if declared not in known:
            raise RecalibrationError(
                f"cube declares an unknown absolute_kind {declared!r}; expected one of "
                + ", ".join(sorted(known))
            )
        return declared
    units = str(factor.attrs.get("units", "")).strip()
    kind = _factor_kind_from_units(units)
    if kind is None:
        raise RecalibrationError(
            "cube does not state which absolute calibration curve it applied: "
            "applied_absolute_calibration_factor carries no absolute_kind attribute and "
            f"its units {units!r} match none of "
            + ", ".join(f"{name} ({unit})" for name, unit in ABSOLUTE_KIND_UNITS)
            + "; re-export the cube so it records absolute_kind"
        )
    return kind


def _cube_detector_grid(ds) -> tuple[np.ndarray, np.ndarray]:
    """Return the cube's own (echelle order, detector pixel) sample addresses."""

    if "detector_pixel" not in ds.coords or "echelle_order" not in ds.coords:
        raise RecalibrationError(
            "absolute-factor recalibration needs aligned detector_pixel and "
            "echelle_order coordinates"
        )
    pixels = np.asarray(ds.coords["detector_pixel"].values, dtype=float)
    orders = np.asarray(ds.coords["echelle_order"].values, dtype=float)
    if not np.all(np.isfinite(pixels)) or not np.all(np.isfinite(orders)):
        raise RecalibrationError("cube detector_pixel/echelle_order coordinates are not finite")
    return orders.astype(int), pixels.astype(int)


def factor_from_snapshot(ds, snapshot: Snapshot) -> np.ndarray:
    """Evaluate a snapshot's absolute factor on the cube's own detector grid.

    A cross-era revision is exactly the case this has to survive: the new
    snapshot carries its own wavelength table, that table moves the order
    borders, and the border move changes which (order, pixel) pairs the new
    snapshot publishes.  Matching the two snapshots' published sample sets is
    therefore the wrong contract -- it fails for precisely the era change the
    revision exists to perform.

    The right grid is the cube's.  Extraction geometry is frozen (the pattern
    digest is checked before anything here runs), so every cube sample names a
    detector column inside a diffraction order that both snapshots cut the same
    way, and the new snapshot's sphere pair measured that column whether or not
    its borders chose to publish it.  ``absolute_on_detector_grid`` re-evaluates
    the snapshot's own factor formula there, so a sample recovered across a
    moved border is a measurement, not an interpolation or an extrapolation of
    the response curve.

    Samples the new snapshot genuinely cannot answer for come back as NaN and
    are dropped and counted upstream: the partial-order pad that leaves the
    sensor, a column of exactly zero net sphere response, and any column outside
    the band the sphere's spectral reference covers.  The sphere is never
    invented where it did not shine.
    """

    from .spectrocube_cli import _snapshot_camera
    from .tools.loader import build_calibration

    orders, pixels = _cube_detector_grid(ds)
    calibration = build_calibration(
        snapshot.root,
        _snapshot_camera(snapshot),
        calibration_files=snapshot.calibration_files(),
    )
    kind = _factor_kind(ds["applied_absolute_calibration_factor"])
    grid = np.asarray(calibration.absolute_on_detector_grid()[kind], dtype=float)
    order_ids = np.asarray(calibration.order_ids, dtype=int)
    rows = {int(order_id): index for index, order_id in enumerate(order_ids)}
    missing = sorted(set(orders.tolist()) - set(rows))
    if missing:
        raise RecalibrationError(
            "new snapshot has no diffraction order for represented order(s): "
            + ", ".join(map(str, missing))
        )
    width = int(grid.shape[1])
    outside = (pixels < 0) | (pixels >= width)
    if np.any(outside):
        raise RecalibrationError(
            "cube detector pixel(s) fall outside the new snapshot's "
            f"{width}-column detector: "
            + ", ".join(map(str, sorted(set(pixels[outside].tolist()))[:8]))
        )
    return grid[np.asarray([rows[int(order)] for order in orders], dtype=int), pixels]


def _replace_factor(revised, original, new_snapshot: Snapshot, new_factor: np.ndarray | None):
    """Apply a new absolute factor, dropping the samples it cannot calibrate.

    Export drops wavelength columns whose absolute factor is not strictly
    positive — one noise-negative sphere-minus-background column must not fail a
    whole file — and recalibration mirrors that tolerance: a sample whose old or
    new factor is non-finite or non-positive is removed and counted in the
    ``dropped_nonpositive_factor_columns`` attribute instead of refusing the
    revision.  Every retained sample keeps the strictly positive guarantee.

    A cross-era revision adds one more way for a sample to be uncalibratable,
    and it deserves its own count.  A non-finite replacement means the new
    snapshot's sphere never answered for that detector column at all, which is a
    coverage statement about the calibration; a finite non-positive one means it
    answered with noise, which is a measurement.  Both are dropped, and the
    coverage half is reported separately as
    ``dropped_uncovered_factor_columns`` so an operator can tell an era whose
    sphere does not reach this cube from an era that merely measured badly.
    """

    if "applied_absolute_calibration_factor" not in revised:
        raise RecalibrationError(
            "absolute-factor recalibration needs applied_absolute_calibration_factor"
        )
    old = revised["applied_absolute_calibration_factor"]
    replacement = (
        np.asarray(new_factor, dtype=float)
        if new_factor is not None
        else factor_from_snapshot(revised, new_snapshot)
    )
    if replacement.shape != old.shape:
        raise RecalibrationError("replacement absolute factor must be wavelength-aligned")
    previous = np.asarray(old.values, dtype=float)
    covered = np.isfinite(replacement)
    usable = covered & (replacement > 0) & np.isfinite(previous) & (previous > 0)
    dropped = int(usable.size - np.count_nonzero(usable))
    uncovered = int(covered.size - np.count_nonzero(covered))
    if not np.any(usable):
        raise RecalibrationError(
            "no wavelength sample has a strictly positive old and replacement absolute factor; "
            f"{uncovered} of {covered.size} samples lie outside the new snapshot's sphere "
            "coverage, which usually means the wrong snapshot or the wrong detector"
        )
    if dropped:
        revised = revised.isel(wavelength=np.flatnonzero(usable))
        replacement = replacement[usable]
        old = revised["applied_absolute_calibration_factor"]
        revised.attrs["dropped_nonpositive_factor_columns"] = (
            int(revised.attrs.get("dropped_nonpositive_factor_columns", 0)) + dropped
        )
    if uncovered:
        revised.attrs["dropped_uncovered_factor_columns"] = (
            int(revised.attrs.get("dropped_uncovered_factor_columns", 0)) + uncovered
        )
    # Rescale by the ratio rather than by dividing out the old factor and
    # multiplying the new one back in.  An unchanged factor then gives a ratio of
    # exactly 1.0 and leaves every stored intensity bit for bit as it was, which
    # is what makes recalibrating a cube onto its own snapshot a true round trip.
    revised["intensity"] = revised["intensity"] * (replacement / old)
    revised["intensity"].attrs.update(original["intensity"].attrs)
    revised["applied_absolute_calibration_factor"] = (
        ("wavelength",),
        replacement,
        dict(old.attrs),
    )
    return revised, dropped, uncovered


def _append_history(
    ds,
    revised,
    new_snapshot: Snapshot,
    changes: list[str],
    *,
    dropped_nonpositive_factor_columns: int = 0,
    dropped_uncovered_factor_columns: int = 0,
) -> dict[str, Any]:
    history = []
    raw_history = ds.attrs.get("recalibration_history_json")
    if raw_history:
        try:
            history = list(json.loads(str(raw_history)))
        except (TypeError, ValueError):
            raise RecalibrationError("cube has malformed recalibration_history_json") from None
    event = {
        "schema": RECALIBRATION_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "old_snapshot_id": str(ds.attrs["snapshot_id"]),
        "old_snapshot_manifest_sha256": str(ds.attrs.get("snapshot_manifest_sha256", "")),
        "old_snapshot_manifest_json": str(ds.attrs["snapshot_manifest_json"]),
        "new_snapshot_id": new_snapshot.snapshot_id,
        "new_snapshot_manifest_sha256": new_snapshot.provenance_attrs()[
            "snapshot_manifest_sha256"
        ],
        "new_snapshot_manifest_json": new_snapshot.provenance_attrs()["snapshot_manifest_json"],
        "changes": changes,
    }
    if dropped_nonpositive_factor_columns:
        event["dropped_nonpositive_factor_columns"] = int(dropped_nonpositive_factor_columns)
    if dropped_uncovered_factor_columns:
        event["dropped_uncovered_factor_columns"] = int(dropped_uncovered_factor_columns)
    history.append(event)
    revised.attrs.update(new_snapshot.provenance_attrs())
    revised.attrs["recalibration_history_json"] = json.dumps(
        history, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    revised.attrs["previous_snapshot_manifest_json"] = str(ds.attrs["snapshot_manifest_json"])
    return event


def _json_records(attrs, name: str) -> list[dict[str, Any]] | None:
    raw = attrs.get(name)
    if not raw:
        return None
    try:
        records = json.loads(str(raw))
    except (TypeError, ValueError):
        raise RecalibrationError(f"cube has malformed {name}") from None
    if not isinstance(records, list):
        raise RecalibrationError(f"cube has malformed {name}")
    return records


def _refresh_wavelength_derived_attrs(revised, solutions: dict[int, list[float]]) -> None:
    """Recompute every attribute the superseded wavelength solution produced.

    A revised cube advertises more than its wavelength coordinate: the export
    also wrote each order's wavelength span and the pre-crop bounds of the whole
    exported axis.  Those are statements about the solution, not about the data,
    so leaving them behind after a revision publishes the old era's numbers under
    the new era's snapshot id -- the exact drift found where an order still
    advertised its old range, off by the recalibration shift.

    Both are rebuilt from the new polynomials over the same detector pixels the
    export used, which the cube still records: ``order_wavelength_ranges_nm_json``
    counts an order's finite pixels from zero (the partial-order pad is a tail),
    and the pre-crop axis is the order-border span intersected with that pad.
    Extraction geometry is frozen, so those pixel sets did not move -- only the
    wavelength each one maps to.
    """

    ranges = _json_records(revised.attrs, "order_wavelength_ranges_nm_json")
    finite_pixels: dict[int, int] = {}
    if ranges is not None:
        rebuilt = []
        for record in ranges:
            try:
                order = int(record["order"])
                n_px = int(record["n_px"])
            except (KeyError, TypeError, ValueError):
                raise RecalibrationError(
                    "cube has malformed order_wavelength_ranges_nm_json"
                ) from None
            finite_pixels[order] = n_px
            coefficients = solutions.get(order)
            if coefficients is None:
                raise RecalibrationError(
                    "new snapshot has no wavelength polynomial for order "
                    f"{order}, which the cube's order_wavelength_ranges_nm_json describes"
                )
            values = np.polyval(coefficients, np.arange(max(n_px, 0), dtype=float))
            rebuilt.append(
                {
                    "order": order,
                    "min_nm": float(np.min(values)) if values.size else float("nan"),
                    "max_nm": float(np.max(values)) if values.size else float("nan"),
                    "n_px": n_px,
                }
            )
        revised.attrs["order_wavelength_ranges_nm_json"] = json.dumps(
            rebuilt, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )

    borders = _json_records(revised.attrs, "order_border_pixel_ranges_json")
    has_bounds = (
        "original_wavelength_min_nm" in revised.attrs
        or "original_wavelength_max_nm" in revised.attrs
    )
    if borders is not None and has_bounds:
        spans = []
        for record in borders:
            try:
                order = int(record["order"])
                start = int(record["start_px"])
                stop = int(record["end_px"]) + 1
            except (KeyError, TypeError, ValueError):
                raise RecalibrationError(
                    "cube has malformed order_border_pixel_ranges_json"
                ) from None
            stop = min(stop, finite_pixels.get(order, stop))
            coefficients = solutions.get(order)
            if coefficients is None or stop <= start:
                continue
            spans.append(np.polyval(coefficients, np.arange(start, stop, dtype=float)))
        if spans:
            joined = np.concatenate(spans)
            revised.attrs["original_wavelength_min_nm"] = float(np.min(joined))
            revised.attrs["original_wavelength_max_nm"] = float(np.max(joined))


def _refresh_wavelength_accuracy(revised, snapshot: Snapshot, payload: dict[str, Any]) -> None:
    """Restate the cube's wavelength accuracy from the new snapshot's alignment.

    The accuracy is the snapshot's alignment RMS carried into nanometres through
    this cube's dispersion, so both halves of it belong to the superseded
    solution.  A new snapshot that records no alignment cannot make the claim at
    all, and the attribute is removed rather than left standing.
    """

    if "wavelength_accuracy_nm" not in revised.attrs:
        return
    from .tools.spectrocube_export import _wavelength_accuracy_nm

    records = {
        int(record["order"]): {"coefficients": list(record["coefficients"])}
        for record in payload["orders"]
    }
    accuracy = _wavelength_accuracy_nm(
        snapshot,
        records,
        np.asarray(revised.coords["detector_pixel"].values, dtype=float),
        np.asarray(revised.coords["echelle_order"].values, dtype=int),
    )
    if accuracy is None:
        revised.attrs.pop("wavelength_accuracy_nm", None)
        revised.attrs.pop("wavelength_accuracy_source", None)
    else:
        revised.attrs["wavelength_accuracy_nm"] = accuracy
        revised.attrs["wavelength_accuracy_source"] = "snapshot alignment rms_px"


def recalibrate_dataset(
    ds,
    new_snapshot: Snapshot,
    *,
    update_wavelength: bool = True,
    update_factor: bool = True,
    new_factor: np.ndarray | None = None,
):
    """Return a revised dataset while preserving unknown aligned fields."""

    old_id = str(ds.attrs.get("snapshot_id", ""))
    if not old_id or "snapshot_manifest_json" not in ds.attrs:
        raise RecalibrationError(
            "cube lacks complete Packet 8 snapshot provenance and cannot be revised safely"
        )
    old_records = _validated_old_records(ds)
    old_pattern = _artifact_digest(old_records, "pattern")
    new_pattern = _new_snapshot_digest(new_snapshot, "pattern")
    if not old_pattern:
        raise RecalibrationError("cube provenance does not identify its extraction pattern digest")
    if old_pattern != new_pattern:
        raise RecalibrationError(
            "extraction-pattern change refused: saved cubes cannot repair extraction geometry; "
            "reprocess the read-only raw SIF data with the new snapshot"
        )
    if not update_wavelength and not update_factor:
        raise RecalibrationError("select wavelength and/or absolute-factor recalibration")
    _require_unchanged_unapplied_roles(
        old_records,
        new_snapshot,
        update_wavelength=update_wavelength,
        update_factor=update_factor,
    )

    revised = ds.copy(deep=True)
    changes: list[str] = []
    polynomial_payload: dict[str, Any] | None = None
    solutions: dict[int, list[float]] = {}
    dropped_nonpositive_factor_columns = 0
    dropped_uncovered_factor_columns = 0
    if update_wavelength:
        wavelength, polynomial_payload, solutions = _evaluate_wavelength(revised, new_snapshot)
        revised = revised.assign_coords(wavelength=("wavelength", wavelength))
        changes.append("wavelength")

    if update_factor:
        revised, dropped_nonpositive_factor_columns, dropped_uncovered_factor_columns = (
            _replace_factor(revised, ds, new_snapshot, new_factor)
        )
        changes.append("absolute-factor")

    event = _append_history(
        ds,
        revised,
        new_snapshot,
        changes,
        dropped_nonpositive_factor_columns=dropped_nonpositive_factor_columns,
        dropped_uncovered_factor_columns=dropped_uncovered_factor_columns,
    )
    if polynomial_payload is not None:
        revised.attrs["wavelength_polynomials_json"] = json.dumps(
            polynomial_payload, sort_keys=True, separators=(",", ":")
        )
        # Every attribute the old solution produced is restated here, after the
        # factor pass has settled which samples survive, so the cube never
        # advertises one era's wavelengths under another era's snapshot id.
        _refresh_wavelength_derived_attrs(revised, solutions)
        _refresh_wavelength_accuracy(revised, new_snapshot, polynomial_payload)
    order = np.argsort(np.asarray(revised.coords["wavelength"].values), kind="stable")
    return revised.isel(wavelength=order), event


def recalibrate_cube(
    input_path: str | Path,
    output_path: str | Path,
    *,
    new_snapshot_path: str | Path,
    update_wavelength: bool = True,
    update_factor: bool = True,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Recalibrate one cube and write an immutable retroactive run manifest."""

    import xarray as xr
    from spectrocube import SpectroCube

    from .catalog import refresh_catalog_row

    source = Path(input_path)
    destination = Path(output_path)
    manifest_path = destination.with_suffix(destination.suffix + ".recalibration.json")
    if not overwrite and (destination.exists() or manifest_path.exists()):
        raise FileExistsError(f"recalibrated output or manifest already exists for {destination}")
    snapshot = load_snapshot(new_snapshot_path)
    with xr.open_dataset(source) as opened:
        revised, event = recalibrate_dataset(
            opened.load(),
            snapshot,
            update_wavelength=update_wavelength,
            update_factor=update_factor,
        )
    report = SpectroCube.from_dataset(revised).validate()
    if not report.ok:
        raise RecalibrationError("recalibrated cube is invalid: " + "; ".join(report.errors))
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        revised.to_netcdf(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    manifest = {
        **event,
        "input": source.name,
        "input_sha256": sha256_file(source),
        "output": destination.name,
        "output_sha256": sha256_file(destination),
    }
    temp_manifest = manifest_path.with_name(f".{manifest_path.name}.tmp")
    temp_manifest.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    os.replace(temp_manifest, manifest_path)
    # A revised cube makes every catalog row describing it stale. Refresh the
    # per-drive catalog beside the output and bump its generated_at, so the next
    # auto or manual merge carries the new digest and snapshot id into the
    # all-years index instead of leaving the old ones presented as current.
    refresh_catalog_row(destination)
    return destination, manifest_path
