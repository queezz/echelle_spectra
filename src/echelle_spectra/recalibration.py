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
            + "; apply both deltas or choose a snapshot whose other roles are unchanged"
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


def _evaluate_wavelength(ds, snapshot: Snapshot) -> tuple[np.ndarray, dict[str, Any]]:
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
    return wavelength, _polynomial_payload(snapshot, represented_solutions)


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


def factor_from_snapshot(ds, snapshot: Snapshot) -> np.ndarray:
    """Load a snapshot's sphere pair and align its factor by raw order/pixel."""

    from .spectrocube_cli import _snapshot_camera
    from .tools.loader import build_calibration

    calibration = build_calibration(
        snapshot.root,
        _snapshot_camera(snapshot),
        calibration_files=snapshot.calibration_files(),
    )
    kind = _factor_kind(ds["applied_absolute_calibration_factor"])
    factor = np.asarray(calibration.absolute[kind], dtype=float)
    wavelength = np.asarray(calibration.wavelength, dtype=float)
    detector_grid = np.broadcast_to(
        np.arange(calibration.DIMW), calibration.order_borders.shape
    )[calibration.order_borders]
    order_grid = np.broadcast_to(
        np.asarray(calibration.order_ids)[:, None], calibration.order_borders.shape
    )[calibration.order_borders]
    # Length trap: Calibrations.wavelength keeps the partial-order NaN pad
    # (N_full) while EchelleImage.wavelength and every absolute factor are built
    # from the already-pruned sphere spectra (N_keep).  The N_full grids above
    # take the keep mask; the factor is already in wavelength[keep] order and
    # must not be masked again.
    keep = np.isfinite(wavelength)
    if factor.size != int(np.count_nonzero(keep)):
        raise RecalibrationError(
            f"snapshot absolute factor has {factor.size} samples for "
            f"{int(np.count_nonzero(keep))} retained wavelength samples"
        )
    lookup = {
        (int(order), int(pixel)): float(value)
        for order, pixel, value in zip(order_grid[keep], detector_grid[keep], factor)
    }
    try:
        return np.asarray(
            [
                lookup[(int(order), int(pixel))]
                for order, pixel in zip(
                    ds.coords["echelle_order"].values,
                    ds.coords["detector_pixel"].values,
                )
            ],
            dtype=float,
        )
    except KeyError as exc:
        raise RecalibrationError(
            f"new snapshot factor has no raw order/pixel sample {exc.args[0]}"
        ) from exc


def _replace_factor(revised, original, new_snapshot: Snapshot, new_factor: np.ndarray | None):
    """Apply a new absolute factor, dropping the samples it cannot calibrate.

    Export drops wavelength columns whose absolute factor is not strictly
    positive — one noise-negative sphere-minus-background column must not fail a
    whole file — and recalibration mirrors that tolerance: a sample whose old or
    new factor is non-finite or non-positive is removed and counted in the
    ``dropped_nonpositive_factor_columns`` attribute instead of refusing the
    revision.  Every retained sample keeps the strictly positive guarantee.
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
    usable = (
        np.isfinite(replacement)
        & (replacement > 0)
        & np.isfinite(previous)
        & (previous > 0)
    )
    dropped = int(usable.size - np.count_nonzero(usable))
    if not np.any(usable):
        raise RecalibrationError(
            "no wavelength sample has a strictly positive old and replacement absolute factor"
        )
    if dropped:
        revised = revised.isel(wavelength=np.flatnonzero(usable))
        replacement = replacement[usable]
        old = revised["applied_absolute_calibration_factor"]
        revised.attrs["dropped_nonpositive_factor_columns"] = (
            int(revised.attrs.get("dropped_nonpositive_factor_columns", 0)) + dropped
        )
    source_signal = revised["intensity"] / old
    revised["intensity"] = source_signal * replacement
    revised["intensity"].attrs.update(original["intensity"].attrs)
    revised["applied_absolute_calibration_factor"] = (
        ("wavelength",),
        replacement,
        dict(old.attrs),
    )
    return revised, dropped


def _append_history(
    ds,
    revised,
    new_snapshot: Snapshot,
    changes: list[str],
    *,
    dropped_nonpositive_factor_columns: int = 0,
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
    history.append(event)
    revised.attrs.update(new_snapshot.provenance_attrs())
    revised.attrs["recalibration_history_json"] = json.dumps(
        history, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    revised.attrs["previous_snapshot_manifest_json"] = str(ds.attrs["snapshot_manifest_json"])
    return event


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
    dropped_nonpositive_factor_columns = 0
    if update_wavelength:
        wavelength, polynomial_payload = _evaluate_wavelength(revised, new_snapshot)
        revised = revised.assign_coords(wavelength=("wavelength", wavelength))
        changes.append("wavelength")

    if update_factor:
        revised, dropped_nonpositive_factor_columns = _replace_factor(
            revised, ds, new_snapshot, new_factor
        )
        changes.append("absolute-factor")

    event = _append_history(
        ds,
        revised,
        new_snapshot,
        changes,
        dropped_nonpositive_factor_columns=dropped_nonpositive_factor_columns,
    )
    if polynomial_payload is not None:
        revised.attrs["wavelength_polynomials_json"] = json.dumps(
            polynomial_payload, sort_keys=True, separators=(",", ":")
        )
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
