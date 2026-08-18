"""Portable per-drive cube catalogs and merged all-years indexes."""

from __future__ import annotations

import json
import os
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .campaign_run import (
    GATE_UNRECORDED,
    RunReceipt,
    read_drive_identity,
    sha256_file,
)

DRIVE_CATALOG_SCHEMA = "echelle-drive-catalog/v1"
MERGED_CATALOG_SCHEMA = "echelle-merged-catalog/v1"
CATALOG_NAME = "echelle-catalog.json"

# One process may merge several drives into one central index concurrently, so
# read-modify-write of that index is serialized here.  Cross-process safety
# still rests on the atomic replace below.
_MERGE_LOCK = threading.Lock()


def _atomic_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _cube_record(path: Path, root: Path, *, gate: str = "") -> dict[str, Any]:
    import xarray as xr

    with xr.open_dataset(path) as ds:
        attrs = dict(ds.attrs)
        wavelength = ds.coords["wavelength"]
        return {
            "path": path.relative_to(root).as_posix(),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
            "shot_number": str(attrs.get("shot_number", path.stem)),
            "year": _record_year(attrs, path),
            "snapshot_id": str(attrs.get("snapshot_id", "unassigned")),
            "calibration_type": str(attrs.get("calibration_type", "unknown")),
            "intensity_units": str(attrs.get("intensity_units", "unknown")),
            "wavelength_min_nm": float(wavelength.min()),
            "wavelength_max_nm": float(wavelength.max()),
            "wavelength_points": int(wavelength.size),
            "source_file": str(attrs.get("source_file", "")),
            "t_start": str(attrs.get("t_start", "")),
            # How the run that produced this cube was authorized to touch its
            # calibration epoch, so a sampled drive never reads as a verified
            # one.  Cubes this catalog's receipt did not publish keep the
            # unrecorded word rather than borrowing the current run's gate.
            "gate": gate or GATE_UNRECORDED,
        }


def _record_year(attrs: dict[str, Any], path: Path) -> int | None:
    for key in ("t_start", "source_file"):
        value = str(attrs.get(key, ""))
        for token in value.replace("\\", "/").split("/"):
            if len(token) >= 4 and token[:4].isdigit() and 1900 <= int(token[:4]) <= 2200:
                return int(token[:4])
    for part in path.parts:
        if len(part) >= 4 and part[:4].isdigit() and 1900 <= int(part[:4]) <= 2200:
            return int(part[:4])
    created = str(attrs.get("created_at", ""))
    if len(created) >= 4 and created[:4].isdigit():
        return int(created[:4])
    return None


def discover_drive_id(root: str | Path) -> str:
    """Return the stable drive id announced at or above *root*, or an empty string.

    Cubes are often written into a folder on the drive rather than at its root,
    so the search walks upwards until a drive announces itself.  A drive that
    announces nothing keys on its volume label instead.
    """

    current = Path(os.path.abspath(root))
    for candidate in (current, *current.parents):
        identity = read_drive_identity(candidate)
        if identity is not None:
            return identity.drive_id
    return ""


def build_drive_catalog(
    cubes_root: str | Path,
    *,
    volume_label: str,
    drive_id: str | None = None,
    receipt_dir: str | Path | None = None,
    output: str | Path | None = None,
) -> Path:
    """Write a deterministic catalog beside one drive's cubes."""

    root = Path(cubes_root).resolve()
    run: dict[str, Any] | None = None
    gate = ""
    exported: set[str] = set()
    if receipt_dir is not None and (Path(receipt_dir) / "run.toml").is_file():
        receipt = RunReceipt.load(Path(receipt_dir))
        gate = receipt.gate
        exported = receipt.exported_outputs()
        run = {
            "id": receipt.directory.name,
            "state": receipt.state,
            "counts": receipt.counts(),
            "snapshot_id": receipt.snapshot_id,
            "gate": receipt.gate,
            "sample": receipt.sample,
        }
        if receipt.drive_warning:
            run["drive_warning"] = receipt.drive_warning
        if receipt.pruned_dirs:
            # Carried only when discovery actually skipped something, the same
            # way the receipt states it: the key's presence means a real skip,
            # and a catalog from a run that pruned nothing is unchanged.  A
            # reader of the merged index is otherwise the one person who cannot
            # tell a complete run from a pruned one.
            run["pruned_dirs"] = list(receipt.pruned_dirs)
    records = []
    errors = []
    for path in sorted(root.rglob("*.nc")):
        relative = path.relative_to(root).as_posix()
        try:
            records.append(
                _cube_record(path, root, gate=gate if relative in exported else "")
            )
        except (OSError, ValueError, KeyError) as exc:
            errors.append({"path": relative, "reason": str(exc)})
    payload = {
        "schema": DRIVE_CATALOG_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="milliseconds"),
        "volume_label": str(volume_label),
        "drive_id": str(discover_drive_id(root) if drive_id is None else drive_id),
        "cube_root": ".",
        "run": run,
        "cubes": records,
        "errors": errors,
    }
    return _atomic_json(Path(output) if output else root / CATALOG_NAME, payload)


def load_catalog(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") not in {DRIVE_CATALOG_SCHEMA, MERGED_CATALOG_SCHEMA}:
        raise ValueError(f"unsupported Echelle catalog schema in {path}")
    return payload


def source_key(source: dict[str, Any]) -> str:
    """Return the identity a merged index files one drive under.

    The stable drive id when the drive announced one; otherwise the volume
    label, which is all an unidentified drive can offer.
    """

    drive_id = str(source.get("drive_id", "")).strip()
    return drive_id or f"label:{str(source.get('volume_label', 'unknown')).strip()}"


def source_catalog_path(source: dict[str, Any]) -> Path:
    """Resolve one merged row's catalog file against its last known drive root."""

    return Path(str(source.get("drive_root", ""))) / str(source.get("catalog_path", ""))


def _source_from_drive_catalog(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    """Describe one per-drive catalog for the merged index.

    The catalog file's own location is that drive's root.  Rows therefore store
    the catalog path *relative* to the root plus the root this machine saw it
    at, so a drive that returns as another letter or under another mount point
    is relocated by re-merging its catalog rather than being marked permanently
    unavailable by an unusable absolute path from a different operating system.
    """

    resolved = path.resolve()
    return {
        "drive_id": str(payload.get("drive_id", "")),
        "volume_label": payload["volume_label"],
        "generated_at": str(payload.get("generated_at", "")),
        "drive_root": resolved.parent.as_posix(),
        "catalog_path": resolved.name,
        "available": resolved.is_file(),
        "run": payload.get("run"),
        "cubes": payload.get("cubes", []),
    }


def merge_catalogs(catalogs: list[str | Path], output: str | Path) -> Path:
    """Merge catalogs while retaining drive identity and last-known location.

    Recency decides, never argument order: for each drive the row with the
    newest ``generated_at`` wins, so merging a stale all-years index after a
    fresh per-drive catalog can no longer revert that drive to older rows.  Rows
    that carry no timestamp at all are older than any that does.
    """

    sources: list[dict[str, Any]] = []
    for raw in catalogs:
        path = Path(raw)
        payload = load_catalog(path)
        if payload["schema"] == MERGED_CATALOG_SCHEMA:
            sources.extend(payload.get("sources", []))
            continue
        sources.append(_source_from_drive_catalog(path, payload))
    by_drive: dict[str, dict[str, Any]] = {}
    for source in sources:
        key = source_key(source)
        current = by_drive.get(key)
        if current is None or str(source.get("generated_at", "")) >= str(
            current.get("generated_at", "")
        ):
            by_drive[key] = source
    for source in by_drive.values():
        source["available"] = source_catalog_path(source).is_file()
    payload = {
        "schema": MERGED_CATALOG_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="milliseconds"),
        "sources": [by_drive[key] for key in sorted(by_drive)],
    }
    return _atomic_json(Path(output), payload)


def merge_into_central_index(catalog: str | Path, central_index: str | Path) -> Path:
    """Fold one per-drive catalog into a durable central index, creating it if new."""

    index = Path(central_index)
    with _MERGE_LOCK:
        inputs: list[str | Path] = [index] if index.is_file() else []
        inputs.append(Path(catalog))
        return merge_catalogs(inputs, index)


def refresh_catalog_row(cube: str | Path) -> Path | None:
    """Update the adjacent per-drive catalog's row for one rewritten cube.

    ``recal-cube`` writes a new cube under a new snapshot id, which makes every
    catalog row describing it stale.  When a per-drive catalog sits beside the
    output — at the cube's folder or any folder above it — its row is rewritten
    with the new digest, size, and snapshot id, and the catalog's
    ``generated_at`` is bumped.  That bump is what a later auto-merge or manual
    ``echelle catalog merge`` needs: the refreshed drive now outranks whatever
    the central index still remembers, so the correction propagates without
    rebuilding the whole index.  Returns the catalog it refreshed, or None when
    the output has no catalog beside it.
    """

    path = Path(cube).resolve()
    for candidate in path.parents:
        catalog_path = candidate / CATALOG_NAME
        if not catalog_path.is_file():
            continue
        payload = load_catalog(catalog_path)
        if payload.get("schema") != DRIVE_CATALOG_SCHEMA:
            continue
        relative = path.relative_to(candidate).as_posix()
        record = _cube_record(path, candidate)
        rows = [row for row in payload.get("cubes", []) if row.get("path") != relative]
        previous = next(
            (row for row in payload.get("cubes", []) if row.get("path") == relative),
            None,
        )
        if previous is not None and previous.get("gate"):
            record["gate"] = str(previous["gate"])
        rows.append(record)
        payload["cubes"] = sorted(rows, key=lambda row: str(row.get("path", "")))
        payload["generated_at"] = datetime.now(timezone.utc).isoformat(timespec="milliseconds")
        return _atomic_json(catalog_path, payload)
    return None
