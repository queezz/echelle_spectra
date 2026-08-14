"""Portable per-drive cube catalogs and merged all-years indexes."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .campaign_run import RunReceipt, sha256_file

DRIVE_CATALOG_SCHEMA = "echelle-drive-catalog/v1"
MERGED_CATALOG_SCHEMA = "echelle-merged-catalog/v1"
CATALOG_NAME = "echelle-catalog.json"


def _atomic_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    os.replace(temporary, path)
    return path


def _cube_record(path: Path, root: Path) -> dict[str, Any]:
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


def build_drive_catalog(
    cubes_root: str | Path,
    *,
    volume_label: str,
    receipt_dir: str | Path | None = None,
    output: str | Path | None = None,
) -> Path:
    """Write a deterministic catalog beside one drive's cubes."""

    root = Path(cubes_root).resolve()
    records = []
    errors = []
    for path in sorted(root.rglob("*.nc")):
        try:
            records.append(_cube_record(path, root))
        except (OSError, ValueError, KeyError) as exc:
            errors.append({"path": path.relative_to(root).as_posix(), "reason": str(exc)})
    run: dict[str, Any] | None = None
    if receipt_dir is not None and (Path(receipt_dir) / "run.toml").is_file():
        receipt = RunReceipt.load(Path(receipt_dir))
        run = {
            "id": receipt.directory.name,
            "state": receipt.state,
            "counts": receipt.counts(),
            "snapshot_id": receipt.snapshot_id,
        }
    payload = {
        "schema": DRIVE_CATALOG_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "volume_label": str(volume_label),
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


def merge_catalogs(catalogs: list[str | Path], output: str | Path) -> Path:
    """Merge catalogs while retaining volume identity and last-known location."""

    sources = []
    for raw in catalogs:
        path = Path(raw)
        payload = load_catalog(path)
        if payload["schema"] == MERGED_CATALOG_SCHEMA:
            sources.extend(payload.get("sources", []))
            continue
        sources.append(
            {
                "volume_label": payload["volume_label"],
                "catalog_path": str(path.resolve()),
                "available": path.is_file(),
                "run": payload.get("run"),
                "cubes": payload.get("cubes", []),
            }
        )
    by_volume = {str(source["volume_label"]): source for source in sources}
    for source in by_volume.values():
        source["available"] = Path(str(source.get("catalog_path", ""))).is_file()
    payload = {
        "schema": MERGED_CATALOG_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": [by_volume[key] for key in sorted(by_volume)],
    }
    return _atomic_json(Path(output), payload)
