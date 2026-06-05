"""Configuration helpers for reproducible SpectroCube export."""

from __future__ import annotations

from pathlib import Path
from typing import Any

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised only on older Python
    import tomli as tomllib


_CALIBRATION_FILE_KEYS = {
    "order_pattern": "orders",
    "wavelength": "wavelength",
    "sphere": "sphr",
    "sphere_background": "bkgr",
    "integral": "integral",
}


def load_toml(path: str | Path) -> dict[str, Any]:
    """Load a TOML file as a plain dictionary."""
    with Path(path).open("rb") as f:
        return tomllib.load(f)


def export_config_from_toml(path: str | Path) -> dict[str, Any]:
    """Return CLI-ready export defaults from a calibration config TOML."""
    raw = load_toml(path)
    calibration = raw.get("calibration", {})
    export = raw.get("export", {})
    metadata = raw.get("metadata", {})

    calibration_files = {
        dst: str(calibration[src])
        for src, dst in _CALIBRATION_FILE_KEYS.items()
        if calibration.get(src)
    }

    settings: dict[str, Any] = {
        "camera": calibration.get("camera"),
        "calibration_dir": calibration.get("calibration_dir"),
        "instrument_id": calibration.get("instrument_id"),
        "wavelength_medium": calibration.get("wavelength_medium"),
        "calibration_files": calibration_files,
        "units": export.get("units"),
        "wavelength_min_nm": export.get("wavelength_min_nm"),
        "wavelength_max_nm": export.get("wavelength_max_nm"),
        "drop_nonfinite_columns": export.get("drop_nonfinite_columns"),
        "calibration_source": export.get("calibration_source"),
        "output_suffix": export.get("output_suffix"),
        "extra_attrs": {},
    }

    if metadata.get("config_id"):
        settings["extra_attrs"]["spectrocube_export_config_id"] = metadata["config_id"]
    if metadata.get("crop_measurement_note"):
        settings["extra_attrs"]["wavelength_crop_note"] = metadata["crop_measurement_note"]
    if metadata.get("crop_measured_at"):
        settings["extra_attrs"]["wavelength_crop_measured_at"] = metadata["crop_measured_at"]

    return settings


def export_plan_from_toml(path: str | Path) -> dict[str, Any]:
    """Load a SpectroCube generation plan TOML."""
    raw = load_toml(path)
    plan = raw.get("plan", {})
    if "config" in plan:
        config_path = Path(plan["config"])
        if not config_path.is_absolute():
            config_path = Path(path).parent / config_path
        plan["config"] = str(config_path)
    return plan

