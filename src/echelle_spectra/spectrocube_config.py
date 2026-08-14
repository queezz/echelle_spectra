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


def _relative_to_config(value: object, config_path: Path) -> str | None:
    if not value:
        return None
    candidate = Path(str(value))
    if not candidate.is_absolute():
        candidate = config_path.parent / candidate
    return str(candidate)


def load_toml(path: str | Path) -> dict[str, Any]:
    """Load a TOML file as a plain dictionary."""
    with Path(path).open("rb") as f:
        return tomllib.load(f)


def export_config_from_toml(path: str | Path) -> dict[str, Any]:
    """Return CLI-ready export defaults from a calibration config TOML."""
    config_path = Path(path)
    raw = load_toml(config_path)
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
        "registry": _relative_to_config(calibration.get("registry"), config_path),
        "calibrations": _relative_to_config(calibration.get("calibrations"), config_path),
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
    if metadata.get("trigger_delay_s") is not None:
        settings["extra_attrs"]["trigger_delay_s"] = float(metadata["trigger_delay_s"])
    if metadata.get("time_axis_reference"):
        settings["extra_attrs"]["time_axis_reference"] = metadata["time_axis_reference"]
    if metadata.get("frame_time_formula"):
        settings["extra_attrs"]["frame_time_formula"] = metadata["frame_time_formula"]
    if metadata.get("trigger_delay_note"):
        settings["extra_attrs"]["trigger_delay_note"] = metadata["trigger_delay_note"]

    return settings


def export_plan_from_toml(path: str | Path) -> dict[str, Any]:
    """Load a SpectroCube generation plan TOML."""
    plan_path = Path(path)
    raw = load_toml(plan_path)
    plan = raw.get("plan", {})
    for key in ("config", "registry", "calibrations"):
        if key not in plan:
            continue
        referenced = Path(plan[key])
        if not referenced.is_absolute():
            referenced = plan_path.parent / referenced
        plan[key] = str(referenced)
    return plan
