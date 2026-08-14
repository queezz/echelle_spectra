"""Canonical provenance-complete LHD text export."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .campaign_run import sha256_file

TEXT_SCHEMA = "echelle-lhd-text/v1"


def _quoted(values: Sequence[object]) -> str:
    return ",".join(repr(str(value)) for value in values)


def render_lhd_header(
    *,
    diagnostic: str,
    shot: str,
    dimension_name: str,
    dimension_size: int,
    dimension_unit: str,
    value_names: Sequence[object],
    value_units: Sequence[object],
    provenance: Mapping[str, object] | None = None,
    created_at: str | None = None,
    comments: Sequence[str] = (),
) -> str:
    """Render the single LHD header used by every text-writing surface."""

    created = created_at or datetime.now(timezone.utc).isoformat(timespec="seconds")
    provenance = provenance or {}
    lines = [
        "# [Parameters]",
        f"# FormatSchema = {TEXT_SCHEMA!r}",
        f"# Name = {diagnostic!r}",
        f"# ShotNo = {str(shot)!r}",
        f"# Date = {created!r}",
        "#",
        "# DimNo = 1",
        f"# DimName = {dimension_name!r}",
        f"# DimSize = {int(dimension_size)}",
        f"# DimUnit = {dimension_unit!r}",
        "#",
        f"# ValNo = {len(value_names)}",
        f"# ValName = {_quoted(value_names)}",
        f"# ValUnit = {_quoted(value_units)}",
        "#",
        "# [Provenance]",
    ]
    for key, value in sorted(provenance.items()):
        if value is None or value == "":
            continue
        rendered = (
            json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            if isinstance(value, (dict, list, tuple))
            else str(value)
        )
        lines.append(f"# {key} = {rendered}")
    lines.extend(["#", "# [Comments]"])
    lines.extend(f"# {line}" for line in comments)
    lines.extend(["#", "# [Data]"])
    return "\n".join(lines)


def write_lhd_text(
    output: str | Path,
    data: np.ndarray,
    *,
    header: str,
    formats: Sequence[str] | str,
    overwrite: bool = False,
) -> Path:
    """Atomically write an LHD table and refuse replacement by default."""

    destination = Path(output)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"text output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        np.savetxt(
            temporary,
            np.asarray(data),
            delimiter=", ",
            header=header.strip(),
            comments="",
            fmt=formats,
        )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def write_cube_text(
    cube_path: str | Path,
    output: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write saved SpectroCube intensity without reopening the raw detector file."""

    import xarray as xr

    source = Path(cube_path)
    with xr.open_dataset(source) as opened:
        ds = opened.load()
    wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
    intensity = ds["intensity"]
    other_dims = [dim for dim in intensity.dims if dim != "wavelength"]
    ordered = intensity.transpose(*other_dims, "wavelength")
    flattened = np.asarray(ordered.values, dtype=float).reshape((-1, wavelength.size))
    if not other_dims:
        flattened = flattened.reshape((1, wavelength.size))
    value_names = []
    if other_dims:
        for index in np.ndindex(*(ordered.sizes[dim] for dim in other_dims)):
            value_names.append("/".join(f"{dim}={value}" for dim, value in zip(other_dims, index)))
    else:
        value_names = ["intensity"]
    units = str(ds.attrs.get("intensity_units", "unknown"))
    attrs: dict[str, Any] = dict(ds.attrs)
    provenance = {
        "cube_path": source.name,
        "cube_sha256": sha256_file(source),
        "source_file": attrs.get("source_file", ""),
        "snapshot_id": attrs.get("snapshot_id", "unassigned"),
        "snapshot_manifest_sha256": attrs.get("snapshot_manifest_sha256", ""),
        "snapshot_manifest_json": attrs.get("snapshot_manifest_json", ""),
        "calibration_file_digests_json": attrs.get("calibration_file_digests_json", ""),
        "calibration_registry_schema": attrs.get("calibration_registry_schema", ""),
        "calibration_registry_sha256": attrs.get("calibration_registry_sha256", ""),
        "calibration_registry_epoch_position": attrs.get(
            "calibration_registry_epoch_position", ""
        ),
        "wavelength_polynomials_json": attrs.get("wavelength_polynomials_json", ""),
        "applied_factor_application": ds.get("applied_absolute_calibration_factor", {}).attrs.get(
            "application", ""
        )
        if "applied_absolute_calibration_factor" in ds
        else "",
        "recalibration_history_json": attrs.get("recalibration_history_json", ""),
    }
    header = render_lhd_header(
        diagnostic=str(attrs.get("instrument_id", "Echelle Spectra")),
        shot=str(attrs.get("shot_number", source.stem)),
        dimension_name="wavelength",
        dimension_size=wavelength.size,
        dimension_unit="nm",
        value_names=value_names,
        value_units=[units] * len(value_names),
        provenance=provenance,
        created_at=str(attrs.get("created_at", "")) or None,
        comments=("Text was derived from the saved cube; raw SIF data was not reopened.",),
    )
    table = np.column_stack((wavelength, flattened.T))
    return write_lhd_text(
        output,
        table,
        header=header,
        formats=["%.9f", *(["%.9e"] * flattened.shape[0])],
        overwrite=overwrite,
    )
