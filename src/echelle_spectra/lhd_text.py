"""LHD text export against a frozen header.

The LHD-side text header is a contract with the machine's data system, so it
stays at the byte shape it had before the Harbor unification (owner ruling,
2026-08-14).  Two legacy dialects survive, and both are rendered here from
templates recovered verbatim from the pre-train tree:

``spec_div1``
    ``resources/header_template.txt`` -- the GUI band save and every
    cube-derived deliverable.  ``DimUnit`` (singular) and the LHD
    viewing-geometry ``[Comments]`` block.

``spectrum``
    ``resources/header_template_spectrum.txt`` -- ``Spectrum.save``.
    ``DimUnits`` (plural), a fixed wavelength dimension, and an ``exposure``
    comment line.

``ShotNo`` is unquoted and ``Date`` is local ``%m/%d/%Y %H:%M`` in both.  No
field may be added outside the templates: provenance rides only as extra
free-text comment lines appended inside the existing ``[Comments]`` block.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime
from functools import cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np

from .campaign_run import sha256_file

TEXT_SCHEMA = "echelle-lhd-text/v1"

#: The GUI/`spec_div1` dialect: LHD's deliverable shape.
SPEC_DIV1 = "spec_div1"
#: The `Spectrum.save` dialect.
SPECTRUM = "spectrum"

#: The recovered legacy template backing each dialect, under ``resources``.
TEMPLATE_FILENAMES = {
    SPEC_DIV1: "header_template.txt",
    SPECTRUM: "header_template_spectrum.txt",
}

#: Timing attributes the frozen header needs, and where a cube gets them.
TIMING_ATTR_SOURCES = {
    "trigger_delay_s": "the export config's [metadata] trigger_delay_s field",
    "frame_interval_s": "the detector's CycleTime, recorded by export_spectrocube",
    "exposure_s": "the detector's ExposureTime, recorded by export_spectrocube",
}

_LEGACY_DATE_FORMAT = "%m/%d/%Y %H:%M"


class LhdTextError(ValueError):
    """A text export that cannot be written honestly against the frozen header."""


@cache
def legacy_template(dialect: str) -> str:
    """Return one frozen header template, exactly as the legacy tree carried it."""

    try:
        name = TEMPLATE_FILENAMES[dialect]
    except KeyError:
        raise LhdTextError(
            f"unknown LHD text dialect {dialect!r}; expected one of "
            f"{', '.join(sorted(TEMPLATE_FILENAMES))}"
        ) from None
    resource = files("echelle_spectra.resources").joinpath(name)
    return resource.read_text(encoding="utf-8").strip("\n")


def _quoted(values: Sequence[object]) -> str:
    """Quote a value list the way both legacy writers did: ``'a','b'``."""

    return ",".join(f"'{value}'" for value in values)


def _append_comments(header: str, comments: Sequence[str]) -> str:
    """Insert free text at the end of the template's own ``[Comments]`` block."""

    if not comments:
        return header
    lines = header.split("\n")
    try:
        data_at = lines.index("# [Data]")
    except ValueError:  # pragma: no cover - templates are frozen and verified
        raise LhdTextError("frozen header template lost its [Data] marker") from None
    if data_at == 0 or lines[data_at - 1] != "#":
        raise LhdTextError(  # pragma: no cover - templates are frozen and verified
            "frozen header template lost the blank comment line before [Data]"
        )
    lines[data_at - 1 : data_at - 1] = [f"# {line}" for line in comments]
    return "\n".join(lines)


def render_lhd_header(
    *,
    shot: object,
    dimension_size: int,
    value_names: Sequence[object],
    value_units: Sequence[object],
    trigger_delay_s: object,
    frame_interval_s: object,
    dialect: str = SPEC_DIV1,
    diagnostic: str = "spec_div1",
    dimension_name: str = "Time",
    dimension_unit: str = "s",
    exposure_s: object | None = None,
    date: str | None = None,
    comments: Sequence[str] = (),
) -> str:
    """Render one frozen LHD header.

    ``diagnostic``, ``dimension_name`` and ``dimension_unit`` apply to the
    ``spec_div1`` dialect only; the ``spectrum`` template fixes its own
    ``Name``, ``DimName`` and ``DimUnits`` and requires ``exposure_s``.  For
    ``spec_div1`` an ``exposure_s`` is written as the first appended comment,
    because that dialect's frozen block has no exposure line of its own.
    """

    template = legacy_template(dialect)
    stamp = date or datetime.now().strftime(_LEGACY_DATE_FORMAT)
    common = {
        "shot": shot,
        "date": stamp,
        "nval": len(value_names),
        "vnames": _quoted(value_names),
        "vunit": _quoted(value_units),
        "trigdelay": trigger_delay_s,
        "cycletime": frame_interval_s,
    }
    extra = list(comments)
    if dialect == SPECTRUM:
        if exposure_s is None:
            raise LhdTextError(
                "the spectrum dialect's frozen header carries an exposure line; "
                "pass exposure_s"
            )
        substitutions: dict[str, object] = {
            **common,
            "size": int(dimension_size),
            "exposure": exposure_s,
        }
    else:
        substitutions = {
            **common,
            "diag_name": diagnostic,
            "dimno": 1,
            "dimname": dimension_name,
            "dimsize": int(dimension_size),
            "dimunits": dimension_unit,
        }
        if exposure_s is not None:
            extra.insert(0, f"exposure = {exposure_s} (s)")
    return _append_comments(template.format(**substitutions), extra)


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


def cube_timing(attrs: Mapping[str, object], *, cube: str) -> dict[str, float]:
    """Return the three timing values the frozen header needs, or refuse by name."""

    timing: dict[str, float] = {}
    for attr, origin in TIMING_ATTR_SOURCES.items():
        value = attrs.get(attr)
        if value is None or value == "":
            raise LhdTextError(
                f"{cube} is missing the timing attribute {attr!r}, which the frozen "
                f"LHD header must state; {origin} supplies it. Re-export the cube "
                "with a config that provides it rather than writing text without it."
            )
        try:
            timing[attr] = float(value)
        except (TypeError, ValueError):
            raise LhdTextError(
                f"{cube} carries a non-numeric timing attribute {attr!r}: {value!r}"
            ) from None
    return timing


def _provenance_comments(source: Path, attrs: Mapping[str, object]) -> list[str]:
    """Provenance as free text, which may live only inside ``[Comments]``."""

    lines = [
        "Text was derived from the saved cube; raw SIF data was not reopened.",
        f"format_schema = {TEXT_SCHEMA}",
        f"cube_file = {source.name}",
        f"cube_sha256 = {sha256_file(source)}",
    ]
    optional = (
        ("cube_created_at", "created_at"),
        ("source_file", "source_file"),
        ("snapshot_id", "snapshot_id"),
        ("snapshot_manifest_sha256", "snapshot_manifest_sha256"),
        ("calibration_registry_sha256", "calibration_registry_sha256"),
        ("calibration_registry_epoch_position", "calibration_registry_epoch_position"),
        ("time_axis_reference", "time_axis_reference"),
        ("frame_time_formula", "frame_time_formula"),
    )
    for label, attr in optional:
        value = attrs.get(attr)
        if value is None or value == "":
            continue
        lines.append(f"{label} = {value}")
    return lines


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
    timing = cube_timing(attrs, cube=source.name)
    header = render_lhd_header(
        diagnostic=str(attrs.get("instrument_id", "Echelle Spectra")),
        shot=str(attrs.get("shot_number", source.stem)),
        dimension_name="wavelength",
        dimension_size=wavelength.size,
        dimension_unit="nm",
        value_names=value_names,
        value_units=[units] * len(value_names),
        trigger_delay_s=timing["trigger_delay_s"],
        frame_interval_s=timing["frame_interval_s"],
        exposure_s=timing["exposure_s"],
        comments=_provenance_comments(source, attrs),
    )
    table = np.column_stack((wavelength, flattened.T))
    return write_lhd_text(
        output,
        table,
        header=header,
        formats=["%.9f", *(["%.9e"] * flattened.shape[0])],
        overwrite=overwrite,
    )
