"""The LHD text header is frozen at its pre-unification legacy byte shape.

Both dialects are diffed line for line against templates recovered verbatim
from the pre-train commit ``859bf86`` and committed under ``tests/golden``.
A rendered header may differ from its golden template only in the substituted
values and in free text appended inside the existing ``[Comments]`` block.
"""

from __future__ import annotations

import os
from datetime import datetime
from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
import xarray as xr

from echelle_spectra.campaign_run import sha256_file
from echelle_spectra.lhd_text import (
    SPEC_DIV1,
    SPECTRUM,
    TEMPLATE_FILENAMES,
    TEXT_SCHEMA,
    LhdTextError,
    legacy_template,
    render_lhd_header,
    write_cube_text,
)

GOLDEN = Path(__file__).resolve().parent / "golden"
DATE = "08/14/2026 09:41"


def _header_block(path: Path) -> list[str]:
    """Every comment line a written table starts with."""

    lines = path.read_text(encoding="utf-8").split("\n")
    return lines[: lines.index("# [Data]") + 1]


def _split_date(lines: list[str]) -> tuple[list[str], str]:
    """Return the header without its Date line, and the Date line's value."""

    index = next(i for i, line in enumerate(lines) if line.startswith("# Date = "))
    value = lines[index].removeprefix("# Date = '").removesuffix("'")
    return lines[:index] + lines[index + 1 :], value


def _assert_written_now(stamp: str) -> None:
    assert datetime.strptime(stamp, "%m/%d/%Y %H:%M")


def test_packaged_templates_are_the_recovered_legacy_bytes() -> None:
    """The writer's templates are the golden bytes, not a re-typed likeness."""

    for dialect, golden in (
        (SPEC_DIV1, "lhd_header_spec_div1.txt"),
        (SPECTRUM, "lhd_header_spectrum.txt"),
    ):
        packaged = Path(
            files("echelle_spectra.resources").joinpath(TEMPLATE_FILENAMES[dialect])
        )
        reference = GOLDEN / golden
        assert packaged.read_bytes() == reference.read_bytes()
        assert legacy_template(dialect) == reference.read_text(
            encoding="utf-8"
        ).strip("\n")
        # `.gitattributes` pins both to LF, so no checkout can re-line-end the
        # frozen bytes and quietly change what LHD receives.
        assert b"\r" not in packaged.read_bytes()


def test_spec_div1_header_matches_the_legacy_template_line_for_line() -> None:
    values = {
        "diag_name": "spec_div1",
        "shot": 193778,
        "date": DATE,
        "dimno": 1,
        "dimname": "Time",
        "dimsize": 41,
        "dimunits": "s",
        "nval": 2,
        "vnames": "'CH','D'",
        "vunit": "'W/m^2','W/m^2'",
        "trigdelay": 2.5,
        "cycletime": 0.5,
    }
    expected = (GOLDEN / "lhd_header_spec_div1.txt").read_text(encoding="utf-8")
    expected = expected.strip("\n").format(**values)
    rendered = render_lhd_header(
        diagnostic="spec_div1",
        shot=193778,
        dimension_name="Time",
        dimension_size=41,
        dimension_unit="s",
        value_names=["CH", "D"],
        value_units=["W/m^2", "W/m^2"],
        trigger_delay_s=2.5,
        frame_interval_s=0.5,
        date=DATE,
    )
    assert rendered.split("\n") == expected.split("\n")
    assert rendered == expected
    # The frozen particulars this packet exists to defend.
    assert "# ShotNo = 193778" in rendered
    assert "# DimUnit = 's'" in rendered
    assert "# Contact: M. Kobayashi (2169)" in rendered
    assert "DimUnits" not in rendered
    assert "FormatSchema" not in rendered
    assert "[Provenance]" not in rendered


def test_spectrum_header_matches_the_legacy_template_line_for_line() -> None:
    values = {
        "shot": "193778",
        "date": DATE,
        "size": 4,
        "nval": 3,
        "vnames": "'0','1','2'",
        "vunit": "'W/(m2 nm)','W/(m2 nm)','W/(m2 nm)'",
        "trigdelay": "none",
        "cycletime": 0.5,
        "exposure": 0.25,
    }
    expected = (GOLDEN / "lhd_header_spectrum.txt").read_text(encoding="utf-8")
    expected = expected.strip("\n").format(**values)
    rendered = render_lhd_header(
        dialect=SPECTRUM,
        shot="193778",
        dimension_size=4,
        value_names=["0", "1", "2"],
        value_units=["W/(m2 nm)"] * 3,
        trigger_delay_s="none",
        frame_interval_s=0.5,
        exposure_s=0.25,
        date=DATE,
    )
    assert rendered.split("\n") == expected.split("\n")
    assert rendered == expected
    assert "# DimUnits = 'nm'" in rendered
    assert "# ShotNo = 193778" in rendered
    assert "# exposure = 0.25 (s)" in rendered


def test_spectrum_dialect_refuses_a_missing_exposure() -> None:
    with pytest.raises(LhdTextError, match="exposure"):
        render_lhd_header(
            dialect=SPECTRUM,
            shot="1",
            dimension_size=1,
            value_names=["0"],
            value_units=["counts"],
            trigger_delay_s=0.0,
            frame_interval_s=0.5,
        )


def test_appended_comments_stay_inside_the_frozen_comments_block() -> None:
    rendered = render_lhd_header(
        shot=1,
        dimension_size=2,
        value_names=["a"],
        value_units=["s"],
        trigger_delay_s=2.5,
        frame_interval_s=0.5,
        exposure_s=0.25,
        date=DATE,
        comments=["snapshot_id = 20250926_cmos"],
    )
    lines = rendered.split("\n")
    opened = lines.index("# [Comments]")
    closed = lines.index("# [Data]")
    assert opened < lines.index("# exposure = 0.25 (s)") < closed
    assert opened < lines.index("# snapshot_id = 20250926_cmos") < closed
    # The block still closes exactly as the legacy template closed it.
    assert lines[closed - 1] == "#"
    # Appending never adds a line above the block: the parameter section keeps
    # the template's own line positions.
    template = legacy_template(SPEC_DIV1).split("\n")
    assert opened == template.index("# [Comments]")


def _timed_cube(path: Path, **overrides: object) -> Path:
    attrs: dict[str, object] = {
        "instrument_id": "spec_div1",
        "intensity_units": "W/m2/nm/sr",
        "calibration_type": "absolute",
        "shot_number": "193778",
        "created_at": "2026-08-14T00:31:00+00:00",
        "snapshot_id": "20250926_cmos",
        "snapshot_manifest_sha256": "a" * 64,
        "trigger_delay_s": 2.5,
        "frame_interval_s": 0.5,
        "exposure_s": 0.25,
        "time_axis_reference": "LHD discharge time",
        "frame_time_formula": "trigger_delay_s + frame * frame_interval_s",
    }
    attrs.update(overrides)
    attrs = {key: value for key, value in attrs.items() if value is not None}
    xr.Dataset(
        {"intensity": (("frame", "wavelength"), np.array([[1.0, 2.0, 3.0]]))},
        coords={"frame": [0], "wavelength": np.array([410.0, 411.0, 412.0])},
        attrs=attrs,
    ).to_netcdf(path)
    return path


def test_cube_text_writes_the_frozen_header_with_cube_timing(tmp_path: Path) -> None:
    cube = _timed_cube(tmp_path / "193778.nc")
    output = write_cube_text(cube, tmp_path / "193778.txt")

    written, stamp = _split_date(_header_block(output))
    _assert_written_now(stamp)
    assert written == [
        "# [Parameters]",
        "# Name = 'spec_div1'",
        "# ShotNo = 193778",
        "#",
        "# DimNo = 1",
        "# DimName = 'wavelength'",
        "# DimSize = 3",
        "# DimUnit = 'nm'",
        "#",
        "# ValNo = 1",
        "# ValName = 'frame=0'",
        "# ValUnit = 'W/m2/nm/sr'",
        "#",
        "# [Comments]",
        "# The Echelle spectroscopy views the divertor region at 3.5L port from 4-O (AL-02).",
        '# For detailes of viewing area, please see "A. Kuzmin et al., PFR vol.13 (2018) 3402058.".',
        "# CH/D-band = 429.0-431.5 nm",
        "# Contact: M. Kobayashi (2169)",
        "# time = 2.5 + frameNo*0.5 (s)",
        "# exposure = 0.25 (s)",
        "# Text was derived from the saved cube; raw SIF data was not reopened.",
        f"# format_schema = {TEXT_SCHEMA}",
        "# cube_file = 193778.nc",
        f"# cube_sha256 = {sha256_file(cube)}",
        "# cube_created_at = 2026-08-14T00:31:00+00:00",
        "# snapshot_id = 20250926_cmos",
        f"# snapshot_manifest_sha256 = {'a' * 64}",
        "# time_axis_reference = LHD discharge time",
        "# frame_time_formula = trigger_delay_s + frame * frame_interval_s",
        "#",
        "# [Data]",
    ]


@pytest.mark.parametrize(
    ("attr", "origin"),
    [
        ("trigger_delay_s", "trigger_delay_s field"),
        ("frame_interval_s", "CycleTime"),
        ("exposure_s", "ExposureTime"),
    ],
)
def test_cube_text_refuses_a_cube_missing_one_timing_attribute(
    tmp_path: Path, attr: str, origin: str
) -> None:
    cube = _timed_cube(tmp_path / f"no-{attr}.nc", **{attr: None})
    with pytest.raises(LhdTextError) as raised:
        write_cube_text(cube, tmp_path / "out.txt")
    message = str(raised.value)
    assert attr in message
    assert origin in message
    assert not (tmp_path / "out.txt").exists()


def test_cube_text_refuses_a_counts_cube_because_lhd_must_be_absolute(
    tmp_path: Path,
) -> None:
    """The LHD deliverable carries an absolute scale or it is not written at all.

    A counts cube has valid timing and would otherwise produce a perfectly
    well-formed header reading ValUnit='counts' -- which is exactly the silent
    wrong-science this refusal exists to stop.
    """

    cube = _timed_cube(
        tmp_path / "counts.nc", intensity_units="counts", calibration_type="counts"
    )
    with pytest.raises(LhdTextError) as raised:
        write_cube_text(cube, tmp_path / "counts.txt")
    message = str(raised.value)
    assert "absolute" in message
    assert "wmsr" in message
    assert "recal-cube" in message
    assert not (tmp_path / "counts.txt").exists()


def test_cube_text_records_the_sphere_standard_the_cube_names(tmp_path: Path) -> None:
    """When the absolute cube names its sphere standard, the deliverable says so."""

    cube = _timed_cube(
        tmp_path / "sourced.nc",
        calibration_source="snapshot 20250926_cmos integrating sphere",
    )
    output = write_cube_text(cube, tmp_path / "sourced.txt")
    text = output.read_text(encoding="utf-8")
    assert "# calibration_source = snapshot 20250926_cmos integrating sphere" in text


class _Value:
    """A Qt spin box / check box as far as the save paths are concerned."""

    def __init__(self, value):
        self._value = value

    def value(self):
        return self._value

    def isChecked(self):
        return bool(self._value)


class _Band:
    def __init__(self, name):
        self.name = name


class _Fit:
    def __init__(self, values):
        self.intensities_combined = np.asarray(values, dtype=float)


class _Em:
    info = {"NumberOfFrames": 3, "CycleTime": 0.5}


def test_gui_band_save_writes_the_legacy_header(tmp_path: Path) -> None:
    """The GUI's own call path, not a re-creation of it."""

    from echelle_spectra.echelle_spectra_gui import EchelleSpectraGUI

    window = SimpleNamespace(
        bands=[_Band("CH"), _Band("D")],
        em=_Em(),
        trigger_delay=_Value(2.5),
        fres={"CH": _Fit([np.nan, 1.0, 2.0]), "D": _Fit([np.nan, 3.0, 4.0])},
        shot_number=_Value(193778),
        output_path=tmp_path,
        config={"diag_name": "spec_div1"},
        overwrite=_Value(True),
    )

    EchelleSpectraGUI.save_intensities(window)

    written, stamp = _split_date(_header_block(tmp_path / "193778_spec_div1.txt"))
    _assert_written_now(stamp)
    assert written == [
        "# [Parameters]",
        "# Name = 'spec_div1'",
        "# ShotNo = 193778",
        "#",
        "# DimNo = 1",
        "# DimName = 'Time'",
        "# DimSize = 2",
        "# DimUnit = 's'",
        "#",
        "# ValNo = 2",
        "# ValName = 'CH','D'",
        "# ValUnit = 'W/m^2','W/m^2'",
        "#",
        "# [Comments]",
        "# The Echelle spectroscopy views the divertor region at 3.5L port from 4-O (AL-02).",
        '# For detailes of viewing area, please see "A. Kuzmin et al., PFR vol.13 (2018) 3402058.".',
        "# CH/D-band = 429.0-431.5 nm",
        "# Contact: M. Kobayashi (2169)",
        "# time = 2.5 + frameNo*0.5 (s)",
        "#",
        "# [Data]",
    ]


def test_spectrum_save_writes_the_legacy_header(tmp_path: Path) -> None:
    """``Spectrum.save``'s own call path, through the plural-``DimUnits`` dialect."""

    from echelle_spectra.tools.echelle import Spectrum

    spectrum = SimpleNamespace(
        shotnumber="193778",
        wavelength=np.array([410.0, 411.0, 412.0, 413.0]),
        info={"NumberOfFrames": 2, "CycleTime": 0.5, "ExposureTime": 0.25},
        trigdelay="none",
        saveunits="wm",
        units_names={"wm": "W/(m2 nm)"},
        spectra_to_save={"wm": np.ones((2, 4))},
        output_path=str(tmp_path),
    )

    Spectrum.save(spectrum)

    written, stamp = _split_date(_header_block(tmp_path / "193778_echelle_spec.txt"))
    _assert_written_now(stamp)
    assert written == [
        "# [Parameters]",
        "# Name = 'Echelle Spectra'",
        "# ShotNo = 193778",
        "#",
        "# DimNo = 1",
        "# DimName = 'wavelength'",
        "# DimSize = 4",
        "# DimUnits = 'nm'",
        "#",
        "# ValNo = 2",
        "# ValName = '0','1'",
        "# ValUnit = 'W/(m2 nm)','W/(m2 nm)'",
        "#",
        "# [Comments]",
        "# time = none + frameNo*0.5 (s)",
        "# exposure = 0.25 (s)",
        "#",
        "# [Data]",
    ]
