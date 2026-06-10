"""Tests for wavelength-calibration QC helpers and CLI."""

from __future__ import annotations

import csv

import pytest

from echelle_spectra.tools.calibration_alignment import CalibrationTableLine
from echelle_spectra.tools.wavelength_calibration_qc import (
    find_focus_orders,
    fit_wavelength_orders,
    run_wavelength_calibration_qc,
)
from echelle_spectra.wavelength_qc_cli import main


def _line(order: int, center: float, wavelength: float, species: str = "NeI"):
    return CalibrationTableLine(
        order_idx=order,
        pixel_from=center - 4.0,
        pixel_to=center + 4.0,
        center_pixel=center,
        wavelength_nm=wavelength,
        species=species,
        comment="ok",
    )


def _write_table(path, rows) -> None:
    text = ["# test table", "# order from to center wavelength band"]
    for row in rows:
        text.append(
            f"{row.order_idx} {row.pixel_from:.3f} {row.pixel_to:.3f} "
            f"{row.center_pixel:.3f} {row.wavelength_nm:.6f} {row.species} # {row.comment}"
        )
    path.write_text("\n".join(text) + "\n", encoding="utf-8")


def test_fit_wavelength_orders_matches_calibration_degree_rule():
    lines = [
        _line(8, 100.0, 630.0),
        _line(8, 200.0, 629.0),
        _line(8, 300.0, 628.0),
        _line(9, 100.0, 615.0),
        _line(9, 200.0, 614.0),
    ]
    fits = fit_wavelength_orders("table", lines, detector_width_px=400)
    assert fits[8].degree == 2
    assert fits[9].degree == 1
    assert fits[8].rms_residual_nm < 1e-10
    assert fits[8].wavelength_min_nm < fits[8].wavelength_max_nm


def test_find_focus_orders_returns_overlapping_orders():
    lines = [
        _line(8, 100.0, 630.0),
        _line(8, 200.0, 629.0),
        _line(8, 300.0, 628.0),
        _line(20, 100.0, 480.0),
        _line(20, 200.0, 479.0),
        _line(20, 300.0, 478.0),
    ]
    fits = fit_wavelength_orders("table", lines, detector_width_px=400)
    assert find_focus_orders(fits, min_nm=620.0, max_nm=640.0) == [8]


def test_run_wavelength_calibration_qc_writes_tables_and_plots(tmp_path):
    table = tmp_path / "wave.txt"
    rows = [
        _line(8, 100.0, 630.0),
        _line(8, 200.0, 629.0),
        _line(8, 300.0, 628.0),
        _line(9, 100.0, 615.0),
        _line(9, 200.0, 614.0),
        _line(9, 300.0, 613.0),
    ]
    _write_table(table, rows)
    out = tmp_path / "out"

    result = run_wavelength_calibration_qc(
        [("test", table)],
        out,
        detector_width_px=400,
        focus_min_nm=620.0,
        focus_max_nm=640.0,
    )

    assert result.tables[0].name == "test"
    assert (out / "order_fit_summary.csv").is_file()
    assert (out / "line_fit_residuals.csv").is_file()
    assert (out / "test_all_orders.png").is_file()
    assert (out / "focus_620_640_orders.png").is_file()
    assert (out / "focus_620_640_line_residuals.png").is_file()

    with (out / "order_fit_summary.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert {row["order"] for row in rows} == {"8", "9"}


def test_wavelength_qc_cli_accepts_named_table(tmp_path):
    cal_dir = tmp_path / "cal"
    cal_dir.mkdir()
    table = cal_dir / "wave.txt"
    _write_table(
        table,
        [
            _line(8, 100.0, 630.0),
            _line(8, 200.0, 629.0),
            _line(8, 300.0, 628.0),
        ],
    )
    out = tmp_path / "out"

    with pytest.raises(SystemExit) as exc:
        main(
            [
                "--calibration-dir",
                str(cal_dir),
                "--table",
                "synthetic=wave.txt",
                "--output-dir",
                str(out),
                "--detector-width-px",
                "400",
            ]
        )

    assert exc.value.code == 0
    assert (out / "synthetic_all_orders.png").is_file()
