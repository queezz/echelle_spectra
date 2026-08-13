from __future__ import annotations

from pathlib import Path

from echelle_spectra import cli, snapshot_cli
from echelle_spectra.snapshot import ROLE_FILENAMES


def _source_args(tmp_path: Path) -> list[str]:
    arguments = []
    option_names = {
        "pattern": "--pattern",
        "wavelength": "--wavelength",
        "sphere": "--sphere",
        "sphere_background": "--sphere-background",
        "integral": "--integral",
    }
    for role in ROLE_FILENAMES:
        path = tmp_path / f"{role}.source"
        path.write_text(role, encoding="utf-8")
        arguments.extend([option_names[role], str(path)])
    return arguments


def test_snapshot_cli_create_validate_and_show(tmp_path: Path, capsys) -> None:
    root = tmp_path / "calibrations"
    arguments = [
        "create",
        "20260813_cmos",
        "--output-root",
        str(root),
        "--detector",
        "cmos",
        "--lamp-used",
        "H2",
        *_source_args(tmp_path),
    ]
    assert snapshot_cli.main(arguments) == 0
    assert "Created calibration snapshot 20260813_cmos" in capsys.readouterr().out

    snapshot_path = root / "20260813_cmos"
    assert snapshot_cli.main(["validate", str(snapshot_path)]) == 0
    assert "all digests verified" in capsys.readouterr().out
    assert snapshot_cli.main(["show", str(snapshot_path)]) == 0
    shown = capsys.readouterr().out
    assert "detector:  cmos" in shown
    assert "lamps:     H2" in shown


def test_umbrella_status_reports_current_snapshot(tmp_path: Path, capsys) -> None:
    root = tmp_path / "calibrations"
    assert (
        snapshot_cli.main(
            [
                "create",
                "20260813_cmos",
                "--output-root",
                str(root),
                "--detector",
                "cmos",
                "--lamp-used",
                "H2",
                *_source_args(tmp_path),
            ]
        )
        == 0
    )
    capsys.readouterr()

    assert cli.main(["status", "--calibrations", str(root)]) == 0
    output = capsys.readouterr().out
    assert "snapshots: 1 valid" in output
    assert "current:   20260813_cmos (H2)" in output
    assert "runs:      none found" in output


def test_bare_umbrella_command_explains_the_next_steps(capsys) -> None:
    assert cli.main([]) == 0
    output = capsys.readouterr().out
    assert "Start with 'echelle status'" in output
    assert "snapshot" in output
    assert "process" in output
