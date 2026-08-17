"""Every campaign command answers a bad path, and never raises one at the operator.

The field failure these pin: a documentation example written with Unix-shaped
placeholder paths was pasted into PowerShell, where ``/data/all-years.json`` is
a perfectly legal path on the current drive.  Nothing refused it near the
command line; the run died several frames inside pathlib and printed a
traceback naming a file the operator had never typed.

So each test here asks the same three questions of one refusal: does it name
the absolute path *this machine* looked at, does it name the flag that supplied
it, and does it name the command that writes such a file -- with no traceback
and a nonzero exit.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from echelle_spectra import campaign_tools_cli, cli

MERGED_CATALOG = {
    "schema": "echelle-merged-catalog/v1",
    "generated_at": "2026-08-17T00:00:00.000+00:00",
    "sources": [],
}


def _shown(capsys) -> str:
    captured = capsys.readouterr()
    return captured.out + captured.err


def _catalog(tmp_path: Path) -> Path:
    path = tmp_path / "all-years.json"
    path.write_text(json.dumps(MERGED_CATALOG), encoding="utf-8")
    return path


def _cube_without_timing(tmp_path: Path) -> Path:
    """A real, openable cube that the frozen LHD header cannot be built from."""

    path = tmp_path / "6001.nc"
    dataset = xr.Dataset(
        {"intensity": (("wavelength",), np.linspace(1.0, 2.0, 8))},
        coords={"wavelength": np.linspace(400.0, 700.0, 8)},
        attrs={"shot_number": "6001", "intensity_units": "counts"},
    )
    dataset.to_netcdf(path)
    return path


# ---------------------------------------------------------------------------
# The reported failure: echelle web
# ---------------------------------------------------------------------------


def test_missing_catalog_names_the_path_the_flag_and_the_producing_command(
    tmp_path: Path, capsys
) -> None:
    missing = tmp_path / "all-years.json"

    assert cli.main(["web", "--catalog", str(missing), "--output", str(tmp_path / "page")]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(missing.resolve()) in shown
    assert "--catalog" in shown
    assert "echelle catalog build" in shown
    assert len([line for line in shown.splitlines() if line.strip()]) == 1


def test_a_relative_catalog_is_answered_with_the_absolute_path_it_looked_at(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """The whole defect: the operator has to be told *where* the process looked."""

    monkeypatch.chdir(tmp_path)

    assert cli.main(["web", "--catalog", "all-years.json", "--output", "page"]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str((tmp_path / "all-years.json").resolve()) in shown


def test_malformed_catalog_json_names_the_file_and_the_parse_error(
    tmp_path: Path, capsys
) -> None:
    broken = tmp_path / "all-years.json"
    broken.write_text('{"schema": "echelle-merged-catalog/v1",', encoding="utf-8")

    assert cli.main(["web", "--catalog", str(broken), "--output", str(tmp_path / "page")]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "not valid JSON" in shown
    assert str(broken.resolve()) in shown
    assert "--catalog" in shown


def test_a_catalog_of_the_wrong_schema_is_refused_by_name(tmp_path: Path, capsys) -> None:
    stranger = tmp_path / "all-years.json"
    stranger.write_text(json.dumps({"schema": "something-else/v1"}), encoding="utf-8")

    assert cli.main(["web", "--catalog", str(stranger), "--output", str(tmp_path / "page")]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "schema" in shown and "--catalog" in shown


@pytest.mark.parametrize(
    ("flag", "value", "expected"),
    [
        ("--registry", "calibration_registry.toml", "epoch registry not found"),
        ("--calibrations", "calibrations", "calibration snapshot root not found"),
        ("--drift", "epoch-drift.json", "drift evidence not found"),
        ("--document", "notes.md", "Markdown document not found"),
    ],
)
def test_web_refuses_every_other_missing_input_by_its_own_flag(
    tmp_path: Path, capsys, flag: str, value: str, expected: str
) -> None:
    missing = tmp_path / value

    code = cli.main(
        [
            "web",
            "--catalog",
            str(_catalog(tmp_path)),
            "--output",
            str(tmp_path / "page"),
            flag,
            str(missing),
        ]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert expected in shown
    assert str(missing.resolve()) in shown
    assert flag in shown


def test_web_refuses_an_output_whose_parent_folder_does_not_exist(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "typo" / "campaign-page"

    code = cli.main(
        ["web", "--catalog", str(_catalog(tmp_path)), "--output", str(output)]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(output.parent.resolve()) in shown
    assert "--output" in shown
    assert not output.exists()


def test_a_malformed_registry_names_the_toml_parse_error(tmp_path: Path, capsys) -> None:
    registry = tmp_path / "calibration_registry.toml"
    registry.write_text("[[epoch]\nsnapshot_id = ", encoding="utf-8")

    code = cli.main(
        [
            "web",
            "--catalog",
            str(_catalog(tmp_path)),
            "--output",
            str(tmp_path / "page"),
            "--registry",
            str(registry),
        ]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "not valid TOML" in shown
    assert str(registry.resolve()) in shown


def test_the_pasted_field_command_leaves_a_real_process_without_a_traceback(
    tmp_path: Path,
) -> None:
    """The verbatim shape of the reported failure, run as its own process."""

    missing = tmp_path / "data" / "all-years.json"
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from echelle_spectra.cli import main; sys.exit(main(sys.argv[1:]))",
            "web",
            "--catalog",
            str(missing),
            "--output",
            str(tmp_path / "data" / "campaign-page"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "Traceback" not in completed.stderr
    assert "Traceback" not in completed.stdout
    assert str(missing.resolve()) in completed.stderr
    assert "--catalog" in completed.stderr


# ---------------------------------------------------------------------------
# The sibling verbs
# ---------------------------------------------------------------------------


def test_txt_refuses_a_missing_cube_by_absolute_path(tmp_path: Path, capsys) -> None:
    missing = tmp_path / "6001.nc"

    assert cli.main(["txt", str(missing), str(tmp_path / "6001.txt")]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(missing.resolve()) in shown
    assert "INPUT" in shown
    assert "echelle process" in shown


def test_txt_delivers_a_cube_refusal_as_one_line_rather_than_a_traceback(
    tmp_path: Path, capsys
) -> None:
    """The message was always good; it just used to arrive wrapped in frames."""

    cube = _cube_without_timing(tmp_path)

    assert cli.main(["txt", str(cube), str(tmp_path / "6001.txt")]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "timing attribute" in shown
    assert not (tmp_path / "6001.txt").exists()


def test_catalog_build_refuses_a_missing_cube_folder(tmp_path: Path, capsys) -> None:
    missing = tmp_path / "cubes"

    code = campaign_tools_cli.catalog_main(
        ["build", str(missing), "--volume-label", "NIFS-A"]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(missing.resolve()) in shown
    assert "echelle process" in shown


def test_catalog_merge_refuses_a_malformed_catalog(tmp_path: Path, capsys) -> None:
    broken = tmp_path / "catalog.json"
    broken.write_text("not json at all", encoding="utf-8")

    code = campaign_tools_cli.catalog_main(
        ["merge", str(broken), "-o", str(tmp_path / "all-years.json")]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "not valid JSON" in shown
    assert str(broken.resolve()) in shown


def test_recal_cube_refuses_a_missing_snapshot_and_says_who_writes_one(
    tmp_path: Path, capsys
) -> None:
    cube = _cube_without_timing(tmp_path)
    missing = tmp_path / "calibrations" / "20260814_cmos"

    code = campaign_tools_cli.recal_main(
        [str(cube), "-o", str(tmp_path / "out.nc"), "--new-snapshot", str(missing)]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(missing.resolve()) in shown
    assert "--new-snapshot" in shown
    assert "echelle snapshot create" in shown


def test_drift_refine_refuses_a_missing_evidence_file(tmp_path: Path, capsys) -> None:
    missing = tmp_path / "epoch-drift.json"
    (tmp_path / "calibrations").mkdir()

    code = campaign_tools_cli.drift_main(
        [
            "refine",
            str(missing),
            "--calibrations",
            str(tmp_path / "calibrations"),
            "--accept-shift",
            "1.5",
        ]
    )

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert str(missing.resolve()) in shown
    assert "echelle drift audit" in shown


def test_status_reports_where_it_looked_not_what_was_typed(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    monkeypatch.chdir(tmp_path)

    assert cli.main(["status", "--calibrations", "calibrations", "--runs", "runs"]) == 0

    shown = _shown(capsys)
    assert str((tmp_path / "calibrations").resolve()) in shown
    assert str((tmp_path / "runs").resolve()) in shown


# ---------------------------------------------------------------------------
# The immutable verdict file: a second audit is the normal path
# ---------------------------------------------------------------------------


def test_drift_audit_refuses_a_taken_evidence_name_and_offers_the_next_one(
    tmp_path: Path, capsys
) -> None:
    cubes = tmp_path / "cubes"
    cubes.mkdir()
    taken = tmp_path / "epoch-drift.json"
    taken.write_text("{}", encoding="utf-8")

    code = campaign_tools_cli.drift_main(["audit", str(cubes), "-o", str(taken)])

    assert code == 1
    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "immutable" in shown
    assert str(taken.resolve()) in shown
    assert str((tmp_path / "epoch-drift-2.json").resolve()) in shown
    # Immutability stands: nothing was written over the existing verdict.
    assert taken.read_text(encoding="utf-8") == "{}"


def test_the_suggested_evidence_name_walks_past_the_names_already_used(
    tmp_path: Path, capsys
) -> None:
    cubes = tmp_path / "cubes"
    cubes.mkdir()
    for name in ("epoch-drift.json", "epoch-drift-2.json", "epoch-drift-3.json"):
        (tmp_path / name).write_text("{}", encoding="utf-8")

    assert (
        campaign_tools_cli.drift_main(
            ["audit", str(cubes), "-o", str(tmp_path / "epoch-drift-2.json")]
        )
        == 1
    )

    assert str((tmp_path / "epoch-drift-4.json").resolve()) in _shown(capsys)


# ---------------------------------------------------------------------------
# The trap that produced the failure in the first place
# ---------------------------------------------------------------------------

# Invented Unix-shaped roots are the dangerous kind of placeholder: `/data/...`
# reads like a real campaign root, and pasted into PowerShell it does not fail
# the way a placeholder should -- it quietly becomes a path on the current
# drive.  Roots that announce themselves as placeholders (`/path/to/...`) and
# real roots inside a labelled POSIX block (`/Volumes/...`) are not the trap.
INVENTED_ROOTS = ("data", "lab")
_FENCE = re.compile(r"```[a-z]*\n(.*?)```", flags=re.S)


def _doc_code_blocks() -> list[tuple[str, str]]:
    docs = Path(__file__).resolve().parents[1] / "docs"
    return [
        (page.name, block)
        for page in sorted(docs.rglob("*.md"))
        for block in _FENCE.findall(page.read_text(encoding="utf-8"))
    ]


def test_no_documented_example_offers_an_invented_unix_root_to_copy() -> None:
    forbidden = re.compile(r"(?<![\w./-])/(" + "|".join(INVENTED_ROOTS) + r")/")
    offenders = [
        f"{name}: {line.strip()}"
        for name, block in _doc_code_blocks()
        for line in block.splitlines()
        if forbidden.search(line)
    ]

    assert not offenders, "copy-paste trap: " + "; ".join(offenders)


def test_no_powershell_block_hands_the_reader_a_posix_absolute_path() -> None:
    docs = Path(__file__).resolve().parents[1] / "docs"
    blocks = [
        (page.name, block)
        for page in sorted(docs.rglob("*.md"))
        for block in re.findall(r"```powershell\n(.*?)```", page.read_text(encoding="utf-8"), re.S)
    ]
    assert blocks, "the pages are supposed to lead with PowerShell"
    offenders = [
        f"{name}: {line.strip()}"
        for name, block in blocks
        for line in block.splitlines()
        if re.search(r'["\s=]/[A-Za-z]', line)
    ]

    assert not offenders, "PowerShell block carries a POSIX path: " + "; ".join(offenders)
