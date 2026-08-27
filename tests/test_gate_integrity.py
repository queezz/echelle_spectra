"""Packet F3 — every registry-backed run is authorized, and its receipt says how.

The gate is exercised through the shipped ``echelle process`` entry point. Only
the SIF reader and the exporter are replaced: the fake exporter writes a real
NetCDF file carrying exactly the attributes the CLI handed it, so the sample
marking is checked on a saved cube rather than on a call argument.
"""

from __future__ import annotations

import json
import re
import sys
from collections.abc import Iterator
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import xarray as xr

from echelle_spectra import campaign_tools_cli, cli
from echelle_spectra.campaign_run import (
    GATE_SAMPLE,
    GATE_UNGATED,
    GATE_UNRECORDED,
    GATE_VERDICT,
    RUN_SCHEMA,
    RunReceipt,
)
from echelle_spectra.drift import (
    DRIFT_SCHEMA,
    DriftError,
    create_refinement_snapshot,
    require_sampled_verdict,
    write_drift_evidence,
)
from echelle_spectra.lhd_text import TIMING_ATTR_SOURCES
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot
from echelle_spectra.spectrocube_cli import ExportResult, main
from echelle_spectra.spectrocube_config import export_config_from_toml

BASE_ID = "20250101_cmos"
REFINED_ID = "20250101_cmos-r1"

DOC_GUIDES = (
    "operator-cheat-sheet.md",
    "campaign-runs.md",
    "calibration-epoch-registry.md",
    # The two routing pages carry copy-paste examples too, and a copied example
    # is exactly how a wrong invocation reaches an operator.
    "usage-overview.md",
    "calibration-to-cube.md",
)
UMBRELLA_COMMANDS = {
    "status",
    "snapshot",
    "process",
    "catalog",
    "txt",
    "recal-cube",
    "drift",
    "inventory",
    "web",
    "historical",
}
NESTED_COMMANDS = {"snapshot", "catalog", "drift"}


@pytest.fixture()
def fake_spectrocube():
    with patch.dict(sys.modules, {"spectrocube": MagicMock()}):
        yield


def _snapshot(root: Path, snapshot_id: str, *, shot_from: int, shot_to: int) -> None:
    sources = root.parent / f"{snapshot_id}-sources"
    sources.mkdir(parents=True)
    files: dict[str, Path] = {}
    for role in ROLE_FILENAMES:
        source = sources / f"{role}.dat"
        if role == "wavelength":
            source.write_text(
                "# order from to center wavelength\n"
                + "".join(
                    f"0\t{index}\t{index + 1}\t{index}\t{410.0 + index:.6f}\tline\n"
                    for index in range(5)
                ),
                encoding="utf-8",
            )
        else:
            source.write_text(f"{snapshot_id}/{role}\n", encoding="utf-8")
        files[role] = source
    create_snapshot(
        root,
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity={"shot_from": shot_from, "shot_to": shot_to},
    )


def _epoch_registry(
    tmp_path: Path, *epochs: tuple[str, int, int]
) -> tuple[Path, Path]:
    """Write one ordered registry and every snapshot it names."""

    snapshots = tmp_path / "calibrations"
    lines = ['schema = "echelle-calibration-registry/v1"']
    for snapshot_id, shot_from, shot_to in epochs or ((BASE_ID, 100000, 299999),):
        _snapshot(snapshots, snapshot_id, shot_from=shot_from, shot_to=shot_to)
        lines.extend(["", "[[epochs]]", f'snapshot_id = "{snapshot_id}"'])
    registry = tmp_path / "calibration_registry.toml"
    registry.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return registry, snapshots


def _shots(root: Path, count: int) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for index in range(count):
        (root / f"20000{index}_Echelle.SIF").write_text(f"raw {index}", encoding="utf-8")
    return root


def _evidence(path: Path, payload: dict) -> Path:
    return write_drift_evidence(path, {"schema": DRIFT_SCHEMA, **payload})


def _cube_writer(calls: list[str]):
    """Write a real cube carrying the attributes the CLI selected for it."""

    def export(sif: Path, nc_out: Path, **kwargs) -> ExportResult:
        calls.append(sif.name)
        xr.Dataset(attrs=dict(kwargs["extra_attrs"])).to_netcdf(nc_out)
        return ExportResult("exported")

    return export


def _process(argv: list[str], calls: list[str] | None = None) -> int:
    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=_cube_writer(calls if calls is not None else []),
        ),
        pytest.raises(SystemExit) as result,
    ):
        main(argv)
    return int(result.value.code)


def _registry_argv(tmp_path: Path, source: Path, registry: Path, snapshots: Path) -> list[str]:
    return [
        str(source),
        "-o",
        str(tmp_path / "cubes"),
        "--runs-dir",
        str(tmp_path / "runs"),
        "--registry",
        str(registry),
        "--calibrations",
        str(snapshots),
    ]


# ---------------------------------------------------------------------------
# Every registry-backed path is gated
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("count", [1, 3])
def test_registry_folder_run_of_any_size_refuses_without_verdict_or_sample(
    tmp_path: Path, fake_spectrocube, capsys, count: int
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", count)

    code = _process(_registry_argv(tmp_path, source, registry, snapshots))

    assert code == 1
    assert not (tmp_path / "runs").exists()
    assert not list(tmp_path.glob("cubes/*.nc"))
    refusal = capsys.readouterr().err
    assert "requires a sampled epoch audit" in refusal
    assert "--sample 5" in refusal
    assert "echelle drift audit" in refusal
    assert "--drift-verdict drift-evidence.json" in refusal


def test_single_file_registry_run_refuses_without_a_verdict(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 1)
    sif = source / "200000_Echelle.SIF"

    code = _process(
        [
            str(sif),
            "-o",
            str(tmp_path / "one.nc"),
            "--registry",
            str(registry),
            "--calibrations",
            str(snapshots),
        ]
    )

    assert code == 1
    assert not (tmp_path / "one.nc").exists()
    refusal = capsys.readouterr().err
    assert "requires a sampled epoch audit" in refusal
    assert "echelle drift audit" in refusal


def test_config_run_without_a_registry_stays_legal_and_records_that_it_was_ungated(
    tmp_path: Path, fake_spectrocube
) -> None:
    config = tmp_path / "export.toml"
    config.write_text('[calibration]\ncamera = "CMOS"\n', encoding="utf-8")
    source = _shots(tmp_path / "shots", 2)

    code = _process(
        [
            str(source),
            "-o",
            str(tmp_path / "cubes"),
            "--runs-dir",
            str(tmp_path / "runs"),
            "--config",
            str(config),
        ]
    )

    assert code == 0
    receipt = RunReceipt.load(next((tmp_path / "runs").iterdir()))
    assert receipt.gate == GATE_UNGATED
    assert receipt.sample is False
    assert receipt.drift_evidence == ""
    manifest = receipt.manifest_path.read_text(encoding="utf-8")
    assert f'gate = "{GATE_UNGATED}"' in manifest


# ---------------------------------------------------------------------------
# The legal first run
# ---------------------------------------------------------------------------


def test_sample_run_marks_receipt_and_cubes_and_still_leaves_bulk_refused(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 4)
    argv = _registry_argv(tmp_path, source, registry, snapshots)
    calls: list[str] = []

    code = _process([*argv, "--sample", "2"], calls)

    assert code == 0
    assert calls == ["200000_Echelle.SIF", "200001_Echelle.SIF"]
    cubes = sorted((tmp_path / "cubes").glob("*.nc"))
    assert [path.name for path in cubes] == [
        "200000_Echelle_spectrocube.nc",
        "200001_Echelle_spectrocube.nc",
    ]
    for cube in cubes:
        with xr.open_dataset(cube) as opened:
            assert opened.attrs["drift_sample"] == 1
            assert opened.attrs["calibration_registry_epoch_position"] == 1

    receipt = RunReceipt.load(next((tmp_path / "runs").iterdir()))
    assert receipt.gate == GATE_SAMPLE
    assert receipt.sample is True
    assert receipt.sample_files == 2
    assert receipt.expected_files == 2
    manifest = receipt.manifest_path.read_text(encoding="utf-8")
    assert f'schema = "{RUN_SCHEMA}"' in manifest
    assert "sample = true" in manifest
    assert "sample_files = 2" in manifest
    capsys.readouterr()

    # A sample authorizes only itself: the whole drive still needs a verdict.
    assert _process(argv) == 1
    assert "requires a sampled epoch audit" in capsys.readouterr().err


def test_sample_cannot_be_combined_with_a_verdict_or_used_without_a_registry(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 2)
    evidence = _evidence(
        tmp_path / "aligned.json", {"verdict": "aligned", "snapshot_ids": [BASE_ID]}
    )
    argv = _registry_argv(tmp_path, source, registry, snapshots)

    assert _process([*argv, "--sample", "2", "--drift-verdict", str(evidence)]) == 2
    assert "cannot be combined" in capsys.readouterr().err

    assert _process([str(source), "-o", str(tmp_path / "cubes"), "--sample", "2"]) == 2
    assert "requires --registry" in capsys.readouterr().err

    assert _process([*argv, "--sample", "0"]) == 2
    assert "at least one file" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# An accepted shift authorizes only the refined snapshot
# ---------------------------------------------------------------------------


def test_accepted_refinement_authorizes_only_the_refined_snapshot(tmp_path: Path) -> None:
    _, snapshots = _epoch_registry(tmp_path)
    evidence = _evidence(
        tmp_path / "drift.json",
        {
            "verdict": "shifted",
            "snapshot_ids": [BASE_ID],
            "summary": {"median_shift_px": 5.0},
        },
    )

    refined, accepted = create_refinement_snapshot(
        evidence, calibrations_root=snapshots, accepted_shift_px=5.0
    )

    assert refined.snapshot_id == REFINED_ID
    record = json.loads(accepted.read_text(encoding="utf-8"))
    assert record["snapshot_ids"] == [REFINED_ID]
    assert record["base_snapshot_ids"] == [BASE_ID]
    with pytest.raises(DriftError, match="update the registry validity"):
        require_sampled_verdict(accepted, {BASE_ID})
    assert require_sampled_verdict(accepted, {REFINED_ID})["accepted"] is True


def test_registry_still_selecting_the_condemned_base_is_refused_with_the_repoint(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 2)
    accepted = _evidence(
        tmp_path / "accepted.json",
        {
            "verdict": "shifted",
            "accepted": True,
            "accepted_snapshot_id": REFINED_ID,
            "snapshot_ids": [REFINED_ID],
            "base_snapshot_ids": [BASE_ID],
        },
    )
    argv = _registry_argv(tmp_path, source, registry, snapshots)

    code = _process([*argv, "--drift-verdict", str(accepted)])

    assert code == 1
    refusal = capsys.readouterr().err
    assert f"registry still selects {BASE_ID}" in refusal
    assert f"the accepted correction produced {REFINED_ID}" in refusal
    assert "update the registry validity to point at the refined snapshot" in refusal


def test_registry_repointed_at_the_refinement_passes_and_the_receipt_names_its_evidence(
    tmp_path: Path, fake_spectrocube
) -> None:
    registry, snapshots = _epoch_registry(tmp_path, (REFINED_ID, 100000, 299999))
    source = _shots(tmp_path / "shots", 2)
    accepted = _evidence(
        tmp_path / "accepted.json",
        {
            "verdict": "shifted",
            "accepted": True,
            "accepted_snapshot_id": REFINED_ID,
            "snapshot_ids": [REFINED_ID],
            "base_snapshot_ids": [BASE_ID],
        },
    )
    calls: list[str] = []

    code = _process(
        [*_registry_argv(tmp_path, source, registry, snapshots), "--drift-verdict", str(accepted)],
        calls,
    )

    assert code == 0
    assert len(calls) == 2
    receipt = RunReceipt.load(next((tmp_path / "runs").iterdir()))
    assert receipt.gate == GATE_VERDICT
    assert receipt.sample is False
    assert receipt.drift_evidence == str(accepted.resolve())
    assert len(receipt.drift_evidence_sha256) == 64
    assert receipt.drift_verdict == "shifted"
    for cube in (tmp_path / "cubes").glob("*.nc"):
        with xr.open_dataset(cube) as opened:
            assert "drift_sample" not in opened.attrs


def test_aligned_verdict_authorizes_the_audited_ids_and_status_shows_the_gate(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 3)
    evidence = _evidence(
        tmp_path / "aligned.json", {"verdict": "aligned", "snapshot_ids": [BASE_ID]}
    )

    code = _process(
        [*_registry_argv(tmp_path, source, registry, snapshots), "--drift-verdict", str(evidence)]
    )

    assert code == 0
    receipt = RunReceipt.load(next((tmp_path / "runs").iterdir()))
    assert receipt.gate == GATE_VERDICT
    assert receipt.drift_verdict == "aligned"
    manifest = receipt.manifest_path.read_text(encoding="utf-8")
    assert 'drift_verdict = "aligned"' in manifest
    capsys.readouterr()

    assert cli.main(["status", "--runs", str(tmp_path / "runs")]) == 0
    shown = capsys.readouterr().out
    assert f"gate:      {GATE_VERDICT} (aligned;" in shown


# ---------------------------------------------------------------------------
# Self-explanation and backward-readable receipts
# ---------------------------------------------------------------------------


def test_bare_catalog_and_drift_print_help_like_bare_snapshot(capsys) -> None:
    assert campaign_tools_cli.catalog_main([]) == 0
    catalog_help = capsys.readouterr().out
    assert "build" in catalog_help and "merge" in catalog_help

    assert campaign_tools_cli.drift_main([]) == 0
    drift_help = capsys.readouterr().out
    assert "audit" in drift_help and "refine" in drift_help

    assert cli.main(["snapshot"]) == 0
    assert cli.main(["catalog"]) == 0
    assert cli.main(["drift"]) == 0
    assert capsys.readouterr().out.count("usage:") == 3


def _documented_commands(text: str) -> Iterator[str]:
    """Yield every ``echelle`` invocation written into one guide's code blocks."""

    for block in re.findall(r"```[a-z]*\n(.*?)```", text, flags=re.S):
        joined = block.replace("`\n", " ").replace("\\\n", " ")
        for raw in joined.splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            match = re.search(r"(?:^|[\\/.])echelle(?:\.ps1)?\s+(\S.*)$", line)
            if match is not None:
                yield match.group(1).strip()


def test_every_documented_command_parses_against_the_shipped_cli(capsys) -> None:
    docs = Path(__file__).resolve().parents[1] / "docs"
    checked = 0
    for guide in DOC_GUIDES:
        for command in _documented_commands((docs / guide).read_text(encoding="utf-8")):
            tokens = command.split()
            if tokens[0].startswith("-"):
                # The umbrella's own flags, answered by `echelle` itself rather
                # than by a subcommand. They are still documented invocations,
                # so they are run rather than skipped.
                assert cli.main(tokens) == 0, f"{guide}: 'echelle {command}' failed"
                capsys.readouterr()
                checked += 1
                continue
            assert tokens[0] in UMBRELLA_COMMANDS, f"{guide}: unknown command '{command}'"
            parts = [tokens[0]]
            if tokens[0] in NESTED_COMMANDS and len(tokens) > 1 and not tokens[1].startswith("-"):
                parts.append(tokens[1])
            with pytest.raises(SystemExit) as shown:
                cli.main([*parts, "--help"])
            assert shown.value.code == 0
            rendered = capsys.readouterr().out
            for flag in re.findall(r"(?<![\w-])--[a-z][a-z-]+", command):
                assert flag in rendered, f"{guide}: 'echelle {' '.join(parts)}' has no {flag}"
            checked += 1
    assert checked >= 15


def test_the_documented_timing_config_reaches_a_cube_through_a_registry_run(
    tmp_path: Path, fake_spectrocube
) -> None:
    """The documented loop can reach `echelle txt`, using only a file it prints.

    ``trigger_delay_s`` is produced by nothing else in the loop -- no snapshot
    carries it, no audit derives it, no processing flag supplies it -- so until
    the guide printed this file the documented trip could not produce a cube
    ``echelle txt`` would accept, and the missing file was one the docs never
    mentioned.
    """

    guide = Path(__file__).resolve().parents[1] / "docs" / "operator-cheat-sheet.md"
    blocks = [
        block
        for block in re.findall(r"```toml\n(.*?)```", guide.read_text(encoding="utf-8"), re.S)
        if "trigger_delay_s" in block
    ]
    assert len(blocks) == 1, "the guide must print exactly one export-timing example"
    config = tmp_path / "export-timing.toml"
    config.write_text(blocks[0], encoding="utf-8")

    # It may name nothing the registry owns, or the run refuses it as a second
    # authority over one export.
    settings = export_config_from_toml(config)
    assert settings["camera"] is None
    assert settings["calibration_dir"] is None
    assert not settings["calibration_files"]
    assert "trigger_delay_s" in settings["extra_attrs"]

    registry, snapshots = _epoch_registry(tmp_path)
    source = _shots(tmp_path / "shots", 1)

    code = _process(
        [
            *_registry_argv(tmp_path, source, registry, snapshots),
            "--config",
            str(config),
            "--sample",
            "1",
        ]
    )

    assert code == 0
    cube = next((tmp_path / "cubes").glob("*.nc"))
    with xr.open_dataset(cube) as opened:
        assert float(opened.attrs["trigger_delay_s"]) > 0.0
        assert opened.attrs["frame_time_formula"] == "trigger_delay_s + frame * frame_interval_s"
    # And this really is the attribute the frozen header has no other source for.
    assert "[metadata] trigger_delay_s" in TIMING_ATTR_SOURCES["trigger_delay_s"]


def test_receipt_written_before_the_gate_still_loads_and_reports_its_silence(
    tmp_path: Path, capsys
) -> None:
    runs = tmp_path / "runs" / "2026-08-13_00-00-00-shots"
    runs.mkdir(parents=True)
    (runs / "run.toml").write_text(
        f"""schema = "{RUN_SCHEMA}"

[run]
id = "2026-08-13_00-00-00-shots"
state = "completed"
created_at = "2026-08-13T00:00:00.000Z"
updated_at = "2026-08-13T00:00:10.000Z"
source_root = "{(tmp_path / 'shots').as_posix()}"
output_root = "{(tmp_path / 'cubes').as_posix()}"
pattern = "*.SIF"
volume_label = "NIFS-A"
snapshot_id = "per-source-registry"
expected_files = 2

[counts]
exported = 2
failed = 0
interrupted = 0
skipped = 0
""",
        encoding="utf-8",
    )

    receipt = RunReceipt.load(runs)
    assert receipt.gate == GATE_UNRECORDED
    assert receipt.sample is False
    assert receipt.expected_files == 2

    assert cli.main(["status", "--runs", str(tmp_path / "runs")]) == 0
    assert "gate:      unrecorded" in capsys.readouterr().out
