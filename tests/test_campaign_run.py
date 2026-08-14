from __future__ import annotations

import json
import sys
import threading
from collections import defaultdict
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from echelle_spectra import cli
from echelle_spectra.campaign_run import RUN_SCHEMA, RunReceipt, sha256_file
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot
from echelle_spectra.spectrocube_cli import ExportResult, _export_one, main


@pytest.fixture()
def fake_spectrocube():
    with patch.dict(sys.modules, {"spectrocube": MagicMock()}):
        yield


def _files(root: Path, *names: str) -> None:
    for name in names:
        (root / name).write_text(f"raw {name}", encoding="utf-8")


def _successful_export(calls: list[str]):
    def export(source: Path, output: Path, **_kwargs) -> ExportResult:
        calls.append(source.name)
        output.write_text(f"cube from {source.name}", encoding="utf-8")
        return ExportResult("exported")

    return export


def _export_kwargs() -> dict:
    return {
        "units": "counts",
        "camera": "CMOS",
        "calibration_dir": None,
        "calibration_files": {},
        "instrument_id": "echelle",
        "wavelength_medium": "air",
        "wavelength_min_nm": None,
        "wavelength_max_nm": None,
        "calibration_source": None,
        "drop_nonfinite_columns": True,
        "extra_attrs": {},
        "overwrite": False,
        "dry_run": False,
        "verbose": False,
    }


def _registry_fixture(tmp_path: Path) -> tuple[Path, Path]:
    snapshots = tmp_path / "calibrations"
    snapshot_ids = ("20240101_cmos", "20250101_cmos")
    bounds = ((100000, 199999), (200000, 299999))
    for snapshot_id, (shot_from, shot_to) in zip(snapshot_ids, bounds):
        sources = tmp_path / f"{snapshot_id}-sources"
        sources.mkdir()
        files = {}
        for role in ROLE_FILENAMES:
            source = sources / f"{role}.dat"
            source.write_text(f"{snapshot_id}/{role}", encoding="utf-8")
            files[role] = source
        create_snapshot(
            snapshots,
            snapshot_id=snapshot_id,
            detector="cmos",
            files=files,
            lamps=("ThAr",),
            validity={"shot_from": shot_from, "shot_to": shot_to},
        )
    registry = tmp_path / "calibration_registry.toml"
    registry.write_text(
        """schema = "echelle-calibration-registry/v1"

[[epochs]]
snapshot_id = "20240101_cmos"

[[epochs]]
snapshot_id = "20250101_cmos"
""",
        encoding="utf-8",
    )
    return registry, snapshots


def test_export_publishes_output_atomically(tmp_path: Path) -> None:
    source = tmp_path / "shot.SIF"
    output = tmp_path / "shot.nc"
    source.write_text("raw", encoding="utf-8")

    def write_cube(_spectrum, path: str, **_kwargs) -> None:
        Path(path).write_text("complete cube", encoding="utf-8")

    with (
        patch("echelle_spectra.tools.loader.load_spectrum", return_value=object()),
        patch(
            "echelle_spectra.tools.spectrocube_export.export_spectrocube", side_effect=write_cube
        ),
    ):
        result = _export_one(source, output, **_export_kwargs())
    assert result.status == "exported"
    assert output.read_text(encoding="utf-8") == "complete cube"
    assert not list(tmp_path.glob(".*.tmp"))


def test_interrupted_export_removes_unpublished_temporary_output(tmp_path: Path) -> None:
    source = tmp_path / "shot.SIF"
    output = tmp_path / "shot.nc"
    source.write_text("raw", encoding="utf-8")

    def interrupt(_spectrum, path: str, **_kwargs) -> None:
        Path(path).write_text("partial cube", encoding="utf-8")
        raise KeyboardInterrupt

    with (
        patch("echelle_spectra.tools.loader.load_spectrum", return_value=object()),
        patch("echelle_spectra.tools.spectrocube_export.export_spectrocube", side_effect=interrupt),
        pytest.raises(KeyboardInterrupt),
    ):
        _export_one(source, output, **_export_kwargs())
    assert not output.exists()
    assert not list(tmp_path.glob(".*.tmp"))


def test_interrupted_batch_resumes_only_unfinished_sources(
    tmp_path: Path, fake_spectrocube
) -> None:
    source = tmp_path / "source"
    output = tmp_path / "cubes"
    runs = tmp_path / "runs"
    source.mkdir()
    _files(source, "a.SIF", "b.SIF", "c.SIF")
    first_calls: list[str] = []

    def interrupted_export(sif: Path, nc_out: Path, **_kwargs) -> ExportResult:
        first_calls.append(sif.name)
        if sif.name == "b.SIF":
            raise KeyboardInterrupt
        nc_out.write_text("completed a", encoding="utf-8")
        return ExportResult("exported")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=interrupted_export),
        pytest.raises(SystemExit) as interrupted,
    ):
        main(
            [
                str(source),
                "-o",
                str(output),
                "--runs-dir",
                str(runs),
                "--snapshot-id",
                "20260813_cmos",
                "--volume-label",
                "NIFS-A",
            ]
        )
    assert interrupted.value.code == 130
    assert first_calls == ["a.SIF", "b.SIF"]
    completed = output / "a_spectrocube.nc"
    completed_digest = sha256_file(completed)

    second_calls: list[str] = []
    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=_successful_export(second_calls),
        ),
        pytest.raises(SystemExit) as resumed,
    ):
        main(
            [
                str(source),
                "-o",
                str(output),
                "--runs-dir",
                str(runs),
                "--snapshot-id",
                "20260813_cmos",
                "--volume-label",
                "NIFS-A",
            ]
        )
    assert resumed.value.code == 0
    assert second_calls == ["b.SIF", "c.SIF"]
    assert sha256_file(completed) == completed_digest

    run_dir = next(runs.iterdir())
    receipt = RunReceipt.load(run_dir)
    assert receipt.state == "completed"
    assert receipt.snapshot_id == "20260813_cmos"
    assert receipt.volume_label == "NIFS-A"
    assert receipt.counts() == {"exported": 2, "failed": 0, "interrupted": 0, "skipped": 1}
    records = [json.loads(line) for line in receipt.records_path.read_text().splitlines()]
    assert [record["status"] for record in records] == [
        "exported",
        "interrupted",
        "skipped",
        "exported",
        "exported",
    ]
    assert records[0]["source_sha256"]
    assert records[0]["output_sha256"]


def test_batch_continues_after_failure_and_accounts_for_every_source(
    tmp_path: Path, fake_spectrocube
) -> None:
    source = tmp_path / "source"
    output = tmp_path / "cubes"
    runs = tmp_path / "runs"
    source.mkdir()
    _files(source, "a.SIF", "b.SIF", "c.SIF")
    calls: list[str] = []

    def export(sif: Path, nc_out: Path, **_kwargs) -> ExportResult:
        calls.append(sif.name)
        if sif.name == "b.SIF":
            return ExportResult("failed", "synthetic detector failure")
        nc_out.write_text(f"cube {sif.name}", encoding="utf-8")
        return ExportResult("exported")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=export),
        pytest.raises(SystemExit) as failed,
    ):
        main([str(source), "-o", str(output), "--runs-dir", str(runs)])
    assert failed.value.code == 1
    assert calls == ["a.SIF", "b.SIF", "c.SIF"]

    receipt = RunReceipt.load(next(runs.iterdir()))
    assert receipt.state == "partial"
    assert receipt.counts() == {"exported": 2, "failed": 1, "interrupted": 0, "skipped": 0}
    failure = next(record for record in receipt._records if record["status"] == "failed")
    assert failure["reason"] == "synthetic detector failure"


def test_multi_drive_workers_run_concurrently_and_reconcile_status(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    drive_a = tmp_path / "drive-a"
    drive_b = tmp_path / "drive-b"
    drive_a.mkdir()
    drive_b.mkdir()
    _files(drive_a, "a.SIF", "b.SIF")
    _files(drive_b, "a.SIF", "b.SIF")
    outputs = tmp_path / "cubes"
    runs = tmp_path / "runs"
    first_reads = threading.Barrier(2)
    calls: dict[str, list[str]] = defaultdict(list)

    def export(sif: Path, nc_out: Path, **_kwargs) -> ExportResult:
        calls[sif.parent.name].append(sif.name)
        if sif.name == "a.SIF":
            first_reads.wait(timeout=2)
        if sif.parent == drive_b and sif.name == "a.SIF":
            return ExportResult("failed", "synthetic drive-b read failure")
        nc_out.write_text(f"cube {sif.parent.name}/{sif.name}", encoding="utf-8")
        return ExportResult("exported")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=export),
        pytest.raises(SystemExit) as result,
    ):
        main(
            [
                str(drive_a),
                str(drive_b),
                "-o",
                str(outputs),
                "--runs-dir",
                str(runs),
                "--volume-label",
                "NIFS-A",
                "--volume-label",
                "NIFS-B",
            ]
        )

    assert result.value.code == 1
    assert calls == {
        "drive-a": ["a.SIF", "b.SIF"],
        "drive-b": ["a.SIF", "b.SIF"],
    }
    assert sorted(path.name for path in (outputs / "drive-a").glob("*.nc")) == [
        "a_spectrocube.nc",
        "b_spectrocube.nc",
    ]
    assert sorted(path.name for path in (outputs / "drive-b").glob("*.nc")) == [
        "b_spectrocube.nc"
    ]

    receipts = [RunReceipt.load(path.parent) for path in runs.rglob("run.toml")]
    assert len(receipts) == 2
    by_volume = {receipt.volume_label: receipt for receipt in receipts}
    assert by_volume["NIFS-A"].state == "completed"
    assert by_volume["NIFS-A"].counts()["exported"] == 2
    assert by_volume["NIFS-B"].state == "partial"
    assert by_volume["NIFS-B"].counts()["exported"] == 1
    assert by_volume["NIFS-B"].counts()["failed"] == 1

    capsys.readouterr()
    assert cli.main(["status", "--runs", str(runs)]) == 0
    shown = capsys.readouterr().out
    assert "targets:   2 independent source(s)" in shown
    assert "combined:  4/4" in shown
    assert "NIFS-A: 2/2 [completed]" in shown
    assert "NIFS-B: 2/2 [partial]" in shown


def test_changed_output_is_not_accepted_as_completed(tmp_path: Path) -> None:
    source_root = tmp_path / "source"
    output_root = tmp_path / "output"
    source_root.mkdir()
    output_root.mkdir()
    source = source_root / "a.SIF"
    output = output_root / "a.nc"
    source.write_text("source", encoding="utf-8")
    output.write_text("cube", encoding="utf-8")
    receipt = RunReceipt.create(
        tmp_path / "run",
        source_root=source_root,
        output_root=output_root,
        pattern="*.SIF",
        volume_label="fixture",
        snapshot_id="20260813_cmos",
        expected_files=1,
    )
    identity = receipt.identity_for(source)
    receipt.append(
        identity,
        output,
        status="exported",
        started_at="2026-08-13T00:00:00Z",
        finished_at="2026-08-13T00:00:01Z",
        duration_s=1.0,
    )
    assert receipt.completed_output_is_valid(identity, output)
    output.write_text("tampered", encoding="utf-8")
    assert not receipt.completed_output_is_valid(identity, output)
    assert receipt.has_export_record(identity)


def test_resume_refuses_changed_completed_output(tmp_path: Path, fake_spectrocube) -> None:
    source_root = tmp_path / "source"
    output_root = tmp_path / "output"
    source_root.mkdir()
    output_root.mkdir()
    source = source_root / "a.SIF"
    output = output_root / "a_spectrocube.nc"
    source.write_text("source", encoding="utf-8")
    output.write_text("original cube", encoding="utf-8")
    receipt = RunReceipt.create(
        tmp_path / "run",
        source_root=source_root,
        output_root=output_root,
        pattern="*.SIF",
        volume_label="fixture",
        snapshot_id="20260813_cmos",
        expected_files=1,
    )
    receipt.append(
        receipt.identity_for(source),
        output,
        status="exported",
        started_at="2026-08-13T00:00:00Z",
        finished_at="2026-08-13T00:00:01Z",
        duration_s=1.0,
    )
    receipt.finish("interrupted")
    output.write_text("changed cube", encoding="utf-8")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one") as export,
        pytest.raises(SystemExit) as failed,
    ):
        main(
            [
                str(source_root),
                "-o",
                str(output_root),
                "--run-dir",
                str(receipt.directory),
            ]
        )
    assert failed.value.code == 1
    export.assert_not_called()
    resumed = RunReceipt.load(receipt.directory)
    assert resumed.state == "partial"
    assert resumed.counts()["failed"] == 1
    assert "output changed" in resumed._records[-1]["reason"]


def test_status_reports_latest_receipt(tmp_path: Path, capsys) -> None:
    runs = tmp_path / "runs"
    receipt = RunReceipt.create(
        runs / "2026-08-13_00-00-fixture",
        source_root=tmp_path / "source",
        output_root=tmp_path / "output",
        pattern="*.SIF",
        volume_label="fixture",
        snapshot_id="20260813_cmos",
        expected_files=3,
    )
    receipt.finish("interrupted")
    assert cli.main(["status", "--runs", str(runs)]) == 0
    shown = capsys.readouterr().out
    assert "1 receipt(s)" in shown
    assert "[interrupted]" in shown
    assert "0/3" in shown
    assert "20260813_cmos" in shown


def test_run_manifest_uses_versioned_schema(tmp_path: Path) -> None:
    receipt = RunReceipt.create(
        tmp_path / "run",
        source_root=tmp_path / "source",
        output_root=tmp_path / "output",
        pattern="*.SIF",
        volume_label="fixture",
        snapshot_id="unassigned",
        expected_files=0,
    )
    assert f'schema = "{RUN_SCHEMA}"' in receipt.manifest_path.read_text(encoding="utf-8")


def test_registry_batch_resolves_and_records_one_snapshot_per_source(
    tmp_path: Path, fake_spectrocube
) -> None:
    registry, snapshots = _registry_fixture(tmp_path)
    source = tmp_path / "source"
    source.mkdir()
    _files(source, "199999_Echelle.SIF", "200000_Echelle.SIF", "200001_Echelle.SIF")
    output = tmp_path / "cubes"
    runs = tmp_path / "runs"
    calibrations = {"20240101_cmos": object(), "20250101_cmos": object()}
    exports: list[tuple[str, str]] = []

    def build(folder: Path, _camera: str, **_kwargs):
        return calibrations[Path(folder).name]

    def export(sif: Path, nc_out: Path, **kwargs) -> ExportResult:
        exports.append((sif.name, kwargs["snapshot"].snapshot_id))
        assert kwargs["calibration"] is calibrations[kwargs["snapshot"].snapshot_id]
        assert len(kwargs["extra_attrs"]["calibration_registry_sha256"]) == 64
        nc_out.write_text("cube", encoding="utf-8")
        return ExportResult("exported")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", side_effect=build) as builder,
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=export),
        pytest.raises(SystemExit) as result,
    ):
        main(
            [
                str(source),
                "-o",
                str(output),
                "--runs-dir",
                str(runs),
                "--registry",
                str(registry),
                "--calibrations",
                str(snapshots),
            ]
        )
    assert result.value.code == 0
    assert exports == [
        ("199999_Echelle.SIF", "20240101_cmos"),
        ("200000_Echelle.SIF", "20250101_cmos"),
        ("200001_Echelle.SIF", "20250101_cmos"),
    ]
    assert builder.call_count == 2
    receipt = RunReceipt.load(next(runs.iterdir()))
    assert receipt.snapshot_id == "per-source-registry"
    assert [record["snapshot_id"] for record in receipt._records] == [
        "20240101_cmos",
        "20250101_cmos",
        "20250101_cmos",
    ]
