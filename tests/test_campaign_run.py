from __future__ import annotations

import json
import os
import sys
import threading
from collections import defaultdict
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from echelle_spectra import cli
from echelle_spectra.campaign_run import (
    RUN_SCHEMA,
    RunReceipt,
    new_run_directory,
    sha256_file,
    target_runs_root,
)
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


def test_two_drives_sharing_a_folder_name_get_receipt_dirs_that_tell_them_apart(
    tmp_path: Path,
) -> None:
    """`\\data` is the folder name every drive uses, so it cannot be the name."""

    runs = tmp_path / "runs"
    first = tmp_path / "drive-a" / "data"
    second = tmp_path / "drive-b" / "data"
    for source in (first, second):
        source.mkdir(parents=True)

    dated_a = new_run_directory(runs, first, volume_label="NIFS-A")
    dated_b = new_run_directory(runs, second, volume_label="NIFS-B")
    assert dated_a.name.endswith("-nifs-a-data")
    assert dated_b.name.endswith("-nifs-b-data")

    tree_a = target_runs_root(runs, first, volume_label="NIFS-A")
    tree_b = target_runs_root(runs, second, volume_label="NIFS-B")
    assert tree_a.name.startswith("nifs-a-data-")
    assert tree_b.name.startswith("nifs-b-data-")
    assert tree_a != tree_b

    # A label that only repeats the folder is not doubled into `data-data`, and
    # a target with no label at all keeps the name it always had.
    assert new_run_directory(runs, first, volume_label="data").name.endswith("-data")
    assert target_runs_root(runs, first).name.startswith("data-")


def test_multi_drive_receipt_trees_are_named_for_the_drive_that_was_processed(
    tmp_path: Path, fake_spectrocube
) -> None:
    drives = [tmp_path / "drive-a" / "data", tmp_path / "drive-b" / "data"]
    for source in drives:
        source.mkdir(parents=True)
        _files(source, "a.SIF")
    runs = tmp_path / "runs"
    calls: list[str] = []

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=_successful_export(calls),
        ),
        pytest.raises(SystemExit) as result,
    ):
        main(
            [
                *(str(source) for source in drives),
                "-o",
                str(tmp_path / "cubes"),
                "--runs-dir",
                str(runs),
                "--volume-label",
                "NIFS-A",
                "--volume-label",
                "NIFS-B",
            ]
        )

    assert result.value.code == 0
    trees = sorted(path.name for path in runs.iterdir())
    assert len(trees) == 2
    assert trees[0].startswith("nifs-a-data-")
    assert trees[1].startswith("nifs-b-data-")


def _exported_receipt(tmp_path: Path) -> tuple[RunReceipt, Path, Path]:
    """Return a receipt holding one export record, with its source and output."""
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
    receipt.append(
        receipt.identity_for(source),
        output,
        status="exported",
        started_at="2026-08-13T00:00:00Z",
        finished_at="2026-08-13T00:00:01Z",
        duration_s=1.0,
    )
    return receipt, source, output


def test_resume_trusts_recorded_identity_by_stat_and_not_by_digest(
    tmp_path: Path, monkeypatch
) -> None:
    """Deliberate replacement for the old always-digest resume check.

    This test used to pin the opposite rule: a completed output counted only
    after the source *and* the output were digested again.  That is what made a
    resume cost a full re-read of every finished file -- 885 MB per shot -- so
    the resume path now trusts the stat identity it recorded and nothing is
    read here at all.  A resized output is still refused; a same-size edit is
    accepted, and re-proving digests belongs to the catalog scrub verb.
    """

    receipt, source, output = _exported_receipt(tmp_path)
    identity = receipt.identity_for(source)

    def refuse(path: Path, **_kwargs) -> str:
        raise AssertionError(f"resume digested {path}")

    monkeypatch.setattr("echelle_spectra.campaign_run.sha256_file", refuse)
    assert receipt.resumable_record(identity, output) is not None

    output.write_text("tampered cube", encoding="utf-8")
    assert receipt.resumable_record(identity, output) is None
    assert receipt.has_export_record(identity)

    output.write_text("CUBE", encoding="utf-8")
    assert receipt.resumable_record(identity, output) is not None


def test_record_without_source_mtime_resumes_on_recorded_size_alone(tmp_path: Path) -> None:
    """Receipts written before the mtime field still resume, by size.

    The sources are read-only campaign mounts, so size alone is the identity
    those older records can offer, and refusing them would mean re-converting
    every drive that was processed before this change.
    """

    receipt, source, output = _exported_receipt(tmp_path)
    for record in receipt._records:
        record.pop("source_mtime_ns")
    identity = receipt.identity_for(source)
    assert receipt.resumable_record(identity, output) is not None

    source.write_text("source grew", encoding="utf-8")
    assert receipt.resumable_record(receipt.identity_for(source), output) is None


def test_changed_source_is_not_resumed_whether_size_or_mtime_moved(tmp_path: Path) -> None:
    """A source that is no longer the recorded file gets converted again."""

    receipt, source, output = _exported_receipt(tmp_path)
    original = source.stat()

    source.write_text("a longer source", encoding="utf-8")
    assert receipt.resumable_record(receipt.identity_for(source), output) is None

    source.write_text("source", encoding="utf-8")
    os.utime(source, ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000_000))
    touched = receipt.identity_for(source)
    assert touched.size_bytes == original.st_size
    assert receipt.resumable_record(touched, output) is None


def test_missing_or_resized_output_is_not_resumed(tmp_path: Path) -> None:
    """An output that is gone, or no longer its recorded size, is not a finished file."""

    receipt, source, output = _exported_receipt(tmp_path)
    identity = receipt.identity_for(source)

    output.write_text("cube with more in it", encoding="utf-8")
    assert receipt.resumable_record(identity, output) is None

    output.unlink()
    assert receipt.resumable_record(identity, output) is None


def test_export_record_carries_stat_identity_and_both_digests(tmp_path: Path) -> None:
    """Recording duty is unchanged: an exported file is still digested on both sides."""

    receipt, source, _output = _exported_receipt(tmp_path)
    record = receipt._records[-1]
    assert record["source_size_bytes"] == source.stat().st_size
    assert record["source_mtime_ns"] == source.stat().st_mtime_ns
    assert len(record["source_sha256"]) == 64
    assert len(record["output_sha256"]) == 64
    assert record["source_sha256"] == sha256_file(source)


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


def test_resumed_batch_skips_finished_sources_without_reading_them(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    """The whole point: a skipped source is never opened, so resume is O(remaining)."""

    source_root = tmp_path / "source"
    output_root = tmp_path / "output"
    source_root.mkdir()
    output_root.mkdir()
    _files(source_root, "a.SIF", "b.SIF")
    done_source = source_root / "a.SIF"
    done_output = output_root / "a_spectrocube.nc"
    done_output.write_text("cube from a.SIF", encoding="utf-8")
    receipt = RunReceipt.create(
        tmp_path / "run",
        source_root=source_root,
        output_root=output_root,
        pattern="*.SIF",
        volume_label="fixture",
        snapshot_id="20260813_cmos",
        expected_files=2,
    )
    receipt.append(
        receipt.identity_for(done_source),
        done_output,
        status="exported",
        started_at="2026-08-13T00:00:00Z",
        finished_at="2026-08-13T00:00:01Z",
        duration_s=1.0,
    )
    receipt.finish("interrupted")
    recorded_digest = receipt._records[-1]["source_sha256"]

    from echelle_spectra import campaign_run

    real_digest = campaign_run.sha256_file

    def digest_unless_skipped(path: Path, **kwargs) -> str:
        if Path(path).resolve() == done_source.resolve():
            raise AssertionError(f"resume re-read the finished source {path}")
        return real_digest(Path(path), **kwargs)

    calls: list[str] = []
    with (
        patch("echelle_spectra.campaign_run.sha256_file", side_effect=digest_unless_skipped),
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=_successful_export(calls),
        ),
        pytest.raises(SystemExit) as result,
    ):
        main([str(source_root), "-o", str(output_root), "--run-dir", str(receipt.directory)])

    assert result.value.code == 0
    assert calls == ["b.SIF"]
    resumed = RunReceipt.load(receipt.directory)
    assert resumed.resume_trust == "stat"
    assert 'resume_trust = "stat"' in resumed.manifest_path.read_text(encoding="utf-8")
    skipped = next(
        record for record in resumed._records if record["status"] == "skipped"
    )
    assert "recorded identity" in skipped["reason"]
    # The digest is carried from the record that proved it, not taken again.
    assert skipped["source_sha256"] == recorded_digest

    shown = capsys.readouterr().out
    assert shown.count("Resumed 1 file(s) on recorded identity (size+mtime)") == 1
    assert "echelle catalog verify" in shown


def test_run_that_resumed_nothing_claims_no_recorded_identity(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    """A run that converted everything it accounted for must not read as a resume."""

    source_root = tmp_path / "source"
    source_root.mkdir()
    _files(source_root, "a.SIF", "b.SIF")
    output_root = tmp_path / "cubes"
    runs = tmp_path / "runs"
    calls: list[str] = []

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=_successful_export(calls),
        ),
        pytest.raises(SystemExit) as result,
    ):
        main([str(source_root), "-o", str(output_root), "--runs-dir", str(runs)])

    assert result.value.code == 0
    assert calls == ["a.SIF", "b.SIF"]
    receipt = RunReceipt.load(next(runs.iterdir()))
    assert receipt.resume_trust == ""
    assert "resume_trust" not in receipt.manifest_path.read_text(encoding="utf-8")
    assert "Resumed" not in capsys.readouterr().out


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
    drift = tmp_path / "drift.json"
    drift.write_text(
        json.dumps(
            {
                "schema": "echelle-drift-evidence/v1",
                "verdict": "aligned",
                "snapshot_ids": ["20240101_cmos", "20250101_cmos"],
            }
        ),
        encoding="utf-8",
    )

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
                "--drift-verdict",
                str(drift),
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
