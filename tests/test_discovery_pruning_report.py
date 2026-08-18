"""A skipped folder is announced, not merely skipped.

Recursive discovery prunes any folder named ``calibrations`` and any folder
holding a ``snapshot.toml``.  That is right -- lamp frames are not science
shots -- but the pruning was only ever described when discovery found nothing
at all.  The bench, however, writes its snapshot into whichever folder it was
launched at, and ``echelle-calib <folder>`` accepts any folder: a day folder
that also holds science SIFs is a legal target.  The next ``echelle process``
over that drive then dropped the whole day, and the console header, the run
receipt, and the campaign page built from it all read as a complete run.

These tests pin the fix: whenever pruning removed a folder that exists, the
batch header names it and the receipt carries it -- and when nothing was
pruned, nothing extra is said.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import tomllib

from echelle_spectra.spectrocube_cli import (
    PRUNED_LISTING_LIMIT,
    ExportResult,
    _discover_science_files,
    main,
)


@pytest.fixture()
def fake_spectrocube():
    with patch.dict(sys.modules, {"spectrocube": MagicMock()}):
        yield


def _shots(folder: Path, *names: str) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        (folder / f"{name}_Echelle.SIF").write_text(f"raw {name}", encoding="utf-8")
    return folder


def _snapshot_marker(folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "snapshot.toml").write_text('id = "20250101_cmos"\n', encoding="utf-8")
    return folder


def _process(argv: list[str]) -> int:
    """Run the shipped entry point with only the exporter replaced."""

    def export(sif: Path, nc_out: Path, **kwargs: object) -> ExportResult:
        # The real exporter writes nothing during a dry run, and neither may
        # this one: the destination folder does not exist yet.
        if kwargs.get("dry_run"):
            return ExportResult("dry-run", "dry run")
        nc_out.write_text(f"cube for {sif.name}", encoding="utf-8")
        return ExportResult("exported")

    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=export),
        pytest.raises(SystemExit) as result,
    ):
        main(argv)
    return int(result.value.code)


def _run_manifest(runs_root: Path) -> dict:
    manifests = sorted(runs_root.rglob("run.toml"))
    assert len(manifests) == 1, f"expected exactly one receipt, found {manifests}"
    with manifests[0].open("rb") as stream:
        return tomllib.load(stream)["run"]


def _batch_argv(tmp_path: Path, source: Path) -> list[str]:
    return [
        str(source),
        "-o",
        str(tmp_path / "cubes"),
        "--runs-dir",
        str(tmp_path / "runs"),
    ]


# ---------------------------------------------------------------------------
# A day folder the bench wrote a snapshot into
# ---------------------------------------------------------------------------


def test_a_day_folder_holding_a_snapshot_is_skipped_and_said_so(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    # The trap exactly as the bench produces it: real shots, and a snapshot.toml
    # written beside them because `echelle-calib` was launched at this folder.
    _shots(_snapshot_marker(source / "20190207"), "193780", "193781")

    discovery = _discover_science_files(source, "*.SIF")

    assert [path.name for path in discovery.files] == ["193778_Echelle.SIF"]
    assert [path.name for path in discovery.pruned] == ["20190207"]

    assert _process(_batch_argv(tmp_path, source)) == 0

    header = capsys.readouterr().out
    assert (
        "Skipped:     1 calibration folder(s) not searched "
        "(named 'calibrations', or holding a snapshot.toml):" in header
    )
    assert "               20190207" in header
    assert _run_manifest(tmp_path / "runs")["pruned_dirs"] == ["20190207"]


def test_a_plain_calibrations_folder_is_skipped_and_said_so(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    _shots(source / "calibrations", "Ne_lamp")

    assert _process(_batch_argv(tmp_path, source)) == 0

    header = capsys.readouterr().out
    assert "Skipped:     1 calibration folder(s) not searched" in header
    assert "               calibrations" in header
    assert _run_manifest(tmp_path / "runs")["pruned_dirs"] == ["calibrations"]


def test_a_drive_with_nothing_pruned_says_nothing_extra(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    _shots(source / "20190207", "193780")

    assert _process(_batch_argv(tmp_path, source)) == 0

    assert "Skipped:" not in capsys.readouterr().out
    # An untouched manifest is the point: a run that pruned nothing carries no
    # pruning key at all, so an auditor reading one is reading a real skip.
    assert "pruned_dirs" not in _run_manifest(tmp_path / "runs")


# ---------------------------------------------------------------------------
# Bounded on the console, complete in the receipt
# ---------------------------------------------------------------------------


def test_a_long_skip_list_is_bounded_on_screen_and_whole_in_the_receipt(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    pruned = [f"2019030{index}" for index in range(1, PRUNED_LISTING_LIMIT + 3)]
    for name in pruned:
        _shots(_snapshot_marker(source / name), f"shot_{name}")

    assert _process(_batch_argv(tmp_path, source)) == 0

    header = capsys.readouterr().out
    assert f"Skipped:     {len(pruned)} calibration folder(s) not searched" in header
    for name in pruned[:PRUNED_LISTING_LIMIT]:
        assert f"               {name}" in header
    for name in pruned[PRUNED_LISTING_LIMIT:]:
        assert f"               {name}" not in header
    assert f"             ... and {len(pruned) - PRUNED_LISTING_LIMIT} more" in header

    assert _run_manifest(tmp_path / "runs")["pruned_dirs"] == pruned


def test_a_nested_skip_is_reported_at_the_folder_that_was_actually_dropped(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    # Everything under a pruned folder went with it; naming the children too
    # would bury the one fact the operator needs.
    snapshot = _snapshot_marker(source / "20190207")
    _shots(snapshot / "sources", "193780")
    _shots(snapshot / "sources" / "deeper", "193781")

    assert _process(_batch_argv(tmp_path, source)) == 0

    header = capsys.readouterr().out
    assert "Skipped:     1 calibration folder(s) not searched" in header
    assert "sources" not in header
    assert _run_manifest(tmp_path / "runs")["pruned_dirs"] == ["20190207"]


def test_a_dry_run_lists_the_files_it_found_without_being_asked_to_be_verbose(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    """The documented step says to inspect the listed sources; so list them.

    ``--dry-run`` printed a count and, behind ``--verbose``, the names. The
    guide's "start with --dry-run and inspect the listed sources" therefore
    described a listing the plain command never produced.
    """

    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    _shots(source / "20190207", "193780")
    _shots(_snapshot_marker(source / "20190208"), "193790")

    assert _process([*_batch_argv(tmp_path, source), "--dry-run"]) == 0

    out = capsys.readouterr().out
    assert "Files:       2 (DRY RUN)" in out
    assert "Would convert:" in out
    assert "  20190206/193778_Echelle.SIF" in out
    assert "  20190207/193780_Echelle.SIF" in out
    # A pruned day is skipped, not silently listed as convertible.
    assert "193790" not in out.split("Would convert:")[1]
    # And a dry run still writes no receipt at all.
    assert not (tmp_path / "runs").exists()


def test_an_operators_own_path_pattern_prunes_nothing_and_reports_nothing(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    source = tmp_path / "drive"
    _shots(source / "calibrations", "Ne_lamp")
    _shots(source / "20190206", "193778")

    assert (
        _process([*_batch_argv(tmp_path, source), "--pattern", "**/*.SIF"]) == 0
    )

    out = capsys.readouterr().out
    assert "Skipped:" not in out
    manifest = _run_manifest(tmp_path / "runs")
    assert "pruned_dirs" not in manifest
    assert manifest["expected_files"] == 2


def test_the_catalog_and_the_page_carry_the_skip_the_receipt_recorded(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    """The campaign page was the last surface where a pruned run read complete.

    The receipt has said it since this packet's first half.  The drive catalog
    now carries it onward, and the page prints it on the drive's own card, so an
    operator reading the index months later is told a day was left out rather
    than left to notice the missing shots.
    """

    import json

    from echelle_spectra.catalog import build_drive_catalog, merge_catalogs
    from echelle_spectra.reading_room import build_reading_room

    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")
    _shots(_snapshot_marker(source / "20190207"), "193780")

    assert _process(_batch_argv(tmp_path, source)) == 0
    capsys.readouterr()

    receipt = sorted((tmp_path / "runs").rglob("run.toml"))[0].parent
    drive_catalog = build_drive_catalog(
        tmp_path / "cubes", volume_label="NIFS-A", drive_id="id-a", receipt_dir=receipt
    )
    payload = json.loads(drive_catalog.read_text(encoding="utf-8"))
    assert payload["run"]["pruned_dirs"] == ["20190207"]

    merged = merge_catalogs([drive_catalog], tmp_path / "all-years.json")
    page = build_reading_room(merged, tmp_path / "web").read_text(encoding="utf-8")
    assert "This run skipped 1 calibration folder(s): 20190207" in page


def test_a_run_that_pruned_nothing_writes_no_key_into_its_catalog(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    """Absence is the statement: no key means no skip, in the catalog too."""

    import json

    from echelle_spectra.catalog import build_drive_catalog

    source = tmp_path / "drive"
    _shots(source / "20190206", "193778")

    assert _process(_batch_argv(tmp_path, source)) == 0
    capsys.readouterr()

    receipt = sorted((tmp_path / "runs").rglob("run.toml"))[0].parent
    drive_catalog = build_drive_catalog(
        tmp_path / "cubes", volume_label="NIFS-A", drive_id="id-a", receipt_dir=receipt
    )
    assert "pruned_dirs" not in json.loads(drive_catalog.read_text(encoding="utf-8"))["run"]
