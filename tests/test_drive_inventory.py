"""What `echelle inventory` says about a drive it has never seen before."""

from __future__ import annotations

import json

import pytest

from echelle_spectra import drive_inventory
from echelle_spectra.campaign_tools_cli import inventory_main
from echelle_spectra.cli import main as echelle_main

_LOGBOOK = "24thcycle_Echelle_setup.xlsx"


def _drive(tmp_path, days=(("20250926", (1, 2, 3)), ("20250927", (10, 11)))):
    """A synthetic campaign drive: day folders, frames, one calibration folder."""

    root = tmp_path / "LHD2025"
    data = root / "DATA"
    for name, shots in days:
        folder = data / name
        folder.mkdir(parents=True)
        for shot in shots:
            (folder / f"{shot:06d}_Echelle.SIF").write_bytes(b"x" * 1024)
        # A Mac on exFAT leaves one of these beside every real frame.
        (folder / f"._{shots[0]:06d}_Echelle.SIF").write_bytes(b"")
    calibration = data / "calibrations" / "20250926_sphere"
    calibration.mkdir(parents=True)
    (calibration / "snapshot.toml").write_text("", encoding="utf-8")
    return root


def test_a_drive_with_no_logbook_is_still_fully_inventoried(tmp_path):
    """Days, shot ranges and frame counts are read off the filenames alone."""

    record = drive_inventory.inventory_drive(_drive(tmp_path))

    assert record["schema"] == "echelle-drive-inventory/v1"
    assert record["volume_name"] == "LHD2025"
    assert record["logbooks"] == []
    assert [day["day"] for day in record["days"]] == ["20250926", "20250927"]
    assert record["days"][0]["sif_count"] == 3
    assert record["days"][0]["shot_range"] == [1, 3]
    assert record["days"][1]["shot_range"] == [10, 11]
    assert record["days"][0]["typical_file_bytes"] == 1024


def test_the_apple_double_siblings_are_not_counted_as_frames(tmp_path):
    """A `._name.SIF` metadata stub is not a frame and never inflates the count."""

    record = drive_inventory.inventory_drive(_drive(tmp_path))

    assert sum(day["sif_count"] for day in record["days"]) == 5


def test_a_calibration_folder_is_named_and_never_counted_as_a_day(tmp_path):
    """Calibration folders belong on their own line, not in the day census."""

    record = drive_inventory.inventory_drive(_drive(tmp_path))

    assert "calibrations" in record["calibration_folders"]
    assert "calibrations/20250926_sphere" in record["calibration_folders"]
    assert all("calib" not in day["day"] for day in record["days"])


def test_the_human_summary_names_the_span_and_the_cube_estimate(tmp_path):
    """The lines an operator reads say the span, the room and the calibrations."""

    record = drive_inventory.inventory_drive(_drive(tmp_path))

    text = "\n".join(drive_inventory.format_inventory(record))
    assert "=== LHD2025 ===" in text
    assert "day folders: 2, SIF files: 5" in text
    assert "span       : 20250926 (shot 1) .. 20250927 (shot 11)" in text
    assert "logbook    : (none on this drive)" in text


def test_a_drive_without_a_data_folder_reads_its_day_named_folders_and_says_so(tmp_path):
    """A DATA-less root counts date-named folders only, and names where it looked."""

    root = tmp_path / "LOOSE"
    for name in ("20250926", "cubes", "campaign-page"):
        (root / name).mkdir(parents=True)
    (root / "20250926" / "000042_Echelle.SIF").write_bytes(b"x")

    record = drive_inventory.inventory_drive(root)

    assert [day["day"] for day in record["days"]] == ["20250926"]
    assert record["data_root"] == str(root)
    assert "no DATA folder here" in "\n".join(drive_inventory.format_inventory(record))


def test_an_unreadable_logbook_is_a_finding_not_a_crash(tmp_path):
    """A workbook that will not parse is reported; the rest of the drive is not lost."""

    pytest.importorskip("openpyxl")
    root = _drive(tmp_path)
    (root / _LOGBOOK).write_text("this is not a workbook", encoding="utf-8")

    record = drive_inventory.inventory_drive(root)

    assert record["logbooks"][0]["file"] == _LOGBOOK
    assert record["logbooks"][0]["error"]
    assert len(record["days"]) == 2
    assert "UNREADABLE" in "\n".join(drive_inventory.format_inventory(record))


def test_without_openpyxl_a_logbook_is_present_but_unread(tmp_path, monkeypatch):
    """No reader installed is said plainly, never mistaken for no logbook."""

    root = _drive(tmp_path)
    (root / _LOGBOOK).write_text("anything", encoding="utf-8")
    monkeypatch.setattr(drive_inventory, "_openpyxl", lambda: None)

    record = drive_inventory.inventory_drive(root)

    assert record["logbooks"] == [
        {"file": _LOGBOOK, "unread": drive_inventory._NO_READER}
    ]
    assert "PRESENT BUT UNREAD" in "\n".join(drive_inventory.format_inventory(record))


def test_a_readable_logbook_reports_its_days_and_odd_trigger_delays(tmp_path):
    """A day taken at another trigger delay is called out by name."""

    openpyxl = pytest.importorskip("openpyxl")
    root = _drive(tmp_path)
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append(["date", "trig. delay(sec)", "exposure(sec)", "frames"])
    sheet.append([20250926, 2.5, 0.5, 40])
    sheet.append([20250927, 3.5, 0.5, 40])
    sheet.append([None, None, None, None])
    book.save(root / _LOGBOOK)

    record = drive_inventory.inventory_drive(root)

    logbook = record["logbooks"][0]
    assert logbook["rows"] == 2
    assert logbook["date_range"] == ["20250926", "20250927"]
    assert logbook["trig_delays_seen"] == [2.5, 3.5]
    assert logbook["days_not_2p5s"] == [
        {"date": "20250927", "trig_delay_s": 3.5, "exposure_s": 0.5, "frames": 40}
    ]
    assert "ODD DAY 20250927" in "\n".join(drive_inventory.format_inventory(record))


def test_the_record_is_written_under_the_output_folder_as_json(tmp_path):
    """The verb's record lands in inventory/<volume>.json and reloads as JSON."""

    root = _drive(tmp_path)
    output = tmp_path / "campaign"
    output.mkdir()

    written = drive_inventory.write_inventory(drive_inventory.inventory_drive(root), output)

    assert written == output / "inventory" / "LHD2025.json"
    reloaded = json.loads(written.read_text(encoding="utf-8"))
    assert reloaded["schema"] == "echelle-drive-inventory/v1"
    assert set(reloaded) >= {
        "generated_at",
        "drive_root",
        "volume_name",
        "size_bytes",
        "free_bytes",
        "logbooks",
        "days",
        "calibration_folders",
    }


def test_the_verb_prints_the_summary_and_leaves_with_zero(tmp_path, capsys):
    """`echelle inventory DRIVE --output DIR` reports and succeeds."""

    root = _drive(tmp_path)
    output = tmp_path / "campaign"
    output.mkdir()

    code = inventory_main([str(root), "--output", str(output)])

    assert code == 0
    assert "=== LHD2025 ===" in capsys.readouterr().out
    assert (output / "inventory" / "LHD2025.json").is_file()


def test_a_missing_drive_is_refused_in_one_line(tmp_path, capsys):
    """A typed path that is not there earns a sentence, not a traceback."""

    code = inventory_main([str(tmp_path / "not-mounted")])

    assert code == 1
    message = capsys.readouterr().err
    assert "drive root not found" in message
    assert "the DRIVE_ROOT argument" in message


def test_inventory_help_exits_zero():
    """`echelle inventory --help` is reachable through the umbrella command."""

    with pytest.raises(SystemExit) as exit_info:
        echelle_main(["inventory", "--help"])

    assert exit_info.value.code == 0
