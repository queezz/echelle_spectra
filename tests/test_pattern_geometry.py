"""What `echelle drift geometry` says about where a drive's light sits."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from echelle_spectra import pattern_survey
from echelle_spectra.campaign_tools_cli import drift_main
from echelle_spectra.cli import main as echelle_main

# The same continuum frame the bench's band-centre tests are built from: four
# order bands, sloping, wider than the strip the extractor sums out of them.
_COLUMNS = 640
_ROWS = 260
_ORDER_ROWS = (40.0, 100.0, 160.0, 220.0)
_SLOPE = 0.01
_BAND_SIGMA = 5.0
_LEVEL = 20000.0
_FLOOR = 40.0


def _traces(shift_rows: float = 0.0) -> np.ndarray:
    columns = np.arange(_COLUMNS, dtype=float)
    return np.column_stack([row + _SLOPE * columns + shift_rows for row in _ORDER_ROWS])


def _image(shift_rows: float = 0.0, *, level: float = _LEVEL) -> np.ndarray:
    rows = np.arange(_ROWS, dtype=float)[:, None]
    traces = _traces(shift_rows)
    image = np.full((_ROWS, _COLUMNS), _FLOOR)
    for order_idx in range(traces.shape[1]):
        centers = traces[:, order_idx][None, :]
        image = image + level * np.exp(-0.5 * ((rows - centers) / _BAND_SIGMA) ** 2)
    return image


def _pattern_file(tmp_path, name="pattern.txt", shift_rows=0.0):
    path = tmp_path / name
    np.savetxt(path, np.rint(_traces(shift_rows)).astype(int), fmt="%d")
    return path


def _drive(tmp_path, day_names=("20250926", "20250927", "20250928")):
    """A synthetic drive: DATA/<day>/<shot>_Echelle.SIF, contents unread."""

    data = tmp_path / "LHD2025" / "DATA"
    for name in day_names:
        folder = data / name
        folder.mkdir(parents=True)
        for shot in range(3):
            (folder / f"{shot:06d}_Echelle.SIF").write_bytes(b"")
    return data.parent


def _reader(shifts):
    """A frame reader answering with a stack whose shift depends on the day."""

    def read(path):
        shift = shifts[path.parent.name] if isinstance(shifts, dict) else shifts
        return np.stack([_image(shift)])

    return read


def _survey(root, patterns, **kwargs):
    lines: list[str] = []
    payload = pattern_survey.survey_geometry(
        root, patterns=patterns, report=lines.append, **kwargs
    )
    return payload, lines


def test_light_on_its_own_pattern_reads_geometry_ok(tmp_path):
    """A campaign sitting on its pattern is inside the alarm and exits zero."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, every=1, reader=_reader(0.0))
    lines, code = pattern_survey.verdict_lines(payload)

    assert payload["verdict"] == "GEOMETRY_OK"
    assert code == 0
    assert lines[-1] == "GEOMETRY_OK"
    assert "inside the 1.5-row alarm" in lines[0]


def test_light_off_the_pattern_raises_the_geometry_alarm(tmp_path):
    """Bands five rows off the pattern refuse the run and exit two."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, every=1, reader=_reader(5.0))
    lines, code = pattern_survey.verdict_lines(payload)

    assert payload["verdict"] == "GEOMETRY_ALARM"
    assert code == 2
    assert lines[-1] == "GEOMETRY_ALARM"
    assert "PAST the 1.5-row alarm" in lines[0]
    assert payload["worst_offset_rows"] == pytest.approx(5.0, abs=0.3)


def test_a_root_with_no_day_folders_is_unmeasured_never_ok(tmp_path):
    """Nothing to read leaves with its own code, so silence is not permission."""

    root = tmp_path / "EMPTY"
    (root / "DATA").mkdir(parents=True)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, reader=_reader(0.0))
    lines, code = pattern_survey.verdict_lines(payload)

    assert code == 3
    assert lines[-1] == "GEOMETRY_ALARM"
    assert "no day folders under" in lines[0]


def test_a_day_folder_with_no_frames_says_no_sif_files(tmp_path):
    """Day folders that hold nothing are a different unmeasured answer."""

    root = tmp_path / "LHD2025"
    (root / "DATA" / "20250926").mkdir(parents=True)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, reader=_reader(0.0))
    lines, code = pattern_survey.verdict_lines(payload)

    assert code == 3
    assert "no SIF files found under" in lines[0]


def test_a_frame_too_dim_to_measure_is_unknown_never_ok(tmp_path):
    """A frame with no bands in it reports geometry UNKNOWN, not zero offset."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    def dark(path):
        return np.full((1, _ROWS, _COLUMNS), 12.0)

    payload, _ = _survey(root, patterns, every=1, reader=dark)
    lines, code = pattern_survey.verdict_lines(payload)

    assert code == 3
    assert payload["worst_offset_rows"] is None
    assert "geometry UNKNOWN" in lines[0]
    assert lines[-1] == "GEOMETRY_ALARM"


def test_an_unreadable_frame_is_one_row_and_the_survey_goes_on(tmp_path):
    """One frame the reader refuses does not cost every day after it."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])
    refused = {"20250927"}

    def sometimes(path):
        if path.parent.name in refused:
            raise OSError("SIF header is truncated")
        return np.stack([_image(0.0)])

    payload, lines = _survey(root, patterns, every=1, reader=sometimes)

    assert len(payload["rows"]) == 3
    bad = [row for row in payload["rows"] if row["day"] == "20250927"]
    assert bad[0]["error"] == "SIF header is truncated"
    assert any("unreadable: SIF header is truncated" in line for line in lines)
    assert pattern_survey.verdict_lines(payload)[1] == 0


def test_the_survey_records_one_row_per_sampled_day(tmp_path):
    """The JSON is the campaign's geometry diary, one row per day it read."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, every=1, reader=_reader(0.0))
    written = pattern_survey.write_survey(payload, tmp_path)
    reloaded = json.loads(written.read_text(encoding="utf-8"))

    assert written == tmp_path / "inventory" / "LHD2025-pattern-survey.json"
    assert reloaded["schema"] == "echelle-pattern-survey/v1"
    assert set(reloaded) >= {
        "generated_at",
        "drive_root",
        "every_nth_day",
        "patterns",
        "rows",
        "verdict",
        "verdict_pattern",
        "alarm_rows",
        "worst_offset_rows",
    }
    assert [row["day"] for row in reloaded["rows"]] == ["20250926", "20250927", "20250928"]
    assert set(reloaded["rows"][0]["offsets_rows"]) == {"pattern.txt"}


def test_a_survey_of_a_data_folder_is_still_filed_under_the_drive(tmp_path):
    """Pointing at DATA must not file every drive's survey under the name DATA."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root / "DATA", patterns, every=1, reader=_reader(0.0))
    written = pattern_survey.write_survey(payload, tmp_path)

    assert written.name == "LHD2025-pattern-survey.json"


def test_every_nth_day_reads_the_days_it_says(tmp_path):
    """--every 2 reads the first, third and fifth day and no others."""

    root = _drive(tmp_path, day_names=("01", "02", "03", "04", "05"))
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, every=2, reader=_reader(0.0))

    assert [row["day"] for row in payload["rows"]] == ["01", "03", "05"]
    assert payload["every_nth_day"] == 2


def test_every_zero_falls_back_to_the_three_frame_spread(tmp_path):
    """Zero asks the older question: first day, middle day, last day."""

    root = _drive(tmp_path, day_names=("01", "02", "03", "04", "05"))
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(root, patterns, every=0, reader=_reader(0.0))

    assert [row["day"] for row in payload["rows"]] == ["01", "03", "05"]


def test_a_single_day_root_is_read_as_that_one_day(tmp_path):
    """A folder of frames handed over on its own is a day, not an empty drive."""

    day = tmp_path / "20250926"
    day.mkdir()
    for shot in range(3):
        (day / f"{shot:06d}_Echelle.SIF").write_bytes(b"")
    patterns = pattern_survey.load_patterns([_pattern_file(tmp_path)])

    payload, _ = _survey(day, patterns, reader=_reader(0.0))

    assert [row["day"] for row in payload["rows"]] == ["20250926"]
    assert pattern_survey.verdict_lines(payload)[1] == 0


def test_the_verdict_follows_the_last_pattern_given(tmp_path):
    """Every pattern is measured; only the last one decides the run."""

    root = _drive(tmp_path)
    patterns = pattern_survey.load_patterns(
        [_pattern_file(tmp_path, "old.txt", shift_rows=-5.0), _pattern_file(tmp_path, "new.txt")]
    )

    payload, _ = _survey(root, patterns, every=1, reader=_reader(0.0))

    assert payload["verdict_pattern"] == "new.txt"
    assert payload["verdict"] == "GEOMETRY_OK"
    row = payload["rows"][0]
    assert row["offsets_rows"]["old.txt"] == pytest.approx(5.0, abs=0.3)
    assert row["offsets_rows"]["new.txt"] == pytest.approx(0.0, abs=0.3)


def test_the_brightest_frame_of_the_stack_is_the_one_measured(tmp_path):
    """Mid-day shots hold several frames; the reading takes the brightest."""

    stack = np.stack([_image(0.0, level=10.0), _image(5.0), _image(0.0, level=10.0)])

    assert pattern_survey.brightest_frame(stack) is not None
    np.testing.assert_allclose(pattern_survey.brightest_frame(stack), stack[1])


def test_the_packaged_eras_are_the_default_patterns():
    """With no --pattern the 2024 and 2025 patterns are measured, newest last."""

    patterns = pattern_survey.load_patterns(None)

    assert [choice.key for choice in patterns] == ["2024", "2025"]
    assert all(choice.rows.ndim == 2 for choice in patterns)


def test_the_verb_prints_the_table_and_ends_on_the_verdict(tmp_path, capsys, monkeypatch):
    """`echelle drift geometry` writes the survey and leaves with the exit code."""

    root = _drive(tmp_path)
    pattern = _pattern_file(tmp_path)
    monkeypatch.setattr(pattern_survey, "frame_reader", lambda rows: _reader(5.0))
    output = tmp_path / "campaign"
    output.mkdir()

    code = drift_main(
        [
            "geometry",
            str(root),
            "--pattern",
            str(pattern),
            "--every",
            "1",
            "--output",
            str(output),
        ]
    )

    printed = capsys.readouterr().out.strip().splitlines()
    assert code == 2
    assert printed[0].startswith("day        shot")
    assert printed[-1] == "GEOMETRY_ALARM"
    assert (output / "inventory" / "LHD2025-pattern-survey.json").is_file()


def test_a_missing_pattern_file_is_refused_in_one_line(tmp_path, capsys):
    """A pattern that is not there earns a sentence, not a traceback."""

    root = _drive(tmp_path)

    code = drift_main(["geometry", str(root), "--pattern", str(tmp_path / "gone.txt")])

    assert code == 1
    assert "order pattern not found" in capsys.readouterr().err


def test_drift_geometry_help_exits_zero():
    """`echelle drift geometry --help` is reachable through the umbrella command."""

    with pytest.raises(SystemExit) as exit_info:
        echelle_main(["drift", "geometry", "--help"])

    assert exit_info.value.code == 0


def test_two_patterns_sharing_a_filename_get_their_own_columns(tmp_path: Path) -> None:
    """Snapshot patterns are all named pattern.txt, so keying columns by bare
    filename collapsed two eras into one silently-identical pair of columns —
    the first field survey against the 2019 and 2024 snapshots printed one
    measurement twice. Colliding names are disambiguated by their folder.
    """

    from echelle_spectra.pattern_survey import load_patterns

    a = tmp_path / "20190529_cmos" / "pattern.txt"
    b = tmp_path / "20240305_cmos" / "pattern.txt"
    for path, base in ((a, 40), (b, 49)):
        path.parent.mkdir()
        rows = np.tile(np.arange(base, base + 3 * 10, 10), (64, 1))
        np.savetxt(path, rows, fmt="%d")
    choices = load_patterns([a, b])
    keys = [choice.key for choice in choices]
    assert keys == ["20190529_cmos/pattern.txt", "20240305_cmos/pattern.txt"]
    assert not np.array_equal(choices[0].rows, choices[1].rows)
