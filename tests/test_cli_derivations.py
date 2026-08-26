"""What the operator no longer has to type, the commands work out themselves.

The campaign page is being narrowed to two questions -- which folder, and which
calibration -- so everything else has to be derived where the answer actually
lives: on disk, at run time.  These tests pin the four derivations that move off
the page and into the CLIs.

The sharpest of them is discovery.  A campaign drive is a tree of date-named day
folders with the drive's own lamp frames sitting beside them, so the flat
``*.SIF`` search found nothing at all, and the hand-written ``**/*.SIF`` an
operator reached for instead quietly exported the calibration lamps as science
shots.

The last section pins ``--practice``: the page learned with empty hands.  The
tool that composes CLI commands must not need a CLI command composed before it
can help, so practice invents its own campaign and builds the same page a real
one would, into a fresh system-temp folder nothing else reads.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from echelle_spectra import campaign_tools_cli
from echelle_spectra.campaign_tools_cli import (
    _derived_every,
    _derived_evidence_path,
    _evidence_home,
    _next_free,
    drift_main,
    web_main,
)
from echelle_spectra.spectrocube_cli import (
    _authorize_run,
    _auto_sample_count,
    _build_parser,
    _discover_input_files,
    main,
)

MERGED_CATALOG = {
    "schema": "echelle-merged-catalog/v1",
    "generated_at": "2026-08-18T00:00:00.000+00:00",
    "sources": [],
}


@pytest.fixture()
def mock_spectrocube():
    """Inject a fake spectrocube module, as the batch tests next door do."""

    fake = MagicMock()
    with patch.dict("sys.modules", {"spectrocube": fake}):
        yield fake


def _drive(tmp_path: Path, *, suffix: str = ".SIF") -> Path:
    """A campaign drive: two day folders, lamp frames, and a saved snapshot."""

    root = tmp_path / "drive"
    for day, shots in (("20190206", ("193778", "193779")), ("20190207", ("193780",))):
        folder = root / day
        folder.mkdir(parents=True)
        for shot in shots:
            (folder / f"{shot}_Echelle{suffix}").touch()

    lamps = root / "calibrations"
    lamps.mkdir()
    for lamp in ("Ne", "Kr"):
        (lamps / f"{lamp}_lamp{suffix}").touch()

    snapshot = root / "20250101_cmos"
    snapshot.mkdir()
    (snapshot / "snapshot.toml").write_text('id = "20250101_cmos"\n', encoding="utf-8")
    (snapshot / f"Ne_source{suffix}").touch()
    nested = snapshot / "sources"
    nested.mkdir()
    (nested / f"Kr_source{suffix}").touch()

    return root


def _names(paths: list[Path], root: Path) -> list[str]:
    return [path.relative_to(root).as_posix() for path in paths]


# ---------------------------------------------------------------------------
# 1. Discovery walks the day folders and refuses to mistake lamps for shots
# ---------------------------------------------------------------------------


def test_discovery_finds_every_day_folder_shot_and_no_calibration_frame(
    tmp_path: Path,
) -> None:
    root = _drive(tmp_path)
    (root / "196000_Echelle.SIF").touch()

    found = _discover_input_files(root, "*.SIF")

    assert _names(found, root) == [
        "196000_Echelle.SIF",
        "20190206/193778_Echelle.SIF",
        "20190206/193779_Echelle.SIF",
        "20190207/193780_Echelle.SIF",
    ]


def test_a_flat_drive_still_works_and_each_file_appears_once(tmp_path: Path) -> None:
    for shot in ("a", "b"):
        (tmp_path / f"{shot}.SIF").touch()

    found = _discover_input_files(tmp_path, "*.SIF")

    assert _names(found, tmp_path) == ["a.SIF", "b.SIF"]
    assert len(found) == len(set(found))


def test_the_lowercase_fallback_survives_the_walk(tmp_path: Path) -> None:
    root = _drive(tmp_path, suffix=".sif")

    found = _discover_input_files(root, "*.SIF")

    assert _names(found, root) == [
        "20190206/193778_Echelle.sif",
        "20190206/193779_Echelle.sif",
        "20190207/193780_Echelle.sif",
    ]
    # On a case-insensitive filesystem `*.SIF` already matched these names; the
    # retry must not hand the same file back a second time.
    assert len(found) == len(set(found))


def test_a_path_shaped_pattern_is_the_operators_own_and_is_used_as_typed(
    tmp_path: Path,
) -> None:
    root = _drive(tmp_path)

    assert _names(_discover_input_files(root, "20190206/*.SIF"), root) == [
        "20190206/193778_Echelle.SIF",
        "20190206/193779_Echelle.SIF",
    ]
    # `**` is the escape hatch, and it escapes the exclusion too.
    assert "calibrations/Ne_lamp.SIF" in _names(_discover_input_files(root, "**/*.SIF"), root)


def test_a_drive_of_only_lamp_frames_is_refused_in_words_that_explain_the_search(
    tmp_path: Path, capsys, mock_spectrocube
) -> None:
    root = tmp_path / "drive"
    lamps = root / "calibrations"
    lamps.mkdir(parents=True)
    (lamps / "Ne_lamp.SIF").touch()

    with pytest.raises(SystemExit) as exit_code:
        main([str(root), "--dry-run"])

    assert exit_code.value.code == 1
    refusal = capsys.readouterr().err
    assert "recursive" in refusal
    assert "calibration folders" in refusal
    assert "snapshot.toml" in refusal
    assert str(root) in refusal


def test_the_batch_run_exports_the_day_folders_it_walked(
    tmp_path: Path, capsys, mock_spectrocube
) -> None:
    root = _drive(tmp_path)

    with pytest.raises(SystemExit) as exit_code:
        main([str(root), "--dry-run", "--verbose", "-o", str(tmp_path / "cubes")])

    assert exit_code.value.code == 0
    shown = capsys.readouterr().out
    assert "Files:       3" in shown
    assert "[1/3] 193778_Echelle.SIF" in shown
    assert "Ne_lamp.SIF" not in shown


# ---------------------------------------------------------------------------
# 2. --sample auto
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("files", "expected"),
    [(0, 0), (1, 1), (4, 4), (5, 5), (30, 5), (125, 5), (200, 8), (750, 30), (5000, 30)],
)
def test_the_auto_sample_size_follows_the_stated_rule(files: int, expected: int) -> None:
    assert _auto_sample_count(files) == expected


def test_sample_takes_auto_or_a_plain_integer_and_refuses_anything_else() -> None:
    assert _build_parser().parse_args(["/data", "--sample", "auto"]).sample == "auto"
    assert _build_parser().parse_args(["/data", "--sample", "12"]).sample == 12
    with pytest.raises(SystemExit):
        _build_parser().parse_args(["/data", "--sample", "five"])


def test_an_auto_sample_says_how_many_of_how_many_it_took(tmp_path: Path, capsys) -> None:
    files = [tmp_path / f"{index:05d}.SIF" for index in range(200)]
    args = SimpleNamespace(sample="auto", drift_verdict=None)

    authorization = _authorize_run(
        args, {}, object(), files, input_path=tmp_path, target_label=None
    )

    assert authorization is not None
    assert len(authorization.files) == 8
    assert authorization.sample is True
    assert "sample auto: 8 of 200 files" in capsys.readouterr().out


def test_a_plain_integer_sample_is_untouched_by_the_new_word(tmp_path: Path) -> None:
    files = [tmp_path / f"{index:05d}.SIF" for index in range(200)]
    args = SimpleNamespace(sample=3, drift_verdict=None)

    authorization = _authorize_run(
        args, {}, object(), files, input_path=tmp_path, target_label=None
    )

    assert authorization is not None
    assert len(authorization.files) == 3


def test_auto_still_needs_a_registry_like_every_other_sample(tmp_path: Path, capsys) -> None:
    with pytest.raises(SystemExit) as exit_code:
        main([str(tmp_path), "--sample", "auto"])

    assert exit_code.value.code == 2
    assert "requires --registry" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# 3. echelle drift audit derives its evidence name and its interval
# ---------------------------------------------------------------------------


def test_the_evidence_is_named_beside_the_cubes_it_describes(tmp_path: Path) -> None:
    cubes = tmp_path / "cubes"
    cubes.mkdir()

    assert _evidence_home([cubes]) == cubes
    assert _derived_evidence_path([cubes]) == cubes / "drift-evidence-001.json"

    (cubes / "drift-evidence-001.json").write_text("{}", encoding="utf-8")
    assert _derived_evidence_path([cubes]) == cubes / "drift-evidence-002.json"


def test_a_taken_evidence_name_is_answered_in_the_same_numbering(tmp_path: Path) -> None:
    """The refusal offers a name from the one series the folder already uses."""

    taken = tmp_path / "drift-evidence.json"
    taken.write_text("{}", encoding="utf-8")
    assert _next_free(taken) == tmp_path / "drift-evidence-001.json"

    (tmp_path / "drift-evidence-001.json").write_text("{}", encoding="utf-8")
    assert _next_free(taken) == tmp_path / "drift-evidence-002.json"
    # An already-numbered name is never offered back to itself.
    assert _next_free(tmp_path / "drift-evidence-001.json") == (
        tmp_path / "drift-evidence-002.json"
    )


def test_cubes_named_one_by_one_share_their_folder(tmp_path: Path) -> None:
    cubes = tmp_path / "cubes"
    (cubes / "day-a").mkdir(parents=True)
    (cubes / "day-b").mkdir(parents=True)
    first = cubes / "day-a" / "1.nc"
    second = cubes / "day-b" / "2.nc"
    first.touch()
    second.touch()

    assert _evidence_home([first, second]) == cubes


def test_the_audit_announces_the_name_it_derived(tmp_path: Path, capsys) -> None:
    cubes = tmp_path / "cubes"
    cubes.mkdir()

    # The empty folder makes the audit itself refuse, which is exactly what an
    # unnamed evidence file must not depend on: the name was derived first.
    assert drift_main(["audit", str(cubes)]) == 1

    shown = capsys.readouterr()
    assert str(cubes / "drift-evidence-001.json") in shown.out
    assert "Traceback" not in shown.out + shown.err


@pytest.mark.parametrize(("cubes", "expected"), [(0, 1), (10, 1), (40, 2), (400, 20)])
def test_the_interval_is_derived_to_measure_about_twenty_cubes(
    tmp_path: Path, cubes: int, expected: int
) -> None:
    folder = tmp_path / "cubes"
    folder.mkdir()
    for index in range(cubes):
        (folder / f"{index:04d}.nc").touch()

    assert _derived_every([folder]) == (expected, cubes)


def _dated_cube(path: Path, day: str) -> Path:
    """One cube carrying nothing but the acquisition date the filter reads."""

    import numpy as np
    import xarray as xr

    path.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(
        {"intensity": ("wavelength", np.zeros(1))},
        coords={"wavelength": np.zeros(1)},
        attrs={"t_start": f"{day}T09:00:00+09:00"},
    ).to_netcdf(path)
    return path


def _dated_drive(folder: Path, inside: int, outside: int) -> Path:
    """A drive whose cubes mostly sit outside the window an operator asks for."""

    for index in range(outside):
        _dated_cube(folder / f"old-{index:04d}.nc", "2019-06-01")
    for index in range(inside):
        _dated_cube(folder / f"new-{index:04d}.nc", "2026-08-14")
    return folder


WINDOW = {"date_from": "2026-08-01", "date_to": "2026-08-31"}


def test_the_interval_is_derived_from_the_cubes_the_dates_actually_keep(
    tmp_path: Path,
) -> None:
    """The gate-integrity case: dividing by the drive would measure one cube.

    ``audit_cubes`` filters by date and *then* takes every Nth survivor.  A
    derivation counted over the whole drive -- 100 here, 40 inside the window --
    would hand ``select_sample_paths`` an interval of 5 for a list of 40 and
    measure 8 of them, under a line promising twenty.  At campaign scale the
    same arithmetic measures exactly one.
    """

    from echelle_spectra.drift import _filter_by_date, resolve_cube_paths, select_sample_paths

    folder = _dated_drive(tmp_path / "cubes", inside=40, outside=60)

    assert _derived_every([folder], **WINDOW) == (2, 40)

    # And the interval, applied the way the audit applies it, measures twenty.
    kept = _filter_by_date(resolve_cube_paths([folder]), **WINDOW)
    assert len(select_sample_paths(kept, every=2)) == 20


def test_an_appledouble_sibling_is_never_resolved_as_a_cube(tmp_path: Path) -> None:
    """A Mac writing cubes to an exFAT trip drive leaves a ``._cube.nc``
    metadata sibling beside every real cube.  The first real field audit
    crashed opening one as a dataset; a directory expansion skips them.
    """

    from echelle_spectra.drift import resolve_cube_paths

    folder = tmp_path / "cubes"
    folder.mkdir()
    real = folder / "196201_Echelle_spectrocube.nc"
    real.write_bytes(b"cube bytes")
    (folder / "._196201_Echelle_spectrocube.nc").write_bytes(b"\x00\x05\x16\x07")

    assert resolve_cube_paths([folder]) == [real]


def test_the_unfiltered_derivation_is_unchanged(tmp_path: Path) -> None:
    folder = _dated_drive(tmp_path / "cubes", inside=40, outside=60)

    assert _derived_every([folder]) == (5, 100)


def test_a_window_holding_fewer_than_twenty_cubes_measures_all_of_them(
    tmp_path: Path,
) -> None:
    from echelle_spectra.drift import _filter_by_date, resolve_cube_paths, select_sample_paths

    folder = _dated_drive(tmp_path / "cubes", inside=12, outside=60)

    assert _derived_every([folder], **WINDOW) == (1, 12)

    kept = _filter_by_date(resolve_cube_paths([folder]), **WINDOW)
    assert len(select_sample_paths(kept, every=1)) == 12


def test_the_derived_interval_says_which_count_it_divided(tmp_path: Path, capsys) -> None:
    folder = _dated_drive(tmp_path / "cubes", inside=40, outside=60)

    with patch("echelle_spectra.drift.audit_cubes") as audited:
        audited.return_value = {"verdict": "aligned"}
        with patch("echelle_spectra.drift.write_drift_evidence") as written:
            written.return_value = tmp_path / "evidence.json"
            assert (
                drift_main(
                    [
                        "audit",
                        str(folder),
                        "--from",
                        WINDOW["date_from"],
                        "--to",
                        WINDOW["date_to"],
                        "-o",
                        str(tmp_path / "evidence.json"),
                    ]
                )
                == 0
            )

    assert audited.call_args.kwargs["every"] == 2
    shown = capsys.readouterr().out
    assert "every: 2 (derived from 40 cubes in --from/--to;" in shown
    # The drive's 100 cubes are not what this number came from, and saying so
    # is the only check an operator can run on it.
    assert "100" not in shown


def test_an_explicit_interval_wins(tmp_path: Path) -> None:
    folder = tmp_path / "cubes"
    folder.mkdir()
    for index in range(400):
        (folder / f"{index:04d}.nc").touch()

    with patch("echelle_spectra.drift.audit_cubes") as audited:
        audited.return_value = {"verdict": "aligned"}
        with patch("echelle_spectra.drift.write_drift_evidence") as written:
            written.return_value = tmp_path / "evidence.json"
            assert drift_main(["audit", str(folder), "--every", "7"]) == 0

    assert audited.call_args.kwargs["every"] == 7


# ---------------------------------------------------------------------------
# 4. echelle web: the campaign home, the absolute page, and --open
# ---------------------------------------------------------------------------


def _catalog(tmp_path: Path, name: str = "all-years.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(MERGED_CATALOG), encoding="utf-8")
    return path


def _home(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "campaign.toml"
    path.write_text(body, encoding="utf-8")
    return path


def test_a_campaign_home_answers_both_questions_the_page_used_to_ask(
    tmp_path: Path, capsys
) -> None:
    _catalog(tmp_path)
    _home(
        tmp_path,
        '# the campaign, written down once\ncatalog = "all-years.json"\noutput = "page"\n',
    )
    (tmp_path / "page").mkdir()

    assert web_main(["--home", str(tmp_path)]) == 0

    shown = capsys.readouterr().out
    assert str(tmp_path / "campaign.toml") in shown
    assert str(tmp_path / "page" / "index.html") in shown
    assert (tmp_path / "page" / "index.html").is_file()


def test_paths_inside_the_home_are_read_relative_to_the_home_itself(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    _catalog(campaign)
    (campaign / "page").mkdir()
    _home(campaign, 'catalog = "all-years.json"\noutput = "page"\n')
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    assert web_main(["--home", str(campaign / "campaign.toml")]) == 0
    assert str(campaign / "page" / "index.html") in capsys.readouterr().out


def test_an_explicit_flag_always_beats_the_home(tmp_path: Path, capsys) -> None:
    _catalog(tmp_path)
    _home(tmp_path, 'catalog = "all-years.json"\noutput = "page"\n')
    (tmp_path / "page").mkdir()
    (tmp_path / "typed").mkdir()

    assert web_main(["--home", str(tmp_path), "--output", str(tmp_path / "typed")]) == 0

    assert (tmp_path / "typed" / "index.html").is_file()
    assert not (tmp_path / "page").exists() or not (tmp_path / "page" / "index.html").exists()


def test_a_campaign_home_in_this_folder_is_used_and_said_out_loud(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _catalog(tmp_path)
    _home(tmp_path, 'catalog = "all-years.json"\noutput = "page"\n')
    (tmp_path / "page").mkdir()
    monkeypatch.chdir(tmp_path)

    assert web_main([]) == 0

    shown = capsys.readouterr().out
    assert "campaign home:" in shown
    assert str(tmp_path / "campaign.toml") in shown


def test_without_a_home_the_required_arguments_still_stand(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    monkeypatch.chdir(tmp_path)

    with pytest.raises(SystemExit) as exit_code:
        web_main([])

    assert exit_code.value.code == 2
    err = capsys.readouterr().err
    assert "the following arguments are required: --catalog, --output" in err
    # The cold start is where an operator first meets the command, so the
    # refusal teaches the campaign home instead of stopping at argparse.
    assert "no campaign.toml was found in this folder" in err
    assert "echelle web --open" in err
    assert 'catalog = "all-years.json"' in err
    assert "echelle status" in err


def test_open_hands_the_written_file_to_the_browser_and_nothing_else(
    tmp_path: Path, monkeypatch
) -> None:
    _catalog(tmp_path)
    (tmp_path / "page").mkdir()
    opened: list[str] = []
    monkeypatch.setattr(
        campaign_tools_cli.webbrowser,
        "open",
        lambda uri: opened.append(uri) or True,
    )

    assert (
        web_main(
            [
                "--catalog",
                str(tmp_path / "all-years.json"),
                "--output",
                str(tmp_path / "page"),
                "--open",
            ]
        )
        == 0
    )

    assert opened == [(tmp_path / "page" / "index.html").resolve().as_uri()]


def test_a_build_without_open_touches_no_browser(tmp_path: Path, monkeypatch) -> None:
    _catalog(tmp_path)
    (tmp_path / "page").mkdir()
    opened: list[str] = []
    monkeypatch.setattr(
        campaign_tools_cli.webbrowser, "open", lambda uri: opened.append(uri) or True
    )

    assert (
        web_main(
            [
                "--catalog",
                str(tmp_path / "all-years.json"),
                "--output",
                str(tmp_path / "page"),
            ]
        )
        == 0
    )

    assert opened == []


# ---------------------------------------------------------------------------
# 5. --practice: the page learned with empty hands
# ---------------------------------------------------------------------------


def _shown(capsys) -> str:
    captured = capsys.readouterr()
    return captured.out + captured.err


def test_practice_builds_an_invented_campaign(capsys) -> None:
    assert web_main(["--practice"]) == 0

    shown = _shown(capsys)
    assert "Traceback" not in shown
    lines = [line for line in shown.splitlines() if line.strip()]
    assert lines[0] == "practice campaign: every fact on this page is invented"

    index_path = Path(lines[1])
    assert index_path.is_absolute()
    assert index_path.is_file()
    assert index_path.name == "index.html"

    text = index_path.read_text(encoding="utf-8")
    assert "PRACTICE-A" in text
    assert "PRACTICE-B" in text
    assert "PRACTICE-C" in text
    assert "PRACTICE-D" in text
    # A real verdict card, and the composer derives the NEXT free evidence
    # name from the invented drift-evidence-001.json it was handed.
    assert "shifted" in text
    assert "drift-evidence-002.json" in text


def test_practice_refuses_combination_with_catalog(capsys) -> None:
    assert web_main(["--practice", "--catalog", "all-years.json"]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "--catalog" in shown
    assert "--practice" in shown
    assert len([line for line in shown.splitlines() if line.strip()]) == 1


def test_practice_refuses_combination_with_home(tmp_path: Path, capsys) -> None:
    assert web_main(["--practice", "--home", str(tmp_path)]) == 1

    shown = _shown(capsys)
    assert "Traceback" not in shown
    assert "--home" in shown
    assert "--practice" in shown
    assert len([line for line in shown.splitlines() if line.strip()]) == 1


def test_practice_writes_nothing_outside_its_temp_folder(tmp_path, monkeypatch, capsys) -> None:
    practice_root = tmp_path / "practice-root"

    def fake_mkdtemp(*args, **kwargs) -> str:
        practice_root.mkdir(parents=True, exist_ok=True)
        return str(practice_root)

    monkeypatch.setattr(campaign_tools_cli.tempfile, "mkdtemp", fake_mkdtemp)

    assert web_main(["--practice"]) == 0

    outside = [
        path
        for path in tmp_path.rglob("*")
        if path != practice_root and practice_root not in path.parents
    ]
    assert outside == [], f"practice wrote outside its temp folder: {outside}"

    shown = _shown(capsys)
    lines = [line for line in shown.splitlines() if line.strip()]
    index_path = Path(lines[1])
    assert practice_root in index_path.parents


def test_the_cold_start_refusal_offers_practice(tmp_path: Path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)

    with pytest.raises(SystemExit):
        web_main([])

    assert "echelle web --practice --open" in capsys.readouterr().err


def test_web_open_with_nothing_serves_the_page_instead_of_lecturing(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """The typed command is somebody asking to USE the page: start it."""

    monkeypatch.chdir(tmp_path)
    served: list[list[str]] = []
    from echelle_spectra import campaign_server

    monkeypatch.setattr(
        campaign_server, "serve_main", lambda argv=None, **_: served.append(argv) or 0
    )

    assert web_main(["--open"]) == 0

    assert served == [["--open"]]
    assert "serving the campaign page instead" in capsys.readouterr().out


def test_web_open_with_a_home_in_reach_still_builds_statically(
    tmp_path: Path, monkeypatch
) -> None:
    _catalog(tmp_path)
    _home(tmp_path, 'catalog = "all-years.json"\noutput = "page"\n')
    (tmp_path / "page").mkdir()
    monkeypatch.chdir(tmp_path)
    opened: list[str] = []
    monkeypatch.setattr(
        campaign_tools_cli.webbrowser, "open", lambda uri: opened.append(uri) or True
    )

    assert web_main(["--open"]) == 0

    assert opened == [(tmp_path / "page" / "index.html").resolve().as_uri()]


# ---------------------------------------------------------------------------
# recal-cube tells the operator when it dropped factor columns
#
# The manifest has always recorded the drop; a cross-era run that quietly loses
# forty columns to the new sphere's narrower coverage needs to say so at the
# terminal too. The two recorded counts are not additive -- nonpositive is the
# TOTAL dropped, uncovered is the subset outside coverage -- so the line must
# name the total and call out the subset, never sum them.


def _write_recal_manifest(path: Path, **counts: int) -> Path:
    record = {"schema": "echelle-recalibration/v1", "changes": ["absolute-factor"]}
    record.update({key: value for key, value in counts.items() if value})
    path.write_text(json.dumps(record), encoding="utf-8")
    return path


def test_recal_reports_the_total_and_the_uncovered_subset(tmp_path: Path, capsys) -> None:
    manifest = _write_recal_manifest(
        tmp_path / "cube.nc.recalibration.json",
        dropped_nonpositive_factor_columns=40,
        dropped_uncovered_factor_columns=37,
    )
    campaign_tools_cli._report_dropped_factor_columns(manifest)
    out = capsys.readouterr().out
    assert "dropped 40" in out
    assert "37 of them outside" in out
    # The counts are a total and its subset; the line never presents their sum.
    assert "77" not in out


def test_recal_names_a_pure_noise_drop_without_a_coverage_clause(
    tmp_path: Path, capsys
) -> None:
    manifest = _write_recal_manifest(
        tmp_path / "cube.nc.recalibration.json", dropped_nonpositive_factor_columns=3
    )
    campaign_tools_cli._report_dropped_factor_columns(manifest)
    out = capsys.readouterr().out
    assert "dropped 3" in out
    assert "outside" not in out


def test_recal_says_nothing_when_no_column_was_dropped(tmp_path: Path, capsys) -> None:
    manifest = _write_recal_manifest(tmp_path / "cube.nc.recalibration.json")
    campaign_tools_cli._report_dropped_factor_columns(manifest)
    assert capsys.readouterr().out == ""
