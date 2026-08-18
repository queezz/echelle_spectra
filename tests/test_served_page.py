"""The served half: the folder picker, and the two states a cold campaign has.

The static build is the contract these tests exist to protect: ``served=False``
writes exactly the bytes it wrote before any of this existed, and every pin in
``test_reading_room_page.py`` still reads that same file.  Everything the server
makes possible — a real picker fed by ``/api/browse``, a setup screen that
writes the campaign home through ``/api/home`` — is appended and asserted here,
never woven into the static half.

The live half of the law (the Perimeter Walk, rail offset-identity sampling in a
real DOM, the dialog's own focus and scroll behaviour) is the owner's
validation work and is deliberately not claimed by these tests.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from echelle_spectra.reading_room import (
    BROWSE_ENDPOINT,
    DATA_FOLDER_MARKER,
    HOME_ENDPOINT,
    STEP_COMMANDS,
    build_reading_room,
    render_empty_campaign_page,
    render_setup_page,
)

#: The build stamps itself with the moment it ran, so two builds a second apart
#: differ in exactly that one place and nowhere else.
_STAMP = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\+00:00")

#: Every ``/api/`` string the served half is allowed to carry.
_ENDPOINTS = (BROWSE_ENDPOINT, HOME_ENDPOINT)


def _catalog(tmp_path: Path) -> Path:
    """One small merged index: one drive that answers, with one cube."""

    drive = tmp_path / "drive-a"
    drive.mkdir(exist_ok=True)
    (drive / "echelle-catalog.json").write_text("{}", encoding="utf-8")
    payload = {
        "schema": "echelle-merged-catalog/v1",
        "generated_at": "2026-08-14T00:00:00.000+00:00",
        "sources": [
            {
                "drive_id": "id-a",
                "volume_label": "NIFS-A",
                "drive_root": drive.as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": {"id": "run-1", "state": "completed", "counts": {"exported": 1},
                        "gate": "verdict"},
                "cubes": [
                    {
                        "path": "a.nc",
                        "shot_number": "193778",
                        "year": 2025,
                        "snapshot_id": "20250926_cmos",
                        "gate": "verdict",
                    }
                ],
            }
        ],
    }
    path = tmp_path / "all-years.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _api_strings(text: str) -> list[str]:
    """Every distinct endpoint-shaped string the page carries."""

    return sorted({match.group(0) for match in re.finditer(r"/api/[\w./?=-]*", text)})


@pytest.fixture()
def static_page(tmp_path: Path) -> str:
    return build_reading_room(_catalog(tmp_path), tmp_path / "web").read_text(encoding="utf-8")


@pytest.fixture()
def served_page(tmp_path: Path) -> str:
    return build_reading_room(
        _catalog(tmp_path), tmp_path / "web-served", served=True
    ).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# The static build does not move
# ---------------------------------------------------------------------------


def test_the_static_build_is_byte_identical_with_the_flag_absent_or_false(
    tmp_path: Path,
) -> None:
    catalog = _catalog(tmp_path)
    absent = build_reading_room(catalog, tmp_path / "web-a").read_bytes()
    explicit = build_reading_room(catalog, tmp_path / "web-b", served=False).read_bytes()
    # The build stamp is the one byte range two builds may legitimately differ
    # in; everything else must match exactly, including whitespace.
    assert _STAMP.sub("STAMP", absent.decode("utf-8")) == _STAMP.sub(
        "STAMP", explicit.decode("utf-8")
    )
    assert len(absent) == len(explicit)


def test_the_static_build_still_reaches_nothing(static_page: str) -> None:
    for forbidden in ("fetch(", "XMLHttpRequest", "http://", "/api/", "Browse"):
        assert forbidden not in static_page, f"the static page must not contain {forbidden!r}"
    # The one https:// the static page may carry is the documentation link,
    # and only as a plain anchor -- the page still fetches nothing.
    for index in range(len(static_page)):
        if static_page.startswith("https://", index):
            assert static_page[max(0, index - 9) : index] == '<a href="'
    assert "fetches nothing" in static_page
    assert 'class="picker"' not in static_page


# ---------------------------------------------------------------------------
# The served build
# ---------------------------------------------------------------------------


def test_the_served_build_puts_browse_beside_the_data_folder_field(served_page: str) -> None:
    assert '<button type="button" class="browse" data-browse="f-input">Browse' in served_page
    # The field itself is unchanged; the button is beside it, in one position.
    assert '<input type="text" id="f-input"' in served_page
    assert 'class="field-row"' in served_page
    # A control that looks pressable: it falls through to the page's own button
    # styling rather than stripping its resting border away.
    assert ".field-row .browse { flex: 0 0 auto; white-space: nowrap; }" in served_page


def test_the_served_build_carries_the_picker_dialog(served_page: str) -> None:
    assert '<div class="picker" id="picker" hidden>' in served_page
    assert 'role="dialog" aria-modal="true"' in served_page
    assert 'id="picker-path"' in served_page  # the header states the current path
    assert ">Choose this folder</button>" in served_page
    assert 'id="picker-body"' in served_page
    # The body scrolls inside itself and the page behind it is locked, with a
    # stable gutter so nothing — no rail — shifts when it opens.
    assert ".picker-body" in served_page
    assert "overscroll-behavior: contain" in served_page
    assert "html.picker-open { overflow: hidden; }" in served_page
    assert "scrollbar-gutter: stable" in served_page


def test_the_dialog_keeps_its_exit_in_reach_and_answers_escape(served_page: str) -> None:
    assert '<button type="button" class="picker-close" id="picker-close"' in served_page
    assert ">Close</button>" in served_page
    # The head holding Close never scrolls: only the body does.
    assert ".picker-head { flex: 0 0 auto;" in served_page
    assert "if (event.key === 'Escape' && !box.hidden) { pickerClose(); }" in served_page


def test_a_row_is_the_press_target_and_its_count_is_a_label(served_page: str) -> None:
    assert "row.className = 'picker-row';" in served_page
    assert "count.className = 'picker-count';" in served_page
    assert "' SIF file(s)'" in served_page
    # The count is a span written into the row, never its own control.
    assert "createElement('span')" in served_page
    assert "count.textContent = note;" in served_page


def test_choosing_writes_the_path_into_the_field_and_derives_from_it(served_page: str) -> None:
    assert "field.value = chosen;" in served_page
    assert "field.dispatchEvent(new Event('input', { bubbles: true }));" in served_page
    assert "compose();" in served_page
    assert "pickerClose();" in served_page


def test_the_picker_takes_a_pasted_path_and_names_what_it_is_reading(served_page: str) -> None:
    # A pasted path is navigation: the go box submits into pickerLoad.
    assert 'id="picker-typed"' in served_page
    assert "or paste a path and press Enter" in served_page
    assert "if (value) { pickerLoad(value); }" in served_page
    # A NAS can take seconds; the head says what is being read, immediately.
    assert "'reading ' + (path || 'the drives')" in served_page


def test_a_browse_error_is_one_sentence_inside_the_dialog(served_page: str) -> None:
    assert "alert(" not in served_page
    assert "confirm(" not in served_page
    assert 'class="picker-error" id="picker-error" hidden' in served_page
    assert (
        "pickerSay(String(answer.payload.error || 'That folder could not be read.'));"
        in served_page
    )
    # A refused folder still leaves somewhere to go: the drive list, whose own
    # empty path cannot refuse in turn, so the retry cannot loop.
    assert "if (path) { pickerLoad('', true); }" in served_page


def test_the_served_build_fetches_only_the_two_endpoints(served_page: str) -> None:
    assert "fetch(" in served_page
    assert BROWSE_ENDPOINT in served_page
    for found in _api_strings(served_page):
        assert found in _ENDPOINTS, f"the served page reaches an endpoint it should not: {found}"
    assert "http://" not in served_page
    # Same narrowed rule as the static page: https:// only inside the one
    # documentation anchor, never as a fetched resource.
    for index in range(len(served_page)):
        if served_page.startswith("https://", index):
            assert served_page[max(0, index - 9) : index] == '<a href="'
    assert "<script src" not in served_page and "<link " not in served_page


def test_the_served_banner_does_not_claim_to_reach_nothing(served_page: str) -> None:
    assert "reaches nothing outside itself" not in served_page
    assert "never executes commands" in served_page
    assert served_page.count("never executes commands") == 1


# ---------------------------------------------------------------------------
# The setup page: the page is the setup
# ---------------------------------------------------------------------------


def test_the_setup_page_names_the_situation_and_carries_the_picker() -> None:
    text = render_setup_page()
    assert "No campaign home yet" in text
    assert '<div class="picker" id="picker" hidden>' in text
    assert ">Choose this folder</button>" in text
    assert 'id="pick-home"' in text
    # One sentence of teaching, not a lecture: the page says what to do once.
    assert text.count("Pick the folder this campaign lives in") == 1


def test_the_setup_page_writes_the_home_and_then_reloads() -> None:
    text = render_setup_page()
    assert f"fetch('{HOME_ENDPOINT}'" in text
    assert "method: 'POST'" in text
    assert "JSON.stringify({ folder: folder })" in text
    assert "window.location.href = '/';" in text
    # A refused folder says so in place, inside the dialog.
    assert "That folder could not become the campaign home." in text
    assert "alert(" not in text
    for found in _api_strings(text):
        assert found in _ENDPOINTS, f"the setup page reaches an endpoint it should not: {found}"


def test_the_setup_page_bakes_in_no_example_drive_path() -> None:
    text = render_setup_page()
    for invented in ("C:\\", "C:/", "D:\\", "D:/", "/mnt/", "/media/", "Users\\", "Users/"):
        assert invented not in text, f"the setup page invents a path: {invented!r}"
    # The drive list is what the server answers with, not a literal in the page.
    assert "PICKER_BROWSE + encodeURIComponent(path || '')" in text


# ---------------------------------------------------------------------------
# The empty campaign page: empty is not broken
# ---------------------------------------------------------------------------


def _home() -> dict[str, str]:
    return {
        "folder": "E:/2026-campaign",
        "registry": "E:/2026-campaign/calibration_registry.toml",
        "calibrations": "E:/2026-campaign/calibrations",
        "catalog": "E:/2026-campaign/echelle-catalog.json",
    }


def test_the_empty_campaign_page_states_empty_rather_than_broken() -> None:
    text = render_empty_campaign_page(_home())
    read = text[text.index('<main class="solo">') : text.index("</main>")]
    assert "This campaign has no catalog yet" in read
    assert "That is empty, not broken." in read
    # Nothing is claimed about cubes or drives, because nothing is known.
    for claim in ("cube(s)", "drive(s)", "state-measured", "verdict-authorized", "aligned"):
        assert claim not in read, f"the empty page claims {claim!r}"
    assert "as soon as that first run writes the catalog" in read


def test_the_empty_campaign_page_composes_the_first_command_in_both_shells() -> None:
    text = render_empty_campaign_page(_home())
    _, template = STEP_COMMANDS["sample"]
    assert template.split()[0] == "echelle"
    assert 'data-shell="powershell"' in text and 'data-shell="posix"' in text
    # The registry, snapshot root and catalog come from the home; the data
    # folder is the one marker, in both shell shapes.
    assert text.count("&lt;data folder&gt;") >= 2
    assert DATA_FOLDER_MARKER == "<data folder>"
    assert "E:\\2026-campaign\\calibration_registry.toml" in text  # PowerShell shape
    assert "E:/2026-campaign/calibration_registry.toml" in text  # POSIX shape
    assert "--sample auto" in text
    # One copy control per shape, carrying the whole command.
    assert text.count('class="copy" data-copy="echelle process') == 2


def test_the_empty_campaign_page_wires_the_controls_it_renders() -> None:
    """A show/hide toggle and a copy button that no script answers are dead
    controls that still look pressable, so the one-screen pages carry their own
    copy of that wiring rather than borrowing the campaign page's."""

    text = render_empty_campaign_page(_home())
    assert 'class="fold-toggle" aria-expanded="false"' in text
    assert "function soloCopy(button)" in text
    assert "event.target.closest('.fold-toggle')" in text
    assert "event.target.closest('.copy')" in text


def test_the_empty_campaign_page_falls_back_to_the_home_folder_alone() -> None:
    text = render_empty_campaign_page({"folder": "E:/2026-campaign"})
    assert "E:/2026-campaign/calibrations" in text
    assert "E:/2026-campaign/echelle-catalog.json" in text
    assert "&lt;data folder&gt;" in text
