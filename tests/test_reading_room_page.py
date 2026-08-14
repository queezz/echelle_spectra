"""Packet F9 — the page flows as the work does, tab by tab and drive by drive.

Every assertion here reads the one static file ``echelle web`` writes.  The
live verification the law also demands — the Perimeter Walk, rail
offset-identity sampling in a real DOM, and the visual review — is the owner's
validation-week work and is deliberately not claimed by these tests.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import pytest

from echelle_spectra.calibration_registry import REGISTRY_SCHEMA
from echelle_spectra.reading_room import _SOURCE_NOTES, build_reading_room, render_markdown
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot

DRIFT_SCHEMA = "echelle-drift-evidence/v1"


def _catalog(tmp_path: Path) -> Path:
    """One merged index holding all four honesty states at once."""

    drive = tmp_path / "drive-a"
    drive.mkdir(exist_ok=True)
    (drive / "echelle-catalog.json").write_text("{}", encoding="utf-8")
    cube = {
        "path": "a.nc",
        "shot_number": "193778",
        "year": 2025,
        "snapshot_id": "20250926_cmos",
        "wavelength_min_nm": 400.0,
        "wavelength_max_nm": 700.0,
    }
    payload = {
        "schema": "echelle-merged-catalog/v1",
        "generated_at": "2026-08-14T00:00:00.000+00:00",
        "sources": [
            {
                "drive_id": "id-a",
                "volume_label": "NIFS-A",
                "drive_root": drive.as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": {
                    "id": "run-1",
                    "state": "completed",
                    "counts": {"exported": 2},
                    "gate": "verdict",
                },
                "cubes": [
                    {**cube, "gate": "verdict"},
                    {**cube, "path": "b.nc", "shot_number": "193779", "gate": "sample"},
                ],
            },
            {
                # Missing drive: the row survives, the files did not answer.
                "drive_id": "id-b",
                "volume_label": "NIFS-B",
                "drive_root": (tmp_path / "gone").as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": {"id": "run-0", "state": "completed", "counts": {}, "gate": "verdict"},
                "cubes": [{**cube, "path": "c.nc", "gate": "unrecorded (pre-gate receipt)"}],
            },
            {
                # Unmeasured: reachable, but no receipt ever described it.
                "drive_id": "id-c",
                "volume_label": "NIFS-C",
                "drive_root": drive.as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": None,
                "cubes": [],
            },
            {
                # Empty: measured, and it published nothing.
                "drive_id": "id-d",
                "volume_label": "NIFS-D",
                "drive_root": drive.as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": {
                    "id": "run-2",
                    "state": "completed",
                    "counts": {},
                    "gate": "ungated (no registry)",
                },
                "cubes": [],
            },
        ],
    }
    path = tmp_path / "all-years.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _shifted_evidence(tmp_path: Path) -> Path:
    payload = {
        "schema": DRIFT_SCHEMA,
        "created_at": "2026-08-14T01:00:00+00:00",
        "verdict": "shifted",
        "snapshot_ids": ["20250926_cmos"],
        "sampled_cubes": [{"cube": "a.nc"}, {"cube": "b.nc"}],
        "skipped_cubes": [
            {"cube": "c.nc", "shot_number": "12", "reason": "no plasma-bright frames"}
        ],
        "interval_warning": "residuals form two groups (shots 1..2 vs 3..4)",
        "thresholds_px": {"alignment_maximum_residual": 0.5, "repair_limit": 25},
        "summary": {
            "median_shift_px": 8.0,
            "maximum_absolute_pixel_residual_px": 8.4,
            "maximum_pixel_deviation_px": 0.4,
            "alignment_tolerance_px": 0.5,
            "quorum": {"satisfied": True, "resolved_lines": 8, "distinct_orders": 4,
                       "coverage_nm": [410.0, 660.0]},
        },
        "per_shot": [
            {
                "shot_number": "193778",
                "cube": "a.nc",
                "date": "2025-09-26",
                "lines": 4,
                "median_shift_px": 8.0,
                "pixel_spread_px": 0.2,
                "group": 1,
            }
        ],
        "order_corrections": [
            {
                "order": 6,
                "reference_pixel": 511.5,
                "dispersion_nm_per_px": 0.0108,
                "predicted_shift_nm": -0.0864,
            }
        ],
        "repair_commands": [
            {
                "shell": "any",
                "purpose": "accept the sampled shift and emit the immutable -rN snapshot",
                "command": 'echelle drift refine "drift.json" --accept-shift 8',
            },
            {"shell": "any", "purpose": "repoint the registry [validity] entry", "command": ""},
            {
                "shell": "powershell",
                "purpose": "revise every cube already exported",
                "command": 'Get-ChildItem "cubes\\*.nc"',
            },
            {
                "shell": "posix",
                "purpose": "revise every cube already exported",
                "command": "for cube in cubes/*.nc; do :; done",
            },
        ],
        "lines": [
            {
                "shot_number": "193778",
                "line": "H-alpha",
                "status": "measured",
                "expected_nm": 656.279,
                "residual_nm": 0.086,
                "pixel_residual_px": 8.0,
            }
        ],
    }
    path = tmp_path / "drift-evidence.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _unknown_evidence(tmp_path: Path) -> Path:
    path = tmp_path / "unknown.json"
    path.write_text(
        json.dumps({"schema": DRIFT_SCHEMA, "verdict": "probably-fine", "lines": []}),
        encoding="utf-8",
    )
    return path


@pytest.fixture()
def page(tmp_path: Path) -> str:
    built = build_reading_room(
        _catalog(tmp_path),
        tmp_path / "web",
        drift_paths=[_shifted_evidence(tmp_path), _unknown_evidence(tmp_path)],
    )
    assert not (built.parent / "reading-room.json").exists(), "the dead sidecar is gone"
    return built.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Helpers that slice one region of the built page
# ---------------------------------------------------------------------------


def _element(text: str, marker: str, tag: str = "section") -> str:
    """Slice the whole element carrying *marker*, nesting included."""

    start = text.rfind(f"<{tag}", 0, text.index(marker) + len(marker))
    depth = 0
    for match in re.finditer(rf"<{tag}\b|</{tag}>", text[start:]):
        depth += -1 if match.group().startswith("</") else 1
        if depth == 0:
            return text[start : start + match.end()]
    raise AssertionError(f"unbalanced <{tag}> around {marker}")


def _escaped(text: str) -> str:
    return html.escape(text, quote=True)


def _rail(text: str, identifier: str) -> str:
    return _element(text, f'id="{identifier}"', "aside")


def _rail_group(text: str, identifier: str, tab: str) -> str:
    rail = _rail(text, identifier)
    start = rail.index(f'class="rail-group" data-tab="{tab}"')
    end = rail.find('class="rail-group"', start + 1)
    return rail[start:] if end == -1 else rail[start:end]


def _view(text: str, tab: str) -> str:
    rail = text[text.index('id="tab-' + tab + '"') :]
    end = rail.find('<section class="tabview"')
    return rail if end == -1 else rail[:end]


def _drive_rows(text: str) -> list[str]:
    return [
        _element(text, match.group(), "article")
        for match in re.finditer(r'<article class="drive-row" id="drive-\d+"', text)
    ]


def _steps(row: str) -> list[tuple[str, str, bool]]:
    """Return (name, state class, primary) for every step box in one flow row."""

    found = []
    for match in re.finditer(
        r'<li class="step ([a-z-]+)( is-primary)?"[^>]*>.*?<span class="step-name">([^<]+)<',
        row,
        flags=re.S,
    ):
        found.append((match.group(3), match.group(1), bool(match.group(2))))
    return found


def _first_not_done(steps: list[tuple[str, str, bool]]) -> tuple[str, str, bool]:
    return next(item for item in steps if item[1] in {"step-ready", "step-blocked"})


# ---------------------------------------------------------------------------
# Fixtures at three pipeline positions
# ---------------------------------------------------------------------------


def _one_drive_catalog(tmp_path: Path, sources: list[dict[str, object]]) -> Path:
    path = tmp_path / "index.json"
    path.write_text(
        json.dumps(
            {
                "schema": "echelle-merged-catalog/v1",
                "generated_at": "2026-08-14T00:00:00.000+00:00",
                "sources": sources,
            }
        ),
        encoding="utf-8",
    )
    return path


def _connected(tmp_path: Path, label: str, **extra: object) -> dict[str, object]:
    drive = tmp_path / label
    drive.mkdir(exist_ok=True)
    (drive / "echelle-catalog.json").write_text("{}", encoding="utf-8")
    source: dict[str, object] = {
        "drive_id": f"id-{label}",
        "volume_label": label,
        "drive_root": drive.as_posix(),
        "catalog_path": "echelle-catalog.json",
        "run": None,
        "cubes": [],
    }
    source.update(extra)
    return source


def _registry(tmp_path: Path, *, covering: bool) -> tuple[Path, Path]:
    """A real registry whose one epoch does or does not cover today."""

    snapshots = tmp_path / "calibrations"
    sources = tmp_path / "snapshot-sources"
    sources.mkdir()
    files = {}
    for role in ROLE_FILENAMES:
        item = sources / f"{role}.dat"
        item.write_text(f"20260812_cmos/{role}\n", encoding="utf-8")
        files[role] = item
    create_snapshot(
        snapshots,
        snapshot_id="20260812_cmos",
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity=(
            {"date_from": "2020-01-01", "date_to": "2099-12-31"}
            if covering
            else {"date_from": "2020-01-01", "date_to": "2020-12-31"}
        ),
    )
    registry = tmp_path / "calibration_registry.toml"
    registry.write_text(
        f'schema = "{REGISTRY_SCHEMA}"\n\n[[epochs]]\nsnapshot_id = "20260812_cmos"\n',
        encoding="utf-8",
    )
    return registry, snapshots


def _aligned_evidence(tmp_path: Path, cube: str) -> Path:
    path = tmp_path / f"aligned-{cube}.json"
    path.write_text(
        json.dumps(
            {
                "schema": DRIFT_SCHEMA,
                "created_at": "2026-08-14T02:00:00+00:00",
                "verdict": "aligned",
                "snapshot_ids": ["20260812_cmos"],
                "sampled_cubes": [{"cube": cube}],
                "per_shot": [{"shot_number": "1", "cube": cube, "lines": 6}],
                "lines": [],
            }
        ),
        encoding="utf-8",
    )
    return path


def _aligned_awaiting_bulk(tmp_path: Path) -> str:
    """One drive whose sample was audited aligned: the product is next."""

    source = _connected(
        tmp_path,
        "NIFS-A",
        run={"id": "run-a", "state": "completed", "counts": {"exported": 1}, "gate": "sample",
             "sample": True},
        cubes=[
            {
                "path": "sample-a.nc",
                "shot_number": "1",
                "year": 2026,
                "snapshot_id": "20260812_cmos",
                "gate": "sample",
            }
        ],
    )
    return build_reading_room(
        _one_drive_catalog(tmp_path, [source]),
        tmp_path / "web",
        drift_paths=[_aligned_evidence(tmp_path, "sample-a.nc")],
    ).read_text(encoding="utf-8")


def _nothing_done(tmp_path: Path) -> str:
    catalog = _one_drive_catalog(tmp_path, [_connected(tmp_path, "NIFS-A")])
    return build_reading_room(catalog, tmp_path / "web").read_text(encoding="utf-8")


def _calibrated_only(tmp_path: Path) -> str:
    catalog = _one_drive_catalog(tmp_path, [_connected(tmp_path, "NIFS-A")])
    registry, snapshots = _registry(tmp_path, covering=True)
    return build_reading_room(
        catalog,
        tmp_path / "web",
        registry_path=registry,
        calibrations_root=snapshots,
    ).read_text(encoding="utf-8")


def _mid_drives(tmp_path: Path) -> str:
    cube = {
        "shot_number": "1",
        "year": 2026,
        "snapshot_id": "20260812_cmos",
        "wavelength_min_nm": 400.0,
        "wavelength_max_nm": 700.0,
        "gate": "sample",
    }
    sampled = _connected(
        tmp_path,
        "NIFS-A",
        run={"id": "run-a", "state": "completed", "counts": {"exported": 1}, "gate": "sample",
             "sample": True},
        cubes=[{**cube, "path": "sample-a.nc"}],
    )
    finished = _connected(
        tmp_path,
        "NIFS-B",
        run={"id": "run-b", "state": "completed", "counts": {"exported": 9}, "gate": "verdict"},
        cubes=[{**cube, "path": "bulk-b.nc", "gate": "verdict"}],
    )
    catalog = _one_drive_catalog(tmp_path, [sampled, finished])
    registry, snapshots = _registry(tmp_path, covering=True)
    return build_reading_room(
        catalog,
        tmp_path / "web",
        drift_paths=[_aligned_evidence(tmp_path, "bulk-b.nc")],
        registry_path=registry,
        calibrations_root=snapshots,
    ).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Tabs, rails and anchors
# ---------------------------------------------------------------------------


def test_four_tabs_by_work_type_with_their_own_views(page: str) -> None:
    for key, title in (
        ("now", "Now"),
        ("drives", "Drives"),
        ("calibration", "Calibration"),
        ("reading", "Reading room"),
    ):
        assert f'class="tab" data-tab="{key}"' in page
        assert f">{title}</button>" in page
        assert f'id="tab-{key}" data-tab="{key}"' in page
    # Now leads; every other view opens closed.
    assert 'id="tab-now" data-tab="now">' in page
    for key in ("drives", "calibration", "reading"):
        assert f'id="tab-{key}" data-tab="{key}" hidden>' in page


def test_pressing_the_active_tab_navigates_home_from_a_sub_view(page: str) -> None:
    # The cookbook's corrected guard: the no-render shortcut survives, with the
    # sub-view named as its exception rather than the guard deleted.
    assert "if (name === state.tab && !subviewOpen(name)) { return; }" in page
    assert "if (name === state.tab) { closeFolds(byId('tab-' + name)); }" in page
    # Every sub-view a tab can hold is what the exception actually tests for.
    assert "view.querySelector('details[open]')" in page
    assert "view.querySelector('.fold-toggle[aria-expanded=\"true\"]')" in page
    # A tab press lands at the top of its destination.
    assert "window.scrollTo(0, 0);" in page


def test_each_rail_carries_only_its_own_tab_cargo(page: str) -> None:
    now_left = _rail_group(page, "rail-left", "now")
    drives_left = _rail_group(page, "rail-left", "drives")
    for field in ("f-input", "f-output", "f-registry", "f-verdict", "f-epoch"):
        assert f'id="{field}"' in now_left
    assert 'id="compose"' in now_left
    for control in ("filter-year", "filter-epoch", "filter-drive", "filter-status"):
        assert f'id="{control}"' in drives_left and f'id="{control}"' not in now_left
    # The composer lives once, in Now; Drives keeps a compact one-press entry.
    assert 'id="send-drive"' in drives_left and 'id="f-input"' not in drives_left
    # Calibration and the reading room own no controls; the rail stays empty
    # rather than inventing navigation for itself.
    for tab in ("calibration", "reading"):
        assert "<article" not in _rail_group(page, "rail-left", tab)
    now_right = _rail_group(page, "rail-right", "now")
    drives_right = _rail_group(page, "rail-right", "drives")
    assert 'id="sectnav-now"' in now_right and "Campaign position" in now_right
    assert 'id="find"' in drives_right and 'id="reset"' in drives_right
    assert 'id="find"' not in now_right
    assert 'id="sectnav-calibration"' in _rail_group(page, "rail-right", "calibration")
    assert 'id="sectnav-reading"' in _rail_group(page, "rail-right", "reading")


def test_two_rails_are_sticky_at_their_own_resting_offset(page: str) -> None:
    assert 'class="rail rail-left" id="rail-left"' in page
    assert 'class="rail rail-right" id="rail-right"' in page
    # The offset is the sum of the page's own header metrics, never a literal.
    assert "top: calc(var(--bar) + var(--gap));" in page
    # The grid parent stretches and only the rail pins: the zero-travel trap.
    assert "align-items: stretch;" in page
    assert "align-self: start;" in page
    assert "@media (max-width: 900px)" in page


def test_every_jump_target_is_a_unique_sec_anchor(page: str) -> None:
    ids = re.findall(r'id="([^"]+)"', page)
    assert len(ids) == len(set(ids)), "duplicate element ids"
    anchors = [item for item in ids if item.startswith("sec-")]
    assert {
        "sec-now-calibrate",
        "sec-now-drives",
        "sec-plan",
        "sec-drives-cards",
        "sec-catalog",
        "sec-cal-epochs",
        "sec-drift",
        "sec-reading-room",
    } <= set(anchors)
    assert '[id^="sec-"] { scroll-margin-top: calc(var(--bar) + 16px); }' in page
    for anchor in anchors:
        assert f'href="#{anchor}"' in page, f"{anchor} has no rail link"
    # The scroll-spy answers for one tab's own view, never the whole column.
    assert "var view = byId('tab-' + state.tab);" in page
    assert "var targets = view.querySelectorAll('[id^=\"sec-\"]');" in page
    assert "document.querySelectorAll('#sectnav-' + state.tab + ' .sectnav-link')" in page


def test_local_find_stays_with_the_catalog_table(page: str) -> None:
    assert 'id="find"' in page and 'id="cube-table"' in page
    assert 'id="cube-table"' in _view(page, "drives")
    assert "function filterCatalog()" in page
    assert "byId('find').addEventListener('input', filterCatalog);" in page
    assert "data-find=" in page


# ---------------------------------------------------------------------------
# The stepper
# ---------------------------------------------------------------------------


def test_nothing_done_asks_for_the_bench_then_the_first_sample(tmp_path: Path) -> None:
    text = _nothing_done(tmp_path)
    calibrate = _element(text, 'id="sec-now-calibrate"')
    steps = _steps(calibrate)
    assert [name for name, _, _ in steps] == [
        "Sphere + lamps",
        "Bench fit",
        "Snapshot saved",
        "Registry epoch",
    ]
    name, state, primary = _first_not_done(steps)
    assert (name, state, primary) == ("Bench fit", "step-ready", True)
    assert "echelle-calib" in calibrate
    # Nothing this page cannot see is claimed as done.
    assert ("Sphere + lamps", "step-unrecorded", False) in steps
    row = _drive_rows(text)[0]
    name, state, primary = _first_not_done(_steps(row))
    assert (name, state, primary) == ("Sample N", "step-ready", True)
    assert "--sample 20" in row


def test_calibrated_only_closes_the_calibrate_stage(tmp_path: Path) -> None:
    text = _calibrated_only(tmp_path)
    calibrate = _element(text, 'id="sec-now-calibrate"')
    steps = _steps(calibrate)
    assert [state for _, state, _ in steps] == ["step-done"] * 4
    assert "is-primary" not in calibrate
    assert "The calibration is in place for this campaign." in calibrate
    assert "20260812_cmos covers today" in calibrate
    # The drive has not moved: its own first step is still the sample.
    row = _drive_rows(text)[0]
    assert _first_not_done(_steps(row))[:2] == ("Sample N", "step-ready")


def test_parallel_drives_each_carry_their_own_position(tmp_path: Path) -> None:
    text = _mid_drives(tmp_path)
    rows = _drive_rows(text)
    assert len(rows) == 2, "one independent stepper row per connected drive"
    sampled, finished = rows
    assert 'data-drive="NIFS-A"' in sampled and 'data-drive="NIFS-B"' in finished
    name, state, primary = _first_not_done(_steps(sampled))
    assert (name, state, primary) == ("Drift audit", "step-ready", True)
    assert "echelle drift audit" in sampled
    # The finished drive has nothing waiting: cubes exist and are catalogued.
    assert not [item for item in _steps(finished) if item[1] in {"step-ready", "step-blocked"}]
    assert "is-primary" not in finished
    assert "This drive is done" in finished
    assert ("Generate cubes", "step-done", False) in _steps(finished)
    # The optional text export is never what makes a drive done.
    assert ("LHD txt", "step-unrecorded", False) in _steps(finished)


def test_a_done_step_links_the_evidence_it_rests_on(page: str) -> None:
    row = next(item for item in _drive_rows(page) if 'data-drive="NIFS-A"' in item)
    assert '<a class="xlink" href="#drift-1" data-tab="calibration"' in row
    # The label is the evidence file's own name; the whole path stays in the title.
    assert ">drift-evidence.json</a>" in row and 'title="' in row
    assert 'id="drift-1"' in _view(page, "calibration")
    assert ("Drift audit", "step-done", False) in _steps(row)


def test_a_ready_step_carries_both_shell_shapes_and_a_full_payload(tmp_path: Path) -> None:
    row = _drive_rows(_nothing_done(tmp_path))[0]
    block = _element(row, 'id="now-d1-sample"', "article")
    assert "PowerShell" in block and "POSIX shell" in block
    payloads = [html.unescape(item) for item in re.findall(r'data-copy="([^"]*)"', block)]
    bodies = [
        html.unescape(item)
        for item in re.findall(r'<pre class="cmd-body"[^>]*>(.*?)</pre>', block, flags=re.S)
    ]
    assert len(payloads) == len(bodies) == 2
    assert payloads == bodies
    assert all(payload.startswith("echelle process ") for payload in payloads)
    assert "\\NIFS-A" in payloads[0] and "/NIFS-A" in payloads[1]


def test_a_shifted_verdict_offers_its_own_repair_step(page: str) -> None:
    row = next(item for item in _drive_rows(page) if 'data-drive="NIFS-A"' in item)
    name, state, primary = _first_not_done(_steps(row))
    assert (name, state, primary) == ("Verdict", "step-ready", True)
    assert "shifted — refine, then repoint the registry" in row
    assert "echelle drift refine" in row


def test_a_remembered_drive_keeps_one_collapsed_line(page: str) -> None:
    assert 'class="drive-absent"' in page
    absent = _element(page, 'class="drive-absent"', "p")
    assert "NIFS-B" in absent and "cube(s) remembered" in absent
    assert "next: Connect + identify" in absent
    # A drive that did not answer never renders a full stepper row.
    assert not any('data-drive="NIFS-B"' in row for row in _drive_rows(page))


def test_the_stepper_is_drawn_with_the_page_own_css(page: str) -> None:
    assert '<ol class="flow">' in page
    assert ".step + .step::before" in page
    assert "mermaid" not in page.lower()


# ---------------------------------------------------------------------------
# The compression pass
# ---------------------------------------------------------------------------


def test_drive_cards_carry_facts_and_the_teaching_lives_once(page: str) -> None:
    cards = _element(page, 'id="sec-drives-cards"')
    for sentence in _SOURCE_NOTES.values():
        assert _escaped(sentence) not in cards, "a drive card is teaching again"
        assert page.count(_escaped(sentence)) == 1, "the teaching must be said exactly once"
    legend = _rail_group(page, "rail-right", "drives")
    assert _escaped(_SOURCE_NOTES["missing-drive"]) in legend
    # What survives on the card is chips, counts, and one truncated path.
    assert 'class="chips"' in cards and "cube(s)" in cards
    assert 'class="chip path" title=' in cards
    assert "<p class=\"note\">" not in cards


def test_the_read_only_sentence_is_said_once_in_the_banner(page: str) -> None:
    assert page.count("never executes commands") == 1
    assert "This page has no code path that runs" not in page
    # The composer preamble is one line, not a paragraph of instruction.
    assert "Pre-filled from this catalog and registry; Compose rewrites text only." in page
    assert "Nothing here runs; Compose only rewrites" not in page


def test_data_reads_at_the_page_base_size(page: str) -> None:
    assert "table { border-collapse: collapse; width: 100%; font-size: 1rem; }" in page
    assert ".muted { color: var(--muted); font-size: .92rem; }" in page


def test_the_product_is_named_as_the_product(tmp_path: Path, page: str) -> None:
    assert "Generate cubes" in page
    # LHD text is the side deliverable, and no file records it.
    assert "the side deliverable: no receipt or catalog field records a txt export" in page
    # An aligned verdict makes the product the very next thing to do.
    row = _drive_rows(_aligned_awaiting_bulk(tmp_path))[0]
    name, state, primary = _first_not_done(_steps(row))
    assert (name, state, primary) == ("Generate cubes", "step-ready", True)
    assert _escaped("Generate the cubes — the campaign's product") in row
    assert "echelle process --plan" in row


# ---------------------------------------------------------------------------
# Surviving F7 contracts
# ---------------------------------------------------------------------------


def test_every_composed_command_carries_both_shell_shapes_and_a_full_payload(
    tmp_path: Path, page: str
) -> None:
    for identifier in ("cmd-process", "cmd-audit", "cmd-bench"):
        block = _element(page, f'id="{identifier}"', "article")
        assert "PowerShell" in block and "POSIX shell" in block
        assert f'id="{identifier}-powershell"' in block and f'id="{identifier}-posix"' in block
        payloads = [html.unescape(item) for item in re.findall(r'data-copy="([^"]*)"', block)]
        bodies = [
            html.unescape(item)
            for item in re.findall(r'<pre class="cmd-body"[^>]*>(.*?)</pre>', block, flags=re.S)
        ]
        assert len(payloads) == len(bodies) == 2
        # The copy button carries the complete command, collapsed or not.
        assert payloads == bodies
        assert all(payload.strip() for payload in payloads)
    audit = page[page.index('id="cmd-audit"') :]
    windows = html.unescape(re.findall(r'data-copy="([^"]*)"', audit)[0])
    posix = html.unescape(re.findall(r'data-copy="([^"]*)"', audit)[1])
    assert windows.startswith("echelle drift audit") and posix.startswith("echelle drift audit")
    assert "\\drive-a" in windows and "/drive-a" in posix
    process = page[page.index('id="cmd-process"') :]
    assert (
        html.unescape(re.findall(r'data-copy="([^"]*)"', process)[0])
        == 'echelle process --plan "campaign-plan.toml"'
    )
    # The composed plan is what the composed process command actually reads.
    plan = html.unescape(
        re.search(r'id="plan-out"[^>]*>(.*?)</textarea>', page, flags=re.S).group(1)
    )
    assert 'drift_verdict = "' in plan and 'central_index = "' in plan
    assert "20250926_cmos" in plan and "{{" not in plan


def test_the_four_honesty_states_and_an_unknown_verdict_stay_distinct(page: str) -> None:
    for state, label in (
        ("state-missing-drive", "missing drive"),
        ("state-unmeasured", "unmeasured"),
        ("state-empty", "empty"),
        ("state-insufficient-data", "insufficient-data"),
    ):
        assert f'<span class="pill {state}">{label}</span>' in page
        assert f".card.{state} {{ border-left" in page
    # An unknown verdict is its own state; it never borrows a known one.
    assert '<span class="pill state-unrecognized">unrecognized: probably-fine</span>' in page
    unknown = page[page.index("unrecognized: probably-fine") - 200 :]
    unknown = unknown[: unknown.index("</article>")]
    for known in ("verdict-aligned", "verdict-shifted", "verdict-beyond-repair", "insufficient"):
        assert known not in unknown
    # A sampled drive reads differently from a verdict-authorized one and from
    # a run the registry never gated.
    for gate, label in (
        ("gate-verdict", "verdict-authorized"),
        ("gate-sample", "unverified sample"),
        ("gate-ungated", "ungated (no registry)"),
        ("gate-unrecorded", "unrecorded (pre-gate receipt)"),
    ):
        assert f'<span class="pill {gate}">{label}</span>' in page


def test_drift_evidence_v2_is_rendered_in_full(page: str) -> None:
    assert "Interval warning" in page and "residuals form two groups" in page
    assert 'class="warn interval-warning"' in page
    assert "<th>Median Δx px</th>" in page and "<th>Spread px</th>" in page
    assert "<td>193778</td>" in page and "<td>2025-09-26</td>" in page
    assert "Per-order correction the refinement would apply" in page
    assert "Thresholds this verdict was judged against" in page
    assert "alignment maximum residual" in page
    assert "Skipped cubes" in page and "no plasma-bright frames" in page
    assert "Quorum satisfied" in page
    # The stored repair sequence becomes meaning-first rows, both shells, and a
    # step with no command says so rather than inventing one.
    assert "Repair sequence composed from this evidence" in page
    assert "accept the sampled shift" in page
    assert "No command: this step is a decision" in page
    assert "Get-ChildItem" in page and "for cube in cubes/*.nc" in page
    # The per-line drill-down sits behind a fold and keeps an exit in reach.
    assert 'class="drill"' in page and "Per-line evidence (1 row(s))" in page
    assert "drill-close" in page and "Escape closes it" in page
    assert "if (event.key === 'Escape') { closeFolds(); }" in page
    # Evidence is read in sequence position, inside the Calibration tab.
    assert 'class="drill"' in _view(page, "calibration")


def test_unmeasured_drift_is_not_rendered_as_aligned(tmp_path: Path) -> None:
    built = build_reading_room(_catalog(tmp_path), tmp_path / "web")
    text = built.read_text(encoding="utf-8")
    assert "Bulk readiness is unmeasured" in text
    assert "aligned</h3>" not in text


def test_markdown_is_rendered_rather_than_escaped() -> None:
    rendered = render_markdown(
        "# Heading\n\nA **bold** word and a [link](procedure.md).\n\n- first\n- second\n\n"
        "Inline `code` here.\n\n```\nraw text\n```\n"
    )
    assert '<h3 id="doc-heading">Heading</h3>' in rendered
    assert "<strong>bold</strong>" in rendered
    assert '<a href="procedure.md">link</a>' in rendered
    assert "<ul><li>first</li><li>second</li></ul>" in rendered
    assert "<code>code</code>" in rendered
    assert "<pre><code>raw text</code></pre>" in rendered
    # None of the source markers survive into the output.
    assert "# Heading" not in rendered and "**" not in rendered
    assert "[link]" not in rendered


def test_document_headings_stay_out_of_the_sec_anchor_namespace() -> None:
    rendered = render_markdown("## Section\n", id_prefix="doc-0-")
    assert 'id="doc-0-section"' in rendered
    assert "sec-" not in rendered


def test_packaged_documents_render_with_no_checkout_in_reach(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    catalog = _catalog(tmp_path)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    text = build_reading_room(catalog, tmp_path / "web").read_text(encoding="utf-8")
    for name in ("Vocabulary", "Procedure", "Provenance"):
        assert f'<h3 class="card-title">{name}</h3>' in text
    assert "the installed package (resources/reading_room/vocabulary.md)" in text
    # The canon is rendered, not escaped into a code block.
    assert "# Vocabulary" not in text
    assert "misaligned-beyond-repair" in text
    # The reading room is its own tab at last.
    assert 'id="tab-reading"' in text
    assert "Vocabulary" in _view(text, "reading")


def test_extra_documents_are_rendered_after_the_packaged_canon(tmp_path: Path) -> None:
    extra = tmp_path / "trip-notes.md"
    extra.write_text("# Trip notes\n\nOne *line*.\n", encoding="utf-8")
    text = build_reading_room(
        _catalog(tmp_path), tmp_path / "web", document_paths=[extra]
    ).read_text(encoding="utf-8")
    assert '<h3 class="card-title">Trip notes</h3>' in text
    assert "<em>line</em>" in text
    assert text.index("Provenance") < text.index("Trip notes")


def test_the_registry_pre_fills_the_composer_or_says_it_was_not_read(tmp_path: Path) -> None:
    text = build_reading_room(_catalog(tmp_path), tmp_path / "web").read_text(encoding="utf-8")
    assert "no registry supplied" in text
    missing = tmp_path / "absent-registry.toml"
    text = build_reading_room(
        _catalog(tmp_path), tmp_path / "web2", registry_path=missing
    ).read_text(encoding="utf-8")
    assert "registry unreadable" in text
    assert "calibration registry not found" in text
    # An unreadable registry blocks the calibrate stage; it never reads as done.
    assert "step step-blocked" in _element(text, 'id="sec-now-calibrate"')


def test_the_page_executes_nothing_and_reaches_nothing(page: str) -> None:
    assert "never executes commands" in page
    for forbidden in (
        "subprocess",
        "fetch(",
        "XMLHttpRequest",
        "http://",
        "https://",
        "<script src",
        "<link ",
        "@import",
    ):
        assert forbidden not in page, f"the page must not contain {forbidden!r}"
    assert "color-scheme: light dark;" in page
    assert "@media (prefers-color-scheme: dark)" in page
    # Every control keeps a resting border and fill.
    assert "button, select, input, textarea {" in page
