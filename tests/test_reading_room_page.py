"""Packet F7 — the reading room page is built to the house web UI law.

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

from echelle_spectra.reading_room import build_reading_room, render_markdown

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


def _rail(text: str, identifier: str) -> str:
    start = text.index(f'id="{identifier}"')
    return text[start : text.index("</aside>", start)]


def test_two_rails_are_sticky_at_their_own_resting_offset(page: str) -> None:
    assert 'class="rail rail-left" id="rail-left"' in page
    assert 'class="rail rail-right" id="rail-right"' in page
    # The offset is the sum of the page's own header metrics, never a literal.
    assert "top: calc(var(--bar) + var(--gap));" in page
    # The grid parent stretches and only the rail pins: the zero-travel trap.
    assert "align-items: stretch;" in page
    assert "align-self: start;" in page
    assert "@media (max-width: 900px)" in page


def test_controls_live_left_and_context_lives_right(page: str) -> None:
    left, right = _rail(page, "rail-left"), _rail(page, "rail-right")
    for control in ("filter-year", "filter-epoch", "filter-drive", "filter-status"):
        assert f'id="{control}"' in left
    for field in ("f-input", "f-output", "f-registry", "f-verdict", "f-epoch"):
        assert f'id="{field}"' in left
    assert 'id="compose"' in left
    assert 'id="find"' in right and 'id="sectnav"' in right
    assert 'id="find"' not in left


def test_every_jump_target_is_a_unique_sec_anchor(page: str) -> None:
    ids = re.findall(r'id="([^"]+)"', page)
    assert len(ids) == len(set(ids)), "duplicate element ids"
    anchors = [item for item in ids if item.startswith("sec-")]
    assert {"sec-catalog", "sec-drift", "sec-plan", "sec-reading-room"} <= set(anchors)
    assert '[id^="sec-"] { scroll-margin-top: calc(var(--bar) + 16px); }' in page
    for anchor in anchors:
        assert f'href="#{anchor}"' in page, f"{anchor} has no rail link"
    # The scroll-spy is scoped to the content column, never the whole page.
    assert "document.querySelectorAll('#content [id^=\"sec-\"]')" in page


def test_local_find_is_present_and_wired_to_the_catalog_filter(page: str) -> None:
    assert 'id="find"' in page and 'id="cube-table"' in page
    assert "function filterCatalog()" in page
    assert "byId('find').addEventListener('input', filterCatalog);" in page
    assert "data-find=" in page


def test_every_composed_command_carries_both_shell_shapes_and_a_full_payload(
    tmp_path: Path, page: str
) -> None:
    for identifier in ("cmd-process", "cmd-audit", "cmd-bench"):
        block = page[page.index(f'id="{identifier}"') : page.index("</article>", page.index(f'id="{identifier}"'))]
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
