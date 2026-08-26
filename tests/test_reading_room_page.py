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
from echelle_spectra.reading_room import (
    _SOURCE_NOTES,
    SNAPSHOT_SCAN_BUDGET,
    SNAPSHOT_SCAN_DEPTH,
    _derived_from_folder,
    _saved_snapshots,
    build_reading_room,
    render_markdown,
    with_answers,
)
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


def test_saved_snapshots_are_seen_without_a_registry(tmp_path: Path) -> None:
    """A snapshot the bench wrote is real before any registry names it.

    The owner's field case, 2026-08-18: a NAS calibrations root holding years
    of snapshots, no registry yet — and the page said "no snapshot is in
    reach", which was false.
    """

    bare = tmp_path / "bare.json"
    bare.write_text(
        json.dumps(
            {
                "schema": "echelle-merged-catalog/v1",
                "generated_at": "2026-08-18T00:00:00.000+00:00",
                "sources": [],
            }
        ),
        encoding="utf-8",
    )
    root = tmp_path / "cal"
    # The owner's real shape: the snapshot sits a couple of levels down.
    nested = root / "calibrations" / "20190314_cmos"
    nested.mkdir(parents=True)
    (nested / "snapshot.toml").write_text('id = "20190314_cmos"\n', encoding="utf-8")

    page = build_reading_room(
        bare, tmp_path / "web-saved", calibrations_root=root
    ).read_text(encoding="utf-8")

    assert "20190314_cmos — saved, not in any registry" in page
    assert "saved snapshot folder(s) on the calibrations root" in page
    # The registry epoch step still honestly wants a registry.
    assert "--registry" in page


def _bare_catalog(tmp_path: Path, name: str = "bare.json") -> Path:
    """A merged index naming no drive: the shape a fresh campaign has."""

    path = tmp_path / name
    path.write_text(
        json.dumps(
            {
                "schema": "echelle-merged-catalog/v1",
                "generated_at": "2026-08-18T00:00:00.000+00:00",
                "sources": [],
            }
        ),
        encoding="utf-8",
    )
    return path


def _snapshot_folder(folder: Path, declared: str | None) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    body = f'id = "{declared}"\n' if declared is not None else "created_utc = \"\"\n"
    (folder / "snapshot.toml").write_text(body, encoding="utf-8")
    return folder


def test_a_saved_snapshot_is_offered_under_the_id_its_manifest_declares(
    tmp_path: Path,
) -> None:
    """The manifest id is the identity a registry and a cube both name.

    A folder renamed on the way to a NAS still binds the same snapshot, so
    offering the folder name would offer a calibration no registry can resolve.
    Two folders binding one id are one choice, not two.
    """

    root = tmp_path / "cal"
    _snapshot_folder(root / "2019-march-run", "20190314_cmos")
    _snapshot_folder(root / "copied-from-the-nas", "20190314_cmos")
    # The one case the folder name is the truth: a binder that declares no id.
    _snapshot_folder(root / "20200101_cmos", None)

    records, truncated = _saved_snapshots(root)
    assert [record["id"] for record in records] == ["20190314_cmos", "20200101_cmos"]
    assert truncated is False
    # Which of two folders binding one id is named must not ride on directory
    # order, so the earlier path wins and the same root always reads the same.
    assert records[0]["root"].endswith("/2019-march-run")

    page = build_reading_room(
        _bare_catalog(tmp_path), tmp_path / "web-ids", calibrations_root=root
    ).read_text(encoding="utf-8")
    assert page.count('<option value="20190314_cmos">') == 1
    assert '<option value="20200101_cmos">' in page
    # A folder name is a place, never an identity: it may be read as where the
    # data is, and it is never offered as a calibration to choose.
    chooser = re.search(r'<select id="f-epoch">(.*?)</select>', page, re.S).group(1)
    for folder in ("2019-march-run", "copied-from-the-nas"):
        assert folder not in chooser, f"the folder name {folder!r} is offered as an identity"
    assert 'class="path-line">' in page
    assert "copied-from-the-nas" not in page, "the second folder binding one id is not a snapshot"


def test_the_snapshot_walk_stops_at_its_budget_and_the_page_says_so(
    tmp_path: Path,
) -> None:
    """The walk runs on the served page's request thread, over what is usually a
    NAS share.  Depth alone bounds nothing: three levels of a wide root is still
    an unbounded number of network stats, so breadth is bounded too — and a list
    cut short must say it was, or it is the same lie as claiming none exist."""

    root = tmp_path / "cal"
    root.mkdir()
    for index in range(SNAPSHOT_SCAN_BUDGET + 40):
        _snapshot_folder(root / f"{index:04d}_cmos", f"{index:04d}_cmos")

    ids, truncated = _saved_snapshots(root)
    assert truncated is True
    assert len(ids) <= SNAPSHOT_SCAN_BUDGET
    assert ids, "the budget returns what it reached, never nothing"

    page = build_reading_room(
        _bare_catalog(tmp_path), tmp_path / "web-wide", calibrations_root=root
    ).read_text(encoding="utf-8")
    assert "this list is what the scan reached and not all of it" in page
    assert "A snapshot missing here is not missing from the root." in page

    # A root the budget covers says nothing of the kind.
    narrow = tmp_path / "narrow"
    _snapshot_folder(narrow / "20190314_cmos", "20190314_cmos")
    quiet = build_reading_room(
        _bare_catalog(tmp_path, "bare-2.json"), tmp_path / "web-narrow", calibrations_root=narrow
    ).read_text(encoding="utf-8")
    assert "what the scan reached" not in quiet


def test_the_composer_never_seeds_unassigned_as_a_calibration(tmp_path: Path) -> None:
    """"unassigned" is a cube's way of saying it had no calibration identity.

    The select has never offered it; seeding it into the composed plan named it
    anyway, as the calibration the registry supposedly selects for these shots.
    """

    payload = {
        "schema": "echelle-merged-catalog/v1",
        "generated_at": "2026-08-18T00:00:00.000+00:00",
        "sources": [
            {
                "drive_id": "id-a",
                "volume_label": "NIFS-A",
                "drive_root": (tmp_path / "drive-a").as_posix(),
                "catalog_path": "echelle-catalog.json",
                "run": None,
                "cubes": [{"path": "a.nc", "shot_number": "1", "snapshot_id": "unassigned"}],
            }
        ],
    }
    catalog = tmp_path / "unassigned.json"
    catalog.write_text(json.dumps(payload), encoding="utf-8")

    page = build_reading_room(catalog, tmp_path / "web-unassigned").read_text(encoding="utf-8")
    data = json.loads(re.search(r"^const DATA=(.*);$", page, re.M).group(1))
    assert data["values"]["epoch"] == ""
    assert "that is unassigned." not in html.unescape(page)
    # The composer's own select does not offer it either.  The catalog FILTER
    # still does, and rightly: a cube that carries no identity is something a
    # reader looks for.
    chooser = re.search(r'<select id="f-epoch">(.*?)</select>', page, re.S).group(1)
    assert "unassigned" not in chooser


def test_the_two_halves_of_the_folder_derivation_agree_at_a_bare_root() -> None:
    """``derivedFrom`` in the page's JavaScript is the other half of
    ``_derived_from_folder``, and the two must answer the same for the same
    folder.  A bare root is where they used to disagree: JavaScript stopped one
    character early and kept the slash Python's ``rstrip`` removes, deriving
    ``//drift-evidence-001.json`` against Python's ``/drift-evidence-001.json``.

    No JavaScript runs here, so the JavaScript half is pinned as the text it
    is; the Python half is pinned by calling it.
    """

    from echelle_spectra import reading_room

    assert _derived_from_folder("/", "drift-evidence-001.json")["verdict"] == (
        "/drift-evidence-001.json"
    )
    source = Path(reading_room.__file__).read_text(encoding="utf-8")
    assert "while (trimmed.length > 0 && trimmed.charAt(trimmed.length - 1) === '/')" in source
    assert "trimmed.length > 1" not in source


#: The campaign home a read-only data drive sends its written things into.
_HOME = "/Users/owner/Echelle-campaigns/HD-LXU3-DATA"


def test_a_drive_that_refuses_writes_derives_into_the_campaign_home() -> None:
    """The one rule the redirect adds, stated where the rule lives.

    Cubes and evidence belong on the drive beside their own data — until that
    drive is an NTFS volume a Mac mounted read-only, where writing them there
    is not a preference but an impossibility.  The data folder is untouched:
    it is still what is READ, and only the written things move.
    """

    derived = _derived_from_folder(
        "/Volumes/HD-LXU3/DATA", "drift-evidence-002.json", redirect=_HOME
    )
    assert derived == {
        "input": "/Volumes/HD-LXU3/DATA",
        "output": f"{_HOME}/cubes/DATA",
        "cubes": f"{_HOME}/cubes/DATA",
        # The label is the data folder's own name either way: it names the
        # shots, not the place their cubes happened to land.
        "label": "DATA",
        "verdict": f"{_HOME}/drift-evidence-002.json",
    }
    # Without the redirect the same folder derives exactly where it always did.
    beside = _derived_from_folder("/Volumes/HD-LXU3/DATA", "drift-evidence-002.json")
    assert beside["output"] == "/Volumes/HD-LXU3/DATA"
    assert beside["verdict"] == "/Volumes/HD-LXU3/DATA/drift-evidence-002.json"
    assert beside["label"] == derived["label"]


def test_an_unanswered_data_folder_redirects_nowhere_at_all() -> None:
    """The unfilled marker is a placeholder, not a folder.

    Redirecting it would put the campaign home into the fields of a page whose
    reader has answered nothing yet — a derived path presented as an answer,
    with the marker still sitting in the field it came from.
    """

    derived = _derived_from_folder("   ", "drift-evidence-001.json", redirect=_HOME)
    assert derived["input"] == "<data folder>"
    assert derived["output"] == "<data folder>"
    assert derived["cubes"] == "<data folder>"
    assert _HOME not in json.dumps(derived)


def test_the_composed_run_carries_the_redirect_into_every_written_value() -> None:
    """``with_answers`` is what the served Run control composes through, so the
    launched command must land where the page's own fields say it will."""

    composed = with_answers(
        {"pattern": "*.SIF", "plan": "campaign-plan.toml", "epoch": "stale"},
        folder="/Volumes/HD-LXU3/DATA",
        epoch="20250926_cmos",
        evidence_name="drift-evidence-002.json",
        redirect=_HOME,
    )
    assert composed["input"] == "/Volumes/HD-LXU3/DATA"
    assert composed["output"] == f"{_HOME}/cubes/DATA"
    # The cubes a step reads are the cubes the previous step wrote.
    assert composed["cubes"] == composed["output"]
    assert composed["verdict"] == f"{_HOME}/drift-evidence-002.json"
    # Everything not derived from the folder is left exactly as it was handed.
    assert composed["pattern"] == "*.SIF" and composed["plan"] == "campaign-plan.toml"
    assert composed["epoch"] == "20250926_cmos"
    # No redirect, no move: the same two answers derive beside the data.
    beside = with_answers(
        {"pattern": "*.SIF"},
        folder="/Volumes/HD-LXU3/DATA",
        epoch="",
        evidence_name="drift-evidence-002.json",
    )
    assert beside["output"] == "/Volumes/HD-LXU3/DATA"


def test_the_two_halves_of_the_redirect_agree_on_a_read_only_drive() -> None:
    """``derivedFrom`` grew the redirect branch that ``_derived_from_folder``
    grew, and the two must answer the same for the same folder and home.

    No JavaScript runs here, so the JavaScript half is pinned as the text it
    is; the Python half is pinned by calling it, and the two shapes are spelled
    side by side so an edit to one without the other reads as the mismatch it
    would be in the browser.
    """

    from echelle_spectra import reading_room

    label, evidence = "DATA", "drift-evidence-002.json"
    derived = _derived_from_folder(f"/Volumes/HD-LXU3/{label}", evidence, redirect=_HOME)
    assert derived["output"] == _HOME + "/cubes/" + label
    assert derived["cubes"] == _HOME + "/cubes/" + label
    assert derived["verdict"] == _HOME + "/" + evidence

    source = Path(reading_room.__file__).read_text(encoding="utf-8")
    # The branch is taken on exactly the same three conditions Python's is: a
    # folder that was actually answered, a measured refusal, and a home.
    assert "if (folder && folderState.readonly && DATA.home) {" in source
    assert "var folderState = { readonly: DATA.readonly === true };" in source
    assert "var base = folderPath(DATA.home);" in source
    # And it spells the same three values.
    assert "output: base + '/cubes/' + label," in source
    assert "cubes: base + '/cubes/' + label," in source
    assert "verdict: base + '/' + DATA.evidence_name" in source
    # The note the page shows follows the same condition, never a fourth rule.
    assert "note.hidden = !(folder && folderState.readonly && DATA.home);" in source


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
    # No Compose button: every field recomposes live on input, and a button
    # whose work already happened is a control that visibly does nothing.
    assert 'id="compose"' not in now_left
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
    # The count is the CLI's to derive from what the folder holds, not the page's.
    assert "--sample auto" in row


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
    assert "Point at the shots and the calibration; the rest is derived." in page
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
# The two-input front door
# ---------------------------------------------------------------------------


def _composer(text: str) -> str:
    return _element(text, 'id="f-input"', "article")


def _numbered_evidence(tmp_path: Path, number: str, cube: str = "sample-a.nc") -> Path:
    path = tmp_path / f"drift-evidence-{number}.json"
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


def test_the_composer_asks_two_questions_and_folds_every_derived_value(page: str) -> None:
    card = _composer(page)
    head, fold = card[: card.index("<details")], card[card.index("<details") :]
    # Two visible controls: where the shots are, and which calibration to use.
    assert head.count("<input") == 1 and head.count("<select") == 1
    assert 'id="f-input"' in head and 'id="f-epoch"' in head
    # Everything else is derived, and the fold that holds it opens closed.
    tag = re.search(r'<details[^>]*id="composer-advanced"[^>]*>', card)
    assert tag is not None and " open" not in tag.group()
    assert "Advanced — every derived value, editable" in fold
    for field in (
        "f-output", "f-label", "f-verdict", "f-registry",
        "f-calibrations", "f-catalog", "f-plan", "f-pattern",
    ):
        assert f'id="{field}"' in fold
    # A card, not a card in a card, and its head is pressable.
    assert '<article class="card' not in fold
    assert ".fold-group > summary { cursor: pointer;" in page
    # The fields this front door removed are gone from the page entirely.
    for gone in ("f-drive", "f-every", "f-bench", "f-sample"):
        assert f'id="{gone}"' not in page


def test_no_example_drive_ships_with_the_page(page: str) -> None:
    for literal in ("T:\\", "D:\\", "E:\\", "F:\\", "T:/", "D:/", "E:/", "F:/"):
        assert literal not in page, f"a baked example drive {literal!r} shipped"
    field = re.search(r'<input type="text" id="f-input"[^>]*>', page)
    assert field is not None
    # Nothing on this machine records the data folder, so the page fills nothing
    # in: a placeholder asks, a value would claim.
    assert "value=" not in field.group()
    assert "placeholder=" in field.group()
    assert 'value="shots"' not in page
    # Until it is answered, every composed line says so in one place.
    assert 'id="unfilled-note"' in page
    assert page.count("Paste the data folder in the rail first") == 1


def test_the_composed_sequence_reads_in_campaign_order(page: str) -> None:
    order = [
        page.index(f'id="{identifier}"')
        for identifier in ("cmd-sample", "cmd-audit", "cmd-process", "cmd-bench")
    ]
    assert order == sorted(order), "sample, then check, then bulk; the bench is the aside"
    sample = html.unescape(_element(page, 'id="cmd-sample"', "article"))
    assert "--sample auto" in sample and "--sample 20" not in sample
    assert "--central-index" in sample
    audit = html.unescape(_element(page, 'id="cmd-audit"', "article"))
    assert "echelle drift audit" in audit
    assert "--every" not in audit, "the CLI derives the interval now"
    # The gate is stated once, on the step it actually gates.
    assert "refused until that evidence exists" in _element(page, 'id="cmd-process"', "article")


def test_the_composed_audit_writes_the_next_free_evidence_name(tmp_path: Path) -> None:
    given = _numbered_evidence(tmp_path, "001")
    text = build_reading_room(
        _one_drive_catalog(tmp_path, [_connected(tmp_path, "NIFS-A")]),
        tmp_path / "web",
        drift_paths=[given],
    ).read_text(encoding="utf-8")
    audit = _element(text, 'id="cmd-audit"', "article")
    payload = html.unescape(re.findall(r'data-copy="([^"]*)"', audit)[1])
    assert payload.endswith('drift-evidence-002.json"')
    # Evidence is immutable: the composed command never names the file it read.
    assert "drift-evidence-001.json" not in payload
    assert "drift-evidence-002.json" in _composer(text)
    # A build that was handed none starts at one.
    fresh = build_reading_room(
        _one_drive_catalog(tmp_path, [_connected(tmp_path, "NIFS-A")]), tmp_path / "web2"
    ).read_text(encoding="utf-8")
    assert "drift-evidence-001.json" in _element(fresh, 'id="cmd-audit"', "article")


def test_the_alignment_verdict_leads_every_drive_card(page: str) -> None:
    leading = re.findall(r'<p class="chips"><span class="pill ([^"]*)">([^<]*)<', page)
    assert leading, "no chip row was rendered"
    assert all("verdict-lead" in classes for classes, _ in leading)
    cards = _element(page, 'id="sec-drives-cards"')
    assert '<span class="pill verdict-shifted verdict-lead">shifted</span>' in cards
    # A drive nothing measured stays unmeasured; it is never dressed as aligned.
    assert 'verdict-lead">alignment unmeasured<' in cards
    assert "aligned</span>" not in cards
    assert ".pill.verdict-lead {" in page
    # The Now tab's drive rows lead with the same fact.
    assert 'verdict-lead">shifted<' in _drive_rows(page)[0]


# ---------------------------------------------------------------------------
# The illustrated how-to
# ---------------------------------------------------------------------------
#
# "why do we need to make a separate explainer? That should be the default UX...
# Or at least the campaign's reading room should have this. Kind how-to with
# illustrations. Not the textbook." (queezz, 2026-08-19.)


_HOWTO_ANCHORS = (
    "sec-howto",
    "sec-howto-eyes",
    "sec-howto-flow",
    "sec-howto-letters",
    "sec-howto-checklist",
)


def _howto(page: str) -> str:
    return _element(page, 'id="sec-howto"')


def test_the_how_to_leads_the_reading_room_with_its_own_anchors(page: str) -> None:
    reading = _view(page, "reading")
    # It leads: a person needs the pictures before the canon, not after it.
    assert reading.index('id="sec-howto"') < reading.index('id="sec-reading-room"')
    for anchor in _HOWTO_ANCHORS:
        assert f'id="{anchor}"' in reading
    # Registered where the room indexes its sections, in the room's own rail
    # pattern -- the page's general pin already demands a link for every sec-
    # anchor; this one demands the order too.
    rail = _rail_group(page, "rail-right", "reading")
    positions = [rail.index(f'href="#{anchor}"') for anchor in _HOWTO_ANCHORS]
    assert positions == sorted(positions)
    assert rail.index('href="#sec-howto"') < rail.index('href="#sec-reading-room"')
    # One bordered region with headed groups inside it, never a card in a card.
    howto = _howto(page)
    assert howto.startswith('<section class="panel howto" id="sec-howto">')
    assert "<article" not in howto
    assert re.findall(r'<h3 id="(sec-howto-[a-z]+)"', howto) == [
        "sec-howto-eyes",
        "sec-howto-flow",
        "sec-howto-letters",
        "sec-howto-checklist",
    ]


def test_the_how_to_teaches_with_figures_rather_than_paragraphs(page: str) -> None:
    howto = _howto(page)
    figures = re.findall(r"<figure class=\"figure\">(.*?)</figure>", howto, re.S)
    assert len(figures) == 3
    for figure in figures:
        assert '<svg viewBox="' in figure, "a figure that does not scale by its own viewBox"
        assert 'role="img"' in figure and 'aria-label="' in figure
        assert "<figcaption>" in figure
        # Hand-drawn primitives on the page's own ink, and the one accent group
        # spent through the page's own variable.
        assert 'stroke="currentColor"' in figure
        assert 'class="accent"' in figure
        assert "<rect" in figure and "<text" in figure and "<path" in figure
        # Wide content scrolls inside its own box: a drawing squeezed to a
        # phone's width would render its labels at a few pixels, and the page
        # body must never gain a sideways scrollbar to avoid that.
        assert figure.startswith('<div class="scroll-x">')
    assert ".figure svg { display: block; width: 100%; min-width: 40rem;" in page
    # The three drawings the owner asked for, named by what they show.
    labels = re.findall(r'aria-label="([^"]+)"', howto)
    assert len(labels) == 3
    assert "echelle status looks from the folder you stand in" in labels[0]
    assert "the drift verdict gates what follows" in labels[1]
    assert "the same folder from every machine" in labels[2]
    assert ".figure .accent { color: var(--accent); }" in page
    # The drawing prints the scan's own numbers, never hand-typed glyphs: a
    # picture is the last place anyone would look for a drifted constant.
    assert f"({SNAPSHOT_SCAN_DEPTH} levels, {SNAPSHOT_SCAN_BUDGET} entries)" in howto
    assert "SCAN_DEPTH" not in howto and "SCAN_BUDGET" not in howto
    # Kind and short: the pictures carry the mechanism, so no part of the
    # how-to turns back into the textbook hum.
    assert "mermaid" not in page.lower()
    assert "<img" not in page and "base64" not in page


def test_the_how_to_checklist_answers_the_calibrations_question(page: str) -> None:
    checklist = _element(_howto(page), 'id="sec-howto-checklist"', "h3")
    howto = _howto(page)
    items = re.findall(r"<li>(.*?)</li>", howto[howto.index('id="sec-howto-checklist"') :], re.S)
    assert len(items) == 4
    joined = html.unescape(" ".join(items))
    # The checklist names the campaign file by its one real name.
    assert "the campaign file (campaign.toml)" in joined
    assert "echelle-campaign.toml" not in joined
    assert "calibrations line" in joined
    assert f"{SNAPSHOT_SCAN_DEPTH} levels down" in joined
    assert f"stops after {SNAPSHOT_SCAN_BUDGET} entries" in joined
    assert "a folder holding snapshot.toml" in joined
    assert "echelle status" in joined
    assert checklist  # the heading itself is the jump target the rail names


def test_the_how_to_bakes_in_no_real_path_only_the_placeholder_convention(
    page: str,
) -> None:
    """Placeholders only: no hostname, no address, no drive letter of anyone's.

    The page's own no-baked-drives pin already refuses a literal drive; this
    one refuses the rest of a real machine and pins what stands in for them.
    """

    howto = _howto(page)
    # A drive letter is the subject of one figure and is still never a literal:
    # the page's own placeholder convention carries it.
    assert "&lt;letter&gt;:\\NIFS" in howto
    assert not re.search(r"[A-Za-z]:\\", html.unescape(howto)), "a real drive letter got in"
    assert r"\\NAS\share\NIFS" in howto
    # No host, no address, no user profile.
    assert not re.search(r"\b\d{1,3}(\.\d{1,3}){3}\b", howto), "an IP address got in"
    for forbidden in ("queezz", "Dropbox", "workdata", "10.10.10"):
        assert forbidden not in howto


def test_inline_doc_code_may_wrap_so_a_long_path_never_widens_the_page(page: str) -> None:
    """A packaged doc's inline code -- a long resource path with no space in it --
    must be allowed to wrap, or it forces a horizontal page scroll at narrow
    widths (measured at 760px). The rule lives on .doc-body code, the same
    overflow-wrap idiom .scan-line already uses for long UNC roots."""

    match = re.search(r"\.doc-body code\s*\{([^}]*)\}", page)
    assert match, ".doc-body code rule is missing"
    assert "overflow-wrap" in match.group(1)


def test_the_reading_room_links_the_documentation_site_once(page: str) -> None:
    link = '<a href="https://queezz.github.io/echelle_spectra">'
    assert page.count(link) == 1
    assert link in _view(page, "reading")
    assert page.count("queezz.github.io/echelle_spectra") == 2, "one anchor, one label"


# ---------------------------------------------------------------------------
# Surviving F7 contracts
# ---------------------------------------------------------------------------


def test_every_composed_command_carries_both_shell_shapes_and_a_full_payload(
    tmp_path: Path, page: str
) -> None:
    for identifier in ("cmd-sample", "cmd-audit", "cmd-process", "cmd-bench"):
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
    # Each shape writes the same path in its own shell's spelling.
    assert "<data folder>\\drift-evidence" in windows
    assert "<data folder>/drift-evidence" in posix
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
        "<script src",
        "<link ",
        "@import",
    ):
        assert forbidden not in page, f"the page must not contain {forbidden!r}"
    # The page still fetches nothing of its own. One outbound address is allowed,
    # and only as a link a reader may choose to follow: never as a request the
    # page makes for itself.
    for match in re.finditer("https://", page):
        prefix = page[max(0, match.start() - len('<a href="')) : match.start()]
        assert prefix == '<a href="', "https:// may appear only as a link target"
    assert "color-scheme: light dark;" in page
    assert "@media (prefers-color-scheme: dark)" in page
    # Every control keeps a resting border and fill.
    assert "button, select, input, textarea {" in page


# ---------------------------------------------------------------------------
# The calibration record — the owner's own reading order
# ---------------------------------------------------------------------------
#
# "Looking at the calibration tab, I'm lost and I don't know what that thing is
# telling me. I need 1. Place, where is the data? 2. What lamps, how about
# spheres? 3. Ok, then show me the data from the calibration bench. Image
# shifted delta x delta y? Theta? Absolute calibration, the snapshot of the
# factors?"  (queezz, 2026-08-18.)  These tests read the built page in exactly
# that order, against the shape his own 20190314_cmos binder has.

_DIGEST = "0" * 64

#: The owner's real binder, trimmed to the fields the tab renders — and one
#: honest hole: Xe is named among the lamps and no frame records it.
_BINDER = f"""\
schema = "echelle-snapshot/v1"
id = "20190314_cmos"
created_utc = "2019-03-14T08:12:00+00:00"
detector = "cmos"
lamps = ["Hg", "Ne", "ThAr", "Xe"]
base_snapshot = "20190207_cmos"

[validity]
date_from = "2019-03-14"

[alignment]
dx_px = 0.42
dy_px = -1.1
rotation_deg = 0.0131
rms_px = 0.031
wavelength_correction_applied = true
pattern_correction_applied = false
pattern_band_offset_note = "the sphere's bands sit above the chosen pattern"
pattern_band_offset_rows = 3.5
pattern_band_offset_orders = 12
vetted_set = "BH paper"
vetted_lineage = ["wavelength.txt", "wavelength-2019.txt"]
science_validation = "measured"
science_lines_validated = 6
science_residual_rms_nm = 0.0082

[qc]
lines_used = 37
worst_residual_px = 0.19
sphere_comparison = "ready"
sphere_comparison_reference = "previous_sphere.sif + previous_sphere_bg.sif"

[[artifacts]]
role = "pattern"
path = "pattern.txt"
sha256 = "{_DIGEST}"
size_bytes = 2048
source_name = "pattern-2019.txt"

[[artifacts]]
role = "wavelength"
path = "wavelength.txt"
sha256 = "{_DIGEST}"
size_bytes = 4096
source_name = "wavelength-2019.txt"

[[artifacts]]
role = "integral"
path = "integral.txt"
sha256 = "{_DIGEST}"
size_bytes = 1024
source_name = "integral-2019.txt"

[[artifacts]]
role = "sphere"
kind = "referenced"
path = "../../sphere-0.1s-x3.sif"
sha256 = "{_DIGEST}"
size_bytes = 398458880
source_name = "sphere-0.1s-x3.sif"

[[artifacts]]
role = "sphere_background"
kind = "referenced"
path = "../../sphere-bg-0.1s-x3.sif"
sha256 = "{_DIGEST}"
size_bytes = 398458880
source_name = "sphere-bg-0.1s-x3.sif"

[[artifacts]]
role = "lamp"
kind = "referenced"
label = "Hg"
path = "../../Hg-30s.sif"
sha256 = "{_DIGEST}"
size_bytes = 1048576
source_name = "Hg-30s.sif"

[[artifacts]]
role = "lamp"
kind = "referenced"
label = "Ne"
path = "../../Ne-30s.sif"
sha256 = "{_DIGEST}"
size_bytes = 1048576
source_name = "Ne-30s.sif"

[[artifacts]]
role = "lamp"
kind = "referenced"
label = "ThAr"
path = "../../ThAr-120s.sif"
sha256 = "{_DIGEST}"
size_bytes = 1048576
source_name = "ThAr-120s.sif"
"""


def _bench_root(tmp_path: Path, binder: str = _BINDER, name: str = "20190314_cmos") -> Path:
    root = tmp_path / "cal"
    folder = root / "calibrations" / name
    folder.mkdir(parents=True)
    (folder / "snapshot.toml").write_text(binder, encoding="utf-8")
    return root


def _record_page(tmp_path: Path, root: Path, out: str = "web-record") -> str:
    return build_reading_room(
        _bare_catalog(tmp_path, f"{out}.json"), tmp_path / out, calibrations_root=root
    ).read_text(encoding="utf-8")


def _snapshot_card(text: str) -> str:
    return _element(text, 'id="sec-cal-snap-20190314-cmos"', "article")


def test_the_calibration_tab_leads_with_the_snapshots_and_reads_in_his_order(
    tmp_path: Path,
) -> None:
    page = _record_page(tmp_path, _bench_root(tmp_path))
    view = _view(page, "calibration")
    # The record leads; the pipeline context follows it, never the other way.
    order = [
        view.index(f'id="{anchor}"')
        for anchor in ("sec-cal-snapshots", "sec-cal-epochs", "sec-drift")
    ]
    assert order == sorted(order)
    card = _snapshot_card(view)
    # Place, then sources, then the bench — his three questions, in his order.
    assert re.findall(r"<h4>([^<]+)</h4>", card) == [
        "Where",
        "Lamps and spheres",
        "Rigid alignment",
        "Wavelength fit",
        "Absolute calibration",
    ]
    # The rail names every snapshot, so a long root is navigable rather than
    # only scrollable.
    rail = _rail_group(page, "rail-right", "calibration")
    assert (
        '<a href="#sec-cal-snap-20190314-cmos" class="sectnav-link">20190314_cmos</a>' in rail
    )
    assert rail.index("sec-cal-snapshots") < rail.index("sec-cal-epochs") < rail.index("sec-drift")
    # And every jump target this tab defines is reachable from that index, in a
    # namespace the rest of the page still does not collide with.
    ids = re.findall(r'id="([^"]+)"', page)
    assert len(ids) == len(set(ids)), "duplicate element ids"
    for anchor in re.findall(r'id="(sec-cal-[^"]+)"', view):
        assert f'href="#{anchor}"' in rail, f"{anchor} has no rail link"


def test_a_snapshot_says_where_it_is(tmp_path: Path) -> None:
    root = _bench_root(tmp_path)
    page = _record_page(tmp_path, root)
    folder = (root / "calibrations" / "20190314_cmos").as_posix()
    assert f'<p class="path-line">{_escaped(folder)}</p>' in _snapshot_card(page)


def test_a_snapshot_names_every_lamp_and_both_spheres(tmp_path: Path) -> None:
    card = _snapshot_card(_record_page(tmp_path, _bench_root(tmp_path)))
    sources = card[card.index("Lamps and spheres") : card.index("Rigid alignment")]
    for species, frame in (("Hg", "Hg-30s.sif"), ("Ne", "Ne-30s.sif"), ("ThAr", "ThAr-120s.sif")):
        assert f"<td>{species}</td>" in sources
        assert f"<td>{frame}</td>" in sources
    # The lamp named in the array with no frame behind it is stated, not dropped
    # and not invented.
    xenon = sources[sources.index("<td>Xe</td>") :]
    assert xenon[: xenon.index("</tr>")].count("<td>not recorded</td>") == 3
    # Both halves of the sphere pair, by the names that carry their exposure.
    assert "<td>sphere-0.1s-x3.sif</td>" in sources
    assert "<td>sphere-bg-0.1s-x3.sif</td>" in sources
    assert "<td>referenced</td>" in sources
    # The binder records no exposure field; that is taught once for the surface.
    assert _record_page(tmp_path, _bench_root(tmp_path / "twin"), out="web-once").count(
        "A binder records no exposure time"
    ) == 1


def test_the_bench_numbers_are_read_off_the_binder(tmp_path: Path) -> None:
    card = _snapshot_card(_record_page(tmp_path, _bench_root(tmp_path)))
    for name, value in (
        ("Δx", "0.42 px"),
        ("Δy", "-1.1 px"),
        ("θ", "0.0131 deg"),
        ("RMS", "0.031 px"),
        ("wavelength table shifted", "yes"),
        ("order pattern shifted", "no"),
        ("lines used", "37"),
        ("worst residual", "0.19 px"),
        ("vetted set", "BH paper"),
        ("vetted lineage", "wavelength.txt → wavelength-2019.txt"),
        ("sphere comparison", "ready"),
        ("compared against", "previous_sphere.sif + previous_sphere_bg.sif"),
        ("factor table", "integral-2019.txt"),
    ):
        pair = (
            f'<span class="facts-name">{_escaped(name)}</span>'
            f'<span class="facts-value">{_escaped(value)}</span>'
        )
        assert pair in card, f"{name} does not read {value!r}"


def test_a_thin_binder_says_not_recorded_rather_than_inventing(tmp_path: Path) -> None:
    """Historical snapshots predate most of these fields, and say so."""

    thin = 'id = "20190314_cmos"\ndetector = "cmos"\nlamps = ["Hg"]\n'
    card = _snapshot_card(_record_page(tmp_path, _bench_root(tmp_path, thin), out="web-thin"))
    for absent in ("Δx", "RMS", "lines used", "worst residual", "vetted set", "sphere comparison"):
        pair = (
            f'<span class="facts-name">{_escaped(absent)}</span>'
            '<span class="facts-value">not recorded</span>'
        )
        assert pair in card, f"{absent} was invented rather than named absent"
    # Nothing is dressed as a measurement it never made.
    assert "0.0" not in card and "0 px" not in card


#: The three sentences an empty snapshot scan may say, verbatim.  They are the
#: contract: the owner spent a day hunting a snapshot the page could see all
#: along because the select said only that it had none, never where it looked.
_UNSET_LINE = "This campaign names no calibrations root, so no folder was scanned for snapshots."


def _missing_line(root: Path) -> str:
    return f"No snapshots under {root.as_posix()} — that folder does not exist on this machine."


def _empty_line(root: Path) -> str:
    return f"{root.as_posix()} holds no snapshot within {SNAPSHOT_SCAN_DEPTH} levels."


def _select_region(page: str) -> str:
    """The composer's calibration select and whatever stands under it."""

    card = _composer(page)
    return card[card.index('<select id="f-epoch"') : card.index("<details")]


def test_no_snapshot_root_and_an_empty_one_are_different_facts(tmp_path: Path) -> None:
    """Three conditions a root can be in, and three different sentences.

    The middle one is the whole point: a root that is simply not on this
    machine used to read exactly like a root that answered and held nothing,
    which is what a mapped drive letter from another box looks like.
    """

    nothing = build_reading_room(
        _bare_catalog(tmp_path, "no-root.json"), tmp_path / "web-no-root"
    ).read_text(encoding="utf-8")
    assert _UNSET_LINE in nothing
    assert "holds no snapshot within" not in nothing
    assert "does not exist on this machine" not in nothing

    empty = tmp_path / "empty-root"
    empty.mkdir()
    swept = build_reading_room(
        _bare_catalog(tmp_path, "empty.json"), tmp_path / "web-empty", calibrations_root=empty
    ).read_text(encoding="utf-8")
    assert _escaped(_empty_line(empty)) in swept
    assert _UNSET_LINE not in swept
    assert "does not exist on this machine" not in swept

    gone = tmp_path / "mapped-elsewhere"
    absent = build_reading_room(
        _bare_catalog(tmp_path, "gone.json"), tmp_path / "web-gone", calibrations_root=gone
    ).read_text(encoding="utf-8")
    assert _escaped(_missing_line(gone)) in absent
    assert _UNSET_LINE not in absent
    assert "holds no snapshot within" not in absent


def test_the_empty_select_names_the_root_it_scanned_and_its_condition(
    tmp_path: Path,
) -> None:
    """The line lives at the select it emptied, not on some other tab.

    "No snapshots" is the answer that cost a day; "no snapshots under THIS
    folder, which is not on this machine" is the answer that ends it.
    """

    gone = tmp_path / "mapped-elsewhere"
    page = build_reading_room(
        _bare_catalog(tmp_path, "beside.json"), tmp_path / "web-beside", calibrations_root=gone
    ).read_text(encoding="utf-8")
    region = _select_region(page)
    assert _escaped(_missing_line(gone)) in region
    # One line, in a state class the page already owns -- no new hue for it.
    assert region.count('class="note scan-line') == 1
    assert 'class="note scan-line state-missing-drive"' in region
    # And the same fact, in the same words, is what the Calibration tab says
    # where it would otherwise list snapshots: one wording, two rooms, and a
    # reader never meets both at once.
    assert page.count(_escaped(_missing_line(gone))) == 2
    assert _escaped(_missing_line(gone)) in _view(page, "calibration")


def test_the_scan_line_is_silent_once_the_scan_finds_something(tmp_path: Path) -> None:
    root = tmp_path / "cal"
    _snapshot_folder(root / "20190314_cmos", "20190314_cmos")
    page = build_reading_room(
        _bare_catalog(tmp_path, "found.json"), tmp_path / "web-found", calibrations_root=root
    ).read_text(encoding="utf-8")
    assert 'class="note scan-line' not in page
    for silent in (_UNSET_LINE, "does not exist on this machine", "holds no snapshot within"):
        assert silent not in page


def test_the_scan_line_says_the_depth_the_walk_actually_uses(tmp_path: Path) -> None:
    """The sentence counts levels, so it reads the constant the walk spends."""

    from echelle_spectra import reading_room

    assert reading_room._saved_snapshots.__kwdefaults__["depth"] == SNAPSHOT_SCAN_DEPTH
    empty = tmp_path / "empty-root"
    empty.mkdir()
    page = build_reading_room(
        _bare_catalog(tmp_path, "depth.json"), tmp_path / "web-depth", calibrations_root=empty
    ).read_text(encoding="utf-8")
    assert f"within {SNAPSHOT_SCAN_DEPTH} levels" in page


def test_an_unreadable_binder_is_its_own_state(tmp_path: Path) -> None:
    root = tmp_path / "cal"
    folder = root / "20190314_cmos"
    folder.mkdir(parents=True)
    (folder / "snapshot.toml").write_text("id = [broken\n", encoding="utf-8")
    page = _record_page(tmp_path, root, out="web-broken")
    assert "This binder could not be read" in page
    assert "Rigid alignment" not in page, "an unreadable binder claims no measurement"


# ---------------------------------------------------------------------------
# What discovery skipped, and the plan sentence that has to stand alone
# ---------------------------------------------------------------------------


def test_a_run_that_pruned_folders_says_which_ones(tmp_path: Path) -> None:
    """The page was the one surface where a pruned run still read complete."""

    source = _connected(
        tmp_path,
        "NIFS-A",
        run={
            "id": "run-a",
            "state": "completed",
            "counts": {"exported": 1},
            "gate": "verdict",
            "pruned_dirs": ["20190207", "calibrations"],
        },
        cubes=[{"path": "a.nc", "shot_number": "1", "snapshot_id": "20260812_cmos"}],
    )
    quiet = _connected(tmp_path, "NIFS-B", run={"id": "run-b", "state": "completed", "counts": {}})
    page = build_reading_room(
        _one_drive_catalog(tmp_path, [source, quiet]), tmp_path / "web-pruned"
    ).read_text(encoding="utf-8")
    line = '<p class="pruned">This run skipped 2 calibration folder(s): 20190207, calibrations</p>'
    # Both surfaces a drive is read on carry it: the Drives card and the Now row.
    assert line in _element(page, 'id="sec-drives-cards"')
    assert line in _drive_rows(page)[0]
    # A run that pruned nothing says nothing; silence is not a claim here.
    assert "skipped 0 calibration folder(s)" not in page
    assert page.count('class="pruned"') == 2


def test_the_plan_sentence_stands_with_and_without_a_calibration(tmp_path: Path) -> None:
    """It used to render "for these shots that is ." with nothing chosen."""

    def _plan(text: str) -> str:
        return html.unescape(
            re.search(r'id="plan-out"[^>]*>(.*?)</textarea>', text, flags=re.S).group(1)
        )

    bare = build_reading_room(
        _bare_catalog(tmp_path, "no-epoch.json"), tmp_path / "web-no-epoch"
    ).read_text(encoding="utf-8")
    assert "# The registry selects the epoch from each shot's own date.\n" in _plan(bare)
    assert "that is ." not in _plan(bare)

    named = build_reading_room(
        _one_drive_catalog(
            tmp_path,
            [
                _connected(
                    tmp_path,
                    "NIFS-A",
                    cubes=[{"path": "a.nc", "shot_number": "1", "snapshot_id": "20260812_cmos"}],
                )
            ],
        ),
        tmp_path / "web-epoch",
    ).read_text(encoding="utf-8")
    assert (
        "# The registry selects the epoch from each shot's own date; for these shots "
        "that is 20260812_cmos.\n"
    ) in _plan(named)
    # The page's own JavaScript is the other half of the same rule.
    assert "function epochClause(epoch)" in named
    assert "values.epoch_clause = epochClause(values.epoch);" in named


def test_the_fifth_verdict_word_is_its_own_loud_state(page: str) -> None:
    # era-misassigned-calibration arrived with the 2026-08-19 auditor round:
    # the lines fit the other isotope at one era shift and the calendar
    # excludes it.  It must render as its own alarm, never as unrecognized.
    from echelle_spectra.reading_room import _verdict_state

    css_class, label = _verdict_state("era-misassigned-calibration")
    assert css_class == "verdict-era-misassigned"
    assert label == "era-misassigned calibration"
    assert ".verdict-era-misassigned { color" in page
    assert ".card.verdict-era-misassigned { border-left" in page
    assert "era-misassigned-calibration" in page
