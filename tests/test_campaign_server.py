"""The served campaign front door answers, browses, and writes exactly once.

Every test here drives a real ``ThreadingHTTPServer`` bound to an ephemeral
loopback port over real HTTP, because the contract the page slice codes against
is the wire contract, not the Python one: the status codes, the JSON shapes and
the ``Cache-Control`` header are the interface.

The one write this server may do gets the most attention.  A campaign home is
created on a cold start and never, ever created twice -- the second POST must
adopt the file byte for byte, because the operator's hand edits live in it.
"""

from __future__ import annotations

import json
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pytest

from echelle_spectra import campaign_server

CATALOG_NAME = "all-years.json"


def _catalog(tmp_path: Path) -> Path:
    """A small merged index in the shape tests/test_reading_room_page.py uses."""

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
            }
        ],
    }
    path = tmp_path / CATALOG_NAME
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class _Client:
    """The wire, and nothing above it."""

    def __init__(self, server: campaign_server.CampaignServer) -> None:
        self.server = server
        self.base = server.url

    def get(self, route: str, **query: str) -> tuple[int, str, dict[str, str]]:
        target = f"{self.base}{route}"
        if query:
            target = f"{target}?{urlencode(query)}"
        with urlopen(target, timeout=10) as response:
            return response.status, response.read().decode("utf-8"), dict(response.headers)

    def json_get(self, route: str, **query: str) -> Any:
        return json.loads(self.get(route, **query)[1])

    def post(self, route: str, payload: dict[str, Any]) -> Any:
        request = Request(
            f"{self.base}{route}",
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urlopen(request, timeout=10) as response:
            return json.loads(response.read().decode("utf-8"))


@contextmanager
def _serve(home: Path | None = None) -> Iterator[_Client]:
    server = campaign_server.make_server(home=home, port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield _Client(server)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


@pytest.fixture
def client() -> Iterator[_Client]:
    with _serve(None) as served:
        yield served


def _error_of(exc: HTTPError) -> dict[str, Any]:
    return json.loads(exc.read().decode("utf-8"))


def _is_html(body: str) -> bool:
    return body.lstrip().startswith("<") and ("<h1" in body.lower() or "<html" in body.lower())


# ---------------------------------------------------------------------------
# Cold start: the setup state
# ---------------------------------------------------------------------------


def test_a_cold_start_serves_a_setup_page_and_says_so(client: _Client) -> None:
    status, body, headers = client.get("/")
    assert status == 200
    assert headers["Content-Type"].startswith("text/html")
    assert _is_html(body)

    assert client.json_get("/api/state") == {"state": "setup", "home": None}


def test_every_answer_forbids_the_browser_cache(client: _Client) -> None:
    """A served campaign page a browser may cache is a page that can lie."""

    assert client.get("/")[2]["Cache-Control"] == "no-store"
    assert client.get("/api/state")[2]["Cache-Control"] == "no-store"
    assert client.get("/api/browse")[2]["Cache-Control"] == "no-store"


def test_an_unknown_route_is_a_json_404(client: _Client) -> None:
    with pytest.raises(HTTPError) as raised:
        client.get("/api/process")
    assert raised.value.code == 404
    assert "no such endpoint" in _error_of(raised.value)["error"]


# ---------------------------------------------------------------------------
# Browsing the owner's folders, read-only
# ---------------------------------------------------------------------------


def test_browse_without_a_path_lists_the_real_drives(client: _Client) -> None:
    payload = client.json_get("/api/browse")
    assert payload["path"] is None
    assert payload["parent"] is None
    assert payload["dirs"] == []
    assert payload["drives"], "a machine with no drives cannot be browsed"
    assert all(isinstance(drive, str) and drive for drive in payload["drives"])


def test_browse_lists_sorted_folders_and_counts_their_sif_files(
    client: _Client, tmp_path: Path
) -> None:
    """The picker's whole point: which of these folders holds the frames?"""

    (tmp_path / "zulu").mkdir()
    (tmp_path / "alpha").mkdir()
    (tmp_path / "mike").mkdir()
    for name in ("a.sif", "b.SIF", "c.sif"):
        (tmp_path / "alpha" / name).write_bytes(b"")
    (tmp_path / "mike" / "one.SIF").write_bytes(b"")
    (tmp_path / "mike" / "notes.txt").write_bytes(b"")
    # A nested frame must not be counted: the count describes this folder only.
    (tmp_path / "zulu" / "deeper").mkdir()
    (tmp_path / "zulu" / "deeper" / "deep.sif").write_bytes(b"")
    # A loose file beside the folders is never listed; directories only.
    (tmp_path / "loose.sif").write_bytes(b"")

    payload = client.json_get("/api/browse", path=str(tmp_path))
    assert payload["drives"] is None
    assert payload["path"] == str(tmp_path.resolve())
    assert payload["parent"] == str(tmp_path.resolve().parent)
    assert [entry["name"] for entry in payload["dirs"]] == ["alpha", "mike", "zulu"]
    assert [entry["sif_count"] for entry in payload["dirs"]] == [3, 1, 0]
    assert Path(payload["dirs"][0]["path"]) == (tmp_path / "alpha").resolve()


def test_browse_refuses_a_path_that_is_not_a_folder(client: _Client, tmp_path: Path) -> None:
    bogus = tmp_path / "no-such-folder"
    with pytest.raises(HTTPError) as raised:
        client.get("/api/browse", path=str(bogus))
    assert raised.value.code == 400
    message = _error_of(raised.value)["error"]
    assert "no-such-folder" in message and "\n" not in message

    a_file = tmp_path / "frame.sif"
    a_file.write_bytes(b"")
    with pytest.raises(HTTPError) as raised:
        client.get("/api/browse", path=str(a_file))
    assert raised.value.code == 400


# ---------------------------------------------------------------------------
# The one write
# ---------------------------------------------------------------------------


def test_choosing_a_folder_writes_one_commented_campaign_toml(
    client: _Client, tmp_path: Path
) -> None:
    home = tmp_path / "campaign"
    home.mkdir()

    written = client.post("/api/home", {"folder": str(home)})
    target = home / "campaign.toml"
    assert Path(written["written"]) == target.resolve()

    text = target.read_text(encoding="utf-8")
    assert text.startswith("#")
    assert "wrote this file" in text
    assert "hand-editable" in text
    assert 'catalog = "all-years.json"' in text
    assert 'output = "campaign-page"' in text
    assert 'registry = "calibration_registry.toml"' in text
    assert 'calibrations = "calibrations"' in text
    # Nothing else was created; the server writes one file and no folders.
    assert sorted(item.name for item in home.iterdir()) == ["campaign.toml"]

    state = client.json_get("/api/state")
    assert state["home"] == str(target.resolve())
    # Written, but nothing catalogued yet: the honest word for that is empty.
    assert state["state"] == "empty"

    status, body, _ = client.get("/")
    assert status == 200
    assert _is_html(body)


def test_a_second_choice_adopts_the_file_and_never_rewrites_it(
    client: _Client, tmp_path: Path
) -> None:
    """The operator's hand edits are the reason this file is never overwritten."""

    home = tmp_path / "campaign"
    home.mkdir()
    client.post("/api/home", {"folder": str(home)})
    target = home / "campaign.toml"
    target.write_text('catalog = "mine.json"\n# edited by hand\n', encoding="utf-8")
    before = target.read_text(encoding="utf-8")
    stamp = target.stat().st_mtime_ns

    adopted = client.post("/api/home", {"folder": str(home)})
    assert Path(adopted["adopted"]) == target.resolve()
    assert "written" not in adopted
    assert target.read_text(encoding="utf-8") == before
    assert target.stat().st_mtime_ns == stamp


def test_choosing_a_folder_that_does_not_exist_is_refused_and_never_created(
    client: _Client, tmp_path: Path
) -> None:
    missing = tmp_path / "not-here"
    with pytest.raises(HTTPError) as raised:
        client.post("/api/home", {"folder": str(missing)})
    assert raised.value.code == 400
    assert "not-here" in _error_of(raised.value)["error"]
    assert not missing.exists(), "the server must never create the folder itself"

    assert client.json_get("/api/state")["state"] == "setup"


def test_a_home_post_without_a_folder_is_refused(client: _Client) -> None:
    with pytest.raises(HTTPError) as raised:
        client.post("/api/home", {})
    assert raised.value.code == 400
    assert "folder" in _error_of(raised.value)["error"]


# ---------------------------------------------------------------------------
# A real campaign: the page itself
# ---------------------------------------------------------------------------


def _ready_home(tmp_path: Path) -> Path:
    _catalog(tmp_path)
    home = tmp_path / "campaign.toml"
    home.write_text(
        f'catalog = "{CATALOG_NAME}"\noutput = "campaign-page"\n', encoding="utf-8"
    )
    return home


def test_a_real_campaign_serves_the_built_page_and_reports_ready(tmp_path: Path) -> None:
    home = _ready_home(tmp_path)
    with _serve(home) as client:
        assert client.json_get("/api/state") == {"state": "ready", "home": str(home.resolve())}
        status, body, headers = client.get("/")
        assert status == 200
        assert headers["Content-Type"].startswith("text/html")
        assert headers["Cache-Control"] == "no-store"
        assert "Echelle campaign" in body
        assert "NIFS-A" in body
        # The page is built, not merely streamed: the output folder now holds it.
        assert (tmp_path / "campaign-page" / "index.html").is_file()


def test_the_page_is_rebuilt_from_the_home_on_every_request(tmp_path: Path) -> None:
    """A page left open overnight must not be able to show yesterday's catalog."""

    home = _ready_home(tmp_path)
    with _serve(home) as client:
        assert "NIFS-A" in client.get("/")[1]
        catalog = json.loads((tmp_path / CATALOG_NAME).read_text(encoding="utf-8"))
        catalog["sources"][0]["volume_label"] = "NIFS-RENAMED"
        (tmp_path / CATALOG_NAME).write_text(json.dumps(catalog), encoding="utf-8")
        assert "NIFS-RENAMED" in client.get("/")[1]


def test_a_home_naming_an_absent_catalog_is_empty_not_broken(tmp_path: Path) -> None:
    home = tmp_path / "campaign.toml"
    home.write_text(f'catalog = "{CATALOG_NAME}"\noutput = "campaign-page"\n', encoding="utf-8")
    with _serve(home) as client:
        assert client.json_get("/api/state")["state"] == "empty"
        status, body, _ = client.get("/")
        assert status == 200
        assert _is_html(body)


# ---------------------------------------------------------------------------
# Resolution and the command surface
# ---------------------------------------------------------------------------


def test_home_resolution_accepts_a_folder_a_file_or_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "campaign.toml"
    home.write_text('catalog = "x.json"\n', encoding="utf-8")

    assert campaign_server.resolve_home(str(tmp_path)) == home.resolve()
    assert campaign_server.resolve_home(str(home)) == home.resolve()

    monkeypatch.chdir(tmp_path)
    assert campaign_server.resolve_home(None) == home.resolve()

    empty = tmp_path / "elsewhere"
    empty.mkdir()
    monkeypatch.chdir(empty)
    # A cold start is not an error; it is the state the setup page exists for.
    assert campaign_server.resolve_home(None) is None
    with pytest.raises(campaign_server.HomeError):
        campaign_server.resolve_home(str(empty))


def test_serve_binds_loopback_only_and_offers_no_host_flag() -> None:
    server = campaign_server.make_server(home=None, port=0)
    try:
        assert server.server_address[0] == "127.0.0.1"
        assert server.url == f"http://127.0.0.1:{server.port}"
    finally:
        server.server_close()

    parser = campaign_server._build_parser("echelle serve")
    flags = {flag for action in parser._actions for flag in action.option_strings}
    assert "--home" in flags and "--port" in flags
    assert "--host" not in flags


def test_serve_main_announces_one_real_url_and_stops_cleanly(capsys) -> None:
    server = campaign_server.make_server(home=None, port=0)
    port = server.port
    finished: list[int] = []
    runner = threading.Thread(target=lambda: finished.append(campaign_server.run_server(server)))
    runner.start()
    # Stop it only once it has genuinely answered, so shutdown cannot outrun the
    # loop it is meant to end.
    for attempt in range(200):
        try:
            with urlopen(f"http://127.0.0.1:{port}/api/state", timeout=5) as answer:
                answer.read()
            break
        except OSError:  # pragma: no cover - only while the socket is still coming up
            assert attempt < 199, "the server never answered"
            threading.Event().wait(0.02)
    server.shutdown()
    runner.join(timeout=10)
    assert finished == [0]
    assert f"serving the campaign page at http://127.0.0.1:{port}" in capsys.readouterr().out


def test_the_cli_dispatches_serve(monkeypatch: pytest.MonkeyPatch) -> None:
    from echelle_spectra import cli

    seen: dict[str, Any] = {}

    def _fake(argv: list[str], *, prog: str) -> int:
        seen["argv"] = argv
        seen["prog"] = prog
        return 0

    monkeypatch.setattr(campaign_server, "serve_main", _fake)
    assert cli.main(["serve", "--port", "0"]) == 0
    assert seen == {"argv": ["--port", "0"], "prog": "echelle serve"}
