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

import http.client
import json
import os
import subprocess
import sys
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pytest

from echelle_spectra import campaign_server, reading_room

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

    def raw(
        self,
        method: str,
        route: str,
        *,
        headers: dict[str, str],
        body: bytes | None = None,
    ) -> tuple[int, str]:
        """One request with the headers spelled by hand, Host included.

        The boundary being tested is made of headers a browser fills in and a
        client library also fills in, so the test writes them itself rather
        than accepting whatever ``urlopen`` would have chosen.
        """

        connection = http.client.HTTPConnection("127.0.0.1", self.server.port, timeout=10)
        try:
            connection.putrequest(method, route, skip_host=True, skip_accept_encoding=True)
            for name, value in headers.items():
                connection.putheader(name, value)
            connection.putheader("Content-Length", str(len(body or b"")))
            connection.endheaders(body or None)
            answer = connection.getresponse()
            return answer.status, answer.read().decode("utf-8")
        finally:
            connection.close()


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
        # Not scanned yet is a state of the FULL page, never a lesser one: the
        # operator begins work before any scan result exists — Browse, the
        # composer and the calibrate stepper are all already there.
        assert "picker-dialog" in body
        assert "Browse" in body
        assert "Compose" in body
        assert "Calibrate" in body
        # The scaffolding catalog never lands in the campaign home; the built
        # page does land in the CONFIGURED output, exactly as a ready build
        # would.
        assert not (tmp_path / CATALOG_NAME).exists()
        assert sorted(item.name for item in tmp_path.iterdir()) == [
            "campaign-page",
            "campaign.toml",
        ]
        assert (tmp_path / "campaign-page" / "index.html").is_file()


def test_the_zero_state_composes_against_the_home_and_never_the_stand_in(
    tmp_path: Path,
) -> None:
    """A fresh operator's very first copied command must land in the campaign.

    The page a not-yet-scanned campaign is served renders over a synthesized
    empty catalog, because zero drives is a state of the full page.  What that
    stand-in must never do is get composed INTO the commands: the first command
    the operator copies is the one that writes the central index, and writing it
    into the system temp leaves the campaign reading an empty file forever.
    """

    home = tmp_path / "campaign.toml"
    home.write_text(f'catalog = "{CATALOG_NAME}"\noutput = "campaign-page"\n', encoding="utf-8")
    with _serve(home) as client:
        body = client.get("/")[1]

    wanted = (tmp_path / CATALOG_NAME).resolve().as_posix()
    windows = wanted.replace("/", "\\")
    # The catalog-writing command names the home's own file, in both shells.
    assert f"--central-index &quot;{wanted}&quot;" in body
    assert f"--central-index &quot;{windows}&quot;" in body
    # So does the field it is composed from, and the plan's central_index.
    assert f'<input type="text" id="f-catalog" value="{wanted}">' in body
    assert f"central_index = &quot;{wanted}&quot;" in body
    # And the drift audit reads the same file it will have been written to.
    assert f"--catalog &quot;{wanted}&quot;" in body
    # Nothing anywhere on the page points at the stand-in.
    assert "empty-catalog.json" not in body
    assert "echelle-serve-" not in body


# ---------------------------------------------------------------------------
# The boundary: loopback is where it starts, not where it ends
# ---------------------------------------------------------------------------


def _post_bytes(payload: dict[str, Any]) -> bytes:
    return json.dumps(payload).encode("utf-8")


def test_a_post_from_a_foreign_origin_is_refused_and_writes_nothing(
    client: _Client, tmp_path: Path
) -> None:
    """Loopback stops a second machine; it does not stop a second tab.

    Any page on any site can post to a loopback port, and a simple request does
    it with no preflight to refuse.  Unguarded, that page could write a
    campaign.toml into any folder on this machine and repoint the server at it.
    """

    target = tmp_path / "campaign"
    target.mkdir()
    status, body = client.raw(
        "POST",
        "/api/home",
        headers={
            "Host": f"127.0.0.1:{client.server.port}",
            "Origin": "http://evil.example",
            "Content-Type": "application/json",
        },
        body=_post_bytes({"folder": str(target)}),
    )
    assert status == 403
    assert "127.0.0.1" in json.loads(body)["error"]
    assert list(target.iterdir()) == [], "a refused request must not have written"
    assert client.server.home is None, "a refused request must not have moved the home"


def test_a_post_addressed_to_a_foreign_host_is_refused(client: _Client, tmp_path: Path) -> None:
    """A name that merely resolves here is the shape DNS rebinding has."""

    target = tmp_path / "campaign"
    target.mkdir()
    status, _ = client.raw(
        "POST",
        "/api/home",
        headers={
            "Host": "campaign.evil.example",
            "Content-Type": "application/json",
        },
        body=_post_bytes({"folder": str(target)}),
    )
    assert status == 403
    assert list(target.iterdir()) == []
    assert client.server.home is None


def test_browsing_the_owner_folders_is_refused_from_a_foreign_host(client: _Client) -> None:
    """The read side is worth as much as the write side: it maps the machine."""

    status, _ = client.raw(
        "GET", "/api/browse", headers={"Host": "campaign.evil.example"}
    )
    assert status == 403
    status, _ = client.raw(
        "GET",
        "/api/browse",
        headers={"Host": f"127.0.0.1:{client.server.port}", "Origin": "null"},
    )
    assert status == 403


def test_a_simple_cross_site_body_shape_is_refused_before_it_is_read(
    client: _Client, tmp_path: Path
) -> None:
    """``text/plain`` is the content type that needs no preflight, so it is the
    one a cross-site form or fetch would arrive as."""

    target = tmp_path / "campaign"
    target.mkdir()
    status, body = client.raw(
        "POST",
        "/api/home",
        headers={
            "Host": f"127.0.0.1:{client.server.port}",
            "Content-Type": "text/plain;charset=UTF-8",
        },
        body=_post_bytes({"folder": str(target)}),
    )
    assert status == 400
    assert "application/json" in json.loads(body)["error"]
    assert list(target.iterdir()) == []


def test_the_ordinary_same_origin_flow_still_works(client: _Client, tmp_path: Path) -> None:
    """The guard must cost the operator nothing: localhost and 127.0.0.1, with
    the page's own Origin or with none at all, all still answer."""

    target = tmp_path / "campaign"
    target.mkdir()
    port = client.server.port
    status, body = client.raw(
        "POST",
        "/api/home",
        headers={
            "Host": f"127.0.0.1:{port}",
            "Origin": f"http://127.0.0.1:{port}",
            "Content-Type": "application/json",
        },
        body=_post_bytes({"folder": str(target)}),
    )
    assert status == 200
    assert Path(json.loads(body)["written"]) == (target / "campaign.toml").resolve()

    assert client.raw("GET", "/api/state", headers={"Host": f"localhost:{port}"})[0] == 200
    assert client.raw("GET", "/api/state", headers={"Host": "127.0.0.1"})[0] == 200
    assert client.get("/api/state")[0] == 200


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


# ---------------------------------------------------------------------------
# Running a campaign verb from the page
# ---------------------------------------------------------------------------
#
# The owner's requirement, in his words: "We can run stuff from WebUI, I want
# that. But I don't want to corrupt data or hang the process and loose the
# progress."  Every test below is one half of that sentence.  What is launched
# is the same CLI verb a terminal would run, derived here and never sent by the
# page; it is launched detached, so nothing this server does can take a night's
# progress with it; and what the page then SAYS about it is read back off the
# files the run leaves, never off anything this process remembers.

#: A stand-in for a campaign verb: it records the argument list it was actually
#: handed, waits for the test to let it finish, and prints on both sides of the
#: wait so the log has something to show.
_STUB = '''\
import json
import pathlib
import sys
import time

root = pathlib.Path(__file__).with_name("stub-out")
root.mkdir(exist_ok=True)
(root / "argv.json").write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
print("stub is running", flush=True)
gate = root / "go.txt"
for _ in range(2400):
    if gate.is_file():
        break
    time.sleep(0.05)
(root / "done.txt").write_text("done", encoding="utf-8")
print("stub has finished", flush=True)
'''


def _campaign(tmp_path: Path) -> Path:
    """A campaign home with a catalog in it, ready to launch from."""

    root = tmp_path / "campaign"
    root.mkdir()
    _catalog(root)
    home = root / "campaign.toml"
    home.write_text(f'catalog = "{CATALOG_NAME}"\noutput = "campaign-page"\n', encoding="utf-8")
    return home


def _shots(tmp_path: Path) -> Path:
    """A data folder that exists, which is all the endpoint asks of one."""

    folder = tmp_path / "shots"
    folder.mkdir()
    (folder / "one.SIF").write_bytes(b"")
    return folder


def _install_stub(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the launcher at the stub and return where the stub reports."""

    stub = tmp_path / "stub_verb.py"
    stub.write_text(_STUB, encoding="utf-8")
    monkeypatch.setattr(campaign_server, "LAUNCH_PREFIX", (sys.executable, str(stub)))
    return stub.with_name("stub-out")


def _wait_for(predicate: Any, seconds: float = 30.0) -> bool:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.05)
    return predicate()


def _expected_argv(home: Path, verb: str, folder: Path, epoch: str = "") -> tuple[list[str], str]:
    """Derive the argument list the page's own composer would, independently."""

    values = dict(campaign_server.home_values(home))
    composition = reading_room.campaign_composition(
        values["catalog"], **campaign_server.campaign_keywords(values)
    )
    composed = reading_room.with_answers(
        composition["values"],
        folder=str(folder.resolve()),
        epoch=epoch,
        evidence_name=composition["evidence_name"],
    )
    return reading_room.run_argv(verb, composed), reading_room.run_command(verb, composed)


def test_a_run_launches_the_derived_argument_list_and_nothing_else(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The whole safety property in one test: the page names a verb, and this
    machine composes the command from the campaign's own files."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)
    out.mkdir(exist_ok=True)
    (out / "go.txt").write_text("go", encoding="utf-8")

    with _serve(home) as client:
        answer = client.post("/api/run", {"verb": "sample", "folder": str(shots), "epoch": ""})
    marker = answer["launched"]

    argv, command = _expected_argv(home, "sample", shots)
    assert marker["argv"] == argv
    assert marker["command"] == command
    assert marker["verb"] == "sample"
    assert marker["args_sha256"] and marker["pid"] > 0 and marker["started_at"]

    # The child really was handed that list -- not a shell string, not a copy.
    assert _wait_for(lambda: (out / "argv.json").is_file())
    assert json.loads((out / "argv.json").read_text(encoding="utf-8")) == argv

    # The marker and the log sit in the campaign's runs area, beside the
    # receipts, and the marker on disk is the one that was answered with.
    launches = campaign_server.launches_root(home)
    assert launches == home.parent / "local" / "runs" / "launches"
    on_disk = campaign_server.read_marker(campaign_server.marker_path(home, "sample", shots))
    assert on_disk == marker
    assert Path(marker["log"]).parent == launches
    assert _wait_for(lambda: "stub is running" in Path(marker["log"]).read_text(encoding="utf-8"))


def test_the_run_endpoint_reads_a_verb_name_and_never_a_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An endpoint that reads a flag list runs whatever another tab decided."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    _install_stub(tmp_path, monkeypatch)

    with _serve(home) as client:
        for payload, expected in (
            ({"verb": "frobnicate", "folder": str(shots)}, "is not a verb this page can run"),
            ({"verb": "", "folder": str(shots)}, "is not a verb this page can run"),
            ({"verb": ["sample"], "folder": str(shots)}, "is not a verb this page can run"),
            (
                {"verb": "sample", "folder": str(shots), "args": ["--force"]},
                "never args",
            ),
            (
                {"verb": "sample", "folder": str(shots), "command": "del /s /q C:\\"},
                "never command",
            ),
            ({"verb": "sample"}, "no data folder was named"),
            ({"verb": "sample", "folder": str(tmp_path / "gone")}, "not a readable folder"),
            (
                {"verb": "sample", "folder": str(shots), "epoch": "never-saved"},
                "knows no calibration named never-saved",
            ),
        ):
            with pytest.raises(HTTPError) as raised:
                client.post("/api/run", payload)
            assert raised.value.code == 400
            message = _error_of(raised.value)["error"]
            assert expected in message and "\n" not in message
        # The refused verb list is stated once, and it is the allowlist itself.
        with pytest.raises(HTTPError) as raised:
            client.post("/api/run", {"verb": "nope", "folder": str(shots)})
        assert "audit, bulk, sample" in _error_of(raised.value)["error"]

    assert not campaign_server.launches_root(home).exists(), "a refused request launched something"


def test_running_needs_a_campaign_home_to_launch_into(
    client: _Client, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    shots = _shots(tmp_path)
    _install_stub(tmp_path, monkeypatch)
    with pytest.raises(HTTPError) as raised:
        client.post("/api/run", {"verb": "sample", "folder": str(shots)})
    assert raised.value.code == 400
    assert "no campaign home yet" in _error_of(raised.value)["error"]


def test_a_run_from_a_foreign_origin_is_refused_and_starts_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The same boundary the write side has: loopback stops a second machine,
    the Host and Origin check stops a second tab."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    _install_stub(tmp_path, monkeypatch)

    with _serve(home) as client:
        status, body = client.raw(
            "POST",
            "/api/run",
            headers={
                "Host": f"127.0.0.1:{client.server.port}",
                "Origin": "http://evil.example",
                "Content-Type": "application/json",
            },
            body=_post_bytes({"verb": "sample", "folder": str(shots)}),
        )
        assert status == 403
        assert "127.0.0.1" in json.loads(body)["error"]
        status, _ = client.raw(
            "POST",
            "/api/run",
            headers={"Host": "campaign.evil.example", "Content-Type": "application/json"},
            body=_post_bytes({"verb": "sample", "folder": str(shots)}),
        )
        assert status == 403
    assert not campaign_server.launches_root(home).exists()


def test_two_clicks_do_not_double_launch_and_the_refusal_names_the_running_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)

    with _serve(home) as client:
        first = client.post("/api/run", {"verb": "sample", "folder": str(shots)})["launched"]
        assert _wait_for(lambda: (out / "argv.json").is_file())

        with pytest.raises(HTTPError) as raised:
            client.post("/api/run", {"verb": "sample", "folder": str(shots)})
        assert raised.value.code == 409
        line = _error_of(raised.value)["error"]
        assert line.startswith("sample is already running on")
        assert f"pid {first['pid']}" in line and "\n" not in line

        # A different verb on the same folder is a different run, not a clash.
        second = client.post("/api/run", {"verb": "audit", "folder": str(shots)})["launched"]
        assert second["pid"] != first["pid"]

        # Once it is gone, the same verb may run again.
        (out / "go.txt").write_text("go", encoding="utf-8")
        assert _wait_for(lambda: not campaign_server.pid_alive(int(first["pid"])))
        again = client.post("/api/run", {"verb": "sample", "folder": str(shots)})["launched"]
        assert again["pid"] != first["pid"]


def test_a_launched_run_outlives_the_server_that_started_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The owner's actual fear, tested: the page and the server going away must
    not take the run -- or the night's progress -- with them."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)

    with _serve(home) as client:
        marker = client.post("/api/run", {"verb": "sample", "folder": str(shots)})["launched"]
        assert _wait_for(lambda: (out / "argv.json").is_file())
        assert campaign_server.pid_alive(int(marker["pid"]))
    # _serve has now shut the server down and closed its socket; the child was
    # still waiting when it did, and is released only afterwards.
    assert campaign_server.pid_alive(int(marker["pid"])), "the server took its run with it"
    (out / "go.txt").write_text("go", encoding="utf-8")
    assert _wait_for(lambda: (out / "done.txt").is_file()), "the run never finished"


def test_a_fresh_server_reads_the_run_state_off_the_files_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No state is held here, so a restarted server -- or a second tab, or a run
    started in a terminal -- reads exactly the same."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)

    with _serve(home) as client:
        marker = client.post("/api/run", {"verb": "sample", "folder": str(shots)})["launched"]
        assert _wait_for(lambda: (out / "argv.json").is_file())

    with _serve(home) as client:
        page = client.get("/")[1]
        assert f"pid {marker['pid']}" in page
        assert 'data-run="sample" disabled' in page
    (out / "go.txt").write_text("go", encoding="utf-8")
    assert _wait_for(lambda: not campaign_server.pid_alive(int(marker["pid"])))


def test_the_card_states_running_finished_and_failed_from_the_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Four states, none inferred from a neighbour: a process that is gone
    without a receipt is failed and shows its log, never a quiet finished."""

    from echelle_spectra.campaign_run import RunReceipt

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)

    idle = campaign_server.run_states(home, plan=home.parent / "campaign-plan.toml")
    assert idle["sample"]["state"] == "idle"
    # The bulk step states the one fact that decides whether it can run at all.
    assert "campaign-plan.toml is not saved" in idle["bulk"]["line"]

    with _serve(home) as client:
        marker = client.post("/api/run", {"verb": "sample", "folder": str(shots)})["launched"]
        assert _wait_for(lambda: (out / "argv.json").is_file())
        running = campaign_server.run_states(home)["sample"]
        assert running["state"] == "running"
        assert str(shots) in running["line"] and f"pid {marker['pid']}" in running["line"]

        (out / "go.txt").write_text("go", encoding="utf-8")
        assert _wait_for(lambda: not campaign_server.pid_alive(int(marker["pid"])))

        # Gone, and it left no receipt: that is failed, with the log to read.
        failed = campaign_server.run_states(home)["sample"]
        assert failed["state"] == "failed"
        assert "left no receipt" in failed["line"]
        assert "stub is running" in failed["tail"]

        # The receipt the real verb would have written turns it into finished.
        receipt = RunReceipt.create(
            campaign_server.runs_root(home) / "2026-08-19_00-00-00-shots",
            source_root=shots,
            output_root=shots,
            pattern="*.SIF",
            volume_label="shots",
            snapshot_id="20250926_cmos",
            expected_files=1,
        )
        receipt.finish("completed")
        finished = campaign_server.run_states(home)["sample"]
        assert finished["state"] == "finished"
        assert "receipt 2026-08-19_00-00-00-shots [completed]" in finished["line"]
        assert finished["tail"] == ""


def test_an_audit_reports_the_immutable_evidence_it_named(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An audit's durable receipt is the verdict file, so that is what is read."""

    home = _campaign(tmp_path)
    shots = _shots(tmp_path)
    out = _install_stub(tmp_path, monkeypatch)
    out.mkdir(exist_ok=True)
    (out / "go.txt").write_text("go", encoding="utf-8")

    with _serve(home) as client:
        marker = client.post("/api/run", {"verb": "audit", "folder": str(shots)})["launched"]
    assert _wait_for(lambda: not campaign_server.pid_alive(int(marker["pid"])))

    evidence = Path(marker["evidence"])
    assert evidence.name == "drift-evidence-001.json"
    assert campaign_server.run_states(home)["audit"]["state"] == "failed"
    evidence.write_text(json.dumps({"verdict": "aligned"}), encoding="utf-8")
    state = campaign_server.run_states(home)["audit"]
    assert state["state"] == "finished"
    assert "verdict aligned in drift-evidence-001.json" in state["line"]


def test_the_launch_is_detached_with_no_console_and_nothing_to_answer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: dict[str, Any] = {}

    class _Recorded:
        def __init__(self, argv: list[str], **keywords: Any) -> None:
            seen["argv"] = argv
            seen["keywords"] = keywords
            self.pid = 4242

    monkeypatch.setattr(campaign_server.subprocess, "Popen", _Recorded)
    log = tmp_path / "local" / "runs" / "launches" / "sample-1234.log"
    assert campaign_server.spawn_detached(["a", "b"], cwd=tmp_path, log=log) == 4242

    keywords = seen["keywords"]
    assert seen["argv"] == ["a", "b"]
    assert keywords["cwd"] == str(tmp_path)
    # A batch run that stops to ask a question nobody can answer is a hang.
    assert keywords["stdin"] == subprocess.DEVNULL
    assert keywords["stderr"] == subprocess.STDOUT
    if os.name == "nt":
        flags = keywords["creationflags"]
        assert flags & subprocess.DETACHED_PROCESS
        assert flags & subprocess.CREATE_NEW_PROCESS_GROUP
    else:  # pragma: no cover - this campaign runs on Windows
        assert keywords["start_new_session"] is True
    # The log is appended to, and says what was launched into it.
    assert "launching: a b" in log.read_text(encoding="utf-8")


def test_asking_whether_a_run_is_alive_never_ends_it() -> None:
    """On Windows ``os.kill`` terminates rather than probes; a liveness check
    that can end the run it asks about would be the exact failure to avoid."""

    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        for _ in range(3):
            assert campaign_server.pid_alive(child.pid)
        assert child.poll() is None, "asking cost the run its life"
    finally:
        child.terminate()
        child.wait(timeout=30)
    assert not campaign_server.pid_alive(child.pid)
    assert not campaign_server.pid_alive(0) and not campaign_server.pid_alive(-1)


def test_the_page_offers_no_way_to_stop_a_run(tmp_path: Path) -> None:
    """Launching is the whole grant: a page that can end a terabyte conversion
    halfway is a page that can lose a night's progress."""

    home = _campaign(tmp_path)
    with _serve(home) as client:
        with pytest.raises(HTTPError) as raised:
            client.post("/api/stop", {"pid": 1})
        assert raised.value.code == 404
    assert not hasattr(campaign_server, "stop_run")
    assert not hasattr(campaign_server, "kill_run")


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
