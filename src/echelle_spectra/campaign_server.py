"""The served campaign front door: ``echelle serve``.

The page is served, not merely written, so that it can browse the owner's real
folders and adopt a campaign home on a cold start.  Everything else about it
stays exactly as strict as the static page it serves:

* it binds ``127.0.0.1`` only, deliberately -- there is no ``--host`` flag,
  because nothing here is meant to answer a second machine;
* it runs no processing, ever.  No endpoint converts a SIF, builds a cube, or
  shells out.  The terabyte conversion is the operator's own command;
* it writes exactly one thing, ever: a ``campaign.toml`` into a folder the
  operator picked, and never over one that already exists.

Every other route reads.  ``GET /`` rebuilds the campaign page fresh from the
campaign home on each request, so a page left open overnight cannot show
yesterday's catalog, and every response carries ``Cache-Control: no-store`` so
the browser cannot show it either.
"""

from __future__ import annotations

import argparse
import json
import os
import string
import sys
import tempfile
import threading
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlsplit

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

HOST = "127.0.0.1"
DEFAULT_PORT = 48917
CAMPAIGN_FILENAME = "campaign.toml"

STATE_SETUP = "setup"
STATE_EMPTY = "empty"
STATE_READY = "ready"

#: The file the page may write, and the only one.  Conventional names, spelled
#: out with the reason each is there, because the operator edits this by hand.
CAMPAIGN_TEMPLATE = """\
# The echelle campaign page wrote this file when you picked this folder as the
# campaign home.  It is an ordinary hand-editable TOML file: rename anything
# here and the page follows it the next time you reload.
#
# Every path is read relative to this file.

# The merged cube index 'echelle catalog' writes.
catalog = "all-years.json"
# Where the built campaign page is written.
output = "campaign-page"
# The ordered epoch registry that says which calibration covers which shots.
registry = "calibration_registry.toml"
# The snapshot root the registry names its epochs from.
calibrations = "calibrations"
"""

SETUP_PLACEHOLDER = "<h1>setup</h1>"
EMPTY_PLACEHOLDER = "<h1>empty campaign</h1>"


class HomeError(ValueError):
    """One plain sentence about a campaign home that cannot be opened."""


# ---------------------------------------------------------------------------
# Campaign home resolution
# ---------------------------------------------------------------------------


def resolve_home(value: str | os.PathLike[str] | None) -> Path | None:
    """Return the ``campaign.toml`` a ``--home`` argument names, or ``None``.

    Without an argument the current folder is consulted; a campaign home found
    there is adopted, and its absence is not an error -- it is the cold start
    the setup page exists for.
    """

    if value is None:
        candidate = Path.cwd() / CAMPAIGN_FILENAME
        return candidate.resolve() if candidate.is_file() else None
    path = Path(value).expanduser()
    if path.is_dir():
        candidate = path / CAMPAIGN_FILENAME
        if not candidate.is_file():
            raise HomeError(f"no {CAMPAIGN_FILENAME} in {path.resolve()}")
        return candidate.resolve()
    if path.is_file():
        return path.resolve()
    raise HomeError(f"no campaign home at {path.resolve()}")


def _relative_to(home: Path, value: Any) -> Path | None:
    """Read one campaign.toml path value, relative to the file that named it."""

    if not isinstance(value, str) or not value.strip():
        return None
    candidate = Path(value.strip()).expanduser()
    if not candidate.is_absolute():
        candidate = home.parent / candidate
    return Path(os.path.normpath(candidate))


def home_values(home: Path) -> dict[str, Any]:
    """Resolve one campaign.toml into the absolute paths the page is built from.

    Unreadable TOML is not fatal here: the home is still the home, and the page
    that reports an empty campaign is a better answer than a traceback.
    """

    try:
        with home.open("rb") as handle:
            document = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError):
        document = {}
    raw_drift = document.get("drift")
    if isinstance(raw_drift, str):
        raw_drift = [raw_drift]
    elif not isinstance(raw_drift, list):
        raw_drift = []
    drift = [resolved for resolved in (_relative_to(home, item) for item in raw_drift) if resolved]
    return {
        "home": home,
        "catalog": _relative_to(home, document.get("catalog", "all-years.json")),
        "output": _relative_to(home, document.get("output", "campaign-page")),
        "registry": _relative_to(home, document.get("registry", "calibration_registry.toml")),
        "calibrations": _relative_to(home, document.get("calibrations", "calibrations")),
        "drift": drift,
    }


def home_state(home: Path | None) -> str:
    """Say which of the three states the server is in right now."""

    if home is None:
        return STATE_SETUP
    catalog = home_values(home)["catalog"]
    return STATE_READY if catalog is not None and catalog.is_file() else STATE_EMPTY


# ---------------------------------------------------------------------------
# Reading the owner's folders
# ---------------------------------------------------------------------------


def list_drives() -> list[str]:
    """Every drive root this machine actually has, newest stdlib first."""

    lister = getattr(os, "listdrives", None)
    if lister is not None:
        try:
            found = [str(drive) for drive in lister()]
        except OSError:
            found = []
        if found:
            return found
    probed = [f"{letter}:\\" for letter in string.ascii_uppercase if os.path.exists(f"{letter}:\\")]
    return probed or [os.path.abspath(os.sep)]


def count_sif(folder: Path) -> int:
    """Count ``*.sif`` files sitting directly in one folder, never recursing.

    Counted case-insensitively but only once per file: a case-folding filesystem
    answers both ``*.sif`` and ``*.SIF`` with the same entries, and a folder of
    forty frames reported as eighty would be a lie in the picker.
    """

    total = 0
    try:
        with os.scandir(folder) as entries:
            for entry in entries:
                if not entry.name.lower().endswith(".sif"):
                    continue
                try:
                    if entry.is_file():
                        total += 1
                except OSError:  # pragma: no cover - vanished mid-walk
                    continue
    except (PermissionError, OSError):
        return 0
    return total


def browse(raw: str | None) -> dict[str, Any]:
    """Answer one read-only picker step: the drives, or one folder's children."""

    if not raw or not raw.strip():
        return {"path": None, "parent": None, "drives": list_drives(), "dirs": []}
    path = Path(raw.strip()).expanduser()
    if not path.is_dir():
        raise HomeError(f"not a readable folder: {path}")
    resolved = path.resolve()
    dirs = []
    try:
        with os.scandir(resolved) as entries:
            for entry in entries:
                try:
                    if not entry.is_dir():
                        continue
                except OSError:  # pragma: no cover - vanished mid-walk
                    continue
                child = Path(entry.path)
                try:
                    has_snapshot = (child / "snapshot.toml").is_file()
                except OSError:  # pragma: no cover - unreadable mount points
                    has_snapshot = False
                dirs.append(
                    {
                        "name": entry.name,
                        "path": str(child),
                        "sif_count": count_sif(child),
                        "has_snapshot": has_snapshot,
                    }
                )
    except (PermissionError, OSError) as exc:
        raise HomeError(f"cannot read {resolved}: {exc.strerror or exc}") from None
    parent = resolved.parent
    return {
        "path": str(resolved),
        "parent": None if parent == resolved else str(parent),
        "dirs": sorted(dirs, key=lambda item: (item["name"].casefold(), item["name"])),
        "drives": None,
    }


def adopt_or_write_home(folder: str | os.PathLike[str]) -> tuple[str, Path]:
    """The one write this server may do, and the guard that keeps it one.

    Returns ``("adopted", path)`` when a campaign home was already there, or
    ``("written", path)`` when this call created it.  The folder itself is never
    created and an existing ``campaign.toml`` is never overwritten.
    """

    root = Path(folder).expanduser()
    if not root.is_dir():
        raise HomeError(f"no such folder: {root}")
    resolved = root.resolve()
    target = resolved / CAMPAIGN_FILENAME
    if target.is_file():
        return "adopted", target
    try:
        # Exclusive creation, so two clicks in two tabs cannot race one file away.
        with open(target, "x", encoding="utf-8", newline="\n") as handle:
            handle.write(CAMPAIGN_TEMPLATE)
    except FileExistsError:
        return "adopted", target
    except OSError as exc:
        raise HomeError(f"cannot write {target}: {exc.strerror or exc}") from None
    return "written", target


# ---------------------------------------------------------------------------
# Page rendering
# ---------------------------------------------------------------------------


def _reading_room_callable(name: str) -> Any:
    """Fetch one reading_room renderer lazily, tolerating its absence."""

    from . import reading_room

    return getattr(reading_room, name, None)


def render_setup(server: CampaignServer) -> str:
    renderer = _reading_room_callable("render_setup_page")
    if renderer is None:
        server.setup_page_rendered = False
        return SETUP_PLACEHOLDER
    server.setup_page_rendered = True
    return str(renderer())


def render_empty(server: CampaignServer, values: dict[str, Any]) -> str:
    renderer = _reading_room_callable("render_empty_campaign_page")
    if renderer is None:
        server.empty_page_rendered = False
        return EMPTY_PLACEHOLDER
    server.empty_page_rendered = True
    return str(renderer(values))


def empty_catalog_path(server: CampaignServer) -> Path:
    """One synthesized zero-sources catalog per server, in the system temp.

    A campaign that has not been scanned yet still gets the FULL page: zero
    drives is a state of the campaign, not a lesser page.  The scaffolding
    catalog lives in the system temp and is never written into the campaign
    home -- the real one appears there when the first processing run writes
    it, and this one is forgotten.
    """

    cached = getattr(server, "empty_catalog", None)
    if cached is not None and cached.is_file():
        return cached
    root = Path(tempfile.mkdtemp(prefix="echelle-serve-"))
    path = root / "empty-catalog.json"
    payload = {
        "schema": "echelle-merged-catalog/v1",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": [],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    server.empty_catalog = path
    return path


def build_page(server: CampaignServer, values: dict[str, Any]) -> bytes:
    """Rebuild the campaign page from the home and return the bytes it wrote."""

    from .reading_room import build_reading_room

    registry = values["registry"]
    calibrations = values["calibrations"]
    keywords: dict[str, Any] = {
        "drift_paths": [path for path in values["drift"] if path.is_file()],
        "document_paths": [],
        # A conventional name the operator has not created yet is not a broken
        # registry; it is simply not supplied.
        "registry_path": registry if registry is not None and registry.is_file() else None,
        "calibrations_root": (
            calibrations if calibrations is not None and calibrations.is_dir() else None
        ),
    }
    with server.build_lock:
        try:
            index = build_reading_room(
                values["catalog"], values["output"], served=True, **keywords
            )
            server.served_keyword = True
        except TypeError:
            index = build_reading_room(values["catalog"], values["output"], **keywords)
            server.served_keyword = False
        return Path(index).read_bytes()


# ---------------------------------------------------------------------------
# The server
# ---------------------------------------------------------------------------


class CampaignServer(ThreadingHTTPServer):
    """A loopback server holding one mutable fact: which campaign home is open."""

    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, address: tuple[str, int], handler: Any, *, home: Path | None) -> None:
        self._home = home
        self._home_lock = threading.Lock()
        self.build_lock = threading.Lock()
        #: Whether ``build_reading_room`` accepted ``served=True`` (None until built).
        self.served_keyword: bool | None = None
        #: Whether reading_room supplied the real setup/empty pages (None until served).
        self.setup_page_rendered: bool | None = None
        self.empty_page_rendered: bool | None = None
        super().__init__(address, handler)

    @property
    def port(self) -> int:
        return int(self.server_address[1])

    @property
    def url(self) -> str:
        return f"http://{HOST}:{self.port}"

    @property
    def home(self) -> Path | None:
        with self._home_lock:
            return self._home

    @home.setter
    def home(self, value: Path | None) -> None:
        with self._home_lock:
            self._home = value


class CampaignHandler(BaseHTTPRequestHandler):
    """Five routes, four of which only read."""

    server_version = "echelle-campaign"
    sys_version = ""
    protocol_version = "HTTP/1.1"

    # The console keeps its one line; per-request chatter would bury it.
    def log_message(self, fmt: str, *args: Any) -> None:
        return

    # -- responses ----------------------------------------------------------

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        # The house freshness contract: a served campaign page is never cached,
        # so what the operator sees is what the folders say now.
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, payload: dict[str, Any], status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self._send(status, body, "application/json")

    def _error(self, message: str, status: int = 400) -> None:
        self._json({"error": message}, status=status)

    def _html(self, markup: str) -> None:
        self._send(200, markup.encode("utf-8"), "text/html; charset=utf-8")

    # -- routing ------------------------------------------------------------

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's name
        parts = urlsplit(self.path)
        route = unquote(parts.path)
        if route == "/":
            self._get_page()
        elif route == "/api/state":
            self._get_state()
        elif route == "/api/browse":
            self._get_browse(parse_qs(parts.query))
        else:
            self._error(f"no such endpoint: {route}", status=404)

    def do_POST(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's name
        route = unquote(urlsplit(self.path).path)
        if route == "/api/home":
            self._post_home()
        else:
            self._error(f"no such endpoint: {route}", status=404)

    # -- endpoints ----------------------------------------------------------

    def _get_page(self) -> None:
        server: CampaignServer = self.server  # type: ignore[assignment]
        home = server.home
        if home is None:
            self._html(render_setup(server))
            return
        values = home_values(home)
        catalog = values["catalog"]
        if catalog is None or not catalog.is_file():
            # Not scanned yet: serve the full page over a synthesized empty
            # catalog, so every control -- Browse, the composer, the calibrate
            # stepper -- exists before the first scan result does.
            values = dict(values)
            values["catalog"] = empty_catalog_path(server)
            if values["output"] is None or not values["output"].parent.is_dir():
                values["output"] = empty_catalog_path(server).parent / "page"
        try:
            self._send(200, build_page(server, values), "text/html; charset=utf-8")
        except (OSError, ValueError) as exc:
            self._error(f"cannot build the campaign page: {exc}", status=500)

    def _get_state(self) -> None:
        server: CampaignServer = self.server  # type: ignore[assignment]
        home = server.home
        self._json({"state": home_state(home), "home": None if home is None else str(home)})

    def _get_browse(self, query: dict[str, list[str]]) -> None:
        values = query.get("path") or []
        try:
            self._json(browse(values[0] if values else None))
        except HomeError as exc:
            self._error(str(exc))

    def _read_body(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            raise HomeError("the request did not say how long its body is") from None
        if length <= 0:
            raise HomeError("the request carried no folder")
        try:
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            raise HomeError("the request body is not JSON") from None
        if not isinstance(payload, dict):
            raise HomeError("the request body is not a JSON object")
        return payload

    def _post_home(self) -> None:
        server: CampaignServer = self.server  # type: ignore[assignment]
        try:
            payload = self._read_body()
            folder = payload.get("folder")
            if not isinstance(folder, str) or not folder.strip():
                raise HomeError("no folder was named")
            outcome, path = adopt_or_write_home(folder.strip())
        except HomeError as exc:
            self._error(str(exc))
            return
        server.home = path
        self._json({outcome: str(path)})


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------


def make_server(*, home: Path | None = None, port: int = DEFAULT_PORT) -> CampaignServer:
    """Bind the loopback server.  Port 0 binds an ephemeral port for tests."""

    return CampaignServer((HOST, port), CampaignHandler, home=home)


def run_server(server: CampaignServer) -> int:
    """Announce one real, absolute URL and serve until Ctrl-C."""

    print(f"serving the campaign page at {server.url} — Ctrl-C stops")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


def _build_parser(prog: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Serve the campaign page on loopback so it can browse your folders "
            "and adopt a campaign home. It never runs any processing."
        ),
    )
    parser.add_argument(
        "--home",
        metavar="DIR-or-TOML",
        help=(
            "Campaign home: a folder holding campaign.toml, or the file itself "
            "(default: campaign.toml in the current folder, if there is one)."
        ),
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        metavar="N",
        help=f"Loopback port to bind (default: {DEFAULT_PORT}; 0 picks a free one).",
    )
    parser.add_argument(
        "--open",
        action="store_true",
        help="Open the served page in the default browser once it is up.",
    )
    return parser


def serve_main(argv: list[str] | None = None, *, prog: str = "echelle serve") -> int:
    """Run ``echelle serve``: bind 127.0.0.1 and serve until Ctrl-C."""

    args = _build_parser(prog).parse_args(argv)
    try:
        home = resolve_home(args.home)
    except HomeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    try:
        server = make_server(home=home, port=args.port)
    except OSError as exc:
        print(f"ERROR: cannot bind {HOST}:{args.port} ({exc.strerror or exc})", file=sys.stderr)
        return 1
    if args.open:
        import webbrowser

        # The server is already bound, so the page is there when the tab lands.
        if not webbrowser.open(server.url):
            print(
                f"WARNING: no default browser answered; open {server.url} by hand",
                file=sys.stderr,
            )
    return run_server(server)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(serve_main())
