"""The served campaign front door: ``echelle serve``.

The page is served, not merely written, so that it can browse the owner's real
folders and adopt a campaign home on a cold start.  Everything else about it
stays exactly as strict as the static page it serves:

* it binds ``127.0.0.1`` only, deliberately -- there is no ``--host`` flag,
  because nothing here is meant to answer a second machine.  Binding loopback
  is not by itself a boundary, though: any page in any tab can post to a
  loopback port, so every request must also arrive addressed to this server
  (:meth:`CampaignHandler._same_origin`) or it is refused unread;
* it processes nothing itself, ever.  No endpoint converts a SIF, builds a
  cube, or holds a worker.  ``POST /api/run`` *launches* one of a fixed few
  campaign verbs as its own detached process and forgets it: the run survives
  the browser tab, this server, and the operator's ``Ctrl-C``, and it leaves
  the same durable receipt a terminal run leaves, because it IS the same
  command.  Nothing crosses that endpoint but a verb NAME from the allowlist
  and the composer's two answers; a command string, a flag or an argument list
  from a request is never read and never run;
* it writes into the campaign in two places: ``POST /api/home`` creates a
  ``campaign.toml`` in a folder the operator picked, never over one that
  already exists and never creating the folder itself, and a launch writes its
  own marker and log into the campaign's runs area beside the receipts.

That is not the same as writing nothing else.  ``GET /`` rebuilds the campaign
page on every request, and building it creates the configured output folder and
writes ``index.html`` into it -- the same build ``echelle web`` does, run on a
request instead of a command.  Nothing else on disk is touched, and no request
body ever decides where that build lands: the campaign home does.

Rebuilding per request is what keeps a page left open overnight from showing
yesterday's catalog, and every response carries ``Cache-Control: no-store`` so
the browser cannot show it either.  It is also the whole of how a run is
reported: what a launch leaves on disk -- its marker, its log, its receipt --
is read back on every build, so this server holds no run state at all and a
restarted server, a second tab, or a run started in a terminal all read alike.

What this server deliberately cannot do is CONTROL a run.  There is no stop, no
pause, no kill: a page that can end a terabyte conversion halfway is a page that
can lose a night's progress, and the owner's whole requirement here was to run
work from the page without risking exactly that.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import string
import subprocess
import sys
import tempfile
import threading
import time
from collections.abc import Sequence
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

#: The host names a request may claim to have been sent to.  Anything else is
#: a name that resolved here without being this server -- the shape a DNS
#: rebinding attack has -- and is refused before it is read.
LOOPBACK_NAMES = frozenset({"127.0.0.1", "localhost", "::1"})

#: How a launch names ``echelle``: this interpreter, running the package this
#: server is itself part of.  RULES.md §14's discipline, and the reason a launch
#: does not depend on which folder's ``PATH`` the server happened to start from.
LAUNCH_PREFIX: tuple[str, ...] = (sys.executable, "-m", "echelle_spectra.cli")

#: The receipt root the campaign CLIs default to, spelled here because a launch
#: runs with the campaign home as its working folder -- so their own default
#: lands in the campaign, and the markers and logs sit beside the receipts they
#: belong to rather than in a place only this server knows about.
RUNS_DIRNAME = ("local", "runs")
LAUNCHES_DIRNAME = "launches"
LAUNCH_SCHEMA = "echelle-launch/v1"

#: How much of a failed run's log the card shows.  Enough to read the error the
#: CLI printed before it stopped; never the whole file, which for a real run is
#: a night's worth of lines.
LOG_TAIL_BYTES = 8192
LOG_TAIL_LINES = 12

RUN_RUNNING = "running"
RUN_FINISHED = "finished"
RUN_FAILED = "failed"
RUN_IDLE = "idle"


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
# Launching a campaign verb, and reading back what it left
# ---------------------------------------------------------------------------
#
# Everything here is deliberately one-directional.  A launch hands a command to
# the operating system, writes down what it handed over, and lets go; nothing in
# this process keeps a handle, a pipe, or a thread attached to it.  That is what
# makes a run survive the tab, the server and the operator's Ctrl-C -- and it is
# also why liveness is read back from the marker plus the process id rather than
# from anything this server remembers.


def runs_root(home: Path) -> Path:
    """The campaign's receipt root: the CLIs' own default, from the home."""

    return home.parent.joinpath(*RUNS_DIRNAME)


def launches_root(home: Path) -> Path:
    """Where launch markers and their logs live -- beside the receipts."""

    return runs_root(home) / LAUNCHES_DIRNAME


def _target_key(target: Path) -> str:
    """One short stable key per data folder, so one marker exists per target."""

    return hashlib.sha256(os.path.normcase(str(target)).encode("utf-8")).hexdigest()[:8]


def _argv_digest(argv: Sequence[str]) -> str:
    """Digest the exact argument list, separated by a byte no argument holds."""

    return hashlib.sha256("\0".join(argv).encode("utf-8")).hexdigest()[:16]


def marker_path(home: Path, verb: str, target: Path) -> Path:
    return launches_root(home) / f"{verb}-{_target_key(target)}.json"


def log_path(home: Path, verb: str, target: Path) -> Path:
    return launches_root(home) / f"{verb}-{_target_key(target)}.log"


def pid_alive(pid: int) -> bool:
    """Whether one recorded process id still names a running process.

    The marker plus this check is the whole liveness story: no handle is held
    anywhere, so a run outlives the server that started it and is still found
    afterwards by any server that reads the same folder.

    ``os.kill`` is never used on Windows: there, ``os.kill`` TERMINATES rather
    than probes, and a liveness check that can end the run it is asking about
    would be the exact failure the owner asked to be spared.

    An operating system that has since handed this id to somebody else would
    read as alive here.  That is why the card reports the RECEIPT once the
    process is gone: the files say whether the work finished, not the id.
    """

    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.OpenProcess.restype = ctypes.c_void_p
        kernel32.OpenProcess.argtypes = (ctypes.c_uint32, ctypes.c_int, ctypes.c_uint32)
        kernel32.CloseHandle.argtypes = (ctypes.c_void_p,)
        kernel32.GetExitCodeProcess.argtypes = (
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_ulong),
        )
        handle = kernel32.OpenProcess(0x1000, False, pid)  # QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        try:
            code = ctypes.c_ulong()
            if not kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return False
            return code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(ctypes.c_void_p(handle))
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:  # pragma: no cover - another account's process
        return True
    return True


def read_marker(path: Path) -> dict[str, Any] | None:
    """One launch marker, or None for anything this version cannot vouch for."""

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict) or payload.get("schema") != LAUNCH_SCHEMA:
        return None
    return payload


def read_markers(home: Path, verb: str) -> list[dict[str, Any]]:
    """Every marker this verb has left, newest launch first."""

    found = []
    try:
        entries = sorted(launches_root(home).glob(f"{verb}-*.json"))
    except OSError:  # pragma: no cover - the runs area vanished mid-read
        return []
    for entry in entries:
        marker = read_marker(entry)
        if marker is not None:
            found.append(marker)
    return sorted(found, key=lambda item: float(item.get("started_epoch") or 0), reverse=True)


def spawn_detached(argv: Sequence[str], *, cwd: Path, log: Path) -> int:
    """Start one command that outlives this process, and return its process id.

    Detached in the strong sense on Windows: its own process group and no
    console, so neither this server exiting nor a Ctrl-C in the console that
    started it reaches the run.  Its input is closed -- a batch run that stops
    to ask a question nobody can answer is a hang -- and both its output streams
    are appended to one log, so a failure explains itself after the fact.
    """

    log.parent.mkdir(parents=True, exist_ok=True)
    handle = open(log, "a", encoding="utf-8", errors="replace", newline="")
    try:
        handle.write(f"\n=== {utc_stamp()} launching: {' '.join(argv)}\n")
        handle.flush()
        keywords: dict[str, Any] = {}
        if os.name == "nt":
            keywords["creationflags"] = (
                getattr(subprocess, "DETACHED_PROCESS", 0x00000008)
                | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x00000200)
            )
        else:  # pragma: no cover - this campaign runs on Windows
            keywords["start_new_session"] = True
        child = subprocess.Popen(  # noqa: S603 - argv is derived, never supplied
            list(argv),
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=handle,
            stderr=subprocess.STDOUT,
            close_fds=True,
            **keywords,
        )
    finally:
        # The child holds its own duplicate of the log; this one is let go
        # immediately, so nothing about the run is attached to this server.
        handle.close()
    return int(child.pid)


def utc_stamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_marker(
    home: Path,
    *,
    verb: str,
    target: Path,
    epoch: str,
    argv: Sequence[str],
    command: str,
    pid: int,
    log: Path,
    evidence: str = "",
) -> dict[str, Any]:
    """Record what was launched, so the page can say it without asking anyone."""

    marker = {
        "schema": LAUNCH_SCHEMA,
        "verb": verb,
        "target": str(target),
        "epoch": epoch,
        "argv": list(argv),
        "command": command,
        "args_sha256": _argv_digest(argv),
        "pid": pid,
        "started_at": utc_stamp(),
        "started_epoch": time.time(),
        "log": str(log),
        # The one durable artifact a verb names on its own command line, so the
        # card can read the outcome back without re-deriving anything.
        "evidence": evidence,
    }
    path = marker_path(home, verb, target)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(marker, indent=2), encoding="utf-8", newline="\n")
    os.replace(temporary, path)
    return marker


def log_tail(path: Path) -> str:
    """The last few lines of one log, read from its end and never whole."""

    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - LOG_TAIL_BYTES))
            chunk = handle.read()
    except OSError:
        return ""
    lines = [line for line in chunk.decode("utf-8", "replace").splitlines() if line.strip()]
    return "\n".join(lines[-LOG_TAIL_LINES:])


def _receipt_since(home: Path, target: Path, since: float) -> dict[str, Any] | None:
    """The run receipt this launch produced, read exactly as ``echelle status``
    reads it: the campaign's own machinery, never a second reader of the files.

    The match is strict -- this target's own receipt, written since this launch
    began.  A bulk run whose plan names some other input folder therefore reads
    as leaving no receipt, which is the honest answer: borrowing the nearest
    receipt would let one verb report another verb's success.
    """

    from .campaign_run import list_run_summaries

    wanted = os.path.normcase(str(target))
    for summary in list_run_summaries(runs_root(home)):
        # One second of slack: a receipt written in the same tick as the launch
        # is this launch's, and clocks are not that precise.
        if float(summary.get("updated_at") or 0) < since - 1.0:
            continue
        if os.path.normcase(str(Path(summary["source_root"]))) == wanted:
            return summary
    return None


def _receipt_line(summary: dict[str, Any]) -> str:
    counts = ", ".join(
        f"{value} {name}" for name, value in sorted((summary.get("counts") or {}).items()) if value
    )
    return (
        f"receipt {summary['id']} [{summary['state']}]"
        + (f", {counts}" if counts else ", no terminal records")
    )


def _verb_state(home: Path, verb: str, *, plan: Path | None = None) -> dict[str, Any]:
    """What the files say about this verb right now: four states, never a fifth.

    Running, finished, failed, and nothing in flight are different facts, and
    none of them is inferred from a neighbour: a process that is gone without a
    receipt is ``failed`` and shows its log, never a quiet ``finished``.
    """

    markers = read_markers(home, verb)
    if not markers:
        line = "Nothing in flight."
        if verb == "bulk" and plan is not None and not plan.is_file():
            line = f"Nothing in flight; {plan.name} is not saved beside the campaign yet."
        return {"state": RUN_IDLE, "line": line, "tail": ""}
    marker = markers[0]
    target = Path(str(marker.get("target") or ""))
    started = str(marker.get("started_at") or "an unrecorded time")
    if pid_alive(int(marker.get("pid") or 0)):
        return {
            "state": RUN_RUNNING,
            "line": f"Running on {target} since {started} (pid {marker.get('pid')}).",
            "tail": "",
        }
    since = float(marker.get("started_epoch") or 0)
    if verb == "audit":
        # An audit's durable receipt is the immutable evidence file it names.
        evidence = Path(str(marker.get("evidence") or ""))
        if str(evidence) and evidence.is_file():
            try:
                verdict = json.loads(evidence.read_text(encoding="utf-8")).get("verdict")
            except (OSError, ValueError):
                verdict = None
            return {
                "state": RUN_FINISHED,
                "line": (
                    f"Finished on {target}: verdict {verdict or 'unrecorded'} "
                    f"in {evidence.name}."
                ),
                "tail": "",
            }
    else:
        summary = _receipt_since(home, target, since)
        if summary is not None:
            return {
                "state": RUN_FINISHED,
                "line": f"Finished on {target}: {_receipt_line(summary)}.",
                "tail": "",
            }
    return {
        "state": RUN_FAILED,
        "line": (
            f"Started {started} on {target}; the process is gone and left no receipt. "
            "The last lines it printed:"
        ),
        "tail": log_tail(Path(str(marker.get("log") or ""))),
    }


def _plan_name() -> str:
    """The plan file's own name, read from the composer rather than repeated."""

    from . import reading_room

    return str(getattr(reading_room, "PLAN_NAME", "campaign-plan.toml"))


def run_states(home: Path, *, plan: Path | None = None) -> dict[str, dict[str, Any]]:
    """One state per launchable verb, read off the campaign's own files."""

    from . import reading_room

    return {verb: _verb_state(home, verb, plan=plan) for verb in reading_room.RUN_ARGV}


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


def campaign_keywords(values: dict[str, Any]) -> dict[str, Any]:
    """The campaign's own files, resolved once for both a build and a launch.

    A conventional name the operator has not created yet is not a broken
    registry; it is simply not supplied.
    """

    registry = values["registry"]
    calibrations = values["calibrations"]
    return {
        "drift_paths": [path for path in values["drift"] if path.is_file()],
        "registry_path": registry if registry is not None and registry.is_file() else None,
        "calibrations_root": (
            calibrations if calibrations is not None and calibrations.is_dir() else None
        ),
    }


def campaign_view(server: CampaignServer, home: Path) -> tuple[dict[str, Any], Path | None]:
    """Resolve one home into the values a build uses, and where to read from.

    The returned ``values["catalog"]`` is always the campaign's own configured
    merged index, because that is the path every command -- composed or
    launched -- must name.  The second element overrides only where the catalog
    is READ from, for the one case that has a stand-in to render.
    """

    values = dict(home_values(home))
    if values["catalog"] is None:
        # A campaign.toml whose 'catalog' is not a usable string still has the
        # conventional name it would have had; naming that beats composing
        # against nothing.
        values["catalog"] = home.parent / "all-years.json"
    if values["catalog"].is_file():
        return values, None
    # Not scanned yet: render the full page over a synthesized empty catalog, so
    # every control -- Browse, the composer, the calibrate stepper -- exists
    # before the first scan result does.  Only the READING moves: every command
    # the page composes or launches still names the campaign's configured
    # catalog, because the first of those commands is the one that creates it.
    load_from = empty_catalog_path(server)
    if values["output"] is None or not values["output"].parent.is_dir():
        values["output"] = load_from.parent / "page"
    return values, load_from


def build_page(
    server: CampaignServer,
    values: dict[str, Any],
    *,
    load_from: Path | None = None,
    launches: dict[str, Any] | None = None,
) -> bytes:
    """Rebuild the campaign page from the home and return the bytes it wrote."""

    from .reading_room import build_reading_room

    keywords: dict[str, Any] = {**campaign_keywords(values), "document_paths": []}
    source = values["catalog"] if load_from is None else load_from
    if load_from is not None:
        keywords["compose_catalog_path"] = values["catalog"]
    with server.build_lock:
        try:
            index = build_reading_room(
                source, values["output"], served=True, launches=launches, **keywords
            )
            server.served_keyword = True
        except TypeError:
            index = build_reading_room(source, values["output"], **keywords)
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
        #: Held only across one launch, so two clicks arriving together cannot
        #: both find no marker and both spawn.  It is mutual exclusion, not
        #: state: nothing about a running child is remembered here.
        self.launch_lock = threading.Lock()
        #: Whether ``build_reading_room`` accepted ``served=True`` (None until built).
        self.served_keyword: bool | None = None
        #: Whether reading_room supplied the real setup page (None until served).
        self.setup_page_rendered: bool | None = None
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

    # -- the boundary -------------------------------------------------------

    def _addresses_this_server(self, authority: str) -> bool:
        """Whether one ``host[:port]`` names this server rather than resolves to it.

        A name that merely resolves to 127.0.0.1 is not this server: that is
        exactly the trick DNS rebinding plays, pointing its own domain at
        loopback so the page it served can read what answers there.
        """

        text = authority.strip()
        if not text:
            return False
        if text.startswith("["):  # an IPv6 literal: [::1] or [::1]:port
            closing = text.find("]")
            if closing < 0:
                return False
            name, rest = text[1:closing], text[closing + 1 :]
            if rest and not rest.startswith(":"):
                return False
            port = rest[1:]
        else:
            name, _, port = text.partition(":")
        if name.casefold() not in LOOPBACK_NAMES:
            return False
        return not port or port == str(self.server.port)  # type: ignore[attr-defined]

    def _same_origin(self) -> bool:
        """Whether this request was addressed to this server by something here.

        Loopback binding stops a second machine, and stops nothing else: any
        web page in any tab can post to a loopback port, and a form or a
        ``text/plain`` fetch does it with no preflight to refuse.  So the two
        headers a browser fills in itself are checked, and both must agree that
        this server is who was being talked to.
        """

        if not self._addresses_this_server(self.headers.get("Host") or ""):
            return False
        origin = (self.headers.get("Origin") or "").strip()
        if not origin:
            # No Origin at all is the shape of a same-origin navigation and of
            # a plain command-line client; neither is a cross-site request.
            return True
        parts = urlsplit(origin)
        if parts.scheme not in ("http", "https"):
            return False  # including the literal "null" a sandboxed frame sends
        return self._addresses_this_server(parts.netloc)

    def _guard(self) -> bool:
        """Refuse an unaddressed request in one line, before reading its body."""

        if self._same_origin():
            return True
        self._error("this server answers only the page it serves on 127.0.0.1", status=403)
        return False

    # -- routing ------------------------------------------------------------

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's name
        if not self._guard():
            return
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
        if not self._guard():
            return
        route = unquote(urlsplit(self.path).path)
        if route == "/api/home":
            self._post_home()
        elif route == "/api/run":
            self._post_run()
        else:
            self._error(f"no such endpoint: {route}", status=404)

    # -- endpoints ----------------------------------------------------------

    def _get_page(self) -> None:
        server: CampaignServer = self.server  # type: ignore[assignment]
        home = server.home
        if home is None:
            self._html(render_setup(server))
            return
        values, load_from = campaign_view(server, home)
        try:
            self._send(
                200,
                build_page(
                    server,
                    values,
                    load_from=load_from,
                    # Every build re-reads the markers, receipts and logs, so
                    # the page never remembers a run: it looks again.
                    launches=run_states(home, plan=home.parent / _plan_name()),
                ),
                "text/html; charset=utf-8",
            )
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
        # The second half of the boundary above: a request that says it carries
        # JSON is one a browser could not have sent cross-site without asking
        # this server first, and this server answers no preflight.
        media = (self.headers.get("Content-Type") or "").split(";", 1)[0].strip().casefold()
        if media != "application/json":
            raise HomeError("this endpoint reads application/json only")
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

    # -- the one endpoint that starts work ----------------------------------

    def _run_request(self, payload: dict[str, Any]) -> tuple[str, Path, str]:
        """Read the three things a run request may carry, and refuse the rest.

        A verb NAME, a data folder, and a calibration.  Not a command, not a
        flag, not an argument list: an endpoint that reads any of those is an
        endpoint that runs whatever a page in another tab decided to send.
        """

        from . import reading_room

        unexpected = sorted(set(payload) - {"verb", "folder", "epoch"})
        if unexpected:
            named = ", ".join(unexpected)
            raise HomeError(
                "this endpoint reads a verb name, a data folder and a calibration only, "
                f"and never {named}"
            )
        verb = payload.get("verb")
        if not isinstance(verb, str) or verb not in reading_room.RUN_ARGV:
            offered = ", ".join(sorted(reading_room.RUN_ARGV))
            spelled = verb if isinstance(verb, str) else ""
            raise HomeError(
                f"{spelled[:40] or 'that'} is not a verb this page can run; it runs {offered}"
            )
        raw_folder = payload.get("folder")
        if not isinstance(raw_folder, str) or not raw_folder.strip():
            raise HomeError("no data folder was named")
        folder = Path(raw_folder.strip()).expanduser()
        if not folder.is_dir():
            raise HomeError(f"not a readable folder: {folder}")
        raw_epoch = payload.get("epoch", "")
        if not isinstance(raw_epoch, str):
            raise HomeError("the calibration must be named as text")
        return verb, folder.resolve(), raw_epoch.strip()

    def _post_run(self) -> None:
        from . import reading_room

        server: CampaignServer = self.server  # type: ignore[assignment]
        home = server.home
        try:
            if home is None:
                raise HomeError("there is no campaign home yet, so there is nothing to run")
            verb, target, epoch = self._run_request(self._read_body())
            values, load_from = campaign_view(server, home)
            # The page's own derivation, run again here over the campaign's own
            # files: what a request decides is the two answers, and the whole
            # argument list follows from them and from this campaign.
            composition = reading_room.campaign_composition(
                values["catalog"] if load_from is None else load_from,
                compose_catalog_path=None if load_from is None else values["catalog"],
                **campaign_keywords(values),
            )
            known = {str(item) for item in composition["epochs"]}
            known |= {str(item) for item in composition["registry"].get("saved") or []}
            if epoch and epoch not in known:
                raise HomeError(f"this campaign knows no calibration named {epoch}")
            composed = reading_room.with_answers(
                composition["values"],
                folder=str(target),
                epoch=epoch,
                evidence_name=composition["evidence_name"],
            )
            argv = reading_room.run_argv(verb, composed)
            command = reading_room.run_command(verb, composed)
        except HomeError as exc:
            self._error(str(exc))
            return

        with server.launch_lock:
            running = read_marker(marker_path(home, verb, target))
            if running is not None and pid_alive(int(running.get("pid") or 0)):
                self._error(
                    f"{verb} is already running on {running.get('target')} "
                    f"(pid {running.get('pid')}, started {running.get('started_at')})",
                    status=409,
                )
                return
            log = log_path(home, verb, target)
            try:
                pid = spawn_detached(
                    [*LAUNCH_PREFIX, *argv], cwd=home.parent, log=log
                )
            except OSError as exc:
                self._error(f"this machine could not start the run: {exc.strerror or exc}")
                return
            marker = write_marker(
                home,
                verb=verb,
                target=target,
                epoch=epoch,
                argv=argv,
                command=command,
                pid=pid,
                log=log,
                evidence=composed["verdict"] if verb == "audit" else "",
            )
        self._json({"launched": marker})


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
