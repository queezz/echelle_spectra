"""Build the read-only campaign flow, catalog, evidence and reading room.

The page is a one-shot static build: ``echelle web`` writes one ``index.html``
that carries its own CSS and JavaScript and fetches nothing.  It is organized
by the work rather than by the data: four tabs — **Now** (the campaign as a
stepper, one independent row per drive), **Drives** (the catalog),
**Calibration** (epochs and drift evidence in sequence position) and
**Reading room** (the packaged canon).  Its structure follows the house web UI
law (fleet's ``WEBUI.md`` and ``WEBUI-COOKBOOK.md``):

* the main column is what a person reads; the rails are what a person presses,
  with controls on the left and context plus the section index on the right,
  each rail carrying only the active tab's own cargo;
* a top tab press always returns that tab to its home state, including when a
  sub-view inside it is open (the view-nesting law's guard exception);
* the teaching lives once, in each tab's own legend and in the reading room;
  cards carry facts — chips, counts, and one truncated path;
* both rails are sticky at an offset derived as the sum of the page's own
  header metrics, and only the rails carry ``align-self: start`` so their grid
  row still spans the content column's full height (the zero-travel trap);
* every jump target's id starts with ``sec-`` so one ``scroll-margin-top`` rule
  covers all of them, and the scroll-spy query is scoped to the content column;
* command rows lead with plain words and keep the literal command behind a
  toggle, in both shell shapes, with the copy button carrying the full payload,
  and they are read in the order the campaign runs them;
* the composer asks two things -- the folder of SIF shots and the calibration
  -- because those are the only two a person holds that no file here records;
  every other value is derived from them and folded away, editable;
* empty, unmeasured, unreachable and judged states are rendered distinctly, and
  a verdict word this page does not know is rendered as unrecognized rather
  than dressed as one it does.

When a local server serves the same page (``served=True``) it gains one power
the file cannot have: it may ask that server which folders exist on this
machine.  That is the whole served half — a folder picker behind a Browse
button on the composer's data-folder field, and the two one-screen pages a cold
campaign is served (:func:`render_setup_page` when no campaign home is recorded
yet, :func:`render_empty_campaign_page` when a home has no catalog yet).  Every
one of those pieces is *appended* to the static build rather than woven into
it, so the one-shot file keeps its "fetches nothing" contract byte for byte.
"""

from __future__ import annotations

import html
import json
import os
import re
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import Any

from .campaign_run import GATE_SAMPLE, GATE_UNGATED, GATE_UNRECORDED, GATE_VERDICT
from .catalog import load_catalog, source_catalog_path

#: The three documents whose canonical text lives inside the installed package.
#: ``echelle web`` renders these from ``importlib.resources``, so a travel kit
#: with no repository checkout still opens a complete reading room.
PACKAGED_DOCUMENTS = ("vocabulary.md", "procedure.md", "provenance.md")


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
#
# A deliberately small, stdlib-only renderer.  Supported: ATX headings, fenced
# and inline code, unordered and ordered lists, pipe tables, horizontal rules,
# links, bold, italic, and paragraphs.  Everything else is rendered as text.

_FENCE = re.compile(r"^```\s*[\w+-]*\s*$")
_HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$")
_RULE = re.compile(r"^(?:-{3,}|\*{3,})$")
_UL_ITEM = re.compile(r"^[-*]\s+(.*)$")
_OL_ITEM = re.compile(r"^\d+\.\s+(.*)$")
_TABLE_RULE = re.compile(r"^\|?\s*:?-{2,}:?\s*(?:\|\s*:?-{2,}:?\s*)*\|?$")
_CODE_SPAN = re.compile(r"`([^`]+)`")
_LINK = re.compile(r"\[([^\]]+)\]\(([^)\s]+)\)")
_STRONG = re.compile(r"\*\*(.+?)\*\*")
_EMPHASIS = re.compile(r"(?<![\w*_])[*_]([^*_\n]+)[*_](?![\w*_])")
_SLUG = re.compile(r"[^a-z0-9]+")


def _e(value: Any) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def _href(target: str) -> str:
    """Return a link target this offline page may keep.

    Anything carrying a scheme is dropped to a fragment: a rendered document's
    own outbound link is a promise about a machine this build never saw, and
    the page's one address is written by the page itself, not by a document.
    """

    cleaned = target.strip()
    return "#" if ":" in cleaned.split("/")[0] else cleaned


def _inline(text: str) -> str:
    """Render one line's inline markup on already-untrusted text."""

    spans: list[str] = []

    def _stash(match: re.Match[str]) -> str:
        spans.append(html.escape(match.group(1)))
        return f"\x00{len(spans) - 1}\x00"

    staged = html.escape(_CODE_SPAN.sub(_stash, text))
    staged = _LINK.sub(
        lambda match: f'<a href="{_e(_href(match.group(2)))}">{match.group(1)}</a>', staged
    )
    staged = _STRONG.sub(r"<strong>\1</strong>", staged)
    staged = _EMPHASIS.sub(r"<em>\1</em>", staged)
    for index, span in enumerate(spans):
        staged = staged.replace(f"\x00{index}\x00", f"<code>{span}</code>")
    return staged


def _unique_id(candidate: str, seen: set[str]) -> str:
    identifier, suffix = candidate, 2
    while identifier in seen:
        identifier = f"{candidate}-{suffix}"
        suffix += 1
    seen.add(identifier)
    return identifier


def _heading(line: str, seen: set[str], prefix: str) -> str:
    match = _HEADING.match(line)
    assert match is not None  # only called after the dispatcher matched
    level = min(len(match.group(1)) + 2, 6)
    text = match.group(2).strip()
    slug = _SLUG.sub("-", text.lower()).strip("-") or "section"
    return f'<h{level} id="{_e(_unique_id(prefix + slug, seen))}">{_inline(text)}</h{level}>'


def _fence_block(lines: list[str], index: int) -> tuple[str, int]:
    body: list[str] = []
    index += 1
    while index < len(lines) and not _FENCE.match(lines[index]):
        body.append(lines[index])
        index += 1
    return f"<pre><code>{html.escape(chr(10).join(body))}</code></pre>", min(index + 1, len(lines))


def _list_block(lines: list[str], index: int) -> tuple[str, int]:
    ordered = _OL_ITEM.match(lines[index]) is not None
    pattern = _OL_ITEM if ordered else _UL_ITEM
    items: list[str] = []
    while index < len(lines):
        match = pattern.match(lines[index])
        if match is not None:
            items.append(match.group(1).strip())
        elif items and lines[index].strip() and lines[index].startswith((" ", "\t")):
            # A wrapped continuation of the item above, not a new block.
            items[-1] = f"{items[-1]} {lines[index].strip()}"
        else:
            break
        index += 1
    tag = "ol" if ordered else "ul"
    body = "".join(f"<li>{_inline(item)}</li>" for item in items)
    return f"<{tag}>{body}</{tag}>", index


def _table_cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def _table_block(lines: list[str], index: int) -> tuple[str, int]:
    head = "".join(f"<th>{_inline(cell)}</th>" for cell in _table_cells(lines[index]))
    index += 2
    rows: list[str] = []
    while index < len(lines) and lines[index].strip() and "|" in lines[index]:
        cells = "".join(f"<td>{_inline(cell)}</td>" for cell in _table_cells(lines[index]))
        rows.append(f"<tr>{cells}</tr>")
        index += 1
    table = f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    return f'<div class="scroll-x">{table}</div>', index


def _starts_block(line: str) -> bool:
    return bool(
        _FENCE.match(line)
        or _HEADING.match(line)
        or _RULE.match(line.strip())
        or _UL_ITEM.match(line)
        or _OL_ITEM.match(line)
    )


def _paragraph(lines: list[str], index: int) -> tuple[str, int]:
    body = [lines[index].strip()]
    index += 1
    while index < len(lines) and lines[index].strip() and not _starts_block(lines[index]):
        body.append(lines[index].strip())
        index += 1
    return f"<p>{_inline(' '.join(body))}</p>", index


def render_markdown(
    text: str, *, heading_ids: set[str] | None = None, id_prefix: str = "doc-"
) -> str:
    """Render the documented Markdown subset to HTML.

    Heading ids are prefixed (``doc-`` by default) rather than ``sec-`` on
    purpose: ``sec-`` is the page's own anchor and scroll-spy namespace, and
    ids slugified from rendered document text must not be able to collide with
    it (the cookbook's scoped-scroll-spy rule).
    """

    seen = heading_ids if heading_ids is not None else set()
    lines = text.replace("\r\n", "\n").replace("\x00", "").split("\n")
    parts: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.strip():
            index += 1
        elif _FENCE.match(line):
            chunk, index = _fence_block(lines, index)
            parts.append(chunk)
        elif _HEADING.match(line):
            parts.append(_heading(line, seen, id_prefix))
            index += 1
        elif _RULE.match(line.strip()):
            parts.append("<hr>")
            index += 1
        elif _UL_ITEM.match(line) or _OL_ITEM.match(line):
            chunk, index = _list_block(lines, index)
            parts.append(chunk)
        elif "|" in line and index + 1 < len(lines) and _TABLE_RULE.match(lines[index + 1].strip()):
            chunk, index = _table_block(lines, index)
            parts.append(chunk)
        else:
            chunk, index = _paragraph(lines, index)
            parts.append(chunk)
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Honesty states
# ---------------------------------------------------------------------------

#: Four different facts about a drive, four different renderings.  Nothing here
#: falls back to a neighbouring state: a drive whose catalog did not answer is
#: never shown as one that answered with nothing.
_SOURCE_STATES = {
    "missing-drive": ("state-missing-drive", "missing drive"),
    "unmeasured": ("state-unmeasured", "unmeasured"),
    "empty": ("state-empty", "empty"),
    "measured": ("state-measured", "measured"),
}

#: What each state means, said once on the page: these sentences are the drive
#: legend's text and appear nowhere else.  A card carries facts, not teaching.
_SOURCE_NOTES = {
    "missing-drive": (
        "This drive's catalog did not answer when the page was built. Its rows are the "
        "last merge's rows; nothing is claimed about the files themselves."
    ),
    "unmeasured": (
        "No run receipt travelled with this catalog, so this page cannot say what was "
        "processed here or under what authorization."
    ),
    "empty": "A run was recorded and it published no cubes. That is a measured zero.",
    "measured": (
        "The catalog answered, a receipt described the run, and cubes were published."
    ),
}

_GATE_STATES = {
    GATE_VERDICT: ("gate-verdict", "verdict-authorized"),
    GATE_SAMPLE: ("gate-sample", "unverified sample"),
    GATE_UNGATED: ("gate-ungated", GATE_UNGATED),
    GATE_UNRECORDED: ("gate-unrecorded", GATE_UNRECORDED),
}

_VERDICT_STATES = {
    "aligned": ("verdict-aligned", "aligned"),
    "shifted": ("verdict-shifted", "shifted"),
    "misaligned-beyond-repair": ("verdict-beyond-repair", "misaligned-beyond-repair"),
    "insufficient-data": ("state-insufficient-data", "insufficient-data"),
}


def source_state(source: dict[str, Any]) -> str:
    """Name one merged catalog row's honesty state."""

    if not source.get("available"):
        return "missing-drive"
    if not source.get("run"):
        return "unmeasured"
    if not source.get("cubes"):
        return "empty"
    return "measured"


def _gate_state(value: str) -> tuple[str, str]:
    return _GATE_STATES.get(str(value), ("gate-unrecognized", f"unrecognized: {value}"))


def _verdict_state(value: Any) -> tuple[str, str]:
    return _VERDICT_STATES.get(str(value), ("state-unrecognized", f"unrecognized: {value}"))


def _refresh_availability(catalog: dict[str, Any]) -> dict[str, Any]:
    if catalog.get("schema") != "echelle-merged-catalog/v1":
        source = {
            "volume_label": catalog["volume_label"],
            "catalog_path": "",
            "available": True,
            "run": catalog.get("run"),
            "cubes": catalog.get("cubes", []),
        }
        return {"schema": "echelle-merged-catalog/v1", "sources": [source]}
    for source in catalog.get("sources", []):
        # Merged rows store the catalog path relative to the drive root, so
        # availability resolves against the root this machine last saw.
        source["available"] = bool(
            source.get("catalog_path") and source_catalog_path(source).is_file()
        )
    return catalog


# ---------------------------------------------------------------------------
# Command composition
# ---------------------------------------------------------------------------
#
# One template per command, filled the same way twice: here, for the page's
# initial state, and in the page's own JavaScript when the reader presses
# Compose.  Keeping the template as the single source is what stops the two
# from drifting apart.

PLAN_TEMPLATE = """# Composed by the Echelle reading room. Nothing here has been run.
# Save this as {{plan}} beside the campaign, read it once, then paste the
# process command below. Relative paths resolve against this file's folder.
# The registry selects the epoch; for these shots that is {{epoch}}.
# Discovery walks input_dir recursively; calibration folders are left out.

[plan]
input_dir = "{{input}}"
output_dir = "{{output}}"
pattern = "{{pattern}}"
registry = "{{registry}}"
calibrations = "{{calibrations}}"
drift_verdict = "{{verdict}}"
central_index = "{{catalog}}"
"""

#: The sequence in campaign order: sample, check the alignment, then convert in
#: bulk.  The order the rows are read in is the order the work runs in, and the
#: bench stays after them as the calibration aside rather than a fourth step.
COMMAND_TEMPLATES = (
    {
        "id": "cmd-sample",
        "title": "1. Sample this data folder",
        "meaning": (
            "Convert a sample of the folder's SIF shots into cubes under the registry, so "
            "there is something to measure the wavelength alignment on. The sample size is "
            "derived from what the folder holds, and every cube it writes is marked an "
            "unverified sample."
        ),
        "powershell": (
            'echelle process "{{input}}" -o "{{output}}" --registry "{{registry}}" '
            '--calibrations "{{calibrations}}" --sample auto --volume-label "{{label}}" '
            '--central-index "{{catalog}}"'
        ),
        "posix": (
            'echelle process "{{input}}" -o "{{output}}" --registry "{{registry}}" '
            '--calibrations "{{calibrations}}" --sample auto --volume-label "{{label}}" '
            '--central-index "{{catalog}}"'
        ),
    },
    {
        "id": "cmd-audit",
        "title": "2. Check the wavelength alignment",
        "meaning": (
            "We do not trust the calibration on this data, we check it: measure Balmer and "
            "Fulcher centroids on the cubes just sampled and write one immutable verdict "
            "file. The name is the next free one, so an audit never overwrites evidence "
            "already on the drive."
        ),
        "powershell": (
            'echelle drift audit "{{output}}" --catalog "{{catalog}}" '
            '--calibrations "{{calibrations}}" -o "{{verdict}}"'
        ),
        "posix": (
            'echelle drift audit "{{output}}" --catalog "{{catalog}}" '
            '--calibrations "{{calibrations}}" -o "{{verdict}}"'
        ),
    },
    {
        "id": "cmd-process",
        "title": "3. Generate cubes from the whole folder",
        "meaning": (
            "The product. This step is unlocked by the verdict the audit above wrote: a "
            "registry-backed bulk run is refused until that evidence exists, and the plan "
            "file names it beside the input, the output, the registry and the snapshot "
            "root, so one resumable receipt records what authorized the run."
        ),
        "powershell": 'echelle process --plan "{{plan}}"',
        "posix": 'echelle process --plan "{{plan}}"',
    },
    {
        "id": "cmd-bench",
        "title": "Calibration aside — check the snapshot by eye",
        "meaning": (
            "Open the live calibration bench on the snapshot root and look at the "
            "current epoch yourself before trusting any verdict this page renders. "
            "The bench reads; it does not change a saved snapshot."
        ),
        "powershell": 'echelle-calib "{{calibrations}}"',
        "posix": 'echelle-calib "{{calibrations}}"',
    },
)

SHELL_NAMES = (("powershell", "PowerShell"), ("posix", "POSIX shell"))


def _shell_path(shell: str, value: str) -> str:
    """Write one value in the path shape the named shell reads."""

    text = str(value)
    return text.replace("/", "\\") if shell == "powershell" else text.replace("\\", "/")


def _fill(template: str, values: dict[str, str], shell: str) -> str:
    filled = template
    for key, value in values.items():
        filled = filled.replace("{{" + key + "}}", _shell_path(shell, value))
    return filled


def _posix(path: str | Path) -> str:
    return str(path).replace("\\", "/")


#: The data folder is the one thing no file on this machine records, so it is
#: the one thing the page asks for.  Until it is answered every value derived
#: from it renders as this marker: visibly unfilled, never an invented path and
#: never a drive letter this build has not seen.
UNFILLED_FOLDER = "<data folder>"

_EVIDENCE_NUMBER = re.compile(r"drift-evidence-(\d+)\.json$", re.IGNORECASE)


def _next_evidence_name(drift: list[dict[str, Any]]) -> str:
    """Name the evidence file the next audit writes, never one already read.

    Drift evidence is immutable: an audit refuses to overwrite one.  So the
    composed command names the next free number rather than the file this build
    was handed -- ``001`` when this build knows of none.
    """

    highest = 0
    for entry in drift:
        match = _EVIDENCE_NUMBER.search(_posix(entry["path"]))
        if match is not None:
            highest = max(highest, int(match.group(1)))
    return f"drift-evidence-{highest + 1:03d}.json"


def _derived_from_folder(folder: str, evidence_name: str) -> dict[str, str]:
    """Everything the data folder decides, derived in one place.

    ``derivedFrom`` in the page's own JavaScript is the other half of this one
    rule; the two must answer the same for the same folder, which is why the
    rule is written once here and mirrored there rather than spelled twice.
    """

    marked = folder.strip() or UNFILLED_FOLDER
    trimmed = _posix(marked).rstrip("/")
    return {
        "input": marked,
        # Cubes and catalogs belong on the drive beside their own data.
        "output": marked,
        "cubes": marked,
        "label": trimmed.rsplit("/", 1)[-1] or marked,
        "verdict": f"{trimmed}/{evidence_name}",
    }


def _composer_values(
    *,
    catalog_path: str | Path,
    registry: dict[str, Any],
    evidence_name: str,
    epochs: list[str],
) -> dict[str, str]:
    """Pre-fill the composer from what the page actually carries.

    Two answers compose a campaign -- the data folder and the calibration --
    and every value here is derived from those two plus this build's own
    inputs.  Nothing is baked in: no example drive ships with the page.
    """

    return {
        **_derived_from_folder("", evidence_name),
        "pattern": "*.SIF",
        "registry": registry.get("path") or "calibration_registry.toml",
        "calibrations": registry.get("calibrations") or "calibrations",
        "catalog": _posix(catalog_path),
        "plan": "campaign-plan.toml",
        "epoch": epochs[0] if epochs else "unassigned",
    }


# ---------------------------------------------------------------------------
# The campaign pipeline
# ---------------------------------------------------------------------------
#
# The Now tab renders the campaign the way it is worked: one calibrate stage
# per campaign, then one row per connected drive, each advancing on its own
# because the work is one worker per drive and several drives at once.
#
# Every step state below is computed at build time from files this build
# already reads -- the merged catalog, the run receipts summarised inside it,
# the registry, and the drift evidence handed to ``--drift``.  Nothing is
# inferred past what those files say: a step whose completion no file records
# is rendered ``not recorded`` rather than guessed from a neighbouring fact,
# and a blocked step names both what is missing and the step that supplies it.

STEP_DONE = "done"
STEP_READY = "ready"
STEP_BLOCKED = "blocked"
STEP_UNRECORDED = "not recorded"

_STEP_CLASSES = {
    STEP_DONE: "step-done",
    STEP_READY: "step-ready",
    STEP_BLOCKED: "step-blocked",
    STEP_UNRECORDED: "step-unrecorded",
}

#: Meaning and command template for every step a reader can act on.  The
#: templates are filled per drive, so a ready step is the command for exactly
#: that drive's next move rather than a generic example.
STEP_COMMANDS = {
    "bench": (
        "Open the live bench on the snapshot root and fit sphere plus every lamp you "
        "measured; the bench writes the snapshot the registry then names.",
        'echelle-calib "{{calibrations}}"',
    ),
    "connect": (
        "Catalog this drive where it is plugged in now, so the index finds its cubes again.",
        'echelle catalog build "{{cubes}}" --volume-label "{{label}}"',
    ),
    "sample": (
        "Process a sample of the data folder as an unverified sample — the legal first "
        "registry run, and the cubes the drift audit then measures. The count is derived "
        "from what the folder holds.",
        'echelle process "{{input}}" -o "{{output}}" --registry "{{registry}}" '
        '--calibrations "{{calibrations}}" --sample auto --volume-label "{{label}}" '
        '--central-index "{{catalog}}"',
    ),
    "audit": (
        "Measure Balmer and Fulcher centroids on this drive's sampled cubes and write one "
        "immutable verdict file under the next free name.",
        'echelle drift audit "{{cubes}}" --catalog "{{catalog}}" '
        '--calibrations "{{calibrations}}" -o "{{verdict}}"',
    ),
    "audit-again": (
        "Audit every cube on this drive: the last sample could not carry a verdict.",
        'echelle drift audit "{{cubes}}" --every 1 --catalog "{{catalog}}" '
        '--calibrations "{{calibrations}}" -o "{{verdict}}"',
    ),
    "cubes": (
        "Generate the cubes — the campaign's product. The plan supplies input, output, "
        "registry, snapshot root and the accepted verdict, so one receipt records the "
        "evidence that authorized the run.",
        'echelle process --plan "{{plan}}"',
    ),
    "txt": (
        "Write LHD text at the frozen legacy header from one cube; repeat per cube LHD asks for.",
        'echelle txt "{{cube}}" "{{txt}}"',
    ),
    "merge": (
        "Fold this drive's catalog into the all-years index so the next build remembers it.",
        'echelle catalog merge "{{drive_catalog}}" -o "{{catalog}}"',
    ),
}


def _step(
    key: str,
    name: str,
    state: str,
    fact: str,
    *,
    evidence: str = "",
    command: str = "",
    values: dict[str, str] | None = None,
) -> dict[str, Any]:
    """One step of the flow: its state, its one fact, and its own command."""

    meaning, template = STEP_COMMANDS[command] if command else ("", "")
    shapes = (
        {shell: _fill(template, values or {}, shell) for shell, _ in SHELL_NAMES}
        if command
        else {}
    )
    return {
        "key": key,
        "name": name,
        "state": state,
        "fact": fact,
        "evidence": evidence,
        "meaning": _fill(meaning, values or {}, "posix") if meaning else "",
        "shapes": shapes,
    }


def _primary_step(steps: list[dict[str, Any]]) -> int:
    """Return the index of the step that answers "what do I do first", or -1.

    The first step a reader can act on wins: ready or blocked.  A step whose
    completion no file records is never the answer — it can never become done
    from this page's side, so it would pin the arrow to itself for the rest of
    the campaign, and a drive is done when its cubes exist and are catalogued
    rather than when its optional LHD text was written.
    """

    for position, step in enumerate(steps):
        if step["state"] in {STEP_READY, STEP_BLOCKED}:
            return position
    return -1


def _epoch_covering(rows: list[dict[str, Any]], today: str) -> dict[str, Any] | None:
    for row in rows:
        start, end = row.get("date_from") or "", row.get("date_to") or ""
        if not start and not end:
            continue
        if (not start or start <= today) and (not end or today <= end):
            return row
    return None


def _calibrate_steps(
    registry: dict[str, Any],
    snapshot_ids: list[str],
    values: dict[str, str],
    today: str,
) -> list[dict[str, Any]]:
    """Compute the campaign's one calibrate stage from registry and catalog."""

    rows = list(registry.get("epoch_rows") or [])
    named = [str(row["snapshot_id"]) for row in rows]
    saved = [item for item in snapshot_ids if item and item != "unassigned"]
    on_disk = [str(item) for item in (registry.get("saved") or []) if item]
    unreadable = registry["status"] == "unreadable"
    absent = registry["status"] == "not supplied"
    # One fact decides the first three steps, because no file records the
    # physical session between them: the snapshot is the only trace it leaves —
    # and a snapshot folder on the calibrations root is exactly such a trace,
    # registry or not.
    trace = named or saved or on_disk
    if trace:
        if named:
            origin = f"{len(named)} registry epoch(s): {', '.join(named)}"
        elif saved:
            origin = f"cubes name snapshot(s) {', '.join(sorted(set(saved)))}"
        else:
            origin = (
                f"saved snapshot folder(s) on the calibrations root: "
                f"{', '.join(on_disk)} — not yet in any registry"
            )
        lamps = _step("lamps", "Sphere + lamps", STEP_DONE, origin)
        fit = _step("fit", "Bench fit", STEP_DONE, "the saved snapshot is the fit")
        snapshot = _step("snapshot", "Snapshot saved", STEP_DONE, origin)
    elif unreadable:
        detail = str(registry.get("detail") or "the registry could not be read")
        lamps = _step("lamps", "Sphere + lamps", STEP_BLOCKED, detail)
        fit = _step("fit", "Bench fit", STEP_BLOCKED, detail)
        snapshot = _step("snapshot", "Snapshot saved", STEP_BLOCKED, detail)
    else:
        missing = (
            "no registry was given to this build, so no snapshot is in reach"
            if absent
            else "no snapshot is in reach"
        )
        lamps = _step("lamps", "Sphere + lamps", STEP_UNRECORDED, missing)
        fit = _step(
            "fit",
            "Bench fit",
            STEP_READY,
            "nothing has been fitted where this page can see it",
            command="bench",
            values=values,
        )
        snapshot = _step("snapshot", "Snapshot saved", STEP_UNRECORDED, missing + "; the bench writes it")
    if named:
        covering = _epoch_covering(rows, today)
        dated = [row for row in rows if row.get("date_from") or row.get("date_to")]
        if covering is not None:
            epoch = _step(
                "epoch",
                "Registry epoch",
                STEP_DONE,
                f"{covering['snapshot_id']} covers today ({today})",
            )
        elif not dated:
            epoch = _step(
                "epoch",
                "Registry epoch",
                STEP_DONE,
                f"{len(named)} epoch(s), all shot-bounded — coverage of today is not recorded",
            )
        else:
            epoch = _step(
                "epoch",
                "Registry epoch",
                STEP_READY,
                f"no epoch covers today ({today}); this campaign needs its own snapshot",
                command="bench",
                values=values,
            )
    elif unreadable:
        epoch = _step("epoch", "Registry epoch", STEP_BLOCKED, str(registry.get("detail") or ""))
    elif saved:
        epoch = _step(
            "epoch",
            "Registry epoch",
            STEP_BLOCKED,
            "no registry names these snapshots — supply one with --registry",
        )
    else:
        epoch = _step(
            "epoch",
            "Registry epoch",
            STEP_UNRECORDED,
            "not supplied to this build; Snapshot saved supplies its entry",
        )
    return [lamps, fit, snapshot, epoch]


def _cube_names(source: dict[str, Any]) -> set[str]:
    return {str(cube.get("path", "")).rsplit("/", 1)[-1] for cube in source.get("cubes", [])}


def _evidence_for_drive(
    source: dict[str, Any], drift: list[dict[str, Any]]
) -> dict[str, Any] | None:
    """Attach one drift verdict to the drive whose cubes it actually measured.

    Evidence names the cubes it sampled by file name and nothing else, so that
    name is the only honest join back to a drive.  Evidence whose cubes match no
    drive stays unattributed and is read in the Calibration tab instead of being
    hung on a drive by a shared epoch, which several drives can share.
    """

    names = _cube_names(source)
    matched = [
        entry
        for entry in drift
        if names
        & {
            str(item.get("cube", ""))
            for key in ("sampled_cubes", "per_shot")
            for item in (entry["evidence"].get(key) or [])
        }
    ]
    if not matched:
        return None
    return max(matched, key=lambda entry: str(entry["evidence"].get("created_at") or ""))


def _drive_values(
    base: dict[str, str], source: dict[str, Any], *, evidence_name: str
) -> dict[str, str]:
    """Fill the command templates for exactly one drive."""

    root = _posix(source.get("drive_root", "")) or base["output"]
    label = str(source.get("volume_label", "unknown"))
    cubes = [str(cube.get("path", "")) for cube in source.get("cubes", [])]
    first = cubes[0] if cubes else ""
    values = dict(base)
    values.update(
        {
            "output": root,
            "cubes": root,
            "label": label,
            # An audit writes; it never overwrites the evidence it was handed,
            # so the name here is the next free one on this drive.
            "verdict": _posix(Path(root) / evidence_name),
            "drive_catalog": _posix(
                Path(root) / (str(source.get("catalog_path") or "echelle-catalog.json"))
            ),
            "cube": _posix(Path(root) / first) if first else f"{root}/<cube>.nc",
            "txt": (
                _posix(Path(root) / first)[:-3] + ".txt"
                if first.endswith(".nc")
                else f"{root}/<cube>.txt"
            ),
        }
    )
    return values


def _evidence_link(entry: dict[str, Any], anchor: str) -> str:
    """Link one drive's step to the evidence card that proves it.

    The label is the file's own name and the whole path stays in the title, the
    same compression every other path on the page gets.
    """

    path = str(entry["path"])
    return (
        f'<a class="xlink" href="#{_e(anchor)}" data-tab="calibration" '
        f'title="{_e(path)}">{_e(path.rsplit("/", 1)[-1])}</a>'
    )


def _verdict_step(
    evidence: dict[str, Any] | None, values: dict[str, str], anchor: str, gate: str
) -> dict[str, Any]:
    """Render the branch this drive's verdict actually points at."""

    if evidence is None:
        if gate == GATE_VERDICT:
            # The receipt proves a verdict authorized this run; the file that
            # carried it was simply not given to this build.
            return _step(
                "verdict",
                "Verdict",
                STEP_UNRECORDED,
                "the receipt was verdict-authorized; no evidence file was given to this build",
            )
        return _step(
            "verdict",
            "Verdict",
            STEP_BLOCKED,
            "no verdict for this drive — Drift audit writes one",
        )
    payload = evidence["evidence"]
    verdict = str(payload.get("verdict", ""))
    link = _evidence_link(evidence, anchor)
    if verdict == "aligned":
        return _step("verdict", "Verdict", STEP_DONE, "aligned — bulk is authorized", evidence=link)
    if verdict == "insufficient-data":
        return _step(
            "verdict",
            "Verdict",
            STEP_READY,
            "insufficient-data — sample more cubes and audit again",
            evidence=link,
            command="audit-again",
            values=values,
        )
    if verdict == "shifted":
        step = _step(
            "verdict",
            "Verdict",
            STEP_READY,
            "shifted — refine, then repoint the registry",
            evidence=link,
        )
        repair = next(
            (item for item in (payload.get("repair_commands") or []) if item.get("command")),
            None,
        )
        if repair is not None:
            command = str(repair["command"])
            shell = str(repair.get("shell", "any"))
            step["meaning"] = str(repair.get("purpose", "accept the sampled shift"))
            step["shapes"] = (
                {"powershell": command, "posix": command}
                if shell == "any"
                else {shell: command}
            )
        return step
    if verdict == "misaligned-beyond-repair":
        return _step(
            "verdict",
            "Verdict",
            STEP_BLOCKED,
            str(payload.get("data_requirement") or "beyond repair — the raw SIF data is needed"),
            evidence=link,
        )
    return _step(
        "verdict",
        "Verdict",
        STEP_UNRECORDED,
        f"unrecognized verdict: {verdict or 'none recorded'}",
        evidence=link,
    )


def _drive_steps(
    source: dict[str, Any],
    *,
    evidence: dict[str, Any] | None,
    anchor: str,
    values: dict[str, str],
    merged: bool,
    catalog_path: str,
) -> list[dict[str, Any]]:
    """Compute one drive's own position in the per-drive loop."""

    available = bool(source.get("available"))
    cubes = source.get("cubes", [])
    run = source.get("run") or {}
    gate = str(run.get("gate", GATE_UNRECORDED)) if run else ""
    steps: list[dict[str, Any]] = []

    drive_id = str(source.get("drive_id") or "")
    steps.append(
        _step(
            "connect",
            "Connect + identify",
            STEP_DONE,
            f"catalog answered · id {drive_id or 'not announced, label only'}",
        )
        if available
        else _step(
            "connect",
            "Connect + identify",
            STEP_READY,
            "this drive's catalog did not answer when the page was built",
            command="connect",
            values=values,
        )
    )

    if cubes:
        steps.append(
            _step(
                "sample",
                "Sample N",
                STEP_DONE,
                f"{len(cubes)} cube(s) catalogued"
                + (" · receipt marks an unverified sample" if run.get("sample") else ""),
            )
        )
    elif available:
        steps.append(
            _step(
                "sample",
                "Sample N",
                STEP_READY,
                "no cube is catalogued for this drive yet",
                command="sample",
                values=values,
            )
        )
    else:
        steps.append(
            _step("sample", "Sample N", STEP_BLOCKED, "needs the drive itself — Connect + identify")
        )

    if evidence is not None:
        steps.append(
            _step(
                "audit",
                "Drift audit",
                STEP_DONE,
                f"{len(evidence['evidence'].get('sampled_cubes') or [])} cube(s) sampled",
                evidence=_evidence_link(evidence, anchor),
            )
        )
    elif gate == GATE_VERDICT:
        steps.append(
            _step(
                "audit",
                "Drift audit",
                STEP_DONE,
                "the receipt records a verdict-authorized run; its evidence file is not "
                "named by this catalog",
            )
        )
    elif cubes and available:
        steps.append(
            _step(
                "audit",
                "Drift audit",
                STEP_READY,
                "cubes are on this drive and no verdict measures them",
                command="audit",
                values=values,
            )
        )
    else:
        steps.append(
            _step("audit", "Drift audit", STEP_BLOCKED, "needs sampled cubes — Sample N supplies them")
        )

    steps.append(_verdict_step(evidence, values, anchor, gate))
    verdict_state = steps[-1]["state"]

    bulk_done = bool(cubes) and bool(run) and not run.get("sample") and run.get("state") == "completed"
    if bulk_done:
        counts = ", ".join(f"{value} {key}" for key, value in sorted((run.get("counts") or {}).items()))
        steps.append(
            _step(
                "cubes",
                "Generate cubes",
                STEP_DONE,
                f"{len(cubes)} cube(s) · receipt {run.get('id')} [{run.get('state')}]"
                + (f" · {counts}" if counts else ""),
            )
        )
    elif cubes and not run:
        steps.append(
            _step(
                "cubes",
                "Generate cubes",
                STEP_UNRECORDED,
                "cubes exist and no receipt says which run made them",
            )
        )
    elif verdict_state == STEP_DONE and available:
        steps.append(
            _step(
                "cubes",
                "Generate cubes",
                STEP_READY,
                "the verdict authorizes the bulk run — this is the product",
                command="cubes",
                values=values,
            )
        )
    else:
        steps.append(
            _step(
                "cubes",
                "Generate cubes",
                STEP_BLOCKED,
                "needs an aligned or refined verdict — Verdict supplies it",
            )
        )

    steps.append(
        _step(
            "txt",
            "LHD txt",
            STEP_UNRECORDED,
            "the side deliverable: no receipt or catalog field records a txt export",
            command="txt",
            values=values,
        )
        if cubes
        else _step("txt", "LHD txt", STEP_BLOCKED, "needs cubes — Generate cubes supplies them")
    )

    steps.append(
        _step("merge", "Catalog merge", STEP_DONE, f"in {catalog_path.rsplit('/', 1)[-1]}")
        if merged and source.get("catalog_path")
        else _step(
            "merge",
            "Catalog merge",
            STEP_READY,
            "this drive is not folded into an all-years index",
            command="merge",
            values=values,
        )
    )

    pending = [step for step in steps if step["state"] in {STEP_READY, STEP_BLOCKED}]
    steps.append(
        _step("done", "Done", STEP_DONE, "cubes exist, catalogued and recalibratable")
        if not pending
        else _step("done", "Done", STEP_BLOCKED, f"waiting on {pending[0]['name']}")
    )
    return steps


# ---------------------------------------------------------------------------
# Page fragments
# ---------------------------------------------------------------------------


def _card(title: str, body: str, *, classes: str = "", identifier: str = "") -> str:
    attribute = f"card {classes}".strip()
    anchor = f' id="{_e(identifier)}"' if identifier else ""
    return (
        f'<article class="{attribute}"{anchor}>'
        f'<h3 class="card-title">{title}</h3>{body}</article>'
    )


def _pill(state: str, label: str) -> str:
    return f'<span class="pill {state}">{_e(label)}</span>'


def _rows_from(sources: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for source in sources:
        run = source.get("run") or {}
        for cube in source.get("cubes", []):
            rows.append(
                {
                    **cube,
                    "volume_label": source.get("volume_label", "unknown"),
                    "drive_id": source.get("drive_id", ""),
                    "available": bool(source.get("available")),
                    "run_state": str(run.get("state") or "unmeasured"),
                    "source_state": source_state(source),
                }
            )
    return rows


def _catalog_table(rows: list[dict[str, Any]]) -> str:
    body = []
    for row in rows:
        gate_class, gate_label = _gate_state(row.get("gate", GATE_UNRECORDED))
        find = " ".join(
            str(row.get(key, ""))
            for key in ("shot_number", "year", "snapshot_id", "volume_label", "path", "source_file")
        ).lower()
        # A row from a drive that did not answer says so where it is read, not
        # only on the drive's own card above.
        drive_note = (
            "" if row["available"] else " " + _pill(*_SOURCE_STATES["missing-drive"])
        )
        coverage = ""
        if row.get("wavelength_min_nm") is not None and row.get("wavelength_max_nm") is not None:
            coverage = f"{float(row['wavelength_min_nm']):.1f}–{float(row['wavelength_max_nm']):.1f}"
        body.append(
            "<tr class=\"cube-row\""
            f' data-year="{_e(row.get("year") or "")}"'
            f' data-epoch="{_e(row.get("snapshot_id") or "")}"'
            f' data-drive="{_e(row.get("volume_label"))}"'
            f' data-status="{_e(row["run_state"])}"'
            f' data-find="{_e(find)}">'
            f'<td>{_e(row.get("year") or "unknown")}</td>'
            f'<td>{_e(row.get("shot_number"))}</td>'
            f'<td>{_e(row.get("snapshot_id"))}</td>'
            f'<td>{_e(row.get("volume_label"))}{drive_note}</td>'
            f'<td>{_e(row["run_state"])}</td>'
            f"<td>{_pill(gate_class, gate_label)}</td>"
            f"<td>{_e(coverage or 'unmeasured')}</td>"
            "</tr>"
        )
    header = (
        "<tr><th>Year</th><th>Shot</th><th>Epoch</th><th>Drive</th><th>Run</th>"
        "<th>Gate</th><th>Coverage nm</th></tr>"
    )
    if not body:
        return '<p class="note state-unmeasured">No catalog rows were merged into this index.</p>'
    return (
        '<div class="scroll-x"><table id="cube-table"><thead>'
        f"{header}</thead><tbody>{''.join(body)}</tbody></table></div>"
    )


def _chip(text: str) -> str:
    return f'<span class="chip">{_e(text)}</span>'


def _path_chip(value: Any) -> str:
    """One truncated path. The whole path stays in the title and in the copy."""

    text = str(value or "")
    if not text:
        return _chip("root unknown")
    shown = text if len(text) <= 36 else "…" + text[-35:]
    return f'<span class="chip path" title="{_e(text)}">{_e(shown)}</span>'


def _drive_verdict(source: dict[str, Any], drift: list[dict[str, Any]]) -> tuple[str, str]:
    """Name this drive's wavelength-alignment verdict.

    A drive no evidence measures is unmeasured, which is its own state and not
    a quiet ``aligned``: nothing on this page turns silence into a pass.
    """

    evidence = _evidence_for_drive(source, drift)
    if evidence is None:
        return ("state-unmeasured", "alignment unmeasured")
    return _verdict_state(evidence["evidence"].get("verdict"))


def _drive_chips(source: dict[str, Any], drift: list[dict[str, Any]]) -> str:
    """Facts only: the states, the counts, and one truncated path.

    The wavelength-alignment verdict leads, because it is the one fact that
    decides whether this drive's data may be converted at all.
    """

    state = source_state(source)
    classes, label = _SOURCE_STATES[state]
    run = source.get("run") or {}
    counts = ", ".join(f"{value} {key}" for key, value in sorted((run.get("counts") or {}).items()))
    epochs = sorted({str(cube.get("snapshot_id") or "") for cube in source.get("cubes", [])} - {""})
    verdict_class, verdict_label = _drive_verdict(source, drift)
    chips = [
        _pill(f"{verdict_class} verdict-lead", verdict_label),
        _pill(classes, label),
        # A drive with no receipt has no gate to report; borrowing the pre-gate
        # word for it would claim a receipt that never existed.
        _pill(*_gate_state(run.get("gate", GATE_UNRECORDED))) if run else "",
        _chip(f"{len(source.get('cubes', []))} cube(s)"),
        _chip(f"id {source.get('drive_id') or 'not announced'}"),
        _chip(f"run {run.get('state')}" + (f" · {counts}" if counts else "")) if run else "",
        _chip(", ".join(epochs[:2]) + ("…" if len(epochs) > 2 else "")) if epochs else "",
        _path_chip(source.get("drive_root")),
    ]
    return f'<p class="chips">{"".join(chips)}</p>'


def _source_cards(sources: list[dict[str, Any]], drift: list[dict[str, Any]]) -> str:
    cards = []
    for source in sources:
        classes, _ = _SOURCE_STATES[source_state(source)]
        run = source.get("run") or {}
        lines = [_drive_chips(source, drift)]
        if run.get("drive_warning"):
            lines.append(f'<p class="warn">{_e(run["drive_warning"])}</p>')
        cards.append(
            _card(_e(source.get("volume_label", "unknown")), "".join(lines), classes=classes)
        )
    if not cards:
        return '<p class="note state-unmeasured">This index names no drives.</p>'
    return f'<div class="card-grid">{"".join(cards)}</div>'


# ---------------------------------------------------------------------------
# The stepper
# ---------------------------------------------------------------------------


def _step_box(step: dict[str, Any], position: int, primary: bool) -> str:
    state = _STEP_CLASSES[step["state"]]
    classes = f"step {state}" + (" is-primary" if primary else "")
    evidence = f'<span class="step-fact">{step["evidence"]}</span>' if step["evidence"] else ""
    current = ' aria-current="step"' if primary else ""
    return (
        f'<li class="{classes}" data-step="{_e(step["key"])}"{current}>'
        f'<span class="step-no">{position}</span>'
        f'<span class="step-name">{_e(step["name"])}</span>'
        f'{_pill(state, step["state"])}'
        f'<span class="step-fact">{_e(step["fact"])}</span>{evidence}</li>'
    )


def _next_block(step: dict[str, Any], identifier: str) -> str:
    head = f'<p class="next-label">Do this next — {_e(step["name"])}</p>'
    if step["shapes"]:
        return head + _command_row(
            identifier, step["name"], step["meaning"] or step["fact"], step["shapes"]
        )
    return head + f'<p class="note">{_e(step["fact"])}</p>'


def _flow(steps: list[dict[str, Any]], identifier: str, *, done_line: str) -> str:
    primary = _primary_step(steps)
    boxes = "".join(
        _step_box(step, position + 1, position == primary)
        for position, step in enumerate(steps)
    )
    if primary < 0:
        detail = f'<p class="next-label">{_e(done_line)}</p>'
    else:
        detail = _next_block(steps[primary], f"{identifier}-{steps[primary]['key']}")
    return f'<ol class="flow">{boxes}</ol><div class="next">{detail}</div>'


def _drive_row(
    source: dict[str, Any],
    steps: list[dict[str, Any]],
    position: int,
    drift: list[dict[str, Any]],
) -> str:
    return (
        f'<article class="drive-row" id="drive-{position}" '
        f'data-drive="{_e(source.get("volume_label", "unknown"))}">'
        f'<h3 class="drive-name">{_e(source.get("volume_label", "unknown"))}</h3>'
        f"{_drive_chips(source, drift)}"
        f'{_flow(steps, f"now-d{position}", done_line="This drive is done: cubes exist, catalogued and recalibratable.")}'
        "</article>"
    )


def _absent_line(source: dict[str, Any], steps: list[dict[str, Any]]) -> str:
    primary = _primary_step(steps)
    following = steps[primary]["name"] if primary >= 0 else "nothing"
    return (
        f'<p class="drive-absent">{_pill(*_SOURCE_STATES["missing-drive"])} '
        f'<strong>{_e(source.get("volume_label", "unknown"))}</strong> '
        f'<span class="muted">{len(source.get("cubes", []))} cube(s) remembered · '
        f"next: {_e(following)}</span></p>"
    )


def _simple_table(headers: tuple[str, ...], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{_e(name)}</th>" for name in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{_e(cell)}</td>" for cell in row) + "</tr>" for row in rows
    )
    return f'<div class="scroll-x"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _number(value: Any, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


def _per_shot_table(evidence: dict[str, Any]) -> str:
    per_shot = evidence.get("per_shot") or []
    if not per_shot:
        return (
            '<p class="note state-unmeasured">No per-shot residuals were resolved, so this '
            "evidence carries no time series.</p>"
        )
    rows = [
        [
            str(item.get("shot_number", "")),
            str(item.get("date", "") or "unknown"),
            str(item.get("lines", 0)),
            _number(item.get("median_shift_px")),
            _number(item.get("pixel_spread_px")),
            str(item.get("group", "one group")),
        ]
        for item in per_shot
    ]
    return _simple_table(
        ("Shot", "Date", "Lines", "Median Δx px", "Spread px", "Group"), rows
    )


def _order_corrections(evidence: dict[str, Any]) -> str:
    corrections = evidence.get("order_corrections") or []
    if not corrections:
        return ""
    rows = [
        [
            str(item.get("order")),
            _number(item.get("reference_pixel"), 1),
            _number(item.get("dispersion_nm_per_px"), 5),
            _number(item.get("predicted_shift_nm"), 4),
            "extrapolated" if item.get("extrapolated") else "within measured pixels",
        ]
        for item in corrections
    ]
    table = _simple_table(
        ("Order", "Reference px", "Dispersion nm/px", "Predicted Δλ nm", "Range"), rows
    )
    return f"<h4>Per-order correction the refinement would apply</h4>{table}"


def _thresholds(evidence: dict[str, Any]) -> str:
    thresholds = evidence.get("thresholds_px") or {}
    if not thresholds:
        return ""
    rows = [[key.replace("_", " "), str(value)] for key, value in sorted(thresholds.items())]
    return f"<h4>Thresholds this verdict was judged against</h4>{_simple_table(('Bound', 'Value'), rows)}"


def _quorum(evidence: dict[str, Any]) -> str:
    quorum = (evidence.get("summary") or {}).get("quorum") or {}
    if not quorum:
        return ""
    if quorum.get("satisfied"):
        return (
            f'<p class="muted">Quorum satisfied: {_e(quorum.get("resolved_lines"))} resolved '
            f'line(s) over {_e(quorum.get("distinct_orders"))} order(s), covering '
            f'{_e(_number((quorum.get("coverage_nm") or [0, 0])[0], 1))}–'
            f'{_e(_number((quorum.get("coverage_nm") or [0, 0])[1], 1))} nm.</p>'
        )
    return (
        '<p class="warn state-insufficient-data">Quorum failed — this interval could not carry '
        f'a verdict: {_e(quorum.get("reason") or "no reason recorded")}.</p>'
    )


def _skipped(evidence: dict[str, Any]) -> str:
    skipped = evidence.get("skipped_cubes") or []
    if not skipped:
        return '<p class="muted">No sampled cube was skipped.</p>'
    rows = [
        [str(item.get("shot_number", "")), str(item.get("cube", "")), str(item.get("reason", ""))]
        for item in skipped
    ]
    return f"<h4>Skipped cubes</h4>{_simple_table(('Shot', 'Cube', 'Reason'), rows)}"


def _command_block(identifier: str, shell: str, name: str, command: str) -> str:
    """One shell shape of one command: meaning above, machinery behind a toggle."""

    body_id = f"{identifier}-{shell}"
    return (
        f'<div class="cmd-shell"><div class="cmd-head"><span class="shell-name">{_e(name)}</span>'
        f'<button type="button" class="fold-toggle" aria-expanded="false" '
        f'aria-controls="{_e(body_id)}" data-show="show">show</button>'
        f'<button type="button" class="copy" data-copy="{_e(command)}">Copy</button></div>'
        f'<pre class="cmd-body" id="{_e(body_id)}" data-shell="{_e(shell)}" hidden>'
        f"{_e(command)}</pre></div>"
    )


def _command_row(identifier: str, title: str, meaning: str, shapes: dict[str, str]) -> str:
    blocks = "".join(
        _command_block(identifier, shell, name, shapes[shell])
        for shell, name in SHELL_NAMES
        if shapes.get(shell)
    )
    return (
        f'<article class="cmd" id="{_e(identifier)}">'
        f'<p class="cmd-meaning"><strong>{_e(title)}</strong> — {_e(meaning)}</p>'
        f"{blocks}</article>"
    )


def _composed_commands(values: dict[str, str]) -> str:
    rows = []
    for template in COMMAND_TEMPLATES:
        shapes = {shell: _fill(template[shell], values, shell) for shell, _ in SHELL_NAMES}
        rows.append(
            _command_row(str(template["id"]), str(template["title"]), str(template["meaning"]), shapes)
        )
    return "".join(rows)


def _repair_rows(evidence: dict[str, Any], identifier: str) -> str:
    """Render the stored repair sequence as meaning-first command rows."""

    steps = evidence.get("repair_commands") or []
    if not steps:
        return ""
    grouped: list[dict[str, Any]] = []
    for step in steps:
        purpose = str(step.get("purpose", ""))
        if grouped and grouped[-1]["purpose"] == purpose:
            grouped[-1]["shapes"][str(step.get("shell", "any"))] = str(step.get("command", ""))
            continue
        grouped.append(
            {
                "purpose": purpose,
                "shapes": {str(step.get("shell", "any")): str(step.get("command", ""))},
            }
        )
    rows = []
    for position, item in enumerate(grouped, start=1):
        shapes = dict(item["shapes"])
        shared = shapes.pop("any", None)
        if shared is not None:
            shapes = {"powershell": shared, "posix": shared}
        title = f"Step {position}"
        if not any(shapes.values()):
            rows.append(
                f'<article class="cmd" id="{_e(identifier)}-repair-{position}">'
                f'<p class="cmd-meaning"><strong>{_e(title)}</strong> — {_e(item["purpose"])}</p>'
                '<p class="note">No command: this step is a decision you make and record '
                "yourself.</p></article>"
            )
            continue
        rows.append(
            _command_row(f"{identifier}-repair-{position}", title, str(item["purpose"]), shapes)
        )
    return f"<h4>Repair sequence composed from this evidence</h4>{''.join(rows)}"


def _line_drilldown(evidence: dict[str, Any], identifier: str) -> str:
    lines = evidence.get("lines") or []
    rows = [
        [
            str(item.get("shot_number", "") or item.get("cube", "")),
            str(item.get("line", "")),
            str(item.get("status", "")),
            _number(item.get("expected_nm"), 4),
            _number(item.get("residual_nm"), 4),
            _number(item.get("pixel_residual_px"), 3),
            str(item.get("reason", "")),
        ]
        for item in lines
    ]
    table = _simple_table(
        ("Shot", "Line", "Status", "Expected nm", "Residual nm", "Δx px", "Reason"), rows
    )
    return (
        f'<details class="drill" id="{_e(identifier)}-lines">'
        f"<summary>Per-line evidence ({len(rows)} row(s)) — Escape closes it</summary>"
        f'<div class="drill-body">{table}'
        '<p class="drill-exit"><button type="button" class="drill-close">Close per-line '
        "evidence</button></p></div></details>"
    )


def _drift_card(entry: dict[str, Any], position: int) -> str:
    evidence = entry["evidence"]
    identifier = f"drift-{position}"
    state, label = _verdict_state(evidence.get("verdict"))
    summary = evidence.get("summary") or {}
    head = [
        f'<p class="state-line">{_pill(state, label)} '
        f'<span class="muted">{_e(len(evidence.get("sampled_cubes") or []))} cube(s) sampled, '
        f'{_e(len(evidence.get("skipped_cubes") or []))} skipped</span></p>',
        f'<p class="muted">epochs {_e(", ".join(evidence.get("snapshot_ids") or []) or "unknown")} '
        f'· written {_e(evidence.get("created_at") or "unknown")} · {_e(entry["path"])}</p>',
    ]
    if evidence.get("interval_warning"):
        head.append(
            f'<p class="warn interval-warning"><strong>Interval warning</strong> — '
            f'{_e(evidence["interval_warning"])}</p>'
        )
    if evidence.get("data_requirement"):
        head.append(f'<p class="note">{_e(evidence["data_requirement"])}</p>')
    if summary:
        head.append(
            '<p class="muted">median shift '
            f'{_e(_number(summary.get("median_shift_px")))} px · largest residual '
            f'{_e(_number(summary.get("maximum_absolute_pixel_residual_px")))} px · spread '
            f'{_e(_number(summary.get("maximum_pixel_deviation_px")))} px · tolerance '
            f'{_e(_number(summary.get("alignment_tolerance_px")))} px</p>'
        )
    body = "".join(
        [
            *head,
            _quorum(evidence),
            "<h4>Per-shot residuals</h4>",
            _per_shot_table(evidence),
            _order_corrections(evidence),
            _thresholds(evidence),
            _skipped(evidence),
            _repair_rows(evidence, identifier),
            _line_drilldown(evidence, identifier),
        ]
    )
    return _card(_e(label), body, classes=f"verdict {state}", identifier=identifier)


def _drift_section(drift: list[dict[str, Any]]) -> str:
    if not drift:
        return (
            '<p class="note state-unmeasured">No drift evidence was supplied to this build. '
            "Bulk readiness is unmeasured — which is not the same as aligned.</p>"
        )
    return "".join(_drift_card(entry, position) for position, entry in enumerate(drift, start=1))


def _bounds(start: Any, end: Any) -> str:
    if start and end:
        return f"{start}–{end}"
    if start:
        return f"from {start}"
    if end:
        return f"to {end}"
    return "—"


def _epoch_table(context: dict[str, Any]) -> str:
    """The registry's epochs, and any epoch the cubes name that it does not."""

    counts: dict[str, int] = {}
    for row in context["rows"]:
        key = str(row.get("snapshot_id") or "unassigned")
        counts[key] = counts.get(key, 0) + 1
    epoch_rows = list(context["registry"].get("epoch_rows") or [])
    named = {str(row["snapshot_id"]) for row in epoch_rows}
    table_rows = [
        [
            str(row["snapshot_id"]),
            _bounds(row.get("shot_from"), row.get("shot_to")),
            _bounds(row.get("date_from"), row.get("date_to")),
            str(counts.get(str(row["snapshot_id"]), 0)),
        ]
        for row in epoch_rows
    ]
    table_rows.extend(
        [key, "not in the registry", "not in the registry", str(value)]
        for key, value in sorted(counts.items())
        if key not in named
    )
    if not table_rows:
        return (
            '<p class="note state-unmeasured">No epoch is named by this build: no registry was '
            "read and no cube claims a snapshot.</p>"
        )
    return _simple_table(("Epoch", "Shots", "Dates", "Cubes here"), table_rows)


def _document_sections(documents: list[dict[str, str]]) -> str:
    parts = []
    for document in documents:
        parts.append(
            f'<article class="card doc" id="{_e(document["anchor"])}">'
            f'<h3 class="card-title">{_e(document["title"])}</h3>'
            f'<p class="muted">Rendered from {_e(document["origin"])}</p>'
            f'<div class="doc-body">{document["html"]}</div></article>'
        )
    return "".join(parts)


# ---------------------------------------------------------------------------
# Rails
# ---------------------------------------------------------------------------


def _select(identifier: str, label: str, options: list[tuple[str, str]], *, first: str) -> str:
    rendered = "".join(
        f'<option value="{_e(value)}">{_e(text)}</option>' for value, text in options
    )
    return (
        f'<label class="field"><span>{_e(label)}</span>'
        f'<select id="{_e(identifier)}"><option value="">{_e(first)}</option>'
        f"{rendered}</select></label>"
    )


def _text_field(
    identifier: str,
    label: str,
    value: str,
    *,
    note: str = "",
    placeholder: str = "",
    browse: bool = False,
) -> str:
    """One text control, and — served only — the Browse button beside it.

    A field with nothing to fill carries a placeholder and no ``value`` at all:
    a baked example path is a claim about a machine this build never saw. The
    Browse button keeps one position whether the picker has ever been opened
    or not, and falls through to the page's base ``button`` styling, so it
    looks pressable at rest rather than only under the pointer.
    """

    hint = f'<small class="muted">{_e(note)}</small>' if note else ""
    attributes = f' value="{_e(value)}"' if value else ""
    if placeholder:
        attributes += f' placeholder="{_e(placeholder)}"'
    control = f'<input type="text" id="{_e(identifier)}"{attributes}>'
    if browse:
        control = (
            f'<span class="field-row">{control}'
            f'<button type="button" class="browse" data-browse="{_e(identifier)}">'
            "Browse…</button></span>"
        )
    return (
        f'<label class="field"><span>{_e(label)}</span>'
        f"{control}{hint}</label>"
    )


def _filter_card(rows: list[dict[str, Any]], sources: list[dict[str, Any]]) -> str:
    def _values(key: str) -> list[tuple[str, str]]:
        return [(item, item) for item in sorted({str(row.get(key) or "") for row in rows} - {""})]

    statuses = sorted({str((source.get("run") or {}).get("state") or "unmeasured") for source in sources})
    body = (
        '<div class="fields">'
        + _select("filter-year", "Year", _values("year"), first="All years")
        + _select("filter-epoch", "Epoch", _values("snapshot_id"), first="All epochs")
        + _select("filter-drive", "Drive", _values("volume_label"), first="All drives")
        + _select(
            "filter-status",
            "Run status",
            [(status, status) for status in statuses],
            first="All run states",
        )
        + "</div>"
    )
    return _card("Filter the catalog", body)


def _epoch_options(epochs: list[str], registry: dict[str, Any]) -> str:
    """The calibrations this build can name — or the one honest reason it cannot.

    The registry maps shots to one immutable snapshot by date, so choosing a
    calibration here is choosing it "for a drive or dates" without saying so
    twice.
    """

    # "unassigned" is a cube's way of saying it had no calibration identity;
    # it is not a calibration anybody can choose.
    known = [item for item in epochs if item and item != "unassigned"]
    saved = [
        item
        for item in (registry.get("saved") or [])
        if item and item not in known
    ]
    options = "".join(f'<option value="{_e(item)}">{_e(item)}</option>' for item in known)
    # A saved snapshot is real without a registry; it is only not yet
    # registered, and the label says exactly that much and no more.
    options += "".join(
        f'<option value="{_e(item)}">{_e(item)} — saved, not in any registry</option>'
        for item in saved
    )
    if options:
        return options
    stated = {
        "not supplied": "no registry supplied to this build",
        "unreadable": "registry unreadable — the commands name it unread",
    }.get(str(registry.get("status")), "this registry names no epoch")
    return f'<option value="">{_e(stated)}</option>'


def _composer_card(
    values: dict[str, str],
    epochs: list[str],
    registry: dict[str, Any],
    *,
    served: bool = False,
) -> str:
    """Two questions, and every answer they decide folded away behind them.

    The data folder and the calibration are the only two facts a person holds
    that this build cannot read off a file; everything else is derived from
    them and stays editable in the fold for the run that needs an exception.
    Served, the data folder also carries the Browse button into the picker.
    """

    body = (
        '<p class="muted">Point at the shots and the calibration; the rest is derived.</p>'
        '<div class="fields">'
        + _text_field(
            "f-input",
            "Data folder (this drive's SIF shots)",
            "",
            placeholder="the folder holding this drive's SIF shots",
            browse=served,
        )
        + '<label class="field"><span>Calibration</span>'
        f'<select id="f-epoch">{_epoch_options(epochs, registry)}</select></label>'
        + "</div>"
        '<details class="fold-group" id="composer-advanced">'
        "<summary>Advanced — every derived value, editable</summary>"
        '<div class="fields">'
        + _text_field(
            "f-output", "Cubes folder", "", placeholder="the data folder, unless you say otherwise"
        )
        + _text_field("f-label", "Volume label", "", placeholder="the data folder's own name")
        + _text_field(
            "f-verdict",
            "Drift evidence to write",
            "",
            placeholder=f"{values['verdict'].rsplit('/', 1)[-1]} in the data folder",
        )
        + _text_field("f-registry", "Registry", values["registry"], browse=served)
        + _text_field("f-calibrations", "Snapshot root", values["calibrations"], browse=served)
        + _text_field("f-catalog", "Merged catalog", values["catalog"])
        + _text_field("f-plan", "Plan file to save", values["plan"])
        + _text_field("f-pattern", "SIF pattern", values["pattern"])
        + "</div></details>"
    )
    return _card("Compose a plan and commands", body, classes="rail-card--growing")


def _scope_card(rows: list[dict[str, Any]]) -> str:
    body = (
        f'<p id="scope-line" class="scope-line">Showing {len(rows)} of {len(rows)} cube(s).</p>'
        '<label class="field"><span>Find in the catalog</span>'
        '<input type="search" id="find" placeholder="shot, epoch, drive, file"></label>'
        '<p class="actions"><button type="button" id="reset">Reset filters and Find</button></p>'
    )
    return _card("Current scope", body)


def _send_to_composer_card(drives: list[dict[str, str]]) -> str:
    """The Drives tab's own composer control: one press, one drive, Now leads it.

    The full composer lives in the Now rail because that is where the ready
    steps' commands are read.  This card is the compact manual entry from the
    catalog: choosing a drive here and pressing the button switches Now to that
    drive's composed commands rather than duplicating the composer's fields.
    """

    options = "".join(
        f'<option value="{_e(drive["root"])}">{_e(drive["label"])}</option>' for drive in drives
    )
    body = (
        '<label class="field"><span>Drive</span>'
        f'<select id="send-drive">{options}</select></label>'
        '<p class="actions"><button type="button" id="send-compose">Compose in Now</button></p>'
    )
    return _card("Compose for one drive", body)


def _position_card(context: dict[str, Any]) -> str:
    body = (
        f'<p>{_e(context["connected_count"])} drive(s) connected, '
        f'{_e(context["absent_count"])} remembered.</p>'
        f'<p>{_e(context["ready_count"])} step(s) ready, '
        f'{_e(context["blocked_count"])} blocked.</p>'
        f'<p class="muted">Built {_e(context["generated_at"])}</p>'
    )
    return _card("Campaign position", body)


def _context_card(context: dict[str, Any]) -> str:
    body = (
        f'<p>{_e(context["source_count"])} drive(s), {_e(context["cube_count"])} cube(s).</p>'
        f'<p class="chips">{_path_chip(context["catalog_path"])}</p>'
        f'<p class="muted">Built {_e(context["generated_at"])}</p>'
    )
    return _card("Context", body)


def _registry_card(context: dict[str, Any]) -> str:
    registry = context["registry"]
    if registry["status"] == "read":
        registry_line = (
            f'<p class="state-line">{_pill("gate-verdict", "registry read")} '
            f'<span class="muted">{_e(len(registry["epochs"]))} epoch(s)</span></p>'
            f'<p class="muted">{_e(registry["path"])}</p>'
        )
    elif registry["status"] == "unreadable":
        registry_line = (
            f'<p class="state-line">{_pill("state-unrecognized", "registry unreadable")}</p>'
            f'<p class="note">{_e(registry["detail"])}</p>'
        )
    else:
        registry_line = (
            f'<p class="state-line">{_pill("state-unmeasured", "no registry supplied")}</p>'
            '<p class="note">This build was given no registry, so the page cannot say which '
            "epoch a run would select.</p>"
        )
    body = (
        f'<p>{_e(len(context["epochs"]))} epoch(s), '
        f'{_e(context["drift_count"])} drift verdict(s).</p>'
        f"{registry_line}"
    )
    return _card("Calibration context", body)


def _legend(title: str, entries: list[tuple[str, str, str]]) -> str:
    """One tab's legend — the only place on the page where a state is taught.

    Folded closed by default: the teaching is on the page exactly once, but it
    sits quietly behind its header until asked (owner correction 2026-08-18,
    "the UI should teach me but not scream").
    """

    rows = "".join(
        f'<li>{_pill(state, label)} <span class="legend-note">{_e(meaning)}</span></li>'
        for state, label, meaning in entries
    )
    return _card(
        title,
        '<details class="fold-group legend-fold"><summary>the words, when you want them'
        f'</summary><ul class="legend">{rows}</ul></details>',
    )


def _step_legend_card() -> str:
    return _legend(
        "What a step state means",
        [
            (_STEP_CLASSES[STEP_DONE], STEP_DONE, "a file records it; the evidence is linked"),
            (_STEP_CLASSES[STEP_READY], STEP_READY, "the composed command for exactly this step"),
            (_STEP_CLASSES[STEP_BLOCKED], STEP_BLOCKED, "what is missing, and the step that supplies it"),
            (
                _STEP_CLASSES[STEP_UNRECORDED],
                STEP_UNRECORDED,
                "no file on this page can say whether it was done",
            ),
        ],
    )


def _drive_legend_card() -> str:
    return _legend(
        "What the drive states mean",
        [
            ("state-missing-drive", "missing drive", _SOURCE_NOTES["missing-drive"]),
            ("state-unmeasured", "unmeasured", _SOURCE_NOTES["unmeasured"]),
            ("state-empty", "empty", _SOURCE_NOTES["empty"]),
            ("state-measured", "measured", _SOURCE_NOTES["measured"]),
            (
                "gate-verdict",
                "verdict-authorized",
                "the run was authorized by drift evidence it names.",
            ),
            (
                "gate-sample",
                "unverified sample",
                "the legal first registry run; its cubes await a verdict.",
            ),
        ],
    )


def _verdict_legend_card() -> str:
    return _legend(
        "What a verdict word means",
        [
            ("verdict-aligned", "aligned", "the sampled residuals sit inside the tolerance."),
            ("verdict-shifted", "shifted", "one rigid detector shift, repairable by refinement."),
            (
                "verdict-beyond-repair",
                "misaligned-beyond-repair",
                "past the repair limit; the raw SIF data is needed.",
            ),
            (
                "state-insufficient-data",
                "insufficient-data",
                "a judged verdict, never aligned: the sample could not carry one.",
            ),
            (
                "state-unmeasured",
                "alignment unmeasured",
                "no evidence in this build measures that drive — never a quiet aligned.",
            ),
            (
                "state-unrecognized",
                "unrecognized",
                "a word this page does not know, never dressed as one it does.",
            ),
        ],
    )


def _index_card(entries: list[tuple[str, str]], identifier: str) -> str:
    links = "".join(
        f'<a href="#{_e(anchor)}" class="sectnav-link">{_e(title)}</a>' for anchor, title in entries
    )
    return _card(
        "On this page",
        f'<nav class="sectnav" id="{_e(identifier)}">{links}</nav>',
        classes="rail-card--growing",
    )


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

_CSS = """
*, *::before, *::after { box-sizing: border-box; }
:root {
  color-scheme: light dark;
  /* The topbar's own height. Under border-box its 1px bottom border is inside
     this number, so the rails' sticky sum adds the wrapper padding and nothing
     else -- adding the border again would pin every rail 1px low. */
  --bar: 56px;
  --gap: 20px;
  --rail-left: 21rem;
  --rail-right: 18rem;
  --bg: #f5f3ee;
  --panel: #fffefa;
  --raised: #efece4;
  --ink: #232019;
  --muted: #66605a;
  --line: #d5cfc3;
  --accent: #2f5d63;
  --miss: #7b4ea8;
  --unmeasured: #5a6b7b;
  --empty: #7a7264;
  --judged: #a86a12;
  --bad: #a8382f;
  --good: #2f6b46;
  --shift: #55519c;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #191b1e;
    --panel: #212429;
    --raised: #272b31;
    --ink: #e7e3da;
    --muted: #9d968c;
    --line: #363b42;
    --accent: #79c0c4;
    --miss: #c3a1e8;
    --unmeasured: #9fb3c6;
    --empty: #b6ab97;
    --judged: #e0a95a;
    --bad: #e4897f;
    --good: #79c295;
    --shift: #a7a2ea;
  }
}
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif;
}
a { color: var(--accent); }
h1, h2, h3, h4 { line-height: 1.25; margin: 0 0 .4rem; }
h1 { font-size: 1.05rem; font-weight: 600; }
h2 { font-size: 1.15rem; }
h3 { font-size: .98rem; }
h4 { font-size: .88rem; margin-top: 1rem; color: var(--muted); text-transform: uppercase;
     letter-spacing: .04em; }
p { margin: .4rem 0; }
.topbar {
  position: sticky; top: 0; z-index: 5;
  height: var(--bar);
  display: flex; align-items: center; gap: .8rem;
  padding: 0 24px;
  background: var(--panel);
  border-bottom: 1px solid var(--line);
}
.topbar .tagline { color: var(--muted); font-size: .82rem; margin-left: auto; }
.tabs { display: flex; gap: .3rem; }
.tab { padding: .3rem .8rem; font-size: .95rem; }
.tab[aria-selected="true"] { border-color: var(--accent); background: var(--bg); font-weight: 600; }
.wrap { padding: var(--gap) 24px 80px; }
/* The rail grid stretches by default so each rail's row spans the content
   column's full height; only the rails themselves pin to the top. Putting
   align-items:start on this parent instead would collapse the row and give a
   sticky rail zero travel -- correct-looking CSS, rail rides the page down. */
.rail-grid {
  display: grid;
  grid-template-columns: var(--rail-left) minmax(0, 1fr) var(--rail-right);
  gap: var(--gap);
  align-items: stretch;
  max-width: 1720px;
  margin: 0 auto;
}
.rail-grid > * { min-width: 0; }
.rail {
  position: sticky;
  top: calc(var(--bar) + var(--gap));
  align-self: start;
  display: flex;
  flex-direction: column;
  gap: 14px;
  max-height: calc(100vh - var(--bar) - var(--gap) - 20px);
  overflow: auto;
}
/* A growing rail card spends only the space the rail can afford, and owns its
   own scroll, so the rail's own overflow stays a backstop rather than the
   working mechanism: the rail never becomes a second page scrollbar. */
[hidden] { display: none !important; }
/* Each rail holds one group per tab and shows only the active tab's group, so
   a rail never carries a neighbouring tab's cargo. */
.rail-group { display: flex; flex: 1 1 auto; flex-direction: column; gap: 14px; min-height: 0; }
/* A tab whose controls are empty gives its column back rather than parking an
   empty rail beside the content. */
body.no-left .rail-grid { grid-template-columns: minmax(0, 1fr) var(--rail-right); }
.rail-card--growing { display: flex; flex: 1 1 auto; flex-direction: column; min-height: 0; }
.rail-card--growing > .sectnav { flex: 1 1 auto; min-height: 0; overflow-y: auto;
  overscroll-behavior: contain; }
/* An open fold spends only the space the rail can afford and scrolls inside
   itself, so opening it never grows the card past the rail. */
.rail-card--growing > .fold-group[open] { flex: 1 1 auto; min-height: 0; overflow-y: auto;
  overscroll-behavior: contain; }
.content { display: flex; flex-direction: column; gap: var(--gap); }
section.panel {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: .55rem;
  padding: 1rem 1.1rem 1.2rem;
}
.card {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: .5rem;
  padding: .7rem .8rem;
}
.rail .card { background: var(--raised); }
.card-title { font-size: .78rem; text-transform: uppercase; letter-spacing: .06em;
  color: var(--muted); }
.card-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(15rem, 1fr));
  gap: .7rem; margin: .6rem 0 1rem; }
/* Data reads at the page's own size. Small type is for genuinely secondary
   metadata -- build stamps, shell names, card headings -- and nothing else. */
.muted { color: var(--muted); font-size: .92rem; }
.note { border-left: .2rem solid var(--line); padding-left: .6rem; color: var(--muted); }
.chips { display: flex; flex-wrap: wrap; gap: .3rem; align-items: center; margin: .35rem 0; }
.chip { display: inline-block; border: 1px solid var(--line); border-radius: .3rem;
  background: var(--raised); padding: 0 .4rem; font-size: .9rem; white-space: nowrap; }
.chip.path { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; max-width: 100%;
  overflow: hidden; text-overflow: ellipsis; }
.legend-note { color: var(--muted); font-size: .9rem; }
.warn { border-left: .3rem solid var(--judged); padding: .4rem .6rem; background: var(--raised);
  border-radius: .3rem; }
.pill { display: inline-block; border: 1px solid currentColor; border-radius: 999px;
  padding: 0 .5rem; font-size: .76rem; white-space: nowrap; }
/* The wavelength-alignment verdict is the loudest fact on a drive: it leads the
   chip row and is read before anything else on the card. */
.pill.verdict-lead { font-size: .98rem; font-weight: 700; border-width: 2px;
  padding: .05rem .6rem; }
.state-missing-drive { color: var(--miss); }
.state-unmeasured { color: var(--unmeasured); }
.state-empty { color: var(--empty); }
.state-insufficient-data { color: var(--judged); }
.state-unrecognized { color: var(--bad); }
.state-measured { color: var(--good); }
.gate-verdict { color: var(--good); }
.gate-sample { color: var(--judged); border-style: dashed; }
.gate-ungated { color: var(--muted); }
.gate-unrecorded { color: var(--unmeasured); border-style: dotted; }
.gate-unrecognized { color: var(--bad); }
.verdict-aligned { color: var(--good); }
.verdict-shifted { color: var(--shift); }
.verdict-beyond-repair { color: var(--bad); }
.card.state-missing-drive { border-left: .35rem solid var(--miss); }
.card.state-unmeasured { border-left: .35rem dashed var(--unmeasured); }
.card.state-empty { border-left: .35rem dotted var(--empty); }
.card.state-insufficient-data { border-left: .35rem solid var(--judged); }
.card.state-unrecognized { border-left: .35rem double var(--bad); }
.card.state-measured { border-left: .35rem solid var(--good); }
.card.verdict-aligned { border-left: .35rem solid var(--good); }
.card.verdict-shifted { border-left: .35rem solid var(--shift); }
.card.verdict-beyond-repair { border-left: .35rem solid var(--bad); }
.card.verdict { margin-bottom: 1rem; color: var(--ink); }
/* The stepper: one horizontal row of boxes per drive, drawn with borders and
   arrows in CSS -- no library, nothing fetched, and it scrolls inside itself
   rather than widening the page. */
.drive-row { border: 1px solid var(--line); border-radius: .5rem; padding: .6rem .7rem;
  margin-bottom: .9rem; background: var(--panel); }
.drive-name { font-size: 1.02rem; }
.drive-absent { display: flex; flex-wrap: wrap; gap: .4rem; align-items: baseline;
  border-top: 1px solid var(--line); padding: .35rem 0; margin: 0; }
/* The row wraps rather than scrolling: the whole flow has to be readable at a
   glance, and a stepper hidden past a horizontal scrollbar is not. */
.flow { display: flex; flex-wrap: wrap; gap: .9rem; list-style: none; margin: .5rem 0 0;
  padding: 0 0 .4rem .9rem; }
.step { position: relative; flex: 1 1 9rem; min-width: 9rem; border: 1px solid var(--line);
  border-radius: .4rem; padding: .4rem .5rem; background: var(--bg); display: flex;
  flex-direction: column; gap: .2rem; }
.step + .step::before { content: "\\2192"; position: absolute; left: -.72rem; top: 50%;
  transform: translateY(-50%); color: var(--muted); }
.step-no { font-size: .78rem; color: var(--muted); }
.step-name { font-weight: 600; }
.step .pill { align-self: flex-start; }
.step-fact { color: var(--muted); font-size: .9rem; overflow-wrap: anywhere; }
.step-done { border-color: var(--good); }
.step-ready { border-color: var(--accent); }
.step-blocked { border-color: var(--judged); border-style: dashed; }
.step-unrecorded { border-color: var(--unmeasured); border-style: dotted; }
.step.is-primary { border-width: 2px; background: var(--raised); }
.pill.step-done { color: var(--good); }
.pill.step-ready { color: var(--accent); }
.pill.step-blocked { color: var(--judged); }
.pill.step-unrecorded { color: var(--unmeasured); }
.next { border-top: 1px solid var(--line); margin-top: .5rem; }
.next-label { font-weight: 600; margin: .5rem 0 0; }
.next .cmd { border-top: none; padding-top: .2rem; }
.scroll-x { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 1rem; }
th, td { border-bottom: 1px solid var(--line); padding: .32rem .45rem; text-align: left;
  vertical-align: top; }
th { color: var(--muted); font-weight: 600; }
/* Every control keeps a border and a fill at rest: an unpressed button must
   never be indistinguishable from plain text. */
button, select, input, textarea {
  font: inherit;
  color: inherit;
  background: var(--raised);
  border: 1px solid var(--line);
  border-radius: .3rem;
  padding: .3rem .55rem;
}
button { cursor: pointer; }
button:hover { border-color: var(--accent); }
:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
.fields { display: flex; flex-direction: column; gap: .5rem; }
.field { display: grid; gap: .2rem; font-size: .8rem; color: var(--muted); }
.field input, .field select { width: 100%; color: var(--ink); }
/* A headed group inside a card, never a second card: one rule above it, a
   pressable head, and the derived fields flat underneath. */
.fold-group { border-top: 1px solid var(--line); margin-top: .6rem; }
.fold-group > summary { cursor: pointer; padding: .4rem .45rem; margin: .35rem 0 0;
  border: 1px solid var(--line); border-radius: .3rem; background: var(--raised);
  font-size: .78rem; text-transform: uppercase; letter-spacing: .05em; color: var(--muted);
  position: sticky; top: 0; }
.fold-group > summary:hover { border-color: var(--accent); }
.fold-group[open] > summary { margin-bottom: .5rem; }
.actions { display: flex; gap: .5rem; flex-wrap: wrap; }
.sectnav { display: flex; flex-direction: column; gap: .15rem; }
.sectnav-link { display: block; padding: .18rem .35rem; border-radius: .25rem;
  text-decoration: none; font-size: .86rem; }
.sectnav-link:hover { background: var(--bg); }
.sectnav-link.current { background: var(--bg); font-weight: 600; }
.legend { list-style: none; margin: 0; padding: 0; display: flex; flex-direction: column;
  gap: .3rem; }
.cmd { border-top: 1px solid var(--line); padding: .7rem 0 .3rem; }
.cmd-meaning { margin: 0 0 .45rem; }
.cmd-shell { margin: .3rem 0; }
.cmd-head { display: flex; align-items: center; gap: .5rem; }
.shell-name { font-size: .76rem; color: var(--muted); text-transform: uppercase;
  letter-spacing: .05em; min-width: 7.5rem; }
.cmd-body { margin: .35rem 0 0; padding: .5rem .6rem; background: var(--raised);
  border: 1px solid var(--line); border-radius: .3rem; white-space: pre-wrap;
  overflow-wrap: anywhere; font-size: .84rem; }
.plan-out { width: 100%; min-height: 12rem; white-space: pre; overflow-wrap: normal;
  overflow-x: auto; }
.drill { border: 1px solid var(--line); border-radius: .4rem; margin-top: .8rem; }
.drill > summary { position: sticky; top: calc(var(--bar) + 4px); background: var(--panel);
  padding: .45rem .6rem; cursor: pointer; border-radius: .4rem; }
.drill-body { padding: 0 .6rem .6rem; }
.drill-exit { position: sticky; bottom: 0; background: var(--panel); padding: .4rem 0; }
.doc-body { font-size: .93rem; }
.doc-body pre { background: var(--raised); border: 1px solid var(--line); border-radius: .3rem;
  padding: .5rem .6rem; overflow-x: auto; }
.doc-body code { background: var(--raised); border-radius: .2rem; padding: 0 .2rem; }
.doc-body pre code { background: none; padding: 0; }
.cube-row.hidden-row { display: none; }
/* One rule covers every jump target on the page, because every jump target's
   id starts with sec-. */
[id^="sec-"] { scroll-margin-top: calc(var(--bar) + 16px); }
@media (max-width: 900px) {
  .rail-grid { grid-template-columns: minmax(0, 1fr); }
  .rail { position: static; max-height: none; overflow: visible; }
  .rail-left { order: 1; }
  .rail-right { order: 2; }
  .content { order: 3; }
}
"""

_JS = """
var state = { tab: 'now' };

function byId(id) { return document.getElementById(id); }

function shellPath(shell, value) {
  var text = String(value === undefined || value === null ? '' : value);
  return shell === 'powershell'
    ? text.split('/').join('\\\\')
    : text.split('\\\\').join('/');
}

function fill(template, values, shell) {
  var filled = template;
  Object.keys(values).forEach(function (key) {
    filled = filled.split('{{' + key + '}}').join(shellPath(shell, values[key]));
  });
  return filled;
}

/* Which derived field the operator has taken over. A value typed in the
   Advanced fold is theirs; the data folder never overwrites it again. */
var derivedEdits = { output: false, label: false, verdict: false };
var DERIVED_FIELDS = { output: 'f-output', label: 'f-label', verdict: 'f-verdict' };

function folderPath(text) {
  var trimmed = String(text === undefined || text === null ? '' : text)
    .split('\\\\').join('/');
  while (trimmed.length > 1 && trimmed.charAt(trimmed.length - 1) === '/') {
    trimmed = trimmed.slice(0, -1);
  }
  return trimmed;
}

/* The other half of _derived_from_folder in reading_room.py: one data folder
   decides these five values, and the two halves must answer the same. */
function derivedFrom(folder) {
  var marked = folder ? folder : DATA.unfilled;
  var trimmed = folderPath(marked);
  return {
    input: marked,
    output: marked,
    cubes: marked,
    label: trimmed.split('/').pop() || marked,
    verdict: trimmed + '/' + DATA.evidence_name
  };
}

function applyDerived() {
  var control = byId('f-input');
  var folder = control ? control.value.trim() : '';
  var derived = derivedFrom(folder);
  Object.keys(DERIVED_FIELDS).forEach(function (key) {
    var field = byId(DERIVED_FIELDS[key]);
    if (!field || derivedEdits[key]) { return; }
    /* No data folder, nothing derived: the field goes back to its placeholder
       rather than showing the marker as if it were an answer. */
    field.value = folder ? derived[key] : '';
  });
  return derived;
}

function composerValues() {
  var values = {};
  Object.keys(DATA.values).forEach(function (key) { values[key] = DATA.values[key]; });
  var derived = applyDerived();
  Object.keys(derived).forEach(function (key) { values[key] = derived[key]; });
  var fields = {
    output: 'f-output', label: 'f-label', verdict: 'f-verdict',
    registry: 'f-registry', calibrations: 'f-calibrations', catalog: 'f-catalog',
    plan: 'f-plan', pattern: 'f-pattern', epoch: 'f-epoch'
  };
  Object.keys(fields).forEach(function (key) {
    var control = byId(fields[key]);
    if (control && control.value !== '') { values[key] = control.value; }
  });
  values.cubes = values.output;
  return values;
}

function compose() {
  var values = composerValues();
  var unfilled = byId('unfilled-note');
  if (unfilled) { unfilled.hidden = values.input !== DATA.unfilled; }
  var plan = byId('plan-out');
  if (plan) { plan.value = fill(DATA.plan_template, values, 'posix'); }
  DATA.commands.forEach(function (command) {
    ['powershell', 'posix'].forEach(function (shell) {
      var body = byId(command.id + '-' + shell);
      if (!body) { return; }
      var text = fill(command[shell], values, shell);
      body.textContent = text;
      var copy = body.parentNode.querySelector('.copy');
      if (copy) { copy.setAttribute('data-copy', text); }
    });
  });
}

function filterCatalog() {
  var year = byId('filter-year').value;
  var epoch = byId('filter-epoch').value;
  var drive = byId('filter-drive').value;
  var status = byId('filter-status').value;
  var needle = byId('find').value.trim().toLowerCase();
  var rows = document.querySelectorAll('#cube-table tbody tr');
  var shown = 0;
  Array.prototype.forEach.call(rows, function (row) {
    var keep =
      (!year || row.getAttribute('data-year') === year) &&
      (!epoch || row.getAttribute('data-epoch') === epoch) &&
      (!drive || row.getAttribute('data-drive') === drive) &&
      (!status || row.getAttribute('data-status') === status) &&
      (!needle || row.getAttribute('data-find').indexOf(needle) !== -1);
    row.classList.toggle('hidden-row', !keep);
    if (keep) { shown += 1; }
  });
  var scope = byId('scope-line');
  if (scope) {
    scope.textContent = 'Showing ' + shown + ' of ' + rows.length + ' cube(s).';
  }
}

function closeFolds(root) {
  var scope = root || document;
  var closed = false;
  Array.prototype.forEach.call(scope.querySelectorAll('details[open]'), function (item) {
    item.open = false;
    closed = true;
  });
  Array.prototype.forEach.call(scope.querySelectorAll('.fold-toggle'), function (button) {
    if (button.getAttribute('aria-expanded') !== 'true') { return; }
    var body = byId(button.getAttribute('aria-controls'));
    if (body) { body.hidden = true; }
    button.setAttribute('aria-expanded', 'false');
    button.textContent = button.getAttribute('data-show');
    closed = true;
  });
  return closed;
}

function copyPayload(button) {
  var text = button.getAttribute('data-copy') || '';
  var done = function () { button.textContent = 'Copied'; };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done, function () { legacyCopy(text, done); });
  } else {
    legacyCopy(text, done);
  }
}

function legacyCopy(text, done) {
  var holder = document.createElement('textarea');
  holder.value = text;
  document.body.appendChild(holder);
  holder.select();
  try { document.execCommand('copy'); done(); } catch (error) { /* nothing to do */ }
  document.body.removeChild(holder);
}

function markCurrentSection() {
  /* Scoped to the active tab's own view: sec- ids in another tab, or anywhere
     else on the page, must never join this query. */
  var view = byId('tab-' + state.tab);
  if (!view) { return; }
  var targets = view.querySelectorAll('[id^="sec-"]');
  var edge = 90;
  var current = '';
  Array.prototype.forEach.call(targets, function (target) {
    if (target.getBoundingClientRect().top <= edge) { current = target.id; }
  });
  var links = document.querySelectorAll('#sectnav-' + state.tab + ' .sectnav-link');
  Array.prototype.forEach.call(links, function (link) {
    link.classList.toggle('current', link.getAttribute('href') === '#' + current);
  });
}

function subviewOpen(name) {
  var view = byId('tab-' + name);
  if (!view) { return false; }
  return !!(view.querySelector('details[open]') ||
    view.querySelector('.fold-toggle[aria-expanded="true"]'));
}

function showTab(name) {
  state.tab = name;
  Array.prototype.forEach.call(document.querySelectorAll('.tabview'), function (view) {
    view.hidden = view.getAttribute('data-tab') !== name;
  });
  Array.prototype.forEach.call(document.querySelectorAll('.rail-group'), function (group) {
    group.hidden = group.getAttribute('data-tab') !== name;
  });
  Array.prototype.forEach.call(document.querySelectorAll('.tab'), function (button) {
    button.setAttribute('aria-selected', button.getAttribute('data-tab') === name
      ? 'true' : 'false');
  });
  var left = byId('rail-left');
  var group = document.querySelector('#rail-left .rail-group[data-tab="' + name + '"]');
  var empty = !group || !group.querySelector('.card');
  document.body.classList.toggle('no-left', empty);
  if (left) { left.hidden = empty; }
  markCurrentSection();
}

function pressTab(name) {
  /* The view-nesting law. The cheap no-render guard is right for a tab with
     nothing open inside it, and silently wrong the moment that tab holds a
     sub-view: pressing the tab you are already on must still navigate home. */
  if (name === state.tab && !subviewOpen(name)) { return; }
  if (name === state.tab) { closeFolds(byId('tab-' + name)); }
  showTab(name);
  /* A tab press lands at the top of its destination, never at the depth the
     previous view was left at. */
  window.scrollTo(0, 0);
}

function openElsewhere(link) {
  var tab = link.getAttribute('data-tab');
  var target = byId((link.getAttribute('href') || '#').slice(1));
  if (tab) { showTab(tab); }
  if (target) { target.scrollIntoView({ block: 'start' }); }
}

function wire() {
  Array.prototype.forEach.call(document.querySelectorAll('.tab'), function (button) {
    button.addEventListener('click', function () {
      pressTab(button.getAttribute('data-tab'));
    });
  });
  var send = byId('send-compose');
  if (send) {
    send.addEventListener('click', function () {
      var drive = byId('send-drive');
      var output = byId('f-output');
      if (drive && drive.value && output) {
        /* Naming a drive here is an answer, so the data folder stops deciding
           where the cubes go. */
        output.value = drive.value;
        derivedEdits.output = true;
      }
      compose();
      pressTab('now');
    });
  }
  ['filter-year', 'filter-epoch', 'filter-drive', 'filter-status'].forEach(function (id) {
    byId(id).addEventListener('change', filterCatalog);
  });
  byId('find').addEventListener('input', filterCatalog);
  byId('reset').addEventListener('click', function () {
    ['filter-year', 'filter-epoch', 'filter-drive', 'filter-status', 'find'].forEach(
      function (id) { byId(id).value = ''; }
    );
    filterCatalog();
  });
  /* Every composer field recomposes live: there is no Compose button,
     because a button whose work already happened on the keystroke is a
     control that visibly does nothing. */
  ['f-input', 'f-output', 'f-label', 'f-verdict', 'f-registry',
   'f-calibrations', 'f-catalog', 'f-plan', 'f-pattern'].forEach(function (id) {
    var field = byId(id);
    if (field) { field.addEventListener('input', compose); }
  });
  var epoch = byId('f-epoch');
  if (epoch) { epoch.addEventListener('change', compose); }
  Object.keys(DERIVED_FIELDS).forEach(function (key) {
    var field = byId(DERIVED_FIELDS[key]);
    if (!field) { return; }
    field.addEventListener('input', function () { derivedEdits[key] = true; });
  });
  document.addEventListener('click', function (event) {
    if (!event.target || !event.target.closest) { return; }
    var jump = event.target.closest('.xlink');
    if (jump) {
      event.preventDefault();
      openElsewhere(jump);
      return;
    }
    var toggle = event.target.closest('.fold-toggle');
    if (toggle) {
      var body = byId(toggle.getAttribute('aria-controls'));
      var open = toggle.getAttribute('aria-expanded') === 'true';
      if (body) { body.hidden = open; }
      toggle.setAttribute('aria-expanded', open ? 'false' : 'true');
      toggle.textContent = open ? toggle.getAttribute('data-show') : 'hide';
      return;
    }
    var copy = event.target.closest('.copy');
    if (copy) { copyPayload(copy); return; }
    var close = event.target.closest('.drill-close');
    if (close) {
      var drill = close.closest('details');
      if (drill) { drill.open = false; drill.scrollIntoView({ block: 'nearest' }); }
    }
  });
  document.addEventListener('keydown', function (event) {
    if (event.key === 'Escape') { closeFolds(); }
  });
  window.addEventListener('scroll', markCurrentSection, { passive: true });
  filterCatalog();
  compose();
  showTab(state.tab);
}

wire();
"""


# ---------------------------------------------------------------------------
# The served half: a folder picker, and the two states a cold campaign has
# ---------------------------------------------------------------------------
#
# Everything below renders only when ``echelle web`` is serving the page from a
# local server, which is the one situation where the page may ask a question:
# the server answers ``/api/browse`` with the folders on this machine and
# ``/api/home`` with the campaign home it wrote.  The static build never reaches
# any of it -- these constants are concatenated into the page only under
# ``served=True``, so the one-shot file keeps its "reaches nothing" contract
# byte for byte.
#
# The dialog obeys the deep-view law: its Close and its Choose sit in a header
# that never scrolls away, Escape dismisses it from any depth, its body owns its
# own scroll (the page behind it is locked, with a stable scrollbar gutter so
# locking moves no rail by a pixel), and every row is one press target with its
# SIF count as a plain label rather than a second door.

#: The one relative endpoint the picker itself reads.
BROWSE_ENDPOINT = "/api/browse?path="

#: The one relative endpoint that writes a campaign home.
HOME_ENDPOINT = "/api/home"

_PICKER_CSS = """
/* Served only. A stable gutter is declared before anything can lock the page,
   so opening the dialog cannot shift the layout -- or a rail -- sideways. */
html { scrollbar-gutter: stable; }
html.picker-open { overflow: hidden; }
.field-row { display: flex; gap: .4rem; align-items: center; }
.field-row input { flex: 1 1 auto; min-width: 0; }
.field-row .browse { flex: 0 0 auto; white-space: nowrap; }
.picker { position: fixed; inset: 0; z-index: 20; display: flex; align-items: center;
  justify-content: center; padding: var(--gap); }
.picker-backdrop { position: absolute; inset: 0; background: rgba(12, 12, 14, .5); }
/* One bordered panel, headed groups inside it as plain rows: no card in card. */
.picker-dialog { position: relative; display: flex; flex-direction: column;
  width: min(46rem, 100%); max-height: min(80vh, 42rem); background: var(--panel);
  border: 1px solid var(--line); border-radius: .55rem; padding: .8rem .9rem .9rem; }
/* The exit stays within reach at any depth: the head never scrolls. */
.picker-head { flex: 0 0 auto; border-bottom: 1px solid var(--line); padding-bottom: .55rem; }
.picker-title { font-size: 1.05rem; margin: 0 0 .15rem; }
.picker-path { margin: 0 0 .5rem; color: var(--muted);
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace; overflow-wrap: anywhere; }
.picker-go { margin: 0 0 .5rem; }
.picker-go input { width: 100%; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
.picker-actions { display: flex; gap: .5rem; flex-wrap: wrap; margin: 0; }
.picker-choose[disabled] { opacity: .55; cursor: default; }
.picker-error { margin: .5rem 0 0; border-left: .2rem solid var(--bad); color: var(--bad);
  padding-left: .6rem; }
.picker-body { flex: 1 1 auto; min-height: 0; overflow-y: auto; overscroll-behavior: contain;
  display: flex; flex-direction: column; gap: .25rem; padding-top: .55rem; }
.picker-row { display: flex; align-items: baseline; gap: .6rem; width: 100%; text-align: left; }
.picker-name { flex: 1 1 auto; overflow-wrap: anywhere; }
.picker-count { flex: 0 0 auto; color: var(--muted); font-size: .9rem; }
.picker-empty { color: var(--muted); margin: .4rem 0 0; }
.solo { max-width: 46rem; margin: 0 auto; padding: 2rem 24px 3rem; }
.solo h2 { font-size: 1.3rem; }
.solo .cmd { border-top: 1px solid var(--line); }
"""

#: The dialog's markup, identical on every served surface so one script drives
#: all of them.
_PICKER_MARKUP = """<div class="picker" id="picker" hidden>
<div class="picker-backdrop" data-picker-close></div>
<section class="picker-dialog" id="picker-dialog" role="dialog" aria-modal="true"
 aria-labelledby="picker-title" tabindex="-1">
<header class="picker-head"><h2 class="picker-title" id="picker-title">Pick a folder</h2>
<p class="picker-path" id="picker-path">Drives on this machine</p>
<form class="picker-go" id="picker-go">
<input type="text" id="picker-typed" placeholder="or paste a path and press Enter">
</form>
<p class="picker-actions">
<button type="button" class="picker-choose" id="picker-choose" disabled>Choose this folder</button>
<button type="button" class="picker-close" id="picker-close" data-picker-close>Close</button>
</p></header>
<p class="picker-error" id="picker-error" hidden></p>
<div class="picker-body" id="picker-body"></div>
</section></div>
"""

_PICKER_JS = (
    """
/* The folder picker. Two states: the drive list (an empty path) and one
   folder's own subfolders. A row is the press target and descends; the count
   beside it is a label, never a second door. */
var PICKER_BROWSE = '"""
    + BROWSE_ENDPOINT
    + """';
var picker = { path: '', choose: null, opener: null };

function pickerEl(id) { return document.getElementById(id); }

function pickerSay(message) {
  var box = pickerEl('picker-error');
  if (!box) { return; }
  box.textContent = message || '';
  box.hidden = !message;
}

function pickerRow(name, path, note) {
  var row = document.createElement('button');
  row.type = 'button';
  row.className = 'picker-row';
  row.setAttribute('data-path', path === null || path === undefined ? '' : String(path));
  var label = document.createElement('span');
  label.className = 'picker-name';
  label.textContent = name;
  row.appendChild(label);
  var count = document.createElement('span');
  count.className = 'picker-count';
  count.textContent = note;
  row.appendChild(count);
  return row;
}

function pickerRender(payload) {
  var body = pickerEl('picker-body');
  var head = pickerEl('picker-path');
  var choose = pickerEl('picker-choose');
  if (!body) { return; }
  picker.path = String(payload.path || '');
  if (head) { head.textContent = picker.path || 'Drives on this machine'; }
  if (choose) { choose.disabled = !picker.path; }
  while (body.firstChild) { body.removeChild(body.firstChild); }
  if (picker.path) {
    var parent = payload.parent ? String(payload.parent) : '';
    body.appendChild(pickerRow('..', parent, parent ? 'one folder up' : 'the drive list'));
  }
  var entries = picker.path ? (payload.dirs || []) : (payload.drives || []);
  Array.prototype.forEach.call(entries, function (entry) {
    var folder = typeof entry === 'string' ? { name: entry, path: entry } : entry;
    var counted = folder.has_snapshot
      ? 'saved snapshot'
      : typeof folder.sif_count === 'number'
        ? folder.sif_count + ' SIF file(s)'
        : 'SIF files not counted';
    body.appendChild(pickerRow(
      String(folder.name || folder.path || ''), String(folder.path || ''), counted
    ));
  });
  if (!entries.length) {
    var empty = document.createElement('p');
    empty.className = 'picker-empty';
    empty.textContent = picker.path
      ? 'This folder holds no subfolders.'
      : 'This machine listed no drives.';
    body.appendChild(empty);
  }
}

function pickerLoad(path, keepError) {
  /* A refused folder says so in one sentence and still leaves somewhere to go:
     the drive list, which is the one listing with no folder to refuse. The
     retry cannot loop, because the drive list's own path is empty. */
  if (!keepError) { pickerSay(''); }
  /* A NAS can take seconds to answer; silence reads as broken. Say what is
     being read the moment the request leaves, and the render replaces it. */
  var head = pickerEl('picker-path');
  if (head) { head.textContent = 'reading ' + (path || 'the drives') + ' …'; }
  fetch(PICKER_BROWSE + encodeURIComponent(path || '')).then(function (response) {
    return response.json().then(function (payload) {
      return { ok: response.ok, payload: payload };
    });
  }).then(function (answer) {
    if (!answer.ok || answer.payload.error) {
      pickerSay(String(answer.payload.error || 'That folder could not be read.'));
      if (path) { pickerLoad('', true); }
      return;
    }
    pickerRender(answer.payload);
  }, function () {
    pickerSay('This machine did not answer, so no folder could be listed.');
  });
}

function pickerOpen(seed, choose, opener) {
  var box = pickerEl('picker');
  if (!box) { return; }
  picker.choose = choose;
  picker.opener = opener || null;
  box.hidden = false;
  document.documentElement.classList.add('picker-open');
  var dialog = pickerEl('picker-dialog');
  if (dialog) { dialog.focus(); }
  pickerLoad(seed || '');
}

function pickerClose() {
  var box = pickerEl('picker');
  if (!box || box.hidden) { return; }
  box.hidden = true;
  document.documentElement.classList.remove('picker-open');
  if (picker.opener && picker.opener.focus) { picker.opener.focus(); }
  picker.opener = null;
}

function pickerWire() {
  var box = pickerEl('picker');
  if (!box) { return; }
  box.addEventListener('click', function (event) {
    if (!event.target || !event.target.closest) { return; }
    if (event.target.closest('[data-picker-close]')) { pickerClose(); return; }
    if (event.target.closest('.picker-choose')) {
      if (picker.path && picker.choose) { picker.choose(picker.path); }
      return;
    }
    var row = event.target.closest('.picker-row');
    if (row) { pickerLoad(row.getAttribute('data-path') || ''); }
  });
  /* Escape dismisses from any depth, including a long folder list. */
  document.addEventListener('keydown', function (event) {
    if (event.key === 'Escape' && !box.hidden) { pickerClose(); }
  });
  /* A pasted path is navigation: Enter goes there, no clicking down. */
  var go = pickerEl('picker-go');
  if (go) {
    go.addEventListener('submit', function (event) {
      event.preventDefault();
      var typed = pickerEl('picker-typed');
      var value = typed ? typed.value.trim() : '';
      if (value) { pickerLoad(value); }
    });
  }
}

pickerWire();
"""
)

#: The served build's own wiring: Browse fills the composer's data folder and
#: derives everything the composer derives from it, exactly as typing would.
_SERVED_JS = """
function browseInto(button) {
  var field = byId(button.getAttribute('data-browse'));
  if (!field) { return; }
  pickerOpen(field.value, function (chosen) {
    field.value = chosen;
    field.dispatchEvent(new Event('input', { bubbles: true }));
    compose();
    pickerClose();
  }, button);
}

document.addEventListener('click', function (event) {
  if (!event.target || !event.target.closest) { return; }
  var browse = event.target.closest('.browse');
  if (browse) { browseInto(browse); }
});
"""

#: The one-screen pages carry the same command rows the campaign page does, so
#: they carry the same wiring: a show/hide toggle that flips one attribute, and
#: a copy button that carries the whole command whether the row is open or shut.
#: Without this the rows would render controls that look pressable and are not.
_SOLO_JS = """
function soloCopy(button) {
  var text = button.getAttribute('data-copy') || '';
  var done = function () { button.textContent = 'Copied'; };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done, function () { done(); });
    return;
  }
  var holder = document.createElement('textarea');
  holder.value = text;
  document.body.appendChild(holder);
  holder.select();
  try { document.execCommand('copy'); done(); } catch (error) { /* nothing to do */ }
  document.body.removeChild(holder);
}

document.addEventListener('click', function (event) {
  if (!event.target || !event.target.closest) { return; }
  var toggle = event.target.closest('.fold-toggle');
  if (toggle) {
    var body = document.getElementById(toggle.getAttribute('aria-controls'));
    var open = toggle.getAttribute('aria-expanded') === 'true';
    if (body) { body.hidden = open; }
    toggle.setAttribute('aria-expanded', open ? 'false' : 'true');
    toggle.textContent = open ? toggle.getAttribute('data-show') : 'hide';
    return;
  }
  var copy = event.target.closest('.copy');
  if (copy) { soloCopy(copy); }
});
"""

#: The setup page's own wiring: the chosen folder becomes the campaign home.
_SETUP_JS = (
    """
function chooseHome(folder) {
  fetch('"""
    + HOME_ENDPOINT
    + """', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ folder: folder })
  }).then(function (response) {
    return response.json().then(function (payload) {
      return { ok: response.ok, payload: payload };
    });
  }).then(function (answer) {
    if (!answer.ok || answer.payload.error) {
      pickerSay(String(answer.payload.error || 'That folder could not become the campaign home.'));
      return;
    }
    window.location.href = '/';
  }, function () {
    pickerSay('This machine did not answer, so nothing was written.');
  });
}

var opener = document.getElementById('pick-home');
if (opener) {
  opener.addEventListener('click', function () { pickerOpen('', chooseHome, opener); });
}
"""
)


#: The four tabs, in work order: what to do now, the drives, the calibration
#: evidence, and the canon last.
TABS = (
    ("now", "Now"),
    ("drives", "Drives"),
    ("calibration", "Calibration"),
    ("reading", "Reading room"),
)


def _tab_bar() -> str:
    buttons = "".join(
        f'<button type="button" class="tab" data-tab="{_e(key)}" '
        f'aria-selected="{"true" if key == TABS[0][0] else "false"}">{_e(title)}</button>'
        for key, title in TABS
    )
    return f'<nav class="tabs" id="tabs" aria-label="Views">{buttons}</nav>'


def _group(key: str, cards: str) -> str:
    hidden = "" if key == TABS[0][0] else " hidden"
    return f'<div class="rail-group" data-tab="{_e(key)}"{hidden}>{cards}</div>'


def _view(key: str, panels: str) -> str:
    hidden = "" if key == TABS[0][0] else " hidden"
    return f'<section class="tabview" id="tab-{_e(key)}" data-tab="{_e(key)}"{hidden}>{panels}</section>'


def _now_view(context: dict[str, Any]) -> str:
    plan_text = _fill(PLAN_TEMPLATE, context["data"]["values"], "posix")
    absent = "".join(
        _absent_line(row["source"], row["steps"]) for row in context["drive_rows"] if row["absent"]
    )
    connected = "".join(
        _drive_row(row["source"], row["steps"], row["position"], context["drift"])
        for row in context["drive_rows"]
        if not row["absent"]
    )
    if not connected:
        connected = (
            '<p class="note state-unmeasured">No drive in this index answered when the page '
            "was built.</p>"
        )
    return "".join(
        [
            '<section class="panel" id="sec-now-calibrate"><h2>Calibrate</h2>',
            "<p>Once per campaign, at the instrument.</p>",
            _flow(
                context["calibrate_steps"],
                "now-cal",
                done_line="The calibration is in place for this campaign.",
            ),
            "</section>",
            '<section class="panel" id="sec-now-drives"><h2>Drives, each on its own step</h2>',
            "<p>One worker per drive; several at once.</p>",
            connected,
            absent,
            "</section>",
            '<section class="panel" id="sec-plan"><h2>Composed plan and commands</h2>',
            "<p>Editable text, composed from the rail.</p>",
            f'<p class="note" id="unfilled-note">Paste the data folder in the rail first: '
            f"every line below carries {_e(UNFILLED_FOLDER)} until you do.</p>",
            "<h3>Plan TOML</h3>",
            f'<textarea class="plan-out" id="plan-out" spellcheck="false">{_e(plan_text)}</textarea>',
            "<h3>Commands, in campaign order</h3>",
            _composed_commands(context["data"]["values"]),
            "</section>",
        ]
    )


def _page(context: dict[str, Any]) -> str:
    encoded = json.dumps(context["data"], ensure_ascii=False).replace("</", "<\\/")
    navigation = {
        "now": [
            ("sec-now-calibrate", "Calibrate"),
            ("sec-now-drives", "Drives, each on its own step"),
            ("sec-plan", "Composed plan and commands"),
        ],
        "drives": [("sec-drives-cards", "Drives"), ("sec-catalog", "Every cube")],
        "calibration": [("sec-cal-epochs", "Epochs"), ("sec-drift", "Drift evidence")],
        "reading": [
            ("sec-reading-room", "Reading room"),
            *[(document["anchor"], document["title"]) for document in context["documents"]],
        ],
    }
    left = "".join(
        [
            _group(
                "now",
                _composer_card(
                    context["data"]["values"],
                    context["epochs"],
                    context["registry"],
                    served=bool(context.get("served")),
                ),
            ),
            _group(
                "drives",
                _filter_card(context["rows"], context["sources"])
                + _send_to_composer_card(context["data"]["drives"]),
            ),
            # Calibration and the reading room own no controls of their own. An
            # empty rail is the honest answer; invented navigation is not.
            _group("calibration", ""),
            _group("reading", ""),
        ]
    )
    right = "".join(
        [
            _group(
                "now",
                _position_card(context)
                + _step_legend_card()
                + _index_card(navigation["now"], "sectnav-now"),
            ),
            _group(
                "drives",
                _scope_card(context["rows"])
                + _context_card(context)
                + _drive_legend_card()
                + _index_card(navigation["drives"], "sectnav-drives"),
            ),
            _group(
                "calibration",
                _registry_card(context)
                + _verdict_legend_card()
                + _index_card(navigation["calibration"], "sectnav-calibration"),
            ),
            _group("reading", _index_card(navigation["reading"], "sectnav-reading")),
        ]
    )
    views = "".join(
        [
            _view("now", _now_view(context)),
            _view(
                "drives",
                '<section class="panel" id="sec-drives-cards"><h2>Drives</h2>'
                + _source_cards(context["sources"], context["drift"])
                + "</section>"
                + '<section class="panel" id="sec-catalog"><h2>Every cube</h2>'
                + _catalog_table(context["rows"])
                + "</section>",
            ),
            _view(
                "calibration",
                '<section class="panel" id="sec-cal-epochs"><h2>Epochs</h2>'
                + _epoch_table(context)
                + "</section>"
                + '<section class="panel" id="sec-drift"><h2>Drift evidence</h2>'
                + _drift_section(context["drift"])
                + "</section>",
            ),
            _view(
                "reading",
                '<section class="panel" id="sec-reading-room"><h2>Reading room</h2>'
                '<p>The canon travels inside this page; the documentation site is '
                '<a href="https://queezz.github.io/echelle_spectra">'
                "queezz.github.io/echelle_spectra</a>.</p>"
                + _document_sections(context["documents"])
                + "</section>",
            ),
        ]
    )
    # The served half is appended, never woven in: with ``served`` false every
    # piece below is the empty string and the file is the static build's own
    # bytes, down to the banner sentence.
    served = bool(context.get("served"))
    picker_css = _PICKER_CSS if served else ""
    picker_markup = _PICKER_MARKUP if served else ""
    served_js = (_PICKER_JS + _SERVED_JS) if served else ""
    tagline = (
        "Served from this machine. This page never executes commands and never starts a "
        "worker; Browse asks this local server which folders exist."
        if served
        else "Read-only. This page never executes commands, never starts a "
        "worker, and fetches nothing."
    )
    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>Echelle campaign</title>\n"
        f"<style>{_CSS}{picker_css}</style></head><body>\n"
        f'<header class="topbar"><h1>Echelle campaign</h1>{_tab_bar()}'
        f'<span class="tagline">{tagline}</span></header>\n'
        '<div class="wrap"><div class="rail-grid">'
        f'<aside class="rail rail-left" id="rail-left" aria-label="Controls">{left}</aside>'
        f'<main class="content" id="content">{views}</main>'
        f'<aside class="rail rail-right" id="rail-right" aria-label="Context">{right}</aside>'
        "</div></div>\n"
        f"{picker_markup}"
        f"<script>\nconst DATA={encoded};\n{_JS}{served_js}</script>\n</body></html>\n"
    )


# ---------------------------------------------------------------------------
# The two one-screen pages a cold campaign is served
# ---------------------------------------------------------------------------

#: The one value neither a campaign home nor a catalog records: where the raw
#: SIF files are.  It is written as a marker so a pasted command fails loudly on
#: the missing folder rather than quietly on a guessed one.
DATA_FOLDER_MARKER = "<data folder>"


def _solo_page(title: str, tagline: str, body: str, script: str) -> str:
    """One self-contained screen in the page's own palette and CSS discipline."""

    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{_e(title)}</title>\n"
        f"<style>{_CSS}{_PICKER_CSS}</style></head><body>\n"
        f'<header class="topbar"><h1>Echelle campaign</h1>'
        f'<span class="tagline">{_e(tagline)}</span></header>\n'
        f'<main class="solo">{body}</main>\n'
        f"{_PICKER_MARKUP}"
        f"<script>\n{_PICKER_JS}{_SOLO_JS}{script}</script>\n</body></html>\n"
    )


def render_setup_page() -> str:
    """The page a cold start is served: picking the folder *is* the setup.

    There is no campaign home yet, so there is nothing to show and one thing to
    do.  The page names that in plain words, teaches it once, and hands over the
    same picker the composer uses — seeded at this machine's own drives, because
    an example path from another machine is a guess dressed as a default.
    """

    body = (
        '<section class="panel">'
        "<h2>No campaign home yet</h2>"
        "<p>Pick the folder this campaign lives in — the one holding, or about to hold, "
        "its calibrations and its catalog. This machine remembers it, so the next time you "
        "open this page it opens on the campaign instead of on this screen.</p>"
        '<p class="actions"><button type="button" id="pick-home">Pick the campaign folder…'
        "</button></p>"
        "</section>"
    )
    return _solo_page(
        "Echelle campaign — setup",
        "Served from this machine. Choosing a folder writes one small campaign file in it "
        "and changes nothing else.",
        body,
        _SETUP_JS,
    )


def _home_value(home: dict[str, Any], keys: tuple[str, ...], default: str = "") -> str:
    for key in keys:
        value = home.get(key)
        if value:
            return _posix(value)
    return default


def render_empty_campaign_page(home: dict[str, Any]) -> str:
    """The page a campaign home with no catalog yet is served.

    Empty is not broken, and the two are rendered as the different facts they
    are: nothing has been processed here, so the file that would list the cubes
    does not exist.  The page therefore claims no cube and no drive; it carries
    the one command that ends the state, composed from the home's own registry,
    snapshot root and catalog path, with the raw data folder left as a marker
    because nothing on disk records it.
    """

    root = _home_value(home, ("folder", "root", "home", "path"))
    named = str(home.get("volume_label") or home.get("label") or "")
    values = {
        "input": DATA_FOLDER_MARKER,
        "output": _home_value(
            home, ("cubes", "output", "cubes_root"), f"{root}/cubes" if root else "cubes"
        ),
        "registry": _home_value(
            home,
            ("registry", "registry_path"),
            f"{root}/calibration_registry.toml" if root else "calibration_registry.toml",
        ),
        "calibrations": _home_value(
            home,
            ("calibrations", "calibrations_root", "snapshots_root"),
            f"{root}/calibrations" if root else "calibrations",
        ),
        "catalog": _home_value(
            home,
            ("catalog", "catalog_path", "central_index"),
            f"{root}/echelle-catalog.json" if root else "echelle-catalog.json",
        ),
        "sample": str(home.get("sample") or "20"),
        # No drive is catalogued here yet, so nothing names one: the home's own
        # folder name is a default to edit, not a drive this page claims to see.
        "label": named or (root.rsplit("/", 1)[-1] if root else "unknown"),
    }
    meaning, template = STEP_COMMANDS["sample"]
    shapes = {shell: _fill(template, values, shell) for shell, _ in SHELL_NAMES}
    body = (
        '<section class="panel">'
        "<h2>This campaign has no catalog yet</h2>"
        "<p>Its home is in place and nothing has been processed here, so the file that would "
        "list the cubes has not been written. That is empty, not broken.</p>"
        f'<p class="muted">Home {_e(root or "not recorded")} · registry {_e(values["registry"])} '
        f'· snapshot root {_e(values["calibrations"])} · catalog to be written at '
        f'{_e(values["catalog"])}</p>'
        "<h3>The first command</h3>"
        + _command_row(
            "first-run", "Process a first sample", _fill(meaning, values, "posix"), shapes
        )
        + f'<p class="note">The folder holding the raw SIF files is the one value no file here '
        f"records: replace {_e(DATA_FOLDER_MARKER)} with it. This page shows the campaign as "
        "soon as that first run writes the catalog.</p>"
        "</section>"
    )
    return _solo_page(
        "Echelle campaign — no catalog yet",
        "Served from this machine. This page never executes commands and never starts a worker.",
        body,
        "",
    )


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def _packaged_documents() -> list[tuple[str, str, str]]:
    """Return the packaged canon as (origin, name, text) triples."""

    root = files("echelle_spectra.resources").joinpath("reading_room")
    loaded = []
    for name in PACKAGED_DOCUMENTS:
        resource = root.joinpath(name)
        loaded.append(
            (
                f"the installed package (resources/reading_room/{name})",
                name,
                resource.read_text(encoding="utf-8"),
            )
        )
    return loaded


def _documents(document_paths: tuple[str | Path, ...] | list[str | Path]) -> list[dict[str, str]]:
    sources = [*_packaged_documents()]
    for path in document_paths:
        candidate = Path(path)
        sources.append((str(candidate), candidate.name, candidate.read_text(encoding="utf-8")))
    seen: set[str] = set()
    documents = []
    for origin, name, text in sources:
        title = name[:-3].replace("-", " ").replace("_", " ").strip().capitalize()
        anchor = _unique_id("sec-doc-" + (_SLUG.sub("-", name.lower()).strip("-") or "doc"), seen)
        documents.append(
            {
                "name": name,
                "title": title,
                "origin": origin,
                "anchor": anchor,
                "html": render_markdown(text, heading_ids=seen, id_prefix=f"doc-{len(documents)}-"),
            }
        )
    return documents


def _saved_snapshots(calibrations_root: str | Path | None, *, depth: int = 3) -> list[str]:
    """The snapshot folders a calibrations root actually holds, registry or not.

    A saved snapshot is real the moment the bench writes it; a page that only
    believes the registry tells an operator with years of snapshots that "no
    snapshot is in reach", which is false (owner, 2026-08-18: "my WEBUI still
    sees no calibration").  The walk is shallow and bounded because snapshots
    sit at most a couple of levels down whatever folder shape the campaign
    grew (``calibrations/<id>/``, ``calibrations/<day>/calibrations/<id>/``).
    """

    if not calibrations_root:
        return []
    found: set[str] = set()

    def walk(folder: Path, budget: int) -> None:
        try:
            entries = [entry for entry in folder.iterdir() if entry.is_dir()]
        except OSError:
            return
        for entry in entries:
            try:
                is_snapshot = (entry / "snapshot.toml").is_file()
            except OSError:  # pragma: no cover - unreadable mount points
                continue
            if is_snapshot:
                found.add(entry.name)
            elif budget > 1:
                walk(entry, budget - 1)

    walk(Path(calibrations_root), depth)
    return sorted(found)


def _registry_context(
    registry_path: str | Path | None, calibrations_root: str | Path | None
) -> dict[str, Any]:
    """Read the registry for the composer, or say honestly that it was not read."""

    if registry_path is None:
        return {
            "status": "not supplied",
            "path": "",
            "detail": "",
            "epochs": [],
            "epoch_rows": [],
            "saved": _saved_snapshots(calibrations_root),
            "calibrations": _posix(calibrations_root) if calibrations_root else "",
        }
    from .calibration_registry import CalibrationRegistryError, load_calibration_registry

    path = Path(registry_path)
    root = Path(calibrations_root) if calibrations_root else path.parent / "calibrations"
    try:
        registry = load_calibration_registry(path, snapshots_root=root)
    except (CalibrationRegistryError, OSError) as exc:
        return {
            "status": "unreadable",
            "path": _posix(path),
            "detail": str(exc),
            "epochs": [],
            "epoch_rows": [],
            "saved": _saved_snapshots(root),
            "calibrations": _posix(root),
        }
    return {
        "status": "read",
        "path": _posix(path),
        "detail": "",
        "saved": _saved_snapshots(root),
        "epochs": [epoch.snapshot_id for epoch in registry.epochs],
        # The bounds each epoch already declares, carried through so the page
        # can say whether one covers today rather than only counting them.
        "epoch_rows": [
            {
                "snapshot_id": epoch.snapshot_id,
                "shot_from": epoch.shot_from,
                "shot_to": epoch.shot_to,
                "date_from": epoch.date_from.isoformat() if epoch.date_from else "",
                "date_to": epoch.date_to.isoformat() if epoch.date_to else "",
            }
            for epoch in registry.epochs
        ],
        "calibrations": _posix(registry.snapshots_root),
    }


def build_reading_room(
    catalog_path: str | Path,
    output_dir: str | Path,
    *,
    drift_paths: tuple[str | Path, ...] | list[str | Path] = (),
    document_paths: tuple[str | Path, ...] | list[str | Path] = (),
    registry_path: str | Path | None = None,
    calibrations_root: str | Path | None = None,
    served: bool = False,
) -> Path:
    """Build one page; no worker or command execution surface exists either way.

    The static build (``served`` false, the default) is the whole artifact: it
    carries its own data inline, so there is no sidecar for it to read and
    nothing for it to fetch, and its bytes are exactly what they were before the
    served half existed.

    With ``served`` true the same page is built for a local server that can
    answer two questions about this machine — which folders exist, and which one
    is the campaign home.  The composer's data-folder field then carries a
    Browse button opening the folder picker, and only that appended block
    fetches anything.
    """

    loaded = load_catalog(catalog_path)
    merged = loaded.get("schema") == "echelle-merged-catalog/v1"
    catalog = _refresh_availability(loaded)
    sources = catalog.get("sources", [])
    drift = [
        {"path": _posix(path), "evidence": json.loads(Path(path).read_text(encoding="utf-8"))}
        for path in drift_paths
    ]
    registry = _registry_context(registry_path, calibrations_root)
    rows = _rows_from(sources)
    epochs = list(
        dict.fromkeys(
            [*registry["epochs"], *[str(row.get("snapshot_id") or "") for row in rows]]
        )
    )
    epochs = [epoch for epoch in epochs if epoch]
    evidence_name = _next_evidence_name(drift)
    values = _composer_values(
        catalog_path=catalog_path,
        registry=registry,
        evidence_name=evidence_name,
        epochs=epochs,
    )
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    drive_rows = []
    for position, source in enumerate(sources, start=1):
        evidence = _evidence_for_drive(source, drift)
        anchor = f"drift-{drift.index(evidence) + 1}" if evidence is not None else ""
        steps = _drive_steps(
            source,
            evidence=evidence,
            anchor=anchor,
            values=_drive_values(values, source, evidence_name=evidence_name),
            merged=merged,
            catalog_path=_posix(catalog_path),
        )
        drive_rows.append(
            {
                "source": source,
                "steps": steps,
                "position": position,
                "absent": not source.get("available"),
            }
        )
    calibrate_steps = _calibrate_steps(
        registry,
        [str(row.get("snapshot_id") or "") for row in rows],
        values,
        generated_at[:10],
    )
    every_step = [*calibrate_steps, *(step for row in drive_rows for step in row["steps"])]
    context = {
        "served": served,
        "catalog_path": _posix(catalog_path),
        "generated_at": generated_at,
        "sources": sources,
        "rows": rows,
        "drift": drift,
        "calibrate_steps": calibrate_steps,
        "drive_rows": drive_rows,
        "connected_count": sum(1 for row in drive_rows if not row["absent"]),
        "absent_count": sum(1 for row in drive_rows if row["absent"]),
        "ready_count": sum(1 for step in every_step if step["state"] == STEP_READY),
        "blocked_count": sum(1 for step in every_step if step["state"] == STEP_BLOCKED),
        "documents": _documents(document_paths),
        "registry": registry,
        "epochs": epochs,
        "source_count": len(sources),
        "cube_count": len(rows),
        "drift_count": len(drift),
        "data": {
            "values": values,
            # The two facts the page's own derivation needs and the templates
            # do not: the unfilled marker and the next free evidence name.
            "unfilled": UNFILLED_FOLDER,
            "evidence_name": evidence_name,
            "drives": [
                {
                    "label": str(source.get("volume_label", "unknown")),
                    "root": _posix(source.get("drive_root", "")),
                }
                for source in sources
            ],
            "plan_template": PLAN_TEMPLATE,
            "commands": [
                {
                    "id": str(template["id"]),
                    "powershell": str(template["powershell"]),
                    "posix": str(template["posix"]),
                }
                for template in COMMAND_TEMPLATES
            ],
        },
    }
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    index_path = root / "index.html"
    temporary = index_path.with_name(".index.html.tmp")
    temporary.write_text(_page(context), encoding="utf-8", newline="\n")
    os.replace(temporary, index_path)
    return index_path
