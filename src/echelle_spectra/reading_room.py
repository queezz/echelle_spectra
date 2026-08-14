"""Build the read-only campaign catalog, composer, and evidence reading room.

The page is a one-shot static build: ``echelle web`` writes one ``index.html``
that carries its own CSS and JavaScript and fetches nothing.  Its structure
follows the house web UI law (fleet's ``WEBUI.md`` and ``WEBUI-COOKBOOK.md``):

* the main column is what a person reads; the rails are what a person presses,
  with controls on the left and context plus the section index on the right;
* both rails are sticky at an offset derived as the sum of the page's own
  header metrics, and only the rails carry ``align-self: start`` so their grid
  row still spans the content column's full height (the zero-travel trap);
* every jump target's id starts with ``sec-`` so one ``scroll-margin-top`` rule
  covers all of them, and the scroll-spy query is scoped to the content column;
* command rows lead with plain words and keep the literal command behind a
  toggle, in both shell shapes, with the copy button carrying the full payload;
* empty, unmeasured, unreachable and judged states are rendered distinctly, and
  a verdict word this page does not know is rendered as unrecognized rather
  than dressed as one it does.
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

    Anything carrying a scheme is dropped to a fragment: the page contract is
    that it reaches nothing outside itself, so a live URL would be a promise
    the page cannot keep.
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
    "measured": "",
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

[plan]
input_dir = "{{input}}"
output_dir = "{{output}}"
pattern = "{{pattern}}"
registry = "{{registry}}"
calibrations = "{{calibrations}}"
drift_verdict = "{{verdict}}"
central_index = "{{catalog}}"
"""

COMMAND_TEMPLATES = (
    {
        "id": "cmd-process",
        "title": "Process this drive in bulk",
        "meaning": (
            "Convert every SIF under the input folder into cubes. The plan file supplies "
            "the input, the output, the registry, the snapshot root and the accepted "
            "drift verdict, so the run is authorized by evidence you have already read, "
            "and one resumable receipt records which evidence that was."
        ),
        "powershell": 'echelle process --plan "{{plan}}"',
        "posix": 'echelle process --plan "{{plan}}"',
    },
    {
        "id": "cmd-audit",
        "title": "Audit this epoch for drift",
        "meaning": (
            "Measure Balmer and Fulcher centroids on a sample of the cubes already on "
            "this drive and write one immutable verdict file. A registry-backed bulk run "
            "is refused until this file exists, and the file is never overwritten."
        ),
        "powershell": (
            'echelle drift audit "{{cubes}}" --every {{every}} --catalog "{{catalog}}" '
            '--calibrations "{{calibrations}}" -o "{{verdict}}"'
        ),
        "posix": (
            'echelle drift audit "{{cubes}}" --every {{every}} --catalog "{{catalog}}" '
            '--calibrations "{{calibrations}}" -o "{{verdict}}"'
        ),
    },
    {
        "id": "cmd-bench",
        "title": "Check the calibration by eye",
        "meaning": (
            "Open the live calibration bench on the calibration folder and look at the "
            "current epoch yourself before trusting any verdict this page renders. "
            "The bench reads; it does not change a saved snapshot."
        ),
        "powershell": 'echelle-calib "{{bench}}"',
        "posix": 'echelle-calib "{{bench}}"',
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


def _composer_values(
    sources: list[dict[str, Any]],
    *,
    catalog_path: str | Path,
    registry: dict[str, Any],
    drift: list[dict[str, Any]],
    epochs: list[str],
) -> dict[str, str]:
    """Pre-fill the composer from what the page actually carries."""

    primary = next(
        (source for source in sources if source.get("available") and source.get("cubes")),
        sources[0] if sources else {},
    )
    cubes = _posix(primary.get("drive_root", "")) or "cubes"
    calibrations = registry.get("calibrations") or "calibrations"
    return {
        # The raw SIF folder is the one thing neither the catalog nor the
        # registry records, so it is the one placeholder here and says so.
        "input": "shots",
        "output": cubes,
        "cubes": cubes,
        "pattern": "*.SIF",
        "registry": registry.get("path") or "calibration_registry.toml",
        "calibrations": calibrations,
        "verdict": (drift[0]["path"] if drift else "drift-evidence.json"),
        "catalog": _posix(catalog_path),
        "plan": "campaign-plan.toml",
        "every": "20",
        "epoch": epochs[0] if epochs else "unassigned",
        "bench": calibrations,
    }


# ---------------------------------------------------------------------------
# Page fragments
# ---------------------------------------------------------------------------


def _card(title: str, body: str, *, classes: str = "") -> str:
    attribute = f"card {classes}".strip()
    return f'<article class="{attribute}"><h3 class="card-title">{title}</h3>{body}</article>'


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


def _source_cards(sources: list[dict[str, Any]]) -> str:
    cards = []
    for source in sources:
        state = source_state(source)
        classes, label = _SOURCE_STATES[state]
        run = source.get("run") or {}
        # A drive with no receipt has no gate to report; borrowing the
        # pre-gate word for it would claim a receipt that never existed.
        gate = _pill(*_gate_state(run.get("gate", GATE_UNRECORDED))) if run else ""
        lines = [
            f'<p class="state-line">{_pill(classes, label)} {gate} '
            f'<span class="muted">{_e(len(source.get("cubes", [])))} cube(s)</span></p>'
        ]
        note = _SOURCE_NOTES[state]
        if note:
            lines.append(f'<p class="note">{_e(note)}</p>')
        if run:
            counts = ", ".join(
                f"{value} {key}" for key, value in sorted((run.get("counts") or {}).items())
            )
            lines.append(
                f'<p class="muted">run {_e(run.get("id"))} [{_e(run.get("state"))}]'
                f"{(' · ' + _e(counts)) if counts else ''}</p>"
            )
            if run.get("drive_warning"):
                lines.append(f'<p class="warn">{_e(run["drive_warning"])}</p>')
        lines.append(
            f'<p class="muted">drive id {_e(source.get("drive_id") or "not announced")} · '
            f'last root {_e(source.get("drive_root") or "unknown")}</p>'
        )
        cards.append(
            _card(_e(source.get("volume_label", "unknown")), "".join(lines), classes=classes)
        )
    if not cards:
        return '<p class="note state-unmeasured">This index names no drives.</p>'
    return f'<div class="card-grid">{"".join(cards)}</div>'


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
    return _card(_e(label), body, classes=f"verdict {state}")


def _drift_section(drift: list[dict[str, Any]]) -> str:
    if not drift:
        return (
            '<p class="note state-unmeasured">No drift evidence was supplied to this build. '
            "Bulk readiness is unmeasured — which is not the same as aligned.</p>"
        )
    return "".join(_drift_card(entry, position) for position, entry in enumerate(drift, start=1))


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


def _text_field(identifier: str, label: str, value: str, *, note: str = "") -> str:
    hint = f'<small class="muted">{_e(note)}</small>' if note else ""
    return (
        f'<label class="field"><span>{_e(label)}</span>'
        f'<input type="text" id="{_e(identifier)}" value="{_e(value)}">{hint}</label>'
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


def _composer_card(
    values: dict[str, str], drives: list[dict[str, str]], epochs: list[str], verdicts: list[str]
) -> str:
    drive_options = "".join(
        f'<option value="{_e(drive["root"])}">{_e(drive["label"])}</option>' for drive in drives
    )
    epoch_options = "".join(f'<option value="{_e(item)}">{_e(item)}</option>' for item in epochs)
    verdict_options = "".join(f'<option value="{_e(item)}">{_e(item)}</option>' for item in verdicts)
    body = (
        '<p class="muted">Pre-filled from this page\'s own catalog and registry. Nothing here '
        "runs; Compose only rewrites the text in the main column.</p>"
        '<div class="fields">'
        f'<label class="field"><span>Drive</span><select id="f-drive">{drive_options}</select></label>'
        f'<label class="field"><span>Epoch</span><select id="f-epoch">{epoch_options}</select></label>'
        f'<label class="field"><span>Drift verdict</span><select id="f-verdict">{verdict_options}'
        "</select></label>"
        + _text_field(
            "f-input",
            "Input folder (raw SIF)",
            values["input"],
            note="Not recorded by any catalog or receipt — the one field this page cannot fill.",
        )
        + _text_field("f-output", "Output folder (cubes)", values["output"])
        + _text_field("f-registry", "Registry", values["registry"])
        + _text_field("f-calibrations", "Snapshot root", values["calibrations"])
        + _text_field("f-catalog", "Merged catalog", values["catalog"])
        + _text_field("f-plan", "Plan file to save", values["plan"])
        + _text_field("f-pattern", "SIF pattern", values["pattern"])
        + _text_field("f-every", "Audit every Nth cube", values["every"])
        + _text_field("f-bench", "Bench folder", values["bench"])
        + "</div>"
        '<p class="actions"><button type="button" id="compose">Compose</button></p>'
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


def _context_card(context: dict[str, Any]) -> str:
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
        f'<p class="muted">Built {_e(context["generated_at"])}</p>'
        f'<p class="muted">{_e(context["source_count"])} drive(s), {_e(context["cube_count"])} '
        f'cube(s), {_e(context["drift_count"])} drift verdict(s)</p>'
        f'<p class="muted">Catalog: {_e(context["catalog_path"])}</p>'
        f"{registry_line}"
    )
    return _card("Context", body)


def _legend_card() -> str:
    entries = [
        ("state-missing-drive", "missing drive", "the catalog did not answer"),
        ("state-unmeasured", "unmeasured", "no receipt; nothing was measured"),
        ("state-empty", "empty", "measured, and it published zero cubes"),
        ("state-insufficient-data", "insufficient-data", "a judged verdict, never aligned"),
        ("state-unrecognized", "unrecognized", "a word this page does not know"),
    ]
    rows = "".join(
        f'<li>{_pill(state, label)} <span class="muted">{_e(meaning)}</span></li>'
        for state, label, meaning in entries
    )
    return _card("What the states mean", f'<ul class="legend">{rows}</ul>')


def _index_card(entries: list[tuple[str, str]]) -> str:
    links = "".join(
        f'<a href="#{_e(anchor)}" class="sectnav-link">{_e(title)}</a>' for anchor, title in entries
    )
    return _card(
        "On this page", f'<nav class="sectnav" id="sectnav">{links}</nav>', classes="rail-card--growing"
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
.topbar .tagline { color: var(--muted); font-size: .82rem; }
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
.rail-card--growing { display: flex; flex: 1 1 auto; flex-direction: column; min-height: 0; }
.rail-card--growing > .fields,
.rail-card--growing > .sectnav { flex: 1 1 auto; min-height: 0; overflow-y: auto;
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
.muted { color: var(--muted); font-size: .85rem; }
.note { border-left: .2rem solid var(--line); padding-left: .6rem; color: var(--muted);
  font-size: .87rem; }
.warn { border-left: .3rem solid var(--judged); padding: .4rem .6rem; background: var(--raised);
  border-radius: .3rem; }
.pill { display: inline-block; border: 1px solid currentColor; border-radius: 999px;
  padding: 0 .5rem; font-size: .76rem; white-space: nowrap; }
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
.scroll-x { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: .87rem; }
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

function composerValues() {
  var values = {};
  Object.keys(DATA.values).forEach(function (key) { values[key] = DATA.values[key]; });
  var fields = {
    input: 'f-input', output: 'f-output', registry: 'f-registry',
    calibrations: 'f-calibrations', catalog: 'f-catalog', plan: 'f-plan',
    pattern: 'f-pattern', every: 'f-every', bench: 'f-bench',
    verdict: 'f-verdict', epoch: 'f-epoch'
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

function closeFolds() {
  var closed = false;
  Array.prototype.forEach.call(document.querySelectorAll('details[open]'), function (item) {
    item.open = false;
    closed = true;
  });
  Array.prototype.forEach.call(document.querySelectorAll('.fold-toggle'), function (button) {
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
  /* Scoped to the content column on purpose: sec- ids elsewhere on the page
     must never join this query. */
  var targets = document.querySelectorAll('#content [id^="sec-"]');
  var edge = 90;
  var current = '';
  Array.prototype.forEach.call(targets, function (target) {
    if (target.getBoundingClientRect().top <= edge) { current = target.id; }
  });
  Array.prototype.forEach.call(document.querySelectorAll('.sectnav-link'), function (link) {
    link.classList.toggle('current', link.getAttribute('href') === '#' + current);
  });
}

function wire() {
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
  byId('compose').addEventListener('click', compose);
  var drive = byId('f-drive');
  if (drive) {
    drive.addEventListener('change', function () {
      if (drive.value) { byId('f-output').value = drive.value; }
      compose();
    });
  }
  document.addEventListener('click', function (event) {
    if (!event.target || !event.target.closest) { return; }
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
  markCurrentSection();
}

wire();
"""


def _page(context: dict[str, Any]) -> str:
    encoded = json.dumps(context["data"], ensure_ascii=False).replace("</", "<\\/")
    index_entries = [
        ("sec-catalog", "Catalog"),
        ("sec-drift", "Drift evidence"),
        ("sec-plan", "Composed plan and commands"),
        ("sec-reading-room", "Reading room"),
        *[(document["anchor"], document["title"]) for document in context["documents"]],
    ]
    left = "".join(
        [
            _filter_card(context["rows"], context["sources"]),
            _composer_card(
                context["data"]["values"],
                context["data"]["drives"],
                context["epochs"],
                context["verdict_paths"],
            ),
        ]
    )
    right = "".join(
        [
            _scope_card(context["rows"]),
            _context_card(context),
            _legend_card(),
            _index_card(index_entries),
        ]
    )
    plan_text = _fill(PLAN_TEMPLATE, context["data"]["values"], "posix")
    main = "".join(
        [
            '<section class="panel" id="sec-catalog"><h2>Catalog</h2>',
            "<p>Every drive this index has ever merged, and every cube it published. A drive that "
            "did not answer keeps its rows and says so.</p>",
            _source_cards(context["sources"]),
            _catalog_table(context["rows"]),
            "</section>",
            '<section class="panel" id="sec-drift"><h2>Drift evidence</h2>',
            _drift_section(context["drift"]),
            "</section>",
            '<section class="panel" id="sec-plan"><h2>Composed plan and commands</h2>',
            "<p>Editable text, composed from the left rail. This page has no code path that runs "
            "a plan, starts a worker, or touches a running batch.</p>",
            "<h3>Plan TOML</h3>",
            f'<textarea class="plan-out" id="plan-out" spellcheck="false">{_e(plan_text)}</textarea>',
            "<h3>Commands</h3>",
            _composed_commands(context["data"]["values"]),
            "</section>",
            '<section class="panel" id="sec-reading-room"><h2>Reading room</h2>',
            "<p>The campaign's own vocabulary, procedure and provenance, rendered from the "
            "documents that ship inside the installed package.</p>",
            _document_sections(context["documents"]),
            "</section>",
        ]
    )
    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>Echelle campaign reading room</title>\n"
        f"<style>{_CSS}</style></head><body>\n"
        '<header class="topbar"><h1>Echelle campaign reading room</h1>'
        '<span class="tagline">Read-only. This page never executes commands, never starts a '
        "worker, and reaches nothing outside itself.</span></header>\n"
        '<div class="wrap"><div class="rail-grid">'
        f'<aside class="rail rail-left" id="rail-left" aria-label="Controls">{left}</aside>'
        f'<main class="content" id="content">{main}</main>'
        f'<aside class="rail rail-right" id="rail-right" aria-label="Context">{right}</aside>'
        "</div></div>\n"
        f"<script>\nconst DATA={encoded};\n{_JS}</script>\n</body></html>\n"
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
            "calibrations": _posix(root),
        }
    return {
        "status": "read",
        "path": _posix(path),
        "detail": "",
        "epochs": [epoch.snapshot_id for epoch in registry.epochs],
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
) -> Path:
    """Build one static page; no worker or command execution surface exists.

    The page is the whole artifact.  It carries its own data inline, so there
    is no sidecar for it to read and nothing for it to fetch.
    """

    catalog = _refresh_availability(load_catalog(catalog_path))
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
    values = _composer_values(
        sources,
        catalog_path=catalog_path,
        registry=registry,
        drift=drift,
        epochs=epochs,
    )
    context = {
        "catalog_path": _posix(catalog_path),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": sources,
        "rows": rows,
        "drift": drift,
        "documents": _documents(document_paths),
        "registry": registry,
        "epochs": epochs,
        "verdict_paths": [entry["path"] for entry in drift] or ["drift-evidence.json"],
        "source_count": len(sources),
        "cube_count": len(rows),
        "drift_count": len(drift),
        "data": {
            "values": values,
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
