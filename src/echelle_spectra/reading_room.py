"""Build the read-only campaign catalog and evidence reading room."""

from __future__ import annotations

import html
import json
import os
from pathlib import Path
from typing import Any

from .catalog import load_catalog, source_catalog_path


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


def _page(payload: dict[str, Any], documents: dict[str, str]) -> str:
    encoded = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    docs_html = "".join(
        f"<details><summary>{html.escape(name)}</summary><pre>{html.escape(text)}</pre></details>"
        for name, text in documents.items()
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Echelle campaign reading room</title>
<style>
:root {{ color-scheme: light dark; font-family: system-ui, sans-serif; }}
body {{ max-width: 1200px; margin: auto; padding: 1.2rem; }}
.filters,.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:.8rem; }}
.card {{ border:1px solid #8888;border-radius:.5rem;padding:.8rem;margin:.6rem 0; }}
.missing,.insufficient-data {{ border-left:.45rem solid #c77b00; }}
.misaligned-beyond-repair {{ border-left:.45rem solid #c33; }}
.aligned {{ border-left:.45rem solid #198754; }} .shifted {{ border-left:.45rem solid #735ac7; }}
label {{ display:grid;gap:.25rem; }} input,select,textarea {{ font:inherit;padding:.4rem; }}
textarea {{ width:100%;min-height:8rem; }} pre {{ white-space:pre-wrap;overflow-wrap:anywhere; }}
table {{ border-collapse:collapse;width:100%; }} td,th {{ border-bottom:1px solid #8885;padding:.35rem;text-align:left; }}
</style></head><body>
<h1>Echelle campaign reading room</h1>
<p>Read-only catalog, run, drift, procedure, vocabulary, and provenance evidence. This page never executes commands or controls workers.</p>
<section><h2>Browse</h2><div class="filters">
<label>Year<select id="year"><option value="">All</option></select></label>
<label>Epoch<select id="epoch"><option value="">All</option></select></label>
<label>Drive<select id="drive"><option value="">All</option></select></label>
<label>Run status<select id="status"><option value="">All</option></select></label>
</div><div id="catalog"></div></section>
<section><h2>Drift evidence</h2><div id="drift"></div></section>
<section><h2>Plan and command composer</h2>
<p>Outputs remain editable text for pasting into a terminal; they are never run here.</p><div class="grid">
<label>Input folder<input id="input" placeholder="/drive/shots"></label>
<label>Output folder<input id="output" placeholder="/drive/cubes"></label>
<label>Registry<input id="registry" value="calibration_registry.toml"></label>
<label>Drift verdict<input id="verdict" placeholder="drift-evidence.json"></label></div>
<button id="compose">Compose</button><h3>Plan TOML</h3><textarea id="plan"></textarea>
<h3>Paste-ready command</h3><textarea id="command"></textarea></section>
<section><h2>Procedure, vocabulary, and provenance</h2>
<article class="card"><h3>Procedure</h3><p>Inspect status and missing drives, audit a representative interval and named shots, accept only a defensible refinement, then paste the reviewed bulk command into a terminal.</p></article>
<article class="card"><h3>Vocabulary</h3><p>A snapshot is an immutable calibration binder; a catalog is a last-known cube index; missing means unreachable rather than deleted; insufficient-data is a real verdict and never means aligned.</p></article>
<article class="card"><h3>Provenance</h3><p>Cube, source, snapshot, registry, calibration, drift, and recalibration identities are displayed from their saved digests and manifests. Local availability never rewrites historical identity.</p></article>
{docs_html}</section>
<script>
const DATA={encoded};
const sources=DATA.catalog.sources||[];
const rows=sources.flatMap(s=>(s.cubes||[]).map(c=>({{...c,volume_label:s.volume_label,available:s.available,run:s.run||{{}}}})));
const byId=id=>document.getElementById(id);
const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}})[c]);
function values(key){{return [...new Set(rows.map(r=>String(r[key]??'')).filter(Boolean))].sort();}}
for(const [id,key] of [['year','year'],['epoch','snapshot_id'],['drive','volume_label']]){{const el=byId(id);for(const v of values(key))el.add(new Option(v,v));}}
const statusControl=byId('status');for(const v of [...new Set(sources.map(s=>s.run?.state||'unknown'))].sort())statusControl.add(new Option(v,v));
function render(){{const filters={{year:byId('year').value,snapshot_id:byId('epoch').value,volume_label:byId('drive').value}};const selected=rows.filter(r=>Object.entries(filters).every(([k,v])=>!v||String(r[k]??'')===v)&&(!statusControl.value||(r.run.state||'unknown')===statusControl.value));
byId('catalog').innerHTML=sources.filter(s=>!s.available).map(s=>`<div class="card missing"><strong>Missing drive: ${{esc(s.volume_label)}}</strong><p>Last catalog remains discoverable; its files are not currently reachable.</p></div>`).join('')+`<p>${{selected.length}} cube(s)</p><table><tr><th>Year</th><th>Shot</th><th>Epoch</th><th>Drive</th><th>Run</th></tr>${{selected.map(r=>`<tr><td>${{esc(r.year??'unknown')}}</td><td>${{esc(r.shot_number)}}</td><td>${{esc(r.snapshot_id)}}</td><td>${{esc(r.volume_label)}}</td><td>${{esc(r.run.state||'unknown')}}</td></tr>`).join('')}}</table>`;}}
for(const id of ['year','epoch','drive','status'])byId(id).addEventListener('change',render);render();
byId('drift').innerHTML=(DATA.drift||[]).map(d=>{{const kind=['aligned','shifted','misaligned-beyond-repair','insufficient-data'].includes(d.verdict)?d.verdict:'insufficient-data';return `<article class="card ${{kind}}"><h3>${{esc(d.verdict)}}</h3><p>Snapshots: ${{esc((d.snapshot_ids||[]).join(', ')||'unknown')}} · sampled ${{(d.sampled_cubes||[]).length}} cube(s)</p><p>${{esc(d.repair_command||'')}}</p><details><summary>Per-line evidence</summary><table><tr><th>Shot/cube</th><th>Line</th><th>Status</th><th>Residual nm</th><th>SNR</th></tr>${{(d.lines||[]).map(x=>`<tr><td>${{esc(x.shot_number||x.cube)}}</td><td>${{esc(x.line)}}</td><td>${{esc(x.status)}}</td><td>${{esc(x.residual_nm??'—')}}</td><td>${{esc(x.snr??'—')}}</td></tr>`).join('')}}</table></details></article>`;}}).join('')||'<div class="card insufficient-data">No drift evidence supplied. Bulk readiness is unknown.</div>';
byId('compose').onclick=()=>{{const input=byId('input').value,output=byId('output').value,registry=byId('registry').value,verdict=byId('verdict').value;byId('plan').value=`# Hand-edit before use. The reading room does not execute this plan.\n[plan]\ninput = ${{JSON.stringify(input)}}\noutput = ${{JSON.stringify(output)}}\nregistry = ${{JSON.stringify(registry)}}\ndrift_verdict = ${{JSON.stringify(verdict)}}\n`;byId('command').value=`echelle process ${{JSON.stringify(input)}} -o ${{JSON.stringify(output)}} --registry ${{JSON.stringify(registry)}} --drift-verdict ${{JSON.stringify(verdict)}}`;}};
</script></body></html>"""


def build_reading_room(
    catalog_path: str | Path,
    output_dir: str | Path,
    *,
    drift_paths: tuple[str | Path, ...] | list[str | Path] = (),
    document_paths: tuple[str | Path, ...] | list[str | Path] = (),
) -> Path:
    """Build static files only; no worker or command execution surface exists."""

    catalog = _refresh_availability(load_catalog(catalog_path))
    drift = [json.loads(Path(path).read_text(encoding="utf-8")) for path in drift_paths]
    documents = {
        Path(path).name: Path(path).read_text(encoding="utf-8") for path in document_paths
    }
    payload = {"catalog": catalog, "drift": drift}
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    data_path = root / "reading-room.json"
    index_path = root / "index.html"
    temporary = data_path.with_name(".reading-room.json.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    os.replace(temporary, data_path)
    temporary = index_path.with_name(".index.html.tmp")
    temporary.write_text(_page(payload, documents), encoding="utf-8", newline="\n")
    os.replace(temporary, index_path)
    return index_path
