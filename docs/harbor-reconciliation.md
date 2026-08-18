# Harbor requirement reconciliation

Status: unreleased implementation candidate at package version 1.6.0. “Fixture”
below means automated synthetic/local evidence, not instrument or field
acceptance.

| Packet / requirement | Code surface | Focused acceptance | Document / evidence status |
| --- | --- | --- | --- |
| 0 immutable digested snapshots | `snapshot.py`, `snapshot_cli.py` | `test_snapshot.py`, `test_snapshot_cli.py` | `calibration-snapshots.md`; previously released evidence |
| 1 durable receipt/resume | `campaign_run.py`, `spectrocube_cli.py` | `test_campaign_run.py` | `campaign-runs.md`; previously released evidence |
| 2 one reader per drive | `spectrocube_cli.py` multi-target executor | `test_campaign_run.py` | `campaign-runs.md`; previously released evidence |
| 3 shared line knowledge | `tools/line_catalog.py`, `tools/line_overlay.py` | `test_line_catalog.py`, `test_line_overlay.py` | `known-line-overlays.md`; previously released evidence |
| 4 live bench core | `calibration_bench.py`, `calibration_bench_gui.py` | bench domain/offscreen tests | `calibration-bench.md`; previously released evidence |
| 5 bench campaign memory | `calibration_campaign.py`, bench GUI | `test_calibration_campaign.py` | `calibration-bench.md`; previously released evidence |
| 6 portable kit | `kit/`, `scripts/build_nifs_kit.py`, reproducible builder | kit/reproducibility tests | `portable-kit.md`; no new rehearsal in this pass |
| 7 neutral optional fields | external SpectroCube 0.2.0 (unchanged) | external 0.2.0 contract suite | external `SPEC.md`; dependency remains pinned at `0b02ac9` |
| 8 registry + complete cube provenance | `calibration_registry.py`, `tools/spectrocube_export.py` | registry/complete-provenance tests | `calibration-epoch-registry.md`; Packet 8 contracts preserved |
| 9 per-drive catalog | `catalog.build_drive_catalog`, automatic batch hook, `campaign_run.ensure_drive_identity`, `campaign_run.resolve_receipt_directory` | `test_catalog_identity.py` drive-id, reconnection, read-only fallback, gate-row, and input/output-seam cases (rehearsal with runs-root paths; refusals) | `campaign-runs.md`, `operator-cheat-sheet.md` step 8 |
| 9 merged missing-drive index | `catalog.merge_catalogs` (recency, drive-relative paths), `merge_into_central_index`, `refresh_catalog_row`, reading-room availability refresh | disconnected catalog fixture; stale-after-fresh, relocation, auto-merge, and recal-refresh cases | `campaign-runs.md` |
| 9 cube-derived text | `lhd_text.py`, `campaign_tools_cli.txt_main` | frozen-header/cube fixture and missing-timing refusal; documented `export-timing.toml` carried through a registry run to the cube | `harbor-candidate.md`, `operator-cheat-sheet.md` "Getting to `echelle txt`" |
| 9 one writer, two frozen dialects | legacy `Spectrum.save` and GUI call `lhd_text`; both legacy templates restored to `resources/` | `test_lhd_text_header.py` golden diffs against `tests/golden/` | `harbor-candidate.md` |
| 10 wavelength/factor deltas | `recalibration.recalibrate_dataset/cube` | A→B wavelength and factor fixture | `harbor-candidate.md` |
| 10 old/new provenance + manifest | recalibration history attrs and adjacent manifest | history event assertions; installed-file gate deferred | `harbor-candidate.md` |
| 10 geometry refusal | pattern digest comparison | changed-pattern refusal names raw SIF | `harbor-candidate.md` |
| 11 interval/shot sampling | `drift.select_sample_paths`, audit CLI | shifted synthetic cube and sample-rule evidence | `harbor-candidate.md` |
| 11 four verdicts + thresholds | `drift.centroid_evidence/verdict_from_evidence` | all four verdict branches | `harbor-candidate.md` |
| 11 repair/refinement | repair command + `create_refinement_snapshot` | immutable `-r1` and accepted evidence fixture | `harbor-candidate.md` |
| 11 pre-bulk gate | `require_sampled_verdict`, registry bulk hook | shifted refusal and aligned acceptance | `harbor-candidate.md` |
| 12 browsing/verdict cards | `reading_room.py` one-file static build: two sticky rails, `sec-` anchors, local Find, v2 drift evidence rendering | `test_reading_room_page.py` rail/anchor/Find, four-honesty-state, unrecognized-verdict and evidence cases; missing-drive fixture | `harbor-candidate.md` |
| 12 composer/no worker control | rail composer asking two inputs (data folder, calibration epoch) and deriving the rest; composed sample/wavelength-alignment-check/bulk-convert commands read in campaign order, the last unlocked by the verdict | both shell shapes and full copy payload per command; page contains no execution or network surface; visual/perimeter gate deferred | `harbor-candidate.md` |
| 12 procedure/vocabulary/provenance | `resources/reading_room/*.md` as the single source, rendered by the stdlib Markdown renderer from the installed package | packaged-resource build from an unrelated working directory; renderer subset case | `harbor-candidate.md`; the docs tree does not restate the packaged canon |
| 13 historical absorption | `historical.py`, three bundled TOMLs, `snapshot import-historical` | all three epochs import, register, and resolve a dated source; missing artifact refused by name | `calibration-epoch-registry.md` |
| 13 connected fixture path | snapshot + receipt + catalog + text + drift/refinement + optional recalibration domain + web | `test_refinement_historical_and_connected_fixture_path` and adjacent focused tests | this reconciliation |

## Deferred validation ledger

| Gate | Candidate state | Required later evidence |
| --- | --- | --- |
| Fable correctness/coverage review | pending | review correctness, missing acceptance, unsafe assumptions, and integration defects |
| Full source suite | 270 passed with established NumPy/netCDF warnings | repeat after review if implementation changes |
| Wheel/sdist, Twine, archive inspection | deliberately omitted | external build outputs, identities, metadata, contents |
| Clean and offline installation | deliberately omitted | fresh external kit/runtime installs without repository state |
| Reproducibility build | deliberately omitted | independent byte-identical artifacts |
| Windows/macOS portable-kit/native rehearsal | deliberately omitted | matching native executions and payload checks |
| Real 2026 snapshot and cube recalibration | unavailable in fixture pass | immutable real snapshot and A→B result |
| Sampled real historical epoch | unavailable in fixture pass | per-shot Balmer/Fulcher evidence and honest verdict |
| Reading-room visual/perimeter/accessibility gate | deliberately omitted | viewport, keyboard, missing/error, and content perimeter evidence |
| Expensive strict documentation/visual gates | deliberately omitted | final strict docs and visual evidence |
| Release versions/tags/push/Harbor shipment | forbidden in this pass | owner-reviewed release certification and explicit owner actions |
