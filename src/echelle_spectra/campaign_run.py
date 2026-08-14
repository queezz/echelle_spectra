"""Durable run receipts for interruptible Echelle campaign processing."""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised only on older Python
    import tomli as tomllib

RUN_SCHEMA = "echelle-run/v1"
RECORD_SCHEMA = "echelle-run-record/v1"
TERMINAL_STATUSES = {"exported", "skipped", "failed", "interrupted"}


def utc_now() -> str:
    """Return a stable UTC timestamp for machine-readable receipts."""
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Digest a file without loading a campaign-sized artifact into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _toml_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "-", value).strip("-").lower()
    return cleaned[:40] or "run"


def _relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def default_volume_label(source_root: Path) -> str:
    """Return a portable best-effort drive identity when no label was supplied."""
    anchor = source_root.resolve().anchor.rstrip("\\/")
    return anchor or source_root.resolve().name or "unknown"


def new_run_directory(runs_root: Path, source_root: Path) -> Path:
    """Reserve a dated run directory without replacing an existing receipt."""
    stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    base = runs_root / f"{stamp}-{_slug(source_root.name)}"
    candidate = base
    revision = 2
    while candidate.exists():
        candidate = runs_root / f"{base.name}-r{revision}"
        revision += 1
    return candidate


def target_runs_root(runs_root: Path, source_root: Path) -> Path:
    """Return a stable isolated receipt root for one multi-target source."""
    resolved = str(source_root.resolve()).encode("utf-8")
    path_digest = hashlib.sha256(resolved).hexdigest()[:8]
    return runs_root / f"{_slug(source_root.name)}-{path_digest}"


@dataclass(frozen=True)
class SourceIdentity:
    """Content identity recorded before processing starts."""

    path: Path
    key: str
    size_bytes: int
    sha256: str


@dataclass
class RunReceipt:
    """Append-only per-source records plus an atomically replaced run summary."""

    directory: Path
    source_root: Path
    output_root: Path
    pattern: str
    volume_label: str
    snapshot_id: str = "unassigned"
    created_at: str = field(default_factory=utc_now)
    expected_files: int = 0
    state: str = "running"
    _records: list[dict[str, Any]] = field(default_factory=list, repr=False)

    @property
    def manifest_path(self) -> Path:
        return self.directory / "run.toml"

    @property
    def records_path(self) -> Path:
        return self.directory / "records.jsonl"

    @classmethod
    def create(
        cls,
        directory: Path,
        *,
        source_root: Path,
        output_root: Path,
        pattern: str,
        volume_label: str,
        snapshot_id: str,
        expected_files: int,
    ) -> RunReceipt:
        if directory.exists():
            raise FileExistsError(f"Run directory already exists: {directory}")
        directory.mkdir(parents=True)
        receipt = cls(
            directory=directory.resolve(),
            source_root=source_root.resolve(),
            output_root=output_root.resolve(),
            pattern=pattern,
            volume_label=volume_label,
            snapshot_id=snapshot_id or "unassigned",
            expected_files=expected_files,
        )
        receipt.write_manifest()
        return receipt

    @classmethod
    def load(cls, directory: Path) -> RunReceipt:
        manifest_path = directory / "run.toml"
        with manifest_path.open("rb") as stream:
            data = tomllib.load(stream)
        if data.get("schema") != RUN_SCHEMA:
            raise ValueError(f"Unsupported run receipt schema in {manifest_path}")
        run = data["run"]
        receipt = cls(
            directory=directory.resolve(),
            source_root=Path(run["source_root"]),
            output_root=Path(run["output_root"]),
            pattern=run["pattern"],
            volume_label=run["volume_label"],
            snapshot_id=run.get("snapshot_id") or "unassigned",
            created_at=run["created_at"],
            expected_files=int(run.get("expected_files", 0)),
            state=run.get("state", "unknown"),
        )
        receipt._records = list(read_records(receipt.records_path))
        return receipt

    def identity_for(self, source: Path) -> SourceIdentity:
        return SourceIdentity(
            path=source.resolve(),
            key=_relative_or_absolute(source, self.source_root),
            size_bytes=source.stat().st_size,
            sha256=sha256_file(source),
        )

    def completed_output_is_valid(
        self,
        source: SourceIdentity,
        output: Path,
        *,
        snapshot_id: str | None = None,
    ) -> bool:
        """Return true only when a prior export still matches source and output."""
        output = output.resolve()
        for record in reversed(self._records):
            if record.get("status") != "exported" or record.get("source") != source.key:
                continue
            if record.get("source_sha256") != source.sha256:
                return False
            if int(record.get("source_size_bytes", -1)) != source.size_bytes:
                return False
            if snapshot_id is not None and record.get("snapshot_id") != snapshot_id:
                return False
            if record.get("output") != _relative_or_absolute(output, self.output_root):
                return False
            if not output.is_file():
                return False
            if int(record.get("output_size_bytes", -1)) != output.stat().st_size:
                return False
            return record.get("output_sha256") == sha256_file(output)
        return False

    def has_export_record(self, source: SourceIdentity, *, snapshot_id: str | None = None) -> bool:
        """Return whether this run ever published an output for the source key."""
        return any(
            record.get("status") == "exported"
            and record.get("source") == source.key
            and (snapshot_id is None or record.get("snapshot_id") == snapshot_id)
            for record in self._records
        )

    def append(
        self,
        source: SourceIdentity,
        output: Path,
        *,
        status: str,
        started_at: str,
        finished_at: str,
        duration_s: float,
        reason: str = "",
        snapshot_id: str | None = None,
    ) -> dict[str, Any]:
        if status not in TERMINAL_STATUSES:
            raise ValueError(f"Unsupported receipt status: {status}")
        output = output.resolve()
        record: dict[str, Any] = {
            "schema": RECORD_SCHEMA,
            "source": source.key,
            "source_path": str(source.path),
            "source_size_bytes": source.size_bytes,
            "source_sha256": source.sha256,
            "volume_label": self.volume_label,
            "snapshot_id": snapshot_id or self.snapshot_id,
            "output": _relative_or_absolute(output, self.output_root),
            "output_path": str(output),
            "status": status,
            "reason": reason,
            "started_at": started_at,
            "finished_at": finished_at,
            "duration_s": round(duration_s, 6),
        }
        if status == "exported" and output.is_file():
            record["output_size_bytes"] = output.stat().st_size
            record["output_sha256"] = sha256_file(output)
        self.directory.mkdir(parents=True, exist_ok=True)
        with self.records_path.open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        self._records.append(record)
        self.write_manifest()
        return record

    def counts(self) -> dict[str, int]:
        latest: dict[str, str] = {}
        for record in self._records:
            latest[str(record["source"])] = str(record["status"])
        return {
            status: sum(value == status for value in latest.values())
            for status in sorted(TERMINAL_STATUSES)
        }

    def finish(self, state: str) -> None:
        self.state = state
        self.write_manifest()

    def write_manifest(self) -> None:
        counts = self.counts()
        lines = [
            f"schema = {_toml_string(RUN_SCHEMA)}",
            "",
            "[run]",
            f"id = {_toml_string(self.directory.name)}",
            f"state = {_toml_string(self.state)}",
            f"created_at = {_toml_string(self.created_at)}",
            f"updated_at = {_toml_string(utc_now())}",
            f"source_root = {_toml_string(str(self.source_root))}",
            f"output_root = {_toml_string(str(self.output_root))}",
            f"pattern = {_toml_string(self.pattern)}",
            f"volume_label = {_toml_string(self.volume_label)}",
            f"snapshot_id = {_toml_string(self.snapshot_id)}",
            f"expected_files = {self.expected_files}",
            "",
            "[counts]",
            *(
                f"{status.replace('-', '_')} = {counts.get(status, 0)}"
                for status in sorted(TERMINAL_STATUSES)
            ),
            "",
        ]
        temporary = self.manifest_path.with_suffix(".toml.tmp")
        temporary.write_text("\n".join(lines), encoding="utf-8", newline="\n")
        os.replace(temporary, self.manifest_path)


def read_records(path: Path) -> Iterable[dict[str, Any]]:
    if not path.is_file():
        return []
    records = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("schema") != RECORD_SCHEMA:
                raise ValueError(f"Unsupported record at {path}:{line_number}")
            records.append(record)
    return records


def find_resumable_run(
    runs_root: Path, source_root: Path, output_root: Path, pattern: str
) -> Path | None:
    """Find the newest unfinished run matching this source/output pair."""
    if not runs_root.is_dir():
        return None
    matches: list[Path] = []
    for manifest in runs_root.glob("*/run.toml"):
        try:
            receipt = RunReceipt.load(manifest.parent)
        except (OSError, ValueError, KeyError):
            continue
        if receipt.state not in {"running", "partial", "interrupted"}:
            continue
        if receipt.source_root.resolve() != source_root.resolve():
            continue
        if receipt.output_root.resolve() != output_root.resolve():
            continue
        if receipt.pattern != pattern:
            continue
        matches.append(manifest.parent)
    return max(matches, key=lambda item: item.name) if matches else None


def list_run_summaries(runs_root: Path) -> list[dict[str, Any]]:
    """Load valid run manifests for the status command, newest first."""
    summaries = []
    if not runs_root.is_dir():
        return summaries
    for manifest in runs_root.rglob("run.toml"):
        try:
            receipt = RunReceipt.load(manifest.parent)
        except (OSError, ValueError, KeyError):
            continue
        summaries.append(
            {
                "id": receipt.directory.name,
                "state": receipt.state,
                "counts": receipt.counts(),
                "expected_files": receipt.expected_files,
                "snapshot_id": receipt.snapshot_id,
                "source_root": str(receipt.source_root),
                "output_root": str(receipt.output_root),
                "pattern": receipt.pattern,
                "volume_label": receipt.volume_label,
                "updated_at": manifest.stat().st_mtime,
            }
        )
    return sorted(summaries, key=lambda item: item["updated_at"], reverse=True)


def latest_run_summaries(runs_root: Path) -> list[dict[str, Any]]:
    """Return the newest receipt for each independent source/output target."""
    latest: dict[tuple[str, str, str], dict[str, Any]] = {}
    for summary in list_run_summaries(runs_root):
        key = (
            os.path.normcase(str(Path(summary["source_root"]).resolve())),
            os.path.normcase(str(Path(summary["output_root"]).resolve())),
            str(summary["pattern"]),
        )
        latest.setdefault(key, summary)
    return list(latest.values())
