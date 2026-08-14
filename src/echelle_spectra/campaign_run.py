"""Durable run receipts for interruptible Echelle campaign processing."""

from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
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

# How a run was authorized to touch a calibration epoch. The authorization keys
# are optional additions to the v1 run schema: receipts written before the gate
# existed stay readable and report GATE_UNRECORDED rather than borrowing an
# authorization they never had.
GATE_VERDICT = "verdict"
GATE_SAMPLE = "sample"
GATE_UNGATED = "ungated (no registry)"
GATE_UNRECORDED = "unrecorded (pre-gate receipt)"

# A drive announces its own identity in a small file at its root, so a USB disk
# that returns as E: on Windows and as /Volumes/... on macOS is still one drive
# in every catalog. The label stays beside it for human display only.
DRIVE_ID_NAME = "echelle-drive-id.toml"
DRIVE_ID_SCHEMA = "echelle-drive-id/v1"
READ_ONLY_DRIVE_WARNING = (
    "drive root is not writable, so no stable drive id was stored; catalogs key "
    f"this drive on its volume label alone until {DRIVE_ID_NAME} can be written"
)


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
    """Return a portable best-effort display label when none was supplied.

    This is a heuristic and never an identity; :func:`ensure_drive_identity`
    owns identity.  What it returns per platform:

    * Windows — the drive letter with its separator removed, such as ``E:``.
      Reconnection under a different letter changes the label, which is exactly
      why the stable drive id exists.
    * POSIX — the mount-like folder the operator named, such as ``NIFS-A`` for
      ``/Volumes/NIFS-A/shots``, because the anchor there is only ``/`` and
      names nothing.  The named folder is used when it has a name, then its
      nearest named parent.

    It never returns an empty label: a root that yields no name at all (``/``
    itself) reports ``unknown`` rather than a silent blank.
    """
    resolved = Path(os.path.abspath(source_root))
    anchor = resolved.anchor.rstrip("\\/")
    if anchor:
        return anchor
    for candidate in (resolved, *resolved.parents):
        if candidate.name:
            return candidate.name
    return "unknown"


@dataclass(frozen=True)
class DriveIdentity:
    """One drive's stable id, its display label, and how both were obtained."""

    drive_id: str
    label: str
    root: Path
    created_at: str = ""
    stored: bool = False
    warning: str = ""


def _drive_id_path(root: Path) -> Path:
    return Path(os.path.abspath(root)) / DRIVE_ID_NAME


def read_drive_identity(root: Path) -> DriveIdentity | None:
    """Return the identity a drive already announces, or None when it has none."""
    path = _drive_id_path(root)
    try:
        with path.open("rb") as stream:
            data = tomllib.load(stream)
    except (FileNotFoundError, NotADirectoryError):
        return None
    except (OSError, tomllib.TOMLDecodeError):
        return None
    drive = data.get("drive")
    if data.get("schema") != DRIVE_ID_SCHEMA or not isinstance(drive, dict):
        return None
    drive_id = str(drive.get("id", "")).strip()
    if not drive_id:
        return None
    return DriveIdentity(
        drive_id=drive_id,
        label=str(drive.get("label", "")).strip() or default_volume_label(root),
        root=path.parent,
        created_at=str(drive.get("created_at", "")),
        stored=True,
    )


def ensure_drive_identity(root: Path, *, label: str | None = None) -> DriveIdentity:
    """Return this drive's stable identity, writing it once on first processing.

    An existing ``echelle-drive-id.toml`` is authoritative for the id and is
    never rewritten: relabelling a drive is a deliberate edit of that file, not
    a side effect of a run.  A supplied *label* is still honoured as this run's
    display name.  When the root cannot be written — a read-only archive mount,
    say — the run continues under the label heuristic and says so in its warning
    so the receipt can record why the drive has no stable id.
    """
    existing = read_drive_identity(root)
    if existing is not None:
        return (
            existing
            if not label or label == existing.label
            else DriveIdentity(
                drive_id=existing.drive_id,
                label=label,
                root=existing.root,
                created_at=existing.created_at,
                stored=True,
            )
        )
    identity = DriveIdentity(
        drive_id=uuid.uuid4().hex,
        label=label or default_volume_label(root),
        root=Path(os.path.abspath(root)),
        created_at=utc_now(),
    )
    path = _drive_id_path(root)
    lines = [
        "# Stable identity for this drive, written once by echelle_spectra.",
        "# Catalogs key on the id, so reconnection under another drive letter or",
        "# mount point keeps one drive's history together. Edit the label freely;",
        "# never edit or copy the id onto a second drive.",
        f"schema = {_toml_string(DRIVE_ID_SCHEMA)}",
        "",
        "[drive]",
        f"id = {_toml_string(identity.drive_id)}",
        f"label = {_toml_string(identity.label)}",
        f"created_at = {_toml_string(identity.created_at)}",
        "",
    ]
    try:
        temporary = path.with_name(f".{path.name}.tmp")
        temporary.write_text("\n".join(lines), encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    except OSError:
        return DriveIdentity(
            drive_id="",
            label=identity.label,
            root=identity.root,
            created_at=identity.created_at,
            stored=False,
            warning=READ_ONLY_DRIVE_WARNING,
        )
    return DriveIdentity(
        drive_id=identity.drive_id,
        label=identity.label,
        root=identity.root,
        created_at=identity.created_at,
        stored=True,
    )


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
    drive_id: str = ""
    drive_warning: str = ""
    gate: str = GATE_UNGATED
    sample: bool = False
    sample_files: int = 0
    drift_evidence: str = ""
    drift_evidence_sha256: str = ""
    drift_verdict: str = ""
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
        drive_id: str = "",
        drive_warning: str = "",
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
            drive_id=drive_id,
            drive_warning=drive_warning,
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
            drive_id=run.get("drive_id", ""),
            drive_warning=run.get("drive_warning", ""),
            gate=run.get("gate") or GATE_UNRECORDED,
            sample=bool(run.get("sample", False)),
            sample_files=int(run.get("sample_files", 0)),
            drift_evidence=run.get("drift_evidence", ""),
            drift_evidence_sha256=run.get("drift_evidence_sha256", ""),
            drift_verdict=run.get("drift_verdict", ""),
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

    def exported_outputs(self) -> set[str]:
        """Return every output key this run published, relative to its output root."""
        return {
            str(record["output"])
            for record in self._records
            if record.get("status") == "exported"
        }

    def counts(self) -> dict[str, int]:
        latest: dict[str, str] = {}
        for record in self._records:
            latest[str(record["source"])] = str(record["status"])
        return {
            status: sum(value == status for value in latest.values())
            for status in sorted(TERMINAL_STATUSES)
        }

    def record_authorization(
        self,
        *,
        gate: str,
        sample: bool = False,
        sample_files: int = 0,
        evidence_path: str = "",
        evidence_sha256: str = "",
        verdict: str = "",
    ) -> None:
        """Record how this run earned the right to process its calibration epoch."""
        self.gate = gate
        self.sample = sample
        self.sample_files = sample_files
        self.drift_evidence = evidence_path
        self.drift_evidence_sha256 = evidence_sha256
        self.drift_verdict = verdict
        self.write_manifest()

    def finish(self, state: str) -> None:
        self.state = state
        self.write_manifest()

    def _authorization_lines(self) -> list[str]:
        lines = [
            f"gate = {_toml_string(self.gate)}",
            f"sample = {'true' if self.sample else 'false'}",
        ]
        if self.sample:
            lines.append(f"sample_files = {self.sample_files}")
        if self.drift_evidence:
            lines.append(f"drift_evidence = {_toml_string(self.drift_evidence)}")
            lines.append(
                f"drift_evidence_sha256 = {_toml_string(self.drift_evidence_sha256)}"
            )
            lines.append(f"drift_verdict = {_toml_string(self.drift_verdict)}")
        return lines

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
            f"drive_id = {_toml_string(self.drive_id)}",
            *(
                [f"drive_warning = {_toml_string(self.drive_warning)}"]
                if self.drive_warning
                else []
            ),
            f"snapshot_id = {_toml_string(self.snapshot_id)}",
            f"expected_files = {self.expected_files}",
            *self._authorization_lines(),
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
                "gate": receipt.gate,
                "sample": receipt.sample,
                "drift_evidence": receipt.drift_evidence,
                "drift_verdict": receipt.drift_verdict,
                "source_root": str(receipt.source_root),
                "output_root": str(receipt.output_root),
                "pattern": receipt.pattern,
                "volume_label": receipt.volume_label,
                "drive_id": receipt.drive_id,
                "drive_warning": receipt.drive_warning,
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
