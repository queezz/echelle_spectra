"""Ordered calibration-epoch registries backed by immutable snapshots.

The registry deliberately stores only snapshot identities.  Shot/date validity
remains in each snapshot's existing ``[validity]`` table, so there is one
authority for epoch boundaries rather than two hand-edited copies.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path, PurePosixPath
from typing import Any

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

from .snapshot import Snapshot, SnapshotError, load_snapshot, validate_snapshot_id

REGISTRY_SCHEMA = "echelle-calibration-registry/v1"
_VALIDITY_KEYS = frozenset({"shot_from", "shot_to", "date_from", "date_to"})
_DATE_PATTERN = re.compile(
    r"(?<!\d)(?P<year>20\d{2})[-_]?(?P<month>0[1-9]|1[0-2])"
    r"[-_]?(?P<day>0[1-9]|[12]\d|3[01])(?!\d)"
)
_NAMED_SHOT_PATTERN = re.compile(r"(?i)(?:^|[^a-z0-9])shot[-_ ]?(\d+)(?!\d)")
# LHD shot numbers run from five to eight digits.  Eight is deliberately inside
# the pattern rather than outside it: a leading ``20240305`` must reach the
# calendar-date guard below and be rejected as a shot, which an upper bound of
# seven digits would silently prevent it from ever doing.  A leading run of
# fewer than five digits is not an LHD shot and is not read as one.
_LHD_SHOT_PATTERN = re.compile(r"^(\d{5,8})(?:\D|$)")


class CalibrationRegistryError(ValueError):
    """The registry or one source-to-epoch selection is invalid."""


@dataclass(frozen=True)
class CalibrationSourceIdentity:
    """Shot/date identity available before loading a source with a calibration."""

    path: Path
    shot_number: int | None = None
    acquisition_date: date | None = None


@dataclass(frozen=True)
class CalibrationEpoch:
    """One ordered registry entry and its already validated snapshot."""

    position: int
    snapshot: Snapshot
    shot_from: int | None
    shot_to: int | None
    date_from: date | None
    date_to: date | None

    @property
    def snapshot_id(self) -> str:
        return self.snapshot.snapshot_id

    @property
    def needs_shot(self) -> bool:
        return self.shot_from is not None or self.shot_to is not None

    @property
    def needs_date(self) -> bool:
        return self.date_from is not None or self.date_to is not None

    def _match_state(self, identity: CalibrationSourceIdentity) -> tuple[bool, set[str]]:
        missing: set[str] = set()
        if self.needs_shot:
            if identity.shot_number is None:
                missing.add("shot number")
            elif not _contains(self.shot_from, self.shot_to, identity.shot_number):
                return False, set()
        if self.needs_date:
            if identity.acquisition_date is None:
                missing.add("acquisition date")
            elif not _contains(self.date_from, self.date_to, identity.acquisition_date):
                return False, set()
        return not missing, missing


@dataclass(frozen=True)
class CalibrationEpochRegistry:
    """A deterministic ordered registry whose entries are validated snapshots."""

    path: Path
    snapshots_root: Path
    epochs: tuple[CalibrationEpoch, ...]

    def resolve(self, identity: CalibrationSourceIdentity) -> CalibrationEpoch:
        """Resolve exactly one epoch or fail without falling back."""

        matches: list[CalibrationEpoch] = []
        missing_by_epoch: list[tuple[CalibrationEpoch, set[str]]] = []
        for epoch in self.epochs:
            matches_epoch, missing = epoch._match_state(identity)
            if matches_epoch:
                matches.append(epoch)
            elif missing:
                missing_by_epoch.append((epoch, missing))

        source = str(identity.path)
        if len(matches) > 1:
            ids = ", ".join(epoch.snapshot_id for epoch in matches)
            raise CalibrationRegistryError(
                f"source {source!r} ambiguously matches calibration snapshots: {ids}"
            )
        if len(matches) == 1:
            if missing_by_epoch:
                possible = ", ".join(
                    f"{epoch.snapshot_id} needs {', '.join(sorted(missing))}"
                    for epoch, missing in missing_by_epoch
                )
                raise CalibrationRegistryError(
                    f"source {source!r} lacks identity needed to exclude another epoch: {possible}"
                )
            return matches[0]
        if missing_by_epoch:
            details = "; ".join(
                f"{epoch.snapshot_id} needs {', '.join(sorted(missing))}"
                for epoch, missing in missing_by_epoch
            )
            raise CalibrationRegistryError(
                f"source {source!r} is missing required calibration identity: {details}"
            )
        raise CalibrationRegistryError(
            "no calibration epoch matches source "
            f"{source!r} (shot={identity.shot_number!r}, "
            f"date={identity.acquisition_date.isoformat() if identity.acquisition_date else None!r})"
        )

    def resolve_source(
        self, source: str | Path, *, root: str | Path | None = None
    ) -> CalibrationEpoch:
        """Parse the supported filename/path identity and resolve one epoch.

        *root* is the source root the operator named.  Supplying it bounds the
        calendar-date scan to that target, so neither a dated mount above it nor
        the current working directory can change which epoch a file selects.
        """

        return self.resolve(source_identity_from_path(source, root=root))


def _contains(start: Any | None, end: Any | None, value: Any) -> bool:
    return (start is None or value >= start) and (end is None or value <= end)


def _intervals_overlap(
    left_start: Any | None,
    left_end: Any | None,
    right_start: Any | None,
    right_end: Any | None,
) -> bool:
    if left_end is not None and right_start is not None and left_end < right_start:
        return False
    if right_end is not None and left_start is not None and right_end < left_start:
        return False
    return True


def _epochs_overlap(left: CalibrationEpoch, right: CalibrationEpoch) -> bool:
    """Return whether two shot/date rectangles can match the same source."""

    return _intervals_overlap(
        left.shot_from, left.shot_to, right.shot_from, right.shot_to
    ) and _intervals_overlap(left.date_from, left.date_to, right.date_from, right.date_to)


def _shot_bound(value: object, *, name: str, snapshot_id: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise CalibrationRegistryError(
            f"snapshot {snapshot_id!r} validity.{name} must be a non-negative integer"
        )
    return value


def _date_bound(value: object, *, name: str, snapshot_id: str) -> date | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise CalibrationRegistryError(
            f"snapshot {snapshot_id!r} validity.{name} must be an ISO YYYY-MM-DD string"
        )
    try:
        parsed = date.fromisoformat(value)
    except ValueError as exc:
        raise CalibrationRegistryError(
            f"snapshot {snapshot_id!r} validity.{name} must be an ISO YYYY-MM-DD string"
        ) from exc
    if parsed.isoformat() != value:
        raise CalibrationRegistryError(
            f"snapshot {snapshot_id!r} validity.{name} must use canonical YYYY-MM-DD"
        )
    return parsed


def _epoch_from_snapshot(position: int, snapshot: Snapshot) -> CalibrationEpoch:
    raw = snapshot.manifest.get("validity")
    if not isinstance(raw, dict):
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} needs a [validity] table for registry use"
        )
    unknown = set(raw) - _VALIDITY_KEYS
    if unknown:
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} validity has unsupported field(s): "
            + ", ".join(sorted(unknown))
        )
    shot_from = _shot_bound(
        raw.get("shot_from"), name="shot_from", snapshot_id=snapshot.snapshot_id
    )
    shot_to = _shot_bound(raw.get("shot_to"), name="shot_to", snapshot_id=snapshot.snapshot_id)
    date_from = _date_bound(
        raw.get("date_from"), name="date_from", snapshot_id=snapshot.snapshot_id
    )
    date_to = _date_bound(raw.get("date_to"), name="date_to", snapshot_id=snapshot.snapshot_id)
    if all(value is None for value in (shot_from, shot_to, date_from, date_to)):
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} validity must declare a shot and/or date bound"
        )
    if shot_from is not None and shot_to is not None and shot_from > shot_to:
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} has shot_from after shot_to"
        )
    if date_from is not None and date_to is not None and date_from > date_to:
        raise CalibrationRegistryError(
            f"snapshot {snapshot.snapshot_id!r} has date_from after date_to"
        )
    return CalibrationEpoch(
        position=position,
        snapshot=snapshot,
        shot_from=shot_from,
        shot_to=shot_to,
        date_from=date_from,
        date_to=date_to,
    )


def _read_registry(path: Path) -> list[object]:
    try:
        with path.open("rb") as stream:
            raw = tomllib.load(stream)
    except FileNotFoundError as exc:
        raise CalibrationRegistryError(f"calibration registry not found: {path}") from exc
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise CalibrationRegistryError(f"cannot read calibration registry {path}: {exc}") from exc
    if raw.get("schema") != REGISTRY_SCHEMA:
        raise CalibrationRegistryError(f"registry schema must be {REGISTRY_SCHEMA!r}")
    raw_epochs = raw.get("epochs")
    if not isinstance(raw_epochs, list) or not raw_epochs:
        raise CalibrationRegistryError("registry must contain a non-empty [[epochs]] list")
    return raw_epochs


def _load_epochs(raw_epochs: list[object], root: Path) -> list[CalibrationEpoch]:
    epochs: list[CalibrationEpoch] = []
    seen: set[str] = set()
    for position, item in enumerate(raw_epochs, start=1):
        if not isinstance(item, dict) or set(item) != {"snapshot_id"}:
            raise CalibrationRegistryError(
                f"registry epochs[{position}] must contain exactly snapshot_id"
            )
        try:
            snapshot_id = validate_snapshot_id(str(item["snapshot_id"]))
            snapshot = load_snapshot(root / snapshot_id)
        except SnapshotError as exc:
            raise CalibrationRegistryError(
                f"registry epochs[{position}] snapshot {item.get('snapshot_id')!r} is invalid: {exc}"
            ) from exc
        if snapshot_id in seen:
            raise CalibrationRegistryError(f"registry repeats snapshot_id {snapshot_id!r}")
        seen.add(snapshot_id)
        epochs.append(_epoch_from_snapshot(position, snapshot))
    return epochs


def _reject_epoch_overlaps(epochs: list[CalibrationEpoch]) -> None:
    for index, left in enumerate(epochs):
        for right in epochs[index + 1 :]:
            if _epochs_overlap(left, right):
                raise CalibrationRegistryError(
                    "calibration epoch overlap is ambiguous between "
                    f"{left.snapshot_id!r} and {right.snapshot_id!r}; boundaries are inclusive"
                )


def load_calibration_registry(
    path: str | Path,
    *,
    snapshots_root: str | Path | None = None,
) -> CalibrationEpochRegistry:
    """Load an ordered registry and verify every referenced snapshot in full."""

    registry_path = Path(path)
    raw_epochs = _read_registry(registry_path)
    root = (
        Path(snapshots_root)
        if snapshots_root is not None
        else registry_path.parent / "calibrations"
    )
    epochs = _load_epochs(raw_epochs, root)
    _reject_epoch_overlaps(epochs)
    return CalibrationEpochRegistry(
        path=registry_path.resolve(),
        snapshots_root=root.resolve(),
        epochs=tuple(epochs),
    )


def _date_scan_text(source: Path, root: Path | None) -> str:
    """Return the path text a source's acquisition date may be read from.

    Without a declared root the path is read exactly as it was written, which is
    all a bare parser call can honestly mean.  With one, the scan is bounded to
    the target the operator named: that root's own name plus the components
    below it.  A dated volume label, home directory, or archive folder above the
    root therefore cannot supply an acquisition date, and one file resolves
    identically whether it was named absolutely or relative to the working
    directory.  A source outside the declared root falls back to its own name.
    """

    if root is None:
        return source.as_posix()
    absolute_source = Path(os.path.normcase(os.path.abspath(source)))
    absolute_root = Path(os.path.normcase(os.path.abspath(root)))
    try:
        relative = absolute_source.relative_to(absolute_root)
    except ValueError:
        return absolute_source.name
    return PurePosixPath(absolute_root.name, *relative.parts).as_posix()


def source_identity_from_path(
    path: str | Path, *, root: str | Path | None = None
) -> CalibrationSourceIdentity:
    """Parse supported LHD shot and calendar-date tokens from a source path.

    The shot number is always read from the file's own stem: a leading five- to
    eight-digit LHD number, or an explicit ``shot_<number>`` token.  An
    eight-digit leading number that is a valid calendar date is a date, not a
    shot.  The acquisition date is read from *root*-bounded path components; see
    :func:`_date_scan_text`.
    """

    source = Path(path)
    text = _date_scan_text(source, Path(root) if root is not None else None)
    dates: set[date] = set()
    for match in _DATE_PATTERN.finditer(text):
        try:
            dates.add(
                date(int(match.group("year")), int(match.group("month")), int(match.group("day")))
            )
        except ValueError:
            continue
    if len(dates) > 1:
        rendered = ", ".join(sorted(value.isoformat() for value in dates))
        raise CalibrationRegistryError(
            f"source {str(source)!r} contains multiple acquisition dates: {rendered}"
        )

    stem = source.stem
    named = _NAMED_SHOT_PATTERN.search(stem)
    lhd = _LHD_SHOT_PATTERN.match(stem)
    shot_number = int(named.group(1)) if named else int(lhd.group(1)) if lhd else None
    if shot_number is not None and len(str(shot_number)) == 8:
        try:
            date.fromisoformat(
                f"{str(shot_number)[:4]}-{str(shot_number)[4:6]}-{str(shot_number)[6:]}"
            )
        except ValueError:
            pass
        else:
            shot_number = None
    return CalibrationSourceIdentity(
        path=source,
        shot_number=shot_number,
        acquisition_date=next(iter(dates), None),
    )
