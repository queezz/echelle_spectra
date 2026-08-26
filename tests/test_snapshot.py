from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from echelle_spectra.snapshot import (
    ROLE_FILENAMES,
    SNAPSHOT_SCHEMA,
    SnapshotError,
    SnapshotValidationError,
    create_snapshot,
    load_snapshot,
)


@pytest.fixture
def calibration_sources(tmp_path: Path) -> dict[str, Path]:
    sources = tmp_path / "sources"
    sources.mkdir()
    answer = {}
    for index, role in enumerate(ROLE_FILENAMES, start=1):
        path = sources / f"source-{role}.dat"
        path.write_bytes((f"{role}-{index}\n" * index).encode())
        answer[role] = path
    return answer


def test_create_snapshot_copies_role_files_and_verifies_digests(
    tmp_path: Path, calibration_sources: dict[str, Path]
) -> None:
    lamp = tmp_path / "h2-source.sif"
    lamp.write_bytes(b"lamp pixels")

    snapshot = create_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20260813_cmos",
        detector="CMOS",
        files=calibration_sources,
        lamps=["H2", "Th-Ar"],
        lamp_files=[("H2", lamp)],
        notes="campaign fixture",
        base_snapshot="20250926_cmos",
        validity={"shot_from": 200000, "date_from": "2026-08-13"},
        alignment={"dx_px": 1.25, "rms_px": 0.42},
        qc={"lines_used": 37},
    )

    assert snapshot.snapshot_id == "20260813_cmos"
    assert snapshot.detector == "cmos"
    assert snapshot.lamps == ("H2", "Th-Ar")
    assert snapshot.manifest["schema"] == SNAPSHOT_SCHEMA
    assert snapshot.manifest["base_snapshot"] == "20250926_cmos"
    assert snapshot.manifest["validity"]["shot_from"] == 200000
    assert {artifact.role for artifact in snapshot.artifacts} == set(ROLE_FILENAMES) | {
        "lamp"
    }
    lamp_artifact = next(artifact for artifact in snapshot.artifacts if artifact.role == "lamp")
    assert lamp_artifact.label == "H2"
    assert {path.name for path in snapshot.root.iterdir()} >= {
        "snapshot.toml",
        *ROLE_FILENAMES.values(),
        "lamps",
    }
    for artifact in snapshot.artifacts:
        target = snapshot.root / artifact.path
        assert hashlib.sha256(target.read_bytes()).hexdigest() == artifact.sha256
        assert target.stat().st_size == artifact.size_bytes


def test_a_destination_that_refuses_file_flags_still_takes_the_save(
    tmp_path: Path, calibration_sources: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A field drive's filesystem may refuse metadata the source carries.

    The real case: Dropbox stamps packaged sources with a BSD file flag
    (UF_TRACKED), and an exFAT trip drive answers ``chflags`` with EINVAL --
    an errno ``shutil.copystat`` does not forgive, so ``copy2`` crashed the
    save after every number was computed.  The save copies bytes first and
    treats unportable metadata as droppable, never fatal.
    """

    import errno
    import shutil as shutil_module

    def refusing_copystat(src, dst, **kwargs):  # noqa: ANN001, ANN003 - test double
        raise OSError(errno.EINVAL, "Invalid argument", str(dst))

    monkeypatch.setattr(shutil_module, "copystat", refusing_copystat)
    snapshot = create_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20260826_cmos",
        detector="CMOS",
        files=calibration_sources,
        lamps=["Ne"],
    )
    for artifact in snapshot.artifacts:
        target = snapshot.root / artifact.path
        assert hashlib.sha256(target.read_bytes()).hexdigest() == artifact.sha256
        source = calibration_sources[artifact.role]
        # The mtime still travels even when the flag copy is refused.
        assert abs(target.stat().st_mtime - source.stat().st_mtime) < 2


def test_changed_artifact_fails_integrity_validation(
    tmp_path: Path, calibration_sources: dict[str, Path]
) -> None:
    snapshot = create_snapshot(
        tmp_path,
        snapshot_id="20260813_cmos",
        detector="cmos",
        files=calibration_sources,
        lamps=["H2"],
    )
    (snapshot.root / "wavelength.txt").write_text("changed", encoding="utf-8")

    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root)

    assert "artifact digest mismatch: wavelength.txt" in caught.value.errors
    assert "artifact size mismatch: wavelength.txt" in caught.value.errors


def test_manifest_path_cannot_escape_snapshot(
    tmp_path: Path, calibration_sources: dict[str, Path]
) -> None:
    snapshot = create_snapshot(
        tmp_path,
        snapshot_id="20260813_cmos",
        detector="cmos",
        files=calibration_sources,
        lamps=["H2"],
    )
    manifest = snapshot.root / "snapshot.toml"
    text = manifest.read_text(encoding="utf-8")
    manifest.write_text(
        text.replace('path = "pattern.txt"', 'path = "../pattern.txt"'),
        encoding="utf-8",
    )

    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root, verify_files=False)

    assert any("path must stay inside the snapshot" in error for error in caught.value.errors)


@pytest.mark.parametrize(
    "snapshot_id",
    ["2026_cmos", "20260813-CMOS", "20260813_cmos_extra", "../20260813_cmos"],
)
def test_invalid_snapshot_ids_are_refused(
    tmp_path: Path, calibration_sources: dict[str, Path], snapshot_id: str
) -> None:
    with pytest.raises(SnapshotError, match="YYYYMMDD"):
        create_snapshot(
            tmp_path,
            snapshot_id=snapshot_id,
            detector="cmos",
            files=calibration_sources,
            lamps=["H2"],
        )


def test_existing_snapshot_is_never_replaced(
    tmp_path: Path, calibration_sources: dict[str, Path]
) -> None:
    create_snapshot(
        tmp_path,
        snapshot_id="20260813_cmos",
        detector="cmos",
        files=calibration_sources,
        lamps=["H2"],
    )

    with pytest.raises(SnapshotError, match="will not be replaced"):
        create_snapshot(
            tmp_path,
            snapshot_id="20260813_cmos",
            detector="cmos",
            files=calibration_sources,
            lamps=["H2"],
        )
