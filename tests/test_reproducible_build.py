"""Tests for the normalized byte-identical distribution contract."""

from __future__ import annotations

import gzip
import io
import tarfile
import zipfile
from pathlib import Path

import pytest

from scripts.companion_wheels import _clean_source
from scripts.nifs_kit import KitError
from scripts.reproducible_build import (
    _validate_clean_source,
    compare_distributions,
    normalize_sdist,
    normalize_wheel,
)

EPOCH = 1786579200


def _raw_wheel(path: Path, *, reverse: bool, timestamp: tuple[int, ...]) -> None:
    entries = [("demo/data.txt", b"data"), ("demo-1.0.dist-info/METADATA", b"Version: 1.0\n")]
    if reverse:
        entries.reverse()
    with zipfile.ZipFile(path, "w") as bundle:
        for name, content in entries:
            info = zipfile.ZipInfo(name, date_time=timestamp)
            info.create_system = 0 if reverse else 3
            info.external_attr = (0o777 if reverse else 0o600) << 16
            bundle.writestr(info, content)


def _raw_sdist(path: Path, *, reverse: bool, mtime: int) -> None:
    entries = [("demo-1.0/a.txt", b"a"), ("demo-1.0/b.txt", b"b")]
    if reverse:
        entries.reverse()
    with path.open("wb") as raw, gzip.GzipFile(
        filename="different-name" if reverse else "name",
        mode="wb",
        fileobj=raw,
        mtime=mtime,
    ) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as bundle:
            for name, content in entries:
                info = tarfile.TarInfo(name)
                info.size = len(content)
                info.mtime = mtime
                info.uid = 1000 if reverse else 2000
                info.gid = 1000 if reverse else 2000
                info.uname = "machine-a" if reverse else "machine-b"
                info.gname = info.uname
                info.mode = 0o777 if reverse else 0o600
                bundle.addfile(info, io.BytesIO(content))


def test_wheel_normalization_removes_order_time_and_permission_drift(tmp_path: Path) -> None:
    first_raw = tmp_path / "first.whl"
    second_raw = tmp_path / "second.whl"
    first = tmp_path / "first-normalized.whl"
    second = tmp_path / "second-normalized.whl"
    _raw_wheel(first_raw, reverse=False, timestamp=(2024, 1, 2, 3, 4, 6))
    _raw_wheel(second_raw, reverse=True, timestamp=(2026, 8, 13, 12, 0, 0))
    normalize_wheel(first_raw, first, EPOCH)
    normalize_wheel(second_raw, second, EPOCH)
    assert first.read_bytes() == second.read_bytes()


def test_sdist_normalization_removes_gzip_tar_and_owner_drift(tmp_path: Path) -> None:
    first_raw = tmp_path / "first.tar.gz"
    second_raw = tmp_path / "second.tar.gz"
    first = tmp_path / "first-normalized.tar.gz"
    second = tmp_path / "second-normalized.tar.gz"
    _raw_sdist(first_raw, reverse=False, mtime=1_700_000_000)
    _raw_sdist(second_raw, reverse=True, mtime=1_780_000_000)
    normalize_sdist(first_raw, first, EPOCH)
    normalize_sdist(second_raw, second, EPOCH)
    assert first.read_bytes() == second.read_bytes()


def test_distribution_comparison_requires_names_and_bytes(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    for root in (first, second):
        (root / "demo.whl").write_bytes(b"wheel")
        (root / "demo.tar.gz").write_bytes(b"sdist")
    assert set(compare_distributions(first, second)) == {"demo.whl", "demo.tar.gz"}
    (second / "demo.whl").write_bytes(b"changed")
    with pytest.raises(KitError, match="not reproducible"):
        compare_distributions(first, second)


def test_live_git_checkout_is_not_accepted_as_clean_source(tmp_path: Path) -> None:
    source = tmp_path / "source"
    (source / "src").mkdir(parents=True)
    (source / ".git").mkdir()
    (source / "pyproject.toml").write_text("[project]\nname='demo'\n", encoding="utf-8")
    (source / "README.md").write_text("demo", encoding="utf-8")
    with pytest.raises(KitError, match="exported clean copy"):
        _validate_clean_source(source)
    with pytest.raises(KitError, match="exported copy"):
        _clean_source(source, "companion")
