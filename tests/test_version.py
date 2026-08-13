from __future__ import annotations

import re
from pathlib import Path

import echelle_spectra


def test_package_version_matches_pyproject() -> None:
    root = Path(__file__).parents[1]
    pyproject = root / "pyproject.toml"
    match = re.search(
        r'^version = "([^"]+)"',
        pyproject.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    assert match is not None
    assert echelle_spectra.__version__ == match.group(1)

    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    heading = re.search(r"^## ([0-9]+\.[0-9]+\.[0-9]+)\b", changelog, re.MULTILINE)
    assert heading is not None
    assert heading.group(1) == match.group(1)
