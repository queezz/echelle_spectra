from __future__ import annotations

import re
from pathlib import Path

import echelle_spectra
from echelle_spectra import cli


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


def test_the_umbrella_command_reports_its_own_version(capsys) -> None:
    """A campaign report starts with which build produced it."""

    assert cli.main(["--version"]) == 0
    assert capsys.readouterr().out.strip() == f"echelle {echelle_spectra.__version__}"

    # And the flag is discoverable rather than merely tolerated.
    assert cli.main(["--help"]) == 0
    assert "--version" in capsys.readouterr().out
