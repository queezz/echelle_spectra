"""Keep pytest scratch outside the Dropbox checkout on every machine."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config: pytest.Config) -> None:
    """Give a default run a fresh OS-local root before pytest builds its factory."""

    if config.option.basetemp is not None:
        return
    scratch = TemporaryDirectory(prefix="echelle-spectra-pytest-")
    config.option.basetemp = Path(scratch.name) / "run"
    config.add_cleanup(scratch.cleanup)
