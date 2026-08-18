"""Keep pytest scratch outside the Dropbox checkout on every machine."""

from __future__ import annotations

import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

# Packet 6 imports repository-local release scripts during tests. Keep their
# bytecode beside neither those scripts nor the synced source package.
sys.dont_write_bytecode = True


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config: pytest.Config) -> None:
    """Give a default run a fresh OS-local root before pytest builds its factory."""

    if config.option.basetemp is not None:
        return
    scratch = TemporaryDirectory(prefix="echelle-spectra-pytest-")
    config.option.basetemp = Path(scratch.name) / "run"
    config.add_cleanup(scratch.cleanup)


@pytest.fixture(autouse=True)
def _folder_dialogs_start_forgetful():
    """No test inherits the folder another test left a dialog standing in.

    The pickers remember where they were last left, for the length of a
    session — and a pytest run is one very long session in which every test
    would otherwise open its dialogs wherever the previous test's tmp_path
    happened to be.
    """

    from echelle_spectra import folder_picker

    folder_picker.forget_folders()
    yield
    folder_picker.forget_folders()
