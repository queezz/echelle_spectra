"""The one folder dialog both GUIs open, and the two habits it has.

The viewer and the calibration bench each ask the same question — *which
folder?* — and each used to ask it with its own hand-built
:class:`QtWidgets.QFileDialog`.  They drifted, as a pair kept in two places
does: the bench learned to show a folder's greyed contents and the viewer
learned it separately, and neither learned what the owner asked for next.  The
dialog lives here now, once, so a habit taught to it is taught to both.

Two habits, beyond Qt's own:

*Typing or pasting a path navigates there.*  Qt's non-native dialog treats its
line edit as the *answer* only — a pasted path fills the field, ``Choose``
returns it, and the folder's contents were never shown.  Owner, 2026-08-18: "I
want a file dialog which allows me to paste in the folder for preview. Now I
can just paste it for selection, not going near."  A path that names a real
directory is therefore read as navigation, and the greyed contents that made
Qt's dialog worth building appear *before* the press rather than never.

*It opens where it was last left.*  A default data directory is right exactly
once per session; every reopening after that is the operator walking back to
where they already were.
"""

from __future__ import annotations

from pathlib import Path

from PyQt5 import QtWidgets

__all__ = [
    "ask_for_folder",
    "build_folder_dialog",
    "configure_folder_dialog",
    "folder_from_typed_text",
    "follow_typed_folder",
    "forget_folders",
    "remember_folder",
    "remembered_folder",
]


def configure_folder_dialog(dialog: QtWidgets.QFileDialog) -> None:
    """Make a folder picker show the files it is not offering to pick.

    Windows' native folder picker hides files outright: a calibration folder
    holding six SIFs reads "No items match your search", so the only way to
    tell one dated folder from another was to open it, look, close it, and open
    the right one — owner, 2026-08-18: "we should show contents, but make them
    gray. Otherwise I have to open the same folder twice."

    Qt's own dialog does exactly the right thing in ``Directory`` mode with
    ``ShowDirsOnly`` switched *off*: the files are listed and greyed, so the
    folder's contents are evidence without ever becoming the answer.  The
    native dialog cannot be told this, which is why this one asks for Qt's by
    name — and asking for Qt's is also what puts a reachable line edit on the
    dialog for :func:`follow_typed_folder` to listen to.  Only folder pickers
    have the defect — a dialog that picks files shows files already — so the
    pattern and previous-pair pickers keep the native ``getOpenFileName`` they
    have always used.
    """

    dialog.setFileMode(QtWidgets.QFileDialog.Directory)
    dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly, False)
    dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
    dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptOpen)


def folder_from_typed_text(text: object) -> Path | None:
    """The existing directory a typed or pasted line names, or ``None``.

    Quoted on purpose: a path copied out of Windows Explorer's address bar or
    out of a shell arrives wrapped in double quotes often enough that refusing
    it would be refusing the exact gesture this exists for.  Anything that is
    not a directory *right now* — a half-typed path, a file, a share that is
    away — is not navigation and is left alone as ordinary text.

    Never ``resolve()``-d: asking the operating system to resolve
    ``\\\\server\\share\\…`` costs a network round trip that can hang while the
    operator is still typing.  A relative path is joined to the working
    directory, which is all it takes to compare it against the dialog's own.
    """

    raw = str(text).strip().strip('"').strip("'")
    if not raw:
        return None
    try:
        candidate = Path(raw).expanduser()
        if not candidate.is_dir():
            return None
    except (OSError, ValueError, RuntimeError):
        # A path the platform will not even parse — a stray ``:`` on Windows,
        # a name past the length limit — is a typo, not a destination.
        return None
    return candidate if candidate.is_absolute() else Path.cwd() / candidate


def follow_typed_folder(dialog: QtWidgets.QFileDialog) -> QtWidgets.QLineEdit | None:
    """Make the dialog's own field navigate, not merely answer.

    Qt's non-native dialog carries exactly one :class:`QtWidgets.QLineEdit` —
    the file-name field — and it is reachable because the dialog is Qt's own
    rather than the platform's.  A directory named there moves the dialog to
    it, so the listing underneath becomes the preview of the folder about to be
    chosen.

    ``setDirectory`` empties the field as a side effect, which would eat the
    text the operator just pasted, so the text is written back inside a guard.
    The guard is what makes this terminate: the write-back re-enters this
    handler with the dialog already standing in that directory, and a directory
    that is already current is not navigated to twice.

    Returns the field it wired, so a caller may test against it, or ``None``
    when the dialog has none — which is what a native dialog looks like from
    here, and is a reason to do nothing rather than to fail.
    """

    edits = dialog.findChildren(QtWidgets.QLineEdit)
    if not edits:
        return None
    edit = edits[0]
    #: A list rather than a ``nonlocal`` flag: this is closed over by two
    #: connections and read in a ``finally``, and a mutable box says so.
    navigating: list[bool] = []

    def navigate(text: object = "") -> None:
        if navigating:
            return
        candidate = folder_from_typed_text(text)
        if candidate is None:
            return
        try:
            standing: Path | None = Path(dialog.directory().absolutePath())
        except (OSError, ValueError):  # pragma: no cover - a dialog without one
            standing = None
        if standing == candidate:
            return
        navigating.append(True)
        try:
            dialog.setDirectory(str(candidate))
            edit.setText(str(text))
        finally:
            navigating.clear()

    edit.textChanged.connect(navigate)
    # Enter on a path that is already the current directory still has to do
    # something: Qt's own accept handles that, and this only guarantees the
    # navigation happened first.
    edit.returnPressed.connect(lambda: navigate(edit.text()))
    # A dialog opened with a path already in its field is the same gesture,
    # made by the caller instead of by the operator.
    navigate(edit.text())
    return edit


#: Where each dialog was last left, keyed by the question it asks.  Session
#: memory only, deliberately: a folder remembered across launches is a folder
#: remembered past the day it was relevant, and the launch argument — or the
#: opened calibration folder — is the better answer at the start of a session.
_LAST_FOLDER: dict[str, str] = {}


def remembered_folder(key: str, fallback: str | Path) -> str:
    """Where this dialog was last left, else what the caller suggested.

    The remembered folder is checked against the filesystem before it is
    offered: a share that has gone away, or a folder deleted since, must not
    open the dialog on nothing.
    """

    remembered = _LAST_FOLDER.get(key)
    if remembered:
        try:
            if Path(remembered).is_dir():
                return remembered
        except OSError:  # pragma: no cover - a share that is away
            pass
    return str(fallback)


def remember_folder(key: str, folder: str | Path) -> None:
    """Record where a dialog was left, for the next one that asks the same thing."""

    chosen = Path(str(folder))
    try:
        home = chosen if chosen.is_dir() else chosen.parent
    except OSError:  # pragma: no cover - a share that went away mid-choice
        return
    _LAST_FOLDER[key] = str(home)


def forget_folders() -> None:
    """Drop every remembered folder — for tests, which must not inherit one."""

    _LAST_FOLDER.clear()


def build_folder_dialog(
    parent, title: str, start_dir: str | Path
) -> QtWidgets.QFileDialog:
    """A configured, path-following folder dialog, not yet shown."""

    dialog = QtWidgets.QFileDialog(parent, title, str(start_dir))
    configure_folder_dialog(dialog)
    follow_typed_folder(dialog)
    return dialog


def ask_for_folder(
    parent, title: str, start_dir: str | Path, *, memory_key: str = ""
) -> str:
    """Ask which folder, and hand back the answer or ``""``.

    One function rather than an inline dialog in each GUI, for the reason both
    GUIs already had one seam apiece: a test must be able to answer this
    without a real modal appearing off-screen and hanging the run.
    """

    opening = remembered_folder(memory_key, start_dir) if memory_key else str(start_dir)
    dialog = build_folder_dialog(parent, title, opening)
    if dialog.exec_() != QtWidgets.QDialog.Accepted:
        return ""
    chosen = dialog.selectedFiles()
    if not chosen:
        return ""
    if memory_key:
        remember_folder(memory_key, chosen[0])
    return chosen[0]
