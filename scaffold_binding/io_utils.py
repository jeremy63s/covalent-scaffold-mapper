"""Path resolution helpers.

Every path the pipeline reads or writes comes through here, so nothing is ever
hardcoded to one machine. A path can be supplied as a command-line argument; if
it was not, and the session is interactive, the user is asked for it.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


class PathError(RuntimeError):
    """A required path was missing, unusable, or could not be asked for."""


def _interactive() -> bool:
    return sys.stdin is not None and sys.stdin.isatty()


def _expand(raw: str) -> Path:
    return Path(os.path.expanduser(raw.strip())).resolve()


def resolve_input_path(
    value: str | None,
    prompt: str,
    *,
    expect_dir: bool = False,
    default: str | Path | None = None,
) -> Path:
    """Return an existing path, asking the user for it if it was not supplied.

    `value` is whatever came from the command line (may be None). `default` is
    offered as a suggestion in the prompt and used when the user hits enter.
    """
    if value:
        path = _expand(value)
    elif default is not None and not _interactive():
        path = _expand(str(default))
    else:
        if not _interactive():
            raise PathError(
                f"{prompt}\nNo value was given and there is no terminal to ask on. "
                "Pass it as a command-line argument instead."
            )
        suffix = f" [{default}]" if default is not None else ""
        raw = input(f"{prompt}{suffix}: ").strip()
        if not raw:
            if default is None:
                raise PathError(f"A path is required: {prompt}")
            raw = str(default)
        path = _expand(raw)

    if not path.exists():
        raise PathError(f"Path does not exist: {path}")
    if expect_dir and not path.is_dir():
        raise PathError(f"Expected a directory, found a file: {path}")
    if not expect_dir and path.is_dir():
        raise PathError(f"Expected a file, found a directory: {path}")
    return path


def resolve_output_dir(
    value: str | None,
    prompt: str,
    *,
    default: str | Path | None = None,
) -> Path:
    """Return a writable directory, creating it if needed and asking if unsupplied."""
    if value:
        path = _expand(value)
    elif default is not None and not _interactive():
        path = _expand(str(default))
    else:
        if not _interactive():
            raise PathError(
                f"{prompt}\nNo value was given and there is no terminal to ask on. "
                "Pass it as a command-line argument instead."
            )
        suffix = f" [{default}]" if default is not None else ""
        raw = input(f"{prompt}{suffix}: ").strip()
        path = _expand(raw or str(default))

    path.mkdir(parents=True, exist_ok=True)
    if not os.access(path, os.W_OK):
        raise PathError(f"Directory is not writable: {path}")
    return path


def confirm(question: str, *, default: bool = False) -> bool:
    """Ask a yes/no question, falling back to `default` when non-interactive."""
    if not _interactive():
        return default
    hint = "Y/n" if default else "y/N"
    raw = input(f"{question} [{hint}]: ").strip().lower()
    if not raw:
        return default
    return raw in ("y", "yes")
