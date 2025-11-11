"""Embedded Helix demo datasets and helpers."""
from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Iterable, List

_DATA_PACKAGE = resources.files(__package__)


def available() -> List[str]:
    """Return all dataset paths relative to the datasets package."""
    paths: List[str] = []
    for entry in _DATA_PACKAGE.rglob("*"):
        if entry.is_file():
            paths.append(str(entry.relative_to(_DATA_PACKAGE)))
    return sorted(paths)


def get_path(relative: str) -> Path:
    """Return an absolute Path to an embedded dataset file."""
    resource = _DATA_PACKAGE.joinpath(relative)
    if not resource.is_file():
        raise FileNotFoundError(f"Dataset '{relative}' not found.")
    return Path(resource)


def read_text(relative: str, encoding: str = "utf-8") -> str:
    """Read the text contents of a dataset file."""
    return get_path(relative).read_text(encoding=encoding)
