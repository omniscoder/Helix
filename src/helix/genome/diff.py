"""Utilities for applying edit events to digital genome sequences."""
from __future__ import annotations

from typing import Iterable

from helix.edit.events import EditEvent


def apply_diffs(sequence: str, events: Iterable[EditEvent]) -> str:
    """
    Apply ordered, non-overlapping edit events to a sequence and return the result.

    For v1 we assume callers provide events sorted by start coordinate and that
    ranges do not overlap. Later versions can enforce/sanitize this.
    """

    events = list(events)
    if not events:
        return sequence

    pieces = []
    cursor = 0
    for event in events:
        if event.start < cursor:
            raise ValueError("EditEvents must be non-overlapping and sorted.")
        pieces.append(sequence[cursor : event.start])
        pieces.append(event.replacement)
        cursor = event.end
    pieces.append(sequence[cursor:])
    return "".join(pieces)
