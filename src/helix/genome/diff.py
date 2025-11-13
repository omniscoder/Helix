"""Utilities for applying edit events to digital genome sequences."""
from __future__ import annotations

from typing import Iterable

from helix.edit.events import EditEvent


def apply_diffs(sequence: str, events: Iterable[EditEvent]) -> str:
    """
    Apply a sequence of edit events onto `sequence` in the order they were recorded.

    Events are interpreted as edits against the *current* sequence, which means
    overlapping or out-of-order coordinates are allowed. This behaviour matches
    how DigitalGenomeView.apply() composes edits during DAG simulations.
    """

    result = sequence
    for event in events:
        start = max(0, min(len(result), event.start))
        end = max(start, min(len(result), event.end))
        result = result[:start] + event.replacement + result[end:]
    return result
