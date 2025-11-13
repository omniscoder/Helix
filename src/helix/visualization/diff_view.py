"""Sequence diff helpers for EditDAG inspection."""
from __future__ import annotations

import difflib
from typing import Iterable

from helix.genome.digital import DigitalGenomeView


def unified_sequence_diff(before: str, after: str, *, context: int = 60) -> str:
    """Return a unified diff between two sequences for CLI/notebook display."""
    before_chunks = [before[i : i + context] for i in range(0, len(before), context)]
    after_chunks = [after[i : i + context] for i in range(0, len(after), context)]
    diff = difflib.unified_diff(
        before_chunks,
        after_chunks,
        fromfile="before",
        tofile="after",
        lineterm="",
    )
    return "\n".join(diff)


def genome_view_diff(
    root_view: DigitalGenomeView,
    node_view: DigitalGenomeView,
    chrom: str,
    *,
    context: int = 80,
) -> str:
    """Compare a chromosome between two genome views and return a unified diff."""
    before = root_view.materialize_chrom(chrom)
    after = node_view.materialize_chrom(chrom)
    return unified_sequence_diff(before, after, context=context)
