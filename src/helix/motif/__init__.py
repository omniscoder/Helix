"""Motif discovery helpers (EM/STEME/online)."""

from __future__ import annotations

from typing import Sequence

from .core import MotifResult
from .em import discover_em
from .steme import discover_steme
from .online import discover_online


def discover_motifs(
    sequences: Sequence[str],
    width: int,
    solver: str = "em",
    **kwargs,
) -> MotifResult:
    solver = solver.lower()
    if solver == "em":
        return discover_em(sequences, width, **kwargs)
    if solver == "steme":
        return discover_steme(sequences, width, **kwargs)
    if solver == "online":
        return discover_online(sequences, width, **kwargs)
    raise ValueError(f"Unknown solver '{solver}'.")


__all__ = ["discover_motifs"]
