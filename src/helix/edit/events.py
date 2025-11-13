"""Edit event primitives."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict


@dataclass(frozen=True)
class EditEvent:
    """
    Local sequence transformation in a digital genome view.

    All coordinates are half-open [start, end) and refer to the view's
    current coordinate system.
    """

    chrom: str
    start: int
    end: int
    replacement: str
    metadata: Dict[str, Any] = field(default_factory=dict)
