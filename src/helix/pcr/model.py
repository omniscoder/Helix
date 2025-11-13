"""PCR simulation data models."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

@dataclass
class Primer:
    """
    Digital primer representation for in-silico PCR simulation.
    """

    name: str
    sequence: str  # canonical 5'->3' string
    max_mismatches: int = 2
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class PrimerPair:
    """Forward + reverse primer pair."""

    name: str
    forward: Primer
    reverse: Primer
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class PCRConfig:
    """Abstract PCR configuration for simulation."""

    cycles: int = 10
    per_cycle_efficiency: float = 0.9
    error_rate: float = 1e-3
    max_amplicon_length: int = 2000
    min_amplicon_length: int = 50
    max_amplicons: int = 256
