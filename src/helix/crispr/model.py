"""
CRISPR simulation models for Helix.

All types in this module are purely computational and intended
for in-silico sequence modeling only, never for wet-lab use.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class CasSystemType(str, Enum):
    """Enumerates supported digital Cas system families."""

    CAS9 = "cas9"
    CAS12A = "cas12a"
    # Extend with additional systems as needed.


@dataclass
class PAMRule:
    """
    Simple pattern-based PAM rule for simulation.

    Example patterns (simulation-only):
      - "NGG" for Cas9-like systems
      - "TTTV" for some Cas12a-like systems

    Pattern syntax is up to the simulator implementation.
    """

    pattern: str
    description: str = ""


@dataclass
class CasSystem:
    """
    Digital representation of a CRISPR effector protein for simulation.

    This does not represent a physical protein or protocol; it only
    encodes parameters for sequence-level modeling.
    """

    name: str
    system_type: CasSystemType
    pam_rules: List[PAMRule]
    cut_offset: int
    """
    Cut offset relative to PAM or guide index, depending on the
    simulator convention. The simulator is responsible for interpreting
    this value consistently.
    """
    max_mismatches: int = 3
    weight_mismatch_penalty: float = 1.0
    weight_pam_penalty: float = 2.0


@dataclass
class GuideRNA:
    """
    Digital guide sequence for CRISPR simulation.

    All sequences here are abstract string representations for
    in-silico modeling only.
    """

    sequence: str  # canonical 5'->3' representation
    pam: Optional[str] = None
    name: Optional[str] = None
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class TargetSite:
    """
    A potential editable site in a digital genome.

    This type holds only coordinate and sequence information for
    simulation; it does not describe any physical construct.
    """

    chrom: str
    start: int  # 0-based start index in the digital genome
    end: int  # end index (exclusive)
    strand: int  # +1 or -1
    sequence: str

    # Optional precomputed scores (simulator may populate these)
    on_target_score: Optional[float] = None
    off_target_score: Optional[float] = None
    pam_match_score: Optional[float] = None


@dataclass
class DigitalGenome:
    """
    Minimal digital genome abstraction for CRISPR simulations.

    This can wrap an existing Helix sequence representation. For now,
    we store only a mapping from chromosome name to a plain string
    sequence. You can later adapt this to Helix's internal types.
    """

    sequences: Dict[str, str]  # chrom -> sequence

    def get_subsequence(self, chrom: str, start: int, end: int) -> str:
        """
        Return a subsequence from the digital genome for simulation.

        Raises KeyError / IndexError if out of range.
        """

        seq = self.sequences[chrom]
        return seq[start:end]
