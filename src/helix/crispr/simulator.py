"""
CRISPR sequence-level simulation engine for Helix.

All functions here operate on abstract digital sequences and do not
describe or imply any wet-lab protocols.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from .model import CasSystem, DigitalGenome, GuideRNA, TargetSite


@dataclass
class CutEvent:
    """
    Represents a simulated cut event in a digital genome.

    This is a purely computational object used to reason about
    sequence changes.
    """

    site: TargetSite
    cut_position: int  # genomic coordinate within [site.start, site.end)
    guide: GuideRNA
    cas: CasSystem
    score: float


def find_candidate_sites(
    genome: DigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_sites: Optional[int] = None,
) -> List[TargetSite]:
    """
    Scan the digital genome for candidate target sites for a given guide.

    This is an in-silico pattern search and scoring operation only.
    The implementation is responsible for:
      - applying PAM rules
      - enforcing mismatch thresholds
      - scoring matches

    Returns an ordered list of candidate sites (e.g., by decreasing score).
    """

    # TODO: implement sequence scanning + scoring.
    raise NotImplementedError("find_candidate_sites is not yet implemented.")


def simulate_cuts(
    genome: DigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_events: Optional[int] = None,
) -> List[CutEvent]:
    """
    Simulate cut events in a digital genome for a given CRISPR system.

    This function does not modify the input genome; instead, it
    returns a list of potential cut events with associated scores.

    The actual modeling of repair outcomes can be layered on top in a
    separate module (e.g., indel outcome distributions).
    """

    # TODO: call find_candidate_sites and map to CutEvent objects.
    raise NotImplementedError("simulate_cuts is not yet implemented.")


def rank_off_targets(
    genome: DigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_candidates: int = 1000,
) -> List[TargetSite]:
    """
    Compute an off-target ranking for a given guide.

    Conceptually:
      - enumerate candidate sites
      - score them with a mismatch/PAM penalty model
      - return a ranked list

    All operations are sequence-based simulations.
    """

    # TODO: reuse find_candidate_sites with different filtering/scoring.
    raise NotImplementedError("rank_off_targets is not yet implemented.")
