"""
CRISPR sequence-level simulation engine for Helix.

All functions here operate on abstract digital sequences and do not
describe or imply any wet-lab protocols.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from ..edit.events import EditEvent
from ..edit.simulate import SimulationContext, build_edit_dag
from .model import CasSystem, DigitalGenome as LegacyDigitalGenome, GuideRNA, TargetSite
from .physics import CRISPRPhysicsBase, create_crispr_physics


@dataclass
class CutEvent:
    """Represents a simulated cut event in a digital genome."""

    site: TargetSite
    cut_position: int
    guide: GuideRNA
    cas: CasSystem
    score: float


def find_candidate_sites(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_sites: Optional[int] = None,
    physics: Optional[CRISPRPhysicsBase] = None,
    use_gpu: bool = False,
) -> List[TargetSite]:
    """Return candidate sites using the new DigitalGenome view primitives."""

    physics_impl = physics or create_crispr_physics(cas, guide, use_gpu=use_gpu)
    results = physics_impl.score_sites(genome.sequences, max_sites=max_sites)
    return [result.site for result in results]


def simulate_cuts(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_events: Optional[int] = None,
    physics: Optional[CRISPRPhysicsBase] = None,
    use_gpu: bool = False,
) -> List[CutEvent]:
    """Simulate cut positions for the highest scoring candidate sites."""

    sites = find_candidate_sites(
        genome,
        cas,
        guide,
        max_sites=max_events,
        physics=physics,
        use_gpu=use_gpu,
    )
    events: List[CutEvent] = []
    for site in sites:
        if site.strand == 1:
            cut = site.start + cas.cut_offset
        else:
            cut = site.end - cas.cut_offset
        event = CutEvent(
            site=site,
            cut_position=max(site.start, min(site.end, cut)),
            guide=guide,
            cas=cas,
            score=site.on_target_score or 0.0,
        )
        events.append(event)
    return events


def rank_off_targets(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_candidates: int = 1000,
    physics: Optional[CRISPRPhysicsBase] = None,
    use_gpu: bool = False,
) -> List[TargetSite]:
    return find_candidate_sites(
        genome,
        cas,
        guide,
        max_sites=max_candidates,
        physics=physics,
        use_gpu=use_gpu,
    )
