"""
CRISPR sequence-level simulation engine for Helix.

This module acts as the public boundary for the CRISPR physics engine:
callers work with DigitalGenome / Guide / Cas dataclasses and rely on
these helpers to invoke the underlying physics backend. The heavy numeric
work happens inside helix.crispr.physics (which can later be swapped for
Numba/C++/Rust) while this module stays pure Python orchestration.

All functions here operate on abstract digital sequences and do not
describe or imply any wet-lab protocols.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

from .. import bioinformatics
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


@dataclass(frozen=True)
class EfficiencyTargetRequest:
    """Batch prediction request for a reference window + guide."""

    target_id: str
    reference_sequence: str
    guide: GuideRNA


@dataclass(frozen=True)
class EfficiencyPrediction:
    """Batch efficiency output."""

    target_id: str
    predicted_score: float
    top_site: TargetSite | None


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


def predict_efficiency_for_targets(
    cas: CasSystem,
    targets: Sequence[EfficiencyTargetRequest],
    *,
    max_sites: int = 1,
    use_gpu: bool = False,
) -> List[EfficiencyPrediction]:
    """
    Batch-oriented helper that scores many (sequence, guide) targets at once.

    The current implementation loops in Python for clarity; the signature is
    intentionally batch-shaped so that future native/vectorized engines can
    drop in without changing callers.
    """

    predictions: List[EfficiencyPrediction] = []
    for request in targets:
        normalized = bioinformatics.normalize_sequence(request.reference_sequence)
        if not normalized:
            predictions.append(
                EfficiencyPrediction(
                    target_id=request.target_id,
                    predicted_score=0.0,
                    top_site=None,
                )
            )
            continue
        genome = LegacyDigitalGenome({"target": normalized})
        sites = find_candidate_sites(
            genome,
            cas,
            request.guide,
            max_sites=max_sites,
            use_gpu=use_gpu,
        )
        top_site = sites[0] if sites else None
        score = float(top_site.on_target_score or 0.0) if top_site else 0.0
        predictions.append(
            EfficiencyPrediction(
                target_id=request.target_id,
                predicted_score=score,
                top_site=top_site,
            )
        )
    return predictions
