"""
CRISPR sequence-level simulation engine for Helix.

All functions here operate on abstract digital sequences and do not
describe or imply any wet-lab protocols.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

from .. import bioinformatics
from ..edit.events import EditEvent
from ..edit.simulate import SimulationContext, build_edit_dag
from .kmm import k_mismatch_positions
from .model import CasSystem, DigitalGenome as LegacyDigitalGenome, GuideRNA, TargetSite
from .physics import CRISPRPhysics


@dataclass
class CutEvent:
    """Represents a simulated cut event in a digital genome."""

    site: TargetSite
    cut_position: int
    guide: GuideRNA
    cas: CasSystem
    score: float


def _normalize_sequence(seq: str) -> str:
    return bioinformatics.normalize_sequence(seq)


def _mismatch_count(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _scan_chromosome(
    chrom: str,
    seq: str,
    guide_seq: str,
    cas: CasSystem,
    physics: CRISPRPhysics,
) -> Iterable[Tuple[TargetSite, float]]:
    seq = seq.upper()
    guide_seq = guide_seq.upper()
    guide_len = len(guide_seq)
    if guide_len == 0 or len(seq) < guide_len:
        return []

    candidates: List[Tuple[TargetSite, float]] = []
    pam_len = len(physics.pam_pattern or "")
    max_mismatch = cas.max_mismatches if cas.max_mismatches is not None else guide_len

    # plus strand: guide precedes PAM
    forward_hits: Sequence[int]
    if len(seq) >= guide_len:
        if cas.max_mismatches is not None:
            forward_hits = k_mismatch_positions(seq, guide_seq, max_mismatch)
        else:
            forward_hits = range(0, len(seq) - guide_len + 1)
        for start in forward_hits:
            pam_start = start + guide_len
            if pam_len and pam_start + pam_len > len(seq):
                continue
            pam_seq = seq[pam_start : pam_start + pam_len] if pam_len else ""
            target = seq[start : start + guide_len]
            mismatches = _mismatch_count(target, guide_seq)
            if cas.max_mismatches is not None and mismatches > cas.max_mismatches:
                continue
            score = physics.score_site(target, pam_seq, strand=1)
            if score <= 0:
                continue
            site = TargetSite(
                chrom=chrom,
                start=start,
                end=start + guide_len,
                strand=1,
                sequence=target,
            )
            candidates.append((site, score))

    # minus strand
    rc_seq = bioinformatics.reverse_complement(seq)
    if len(rc_seq) >= guide_len:
        if cas.max_mismatches is not None:
            reverse_hits = k_mismatch_positions(rc_seq, guide_seq, max_mismatch)
        else:
            reverse_hits = range(0, len(rc_seq) - guide_len + 1)
        for start in reverse_hits:
            pam_start = start + guide_len
            if pam_len and pam_start + pam_len > len(rc_seq):
                continue
            pam_rc = rc_seq[pam_start : pam_start + pam_len] if pam_len else ""
            target_rc = rc_seq[start : start + guide_len]
            mismatches = _mismatch_count(target_rc, guide_seq)
            if cas.max_mismatches is not None and mismatches > cas.max_mismatches:
                continue
            score = physics.score_site(target_rc, pam_rc, strand=-1)
            if score <= 0:
                continue
            rc_start = start
            rc_end = start + guide_len
            orig_start = len(seq) - rc_end
            orig_end = len(seq) - rc_start
            site = TargetSite(
                chrom=chrom,
                start=orig_start,
                end=orig_end,
                strand=-1,
                sequence=bioinformatics.reverse_complement(target_rc),
            )
            candidates.append((site, score))

    return candidates


def find_candidate_sites(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_sites: Optional[int] = None,
    physics: Optional[CRISPRPhysics] = None,
) -> List[TargetSite]:
    """Return candidate sites using the new DigitalGenome view primitives."""

    guide_seq = _normalize_sequence(guide.sequence)
    physics = physics or CRISPRPhysics.from_system(cas, guide)
    scored: List[Tuple[float, TargetSite]] = []

    for chrom, sequence in genome.sequences.items():
        for site, score in _scan_chromosome(chrom, sequence, guide_seq, cas, physics):
            site.on_target_score = score
            scored.append((score, site))

    scored.sort(key=lambda pair: pair[0], reverse=True)
    if max_sites is not None:
        scored = scored[:max_sites]
    return [site for _, site in scored]


def simulate_cuts(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_events: Optional[int] = None,
    physics: Optional[CRISPRPhysics] = None,
) -> List[CutEvent]:
    """Simulate cut positions for the highest scoring candidate sites."""

    sites = find_candidate_sites(genome, cas, guide, max_sites=max_events, physics=physics)
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
    physics: Optional[CRISPRPhysics] = None,
) -> List[TargetSite]:
    return find_candidate_sites(genome, cas, guide, max_sites=max_candidates, physics=physics)
