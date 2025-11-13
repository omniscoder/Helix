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
from ..genome.digital import DigitalGenome, DigitalGenomeView
from .model import CasSystem, DigitalGenome as LegacyDigitalGenome, GuideRNA, TargetSite


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
    pam_pattern: str,
    cas: CasSystem,
) -> Iterable[Tuple[TargetSite, int]]:
    seq = seq.upper()
    guide_len = len(guide_seq)
    if guide_len == 0 or len(seq) < guide_len:
        return []

    candidates: List[Tuple[TargetSite, int]] = []

    # plus strand: guide immediately precedes PAM
    pam_len = len(pam_pattern)
    if pam_len > 0 and len(seq) >= guide_len + pam_len:
        for pos in range(guide_len, len(seq) - pam_len + 1):
            target = seq[pos - guide_len : pos]
            mismatches = _mismatch_count(target, guide_seq)
            if mismatches > cas.max_mismatches:
                continue
            site = TargetSite(
                chrom=chrom,
                start=pos - guide_len,
                end=pos,
                strand=1,
                sequence=target,
            )
            candidates.append((site, mismatches))

    # minus strand: scan reverse complement and map coordinates
    rc_seq = bioinformatics.reverse_complement(seq)
    for pos in range(guide_len, len(rc_seq) - pam_len + 1):
        target = rc_seq[pos - guide_len : pos]
        mismatches = _mismatch_count(target, guide_seq)
        if mismatches > cas.max_mismatches:
            continue
        # map back to original coordinates
        rc_start = pos - guide_len
        rc_end = pos
        start = len(seq) - rc_end
        end = len(seq) - rc_start
        site = TargetSite(
            chrom=chrom,
            start=start,
            end=end,
            strand=-1,
            sequence=bioinformatics.reverse_complement(target),
        )
        candidates.append((site, mismatches))

    return candidates


def find_candidate_sites(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_sites: Optional[int] = None,
) -> List[TargetSite]:
    """Return candidate sites using the new DigitalGenome view primitives."""

    guide_seq = _normalize_sequence(guide.sequence)
    pam_pattern = cas.pam_rules[0].pattern.upper() if cas.pam_rules else "NGG"
    digital = DigitalGenome(sequences=genome.sequences)
    view = digital.view()
    scored: List[Tuple[float, TargetSite]] = []

    for chrom, sequence in genome.sequences.items():
        for site, mismatches in _scan_chromosome(chrom, sequence, guide_seq, pam_pattern, cas):
            score = max(0.0, 1.0 - (mismatches / len(guide_seq)))
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
) -> List[CutEvent]:
    """Simulate cut positions for the highest scoring candidate sites."""

    sites = find_candidate_sites(genome, cas, guide, max_sites=max_events)
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
) -> List[TargetSite]:
    return find_candidate_sites(genome, cas, guide, max_sites=max_candidates)
