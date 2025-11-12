"""
CRISPR sequence-level simulation engine for Helix.

All functions here operate on abstract digital sequences and do not
describe or imply any wet-lab protocols.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

from .. import bioinformatics
from .model import CasSystem, DigitalGenome, GuideRNA, PAMRule, TargetSite
from .pam import match_pam, reverse_complement_pattern


@dataclass
class CutEvent:
    """
    Represents a simulated cut event in a digital genome.

    This is a purely computational object used to reason about
    sequence changes.
    """

    site: TargetSite
    cut_position: int  # genomic coordinate within the chromosome
    guide: GuideRNA
    cas: CasSystem
    score: float


def _normalize_guide_sequence(guide: GuideRNA) -> str:
    """Return an uppercase guide sequence suitable for digital matching."""

    if not guide.sequence:
        raise ValueError("GuideRNA.sequence is required for simulation.")
    return bioinformatics.normalize_sequence(guide.sequence)


def _mismatch_count(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Sequences must have equal length for mismatch counting.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _pam_penalty(pattern: str, cas: CasSystem) -> float:
    pattern = pattern.upper()
    if not pattern:
        return 0.0
    ambiguous = sum(1 for char in pattern if char == "N")
    degeneracy = ambiguous / len(pattern)
    if degeneracy == 0:
        return 0.0
    # Scale by the configured penalty to keep the score within [0, 1].
    return degeneracy * (cas.weight_pam_penalty * 0.1)


def _score_site(mismatches: int, guide_len: int, cas: CasSystem, pam_penalty: float) -> float:
    mismatch_fraction = mismatches / guide_len if guide_len else 1.0
    penalty = mismatch_fraction * cas.weight_mismatch_penalty + pam_penalty
    return max(0.0, 1.0 - penalty)


def _select_pam_rules(rules: Sequence[PAMRule], guide: GuideRNA) -> Sequence[PAMRule]:
    if guide.pam:
        target = guide.pam.upper()
        filtered = [rule for rule in rules if rule.pattern.upper() == target]
        if filtered:
            return filtered
    return rules


def _scan_plus_strand(
    chrom: str,
    sequence: str,
    pattern: str,
    guide_seq: str,
    cas: CasSystem,
) -> Iterable[Tuple[TargetSite, int]]:
    pam_len = len(pattern)
    guide_len = len(guide_seq)
    if guide_len == 0 or pam_len == 0 or len(sequence) < guide_len + pam_len:
        return []
    pam = {"pattern": pattern}
    for pam_start in range(guide_len, len(sequence) - pam_len + 1):
        if not match_pam(sequence, pam, pam_start):
            continue
        guide_start = pam_start - guide_len
        target_seq = sequence[guide_start:pam_start]
        mismatches = _mismatch_count(target_seq, guide_seq)
        if mismatches > cas.max_mismatches:
            continue
        site = TargetSite(
            chrom=chrom,
            start=guide_start,
            end=pam_start,
            strand=1,
            sequence=target_seq,
        )
        yield site, mismatches


def _scan_minus_strand(
    chrom: str,
    sequence: str,
    pattern: str,
    guide_seq: str,
    cas: CasSystem,
) -> Iterable[Tuple[TargetSite, int]]:
    pam_len = len(pattern)
    guide_len = len(guide_seq)
    if guide_len == 0 or pam_len == 0 or len(sequence) < guide_len + pam_len:
        return []
    rc_pattern = reverse_complement_pattern(pattern)
    pam = {"pattern": rc_pattern}
    max_start = len(sequence) - (pam_len + guide_len)
    if max_start < 0:
        return []
    for pam_start in range(0, max_start + 1):
        if not match_pam(sequence, pam, pam_start):
            continue
        guide_start = pam_start + pam_len
        guide_end = guide_start + guide_len
        raw_seq = sequence[guide_start:guide_end]
        guide_oriented = bioinformatics.reverse_complement(raw_seq)
        mismatches = _mismatch_count(guide_oriented, guide_seq)
        if mismatches > cas.max_mismatches:
            continue
        site = TargetSite(
            chrom=chrom,
            start=guide_start,
            end=guide_end,
            strand=-1,
            sequence=guide_oriented,
        )
        yield site, mismatches


def _score_and_annotate_site(
    site: TargetSite,
    mismatches: int,
    guide_len: int,
    cas: CasSystem,
    pam_penalty: float,
) -> TargetSite:
    score = _score_site(mismatches, guide_len, cas, pam_penalty)
    site.on_target_score = round(score, 4)
    site.off_target_score = round(max(0.0, 1.0 - mismatches / guide_len), 4)
    site.pam_match_score = round(max(0.0, 1.0 - pam_penalty), 4)
    return site


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
    Returns an ordered list of candidate sites (highest on-target score first).
    """

    if max_sites is not None and max_sites <= 0:
        return []
    guide_seq = _normalize_guide_sequence(guide)
    pam_rules = _select_pam_rules(cas.pam_rules, guide)
    if not pam_rules:
        return []

    scored: List[Tuple[float, int, TargetSite]] = []
    for chrom, raw_seq in genome.sequences.items():
        sequence = bioinformatics.normalize_sequence(raw_seq)
        if not sequence:
            continue
        for rule in pam_rules:
            pattern = rule.pattern.upper()
            if not pattern:
                continue
            pam_penalty = _pam_penalty(pattern, cas)
            for site, mismatches in _scan_plus_strand(chrom, sequence, pattern, guide_seq, cas):
                annotated = _score_and_annotate_site(site, mismatches, len(guide_seq), cas, pam_penalty)
                scored.append((annotated.on_target_score or 0.0, mismatches, annotated))
            for site, mismatches in _scan_minus_strand(chrom, sequence, pattern, guide_seq, cas):
                annotated = _score_and_annotate_site(site, mismatches, len(guide_seq), cas, pam_penalty)
                scored.append((annotated.on_target_score or 0.0, mismatches, annotated))

    scored.sort(key=lambda item: (-item[0], item[1], item[2].chrom, item[2].start))
    if max_sites is not None:
        scored = scored[:max_sites]
    return [entry[2] for entry in scored]


def _chrom_length(genome: DigitalGenome, chrom: str) -> int:
    seq = genome.sequences.get(chrom, "")
    return len(seq)


def _compute_cut_position(site: TargetSite, cas: CasSystem, genome: DigitalGenome) -> int:
    chrom_len = _chrom_length(genome, site.chrom)
    boundary = site.end if site.strand == 1 else site.start
    offset = int(cas.cut_offset)
    if site.strand == 1:
        cut = boundary - offset
    else:
        cut = boundary + offset
    return max(0, min(chrom_len, cut))


def simulate_cuts(
    genome: DigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_events: Optional[int] = None,
) -> List[CutEvent]:
    """
    Simulate cut events in a digital genome for a given CRISPR system.

    Returns a list of potential cut events (sorted by decreasing score).
    """

    if max_events is not None and max_events <= 0:
        return []
    sites = find_candidate_sites(
        genome,
        cas,
        guide,
        max_sites=max_events,
    )
    events: List[CutEvent] = []
    for site in sites:
        score = site.on_target_score or 0.0
        cut_position = _compute_cut_position(site, cas, genome)
        events.append(
            CutEvent(
                site=site,
                cut_position=cut_position,
                guide=guide,
                cas=cas,
                score=score,
            )
        )
    return events


def rank_off_targets(
    genome: DigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    max_candidates: int = 1000,
) -> List[TargetSite]:
    """
    Compute an off-target ranking for a given guide.

    Sites are sorted by the same score used for simulate_cuts, favoring
    lower mismatch counts.
    """

    return find_candidate_sites(
        genome,
        cas,
        guide,
        max_sites=max_candidates,
    )
