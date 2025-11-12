"""
Prime editing sequence-level simulation for Helix.

All operations here work on digital sequences only and are not
wet-lab protocols.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

from helix.crispr.model import DigitalGenome, TargetSite
from helix import bioinformatics

from .model import PegRNA, PrimeEditOutcome, PrimeEditor


def _normalize_sequence(value: str, *, allow_ambiguous: bool = True) -> str:
    return bioinformatics.normalize_sequence(value, allow_ambiguous=allow_ambiguous)


def _find_exact_match(sequence: str, query: str) -> List[int]:
    positions: List[int] = []
    start = 0
    while True:
        idx = sequence.find(query, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def _best_approximate_match(sequence: str, query: str, limit: int) -> Optional[Tuple[int, int]]:
    best_idx: Optional[int] = None
    best_mismatches: Optional[int] = None
    max_start = max(0, min(len(sequence) - len(query), limit))
    for idx in range(max_start + 1):
        window = sequence[idx : idx + len(query)]
        if len(window) != len(query):
            continue
        mismatches = sum(a != b for a, b in zip(window, query))
        if best_mismatches is None or mismatches < best_mismatches:
            best_idx = idx
            best_mismatches = mismatches
            if mismatches == 0:
                break
    if best_idx is None or best_mismatches is None:
        return None
    return best_idx, best_mismatches


def locate_prime_target_site(
    genome: DigitalGenome,
    peg: PegRNA,
    *,
    search_window: int = 200,
) -> Optional[TargetSite]:
    """
    Identify a primary digital target site for a pegRNA.

    Returns None if no plausible site is found.
    """

    spacer = _normalize_sequence(peg.spacer)
    if not spacer:
        return None
    rc_spacer = bioinformatics.reverse_complement(spacer)

    for chrom, raw_seq in genome.sequences.items():
        sequence = _normalize_sequence(raw_seq)
        if not sequence:
            continue
        for idx in _find_exact_match(sequence, spacer):
            site = TargetSite(
                chrom=chrom,
                start=idx,
                end=idx + len(spacer),
                strand=1,
                sequence=spacer,
                on_target_score=1.0,
            )
            return site
        for idx in _find_exact_match(sequence, rc_spacer):
            site = TargetSite(
                chrom=chrom,
                start=idx,
                end=idx + len(spacer),
                strand=-1,
                sequence=spacer,
                on_target_score=1.0,
            )
            return site

    # Fallback: approximate match within the provided window.
    limit = max(search_window, 0)
    best_site: Optional[TargetSite] = None
    best_mismatches: Optional[int] = None
    for chrom, raw_seq in genome.sequences.items():
        sequence = _normalize_sequence(raw_seq)
        if not sequence:
            continue
        approx = _best_approximate_match(sequence, spacer, limit)
        if approx:
            idx, mismatches = approx
            if best_mismatches is None or mismatches < best_mismatches:
                best_mismatches = mismatches
                best_site = TargetSite(
                    chrom=chrom,
                    start=idx,
                    end=idx + len(spacer),
                    strand=1,
                    sequence=sequence[idx : idx + len(spacer)],
                    on_target_score=max(0.0, 1.0 - mismatches / len(spacer)),
                )
        approx_rc = _best_approximate_match(sequence, rc_spacer, limit)
        if approx_rc:
            idx, mismatches = approx_rc
            if best_mismatches is None or mismatches < best_mismatches:
                seq_slice = sequence[idx : idx + len(spacer)]
                best_mismatches = mismatches
                best_site = TargetSite(
                    chrom=chrom,
                    start=idx,
                    end=idx + len(spacer),
                    strand=-1,
                    sequence=bioinformatics.reverse_complement(seq_slice),
                    on_target_score=max(0.0, 1.0 - mismatches / len(spacer)),
                )
    return best_site


def _extract_window(sequence: str, start: int, end: int) -> str:
    return sequence[start:end]


def _apply_rtt(reference: str, offset: int, template: str) -> str:
    if not reference:
        return template
    offset = max(0, min(len(reference), offset))
    window = list(reference)
    for idx, base in enumerate(template):
        pos = offset + idx
        if pos >= len(window):
            break
        window[pos] = base
    return "".join(window)


def simulate_prime_edit(
    genome: DigitalGenome,
    editor: PrimeEditor,
    peg: PegRNA,
    *,
    max_outcomes: int = 16,
) -> List[PrimeEditOutcome]:
    """
    Simulate prime editing outcomes in a digital genome.

    Returns hypothetical outcomes with heuristic logit scores.
    """

    if max_outcomes <= 0:
        return []
    site = locate_prime_target_site(genome, peg)
    if site is None:
        return []

    chrom_seq = _normalize_sequence(genome.sequences[site.chrom])
    flank = max(len(peg.pbs), len(peg.rtt))
    window_start = max(0, site.start - flank)
    window_end = min(len(chrom_seq), site.end + flank)
    reference_window = _extract_window(chrom_seq, window_start, window_end)

    offset = (site.start - window_start) + editor.nick_to_edit_offset
    intended_sequence = _apply_rtt(reference_window, offset, _normalize_sequence(peg.rtt))

    spacer_seq = _normalize_sequence(peg.spacer)
    target_sequence = site.sequence or spacer_seq
    mismatch_count = sum(a != b for a, b in zip(target_sequence, spacer_seq))
    efficiency = max(0.0, editor.efficiency_scale - 0.1 * mismatch_count)
    if mismatch_count > editor.mismatch_tolerance:
        efficiency = 0.0
    no_edit_score = max(0.0, 1.0 - efficiency)
    indel_score = max(0.0, editor.indel_bias)

    outcomes: List[PrimeEditOutcome] = []
    if efficiency > 0:
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=intended_sequence,
                logit_score=efficiency,
                description="intended_edit",
            )
        )
    if indel_score > 0:
        deletion_end = min(len(reference_window), offset + len(peg.pbs))
        indel_sequence = reference_window[:offset] + reference_window[deletion_end:]
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=indel_sequence,
                logit_score=indel_score,
                description="indel_loss",
            )
        )
    if no_edit_score > 0:
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=reference_window,
                logit_score=no_edit_score,
                description="no_edit",
            )
        )

    outcomes.sort(key=lambda outcome: outcome.logit_score, reverse=True)
    return outcomes[:max(1, max_outcomes)]
