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
from .physics import PrimePhysics, PrimeBranchDistribution


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


def apply_rtt_edit(seq: str, site_start: int, site_end: int, rtt: str) -> str:
    """Return a copy of `seq` with [site_start:site_end] replaced by rtt."""
    return seq[:site_start] + rtt + seq[site_end:]


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
    physics: Optional[PrimePhysics] = None,
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
    physics = physics or PrimePhysics.from_config(editor, peg)
    flank = max(len(peg.pbs), len(peg.rtt))
    window_start = max(0, site.start - flank)
    window_end = min(len(chrom_seq), site.end + flank)
    reference_window = _extract_window(chrom_seq, window_start, window_end)

    offset = (site.start - window_start) + editor.nick_to_edit_offset
    rtt_norm = _normalize_sequence(peg.rtt)
    intended_sequence = _apply_rtt(reference_window, offset, rtt_norm)

    spacer_seq = _normalize_sequence(peg.spacer)
    target_sequence = site.sequence or spacer_seq
    mismatch_count = physics.mismatch_count(target_sequence)
    window_gc = bioinformatics.gc_content(reference_window) if reference_window else None
    if editor.mismatch_tolerance is not None and mismatch_count > editor.mismatch_tolerance:
        distribution = PrimeBranchDistribution(
            efficiency=0.0,
            left=0.0,
            right=0.0,
            reanneal=0.0,
            indel=max(physics.min_score, editor.indel_bias),
            no_edit=1.0,
        )
    else:
        distribution = physics.branch_distribution(mismatch_count, local_gc=window_gc)

    outcomes: List[PrimeEditOutcome] = []
    efficiency_str = f"{distribution.efficiency:.4f}"
    left_logit = distribution.left
    right_logit = distribution.right
    reanneal_logit = distribution.reanneal
    indel_score = distribution.indel
    no_edit_score = distribution.no_edit

    if left_logit > 0:
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=intended_sequence,
                logit_score=left_logit,
                description="prime_left_flap",
                stage="flap",
                metadata={"branch": "left", "offset": str(offset), "physics_efficiency": efficiency_str},
            )
        )
    if right_logit > 0:
        keep = max(1, len(rtt_norm) // 2)
        partial = rtt_norm[:keep]
        right_sequence = reference_window[:offset] + partial + reference_window[offset + keep :]
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=right_sequence,
                logit_score=right_logit,
                description="prime_right_flap",
                stage="flap",
                metadata={
                    "branch": "right",
                    "partial_len": str(keep),
                    "physics_efficiency": efficiency_str,
                },
            )
        )
    if reanneal_logit > 0:
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=reference_window,
                logit_score=reanneal_logit,
                description="prime_reanneal",
                stage="reanneal",
                metadata={"branch": "reanneal", "physics_efficiency": efficiency_str},
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
                stage="error",
                metadata={
                    "branch": "deletion",
                    "deleted_len": str(deletion_end - offset),
                    "physics_efficiency": efficiency_str,
                },
            )
        )
        if len(reference_window) > offset:
            micro_del = reference_window[:offset] + reference_window[offset + 1 :]
            outcomes.append(
                PrimeEditOutcome(
                    site=site,
                    edited_sequence=micro_del,
                    logit_score=indel_score * 0.5,
                    description="micro_del",
                    stage="error",
                    metadata={"branch": "micro_del", "physics_efficiency": efficiency_str},
                )
            )
        micro_insert_base = rtt_norm[:1] or reference_window[offset : offset + 1]
        if micro_insert_base:
            micro_ins = reference_window[:offset] + micro_insert_base + reference_window[offset:]
            outcomes.append(
                PrimeEditOutcome(
                    site=site,
                    edited_sequence=micro_ins,
                    logit_score=indel_score * 0.5,
                    description="micro_ins",
                    stage="error",
                    metadata={"branch": "micro_ins", "physics_efficiency": efficiency_str},
                )
            )
    if no_edit_score > 0:
        outcomes.append(
            PrimeEditOutcome(
                site=site,
                edited_sequence=reference_window,
                logit_score=no_edit_score,
                description="no_edit",
                stage="root",
                metadata={"branch": "no_edit", "physics_efficiency": efficiency_str},
            )
        )

    outcomes.sort(key=lambda outcome: outcome.logit_score, reverse=True)
    return outcomes[:max(1, max_outcomes)]
