"""Simple deterministic PCR simulation engine for Helix Studio."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .session import PCRConfig, PCRCycleStats, PCRResult

BP_MASS_DA = 650.0
DA_TO_NG = 1.66053906660e-15  # 1 dalton in nanograms


def _revcomp(seq: str) -> str:
    complement = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(complement)[::-1]


@dataclass(frozen=True)
class AmpliconInfo:
    start: int
    end: int
    length: int


def find_amplicon(template_seq: str, fwd_primer: str, rev_primer: str) -> AmpliconInfo:
    template = template_seq.upper()
    fwd = fwd_primer.upper()
    rev_rc = _revcomp(rev_primer.upper())

    start = template.find(fwd)
    end_rc = template.rfind(rev_rc)
    if start == -1 or end_rc == -1:
        raise ValueError("Primers do not bind template to produce an amplicon")
    end = end_rc + len(rev_rc)
    if end <= start:
        raise ValueError("Reverse primer binds upstream of forward primer")
    return AmpliconInfo(start=start, end=end, length=end - start)


def _copies_to_mass_ng(copies: float, length_bp: int) -> float:
    mol_mass_da = length_bp * BP_MASS_DA
    return copies * mol_mass_da * DA_TO_NG


def simulate_pcr(config: PCRConfig, rng: Optional[object] = None) -> PCRResult:
    amplicon = find_amplicon(config.template_seq, config.fwd_primer, config.rev_primer)
    length = amplicon.length
    result = PCRResult(
        config=config,
        amplicon_length=length,
        amplicon_start=amplicon.start,
        amplicon_end=amplicon.end,
    )
    copies = float(max(config.initial_copies, 1))

    for cycle in range(1, config.cycles + 1):
        copies *= (1.0 + max(config.efficiency, 0.0))

        per_base = max(config.base_error_rate, 0.0)
        per_molecule = 1.0 - (1.0 - per_base) ** length
        mutated_fraction = 1.0 - (1.0 - per_molecule) ** cycle

        sub_rate = mutated_fraction * (1.0 - config.indel_fraction)
        indel_rate = mutated_fraction * config.indel_fraction
        mass_ng = _copies_to_mass_ng(copies, length)

        result.cycles.append(
            PCRCycleStats(
                cycle_index=cycle,
                amplicon_copies=copies,
                amplicon_mass_ng=mass_ng,
                cumulative_mutation_rate=mutated_fraction,
                sub_rate=sub_rate,
                indel_rate=indel_rate,
            )
        )

    if result.cycles:
        final_cycle = result.cycles[-1]
        result.final_mutation_rate = final_cycle.cumulative_mutation_rate
        result.final_amplicon_copies = final_cycle.amplicon_copies
        result.final_amplicon_mass_ng = final_cycle.amplicon_mass_ng
    return result


__all__ = ["simulate_pcr", "find_amplicon", "AmpliconInfo"]
