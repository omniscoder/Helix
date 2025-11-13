"""PCR edit-dag rules ("physics")."""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

from helix.edit.dag import EditNode
from helix.edit.events import EditEvent
from helix.edit.physics import edit_rule
from helix.edit.simulate import SimulationContext
from helix.genome.digital import DigitalGenome

from .model import PCRConfig, PrimerPair


def _find_k_mismatch_hits(seq: str, query: str, k: int) -> List[int]:
    seq = seq.upper()
    query = query.upper()
    hits: List[int] = []
    n, m = len(seq), len(query)
    if m == 0 or n < m:
        return hits
    for i in range(0, n - m + 1):
        mismatches = 0
        for j in range(m):
            if seq[i + j] != query[j]:
                mismatches += 1
                if mismatches > k:
                    break
        if mismatches <= k:
            hits.append(i)
    return hits


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(comp)[::-1]


def _pcr_context(ctx: SimulationContext) -> tuple[DigitalGenome, PrimerPair, PCRConfig]:
    genome: DigitalGenome = ctx.extra["core_genome"]  # type: ignore[assignment]
    pair: PrimerPair = ctx.extra["primer_pair"]  # type: ignore[assignment]
    cfg: PCRConfig = ctx.extra["config"]  # type: ignore[assignment]
    return genome, pair, cfg


def _limit_events(events: List[Tuple[EditEvent, float, Dict]], limit: int) -> List[Tuple[EditEvent, float, Dict]]:
    if len(events) <= limit:
        return events
    # Prefer amplicons near the midpoint of allowable lengths
    def _score(ev: EditEvent) -> float:
        length = ev.metadata.get("amplicon_length", 0)
        return abs(length - 200.0)

    return sorted(events, key=lambda entry: _score(entry[0]))[:limit]


@edit_rule("pcr.primer_binding")
def pcr_primer_binding(node: EditNode, context: SimulationContext):
    """
    Identify candidate binding sites for a primer pair.
    """

    stage = node.metadata.get("stage", "root")
    if stage not in ("root", "no_product"):
        return []

    genome, pair, cfg = _pcr_context(context)
    events: List[Tuple[EditEvent, float, Dict]] = []
    for chrom, seq in genome.sequences.items():
        forward_hits = _find_k_mismatch_hits(seq, pair.forward.sequence, pair.forward.max_mismatches)
        rc_seq = _revcomp(seq)
        reverse_hits_rc = _find_k_mismatch_hits(rc_seq, pair.reverse.sequence, pair.reverse.max_mismatches)
        reverse_hits: List[int] = []
        rev_len = len(pair.reverse.sequence)
        for hit in reverse_hits_rc:
            start = len(seq) - (hit + rev_len)
            reverse_hits.append(start)

        for f_start in forward_hits:
            for r_start in reverse_hits:
                if r_start <= f_start:
                    continue
                amplicon_end = r_start + rev_len
                length = amplicon_end - f_start
                if not (cfg.min_amplicon_length <= length <= cfg.max_amplicon_length):
                    continue
                logp = -math.log(1.0 + abs(length - 200.0) / 200.0)
                ev = EditEvent(
                    chrom=chrom,
                    start=f_start,
                    end=f_start,
                    replacement="",
                    metadata={
                        "label": "amplicon_template",
                        "mechanism": "pcr",
                        "stage": "binding",
                        "amplicon_chrom": chrom,
                        "amplicon_start": f_start,
                        "amplicon_end": amplicon_end,
                        "amplicon_length": length,
                        "primer_pair": pair.name,
                    },
                )
                events.append((ev, logp, {"stage": "binding"}))

    if not events:
        ev = EditEvent(
            chrom="__none__",
            start=0,
            end=0,
            replacement="",
            metadata={"label": "no_product", "mechanism": "pcr", "stage": "no_product"},
        )
        return [(ev, math.log(1e-12), {"stage": "no_product"})]

    events = _limit_events(events, cfg.max_amplicons)
    base_logp = -math.log(len(events))
    normalized = []
    for ev, logp, meta in events:
        normalized.append((ev, logp + base_logp, meta))
    return normalized


@edit_rule("pcr.amplify_cycle")
def pcr_amplify_cycle(node: EditNode, context: SimulationContext):
    """
    Model one amplification cycle by boosting abundance probability.
    """

    if node.metadata.get("stage") not in ("binding", "cycled"):
        return []

    cfg: PCRConfig = context.extra["config"]  # type: ignore[assignment]
    eff = max(0.0, min(1.0, cfg.per_cycle_efficiency))
    logp_delta = math.log(1.0 + eff)
    cycle = int(node.metadata.get("cycle", 0)) + 1
    ev = EditEvent(
        chrom="__none__",
        start=0,
        end=0,
        replacement="",
        metadata={"label": "amplify_cycle", "mechanism": "pcr", "stage": "cycled"},
    )
    return [(ev, logp_delta, {"stage": "cycled", "cycle": cycle})]


@edit_rule("pcr.error_branch")
def pcr_error_branch(node: EditNode, context: SimulationContext):
    """
    Introduce a substitution error within the amplicon.
    """

    stage = node.metadata.get("stage")
    if stage not in ("binding", "cycled", "amplicon"):
        return []

    chrom = node.metadata.get("amplicon_chrom")
    start = node.metadata.get("amplicon_start")
    end = node.metadata.get("amplicon_end")
    if chrom is None or start is None or end is None:
        return []
    start = int(start)
    end = int(end)
    if end <= start:
        return []

    seq = node.genome_view.materialize_chrom(chrom)
    if end > len(seq):
        end = len(seq)
    if end <= start:
        return []
    pos = context.rng.randint(start, end - 1)
    original = seq[pos].upper()
    bases = [base for base in "ACGT" if base != original]
    if not bases:
        return []
    replacement = context.rng.choice(bases)

    cfg: PCRConfig = context.extra["config"]  # type: ignore[assignment]
    error_rate = max(cfg.error_rate * max(cfg.cycles, 1), 1e-9)
    logp_delta = math.log(error_rate) - math.log(max(1, end - start))

    ev = EditEvent(
        chrom=chrom,
        start=pos,
        end=pos + 1,
        replacement=replacement,
        metadata={
            "label": "pcr_error",
            "mechanism": "pcr",
            "stage": "error",
            "amplicon_chrom": chrom,
            "amplicon_start": start,
            "amplicon_end": end,
            "amplicon_length": end - start,
            "error_position": pos,
        },
    )
    return [(ev, logp_delta, {"stage": "error"})]
