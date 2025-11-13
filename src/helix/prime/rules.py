"""Prime editing rules for the edit DAG runtime."""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Tuple

from helix import bioinformatics
from helix.edit.dag import EditNode
from helix.edit.events import EditEvent
from helix.edit.physics import edit_rule
from helix.edit.simulate import SimulationContext

from helix.crispr.model import DigitalGenome
from .model import PegRNA, PrimeEditor
from .physics import PrimeBranchDistribution, PrimePhysics
from .simulator import apply_rtt_edit, locate_prime_target_site


def _get_prime_objects(ctx: SimulationContext) -> Tuple[DigitalGenome, PrimeEditor, PegRNA]:
    genome: DigitalGenome = ctx.extra["legacy_genome"]
    editor: PrimeEditor = ctx.extra["editor"]
    peg: PegRNA = ctx.extra["peg"]
    return genome, editor, peg


def _get_prime_physics(ctx: SimulationContext) -> PrimePhysics | None:
    physics = ctx.extra.get("prime_physics")
    if isinstance(physics, PrimePhysics):
        return physics
    return None


def _local_gc_fraction(sequence: str, start: int, end: int, peg: PegRNA) -> float | None:
    flank = max(len(peg.pbs), len(peg.rtt))
    window_start = max(0, start - flank)
    window_end = min(len(sequence), end + flank)
    if window_start >= window_end:
        return None
    window = sequence[window_start:window_end]
    if not window:
        return None
    return bioinformatics.gc_content(window)


@edit_rule("prime.rtt_clean")
def prime_rtt_clean(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """First-stage prime editing rule: apply RTT at spacer match."""
    if node.metadata.get("stage", "root") != "root":
        return []

    genome, editor, peg = _get_prime_objects(ctx)
    physics = _get_prime_physics(ctx)
    site = locate_prime_target_site(genome, peg)
    if site is None:
        event = EditEvent(chrom="__none__", start=0, end=0, replacement="", metadata={"label": "prime_no_target"})
        return [(event, 0.0, {"stage": "no_target"})]

    chrom = site.chrom
    seq = node.genome_view.materialize_chrom(chrom)
    edited_seq = apply_rtt_edit(seq, site.start, site.end, peg.rtt)
    replacement = peg.rtt
    mismatch_count = 0
    local_gc = _local_gc_fraction(seq, site.start, site.end, peg)
    if physics:
        target_slice = seq[site.start : site.end]
        if site.strand == -1:
            target_slice = bioinformatics.reverse_complement(target_slice)
        mismatch_count = physics.mismatch_count(target_slice)
        log_prob = physics.log_efficiency(mismatch_count, local_gc=local_gc)
    else:
        log_prob = math.log(max(editor.efficiency_scale or 1e-6, 1e-6))

    event = EditEvent(
        chrom=chrom,
        start=site.start,
        end=site.end,
        replacement=replacement,
        metadata={
            "label": "prime_rtt_clean",
            "mechanism": "prime",
            "stage": "prime_rtt",
            "site_start": site.start,
            "site_end": site.end,
        },
    )
    return [
        (
            event,
            log_prob,
            {
                "stage": "prime_rtt",
                "site_chrom": site.chrom,
                "site_start": site.start,
                "site_end": site.end,
                "physics_mismatch_count": mismatch_count,
                "physics_local_gc": local_gc if local_gc is not None else None,
            },
        )
    ]


@edit_rule("prime.no_edit")
def prime_no_edit(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """No-op branch for prime editing DAGs."""
    if node.metadata.get("stage", "root") != "root":
        return []
    event = EditEvent(
        chrom="__none__",
        start=0,
        end=0,
        replacement="",
        metadata={"label": "prime_no_edit", "mechanism": "prime"},
    )
    stage = node.metadata.get("stage", "root")
    return [(event, -1.0, {"stage": stage})]


@edit_rule("prime.flap_resolution")
def prime_flap_resolution(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """Resolve flap outcomes after an RTT application."""
    if node.metadata.get("stage") != "prime_rtt":
        return []
    chrom = node.metadata.get("site_chrom")
    start = node.metadata.get("site_start")
    end = node.metadata.get("site_end")
    if chrom is None or start is None or end is None:
        return []
    chrom = str(chrom)
    start = int(start)
    end = int(end)
    seq = node.genome_view.materialize_chrom(chrom)
    physics = _get_prime_physics(ctx)
    mismatch_count = int(node.metadata.get("physics_mismatch_count", 0))
    local_gc_raw = node.metadata.get("physics_local_gc")
    local_gc = float(local_gc_raw) if isinstance(local_gc_raw, (int, float)) else None
    if physics:
        distribution = physics.branch_distribution(mismatch_count, local_gc=local_gc)
    else:
        distribution = PrimeBranchDistribution(
            efficiency=0.5,
            left=0.25,
            right=0.25,
            reanneal=0.1,
            indel=0.1,
            no_edit=0.3,
        )

    proposals: List[Tuple[EditEvent, float, Dict]] = []

    def _log_weight(value: float) -> float:
        return math.log(max(value, 1e-6))

    if distribution.left > 0:
        left_event = EditEvent(
            chrom=chrom,
            start=max(0, start - 1),
            end=end,
            replacement=seq[max(0, start - 1) : end],
            metadata={"label": "prime_flap_left", "mechanism": "prime"},
        )
        proposals.append(
            (
                left_event,
                _log_weight(distribution.left),
                {
                    "stage": "repaired",
                    "branch": "flap_left",
                    "physics_efficiency": distribution.efficiency,
                },
            )
        )

    if distribution.right > 0:
        right_event = EditEvent(
            chrom=chrom,
            start=start,
            end=min(len(seq), end + 1),
            replacement=seq[start : min(len(seq), end + 1)],
            metadata={"label": "prime_flap_right", "mechanism": "prime"},
        )
        proposals.append(
            (
                right_event,
                _log_weight(distribution.right),
                {
                    "stage": "repaired",
                    "branch": "flap_right",
                    "physics_efficiency": distribution.efficiency,
                },
            )
        )

    if distribution.reanneal > 0:
        reanneal_event = EditEvent(
            chrom=chrom,
            start=start,
            end=end,
            replacement=seq[start:end],
            metadata={"label": "prime_reanneal", "mechanism": "prime"},
        )
        proposals.append(
            (
                reanneal_event,
                _log_weight(distribution.reanneal),
                {
                    "stage": "repaired",
                    "branch": "reanneal",
                    "physics_efficiency": distribution.efficiency,
                },
            )
        )

    if distribution.indel > 0:
        deletion_end = min(len(seq), end + len(seq[start:end]) // 2)
        indel_event = EditEvent(
            chrom=chrom,
            start=start,
            end=deletion_end,
            replacement="",
            metadata={"label": "prime_indel", "mechanism": "prime"},
        )
        proposals.append(
            (
                indel_event,
                _log_weight(distribution.indel),
                {"stage": "error", "physics_efficiency": distribution.efficiency},
            )
        )

    if distribution.no_edit > 0:
        no_edit_event = EditEvent(
            chrom="__none__",
            start=0,
            end=0,
            replacement="",
            metadata={"label": "prime_no_edit", "mechanism": "prime"},
        )
        proposals.append(
            (
                no_edit_event,
                _log_weight(distribution.no_edit),
                {"stage": "no_edit", "physics_efficiency": distribution.efficiency},
            )
        )

    return proposals
