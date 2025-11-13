"""Prime editing rules for the edit DAG runtime."""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Tuple

from helix.edit.dag import EditNode
from helix.edit.events import EditEvent
from helix.edit.physics import edit_rule
from helix.edit.simulate import SimulationContext

from helix.crispr.model import DigitalGenome
from .model import PegRNA, PrimeEditor
from .simulator import apply_rtt_edit, locate_prime_target_site


def _get_prime_objects(ctx: SimulationContext) -> Tuple[DigitalGenome, PrimeEditor, PegRNA]:
    genome: DigitalGenome = ctx.extra["legacy_genome"]
    editor: PrimeEditor = ctx.extra["editor"]
    peg: PegRNA = ctx.extra["peg"]
    return genome, editor, peg


@edit_rule("prime.rtt_clean")
def prime_rtt_clean(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """First-stage prime editing rule: apply RTT at spacer match."""
    if node.metadata.get("stage", "root") != "root":
        return []

    genome, editor, peg = _get_prime_objects(ctx)
    site = locate_prime_target_site(genome, peg)
    if site is None:
        event = EditEvent(chrom="__none__", start=0, end=0, replacement="", metadata={"label": "prime_no_target"})
        return [(event, 0.0, {"stage": "no_target"})]

    chrom = site.chrom
    seq = node.genome_view.materialize_chrom(chrom)
    edited_seq = apply_rtt_edit(seq, site.start, site.end, peg.rtt)
    replacement = peg.rtt
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
            },
        )
    ]


@edit_rule("prime.no_edit")
def prime_no_edit(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """No-op branch for prime editing DAGs."""
    event = EditEvent(
        chrom="__none__",
        start=0,
        end=0,
        replacement="",
        metadata={"label": "prime_no_edit", "mechanism": "prime"},
    )
    stage = node.metadata.get("stage", "root")
    return [(event, -1.0, {"stage": stage})]
