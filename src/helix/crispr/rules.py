"""CRISPR edit rules plugged into the edit DAG runtime."""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Tuple

from helix.edit.dag import EditNode
from helix.edit.events import EditEvent
from helix.edit.physics import edit_rule
from helix.edit.simulate import SimulationContext

from .model import CasSystem, DigitalGenome, GuideRNA
from .physics import CRISPRPhysics
from .simulator import find_candidate_sites


def _get_context_objects(ctx: SimulationContext) -> Tuple[DigitalGenome, CasSystem, GuideRNA]:
    genome: DigitalGenome = ctx.extra["legacy_genome"]
    cas: CasSystem = ctx.extra["cas"]
    guide: GuideRNA = ctx.extra["guide"]
    return genome, cas, guide


def _get_physics(ctx: SimulationContext) -> CRISPRPhysics | None:
    physics = ctx.extra.get("physics")
    if isinstance(physics, CRISPRPhysics):
        return physics
    return None


@edit_rule("crispr.clean_cut")
def crispr_clean_cut(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """First-stage CRISPR rule: propose clean-cut events."""
    if node.metadata.get("stage", "root") != "root":
        return []

    genome, cas, guide = _get_context_objects(ctx)
    max_sites = ctx.extra.get("max_sites", 3) or 3
    no_edit_prob = ctx.extra.get("no_edit_prob", 0.1)

    physics = _get_physics(ctx)
    candidates = find_candidate_sites(genome, cas, guide, max_sites=max_sites, physics=physics)
    proposals: List[Tuple[EditEvent, float, Dict]] = []

    if not candidates:
        event = EditEvent(chrom="__none__", start=0, end=0, replacement="", metadata={"label": "no_target"})
        proposals.append((event, 0.0, {"stage": "no_target"}))
        return proposals

    for site in candidates:
        cut_pos = site.start + cas.cut_offset
        cut_pos = max(site.start, min(site.end, cut_pos))
        log_prob = math.log(max(site.on_target_score or 1e-6, 1e-6))
        event = EditEvent(
            chrom=site.chrom,
            start=cut_pos,
            end=cut_pos,
            replacement="",
            metadata={
                "label": "clean_cut",
                "mechanism": "crispr",
                "site_chrom": site.chrom,
                "site_start": site.start,
                "site_end": site.end,
                "cut_position": cut_pos,
                "strand": site.strand,
                "stage": "cut",
            },
        )
        proposals.append(
            (
                event,
                log_prob,
                {
                    "stage": "cut",
                    "site_chrom": site.chrom,
                    "site_start": site.start,
                    "site_end": site.end,
                    "cut_position": cut_pos,
                    "strand": site.strand,
                },
            )
        )

    if no_edit_prob > 0:
        event = EditEvent(
            chrom="__none__",
            start=0,
            end=0,
            replacement="",
            metadata={"label": "no_edit", "mechanism": "crispr", "stage": "no_edit"},
        )
        proposals.append((event, math.log(max(no_edit_prob, 1e-6)), {"stage": "no_edit"}))

    return proposals


@edit_rule("crispr.indel_branch")
def crispr_indel_branch(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """Second-stage CRISPR rule: branch into indel/intended repair outcomes."""
    if node.metadata.get("stage") != "cut":
        return []

    chrom = node.metadata.get("site_chrom")
    cut_pos = node.metadata.get("cut_position")
    if chrom is None or cut_pos is None:
        return []

    seq = node.genome_view.materialize_chrom(chrom)
    window = ctx.extra.get("indel_window", 3)
    start = max(0, cut_pos - window)
    end = min(len(seq), cut_pos + window)
    if start > end:
        start, end = end, start

    original = seq[start:end]
    proposals: List[Tuple[EditEvent, float, Dict]] = []

    intended = original[::-1] if original else ""
    intended_event = EditEvent(
        chrom=chrom,
        start=start,
        end=end,
        replacement=intended,
        metadata={"label": "intended_edit", "mechanism": "crispr", "stage": "repaired"},
    )
    proposals.append((intended_event, 0.0, {"stage": "repaired"}))

    del_start = cut_pos
    del_end = min(cut_pos + 1, len(seq))
    indel_event = EditEvent(
        chrom=chrom,
        start=del_start,
        end=del_end,
        replacement="",
        metadata={"label": "indel", "mechanism": "crispr", "stage": "repaired"},
    )
    proposals.append((indel_event, -0.5, {"stage": "repaired"}))

    return proposals


@edit_rule("crispr.no_edit")
def crispr_no_edit(node: EditNode, ctx: SimulationContext) -> Iterable[Tuple[EditEvent, float, Dict]]:
    """Generic no-op branch to keep 'no edit' paths alive."""
    event = EditEvent(
        chrom="__none__",
        start=0,
        end=0,
        replacement="",
        metadata={"label": "no_edit_branch", "mechanism": "crispr"},
    )
    stage = node.metadata.get("stage", "root")
    return [(event, -1.0, {"stage": stage})]
