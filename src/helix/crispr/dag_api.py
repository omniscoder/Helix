"""CRISPR-specific helpers for building edit DAGs."""
from __future__ import annotations

import math
import random
from typing import Optional

from helix.edit.dag import EditDAG
from helix.edit.simulate import SimulationContext, build_edit_dag
from helix.edit.post import dedupe_terminal_nodes
from helix.genome.digital import DigitalGenome as CoreDigitalGenome

from . import rules  # noqa: F401
from .model import CasSystem, DigitalGenome as LegacyDigitalGenome, GuideRNA
from .physics import CRISPRPhysics


def build_crispr_edit_dag(
    genome: LegacyDigitalGenome,
    cas: CasSystem,
    guide: GuideRNA,
    *,
    rng_seed: int = 0,
    max_depth: int = 1,
    min_prob: float = 1e-4,
    max_sites: Optional[int] = 5,
) -> EditDAG:
    """
    Construct a CRISPR edit DAG using registered edit rules.
    """

    core_genome = CoreDigitalGenome(sequences=dict(genome.sequences))
    min_prob = max(min_prob, 1e-12)
    context = SimulationContext(
        rng=random.Random(rng_seed),
        max_depth=max_depth,
        min_log_prob=math.log(min_prob),
        rules=("crispr.clean_cut", "crispr.indel_branch", "crispr.no_edit"),
        extra={
            "legacy_genome": genome,
            "core_genome": core_genome,
            "cas": cas,
            "guide": guide,
            "max_sites": max_sites,
            "no_edit_prob": 0.1,
            "indel_window": 3,
            "physics": CRISPRPhysics.from_system(cas, guide),
        },
    )
    dag = build_edit_dag(core_genome.view(), context)
    return dedupe_terminal_nodes(dag)
