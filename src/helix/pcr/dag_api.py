"""PCR DAG construction helpers."""
from __future__ import annotations

import math
import random

from helix.edit.dag import EditDAG
from helix.edit.simulate import SimulationContext, build_edit_dag
from helix.edit.post import dedupe_terminal_nodes
from helix.genome.digital import DigitalGenome as CoreDigitalGenome
from helix.crispr.model import DigitalGenome as LegacyDigitalGenome

from .model import PCRConfig, PrimerPair

# Ensure rules are registered
from . import rules  # noqa: F401


def pcr_edit_dag(
    genome: LegacyDigitalGenome,
    primer_pair: PrimerPair,
    config: PCRConfig,
    *,
    rng_seed: int = 0,
    min_prob: float = 1e-6,
) -> EditDAG:
    """
    Build an EditDAG describing PCR amplification dynamics.
    """

    min_prob = max(min_prob, 1e-12)
    core_genome = CoreDigitalGenome(sequences=dict(genome.sequences))
    rules = ["pcr.primer_binding", "pcr.amplify_cycle"]
    if config.error_rate > 0:
        rules.append("pcr.error_branch")

    context = SimulationContext(
        rng=random.Random(rng_seed),
        max_depth=max(1, config.cycles + 1),
        min_log_prob=math.log(min_prob),
        rules=rules,
        extra={
            "core_genome": core_genome,
            "legacy_genome": genome,
            "primer_pair": primer_pair,
            "config": config,
        },
    )
    dag = build_edit_dag(core_genome.view(), context)
    root = dag.nodes.get(dag.root_id)
    if root:
        root.metadata.setdefault("stage", "root")
        root.metadata.setdefault("time_step", 0)
    return dedupe_terminal_nodes(dag)
