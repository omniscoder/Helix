"""Prime editing DAG helpers."""
from __future__ import annotations

import math
import random

from typing import Callable, Optional

from helix.edit.dag import EditDAG, EditDAGFrame
from helix.edit.simulate import SimulationContext, build_edit_dag
from helix.edit.post import dedupe_terminal_nodes
from helix.genome.digital import DigitalGenome as CoreDigitalGenome

from . import rules  # noqa: F401
from .model import PegRNA, PrimeEditor
from .physics import PrimePhysics
from helix.crispr.model import DigitalGenome as LegacyDigitalGenome


def build_prime_edit_dag(
    genome: LegacyDigitalGenome,
    editor: PrimeEditor,
    peg: PegRNA,
    *,
    rng_seed: int = 0,
    max_depth: int = 1,
    min_prob: float = 1e-4,
    frame_consumer: Optional[Callable[[EditDAGFrame], None]] = None,
) -> EditDAG:
    min_prob = max(min_prob, 1e-12)
    core_genome = CoreDigitalGenome(sequences=dict(genome.sequences))
    context = SimulationContext(
        rng=random.Random(rng_seed),
        max_depth=max_depth,
        min_log_prob=math.log(min_prob),
        rules=("prime.rtt_clean", "prime.flap_resolution", "prime.no_edit"),
        extra={
            "legacy_genome": genome,
            "core_genome": core_genome,
            "peg": peg,
            "editor": editor,
            "prime_physics": PrimePhysics.from_config(editor, peg),
        },
    )
    dag = build_edit_dag(core_genome.view(), context, frame_consumer=frame_consumer)
    return dedupe_terminal_nodes(dag)
