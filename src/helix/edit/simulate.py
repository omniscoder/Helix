"""Generic edit DAG simulation runtime."""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Sequence

from helix.genome.digital import DigitalGenomeView

from .dag import EditDAG, EditEdge, EditNode
from .events import EditEvent
from .physics import get_rule


@dataclass
class SimulationContext:
    """
    Runtime configuration for building an edit DAG.

    Attributes
    ----------
    rng: random.Random
        RNG for stochastic rule behaviour.
    max_depth: int
        Maximum number of expansion layers (depth-first).
    min_log_prob: float
        Minimum log probability threshold for nodes (prunes very unlikely paths).
    rules: Sequence[str]
        Ordered list of rule names to apply at each depth.
    extra: dict
        Free-form metadata/rule configuration bag.
    """

    rng: random.Random
    max_depth: int = 1
    min_log_prob: float = -math.inf
    rules: Sequence[str] = field(default_factory=list)
    extra: Dict[str, object] = field(default_factory=dict)


def build_edit_dag(
    root_view: DigitalGenomeView,
    ctx: SimulationContext,
) -> EditDAG:
    """
    Construct an edit DAG by applying registered rules up to `max_depth`.
    """

    nodes: Dict[str, EditNode] = {}
    edges: List[EditEdge] = []

    root = EditNode(
        id="node_0",
        genome_view=root_view,
        log_prob=0.0,
        metadata={"time_step": 0, "stage": "root"},
    )
    nodes[root.id] = root
    frontier: List[EditNode] = [root]

    next_id = 1

    for depth in range(ctx.max_depth):
        if not frontier:
            break
        new_frontier: List[EditNode] = []
        for node in frontier:
            for rule_name in ctx.rules:
                rule = get_rule(rule_name)
                proposals = list(rule.propose(node, ctx))
                for event, logp_delta, metadata in proposals:
                    new_log_prob = node.log_prob + logp_delta
                    if new_log_prob < ctx.min_log_prob:
                        continue
                    new_view = node.genome_view.apply(event)
                    parent_time = node.metadata.get("time_step", 0)
                    parent_stage = node.metadata.get("stage", "root")
                    new_time = parent_time + 1
                    new_stage = metadata.get("stage", parent_stage)
                    node_id = f"node_{next_id}"
                    next_id += 1
                    new_node = EditNode(
                        id=node_id,
                        genome_view=new_view,
                        log_prob=new_log_prob,
                        metadata={
                            **node.metadata,
                            **metadata,
                            "time_step": new_time,
                            "stage": new_stage,
                        },
                    )
                    nodes[new_node.id] = new_node
                    edges.append(
                        EditEdge(
                            source=node.id,
                            target=new_node.id,
                            rule_name=rule_name,
                            event=event,
                            metadata=metadata,
                        )
                    )
                    new_frontier.append(new_node)
        frontier = new_frontier

    return EditDAG(nodes=nodes, edges=edges, root_id=root.id)
