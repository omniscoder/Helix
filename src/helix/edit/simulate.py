"""Generic edit DAG simulation runtime."""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterator, List, Optional, Sequence

from helix.genome.digital import DigitalGenomeView

from .dag import EditDAG, EditEdge, EditNode, EditDAGFrame, _compute_seq_hashes
from .events import EditEvent
from .hashing import hash_event, hash_node_id
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


def iter_edit_dag_frames(
    root_view: DigitalGenomeView,
    ctx: SimulationContext,
) -> Iterator[EditDAGFrame]:
    nodes: Dict[str, EditNode] = {}
    root_id = hash_node_id((), "", "root", 0)
    root = EditNode(
        id=root_id,
        genome_view=root_view,
        log_prob=0.0,
        metadata={"time_step": 0, "stage": "root"},
        parents=(),
        seq_hashes=_compute_seq_hashes(root_view),
        diffs=(),
    )
    nodes[root.id] = root
    frontier: List[EditNode] = [root]
    yield EditDAGFrame(step=0, new_nodes={root.id: root}, new_edges=[])
    frame_index = 1
    for depth in range(ctx.max_depth):
        if not frontier:
            break
        new_frontier: List[EditNode] = []
        new_nodes: Dict[str, EditNode] = {}
        new_edges: List[EditEdge] = []
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
                    parents = (node.id,)
                    ev_hash = hash_event(event.chrom, event.start, event.end, event.replacement, str(new_stage))
                    node_id = hash_node_id(parents, ev_hash, str(new_stage), new_time)
                    seq_hashes = _compute_seq_hashes(new_view)
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
                        parents=parents,
                        seq_hashes=seq_hashes,
                        diffs=(event,),
                    )
                    new_edge = EditEdge(
                        source=node.id,
                        target=new_node.id,
                        rule_name=rule_name,
                        event=event,
                        metadata=metadata,
                    )
                    new_nodes[new_node.id] = new_node
                    new_edges.append(new_edge)
                    new_frontier.append(new_node)
        frontier = new_frontier
        if new_nodes or new_edges:
            yield EditDAGFrame(step=frame_index, new_nodes=new_nodes, new_edges=new_edges)
            frame_index += 1


def build_edit_dag(
    root_view: DigitalGenomeView,
    ctx: SimulationContext,
    *,
    frame_consumer: Optional[Callable[[EditDAGFrame], None]] = None,
) -> EditDAG:
    """
    Construct an edit DAG by applying registered rules up to `max_depth`.
    Optionally stream intermediate frames via `frame_consumer`.
    """

    nodes: Dict[str, EditNode] = {}
    edges: List[EditEdge] = []
    root_id: Optional[str] = None

    for frame in iter_edit_dag_frames(root_view, ctx):
        if frame_consumer:
            frame_consumer(frame)
        nodes.update(frame.new_nodes)
        edges.extend(frame.new_edges)
        if root_id is None and frame.new_nodes:
            root_id = next(iter(frame.new_nodes))

    if root_id is None:
        raise RuntimeError("Failed to build edit DAG: no root node generated.")

    return EditDAG(nodes=nodes, edges=edges, root_id=root_id)
