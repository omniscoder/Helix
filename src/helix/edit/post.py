"""Post-processing helpers for edit DAGs."""
from __future__ import annotations

import hashlib
import math
from collections import defaultdict

from .dag import EditDAG


def _seq_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def _log_sum_exp(values):
    if not values:
        return float("-inf")
    m = max(values)
    return m + math.log(sum(math.exp(v - m) for v in values))


def dedupe_terminal_nodes(dag: EditDAG) -> EditDAG:
    """
    Merge terminal nodes with identical sequences.

    The log_prob of the canonical node is updated via log-sum-exp; duplicates are removed.
    """
    if not dag.nodes:
        return dag
    root_view = dag.nodes[dag.root_id].genome_view
    chroms = list(root_view.base.sequences.keys())
    groups = defaultdict(list)
    for node in dag.terminal_nodes():
        key = []
        for chrom in chroms:
            seq = node.genome_view.materialize_chrom(chrom)
            key.append(_seq_hash(seq))
        groups[tuple(key)].append(node.id)

    drop_ids = set()
    for ids in groups.values():
        if len(ids) < 2:
            continue
        kept = ids[0]
        values = [dag.nodes[node_id].log_prob for node_id in ids]
        dag.nodes[kept].log_prob = _log_sum_exp(values)
        drop_ids.update(ids[1:])

    if not drop_ids:
        return dag

    dag.nodes = {node_id: node for node_id, node in dag.nodes.items() if node_id not in drop_ids}
    dag.edges = [
        edge for edge in dag.edges if edge.source not in drop_ids and edge.target not in drop_ids
    ]
    return dag
