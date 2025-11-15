"""Lean check summary generation and validation helpers."""
from __future__ import annotations

import json
import math
from typing import Any, Dict, Iterable, Mapping, Sequence


class LeanCheckError(RuntimeError):
    """Raised when a lean-check payload fails validation."""


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        return float("-inf")
    top = max(values)
    if math.isinf(top):
        return top
    return top + math.log(sum(math.exp(v - top) for v in values))


def _ensure(condition: bool, message: str) -> None:
    if not condition:
        raise LeanCheckError(message)


def build_lean_check(
    payload: Mapping[str, Any],
    dag_name: str,
    *,
    schema_version: str = "1.0",
) -> Dict[str, Any]:
    """Construct a lean-check summary from an edit DAG payload."""

    nodes = payload.get("nodes") or {}
    edges = payload.get("edges") or []
    if not isinstance(nodes, Mapping) or not isinstance(edges, Sequence):
        raise ValueError("Payload missing 'nodes' or 'edges'.")
    root_id = payload.get("root_id") or (nodes and next(iter(nodes))) or ""
    node_entries: Dict[str, Dict[str, Any]] = {}
    for node_id, node in nodes.items():
        node_entries[node_id] = {
            "log_prob": float(node.get("log_prob", 0.0)),
            "metadata": dict(node.get("metadata") or {}),
            "parent_ids": list(node.get("parent_ids") or node.get("parents") or []),
            "seq_hashes": dict(node.get("seq_hashes") or {}),
        }
    edge_entries = []
    for edge in edges:
        edge_entries.append(
            {
                "source": edge.get("source"),
                "target": edge.get("target"),
                "rule": edge.get("rule") or edge.get("rule_name") or "",
            }
        )
    check = {
        "schema": {"kind": "helix.veribiota.lean_check", "version": schema_version},
        "dag_name": dag_name,
        "root_id": root_id,
        "node_count": len(node_entries),
        "edge_count": len(edge_entries),
        "nodes": node_entries,
        "edges": edge_entries,
    }
    # derived stats
    stats = _compute_probability_stats(node_entries, edge_entries, root_id)
    check.update(stats)
    check["meta"] = payload.get("meta") or {}
    return check


def _compute_probability_stats(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: Sequence[Mapping[str, Any]],
    root_id: str,
) -> Dict[str, Any]:
    adjacency: Dict[str, list[str]] = {}
    for edge in edges:
        source = edge.get("source")
        target = edge.get("target")
        if source is None or target is None:
            continue
        adjacency.setdefault(str(source), []).append(str(target))
    root_log = float(nodes.get(root_id, {}).get("log_prob", 0.0))
    terminal_nodes = [node_id for node_id in nodes if node_id not in adjacency]
    terminal_sum = sum(math.exp(float(nodes[node_id]["log_prob"]) - root_log) for node_id in terminal_nodes)
    return {
        "terminal_nodes": terminal_nodes,
        "terminal_prob_sum": terminal_sum,
    }


def validate_lean_check(
    check_payload: Mapping[str, Any],
    *,
    prob_tolerance: float = 5e-3,
) -> None:
    """Validate structural and probability invariants encoded in a lean-check payload."""

    nodes = check_payload.get("nodes")
    edges = check_payload.get("edges")
    root_id = check_payload.get("root_id")
    _ensure(isinstance(nodes, Mapping) and nodes, "Lean check missing nodes.")
    _ensure(isinstance(edges, Sequence), "Lean check missing edges.")
    node_count = check_payload.get("node_count")
    edge_count = check_payload.get("edge_count")
    _ensure(len(nodes) == int(node_count or 0), "Node count mismatch.")
    _ensure(len(edges) == int(edge_count or 0), "Edge count mismatch.")
    _ensure(root_id in nodes, "Root node missing from lean check.")

    _validate_graph_structure(nodes, edges)
    recorded_sum = float(check_payload.get("terminal_prob_sum", 0.0))
    _validate_probabilities(nodes, edges, root_id, recorded_sum, prob_tolerance)


def _validate_graph_structure(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: Sequence[Mapping[str, Any]],
) -> None:
    indegree = {node_id: 0 for node_id in nodes}
    adjacency: Dict[str, list[str]] = {}
    for edge in edges:
        source = edge.get("source")
        target = edge.get("target")
        _ensure(source in nodes, f"Edge source '{source}' missing from nodes.")
        _ensure(target in nodes, f"Edge target '{target}' missing from nodes.")
        adjacency.setdefault(source, []).append(target)
        indegree[target] += 1

    # Kahn's algorithm to detect cycles
    queue = [node_id for node_id, deg in indegree.items() if deg == 0]
    visited = 0
    while queue:
        node = queue.pop()
        visited += 1
        for target in adjacency.get(node, []):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
    _ensure(visited == len(nodes), "Lean check graph contains a cycle.")


def _validate_probabilities(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: Sequence[Mapping[str, Any]],
    root_id: str,
    recorded_sum: float,
    prob_tolerance: float,
) -> None:
    by_source: Dict[str, list[str]] = {}
    for edge in edges:
        by_source.setdefault(edge["source"], []).append(edge["target"])

    for source, targets in by_source.items():
        parent_log = float(nodes[source]["log_prob"])
        probs = [math.exp(float(nodes[target]["log_prob"]) - parent_log) for target in targets]
        total = sum(probs)
        _ensure(
            math.isfinite(total) and abs(total - 1.0) <= prob_tolerance,
            f"Outgoing probabilities from '{source}' sum to {total}, expected 1.",
        )

    root_log = float(nodes[root_id]["log_prob"])
    terminals = [node_id for node_id in nodes if node_id not in by_source]
    terminal_sum = sum(math.exp(float(nodes[node_id]["log_prob"]) - root_log) for node_id in terminals)
    _ensure(abs(terminal_sum - recorded_sum) <= prob_tolerance, "Terminal probability sum mismatch.")
    _ensure(
        abs(terminal_sum - 1.0) <= prob_tolerance,
        f"Terminal probability sum {terminal_sum} deviates from 1.0",
    )
