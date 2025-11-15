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

    # derive per-edge conditional probabilities from log_prob deltas
    node_logs: Dict[str, float] = {node_id: entry["log_prob"] for node_id, entry in node_entries.items()}
    edges_by_source: Dict[str, list[Mapping[str, Any]]] = {}
    for edge in edges:
        source = str(edge.get("source"))
        edges_by_source.setdefault(source, []).append(edge)

    edge_entries: list[Dict[str, Any]] = []
    for source, edge_list in edges_by_source.items():
        parent_log = node_logs.get(source, 0.0)
        weights: list[float] = []
        for edge in edge_list:
            target = edge.get("target")
            child_log = node_logs.get(target, 0.0)
            delta = child_log - parent_log
            if math.isinf(delta):
                weights.append(0.0)
            else:
                weights.append(math.exp(delta))
        total = sum(weights)
        if not math.isfinite(total) or total <= 0:
            total = 1.0
        for edge, weight in zip(edge_list, weights):
            target = edge.get("target")
            edge_entries.append(
                {
                    "source": source,
                    "target": target,
                    "rule": edge.get("rule") or edge.get("rule_name") or "",
                    "prob": weight / total,
                }
            )

    # propagate global masses from root using conditional edge probabilities
    node_masses = _compute_prob_masses(node_entries, edge_entries, root_id)
    for node_id, mass in node_masses.items():
        node_entries[node_id]["prob_mass"] = mass
    check = {
        "schema": {"kind": "helix.veribiota.lean_check", "version": schema_version},
        "dag_name": dag_name,
        "root_id": root_id,
        "node_count": len(node_entries),
        "edge_count": len(edge_entries),
        "nodes": node_entries,
        "edges": edge_entries,
    }
    stats = _compute_probability_stats(node_masses, edge_entries)
    check.update(stats)
    check["meta"] = payload.get("meta") or {}
    return check


def _safe_prob_mass(log_prob: float, root_log: float) -> float:
    if math.isinf(log_prob):
        return 0.0
    return math.exp(log_prob - root_log)


def _compute_prob_masses(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: Sequence[Mapping[str, Any]],
    root_id: str,
) -> Dict[str, float]:
    """Compute global probability mass for each node given conditional edges."""

    masses: Dict[str, float] = {node_id: 0.0 for node_id in nodes}
    if root_id not in masses:
        return masses
    masses[root_id] = 1.0

    by_source: Dict[str, list[Mapping[str, Any]]] = {}
    for edge in edges:
        source = str(edge.get("source"))
        by_source.setdefault(source, []).append(edge)

    # rely on monotonic time_step metadata for a topological order
    ordered_nodes = sorted(
        nodes.items(),
        key=lambda item: int(item[1].get("metadata", {}).get("time_step", 0)),
    )
    for node_id, entry in ordered_nodes:
        mass = masses.get(node_id, 0.0)
        if mass <= 0.0:
            continue
        for edge in by_source.get(node_id, []):
            target = edge.get("target")
            prob = float(edge.get("prob", 0.0))
            if target is None or prob <= 0.0:
                continue
            masses[target] = masses.get(target, 0.0) + mass * prob
    return masses


def _compute_probability_stats(
    node_masses: Mapping[str, float],
    edges: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    adjacency: Dict[str, list[str]] = {}
    for edge in edges:
        source = edge.get("source")
        target = edge.get("target")
        if source is None or target is None:
            continue
        adjacency.setdefault(str(source), []).append(str(target))
    terminal_nodes = [node_id for node_id in node_masses if node_id not in adjacency]
    terminal_sum = sum(node_masses.get(node_id, 0.0) for node_id in terminal_nodes)
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

    node_masses = _extract_prob_masses(nodes, root_id)
    _validate_graph_structure(nodes, edges)
    recorded_sum = float(check_payload.get("terminal_prob_sum", 0.0))
    _validate_probabilities(node_masses, edges, recorded_sum, prob_tolerance)


def _extract_prob_masses(nodes: Mapping[str, Mapping[str, Any]], root_id: str) -> Dict[str, float]:
    root_log = float(nodes[root_id].get("log_prob", 0.0))
    masses: Dict[str, float] = {}
    for node_id, entry in nodes.items():
        if "prob_mass" in entry:
            masses[node_id] = float(entry["prob_mass"])
        else:
            masses[node_id] = _safe_prob_mass(float(entry.get("log_prob", 0.0)), root_log)
    return masses


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
    node_masses: Mapping[str, float],
    edges: Sequence[Mapping[str, Any]],
    recorded_sum: float,
    prob_tolerance: float,
) -> None:
    by_source: Dict[str, list[Mapping[str, Any]]] = {}
    terminals = set(node_masses.keys())
    for edge in edges:
        source = str(edge.get("source"))
        by_source.setdefault(source, []).append(edge)
        terminals.discard(source)

    for source, edge_list in by_source.items():
        probs = [float(edge.get("prob", 0.0)) for edge in edge_list]
        total = sum(probs)
        _ensure(
            math.isfinite(total) and abs(total - 1.0) <= prob_tolerance,
            f"Outgoing probabilities from '{source}' sum to {total}, expected 1.",
        )

    terminal_sum = sum(node_masses.get(node_id, 0.0) for node_id in terminals)
    _ensure(abs(terminal_sum - recorded_sum) <= prob_tolerance, "Terminal probability sum mismatch.")
    _ensure(
        abs(terminal_sum - 1.0) <= prob_tolerance,
        f"Terminal probability sum {terminal_sum} deviates from 1.0",
    )
