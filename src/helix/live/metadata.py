"""Metadata helpers for realtime delta payloads."""

from __future__ import annotations

from typing import Dict

from helix.core.graph import GraphIR, GraphNode


def describe_graph_nodes(graph: GraphIR) -> Dict[str, Dict[str, object]]:
    """
    Build a serializable node metadata map for the realtime delta schema.

    Each entry captures the node kind and the available ports (direction/unit/positivity)
    so downstream visualizers can attach semantics without importing Helix types.
    """

    description: Dict[str, Dict[str, object]] = {}
    for node in graph.iter_nodes():
        description[node.name] = _node_entry(node)
    return description


def _node_entry(node: GraphNode) -> Dict[str, object]:
    port_entries = {}
    for name, port in node.ports.items():
        port_entries[name] = {
            "direction": port.direction,
            "unit": port.unit,
            "positive": bool(port.positive),
        }
    return {"kind": node.kind, "ports": port_entries}
