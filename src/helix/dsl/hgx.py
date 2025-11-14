"""YAML based LiveGraph description language (HGX)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import yaml

from ..core.graph import GraphIR, build_graph
from ..nodes import ABMNode, CouplerNode, FieldNode, GRNNode, MechNode, ObserverNode, ProteinNode, RuleNetNode


NODE_REGISTRY = {
    "ProteinNode": ProteinNode,
    "RuleNetNode": RuleNetNode,
    "GRNNode": GRNNode,
    "FieldNode": FieldNode,
    "ABMNode": ABMNode,
    "MechNode": MechNode,
    "CouplerNode": CouplerNode,
    "ObserverNode": ObserverNode,
}


@dataclass
class HGXModel:
    name: str
    nodes: List[Dict[str, Any]]
    edges: List[Dict[str, Any]]


def load_hgx(path: str | Path) -> HGXModel:
    data = yaml.safe_load(Path(path).read_text())
    graph = data.get("graph", data)
    return HGXModel(name=graph["name"], nodes=graph.get("nodes", []), edges=graph.get("edges", []))


def dump_hgx(model: HGXModel, path: str | Path) -> None:
    payload = {"graph": {"name": model.name, "nodes": model.nodes, "edges": model.edges}}
    Path(path).write_text(yaml.safe_dump(payload, sort_keys=False))


def build_graph_from_hgx(model: HGXModel) -> GraphIR:
    nodes = []
    for node_spec in model.nodes:
        node_type = node_spec["type"]
        cls = NODE_REGISTRY.get(node_type)
        if cls is None:
            raise KeyError(f"Unknown node type {node_type}")
        params = node_spec.get("params", {})
        nodes.append(cls(name=node_spec["name"], **params))
    edges: List[Tuple[str, str, str, str]] = []
    for edge in model.edges:
        src = edge["source"]
        dst = edge["target"]
        edges.append((src["node"], src["port"], dst["node"], dst["port"]))
    return build_graph(nodes, edges)
