"""De Bruijn graph utilities."""
from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Set


@dataclass
class DBGNode:
    kmer: str
    out_edges: Set[str] = field(default_factory=set)
    in_edges: Set[str] = field(default_factory=set)
    coverage: int = 0


@dataclass
class DeBruijnGraph:
    k: int
    nodes: Dict[str, DBGNode] = field(default_factory=dict)


def build_dbg(reads: Iterable[str], k: int) -> DeBruijnGraph:
    graph = DeBruijnGraph(k=k)
    for read in reads:
        seq = read.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            node_prefix = graph.nodes.setdefault(prefix, DBGNode(prefix))
            node_suffix = graph.nodes.setdefault(suffix, DBGNode(suffix))
            node_prefix.out_edges.add(suffix)
            node_suffix.in_edges.add(prefix)
            node_suffix.coverage += 1
            node_prefix.coverage += 1
    return graph


def clean_dbg(graph: DeBruijnGraph, tips: bool = True, bubbles: bool = True, tip_length: int = 2) -> DeBruijnGraph:
    if tips:
        remove_tips(graph, tip_length)
    if bubbles:
        remove_simple_bubbles(graph)
    return graph


def remove_tips(graph: DeBruijnGraph, tip_length: int) -> None:
    to_remove = set()
    for node_id, node in graph.nodes.items():
        if len(node.in_edges) == 0 and len(node.out_edges) == 1:
            if tip_walk(graph, node_id, tip_length):
                to_remove.add(node_id)
    for node_id in to_remove:
        delete_node(graph, node_id)


def tip_walk(graph: DeBruijnGraph, start: str, tip_length: int) -> bool:
    current = start
    for _ in range(tip_length):
        node = graph.nodes.get(current)
        if not node or len(node.out_edges) != 1:
            return False
        nxt = next(iter(node.out_edges))
        if len(graph.nodes[nxt].in_edges) > 1:
            return False
        current = nxt
    return True


def remove_simple_bubbles(graph: DeBruijnGraph) -> None:
    for node_id, node in list(graph.nodes.items()):
        if len(node.out_edges) == 2:
            paths = list(node.out_edges)
            targets = []
            for path in paths:
                walk = path
                visited = set()
                while True:
                    if walk in visited or walk not in graph.nodes:
                        break
                    visited.add(walk)
                    successors = graph.nodes[walk].out_edges
                    if len(successors) != 1:
                        break
                    walk = next(iter(successors))
                    if walk in graph.nodes[node_id].in_edges:
                        targets.append(walk)
                        break
            if len(targets) == 2 and targets[0] == targets[1]:
                remove_path(graph, paths[0], targets[0])


def remove_path(graph: DeBruijnGraph, start: str, end: str) -> None:
    current = start
    while current != end and current in graph.nodes:
        nxt_nodes = list(graph.nodes[current].out_edges)
        delete_node(graph, current)
        if not nxt_nodes:
            break
        current = nxt_nodes[0]


def delete_node(graph: DeBruijnGraph, node_id: str) -> None:
    node = graph.nodes.pop(node_id, None)
    if not node:
        return
    for pred in node.in_edges:
        if pred in graph.nodes:
            graph.nodes[pred].out_edges.discard(node_id)
    for succ in node.out_edges:
        if succ in graph.nodes:
            graph.nodes[succ].in_edges.discard(node_id)


def serialize_graph(graph: DeBruijnGraph) -> Dict[str, object]:
    return {
        "k": graph.k,
        "nodes": {
            node_id: {
                "out": sorted(node.out_edges),
                "in": sorted(node.in_edges),
                "coverage": node.coverage,
            }
            for node_id, node in graph.nodes.items()
        },
    }


def deserialize_graph(payload: Dict[str, object]) -> DeBruijnGraph:
    graph = DeBruijnGraph(k=payload["k"])
    for node_id, data in payload["nodes"].items():
        node = DBGNode(kmer=node_id, out_edges=set(data["out"]), in_edges=set(data["in"]), coverage=data["coverage"])
        graph.nodes[node_id] = node
    return graph


def export_graphml(graph: DeBruijnGraph) -> str:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<graphml>",
        '  <graph edgedefault="directed">',
    ]
    for node_id, node in graph.nodes.items():
        lines.append(f'    <node id="{node_id}" />')
        for succ in node.out_edges:
            lines.append(f'    <edge source="{node_id}" target="{succ}" />')
    lines.extend(["  </graph>", "</graphml>"])
    return "\n".join(lines)
