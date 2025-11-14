"""Graph IR used by the Helix runtime."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Set, Tuple


@dataclass(frozen=True)
class Port:
    """Typed port definition (unit + positivity)."""

    name: str
    direction: str  # "in" or "out"
    unit: str = "arb"
    positive: bool = True


@dataclass(frozen=True)
class Edge:
    """Connection between two ports."""

    source: Tuple[str, str]
    target: Tuple[str, str]


@dataclass
class GraphNode:
    """Base class for all runtime nodes."""

    name: str
    kind: str
    ports: Dict[str, Port]
    state: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def step(self, t: float, dt: float, inputs: Mapping[str, Any]) -> Dict[str, Any]:
        """
        Advance the node one time step.

        Default implementation simply echos the provided inputs allowing the
        scheduler to wire pass-through nodes or placeholders.
        """

        return dict(inputs)

    def signature(self) -> str:
        """Stable identifier used by caches."""

        return self.metadata.get("hash") or f"{self.kind}:{self.name}"


class GraphIR:
    """Graph container with helper views for SCC + scheduler wiring."""

    def __init__(self) -> None:
        self.nodes: Dict[str, GraphNode] = {}
        self.edges: Set[Edge] = set()
        self._adj: MutableMapping[str, Set[str]] = {}
        self._incoming: MutableMapping[str, Set[Edge]] = {}
        self._outgoing: MutableMapping[str, Set[Edge]] = {}

    def add_node(self, node: GraphNode) -> None:
        if node.name in self.nodes:
            raise ValueError(f"node {node.name} already exists")
        self.nodes[node.name] = node
        self._adj.setdefault(node.name, set())
        self._incoming.setdefault(node.name, set())
        self._outgoing.setdefault(node.name, set())

    def connect(self, source_node: str, source_port: str, target_node: str, target_port: str) -> None:
        src = self.nodes[source_node]
        dst = self.nodes[target_node]
        if source_port not in src.ports:
            raise KeyError(f"{source_node} lacks port {source_port}")
        if target_port not in dst.ports:
            raise KeyError(f"{target_node} lacks port {target_port}")
        if src.ports[source_port].direction != "out":
            raise ValueError(f"{source_node}.{source_port} must be an output")
        if dst.ports[target_port].direction != "in":
            raise ValueError(f"{target_node}.{target_port} must be an input")
        edge = Edge(source=(source_node, source_port), target=(target_node, target_port))
        self.edges.add(edge)
        self._adj.setdefault(source_node, set()).add(target_node)
        self._incoming.setdefault(target_node, set()).add(edge)
        self._outgoing.setdefault(source_node, set()).add(edge)

    def adjacency(self) -> Mapping[str, Set[str]]:
        return self._adj

    def incoming_edges(self, node_name: str) -> Iterable[Edge]:
        return self._incoming.get(node_name, set())

    def outgoing_edges(self, node_name: str) -> Iterable[Edge]:
        return self._outgoing.get(node_name, set())

    def iter_nodes(self) -> Iterator[GraphNode]:
        return iter(self.nodes.values())

    def node(self, name: str) -> GraphNode:
        return self.nodes[name]

    def upstream(self, node_name: str) -> Set[str]:
        return {edge.source[0] for edge in self._incoming.get(node_name, set())}

    def downstream(self, node_name: str) -> Set[str]:
        return {edge.target[0] for edge in self._outgoing.get(node_name, set())}

    def strongly_connected_components(self) -> List[Tuple[str, ...]]:
        from .scc import Graph as SCCGraph, tarjan_scc

        g = SCCGraph()
        for src, targets in self._adj.items():
            for dst in targets:
                g.add_edge(src, dst)
        return tarjan_scc(g)

    def condensation_dag(self) -> Tuple[List[Tuple[str, ...]], "GraphIR"]:
        from .scc import Graph as SCCGraph, condensation_dag

        g = SCCGraph()
        for src, targets in self._adj.items():
            for dst in targets:
                g.add_edge(src, dst)
        comps, dag = condensation_dag(g)
        condensation = GraphIR()
        for idx, comp in enumerate(comps):
            node = GraphNode(name=f"island_{idx}", kind="island", ports={})
            condensation.add_node(node)
        for src, targets in dag.adj.items():
            for dst in targets:
                condensation._adj.setdefault(src, set()).add(dst)
        return comps, condensation


def build_graph(nodes: Iterable[GraphNode], edges: Iterable[Tuple[str, str, str, str]]) -> GraphIR:
    """Convenience helper to build a graph from iterables."""

    graph = GraphIR()
    for node in nodes:
        graph.add_node(node)
    for src_node, src_port, dst_node, dst_port in edges:
        graph.connect(src_node, src_port, dst_node, dst_port)
    return graph
