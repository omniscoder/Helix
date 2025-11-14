import pytest

from helix.core.graph import GraphIR, GraphNode, Port
from helix.core import invariants


def build_simple_graph(unit_b: str = "arb") -> GraphIR:
    ports_a = {"out": Port(name="out", direction="out", unit="rate")}
    ports_b = {"inp": Port(name="inp", direction="in", unit=unit_b)}
    node_a = GraphNode(name="a", kind="source", ports=ports_a)
    node_b = GraphNode(name="b", kind="sink", ports=ports_b)
    graph = GraphIR()
    graph.add_node(node_a)
    graph.add_node(node_b)
    graph.connect("a", "out", "b", "inp")
    return graph


def test_unit_validation_passes():
    graph = build_simple_graph(unit_b="rate")
    invariants.validate_edge_units(graph)


def test_unit_validation_fails_on_mismatch():
    graph = build_simple_graph(unit_b="concentration")
    with pytest.raises(invariants.UnitMismatch):
        invariants.validate_edge_units(graph)


def test_snapshot_positivity():
    graph = build_simple_graph(unit_b="rate")
    snapshot = {"a": {"out": 0.1}}
    invariants.validate_snapshot(graph, snapshot)
    snapshot["a"]["out"] = -0.1
    with pytest.raises(invariants.NegativeSignal):
        invariants.validate_snapshot(graph, snapshot)
