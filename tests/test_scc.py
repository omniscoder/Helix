from helix.core.graph import GraphIR, GraphNode, Port, build_graph
from helix.core import scc


def make_node(name: str) -> GraphNode:
    ports = {
        "out": Port(name="out", direction="out"),
        "inp": Port(name="inp", direction="in"),
    }
    return GraphNode(name=name, kind="stub", ports=ports)


def test_tarjan_scc_identifies_cycles():
    graph = build_graph(
        [make_node("a"), make_node("b"), make_node("c")],
        [("a", "out", "b", "inp"), ("b", "out", "a", "inp"), ("b", "out", "c", "inp")],
    )

    comps = graph.strongly_connected_components()
    assert any(set(comp) == {"a", "b"} for comp in comps)


def test_condensation_dag_shape():
    g = scc.Graph()
    g.add_edge("x", "y")
    g.add_edge("y", "z")
    comps, dag = scc.condensation_dag(g)
    assert len(comps) == 3
    assert dag.adj[0] or dag.adj[1] or dag.adj[2]  # ensures edges exist without assuming ordering
