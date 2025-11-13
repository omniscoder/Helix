"""Edit DAG visualization helpers for Helix."""
from __future__ import annotations

import math
from typing import Mapping, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx

from helix.edit.dag import EditDAG, EditEdge, EditNode


def _node_probability(node: EditNode) -> float:
    """Convert log_prob to a probability in [0, 1]."""
    return math.exp(node.log_prob)


def _build_nx_graph(dag: EditDAG) -> Tuple[nx.DiGraph, Mapping[str, float]]:
    """Convert an EditDAG into a networkx DiGraph."""
    graph = nx.DiGraph()
    probs = {}
    for node_id, node in dag.nodes.items():
        prob = _node_probability(node)
        probs[node_id] = prob
        stage = node.metadata.get("stage", "unknown")
        time_step = node.metadata.get("time_step")
        label = f"{node_id}\\nstage={stage}"
        if time_step is not None:
            label += f"\\nt={time_step}"
        graph.add_node(
            node_id,
            label=label,
            prob=prob,
            stage=stage,
            time_step=time_step,
        )
    for edge in dag.edges:
        rule = edge.rule_name
        graph.add_edge(edge.source, edge.target, rule=rule)
    return graph, probs


def plot_edit_dag(
    dag: EditDAG,
    *,
    figsize: Tuple[int, int] = (10, 8),
    node_size: int = 800,
    with_labels: bool = True,
    cmap_name: str = "viridis",
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Render an EditDAG using networkx + matplotlib.

    Node colors represent branch probability (log_prob → probability).
    Edge labels show rule names.
    """
    graph, probs = _build_nx_graph(dag)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    pos = nx.spring_layout(graph, seed=0)
    probability_values = list(probs.values())
    if probability_values:
        min_prob = min(probability_values)
        max_prob = max(probability_values)
    else:
        min_prob = max_prob = 1.0
    span = max(max_prob - min_prob, 1e-9)
    normalized = [(probs[node] - min_prob) / span for node in graph.nodes]
    cmap = plt.get_cmap(cmap_name)
    node_colors = [cmap(value) for value in normalized]

    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=node_size,
        node_color=node_colors,
        edgecolors="black",
    )
    nx.draw_networkx_edges(graph, pos, ax=ax, arrows=True, arrowstyle="->")

    if with_labels:
        labels = {node: graph.nodes[node]["label"] for node in graph.nodes}
        nx.draw_networkx_labels(graph, pos, labels, font_size=8, ax=ax)

    edge_labels = {(edge[0], edge[1]): graph.edges[edge]["rule"] for edge in graph.edges}
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_size=7, ax=ax)

    ax.axis("off")
    ax.set_title("Helix Edit DAG (node color = probability)")
    return ax


def save_edit_dag_png(
    dag: EditDAG,
    out_path: str,
    *,
    figsize: Tuple[int, int] = (10, 8),
) -> None:
    """Render an EditDAG and save it as a PNG."""
    fig, ax = plt.subplots(figsize=figsize)
    plot_edit_dag(dag, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
