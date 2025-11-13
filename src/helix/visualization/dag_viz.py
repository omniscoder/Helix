"""Edit DAG visualization helpers for Helix."""
from __future__ import annotations

import math
from typing import Mapping, Optional, Tuple

from matplotlib import colors as mcolors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import networkx as nx

from helix.edit.dag import EditDAG, EditNode
from ._dag_utils import compute_layout, stage_palette, format_stage_label


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
        stage_label = format_stage_label(stage)
        primary = stage_label
        if time_step is not None:
            primary += f" · t={time_step}"
        label = f"{primary}\n{node_id}"
        graph.add_node(
            node_id,
            label=label,
            display_label=label,
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
    node_size: int = 900,
    with_labels: bool = True,
    cmap_name: str = "viridis",
    stage_cmap: str = "tab10",
    layout: str = "spring",
    show_colorbar: bool = True,
    background_color: str = "#0b0e14",
    seed: int = 0,
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

    pos = compute_layout(graph, layout, seed=seed)
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
    node_sizes = [node_size * (0.6 + 0.4 * value) for value in normalized]
    node_norm_map = {node: value for node, value in zip(graph.nodes, normalized)}

    stage_values = [graph.nodes[node].get("stage", "unknown") for node in graph.nodes]
    stage_colors = stage_palette(stage_values, stage_cmap)
    border_colors = [stage_colors.get(stage, "#FFFFFF") for stage in stage_values]

    ax.set_facecolor(background_color)
    ax.figure.set_facecolor(background_color)

    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors=border_colors,
        linewidths=2.0,
        cmap=cmap,
    )
    edge_widths = [1.1 + 3.2 * node_norm_map.get(edge[1], 0.5) for edge in graph.edges]
    nx.draw_networkx_edges(
        graph,
        pos,
        ax=ax,
        arrows=True,
        arrowstyle="->",
        edge_color="#a7b0c2",
        width=edge_widths,
        alpha=0.85,
        connectionstyle="arc3,rad=0.05",
    )

    if with_labels:
        labels = {node: graph.nodes[node]["display_label"] for node in graph.nodes}
        nx.draw_networkx_labels(
            graph,
            pos,
            labels,
            font_size=8,
            font_color="#f0f4ff",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#00000055", edgecolor="none"),
            ax=ax,
        )

    edge_labels = {(edge[0], edge[1]): graph.edges[edge]["rule"] for edge in graph.edges}
    nx.draw_networkx_edge_labels(
        graph,
        pos,
        edge_labels=edge_labels,
        font_size=7,
        font_color="#d8dee9",
        alpha=0.9,
        bbox=dict(boxstyle="round,pad=0.15", facecolor="#1f2430cc", edgecolor="none"),
        ax=ax,
    )

    if show_colorbar:
        norm = mcolors.Normalize(vmin=min_prob, vmax=max_prob)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        colorbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
        colorbar.set_label("Branch probability", color="#e6e1cf")
        colorbar.ax.yaxis.set_tick_params(color="#c0c5d2")
        plt.setp(plt.getp(colorbar.ax.axes, "yticklabels"), color="#c0c5d2")

    if stage_colors:
        handles = [
            Patch(edgecolor=color, facecolor="none", linewidth=2.0, label=stage)
            for stage, color in stage_colors.items()
        ]
        ax.legend(handles=handles, title="Stage", loc="upper left", frameon=False, labelcolor="#e6e1cf")

    ax.annotate(
        "",
        xy=(0.85, 0.92),
        xytext=(0.55, 0.92),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops=dict(arrowstyle="->", color="#7aa2f7", linewidth=1.4, alpha=0.35),
    )
    ax.text(
        0.70,
        0.94,
        "Process flow",
        transform=ax.transAxes,
        color="#c0c5d2",
        fontsize=8,
        ha="center",
        va="center",
        alpha=0.7,
    )

    leaves = [node for node, out_degree in graph.out_degree() if out_degree == 0]
    if leaves:
        top_leaves = sorted(leaves, key=lambda node: graph.nodes[node]["prob"], reverse=True)[:5]
        lines = [
            f"{format_stage_label(graph.nodes[node].get('stage', 'Outcome'))}: {graph.nodes[node]['prob']:.1%}"
            for node in top_leaves
        ]
        overlay_text = "Top outcomes\n" + "\n".join(lines)
        ax.text(
            0.98,
            0.02,
            overlay_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#f4f6ff",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#05070dcc", edgecolor="#3a4154", linewidth=0.6),
        )

    ax.axis("off")
    ax.set_title("Helix Edit DAG", color="#e6e1cf")
    return ax


def save_edit_dag_png(
    dag: EditDAG,
    out_path: str,
    *,
    figsize: Tuple[int, int] = (10, 8),
    **plot_kwargs,
) -> None:
    """Render an EditDAG and save it as a PNG."""
    fig, ax = plt.subplots(figsize=figsize)
    plot_edit_dag(dag, ax=ax, figsize=figsize, **plot_kwargs)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
