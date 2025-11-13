"""Edit DAG visualization helpers for Helix."""
from __future__ import annotations

import math
from typing import Dict, Mapping, Optional, Tuple

import numpy as np
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import networkx as nx

from helix.edit.dag import EditDAG, EditNode
from ._dag_utils import compute_layout, stage_palette, format_stage_label

COMPARE_COLORS = {
    "a": "#7dd3fc",
    "b": "#f4a261",
    "shared": "#f0f4ff",
}


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
        time_value = int(time_step) if isinstance(time_step, (int, float)) else 0
        stage_label = format_stage_label(stage)
        icon = {
            "root": "●",
            "binding": "✂",
            "cut": "✂",
            "cycled": "↺",
            "amplicon": "🧬",
            "error": "⚠",
            "no_product": "Ø",
        }.get(stage, "●")
        label = f"{stage_label} {icon}\n(t={time_value}, n={node_id})"
        graph.add_node(
            node_id,
            label=label,
            display_label=label,
            prob=prob,
            stage=stage,
            time_step=time_step,
            subset=time_value,
        )
    for edge in dag.edges:
        rule = edge.rule_name
        graph.add_edge(edge.source, edge.target, rule=rule)
    return graph, probs


def _normalize_probabilities(probs: Mapping[str, float]) -> Dict[str, float]:
    total = sum(probs.values()) or 1.0
    return {node: value / total for node, value in probs.items()}


def plot_edit_dag(
    dag: EditDAG,
    *,
    figsize: Tuple[int, int] = (10, 8),
    node_size: int = 900,
    with_labels: bool = True,
    cmap_name: str = "viridis",
    stage_cmap: str = "tab10",
    layout: str = "timeline",
    show_colorbar: bool = True,
    background_color: str = "#0b0e14",
    seed: int = 0,
    min_prob_filter: float = 0.0,
    max_time_filter: Optional[int] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Render an EditDAG using networkx + matplotlib.

    Node colors represent branch probability (log_prob → probability).
    Edge labels show rule names.
    """
    graph, probs = _build_nx_graph(dag)
    normalized_probs = _normalize_probabilities(probs)
    prob_values_full = list(probs.values())
    total_prob = sum(prob_values_full) or 1.0
    normalized_probs = {node: probs[node] / total_prob for node in graph.nodes}
    prob_threshold = max(0.0, min_prob_filter)

    def _visible(node: str) -> bool:
        if normalized_probs.get(node, 0.0) < prob_threshold:
            return False
        node_time = graph.nodes[node].get("subset", 0) or 0
        if max_time_filter is not None and node_time > max_time_filter:
            return False
        return True

    visible_nodes = [node for node in graph.nodes if _visible(node)]
    if not visible_nodes:
        raise ValueError("No nodes satisfy the filtering criteria; relax --min-prob/--max-time.")
    graph = graph.subgraph(visible_nodes).copy()
    probs = {node: probs[node] for node in graph.nodes}
    normalized_probs = {node: normalized_probs[node] for node in graph.nodes}

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
    time_positions: Mapping[int, list[float]] = {}
    for node, (x, y) in pos.items():
        time_step = graph.nodes[node].get("subset", 0)
        bucket = time_positions.setdefault(time_step, [])
        bucket.append(x)

    stage_values = [graph.nodes[node].get("stage", "unknown") for node in graph.nodes]
    stage_colors = stage_palette(stage_values, stage_cmap)
    border_colors = [stage_colors.get(stage, "#FFFFFF") for stage in stage_values]

    ax.set_facecolor(background_color)
    ax.figure.set_facecolor(background_color)

    # DNA watermark (subtle double helix)
    helix_x = np.linspace(-np.pi, np.pi, 200)
    for phase, alpha in [(0, 0.03), (np.pi, 0.02)]:
        helix_y = 0.1 * np.sin(helix_x * 2 + phase)
        ax.plot(
            helix_x,
            helix_y,
            color="#6f7dbf",
            alpha=alpha,
            linewidth=2,
            transform=ax.transAxes,
        )

    # Temporal bands
    if time_positions:
        for t_value, xs in sorted(time_positions.items()):
            center = sum(xs) / len(xs)
            ax.axvline(center, color="#ffffff10", linewidth=1.0, linestyle="--")
            ax.text(
                center,
                1.02,
                f"t={t_value}",
                transform=ax.transData,
                color="#cdd2e3",
                fontsize=8,
                ha="center",
                va="bottom",
                alpha=0.7,
            )

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
    prob_threshold = 0.25
    halo_nodes = [node for node, value in node_norm_map.items() if value <= prob_threshold]
    if halo_nodes:
        xs = [pos[node][0] for node in halo_nodes]
        ys = [pos[node][1] for node in halo_nodes]
        halo_sizes = [node_size * 1.8 for _ in halo_nodes]
        ax.scatter(
            xs,
            ys,
            s=halo_sizes,
            color="#f5a97f",
            alpha=0.12,
            linewidths=0,
        )
    def _connection_style(rule: str) -> str:
        if "clean" in (rule or "").lower():
            return "arc3,rad=0.0"
        if "indel" in (rule or "").lower():
            return "arc3,rad=0.25"
        if "error" in (rule or "").lower():
            return "arc3,rad=-0.2"
        return "arc3,rad=0.05"

    margin_args: dict[str, float] = {}
    try:
        nx.draw_networkx_edges(graph, pos, edgelist=[], ax=ax, min_source_margin=0, min_target_margin=0)
        margin_args = {"min_source_margin": 10, "min_target_margin": 18}
    except TypeError:
        margin_args = {}

    edge_labels = {}
    for edge in graph.edges:
        width = 1.1 + 3.2 * node_norm_map.get(edge[1], 0.5)
        rule = graph.edges[edge].get("rule")
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[edge],
            ax=ax,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=18,
            edge_color="#a7b0c2",
            width=width,
            alpha=0.9,
            connectionstyle=_connection_style(rule),
            **margin_args,
        )
        edge_labels[(edge[0], edge[1])] = f"{rule or 'process'} p={normalized_probs.get(edge[1], 0.0):.2f}"

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
        legend_loc = "lower left" if layout == "timeline" else "upper left"
        ax.legend(
            handles=handles,
            title="Edit stage",
            loc=legend_loc,
            frameon=False,
            labelcolor="#e6e1cf",
        )

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
        top_leaves = sorted(leaves, key=lambda node: normalized_probs.get(node, 0.0), reverse=True)[:5]
        labels = [format_stage_label(graph.nodes[node].get("stage", "Outcome")) for node in top_leaves]
        values = [normalized_probs.get(node, 0.0) for node in top_leaves]
        inset = inset_axes(ax, width="25%", height="30%", loc="lower right", borderpad=1.0)
        inset.barh(range(len(values)), values, color="#7aa2f7")
        inset.set_yticks(range(len(values)))
        inset.set_yticklabels(labels, fontsize=7, color="#e6e1cf")
        inset.set_xlim(0, max(values) * 1.1 if values else 1)
        inset.set_xticks([])
        inset.set_title("Outcome probability", fontsize=8, color="#e6e1cf")
        inset.set_facecolor("#0b0e14")
        for spine in inset.spines.values():
            spine.set_color("#2a2f3a")

    ax.axis("off")
    ax.set_title(
        "Helix CRISPR Repair Simulation\nStochastic lineage graph of possible edit outcomes",
        color="#e6e1cf",
        fontsize=13,
        loc="left",
    )
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


def plot_edit_dag_compare(
    dag_a: EditDAG,
    dag_b: EditDAG,
    *,
    label_a: str = "A",
    label_b: str = "B",
    layout: str = "timeline",
    min_prob_filter: float = 0.0,
    max_time_filter: Optional[int] = None,
    figsize: Tuple[int, int] = (10, 8),
    background_color: str = "#0b0e14",
    seed: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Overlay two DAGs and highlight unique/shared branches."""
    graph_a, probs_a = _build_nx_graph(dag_a)
    graph_b, probs_b = _build_nx_graph(dag_b)
    norm_a = _normalize_probabilities(probs_a)
    norm_b = _normalize_probabilities(probs_b)

    combined = nx.DiGraph()
    nodes_a = set(graph_a.nodes)
    nodes_b = set(graph_b.nodes)
    for node in nodes_a | nodes_b:
        attrs = {}
        if node in graph_a:
            attrs.update(graph_a.nodes[node])
        if node in graph_b:
            attrs.setdefault("stage", graph_b.nodes[node].get("stage"))
            attrs.setdefault("subset", graph_b.nodes[node].get("subset"))
        combined.add_node(node, **attrs)

    edges_union = set(graph_a.edges) | set(graph_b.edges)
    for edge in edges_union:
        in_a = edge in graph_a.edges
        in_b = edge in graph_b.edges
        rule = None
        if in_a:
            rule = graph_a.edges[edge].get("rule")
        if rule is None and in_b:
            rule = graph_b.edges[edge].get("rule")
        combined.add_edge(
            edge[0],
            edge[1],
            in_a=in_a,
            in_b=in_b,
            rule=rule,
        )

    def _visible(node: str) -> bool:
        prob = max(norm_a.get(node, 0.0), norm_b.get(node, 0.0))
        if prob < min_prob_filter:
            return False
        time_value = combined.nodes[node].get("subset", 0) or 0
        if max_time_filter is not None and time_value > max_time_filter:
            return False
        return True

    allowed_nodes = [node for node in combined.nodes if _visible(node)]
    if not allowed_nodes:
        allowed_nodes = list(combined.nodes)
    combined = combined.subgraph(allowed_nodes).copy()
    graph = combined

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    pos = compute_layout(graph, layout, seed=seed)
    ax.set_facecolor(background_color)
    ax.figure.set_facecolor(background_color)

    stages = [graph.nodes[node].get("stage", "unknown") for node in graph.nodes]
    stage_colors = stage_palette(stages, "tab10")

    node_sizes = []
    node_colors = []
    for node in graph.nodes:
        weight = max(norm_a.get(node, 0.0), norm_b.get(node, 0.0))
        base_size = 800 + 600 * weight
        node_sizes.append(base_size)
        in_a = node in nodes_a
        in_b = node in nodes_b
        if in_a and in_b:
            node_colors.append(COMPARE_COLORS["shared"])
        elif in_a:
            node_colors.append(COMPARE_COLORS["a"])
        else:
            node_colors.append(COMPARE_COLORS["b"])

    border_colors = [
        stage_colors.get(graph.nodes[node].get("stage", "unknown"), "#FFFFFF")
        for node in graph.nodes
    ]

    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors=border_colors,
        linewidths=2.0,
    )

    for edge in graph.edges:
        in_a = graph.edges[edge].get("in_a")
        in_b = graph.edges[edge].get("in_b")
        if in_a and in_b:
            color = "#d8dee9"
        elif in_a:
            color = COMPARE_COLORS["a"]
        else:
            color = COMPARE_COLORS["b"]
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[edge],
            ax=ax,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=16,
            edge_color=color,
            width=2.0,
            alpha=0.9,
        )

    labels = {node: graph.nodes[node].get("label", node) for node in graph.nodes}
    nx.draw_networkx_labels(
        graph,
        pos,
        labels,
        font_size=7,
        font_color="#f0f4ff",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#00000055", edgecolor="none"),
        ax=ax,
    )

    legend_handles = [
        Patch(facecolor=COMPARE_COLORS["a"], edgecolor="none", label=f"{label_a} only"),
        Patch(facecolor=COMPARE_COLORS["b"], edgecolor="none", label=f"{label_b} only"),
        Patch(facecolor=COMPARE_COLORS["shared"], edgecolor="none", label="Shared"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=8, labelcolor="#e6e1cf")

    leaves_a = sorted(
        (node for node in graph.nodes if graph.out_degree(node) == 0 and node in nodes_a),
        key=lambda node: norm_a.get(node, 0.0),
        reverse=True,
    )[:3]
    leaves_b = sorted(
        (node for node in graph.nodes if graph.out_degree(node) == 0 and node in nodes_b),
        key=lambda node: norm_b.get(node, 0.0),
        reverse=True,
    )[:3]
    if leaves_a or leaves_b:
        inset = inset_axes(ax, width="30%", height="35%", loc="lower right", borderpad=1.0)
        y = list(range(len(leaves_a) + len(leaves_b)))
        bars = []
        for node in leaves_a:
            bars.append((format_stage_label(graph.nodes[node].get("stage", "Outcome")), norm_a.get(node, 0.0), label_a))
        for node in leaves_b:
            bars.append((format_stage_label(graph.nodes[node].get("stage", "Outcome")), norm_b.get(node, 0.0), label_b))
        bar_labels = [f"{entry[0]} ({entry[2]})" for entry in bars]
        values = [entry[1] for entry in bars]
        colors = [
            COMPARE_COLORS["a"] if entry[2] == label_a else COMPARE_COLORS["b"] for entry in bars
        ]
        inset.barh(range(len(bars)), values, color=colors)
        inset.set_yticks(range(len(bars)))
        inset.set_yticklabels(bar_labels, fontsize=7, color="#e6e1cf")
        inset.set_xlim(0, max(values) * 1.1 if values else 1)
        inset.set_xticks([])
        inset.set_facecolor("#0b0e14")
        inset.set_title("Top outcomes", fontsize=8, color="#e6e1cf")
        for spine in inset.spines.values():
            spine.set_color("#2a2f3a")

    ax.axis("off")
    ax.set_title(
        f"Helix Edit DAG Comparison\n{label_a} vs {label_b}",
        color="#e6e1cf",
        fontsize=13,
        loc="left",
    )
    return ax


def save_edit_dag_compare_png(
    dag_a: EditDAG,
    dag_b: EditDAG,
    out_path: str,
    *,
    figsize: Tuple[int, int] = (10, 8),
    **plot_kwargs,
) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    plot_edit_dag_compare(dag_a, dag_b, ax=ax, figsize=figsize, **plot_kwargs)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
