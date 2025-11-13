"""Animated Edit DAG rendering helpers."""
from __future__ import annotations

import io
import math
from typing import List, Tuple, Optional

import imageio.v2 as iio
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import networkx as nx

from helix.edit.dag import EditDAG
from ._dag_utils import compute_layout, stage_palette, format_stage_label


def animate_edit_dag(
    dag: EditDAG,
    out_gif: str,
    *,
    fps: int = 3,
    layout: str = "timeline",
    figsize: Tuple[int, int] = (8, 6),
    background_color: str = "#0b0e14",
    prob_cmap: str = "viridis",
    stage_cmap: str = "tab10",
    highlight_color: str = "#f5a97f",
    dpi: int = 150,
    min_prob: float = 0.0,
    max_time_filter: Optional[int] = None,
) -> None:
    """
    Animate EditDAG growth over time and export to a GIF.

    Each frame reveals nodes/edges up to the current time_step (metadata attribute).
    """
    graph = nx.DiGraph()
    time_by_node: dict[str, int] = {}
    prob_by_node: dict[str, float] = {}
    stage_by_node: dict[str, str] = {}
    display_labels: dict[str, str] = {}
    max_time_seen = 0

    for node_id, node in dag.nodes.items():
        time_step = int(node.metadata.get("time_step", 0))
        stage = node.metadata.get("stage", "unknown")
        max_time_seen = max(max_time_seen, time_step)
        probability = math.exp(node.log_prob)
        time_by_node[node_id] = time_step
        prob_by_node[node_id] = probability
        stage_by_node[node_id] = stage
        stage_label = format_stage_label(stage)
        label = stage_label
        if time_step is not None:
            label += f" · t={time_step}"
        label += f"\n{node_id}"
        display_labels[node_id] = label
        graph.add_node(node_id)

    for edge in dag.edges:
        graph.add_edge(edge.source, edge.target)

    total_prob = sum(prob_by_node.values()) or 1.0
    normalized_probs = {node: prob_by_node[node] / total_prob for node in graph.nodes}
    prob_threshold = max(0.0, min_prob)
    allowed_nodes = [
        node
        for node in graph.nodes
        if normalized_probs.get(node, 0.0) >= prob_threshold
        and (max_time_filter is None or time_by_node.get(node, 0) <= max_time_filter)
    ]
    if not allowed_nodes:
        allowed_nodes = list(graph.nodes)
    graph = graph.subgraph(allowed_nodes).copy()
    time_by_node = {node: time_by_node[node] for node in graph.nodes}
    stage_by_node = {node: stage_by_node[node] for node in graph.nodes}
    display_labels = {node: display_labels[node] for node in graph.nodes}
    prob_by_node = {node: prob_by_node[node] for node in graph.nodes}
    normalized_probs = {node: normalized_probs[node] for node in graph.nodes}
    max_time_seen = max((time_by_node.get(node, 0) for node in graph.nodes), default=0)

    position = compute_layout(graph, layout, seed=0)
    frames: List = []

    prob_values = list(prob_by_node.values()) or [1.0]
    min_prob = min(prob_values)
    max_prob = max(prob_values)
    if math.isclose(max_prob, min_prob):
        max_prob = min_prob + 1e-9
    prob_norm = mcolors.Normalize(vmin=min_prob, vmax=max_prob)
    cmap = plt.get_cmap(prob_cmap)
    node_color_map = {node: cmap(prob_norm(prob_by_node[node])) for node in graph.nodes}
    node_size_map = {node: 600 + 900 * prob_norm(prob_by_node[node]) for node in graph.nodes}
    node_norm_map = {node: prob_norm(prob_by_node[node]) for node in graph.nodes}
    border_colors = stage_palette(stage_by_node.values(), stage_cmap)
    edge_width_map = {
        (source, target): 1.1 + 3.0 * node_norm_map.get(target, 0.5) for (source, target) in graph.edges
    }
    time_positions: dict[int, list[float]] = {}
    for node, (x, _) in position.items():
        bucket = time_positions.setdefault(time_by_node.get(node, 0), [])
        bucket.append(x)

    leaves = [node for node, out_degree in graph.out_degree() if out_degree == 0]
    if leaves:
        top_leaves = sorted(leaves, key=lambda node: normalized_probs[node], reverse=True)[:5]
        top_outcomes_lines = [
            f"{format_stage_label(stage_by_node[node])}: {normalized_probs[node]:.1%}" for node in top_leaves
        ]
        top_outcomes_text = "Top outcomes\n" + "\n".join(top_outcomes_lines)
    else:
        top_outcomes_text = None

    def _draw_nodes(ax, nodes, alpha: float = 1.0):
        if not nodes:
            return
        nx.draw_networkx_nodes(
            graph,
            position,
            nodelist=nodes,
            node_size=[node_size_map[node] for node in nodes],
            node_color=[node_color_map[node] for node in nodes],
            edgecolors=[border_colors.get(stage_by_node[node], "#FFFFFF") for node in nodes],
            linewidths=2.0,
            alpha=alpha,
            ax=ax,
        )

    for current_time in range(0, max_time_seen + 1):
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_facecolor(background_color)
        fig.patch.set_facecolor(background_color)

        if time_positions:
            for t_value, xs in sorted(time_positions.items()):
                center = sum(xs) / len(xs)
                ax.axvline(center, color="#ffffff10", linewidth=0.8, linestyle="--")
                ax.annotate(
                    f"t={t_value}",
                    xy=(center, 1.02),
                    xycoords=("data", "axes fraction"),
                    color="#cdd2e3",
                    fontsize=7,
                    ha="center",
                    va="bottom",
                    alpha=0.6,
                )

        visible_nodes = [node for node in graph.nodes if time_by_node.get(node, 0) <= current_time]
        fresh_nodes = [node for node in visible_nodes if time_by_node.get(node, 0) == current_time]
        prior_nodes = [node for node in visible_nodes if time_by_node.get(node, 0) < current_time]
        visible_edges = [
            (source, target)
            for (source, target) in graph.edges
            if time_by_node.get(source, 0) <= current_time and time_by_node.get(target, 0) <= current_time
        ]
        fresh_edges = [
            (source, target)
            for (source, target) in visible_edges
            if time_by_node.get(target, 0) == current_time
        ]

        def _connection_style(rule: str | None) -> str:
            if not rule:
                return "arc3,rad=0.05"
            rule_lower = rule.lower()
            if "clean" in rule_lower:
                return "arc3,rad=0.0"
            if "indel" in rule_lower:
                return "arc3,rad=0.25"
            if "error" in rule_lower:
                return "arc3,rad=-0.2"
            return "arc3,rad=0.05"

        margin_args = {}
        try:
            nx.draw_networkx_edges(graph, position, edgelist=[], ax=ax, min_source_margin=0, min_target_margin=0)
            margin_args = {"min_source_margin": 10, "min_target_margin": 18}
        except TypeError:
            margin_args = {}

        for edge in visible_edges:
            rule = graph.edges[edge].get("rule")
            nx.draw_networkx_edges(
                graph,
                position,
                edgelist=[edge],
                ax=ax,
                arrows=True,
                arrowstyle="-|>",
                arrowsize=18,
                edge_color="#6b7085",
                width=edge_width_map[edge],
                alpha=0.65,
                connectionstyle=_connection_style(rule),
                **margin_args,
            )
        for edge in fresh_edges:
            rule = graph.edges[edge].get("rule")
            nx.draw_networkx_edges(
                graph,
                position,
                edgelist=[edge],
                ax=ax,
                arrows=True,
                arrowstyle="-|>",
                arrowsize=20,
                edge_color=highlight_color,
                width=edge_width_map[edge] * 1.15,
                alpha=0.95,
                connectionstyle=_connection_style(rule),
                **margin_args,
            )

        _draw_nodes(ax, prior_nodes, alpha=0.85)
        _draw_nodes(ax, fresh_nodes, alpha=1.0)

        if fresh_nodes:
            xs = [position[node][0] for node in fresh_nodes]
            ys = [position[node][1] for node in fresh_nodes]
            sizes = [node_size_map[node] * 1.8 for node in fresh_nodes]
            ax.scatter(xs, ys, s=sizes, color=highlight_color, alpha=0.25, linewidths=0)

        nx.draw_networkx_labels(
            graph,
            position,
            labels={node: display_labels[node] for node in visible_nodes},
            font_size=8,
            font_color="#f0f4ff",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#00000055", edgecolor="none"),
            ax=ax,
        )

        if border_colors:
            handles = [
                Patch(edgecolor=color, facecolor="none", linewidth=2.0, label=stage)
                for stage, color in border_colors.items()
            ]
            legend_loc = "lower left" if layout == "timeline" else "upper left"
            legend = ax.legend(
                handles=handles,
                title="Stage",
                loc=legend_loc,
                frameon=False,
                fontsize=8,
            )
            if legend:
                legend.get_title().set_color("#e6e1cf")
                for text in legend.get_texts():
                    text.set_color("#e6e1cf")

        ax.text(
            0.02,
            0.98,
            f"time_step ≤ {current_time}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color="#e6e1cf",
            fontsize=12,
            fontweight="bold",
        )
        ax.annotate(
            "",
            xy=(0.85, 0.9),
            xytext=(0.55, 0.9),
            xycoords="axes fraction",
            textcoords="axes fraction",
            arrowprops=dict(arrowstyle="->", color="#7aa2f7", linewidth=1.3, alpha=0.35),
        )
        ax.text(
            0.7,
            0.92,
            "Process flow",
            transform=ax.transAxes,
            color="#c0c5d2",
            fontsize=8,
            ha="center",
            va="center",
            alpha=0.7,
        )
        if top_outcomes_text:
            ax.text(
                0.98,
                0.02,
                top_outcomes_text,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="#f4f6ff",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#05070dcc", edgecolor="#3a4154", linewidth=0.6),
            )
        ax.axis("off")
        ax.set_title(
            "Helix CRISPR Repair Simulation\nStochastic lineage graph of possible edit outcomes",
            color="#e6e1cf",
            fontsize=12,
            loc="left",
            pad=10,
        )
        fig.tight_layout()

        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=dpi)
        plt.close(fig)
        buffer.seek(0)
        frames.append(iio.imread(buffer))

    iio.mimsave(out_gif, frames, duration=1.0 / float(max(fps, 1)))
