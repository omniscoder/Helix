"""Adapters that turn EditDAG payloads into ModernGL visualization specs.

These are the bridge between the \"edit DAG oracle\" (helix.*.edit_dag.v1.* /
helix.edit_dag.frame.v1) and the existing 2.5D viz builders used by the Qt GUI.

Nothing in the GUI calls these yet; they're scaffolding for a future DAG-first
path so Qt/Playground can share the same semantics.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from math import exp
from typing import Any, Dict, List, Mapping, Optional, Tuple

from helix.gui.modern.spec import EditVisualizationSpec
from helix.gui.modern.builders import (
    build_crispr_viz_spec,
    build_prime_viz_spec,
)
from helix.edit.dag import dag_from_payload
from helix import bioinformatics
from helix.crispr.model import DigitalGenome
from helix.prime.model import PegRNA, PrimeEditor


@dataclass(frozen=True)
class CrisprVizContext:
    """Minimal inputs needed to feed the existing CRISPR viz builder."""

    sequence: str
    guide: Mapping[str, Any]
    sim_payload: Mapping[str, Any]


def _extract_root_and_leaves(dag) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """Return (root_node, leaf_nodes) in a simple dict form."""
    root = dag.nodes[dag.root_id]
    leaves: List[Dict[str, Any]] = []
    outgoing = {edge.source for edge in dag.edges}
    for node_id, node in dag.nodes.items():
        if node_id not in outgoing:
            leaves.append(
                {
                    "id": node.id,
                    "log_prob": float(node.log_prob),
                    "metadata": dict(node.metadata),
                    "sequences": node.genome_view.materialize_all(),
                }
            )
    root_entry = {
        "id": root.id,
        "log_prob": float(root.log_prob),
        "metadata": dict(root.metadata),
        "sequences": root.genome_view.materialize_all(),
    }
    return root_entry, leaves


def _prob_from_log_prob(log_p: float, root_log_p: float) -> float:
    """Convert node log_prob to probability relative to root."""
    delta = log_p - root_log_p
    if math.isinf(delta):
        return 0.0
    return float(exp(delta))


def _build_fake_crispr_sim_payload(root: Mapping[str, Any], leaves: List[Mapping[str, Any]]) -> Dict[str, Any]:
    """Construct a crispr.sim-like payload from DAG leaf nodes."""
    root_log_p = float(root.get("log_prob", 0.0))

    outcomes: List[Dict[str, Any]] = []
    for leaf in leaves:
        meta = leaf.get("metadata") or {}
        label = meta.get("label") or leaf.get("id") or "outcome"
        log_p = float(leaf.get("log_prob", 0.0))
        prob = _prob_from_log_prob(log_p, root_log_p)
        diff = meta.get("diff") or meta.get("edit") or None

        outcomes.append(
            {
                "label": label,
                "count": 0,
                "probability": prob,
                "diff": diff,
            }
        )

    total = sum(o["probability"] for o in outcomes) or 1.0
    for o in outcomes:
        o["probability"] = float(o["probability"]) / total

    sequences = root.get("sequences") or {}
    sequence = ""
    if isinstance(sequences, Mapping) and sequences:
        # Prefer first chromosome sequence as the \"site\" window.
        sequence = next(iter(sequences.values()))

    payload: Dict[str, Any] = {
        "schema": {"kind": "crispr.sim.from_dag", "spec_version": "1.0"},
        "meta": {"source": "dag_adapter"},
        "site": {"length": len(sequence)},
        "guide": {},
        "priors": {},
        "draws": 1,
        "outcomes": outcomes,
    }
    return payload


def _extract_crispr_viz_context_from_dag_payload(dag_payload: Mapping[str, Any]) -> Optional[CrisprVizContext]:
    """Best-effort extraction of CRISPR viz inputs from an edit_dag payload."""
    dag = dag_from_payload(dag_payload)
    root_node, leaf_nodes = _extract_root_and_leaves(dag)

    sequences = root_node.get("sequences") or {}
    sequence = ""
    if isinstance(sequences, Mapping) and sequences:
        sequence = next(iter(sequences.values()))
    sequence = bioinformatics.normalize_sequence(sequence)
    if not sequence:
        return None

    meta = dag_payload.get("meta") or {}
    guide_meta = meta.get("guide") or {}
    guide: Dict[str, Any] = dict(guide_meta)
    guide.setdefault("start", guide_meta.get("start", 0))
    guide.setdefault("end", guide_meta.get("end", len(sequence)))
    guide.setdefault("strand", guide_meta.get("strand", "+"))
    guide.setdefault("id", guide_meta.get("id", "guide"))
    guide.setdefault("gc_content", guide_meta.get("gc_content", 0.5))

    sim_payload = _build_fake_crispr_sim_payload(root_node, leaf_nodes)
    sim_payload["guide"] = guide

    return CrisprVizContext(sequence=sequence, guide=guide, sim_payload=sim_payload)


def crispr_dag_to_viz_spec(dag_payload: Mapping[str, Any]) -> Optional[EditVisualizationSpec]:
    """Turn a helix.crispr.edit_dag.v1.* payload into an EditVisualizationSpec.

    This is a thin adapter around build_crispr_viz_spec. It does NOT change any
    behavior yet; nothing calls this from the GUI until you wire it in.
    """
    ctx = _extract_crispr_viz_context_from_dag_payload(dag_payload)
    if ctx is None:
        return None
    return build_crispr_viz_spec(
        sequence=ctx.sequence,
        guide=ctx.guide,
        sim_payload=ctx.sim_payload,
    )


def prime_dag_to_viz_spec(_: Mapping[str, Any]) -> Optional[EditVisualizationSpec]:
    """Placeholder for prime DAG → EditVisualizationSpec adapter.

    Implement once the prime edit_dag schema carries enough metadata about
    PegRNA/PrimeEditor and site windows to reconstruct a sim-like payload.
    """
    return None

