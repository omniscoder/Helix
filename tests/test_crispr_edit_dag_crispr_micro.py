from __future__ import annotations

import math
from pathlib import Path

from helix.cli import _edit_dag_to_payload, _frame_to_payload, _frames_to_artifact_payload, _load_digital_genome, _resolve_cas_system, _build_guide
from helix.crispr.dag_api import build_crispr_edit_dag
from helix.edit.dag import EditDAGFrame, dag_from_payload
from helix.edit.post import dedupe_terminal_nodes
PROB_TOL = 1e-6


def _collect_frames(genome_fasta: Path, guide_seq: str, *, seed: int = 0) -> tuple[dict, list[dict]]:
    genome = _load_digital_genome(genome_fasta)
    cas = _resolve_cas_system("cas9", None)
    guide = _build_guide(guide_seq, name="micro", pam=None, metadata=None)

    frames: list[EditDAGFrame] = []

    def consumer(frame: EditDAGFrame) -> None:
        frames.append(frame)

    dag = build_crispr_edit_dag(
        genome,
        cas,
        guide,
        rng_seed=seed,
        max_depth=2,
        min_prob=1e-4,
        max_sites=50,
        use_gpu=False,
        frame_consumer=consumer,
    )
    dag = dedupe_terminal_nodes(dag)
    static_payload = _edit_dag_to_payload(
        dag,
        artifact="helix.crispr.edit_dag.v1.1",
        metadata={"test": "crispr_micro"},
    )
    frame_payloads = [_frame_to_payload(frame) for frame in frames]
    frames_artifact = _frames_to_artifact_payload(frame_payloads)
    return static_payload, frames_artifact


def _dag_node_index(payload: dict) -> dict[str, dict]:
    nodes = payload.get("nodes") or {}
    return {node_id: node for node_id, node in nodes.items()}


def _dag_edge_multiset(payload: dict) -> dict[tuple[str, str], dict]:
    edges = payload.get("edges") or []
    result: dict[tuple[str, str], dict] = {}
    for edge in edges:
        key = (edge["source"], edge["target"])
        result[key] = {
            "rule": edge.get("rule"),
        }
    return result


def test_crispr_micro_static_vs_frames_roundtrip() -> None:
    genome_fasta = Path("tests/data/crispr_micro.fna")
    guide_seq = "ACGTACGTACGTACGTACGT"

    static_payload, frames_payload = _collect_frames(genome_fasta, guide_seq, seed=0)

    assert static_payload["artifact"] == "helix.crispr.edit_dag.v1.1"

    # Rebuild an EditDAG from frames, apply the same dedupe pass, and convert
    # back to a static payload so we compare like-for-like.
    frames_dag = dag_from_payload(frames_payload)
    frames_dag = dedupe_terminal_nodes(frames_dag)
    frames_static = _edit_dag_to_payload(
        frames_dag,
        artifact="helix.crispr.edit_dag.v1.1",
        metadata={"test": "crispr_micro_frames"},
    )

    assert frames_static["root_id"] == static_payload["root_id"]

    static_nodes = _dag_node_index(static_payload)
    frame_nodes = _dag_node_index(frames_static)

    # Frame-derived DAG may carry extra internal nodes (e.g. non-deduped paths),
    # but every static node must be present.
    assert set(static_nodes.keys()).issubset(set(frame_nodes.keys()))

    for nid, node_s in static_nodes.items():
        node_f = frame_nodes[nid]
        lp_s = float(node_s.get("log_prob", 0.0))
        lp_f = float(node_f.get("log_prob", 0.0))
        assert math.isclose(lp_s, lp_f, rel_tol=0.0, abs_tol=PROB_TOL)
        seq_s = node_s.get("sequences") or {}
        seq_f = node_f.get("sequences") or {}
        assert seq_s == seq_f

    edges_s = _dag_edge_multiset(static_payload)
    edges_f = _dag_edge_multiset(frames_static)

    assert set(edges_s.keys()).issubset(set(edges_f.keys()))
    for key, meta_s in edges_s.items():
        meta_f = edges_f[key]
        assert meta_s["rule"] == meta_f["rule"]
