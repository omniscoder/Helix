from __future__ import annotations

from pathlib import Path

from helix.edit.dag import EditDAG, EditDAGFrame, EditEdge, EditNode
from helix.genome.digital import DigitalGenome
from helix.cli import _edit_dag_to_payload, _frame_to_payload, _frames_to_artifact_payload


def _build_toy_dag() -> tuple[EditDAG, list[EditDAGFrame]]:
    genome = DigitalGenome({"chr": "ACGT"})
    root_view = genome.view()

    root = EditNode(
        id="n_root",
        genome_view=root_view,
        log_prob=0.0,
        metadata={"time_step": 0, "stage": "root"},
        parents=(),
        seq_hashes={},
        diffs=(),
    )
    child = EditNode(
        id="n_child",
        genome_view=root_view,
        log_prob=-0.5,
        metadata={"time_step": 1, "stage": "repaired"},
        parents=("n_root",),
        seq_hashes={},
        diffs=(),
    )
    # Minimal no-op edit event to satisfy the payload conversion helpers.
    from helix.edit.events import EditEvent

    event = EditEvent(chrom="chr", start=0, end=0, replacement="", metadata={})
    edge = EditEdge(
        source="n_root",
        target="n_child",
        rule_name="toy.rule",
        event=event,
        metadata={},
    )
    dag = EditDAG(nodes={"n_root": root, "n_child": child}, edges=[edge], root_id="n_root")

    frames = [
        EditDAGFrame(step=0, new_nodes={"n_root": root}, new_edges=[]),
        EditDAGFrame(step=1, new_nodes={"n_child": child}, new_edges=[edge]),
    ]
    return dag, frames


def test_frames_roundtrip_to_artifact_payload() -> None:
    dag, frames = _build_toy_dag()

    # Static payload from the DAG itself.
    static_payload = _edit_dag_to_payload(
        dag,
        artifact="helix.crispr.edit_dag.v1.1",
        metadata={"test": "toy"},
    )

    # Frame payloads as emitted by the CLI.
    frame_payloads = [_frame_to_payload(frame) for frame in frames]
    reconstructed = _frames_to_artifact_payload(frame_payloads)

    # Root id must match.
    assert reconstructed["root_id"] == static_payload["root_id"]

    nodes_static = static_payload["nodes"]
    nodes_recon = reconstructed["nodes"]
    assert set(nodes_static.keys()) == set(nodes_recon.keys())
    for node_id, node in nodes_static.items():
        recon = nodes_recon[node_id]
        assert node["log_prob"] == recon["log_prob"]
        assert node["metadata"] == recon["metadata"]
        assert node["sequences"] == recon["sequences"]

    edges_static = static_payload["edges"]
    edges_recon = reconstructed["edges"]
    assert len(edges_static) == len(edges_recon)
    for edge_s, edge_r in zip(edges_static, edges_recon):
        assert edge_s["source"] == edge_r["source"]
        assert edge_s["target"] == edge_r["target"]
        assert edge_s["rule"] == edge_r["rule"]
