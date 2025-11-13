import json
from pathlib import Path
from typing import Dict

from helix.crispr.dag_api import build_crispr_edit_dag
from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, GuideRNA, PAMRule
from helix.edit.dag import EditDAG, EditEdge, EditNode
from helix.edit.events import EditEvent
from helix.edit.post import dedupe_terminal_nodes

try:  # pragma: no cover
    from tests.test_helix_cli import run_cli
except ModuleNotFoundError:  # pragma: no cover
    from test_helix_cli import run_cli


def _demo_genome() -> tuple[DigitalGenome, str]:
    guide = "ACCCAGGAAACCCGGGTTTT"
    plus = f"TTT{guide}AGGTTT"
    genome = DigitalGenome({"chr": plus})
    return genome, guide


def _demo_cas() -> CasSystem:
    return CasSystem(
        name="demo-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
        max_mismatches=3,
    )


def _deterministic_crispr_dag() -> EditDAG:
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    return build_crispr_edit_dag(genome, cas, guide, rng_seed=42, max_depth=1)


def test_dag_deterministic_ids():
    dag1 = _deterministic_crispr_dag()
    dag2 = _deterministic_crispr_dag()
    assert set(dag1.nodes.keys()) == set(dag2.nodes.keys())
    for node_id in dag1.nodes:
        assert dag1.nodes[node_id].seq_hashes == dag2.nodes[node_id].seq_hashes


def test_terminal_dedupe_logsumexp():
    dag = _deterministic_crispr_dag()
    terminal = dag.terminal_nodes()[0]
    parent_edge = next(edge for edge in dag.edges if edge.target == terminal.id)
    duplicate_id = f"{terminal.id}_dup"
    duplicate_node = EditNode(
        id=duplicate_id,
        genome_view=terminal.genome_view,
        log_prob=terminal.log_prob - 1.0,
        metadata=dict(terminal.metadata),
        parents=(parent_edge.source,),
        seq_hashes=dict(terminal.seq_hashes),
        diffs=terminal.diffs,
    )
    dag.nodes[duplicate_id] = duplicate_node
    dag.edges.append(
        EditEdge(
            source=parent_edge.source,
            target=duplicate_id,
            rule_name=parent_edge.rule_name,
            event=parent_edge.event,
            metadata=parent_edge.metadata,
        )
    )
    combined = dedupe_terminal_nodes(dag)
    assert duplicate_id not in combined.nodes


def test_edit_dag_dataset_cli(tmp_path: Path):
    out = tmp_path / "dataset.jsonl"
    run_cli(
        "edit-dag",
        "generate-dataset",
        "--mechanism",
        "crispr",
        "--n",
        "2",
        "--out",
        str(out),
        "--seed",
        "7",
    )
    lines = [line for line in out.read_text().splitlines() if line.strip()]
    assert len(lines) == 2
    record: Dict[str, object] = json.loads(lines[0])
    assert "artifact" in record
    assert record["mechanism"] == "crispr"
