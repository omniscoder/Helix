from __future__ import annotations

from helix.crispr.dag_api import build_crispr_edit_dag
from helix.crispr.model import GuideRNA

try:  # pragma: no cover - import path differs under pytest
    from tests.test_crispr_simulator import _demo_cas, _demo_genome
except ModuleNotFoundError:  # pragma: no cover - fallback for direct execution
    from test_crispr_simulator import _demo_cas, _demo_genome


def test_crispr_dag_nodes_carry_edit_class_and_mechanism() -> None:
    """Ensure CRISPR edit DAG nodes/edges expose edit_class and mechanism metadata.

    This metadata is used by DAG clients (Qt / Playground) to color outcomes
    consistently with the core geometry builders.
    """

    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)

    dag = build_crispr_edit_dag(
        genome,
        cas,
        guide,
        rng_seed=123,
        max_depth=2,
        min_prob=1e-4,
        max_sites=10,
        use_gpu=False,
        frame_consumer=None,
    )

    classes = set()
    for node in dag.nodes.values():
        stage = node.metadata.get("stage")
        edit_class = node.metadata.get("edit_class")
        if stage in {"repaired", "error", "no_edit"}:
            assert isinstance(edit_class, str)
            classes.add(edit_class)

    # We expect at least a deletion-like class and a non-cut / substitution-like class.
    assert "deletion" in classes
    assert any(cls in {"substitution", "no_cut"} for cls in classes)

    mechanisms = {edge.metadata.get("mechanism") for edge in dag.edges}
    # All edges should either carry a mechanism label or explicitly opt out.
    assert any(m for m in mechanisms)

