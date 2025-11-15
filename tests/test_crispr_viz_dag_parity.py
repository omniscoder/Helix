from __future__ import annotations

from helix.crispr.pam import build_crispr_pam_mask
from helix.crispr.simulate import resolve_crispr_priors, simulate_cut_repair
from helix.gui.modern.builders import build_crispr_viz_spec
from helix.gui.modern.dag_adapters import crispr_dag_to_viz_spec
from helix.crispr.dag_api import build_crispr_edit_dag
from helix.crispr.model import GuideRNA
from helix.cli import _edit_dag_to_payload

try:  # pragma: no cover - import path differs under pytest
    from tests.test_crispr_simulator import _demo_cas, _demo_genome
except ModuleNotFoundError:  # pragma: no cover - fallback for direct execution
    from test_crispr_simulator import _demo_cas, _demo_genome


def test_crispr_viz_spec_legacy_vs_dag_mode() -> None:
    """CRISPR viz spec invariants should match between legacy and DAG modes.

    This focuses on core geometric/semantic fields that drive the GUI:
      - sequence window
      - guide_range
      - pam_index
      - event type sequence
    The underlying probability model may differ between simulate_cut_repair
    and the edit DAG engine, so we intentionally do not assert on outcome
    labels or probabilities here.
    """

    # Shared CRISPR micro-universe taken from the simulator tests.
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    sequence = genome.sequences["chr_plus"]

    guide_start = 3  # see _demo_genome: "TTT{guide}AGGTTT"
    guide_end = guide_start + len(guide_seq)
    guide_dict = {
        "id": "demo-guide",
        "sequence": guide_seq,
        "start": guide_start,
        "end": guide_end,
        "strand": "+",
        "gc_content": 0.5,
    }

    # Legacy path: crispr.sim payload -> viz spec.
    priors = resolve_crispr_priors("default_indel")
    pam_profile = "SpCas9_NGG"
    pam_mask = build_crispr_pam_mask(sequence, guide_dict, pam_profile, 1.0)
    sim_payload = simulate_cut_repair(
        site_seq=sequence,
        guide=guide_dict,
        priors=priors,
        draws=200,
        seed=123,
        emit_sequence=True,
        pam_mask=pam_mask,
    )
    spec_legacy = build_crispr_viz_spec(sequence, guide_dict, sim_payload)
    assert spec_legacy is not None

    # DAG path: CRISPR edit DAG -> artifact payload -> viz spec via adapter.
    guide_obj = GuideRNA(sequence=guide_seq, pam="NGG", name="demo-guide", metadata={"strand": "+"})
    dag = build_crispr_edit_dag(
        genome,
        cas,
        guide_obj,
        rng_seed=123,
        max_depth=2,
        min_prob=1e-4,
        max_sites=50,
        use_gpu=False,
        frame_consumer=None,
    )
    dag_payload = _edit_dag_to_payload(
        dag,
        artifact="helix.crispr.edit_dag.v1.1",
        metadata={
            "source": "test.crispr_viz_dag_parity",
            "guide": guide_dict,
            "site_sequence": sequence,
        },
    )
    spec_dag = crispr_dag_to_viz_spec(dag_payload)
    assert spec_dag is not None

    # Core invariants: sequence, guide range, PAM index, and event types
    # should align between the two visualization paths.
    assert spec_legacy.sequence == spec_dag.sequence
    assert spec_legacy.guide_range == spec_dag.guide_range
    assert spec_legacy.pam_index == spec_dag.pam_index

    legacy_types = [event.type for event in spec_legacy.events]
    dag_types = [event.type for event in spec_dag.events]
    assert legacy_types == dag_types

    # The human-facing edit_type label should at least share the same prefix
    # (CRISPR – <guide_id>); outcome labels may legitimately differ.
    def _prefix(label: str) -> str:
        return label.split("(", 1)[0].strip()

    assert _prefix(spec_legacy.edit_type) == _prefix(spec_dag.edit_type)
