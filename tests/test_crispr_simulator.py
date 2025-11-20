import json
from pathlib import Path

import numpy as np
import pytest

from helix import bioinformatics
from helix.crispr.model import (
    CasSystem,
    CasSystemType,
    DigitalGenome,
    GuideRNA,
    PAMRule,
)
from helix.crispr.simulator import (
    EfficiencyTargetRequest,
    find_candidate_sites,
    predict_efficiency_for_targets,
    rank_off_targets,
    simulate_cuts,
)
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.simulator import locate_prime_target_site, simulate_prime_edit
from helix.crispr.pam import build_prime_pam_mask
from helix.crispr.physics import (
    _encode_sequence_to_uint8,
    create_crispr_physics,
    score_pairs_encoded,
    score_pairs_encoded_multi,
)
from helix.crispr.dag_api import build_crispr_edit_dag
from helix.prime.dag_api import build_prime_edit_dag

try:  # pragma: no cover
    from tests.test_helix_cli import run_cli
except ModuleNotFoundError:  # pragma: no cover
    from test_helix_cli import run_cli


def _demo_cas() -> CasSystem:
    return CasSystem(
        name="demo-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG", description="demo")],
        cut_offset=3,
    )


def _demo_genome() -> tuple[DigitalGenome, str]:
    guide = "ACCCAGGAAACCCGGGTTTT"
    plus = f"TTT{guide}AGGTTT"
    rc = bioinformatics.reverse_complement(guide)
    minus = f"TTTCCG{rc}TTT"
    genome = DigitalGenome({"chr_plus": plus, "chr_minus": minus})
    return genome, guide


def _degenerate_cas() -> CasSystem:
    return CasSystem(
        name="demo-degenerate",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NRG", description="degenerate-demo")],
        cut_offset=3,
    )


def test_find_candidate_sites_detects_both_strands():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    sites = find_candidate_sites(genome, cas, guide)
    assert sites
    strands = {site.strand for site in sites}
    assert 1 in strands and -1 in strands
    assert all(site.on_target_score is not None for site in sites)


def test_simulate_cuts_returns_events():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    events = simulate_cuts(genome, cas, guide, max_events=1)
    assert len(events) == 1
    event = events[0]
    assert event.cut_position >= 0
    assert event.score == event.site.on_target_score


def test_rank_off_targets_limits_candidates():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    ranked = rank_off_targets(genome, cas, guide, max_candidates=1)
    assert len(ranked) == 1


def test_pam_enforcement_rejects_invalid_sites():
    genome, guide_seq = _demo_genome()
    cas = CasSystem(
        name="pam-aaa",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="AAA", description="strict")],
        cut_offset=3,
    )
    guide = GuideRNA(sequence=guide_seq)
    sites = find_candidate_sites(genome, cas, guide)
    assert not sites


def test_iupac_pam_patterns_supported():
    genome, guide_seq = _demo_genome()
    cas = _degenerate_cas()
    guide = GuideRNA(sequence=guide_seq)
    sites = find_candidate_sites(genome, cas, guide)
    assert sites  # AGG satisfies NRG (N=any, R=AG)


def test_seed_weight_penalizes_pam_proximal_mismatch():
    guide_seq = "ACCCAGGAAACCCGGGTTTT"
    pam = "AGG"
    near = guide_seq[:-1] + ("A" if guide_seq[-1] != "A" else "C")
    far = ("A" if guide_seq[0] != "A" else "C") + guide_seq[1:]
    genome = DigitalGenome(
        {
            "near_chr": f"TTT{near}{pam}TTT",
            "far_chr": f"TTT{far}{pam}TTT",
        }
    )
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    sites = find_candidate_sites(genome, cas, guide)
    scores = {site.chrom: site.on_target_score for site in sites}
    assert {"near_chr", "far_chr"} <= scores.keys()
    assert scores["far_chr"] > scores["near_chr"]


def test_predict_efficiency_batch_matches_single():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    requests = [
        EfficiencyTargetRequest(
            target_id=chrom,
            reference_sequence=sequence,
            guide=GuideRNA(sequence=guide_seq, name=f"{chrom}_guide"),
        )
        for chrom, sequence in genome.sequences.items()
    ]
    predictions = predict_efficiency_for_targets(cas, requests)
    lookup = {pred.target_id: pred.predicted_score for pred in predictions}
    for chrom, sequence in genome.sequences.items():
        manual_genome = DigitalGenome({"target": sequence})
        manual_sites = find_candidate_sites(manual_genome, cas, GuideRNA(sequence=guide_seq), max_sites=1)
        expected = float(manual_sites[0].on_target_score or 0.0) if manual_sites else 0.0
        assert chrom in lookup
        assert lookup[chrom] == expected


def test_crispr_physics_encoded_scoring_matches_sites():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    physics = create_crispr_physics(cas, guide)
    sites = find_candidate_sites(genome, cas, guide, max_sites=1)
    assert sites
    encoded = _encode_sequence_to_uint8(sites[0].sequence)
    direct = physics.on_target_score_encoded(encoded)
    assert direct == sites[0].on_target_score


def test_crispr_physics_instances_isolated_caches():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    sites = find_candidate_sites(genome, cas, guide, max_sites=1)
    assert sites
    encoded = _encode_sequence_to_uint8(sites[0].sequence)
    physics_a = create_crispr_physics(cas, guide)
    physics_b = create_crispr_physics(cas, guide)
    score_a1 = physics_a.on_target_score_encoded(encoded)
    score_a2 = physics_a.on_target_score_encoded(encoded)
    score_b = physics_b.on_target_score_encoded(encoded)
    assert score_a1 == score_a2 == score_b


def test_score_pairs_encoded_matches_single_site():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    physics = create_crispr_physics(cas, guide)
    sites = find_candidate_sites(genome, cas, guide, max_sites=3)
    assert sites
    windows = np.stack([_encode_sequence_to_uint8(site.sequence) for site in sites])
    guides = np.expand_dims(_encode_sequence_to_uint8(guide_seq), axis=0)
    matrix = score_pairs_encoded(guides, windows, physics)
    assert matrix.shape == (1, len(windows))
    for idx, window in enumerate(windows):
        expected = physics.on_target_score_encoded(window)
        assert matrix[0, idx] == pytest.approx(expected, rel=0, abs=1e-9)


def test_score_pairs_encoded_multiple_guides():
    guides = ["ACCCAGGAAACCCGGGTTTT", "TGGGGAAACCCGGGTTTACC"]
    windows = ["ACCCAGGAAACCCGGGTTTT", "ACCCAGGAAACCCGGGTTTA", "TGGGGAAACCCGGGTTTACC"]
    cas = _demo_cas()
    physics_list = [create_crispr_physics(cas, GuideRNA(sequence=seq)) for seq in guides]
    guides_enc = np.stack([_encode_sequence_to_uint8(seq) for seq in guides])
    windows_enc = np.stack([_encode_sequence_to_uint8(seq) for seq in windows])
    pam = np.zeros((len(guides), len(windows)))
    matrix = score_pairs_encoded_multi(physics_list, guides_enc, windows_enc, pam_penalties=pam)
    assert matrix.shape == (len(guides), len(windows))
    for idx, phys in enumerate(physics_list):
        solo = score_pairs_encoded(
            guides_enc[idx : idx + 1],
            windows_enc,
            phys,
            pam_penalties=pam[idx : idx + 1],
        )
        assert np.allclose(matrix[idx : idx + 1], solo)


def test_native_backend_requires_extension():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    try:
        from helix_engine import native as native_engine
    except Exception:
        native_engine = None
    if native_engine and native_engine.is_available():
        physics = create_crispr_physics(cas, guide, backend="native-cpu")
        assert physics.on_target_score_encoded(_encode_sequence_to_uint8(guide_seq)) >= 0.0
    else:
        with pytest.raises(RuntimeError):
            create_crispr_physics(cas, guide, backend="native-cpu")


def test_prime_locate_and_simulate_outcomes():
    genome, guide_seq = _demo_genome()
    peg = PegRNA(spacer=guide_seq[:15], pbs="GAAAC", rtt="TTTTAA")
    site = locate_prime_target_site(genome, peg)
    assert site is not None
    site_seq = genome.sequences[site.chrom]
    pam_mask = build_prime_pam_mask(site_seq, peg, "SpCas9_NGG", 1.0)
    payload = simulate_prime_edit(site_seq, peg, draws=200, seed=123, pam_mask=pam_mask, emit_sequence=True)
    assert payload["schema"]["kind"] == "prime.edit_sim"
    assert payload["outcomes"]
    assert payload["draws"] == 200


def test_cli_crispr_genome_sim(tmp_path: Path):
    genome, guide_seq = _demo_genome()
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">chr\n" + genome.sequences["chr_plus"] + "\n", encoding="utf-8")
    out_path = tmp_path / "cuts.json"
    run_cli(
        "crispr",
        "genome-sim",
        "--genome",
        str(fasta),
        "--guide-sequence",
        guide_seq,
        "--json",
        str(out_path),
    )
    payload = json.loads(out_path.read_text())
    assert payload["schema"]["kind"] == "crispr.cut_events"
    assert payload["events"]
    assert payload["meta"]["crispr_scoring_version"]
    assert payload["meta"]["crispr_engine_backend"]


def test_cli_prime_simulate(tmp_path: Path):
    genome, guide_seq = _demo_genome()
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">chr\n" + genome.sequences["chr_plus"] + "\n", encoding="utf-8")

    peg_config = tmp_path / "peg.json"
    peg_config.write_text(
        json.dumps(
            {
                "name": "peg-demo",
                "spacer": guide_seq[:15],
                "pbs": "GAAAC",
                "rtt": "TTTTAA",
            }
        ),
        encoding="utf-8",
    )

    editor_config = tmp_path / "editor.json"
    editor_config.write_text(
        json.dumps(
            {
                "name": "pe-demo",
                "nick_to_edit_offset": 1,
                "efficiency_scale": 0.8,
                "indel_bias": 0.2,
                "mismatch_tolerance": 3,
                "cas": {
                    "name": "demo-cas9",
                    "system_type": "cas9",
                    "pam_rules": [{"pattern": "NGG"}],
                    "cut_offset": 3,
                },
            }
        ),
        encoding="utf-8",
    )

    out_path = tmp_path / "prime.json"
    viz_spec = tmp_path / "prime_spec.json"
    run_cli(
        "prime",
        "simulate",
        "--genome",
        str(fasta),
        "--peg-config",
        str(peg_config),
        "--editor-config",
        str(editor_config),
        "--draws",
        "200",
        "--json",
        str(out_path),
        "--viz-spec",
        str(viz_spec),
    )
    payload = json.loads(out_path.read_text())
    assert payload["schema"]["kind"] == "prime.edit_sim"
    assert payload["outcomes"]
    assert payload["draws"] == 200
    spec_payload = json.loads(viz_spec.read_text())
    assert spec_payload["edit_events"]
    assert payload["meta"]["prime_scoring_version"]
    assert payload["meta"]["prime_engine_backend"]


def test_cli_crispr_dag(tmp_path: Path):
    genome, guide_seq = _demo_genome()
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">chr\n" + genome.sequences["chr_plus"] + "\n", encoding="utf-8")
    out_path = tmp_path / "dag.json"
    run_cli(
        "crispr",
        "dag",
        "--genome",
        str(fasta),
        "--guide-sequence",
        guide_seq,
        "--json",
        str(out_path),
    )
    payload = json.loads(out_path.read_text())
    assert payload["artifact"] == "helix.crispr.edit_dag.v1.1"
    assert payload["nodes"]


def test_cli_prime_dag(tmp_path: Path):
    genome, guide_seq = _demo_genome()
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">chr\n" + genome.sequences["chr_plus"] + "\n", encoding="utf-8")
    peg_config = tmp_path / "peg.json"
    peg_config.write_text(
        json.dumps({"spacer": guide_seq[:15], "pbs": "GAAAC", "rtt": "TTTTAA"}),
        encoding="utf-8",
    )
    editor_config = tmp_path / "editor.json"
    editor_config.write_text(
        json.dumps(
            {
                "name": "pe-demo",
                "nick_to_edit_offset": 0,
                "efficiency_scale": 0.5,
                "cas": {
                    "name": "demo-cas9",
                    "system_type": "cas9",
                    "pam_rules": [{"pattern": "NGG"}],
                    "cut_offset": 3,
                },
            }
        ),
        encoding="utf-8",
    )
    out_path = tmp_path / "prime_dag.json"
    run_cli(
        "prime",
        "dag",
        "--genome",
        str(fasta),
        "--peg-config",
        str(peg_config),
        "--editor-config",
        str(editor_config),
        "--json",
        str(out_path),
    )
    payload = json.loads(out_path.read_text())
    assert payload["artifact"] == "helix.prime.edit_dag.v1.1"
    assert payload["nodes"]


def test_prime_dag_contains_flap_stage():
    genome, guide_seq = _demo_genome()
    peg = PegRNA(spacer=guide_seq[:15], pbs="GAAAC", rtt="TTTTAA")
    editor = PrimeEditor(
        name="pe-demo",
        cas=_demo_cas(),
        nick_to_edit_offset=1,
        efficiency_scale=0.8,
        indel_bias=0.2,
    )
    dag = build_prime_edit_dag(genome, editor, peg, max_depth=2)
    stages = {node.metadata.get("stage") for node in dag.nodes.values()}
    assert "repaired" in stages


def test_crispr_edit_dag_builder():
    genome, guide_seq = _demo_genome()
    cas = _demo_cas()
    guide = GuideRNA(sequence=guide_seq)
    dag = build_crispr_edit_dag(genome, cas, guide, max_depth=1, max_sites=1)
    assert dag.nodes
    assert dag.root_id in dag.nodes


def test_prime_edit_dag_builder():
    genome, guide_seq = _demo_genome()
    peg = PegRNA(spacer=guide_seq[:15], pbs="GAAAC", rtt="TTTTAA")
    editor = PrimeEditor(name="pe-demo", cas=_demo_cas(), efficiency_scale=0.5)
    dag = build_prime_edit_dag(genome, editor, peg, max_depth=1)
    assert dag.nodes
