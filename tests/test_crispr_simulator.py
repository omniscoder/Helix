import json
from pathlib import Path

from helix import bioinformatics
from helix.crispr.model import (
    CasSystem,
    CasSystemType,
    DigitalGenome,
    GuideRNA,
    PAMRule,
)
from helix.crispr.simulator import find_candidate_sites, rank_off_targets, simulate_cuts
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.simulator import locate_prime_target_site, simulate_prime_edit
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


def test_prime_locate_and_simulate_outcomes():
    genome, guide_seq = _demo_genome()
    peg = PegRNA(spacer=guide_seq[:15], pbs="GAAAC", rtt="TTTTAA")
    editor = PrimeEditor(
        name="pe-demo",
        cas=_demo_cas(),
        nick_to_edit_offset=1,
        efficiency_scale=0.8,
        indel_bias=0.2,
    )
    site = locate_prime_target_site(genome, peg)
    assert site is not None
    outcomes = simulate_prime_edit(genome, editor, peg, max_outcomes=3)
    assert outcomes
    assert 0 < len(outcomes) <= 3
    assert all(outcome.logit_score >= 0 for outcome in outcomes)


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
    run_cli(
        "prime",
        "simulate",
        "--genome",
        str(fasta),
        "--peg-config",
        str(peg_config),
        "--editor-config",
        str(editor_config),
        "--max-outcomes",
        "2",
        "--json",
        str(out_path),
    )
    payload = json.loads(out_path.read_text())
    assert payload["schema"]["kind"] == "prime.edit_sim"
    assert payload["outcomes"]


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
    assert payload["artifact"] == "helix.crispr.edit_dag.v1"
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
    assert payload["artifact"] == "helix.prime.edit_dag.v1"
    assert payload["nodes"]


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
