import json
from pathlib import Path

from helix.crispr.guide import find_guides
from helix.crispr.pam import get_pam
from helix.crispr import score as crispr_score

try:  # pragma: no cover - support direct module invocation
    from tests.test_helix_cli import run_cli
except ModuleNotFoundError:  # pragma: no cover
    from test_helix_cli import run_cli


def test_find_guides_plus_and_minus():
    seq = "TTTTTACCCAGGAAACCCGGGTTTT"
    pam = get_pam("SpCas9-NGG")
    guides = find_guides(seq, pam, guide_len=4, strand="both")
    strands = {guide["strand"] for guide in guides}
    assert "+" in strands
    assert "-" in strands
    plus = next(guide for guide in guides if guide["strand"] == "+")
    assert plus["pam_site"]["end"] - plus["pam_site"]["start"] == len(pam["pattern"])
    minus = next(guide for guide in guides if guide["strand"] == "-")
    assert minus["start"] < minus["end"]
    assert minus["gc_content"] <= 1.0


def test_find_guides_window_limits_results():
    seq = "TTTTTACCCAGGAAACCCGGGTTTT"
    pam = get_pam("SpCas9-NGG")
    guides = find_guides(seq, pam, guide_len=4, strand="both", window=(0, 15))
    assert all(guide["end"] <= 15 for guide in guides)


def test_cli_crispr_find_guides(tmp_path: Path):
    fasta = tmp_path / "target.fna"
    fasta.write_text(">seq\nTTTTTACCCAGGAAACCCGGGTTTT\n", encoding="utf-8")
    json_path = tmp_path / "guides.json"
    run_cli(
        "crispr",
        "find-guides",
        "--fasta",
        str(fasta),
        "--pam",
        "SpCas9-NGG",
        "--guide-len",
        "4",
        "--json",
        str(json_path),
    )
    payload = json.loads(json_path.read_text())
    assert payload["schema"]["kind"] == "crispr.guides"
    assert payload["guides"]
    assert payload["guides"][0]["sequence"] is None


def test_enumerate_off_targets_detects_plus_and_minus():
    genome = "TTTTTACCCAGGAAACCCGGGTTTT"
    pam = get_pam("SpCas9-NGG")
    guide = {"id": "g1", "start": 5, "end": 9, "strand": "+", "sequence": "ACCC"}
    hits = crispr_score.enumerate_off_targets(genome, guide, pam, max_mm=1)
    assert any(hit["strand"] == "+" for hit in hits)
    assert all(hit["distance"] <= 1 for hit in hits)


def test_cli_crispr_offtargets_and_score(tmp_path: Path):
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">seq\nTTTTTACCCAGGAAACCCGGGTTTT\n", encoding="utf-8")
    guides_json = tmp_path / "guides.json"
    run_cli(
        "crispr",
        "find-guides",
        "--fasta",
        str(fasta),
        "--pam",
        "SpCas9-NGG",
        "--guide-len",
        "4",
        "--emit-sequences",
        "--json",
        str(guides_json),
    )

    hits_json = tmp_path / "hits.json"
    run_cli(
        "crispr",
        "offtargets",
        "--fasta",
        str(fasta),
        "--guides",
        str(guides_json),
        "--pam",
        "SpCas9-NGG",
        "--max-mm",
        "1",
        "--json",
        str(hits_json),
    )
    payload = json.loads(hits_json.read_text())
    assert payload["schema"]["kind"] == "crispr.offtargets"
    assert payload["guides"]
    assert payload["guides"][0]["hits"]

    weights = tmp_path / "weights.json"
    weights.write_text(
        json.dumps(
            {
                "on_target": {"gc_optimum": 0.5, "gc_weight": 0.5},
                "off_target": {"default_mismatch": 0.7, "pam_penalty": 0.6},
            }
        ),
        encoding="utf-8",
    )
    scored_json = tmp_path / "scored.json"
    run_cli(
        "crispr",
        "score",
        "--guides",
        str(guides_json),
        "--hits",
        str(hits_json),
        "--weights",
        str(weights),
        "--json",
        str(scored_json),
    )
    scored = json.loads(scored_json.read_text())
    guide_entry = scored["guides"][0]
    assert "on_target_score" in guide_entry
    assert "score" in guide_entry["hits"][0]


def test_cli_crispr_simulate_and_viz(tmp_path: Path):
    import pytest

    pytest.importorskip("matplotlib")
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">seq\nTTTTTACCCAGGAAACCCGGGTTTT\n", encoding="utf-8")
    guides_json = tmp_path / "guides.json"
    run_cli(
        "crispr",
        "find-guides",
        "--fasta",
        str(fasta),
        "--pam",
        "SpCas9-NGG",
        "--guide-len",
        "4",
        "--json",
        str(guides_json),
    )

    sim_json = tmp_path / "sim.json"
    guides_payload = json.loads(guides_json.read_text())
    guide_id = guides_payload["guides"][0]["id"]
    spec_path = tmp_path / "sim_spec.json"
    run_cli(
        "crispr",
        "simulate",
        "--fasta",
        str(fasta),
        "--guides",
        str(guides_json),
        "--guide-id",
        guide_id,
        "--draws",
        "200",
        "--seed",
        "123",
        "--json",
        str(sim_json),
        "--viz-spec",
        str(spec_path),
    )
    sim_payload = json.loads(sim_json.read_text())
    assert sim_payload["schema"]["kind"] == "crispr.sim"
    assert sim_payload["outcomes"]
    spec_payload = json.loads(spec_path.read_text())
    assert spec_payload["sequence"]
    assert spec_payload["edit_events"]

    output = tmp_path / "track.png"
    run_cli(
        "viz",
        "crispr-track",
        "--input",
        str(sim_json),
        "--save",
        str(output),
    )
    assert output.exists()
    assert output.with_suffix(".viz.json").exists()
