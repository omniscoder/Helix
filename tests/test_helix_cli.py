import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from helix import datasets
from helix.schema import SPEC_VERSION

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    path_entries = [str(SRC)]
    if existing:
        path_entries.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(path_entries)
    result = subprocess.run(
        [sys.executable, "-m", "helix.cli", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        env=env,
    )
    if check and result.returncode != 0:
        raise AssertionError(f"Command failed: {result.stderr}")
    return result


def test_cli_dna_smoke():
    result = run_cli("dna", "--sequence", "ACGTACGT", "--window", "4", "--step", "2", "--k", "2", "--max-diff", "0")
    assert "GC content" in result.stdout
    assert "Top" in result.stdout


def test_cli_spectrum_leaderboard():
    spectrum = "0,113,114,128,227,242,242,355,356,370,371,484"
    result = run_cli(
        "spectrum",
        "--peptide",
        "NQEL",
        "--spectrum",
        spectrum,
        "--leaderboard",
        "3",
    )
    assert "Score vs provided spectrum" in result.stdout
    assert "Leaderboard" in result.stdout


def test_cli_rna_mfe(tmp_path: Path):
    fasta = tmp_path / "test.fna"
    fasta.write_text(">seq1\nGGGAAACCC\n", encoding="utf-8")
    result = run_cli("rna", "mfe", "--fasta", str(fasta), "--dotbracket", str(tmp_path / "out.dbn"))
    assert "dotbracket" in result.stdout


def test_cli_triage_json(tmp_path: Path):
    json_path = tmp_path / "report.json"
    result = run_cli(
        "triage",
        "--sequence",
        "AUGGCCUUUUAA",
        "--k",
        "3",
        "--max-diff",
        "0",
        "--min-orf-length",
        "9",
        "--top",
        "1",
        "--json",
        str(json_path),
    )
    assert json_path.exists()
    assert "JSON report saved" in result.stdout
    payload = json.loads(json_path.read_text())
    assert payload["orfs"]


def test_cli_viz_triage(tmp_path: Path):
    json_path = tmp_path / "triage.json"
    json_path.write_text(
        json.dumps(
            {
                "sequence": "AUGGCCUUUUAA",
                "skew": [0, 1, 0],
                "clusters": [{"canonical": "AUG", "count": 2, "patterns": ["AUG"], "positions": [0, 3]}],
                "orfs": [
                    {"start": 0, "end": 12, "strand": "+", "frame": 0, "length_nt": 12, "length_aa": 4, "peptide": "MAF"}
                ],
            }
        ),
        encoding="utf-8",
    )
    output_path = tmp_path / "triage.png"
    result = run_cli("viz", "triage", "--json", str(json_path), "--output", str(output_path), "--top", "1")
    assert output_path.exists()
    assert "Triage visualization saved" in result.stdout


def test_cli_viz_schema_command():
    result = run_cli("viz", "schema")
    assert "viz_minimizers" in result.stdout
    result = run_cli("viz", "schema", "--kind", "viz_minimizers")
    assert "properties" in result.stdout


def test_cli_workflows_runner(tmp_path: Path):
    pytest.importorskip("yaml")
    config = tmp_path / "workflow.yaml"
    config.write_text(
        """
workflows:
  - name: test_run
    steps:
      - command: dna
        args:
          sequence: ACGTACGT
          k: 2
          max_diff: 0
        stdout: dna.txt
      - command: triage
        args:
          sequence: AUGGCCUUUUAA
          k: 3
          max_diff: 0
          min_orf_length: 9
        stdout: triage.txt
""",
        encoding="utf-8",
    )
    output_dir = tmp_path / "runs"
    result = run_cli("workflows", "--config", str(config), "--output-dir", str(output_dir))
    assert "Workflow 'test_run' completed" in result.stdout
    assert (output_dir / "test_run" / "dna.txt").exists()
    assert (output_dir / "test_run" / "triage.txt").exists()


def test_cli_workflow_schema_validation(tmp_path: Path):
    pytest.importorskip("yaml")
    ref = tmp_path / "ref.fna"
    reads = tmp_path / "reads.fna"
    ref.write_text(">ref\nACGTACGTACGT\n", encoding="utf-8")
    reads.write_text(">read1\nACGTACGTACGT\n", encoding="utf-8")
    config = tmp_path / "workflow_schema.yaml"
    config.write_text(
        f"""
workflows:
  - name: schema_run
    steps:
      - command: ["seed", "map"]
        args:
          ref: {ref}
          reads: {reads}
          k: 3
          window: 2
          json: map.json
        schema:
          kind: viz_alignment_ribbon
          output: map.json
""",
        encoding="utf-8",
    )
    output_dir = tmp_path / "runs"
    result = run_cli(
        "workflows",
        "--config",
        str(config),
        "--output-dir",
        str(output_dir),
        "--with-schema",
        "--as-json",
    )
    stdout = result.stdout
    start = stdout.find("[\n")
    assert start != -1
    schema_info = json.loads(stdout[start:])
    assert schema_info[0]["workflow"] == "schema_run"
    assert schema_info[0]["steps"][0]["schema_kind"] == "viz_alignment_ribbon"
    artifact = output_dir / "schema_run" / "map.json"
    payload = json.loads(artifact.read_text())
    assert payload["meta"]["spec_version"] == SPEC_VERSION


def test_cli_demo_viz(tmp_path: Path):
    output_dir = tmp_path / "demo"
    result = run_cli("demo", "viz", "--output", str(output_dir))
    assert "Demo visualizations written" in result.stdout
    assert (output_dir / "minimizers.png").exists()
    assert (output_dir / "minimizers.viz.json").exists()


def test_cli_viz_hydropathy(tmp_path: Path):
    pytest.importorskip("Bio")
    output = tmp_path / "hydro.png"
    protein_path = datasets.get_path("protein/demo_protein.faa")
    result = run_cli(
        "viz",
        "hydropathy",
        "--input",
        str(protein_path),
        "--window",
        "11",
        "--output",
        str(output),
    )
    assert output.exists()
    assert "Hydropathy chart saved" in result.stdout


def test_cli_schema_manifest(tmp_path: Path):
    out = tmp_path / "schemas.json"
    result = run_cli("schema", "manifest", "--out", str(out))
    assert out.exists()
    payload = json.loads(out.read_text())
    assert payload["spec_version"] == SPEC_VERSION
    assert "viz_minimizers" in payload["schemas"]


def test_cli_schema_diff_identity(tmp_path: Path):
    manifest_path = tmp_path / "manifest.json"
    manifest_data = run_cli("schema", "manifest").stdout
    manifest_path.write_text(manifest_data, encoding="utf-8")
    result = run_cli("schema", "diff", "--base", str(manifest_path))
    assert "No schema changes detected." in result.stdout


def test_cli_schema_diff_changes(tmp_path: Path):
    manifest_data = json.loads(run_cli("schema", "manifest").stdout)
    schemas = manifest_data["schemas"]
    removed_key = next(iter(schemas.keys()))
    schemas.pop(removed_key)
    base_path = tmp_path / "base.json"
    base_path.write_text(json.dumps(manifest_data, indent=2), encoding="utf-8")
    result = run_cli("schema", "diff", "--base", str(base_path))
    assert "Added schemas" in result.stdout


def test_cli_viz_schema_flag():
    result = run_cli("viz", "alignment-ribbon", "--schema")
    assert "Sample payload" in result.stdout
    assert "demo_read" in result.stdout


def test_cli_rna_ensemble(tmp_path: Path):
    fasta = tmp_path / "test.fna"
    fasta.write_text(">seq1\nAUGCUA\n", encoding="utf-8")
    output = tmp_path / "mea.json"
    result = run_cli(
        "rna",
        "ensemble",
        "--fasta",
        str(fasta),
        "--gamma",
        "1.0",
        "--json",
        str(output),
    )
    assert output.exists()
    assert "structure" in result.stdout


def test_cli_pcr_dag(tmp_path: Path):
    genome = "src/helix/datasets/dna/plasmid_demo.fna"
    out_path = tmp_path / "pcr_dag.json"
    result = run_cli(
        "pcr",
        "dag",
        "--genome",
        genome,
        "--primer-config",
        "examples/pcr_primers.json",
        "--pcr-config",
        "examples/pcr_config.json",
        "--out",
        str(out_path),
    )
    assert out_path.exists()
    assert out_path.exists()
    payload = json.loads(out_path.read_text())
    assert payload["artifact"] == "helix.pcr.amplicon_dag.v1"
    assert payload["nodes"]
