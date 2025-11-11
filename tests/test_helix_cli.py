import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [sys.executable, "helix_cli.py", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
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


def test_cli_rna_nussinov():
    result = run_cli("rna", "--sequence", "GGGAAACCC", "--min-loop", "0")
    assert "Dot-bracket" in result.stdout
    assert "(" in result.stdout


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


def test_cli_viz_hydropathy(tmp_path: Path):
    pytest.importorskip("Bio")
    output = tmp_path / "hydro.png"
    result = run_cli(
        "viz",
        "hydropathy",
        "--input",
        "input/protein/demo_protein.faa",
        "--window",
        "11",
        "--output",
        str(output),
    )
    assert output.exists()
    assert "Hydropathy chart saved" in result.stdout
