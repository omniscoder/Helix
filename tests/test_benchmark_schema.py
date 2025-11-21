"""Benchmark CLI JSON schema tests."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests.test_helix_cli import run_cli


def test_engine_benchmark_json_schema(tmp_path: Path):
    out_path = tmp_path / "bench.json"
    run_cli(
        "engine",
        "benchmark",
        "--backends",
        "cpu-reference",
        "--crispr-shapes",
        "1x16x10",
        "--prime-workloads",
        "4x200x12",
        "--seed",
        "2",
        "--json",
        str(out_path),
    )
    data = json.loads(out_path.read_text())
    assert data["helix_version"]
    assert isinstance(data["seed"], int)
    scoring = data["scoring_versions"]
    assert scoring["crispr"]
    assert scoring["prime"]
    env = data["env"]
    assert "platform" in env and env["platform"]
    assert "python_version" in env
    benches = data["benchmarks"]
    crispr_entries = benches["crispr"]
    prime_entries = benches["prime"]
    assert crispr_entries and prime_entries
    crispr_entry = crispr_entries[0]
    assert crispr_entry["backend_requested"]
    assert "shape" in crispr_entry
    assert isinstance(crispr_entry["mpairs_per_s"], float)
    assert crispr_entry["backend_used"] in {"cpu-reference", "native-cpu", "gpu"}
    prime_entry = prime_entries[0]
    assert prime_entry["backend_used"]
    assert prime_entry["workload"]
    assert isinstance(prime_entry["predictions_per_s"], float)
