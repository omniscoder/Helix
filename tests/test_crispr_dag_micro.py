"""Ensure the micro-universe DAG matches the brute-force reference."""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.verify_crispr_micro import (  # noqa: E402
    GENOME_PATH,
    PROB_TOL,
    assert_probability_invariants,
    assert_structural_invariants,
    brute_force_crispr_micro,
    collect_terminal_logprobs_from_brute,
    collect_terminal_logprobs_from_payload,
    normalized_from_logs,
    read_fasta,
    run_helix_dag,
)


def test_crispr_micro_dag_matches_bruteforce(tmp_path: Path) -> None:
    payload = run_helix_dag(tmp_path)
    genome = read_fasta(GENOME_PATH)
    brute_dag = brute_force_crispr_micro(genome)

    assert_structural_invariants(payload)
    assert_probability_invariants(payload)

    helix_terminal_logs = collect_terminal_logprobs_from_payload(payload)
    brute_terminal_logs = collect_terminal_logprobs_from_brute(brute_dag)
    assert set(helix_terminal_logs) == set(brute_terminal_logs), "Terminal states diverged."

    helix_probs = normalized_from_logs(helix_terminal_logs)
    brute_probs = normalized_from_logs(brute_terminal_logs)
    for signature, helix_prob in helix_probs.items():
        brute_prob = brute_probs[signature]
        assert math.isclose(helix_prob, brute_prob, rel_tol=PROB_TOL, abs_tol=PROB_TOL)
