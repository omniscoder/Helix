from __future__ import annotations

import pytest

from helix.veribiota.checks import build_lean_check, validate_lean_check, LeanCheckError


@pytest.fixture()
def sample_payload() -> dict:
    return {
        "root_id": "root",
        "nodes": {
            "root": {
                "log_prob": 0.0,
                "metadata": {"time_step": 0},
                "seq_hashes": {"chr1": "aaa"},
            },
            "child": {
                "log_prob": 0.0,
                "metadata": {"time_step": 1},
                "parent_ids": ["root"],
                "seq_hashes": {"chr1": "bbb"},
            },
        },
        "edges": [
            {"source": "root", "target": "child", "rule": "noop"},
        ],
    }


def test_build_and_validate_lean_check(sample_payload: dict) -> None:
    summary = build_lean_check(sample_payload, dag_name="sample")
    assert summary["node_count"] == 2
    assert summary["edge_count"] == 1
    assert summary["terminal_nodes"] == ["child"]
    validate_lean_check(summary)


def test_validate_lean_check_detects_cycle(sample_payload: dict) -> None:
    summary = build_lean_check(sample_payload, dag_name="bad")
    # introduce a cycle manually
    summary["edges"].append({"source": "child", "target": "root"})
    with pytest.raises(LeanCheckError):
        validate_lean_check(summary)


def test_validate_detects_probability_mismatch(sample_payload: dict) -> None:
    summary = build_lean_check(sample_payload, dag_name="bad_probs")
    summary["nodes"]["child"]["log_prob"] = 1.0
    with pytest.raises(LeanCheckError):
        validate_lean_check(summary)
