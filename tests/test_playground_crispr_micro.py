from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Tuple

import pytest

from helix.cli import _frames_to_artifact_payload
from helix.edit.dag import dag_from_payload

PROB_TOL = 1e-6


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _dag_node_index(payload: dict) -> Dict[str, dict]:
    nodes = payload.get("nodes") or {}
    return {node_id: node for node_id, node in nodes.items()}


def _dag_edge_multiset(payload: dict) -> Dict[Tuple[str, str], str]:
    edges = payload.get("edges") or []
    result: Dict[Tuple[str, str], str] = {}
    for edge in edges:
        key = (edge["source"], edge["target"])
        result[key] = edge.get("rule", "")
    return result


@pytest.mark.parametrize("fixture_dir", ["docs/playground/data"])
def test_playground_crispr_micro_fixtures_match_engine(fixture_dir: str) -> None:
    """Ensure the Playground CRISPR micro fixtures stay in sync with the engine.

    Assumes the Playground uses:
      - crispr_micro.edit_dag.json   (static artifact)
      - crispr_micro.frames.json     (list or dict {frames: [...]})
    """

    base = Path(fixture_dir)
    static_path = base / "crispr_micro.edit_dag.json"
    frames_path = base / "crispr_micro.frames.json"

    assert static_path.is_file(), f"Missing Playground static fixture: {static_path}"
    assert frames_path.is_file(), f"Missing Playground frames fixture: {frames_path}"

    static_payload = _load_json(static_path)
    frames_payload = _load_json(frames_path)

    if isinstance(frames_payload, dict) and "frames" in frames_payload:
        frames = frames_payload["frames"]
    else:
        frames = frames_payload

    assert isinstance(frames, list) and frames, "Playground frames fixture is empty"

    frame_artifact = _frames_to_artifact_payload(frames)

    assert static_payload.get("root_id") == frame_artifact.get("root_id")
    assert static_payload.get("artifact") == frame_artifact.get("artifact")

    static_nodes = _dag_node_index(static_payload)
    frame_nodes = _dag_node_index(frame_artifact)

    # Frame-derived artifact may carry extra nodes; static nodes must be a subset.
    assert set(static_nodes.keys()).issubset(set(frame_nodes.keys()))

    for nid, node_s in static_nodes.items():
        node_f = frame_nodes[nid]
        lp_s = float(node_s.get("log_prob", 0.0))
        lp_f = float(node_f.get("log_prob", 0.0))
        assert math.isclose(lp_s, lp_f, rel_tol=0.0, abs_tol=PROB_TOL)

    edges_s = _dag_edge_multiset(static_payload)
    edges_f = _dag_edge_multiset(frame_artifact)

    assert set(edges_s.keys()).issubset(set(edges_f.keys()))
    for key, rule_s in edges_s.items():
        assert rule_s == edges_f[key]
