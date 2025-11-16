from __future__ import annotations

from helix.studio.gl_panel import build_overlay_lines
from helix.studio.run_metrics import compute_run_metrics


def test_build_overlay_lines_contains_run_stats() -> None:
    snapshot = {
        "state": {
            "run_kind": "CRISPR",
            "run_id": 1,
            "viz_dirty": False,
            "config": {"run_config": {"draws": 1000}},
            "genome_source": "inline",
            "outcomes": [
                {"label": "intended", "tags": ["intent"], "diff": {"kind": "none"}, "probability": 0.6},
                {"label": "off-target", "diff": {"kind": "del"}, "probability": 0.2},
            ],
            "dag_from_runtime": {
                "nodes": [
                    {"children": [1, 2]},
                    {"children": []},
                    {"children": []},
                ],
                "max_depth": 3,
            },
        }
    }
    metrics = compute_run_metrics(snapshot)
    lines = build_overlay_lines(metrics)
    combined = "\n".join(lines)
    assert "Crispr" in combined
    assert "1000" in combined
    assert "branches" in combined
    assert "Intended: 1" in combined
    assert "Prob mass" in combined
