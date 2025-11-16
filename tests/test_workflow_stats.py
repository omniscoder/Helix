from __future__ import annotations

from helix.studio.workflow_stats import summarize_workflow_meta


def test_workflow_stats_summary_contains_expected_lines() -> None:
    meta = {
        "inst_events": [1, 2, 3],
        "heat_values": [0.1, 0.2],
        "captured_mass": 25.0,
        "total_mass": 100.0,
        "x_extent": 42.5,
    }
    summary = summarize_workflow_meta(meta)
    assert any("Bars" in line for line in summary)
    assert any("Heat" in line for line in summary)
    assert any("Captured mass" in line for line in summary)
