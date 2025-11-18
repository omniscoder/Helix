from __future__ import annotations

import pytest

from helix.studio.pcr_engine import simulate_pcr, find_amplicon
from helix.studio.run_metrics import compute_run_metrics, build_report, render_markdown_report
from helix.studio.session import (
    PCRConfig,
    _serialize_pcr_config,
    _serialize_pcr_result,
)


def _basic_config(**overrides):
    cfg = PCRConfig(
        template_seq="AAACCCGGGTTTAAACCCGGGTTT",
        fwd_primer="AAACCC",
        rev_primer="CCCGGG",
        cycles=5,
        polymerase_name="TestPol",
        base_error_rate=1e-5,
        indel_fraction=0.1,
        efficiency=0.9,
        initial_copies=1,
    )
    return cfg if not overrides else cfg.__class__(**{**cfg.__dict__, **overrides})


def test_pcr_generates_amplicon_and_growth() -> None:
    cfg = _basic_config()
    result = simulate_pcr(cfg)
    assert result.amplicon_length > 0
    assert len(result.cycles) == cfg.cycles
    assert result.final_amplicon_copies > cfg.initial_copies
    assert result.final_mutation_rate > 0.0


def test_pcr_invalid_primers_raise() -> None:
    cfg = _basic_config(fwd_primer="AAAAAA", rev_primer="TTTTTT")
    with pytest.raises(ValueError):
        simulate_pcr(cfg)


def test_pcr_metrics_roundtrip() -> None:
    cfg = _basic_config()
    result = simulate_pcr(cfg)
    snapshot = {
        "state": {
            "run_kind": "PCR",
            "run_id": 1,
            "viz_dirty": False,
            "config": {},
            "pcr_config": _serialize_pcr_config(cfg),
            "pcr_result": _serialize_pcr_result(result),
        }
    }
    metrics = compute_run_metrics(snapshot)
    assert metrics.is_pcr is True
    assert metrics.pcr_amplicon_length == result.amplicon_length
    assert metrics.pcr_mass_curve is not None
    assert len(metrics.pcr_mass_curve or []) == cfg.cycles
    assert metrics.pcr_final_mass_ng == pytest.approx(result.final_amplicon_mass_ng)
    assert metrics.pcr_amplicon_start == result.amplicon_start
    assert metrics.pcr_forward_primer is not None


def test_find_amplicon_locations() -> None:
    template = "AAACCCGGGTTTAAACCCGGGTTT"
    info = find_amplicon(template, "AAACCC", "AAACCC")
    assert info.start == 0
    assert info.length > 0


def test_pcr_report_markdown_contains_section() -> None:
    cfg = _basic_config()
    result = simulate_pcr(cfg)
    snapshot = {
        "state": {
            "run_kind": "PCR",
            "run_id": 1,
            "viz_dirty": False,
            "config": {},
            "pcr_config": _serialize_pcr_config(cfg),
            "pcr_result": _serialize_pcr_result(result),
        }
    }
    report = build_report(snapshot)
    assert report["pcr"]["amplicon_length"] == result.amplicon_length
    md = render_markdown_report(report)
    assert "PCR Simulation" in md
    assert str(result.amplicon_length) in md
