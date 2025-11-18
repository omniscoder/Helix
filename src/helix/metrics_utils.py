from __future__ import annotations

from typing import Any, Dict, Mapping, Sequence

from helix.studio.run_metrics import RunMetrics


def run_metrics_to_summary(metrics: RunMetrics) -> Dict[str, Any]:
    counts = {
        "intended": metrics.intended_count,
        "off_target": metrics.off_target_count,
        "neutral": metrics.neutral_count,
    }
    probability_mass = {
        "intended": metrics.intended_mass,
        "off_target": metrics.off_target_mass,
        "neutral": metrics.neutral_mass,
    }
    summary = {
        "summary": dict(metrics.summary),
        "counts": counts,
        "probability_mass": probability_mass,
        "dag": dict(metrics.dag_stats),
        "top_outcomes": _normalize_top(metrics.top_outcomes),
        "perf": dict(metrics.perf),
        "config": dict(metrics.config),
    }
    if metrics.is_pcr:
        summary["pcr"] = {
            "amplicon_length": metrics.pcr_amplicon_length,
            "final_mass_ng": metrics.pcr_final_mass_ng,
            "final_copies": metrics.pcr_final_copies,
            "final_mutation_rate": metrics.pcr_final_mutation_rate,
            "cycles": metrics.pcr_cycles,
        }
    return summary


def _normalize_top(entries: Sequence[tuple[str, float]]) -> Sequence[tuple[str, float]]:
    return [(label, prob) for label, prob in entries]
