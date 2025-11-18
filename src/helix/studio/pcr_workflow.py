from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Mapping

from .session import PCRResult


@dataclass(frozen=True)
class PCRWorkflowBar:
    cycle_index: int
    mass_norm: float
    mutation_norm: float
    mass_ng: float
    mutation_rate: float
    copies: float


@dataclass(frozen=True)
class PCRWorkflowSpec:
    bars: list[PCRWorkflowBar]
    max_mass_ng: float
    max_mutation_rate: float


def build_pcr_workflow_spec(result: PCRResult) -> PCRWorkflowSpec:
    cycles = result.cycles or []
    if not cycles:
        return PCRWorkflowSpec(bars=[], max_mass_ng=0.0, max_mutation_rate=0.0)
    max_mass = max(c.amplicon_mass_ng for c in cycles) or 1e-9
    max_mut = max(c.cumulative_mutation_rate for c in cycles) or 1e-9
    bars = [
        PCRWorkflowBar(
            cycle_index=c.cycle_index,
            mass_norm=c.amplicon_mass_ng / max_mass,
            mutation_norm=c.cumulative_mutation_rate / max_mut,
            mass_ng=c.amplicon_mass_ng,
            mutation_rate=c.cumulative_mutation_rate,
            copies=c.amplicon_copies,
        )
        for c in cycles
    ]
    return PCRWorkflowSpec(bars=bars, max_mass_ng=max_mass, max_mutation_rate=max_mut)


def build_pcr_workflow_meta(result: PCRResult) -> dict[str, Any]:
    return {
        "kind": "PCR",
        "amplicon_length": result.amplicon_length,
        "amplicon_start": result.amplicon_start,
        "amplicon_end": result.amplicon_end,
        "final_mass_ng": result.final_amplicon_mass_ng,
        "final_mutation_rate": result.final_mutation_rate,
        "final_copies": result.final_amplicon_copies,
        "cycles": [
            {
                "cycle": c.cycle_index,
                "mass_ng": c.amplicon_mass_ng,
                "mutation_rate": c.cumulative_mutation_rate,
                "copies": c.amplicon_copies,
            }
            for c in (result.cycles or [])
        ],
    }


def summarize_pcr_workflow_meta(meta: Mapping[str, Any]) -> list[str]:
    cycles = meta.get("cycles") or []
    final_mass = meta.get("final_mass_ng", 0.0)
    final_mut = meta.get("final_mutation_rate", 0.0)
    lines = [
        f"Amplicon length: {meta.get('amplicon_length', 0)} bp",
        f"Amplicon region: {meta.get('amplicon_start', 0)}–{meta.get('amplicon_end', 0)}",
        f"Final mass: {float(final_mass):.2f} ng",
        f"Final mutation rate: {float(final_mut):.4f}",
        f"Cycles captured: {len(cycles)}",
    ]
    if cycles:
        first = cycles[0]
        last = cycles[-1]
        try:
            lines.append(
                f"Cycle 1 mass: {float(first.get('mass_ng', 0.0)):.4f} ng · mutations {float(first.get('mutation_rate', 0.0)):.4f}"
            )
            lines.append(
                f"Cycle {last.get('cycle')} mass: {float(last.get('mass_ng', 0.0)):.4f} ng · mutations {float(last.get('mutation_rate', 0.0)):.4f}"
            )
        except Exception:
            pass
    return lines


__all__ = [
    "build_pcr_workflow_spec",
    "build_pcr_workflow_meta",
    "summarize_pcr_workflow_meta",
]
