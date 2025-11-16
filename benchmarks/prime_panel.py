"""Benchmark prime editing outcome spectra against synthetic counts.

This harness treats each entry as a prime editing experiment over a digital
reference window. For each target it:

  * Runs helix.prime.simulator.simulate_prime_edit on the reference window.
  * Aggregates outcome probabilities by label.
  * Optionally loads a counts table and compares distributions.

All sequences are purely computational DNA strings; no wet-lab conditions
are modeled here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from helix import bioinformatics
from helix.prime.model import PegRNA
from helix.prime.priors import resolve_prime_priors
from helix.prime.simulator import simulate_prime_edit

PANEL_KIND = "helix.prime.editing_panel.v1"


@dataclass(frozen=True)
class PrimeSimulationSpec:
    draws: int
    priors_profile: str
    seed: Optional[int]


@dataclass(frozen=True)
class PrimeMeasurementSpec:
    counts_file: Path
    label_column: str
    count_column: str
    timepoint_h: Optional[float]


@dataclass(frozen=True)
class PrimePegSpec:
    id: str
    spacer: str
    pbs: str
    rtt: str
    simulation: PrimeSimulationSpec
    measurements: Optional[PrimeMeasurementSpec]


@dataclass(frozen=True)
class PrimeTargetSpec:
    id: str
    reference_sequence: str
    peg: PrimePegSpec


@dataclass(frozen=True)
class PrimePanelSpec:
    name: str
    description: Optional[str]
    kind: str
    defaults: PrimeSimulationSpec
    targets: List[PrimeTargetSpec]


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def _parse_simulation_defaults(raw: Mapping[str, Any]) -> PrimeSimulationSpec:
    draws = int(raw.get("draws", 2000))
    priors_profile = str(raw.get("priors_profile") or "default_indel").strip() or "default_indel"
    seed_raw = raw.get("seed")
    seed: Optional[int]
    if seed_raw is None or str(seed_raw).strip() == "":
        seed = None
    else:
        seed = int(seed_raw)
    return PrimeSimulationSpec(draws=draws, priors_profile=priors_profile, seed=seed)


def _parse_measurements(base_dir: Path, raw: Mapping[str, Any]) -> PrimeMeasurementSpec:
    counts_file = _resolve_path(base_dir, str(raw.get("counts_file")))
    label_column = str(raw.get("label_column") or "label")
    count_column = str(raw.get("count_column") or "count")
    t_raw = raw.get("timepoint_h")
    timepoint_h = float(t_raw) if t_raw is not None else None
    return PrimeMeasurementSpec(
        counts_file=counts_file,
        label_column=label_column,
        count_column=count_column,
        timepoint_h=timepoint_h,
    )


def _parse_peg(
    base_dir: Path,
    defaults: PrimeSimulationSpec,
    raw: Mapping[str, Any],
) -> PrimePegSpec:
    peg_id = str(raw.get("id") or "").strip() or "peg"
    spacer = bioinformatics.normalize_sequence(str(raw.get("spacer") or ""))
    pbs = bioinformatics.normalize_sequence(str(raw.get("pbs") or ""))
    rtt = bioinformatics.normalize_sequence(str(raw.get("rtt") or ""))
    if not spacer or not pbs or not rtt:
        raise ValueError(f"Peg '{peg_id}' requires non-empty spacer/pbs/rtt.")

    sim_raw = raw.get("simulation") or {}
    if not isinstance(sim_raw, Mapping):
        sim_raw = {}
    sim_spec = _parse_simulation_defaults(
        {
            "draws": sim_raw.get("draws", defaults.draws),
            "priors_profile": sim_raw.get("priors_profile", defaults.priors_profile),
            "seed": sim_raw.get("seed", defaults.seed),
        }
    )

    meas_raw = raw.get("measurements")
    measurements: Optional[PrimeMeasurementSpec]
    if isinstance(meas_raw, Mapping) and meas_raw.get("counts_file"):
        measurements = _parse_measurements(base_dir, meas_raw)
    else:
        measurements = None

    return PrimePegSpec(
        id=peg_id,
        spacer=spacer,
        pbs=pbs,
        rtt=rtt,
        simulation=sim_spec,
        measurements=measurements,
    )


def _parse_target(base_dir: Path, defaults: PrimeSimulationSpec, raw: Mapping[str, Any]) -> PrimeTargetSpec:
    target_id = str(raw.get("id") or "").strip() or "target"
    ref_seq = str(raw.get("reference_sequence") or "").strip()
    reference_sequence = bioinformatics.normalize_sequence(ref_seq)
    if not reference_sequence:
        raise ValueError(f"Target '{target_id}' requires a non-empty 'reference_sequence'.")

    peg_raw = raw.get("peg") or {}
    if not isinstance(peg_raw, Mapping):
        raise ValueError(f"Target '{target_id}' requires a 'peg' mapping.")
    peg_spec = _parse_peg(base_dir, defaults, peg_raw)

    return PrimeTargetSpec(id=target_id, reference_sequence=reference_sequence, peg=peg_spec)


def load_panel_spec(path: Path) -> PrimePanelSpec:
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Prime editing panel config '{cfg_path}' not found.")
    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError("Prime editing panel config must be a YAML mapping.")

    kind = str(data.get("kind") or "").strip()
    if kind and kind != PANEL_KIND:
        raise ValueError(f"Unexpected panel kind '{kind}'. Expected '{PANEL_KIND}'.")

    name = str(data.get("name") or "unnamed_prime_panel")
    desc_raw = data.get("description")
    description = str(desc_raw).strip() or None if desc_raw is not None else None

    defaults_raw = data.get("defaults") or {}
    if not isinstance(defaults_raw, Mapping):
        defaults_raw = {}
    defaults = _parse_simulation_defaults(defaults_raw)

    targets_raw = data.get("targets") or []
    if not isinstance(targets_raw, list) or not targets_raw:
        raise ValueError("Prime editing panel config requires a non-empty 'targets' list.")
    base_dir = cfg_path.parent
    targets = [_parse_target(base_dir, defaults, entry) for entry in targets_raw if isinstance(entry, Mapping)]

    return PrimePanelSpec(
        name=name,
        description=description,
        kind=PANEL_KIND,
        defaults=defaults,
        targets=targets,
    )


def _normalize_counts(counts: Mapping[str, int]) -> Dict[str, float]:
    total = float(sum(max(0, int(v)) for v in counts.values()))
    if total <= 0.0:
        return {label: 0.0 for label in counts}
    return {label: int(v) / total for label, v in counts.items()}


def _load_counts(spec: PrimeMeasurementSpec) -> Optional[Dict[str, float]]:
    if not spec.counts_file.exists():
        return None
    counts: Dict[str, int] = {}
    with spec.counts_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return None
        for row in reader:
            label = str(row.get(spec.label_column) or "").strip()
            if not label:
                continue
            raw_count = row.get(spec.count_column)
            if raw_count is None or str(raw_count).strip() == "":
                continue
            try:
                count = int(raw_count)
            except ValueError:
                continue
            counts[label] = counts.get(label, 0) + max(0, count)
    if not counts:
        return None
    return _normalize_counts(counts)


def _aligned_distributions(
    sim_probs: Mapping[str, float],
    obs_probs: Mapping[str, float],
) -> Tuple[List[float], List[float], List[str]]:
    labels = sorted(set(sim_probs) | set(obs_probs))
    p: List[float] = []
    q: List[float] = []
    for label in labels:
        p.append(float(sim_probs.get(label, 0.0)))
        q.append(float(obs_probs.get(label, 0.0)))
    return p, q, labels


def _l1_distance(p: Sequence[float], q: Sequence[float]) -> float:
    return float(sum(abs(a - b) for a, b in zip(p, q)))


def _js_divergence(p: Sequence[float], q: Sequence[float], eps: float = 1e-12) -> float:
    def _kl(a: Sequence[float], b: Sequence[float]) -> float:
        total = 0.0
        for ai, bi in zip(a, b):
            ai = max(eps, float(ai))
            bi = max(eps, float(bi))
            total += ai * math.log(ai / bi)
        return total

    m = [(float(ai) + float(bi)) * 0.5 for ai, bi in zip(p, q)]
    return 0.5 * (_kl(p, m) + _kl(q, m))


def _simulate_prime_distribution(reference_sequence: str, peg: PrimePegSpec) -> Dict[str, float]:
    seq = bioinformatics.normalize_sequence(reference_sequence)
    if not seq:
        raise ValueError("Reference sequence is empty after normalization.")

    priors = resolve_prime_priors(peg.simulation.priors_profile) if peg.simulation.priors_profile else None
    peg_obj = PegRNA(spacer=peg.spacer, pbs=peg.pbs, rtt=peg.rtt, name=peg.id)
    payload = simulate_prime_edit(
        seq,
        peg_obj,
        priors=priors,
        draws=peg.simulation.draws,
        seed=peg.simulation.seed,
        emit_sequence=False,
        pam_mask=None,
    )
    outcomes = payload.get("outcomes") or []
    sim_probs: Dict[str, float] = {}
    for entry in outcomes:
        label = str(entry.get("label") or "").strip() or "outcome"
        prob = float(entry.get("probability") or 0.0)
        sim_probs[label] = sim_probs.get(label, 0.0) + prob
    total = sum(sim_probs.values())
    if total > 0.0:
        sim_probs = {label: value / total for label, value in sim_probs.items()}
    return sim_probs


def run_panel(panel: PrimePanelSpec) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    aggregate_l1: List[float] = []
    aggregate_js: List[float] = []

    for target in panel.targets:
        peg = target.peg
        entry: Dict[str, Any] = {
            "target_id": target.id,
            "peg_id": peg.id,
            "draws": peg.simulation.draws,
            "priors_profile": peg.simulation.priors_profile,
            "seed": peg.simulation.seed,
        }
        status = "ok"
        metrics: Dict[str, Any] = {}
        sim_probs: Dict[str, float] = {}
        obs_probs: Optional[Dict[str, float]] = None

        try:
            sim_probs = _simulate_prime_distribution(target.reference_sequence, peg)
        except Exception as exc:  # pragma: no cover - defensive
            status = "simulation_error"
            metrics["error"] = str(exc)

        if peg.measurements is not None:
            obs_probs = _load_counts(peg.measurements)
            if obs_probs is None:
                if status == "ok":
                    status = "missing_counts"
            elif status == "ok":
                p, q, _labels = _aligned_distributions(sim_probs, obs_probs)
                l1 = _l1_distance(p, q)
                js = _js_divergence(p, q)
                metrics["l1_distance"] = l1
                metrics["js_divergence"] = js
                aggregate_l1.append(l1)
                aggregate_js.append(js)
                if peg.measurements.timepoint_h is not None:
                    metrics["timepoint_h"] = peg.measurements.timepoint_h

        entry["status"] = status
        entry["simulated"] = sim_probs
        entry["observed"] = obs_probs
        entry["metrics"] = metrics
        results.append(entry)

    summary: Dict[str, Any] = {
        "panel_name": panel.name,
        "panel_description": panel.description,
        "targets": len(panel.targets),
        "evaluated": len(results),
    }
    if aggregate_l1:
        summary["l1_mean"] = float(sum(aggregate_l1) / len(aggregate_l1))
    if aggregate_js:
        summary["js_mean"] = float(sum(aggregate_js) / len(aggregate_js))

    payload: Dict[str, Any] = {
        "schema": {"kind": "prime_editing_panel", "spec_version": "1.0"},
        "panel": {
            "name": panel.name,
            "description": panel.description,
            "kind": panel.kind,
        },
        "summary": summary,
        "results": results,
    }
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare Helix prime editing spectra to outcome counts."
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Prime editing panel YAML config (see templates/prime_editing_panel.helix.yml).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output JSON path for benchmark results.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    panel = load_panel_spec(args.config)
    payload = run_panel(panel)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Prime editing panel results written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

