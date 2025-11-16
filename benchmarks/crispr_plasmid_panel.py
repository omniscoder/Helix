"""Benchmark CRISPR cut/repair spectra against plasmid library counts.

This harness is deliberately pragmatic:

  * Input: a YAML panel spec (see templates/lab_crispr_plasmid_bench.helix.yml)
  * For each (plasmid target, guide):
      - run helix.crispr.simulate.simulate_cut_repair on the reference window
      - optionally load an amplicon counts table and normalize it
      - compute simple distribution metrics (L1 distance, Jensen–Shannon)
  * Output: a JSON blob summarizing per-guide and aggregate metrics.

All sequences are treated as digital DNA strings only; no wet-lab protocols
or experimental conditions are modeled here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import yaml

from helix import bioinformatics
from helix.crispr.simulate import resolve_crispr_priors, simulate_cut_repair
from helix.prime.model import PegRNA
from helix.prime.priors import resolve_prime_priors
from helix.prime.simulator import simulate_prime_edit

PANEL_KIND = "helix.crispr.plasmid_panel.v1"


@dataclass(frozen=True)
class PlasmidGuideSimulationSpec:
    draws: int
    priors_profile: str
    seed: Optional[int]


@dataclass(frozen=True)
class PlasmidPrimeSimulationSpec:
    draws: int
    priors_profile: str
    seed: Optional[int]


@dataclass(frozen=True)
class PlasmidGuideMeasurementSpec:
    counts_file: Path
    label_column: str = "label"
    count_column: str = "count"
    timepoint_h: Optional[float] = None


@dataclass(frozen=True)
class PlasmidGuideSpec:
    id: str
    sequence: str
    pam_profile: str
    simulation: PlasmidGuideSimulationSpec
    measurements: Optional[PlasmidGuideMeasurementSpec]


@dataclass(frozen=True)
class PlasmidTargetSpec:
    id: str
    reference_sequence: str
    guides: List[PlasmidGuideSpec]
    prime_pegs: List["PlasmidPegSpec"]


@dataclass(frozen=True)
class PlasmidPrimeMeasurementSpec:
    counts_file: Path
    label_column: str = "label"
    count_column: str = "count"
    timepoint_h: Optional[float] = None


@dataclass(frozen=True)
class PlasmidPegSpec:
    id: str
    spacer: str
    pbs: str
    rtt: str
    simulation: PlasmidPrimeSimulationSpec
    measurements: Optional[PlasmidPrimeMeasurementSpec]


def _parse_prime_peg(base_dir: Path, raw: Mapping[str, Any]) -> PlasmidPegSpec:
    peg_id = str(raw.get("id") or "").strip() or "peg"
    spacer = bioinformatics.normalize_sequence(str(raw.get("spacer") or ""))
    pbs = bioinformatics.normalize_sequence(str(raw.get("pbs") or ""))
    rtt = bioinformatics.normalize_sequence(str(raw.get("rtt") or ""))
    if not spacer or not pbs or not rtt:
        raise ValueError(f"Prime peg '{peg_id}' requires non-empty spacer/pbs/rtt.")

    sim_section = raw.get("simulation") or {}
    if not isinstance(sim_section, Mapping):
        sim_section = {}
    sim_spec = _parse_prime_simulation_defaults(sim_section)

    measurements_raw = raw.get("measurements")
    measurements: Optional[PlasmidPrimeMeasurementSpec]
    if isinstance(measurements_raw, Mapping) and measurements_raw.get("counts_file"):
        counts_path = _resolve_path(base_dir, str(measurements_raw["counts_file"]))
        label_col = str(measurements_raw.get("label_column") or "label")
        count_col = str(measurements_raw.get("count_column") or "count")
        t_raw = measurements_raw.get("timepoint_h")
        timepoint_h = float(t_raw) if t_raw is not None else None
        measurements = PlasmidPrimeMeasurementSpec(
            counts_file=counts_path,
            label_column=label_col,
            count_column=count_col,
            timepoint_h=timepoint_h,
        )
    else:
        measurements = None

    return PlasmidPegSpec(
        id=peg_id,
        spacer=spacer,
        pbs=pbs,
        rtt=rtt,
        simulation=sim_spec,
        measurements=measurements,
    )


@dataclass(frozen=True)
class PlasmidPanelSpec:
    name: str
    description: Optional[str]
    kind: str
    defaults: PlasmidGuideSimulationSpec
    targets: List[PlasmidTargetSpec]


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def _load_fasta_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    seq_lines = [line for line in lines if not line.startswith(">")]
    return bioinformatics.normalize_sequence("".join(seq_lines))


def _parse_simulation_defaults(raw: Mapping[str, Any]) -> PlasmidGuideSimulationSpec:
    draws = int(raw.get("draws", 2000))
    priors_profile = str(raw.get("priors_profile") or "default_indel").strip() or "default_indel"
    seed_raw = raw.get("seed")
    seed: Optional[int]
    if seed_raw is None or seed_raw == "":
        seed = None
    else:
        seed = int(seed_raw)
    return PlasmidGuideSimulationSpec(draws=draws, priors_profile=priors_profile, seed=seed)


def _parse_prime_simulation_defaults(raw: Mapping[str, Any]) -> PlasmidPrimeSimulationSpec:
    draws = int(raw.get("draws", 2000))
    priors_profile = str(raw.get("priors_profile") or "default_indel").strip() or "default_indel"
    seed_raw = raw.get("seed")
    seed: Optional[int]
    if seed_raw is None or seed_raw == "":
        seed = None
    else:
        seed = int(seed_raw)
    return PlasmidPrimeSimulationSpec(draws=draws, priors_profile=priors_profile, seed=seed)


def _parse_guide(
    base_dir: Path,
    target_defaults: PlasmidGuideSimulationSpec,
    raw: Mapping[str, Any],
) -> PlasmidGuideSpec:
    guide_id = str(raw.get("id") or "").strip() or "guide"
    sequence = bioinformatics.normalize_sequence(str(raw.get("sequence") or ""))
    if not sequence:
        raise ValueError(f"Guide '{guide_id}' is missing a valid 'sequence'.")
    pam_profile = str(raw.get("pam_profile") or "SpCas9_NGG").strip() or "SpCas9_NGG"

    sim_section = raw.get("simulation") or {}
    if not isinstance(sim_section, Mapping):
        sim_section = {}
    sim_spec = _parse_simulation_defaults(
        {
            "draws": sim_section.get("draws", target_defaults.draws),
            "priors_profile": sim_section.get("priors_profile", target_defaults.priors_profile),
            "seed": sim_section.get("seed", target_defaults.seed),
        }
    )

    measurements_raw = raw.get("measurements")
    measurements: Optional[PlasmidGuideMeasurementSpec]
    if isinstance(measurements_raw, Mapping) and measurements_raw.get("counts_file"):
        counts_path = _resolve_path(base_dir, str(measurements_raw["counts_file"]))
        label_col = str(measurements_raw.get("label_column") or "label")
        count_col = str(measurements_raw.get("count_column") or "count")
        t_raw = measurements_raw.get("timepoint_h")
        timepoint_h = float(t_raw) if t_raw is not None else None
        measurements = PlasmidGuideMeasurementSpec(
            counts_file=counts_path,
            label_column=label_col,
            count_column=count_col,
            timepoint_h=timepoint_h,
        )
    else:
        measurements = None

    return PlasmidGuideSpec(
        id=guide_id,
        sequence=sequence,
        pam_profile=pam_profile,
        simulation=sim_spec,
        measurements=measurements,
    )


def _parse_target(base_dir: Path, defaults: PlasmidGuideSimulationSpec, raw: Mapping[str, Any]) -> PlasmidTargetSpec:
    target_id = str(raw.get("id") or "").strip() or "target"

    ref_seq = str(raw.get("reference_sequence") or "").strip()
    ref_fasta = raw.get("reference_fasta")
    if ref_seq:
        reference_sequence = bioinformatics.normalize_sequence(ref_seq)
    elif ref_fasta:
        ref_path = _resolve_path(base_dir, str(ref_fasta))
        if not ref_path.exists():
            raise FileNotFoundError(f"reference_fasta for target '{target_id}' not found: {ref_path}")
        reference_sequence = _load_fasta_sequence(ref_path)
    else:
        raise ValueError(
            f"Target '{target_id}' requires either 'reference_sequence' or 'reference_fasta'."
        )
    if not reference_sequence:
        raise ValueError(f"Target '{target_id}' has an empty reference sequence.")

    guides_raw = raw.get("guides") or []
    if not isinstance(guides_raw, list) or not guides_raw:
        raise ValueError(f"Target '{target_id}' requires at least one guide.")
    guides = [_parse_guide(base_dir, defaults, entry) for entry in guides_raw if isinstance(entry, Mapping)]

    prime_pegs_raw = raw.get("prime_pegs") or []
    prime_pegs: List[PlasmidPegSpec] = []
    if isinstance(prime_pegs_raw, list):
        for entry in prime_pegs_raw:
            if not isinstance(entry, Mapping):
                continue
            prime_pegs.append(_parse_prime_peg(base_dir, entry))

    return PlasmidTargetSpec(
        id=target_id,
        reference_sequence=reference_sequence,
        guides=guides,
        prime_pegs=prime_pegs,
    )


def load_panel_spec(path: Path) -> PlasmidPanelSpec:
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Panel config '{cfg_path}' not found.")
    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError("Plasmid panel config must be a YAML mapping.")

    kind = str(data.get("kind") or "").strip()
    if kind and kind != PANEL_KIND:
        raise ValueError(f"Unexpected panel kind '{kind}'. Expected '{PANEL_KIND}'.")

    name = str(data.get("name") or "unnamed_plasmid_panel")
    description_raw = data.get("description")
    description = str(description_raw).strip() or None if description_raw is not None else None

    defaults_raw = data.get("defaults") or {}
    if not isinstance(defaults_raw, Mapping):
        defaults_raw = {}
    defaults = _parse_simulation_defaults(defaults_raw)

    targets_raw = data.get("targets") or []
    if not isinstance(targets_raw, list) or not targets_raw:
        raise ValueError("Plasmid panel config requires a non-empty 'targets' list.")

    base_dir = cfg_path.parent
    targets = [_parse_target(base_dir, defaults, entry) for entry in targets_raw if isinstance(entry, Mapping)]

    return PlasmidPanelSpec(
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


def _load_counts(meas: PlasmidGuideMeasurementSpec) -> Optional[Dict[str, float]]:
    if not meas.counts_file.exists():
        return None
    counts: Dict[str, int] = {}
    with meas.counts_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return None
        for row in reader:
            label = str(row.get(meas.label_column) or "").strip()
            if not label:
                continue
            raw_count = row.get(meas.count_column)
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


def _simulate_crispr_distribution(
    reference_sequence: str,
    guide: PlasmidGuideSpec,
) -> Dict[str, float]:
    seq = bioinformatics.normalize_sequence(reference_sequence)
    if not seq:
        raise ValueError("Reference sequence is empty after normalization.")

    guide_len = len(guide.sequence)
    if guide_len == 0 or guide_len > len(seq):
        raise ValueError(
            f"Guide '{guide.id}' length {guide_len} invalid for reference length {len(seq)}."
        )

    # For plasmid windows we anchor the guide at the left of the window by default.
    guide_start = 0
    guide_end = guide_start + guide_len

    guide_dict: Dict[str, Any] = {
        "id": guide.id,
        "sequence": guide.sequence,
        "start": guide_start,
        "end": guide_end,
        "strand": "+",
        "gc_content": bioinformatics.gc_content(guide.sequence),
    }

    priors = resolve_crispr_priors(guide.simulation.priors_profile)
    payload = simulate_cut_repair(
        seq,
        guide_dict,
        priors,
        draws=guide.simulation.draws,
        seed=guide.simulation.seed,
        emit_sequence=False,
        pam_mask=None,
    )
    outcomes = payload.get("outcomes") or []
    sim_probs: Dict[str, float] = {}
    for entry in outcomes:
        label = str(entry.get("label") or "").strip() or "outcome"
        prob = float(entry.get("probability") or 0.0)
        sim_probs[label] = sim_probs.get(label, 0.0) + prob
    # Normalize defensively.
    total = sum(sim_probs.values())
    if total > 0.0:
        sim_probs = {label: value / total for label, value in sim_probs.items()}
    return sim_probs


def _simulate_prime_distribution(
    reference_sequence: str,
    peg: PlasmidPegSpec,
) -> Dict[str, float]:
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


def run_panel(panel: PlasmidPanelSpec) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    aggregate_l1: List[float] = []
    aggregate_js: List[float] = []

    for target in panel.targets:
        # CRISPR guides.
        for guide in target.guides:
            entry: Dict[str, Any] = {
                "mechanism": "crispr",
                "target_id": target.id,
                "guide_id": guide.id,
                "pam_profile": guide.pam_profile,
                "draws": guide.simulation.draws,
                "priors_profile": guide.simulation.priors_profile,
                "seed": guide.simulation.seed,
            }
            status = "ok"
            metrics: Dict[str, Any] = {}
            sim_probs: Dict[str, float] = {}
            obs_probs: Optional[Dict[str, float]] = None

            try:
                sim_probs = _simulate_crispr_distribution(target.reference_sequence, guide)
            except Exception as exc:  # pragma: no cover - defensive
                status = "simulation_error"
                metrics["error"] = str(exc)

            if guide.measurements is not None:
                obs_probs = _load_counts(guide.measurements)
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
                    if guide.measurements.timepoint_h is not None:
                        metrics["timepoint_h"] = guide.measurements.timepoint_h

            entry["status"] = status
            entry["simulated"] = sim_probs
            entry["observed"] = obs_probs
            entry["metrics"] = metrics
            results.append(entry)

        # Prime editing pegs.
        for peg in target.prime_pegs:
            entry = {
                "mechanism": "prime",
                "target_id": target.id,
                "peg_id": peg.id,
                "draws": peg.simulation.draws,
                "priors_profile": peg.simulation.priors_profile,
                "seed": peg.simulation.seed,
            }
            status = "ok"
            metrics = {}
            sim_probs: Dict[str, float] = {}
            obs_probs = None

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
        "pairs": len(results),
    }
    if aggregate_l1:
        summary["l1_mean"] = float(sum(aggregate_l1) / len(aggregate_l1))
    if aggregate_js:
        summary["js_mean"] = float(sum(aggregate_js) / len(aggregate_js))

    payload: Dict[str, Any] = {
        "schema": {"kind": "crispr_plasmid_panel", "spec_version": "1.0"},
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
        description="Compare Helix CRISPR cut/repair spectra to plasmid library counts."
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Plasmid panel YAML config (see templates/lab_crispr_plasmid_bench.helix.yml).",
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
    print(f"[OK] CRISPR plasmid panel results written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
