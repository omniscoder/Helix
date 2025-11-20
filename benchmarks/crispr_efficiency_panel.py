"""Benchmark CRISPR on-target efficiency scores against edited fractions.

This harness treats each panel entry as an endogenous-like target with a
single guide and an observed edited fraction. It:

  * Computes an on-target score through the CRISPR simulator API
    (which calls the underlying physics engine).
  * Compares predicted scores to observed edited fractions.
  * Reports basic correlation / error metrics in a JSON payload.

All operations are purely computational on digital DNA sequences.
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
from helix.crispr.model import CasSystem, CasSystemType, GuideRNA
from helix.crispr.simulator import (
    EfficiencyTargetRequest as EngineEfficiencyTargetRequest,
    predict_efficiency_for_targets,
)

PANEL_KIND = "helix.crispr.efficiency_panel.v1"


@dataclass(frozen=True)
class CasConfigSpec:
    name: str
    system_type: CasSystemType
    pam_pattern: str
    cut_offset: int
    max_mismatches: int
    weight_mismatch_penalty: float
    weight_pam_penalty: float


@dataclass(frozen=True)
class EfficiencyMeasurementSpec:
    edited_fraction_value: Optional[float]
    edited_fraction_file: Optional[Path]
    edited_fraction_column: str
    chromatin_features: Dict[str, float]


@dataclass(frozen=True)
class EfficiencyTargetSpec:
    id: str
    reference_sequence: str
    guide_id: str
    guide_sequence: str
    guide_strand: str
    measurements: EfficiencyMeasurementSpec


@dataclass(frozen=True)
class EfficiencyPanelSpec:
    name: str
    description: Optional[str]
    kind: str
    cas: CasConfigSpec
    edited_fraction_column_default: str
    targets: List[EfficiencyTargetSpec]


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def _parse_cas_config(raw: Mapping[str, Any]) -> CasConfigSpec:
    name = str(raw.get("name") or "cas_demo")
    system_type_str = str(raw.get("system_type") or "cas9").lower()
    try:
        system_type = CasSystemType(system_type_str)
    except ValueError:
        system_type = CasSystemType.CAS9
    pam_pattern = str(raw.get("pam_pattern") or "NGG").strip() or "NGG"
    cut_offset = int(raw.get("cut_offset", 3))
    max_mismatches = int(raw.get("max_mismatches", 3))
    weight_mismatch_penalty = float(raw.get("weight_mismatch_penalty", 1.0))
    weight_pam_penalty = float(raw.get("weight_pam_penalty", 2.0))
    return CasConfigSpec(
        name=name,
        system_type=system_type,
        pam_pattern=pam_pattern,
        cut_offset=cut_offset,
        max_mismatches=max_mismatches,
        weight_mismatch_penalty=weight_mismatch_penalty,
        weight_pam_penalty=weight_pam_penalty,
    )


def _parse_measurements(
    base_dir: Path,
    defaults: Mapping[str, Any],
    raw: Mapping[str, Any],
) -> EfficiencyMeasurementSpec:
    edited_raw = raw.get("edited_fraction") or {}
    if not isinstance(edited_raw, Mapping):
        edited_raw = {}
    value_raw = edited_raw.get("value")
    edited_value: Optional[float]
    if value_raw is None or str(value_raw).strip() == "":
        edited_value = None
    else:
        edited_value = float(value_raw)

    file_raw = edited_raw.get("file")
    edited_file: Optional[Path]
    if file_raw:
        edited_file = _resolve_path(base_dir, str(file_raw))
    else:
        edited_file = None

    col_default = str(defaults.get("edited_fraction_column") or "edited_fraction")
    col = str(edited_raw.get("column") or col_default)

    chromatin_raw = raw.get("chromatin_features") or {}
    chromatin: Dict[str, float] = {}
    if isinstance(chromatin_raw, Mapping):
        for key, val in chromatin_raw.items():
            try:
                chromatin[str(key)] = float(val)
            except (TypeError, ValueError):
                continue

    return EfficiencyMeasurementSpec(
        edited_fraction_value=edited_value,
        edited_fraction_file=edited_file,
        edited_fraction_column=col,
        chromatin_features=chromatin,
    )


def _parse_target(
    base_dir: Path,
    defaults: Mapping[str, Any],
    raw: Mapping[str, Any],
) -> EfficiencyTargetSpec:
    target_id = str(raw.get("id") or "").strip() or "target"
    ref_seq = str(raw.get("reference_sequence") or "").strip()
    reference_sequence = bioinformatics.normalize_sequence(ref_seq)
    if not reference_sequence:
        raise ValueError(f"Target '{target_id}' requires a non-empty 'reference_sequence'.")

    guide_raw = raw.get("guide") or {}
    if not isinstance(guide_raw, Mapping):
        raise ValueError(f"Target '{target_id}' requires a 'guide' mapping.")
    guide_id = str(guide_raw.get("id") or "guide")
    guide_seq = bioinformatics.normalize_sequence(str(guide_raw.get("sequence") or ""))
    if not guide_seq:
        raise ValueError(f"Guide '{guide_id}' for target '{target_id}' is missing a valid 'sequence'.")
    strand = str(guide_raw.get("strand") or "+")
    if strand not in {"+", "-"}:
        strand = "+"

    meas_raw = raw.get("measurements") or {}
    if not isinstance(meas_raw, Mapping):
        meas_raw = {}
    meas_spec = _parse_measurements(base_dir, defaults, meas_raw)

    return EfficiencyTargetSpec(
        id=target_id,
        reference_sequence=reference_sequence,
        guide_id=guide_id,
        guide_sequence=guide_seq,
        guide_strand=strand,
        measurements=meas_spec,
    )


def load_panel_spec(path: Path) -> EfficiencyPanelSpec:
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Efficiency panel config '{cfg_path}' not found.")
    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError("Efficiency panel config must be a YAML mapping.")

    kind = str(data.get("kind") or "").strip()
    if kind and kind != PANEL_KIND:
        raise ValueError(f"Unexpected panel kind '{kind}'. Expected '{PANEL_KIND}'.")

    name = str(data.get("name") or "unnamed_efficiency_panel")
    desc_raw = data.get("description")
    description = str(desc_raw).strip() or None if desc_raw is not None else None

    cas_raw = data.get("cas") or {}
    if not isinstance(cas_raw, Mapping):
        cas_raw = {}
    cas_spec = _parse_cas_config(cas_raw)

    defaults_raw = data.get("defaults") or {}
    if not isinstance(defaults_raw, Mapping):
        defaults_raw = {}
    edited_fraction_column_default = str(defaults_raw.get("edited_fraction_column") or "edited_fraction")

    targets_raw = data.get("targets") or []
    if not isinstance(targets_raw, list) or not targets_raw:
        raise ValueError("Efficiency panel config requires a non-empty 'targets' list.")
    base_dir = cfg_path.parent
    targets = [_parse_target(base_dir, defaults_raw, entry) for entry in targets_raw if isinstance(entry, Mapping)]

    return EfficiencyPanelSpec(
        name=name,
        description=description,
        kind=PANEL_KIND,
        cas=cas_spec,
        edited_fraction_column_default=edited_fraction_column_default,
        targets=targets,
    )


def _make_cas_system(spec: CasConfigSpec) -> CasSystem:
    return CasSystem(
        name=spec.name,
        system_type=spec.system_type,
        pam_rules=[],
        cut_offset=spec.cut_offset,
        max_mismatches=spec.max_mismatches,
        weight_mismatch_penalty=spec.weight_mismatch_penalty,
        weight_pam_penalty=spec.weight_pam_penalty,
    )


def _load_edited_fraction_from_file(path: Path, column: str) -> Optional[float]:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or column not in reader.fieldnames:
            return None
        values: List[float] = []
        for row in reader:
            raw = row.get(column)
            if raw is None or str(raw).strip() == "":
                continue
            try:
                values.append(float(raw))
            except ValueError:
                continue
        if not values:
            return None
        return sum(values) / len(values)


def _pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) != len(y) or len(x) < 2:
        return None
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    num = 0.0
    den_x = 0.0
    den_y = 0.0
    for xi, yi in zip(x, y):
        dx = xi - mean_x
        dy = yi - mean_y
        num += dx * dy
        den_x += dx * dx
        den_y += dy * dy
    if den_x <= 0.0 or den_y <= 0.0:
        return None
    return num / math.sqrt(den_x * den_y)


def run_panel(panel: EfficiencyPanelSpec) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    preds: List[float] = []
    obs: List[float] = []
    cas_system = _make_cas_system(panel.cas)
    efficiency_targets = [
        EngineEfficiencyTargetRequest(
            target_id=target.id,
            reference_sequence=target.reference_sequence,
            guide=GuideRNA(sequence=target.guide_sequence, name=target.guide_id),
        )
        for target in panel.targets
    ]
    prediction_lookup = {
        prediction.target_id: prediction.predicted_score
        for prediction in predict_efficiency_for_targets(cas_system, efficiency_targets)
    }

    for target in panel.targets:
        predicted = float(prediction_lookup.get(target.id, 0.0))
        meas = target.measurements

        if meas.edited_fraction_value is not None:
            observed = meas.edited_fraction_value
        elif meas.edited_fraction_file is not None:
            observed = _load_edited_fraction_from_file(
                meas.edited_fraction_file,
                meas.edited_fraction_column,
            )
        else:
            observed = None

        entry: Dict[str, Any] = {
            "target_id": target.id,
            "guide_id": target.guide_id,
            "predicted_efficiency": predicted,
            "observed_edited_fraction": observed,
            "chromatin_features": target.measurements.chromatin_features or {},
        }
        if observed is None:
            entry["status"] = "missing_observation"
        else:
            entry["status"] = "ok"
            preds.append(predicted)
            obs.append(observed)
            entry["abs_error"] = abs(predicted - observed)

        results.append(entry)

    summary: Dict[str, Any] = {
        "panel_name": panel.name,
        "panel_description": panel.description,
        "targets": len(panel.targets),
        "evaluated": len(preds),
    }
    if preds and obs:
        abs_errors = [abs(p - o) for p, o in zip(preds, obs)]
        summary["mae"] = sum(abs_errors) / len(abs_errors)
        r = _pearson(preds, obs)
        if r is not None:
            summary["pearson_r"] = r

    payload: Dict[str, Any] = {
        "schema": {"kind": "crispr_efficiency_panel", "spec_version": "1.0"},
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
        description="Compare Helix CRISPR efficiency scores to edited fractions."
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Efficiency panel YAML config (see templates/crispr_efficiency_panel.helix.yml).",
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
    print(f"[OK] CRISPR efficiency panel results written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
