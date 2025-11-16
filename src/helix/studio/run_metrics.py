from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Mapping, Sequence

Snapshot = Mapping[str, Any]


@dataclass(frozen=True)
class RunMetrics:
    summary: Mapping[str, str]
    intended_count: int
    off_target_count: int
    neutral_count: int
    intended_mass: float
    off_target_mass: float
    neutral_mass: float
    dag_stats: Mapping[str, int]
    top_outcomes: Sequence[tuple[str, float]]
    perf: Mapping[str, float]
    config: Mapping[str, Any]


@dataclass(frozen=True)
class RunVerdict:
    label: str
    details: str
    intended_delta: float = 0.0
    off_target_delta: float = 0.0


def ensure_state_payload(snapshot: Snapshot) -> Mapping[str, Any] | None:
    state_payload = snapshot.get("state") or snapshot
    if isinstance(state_payload, str):
        try:
            state_payload = json.loads(state_payload)
        except Exception:
            return None
    if not isinstance(state_payload, Mapping):
        return None
    return state_payload


def summarize_snapshot(snapshot: Snapshot) -> dict[str, str]:
    state_payload = ensure_state_payload(snapshot) or {}
    timestamp = str(snapshot.get("timestamp") or state_payload.get("timestamp") or "")
    run_kind_key = str(state_payload.get("run_kind", "NONE"))
    run_kind_label = run_kind_key.title()
    run_id = str(state_payload.get("run_id", "—"))
    guide = _guide_snippet(state_payload)
    genome = _genome_label(state_payload)
    label = str(state_payload.get("label", ""))
    config = state_payload.get("config") or {}
    run_cfg = config.get("run_config") or {}
    draws = ""
    if isinstance(run_cfg, Mapping):
        draws = str(run_cfg.get("draws") or "")
    outcomes = state_payload.get("outcomes") or []
    intended_count = sum(1 for outcome in outcomes if _is_intended_outcome(outcome))
    profile_key = (run_kind_key, genome, guide)
    profile_label = f"{run_kind_label} | {genome} | {guide}"
    return {
        "timestamp": timestamp or "—",
        "run_kind": run_kind_label,
        "run_kind_key": run_kind_key,
        "run_id": run_id,
        "guide": guide,
        "genome": genome,
        "label": label,
        "draws": draws or "—",
        "outcomes": str(len(outcomes)),
        "has_intended": intended_count > 0,
        "profile_key": profile_key,
        "profile_label": profile_label,
    }


def compute_run_metrics(snapshot: Snapshot) -> RunMetrics:
    summary = summarize_snapshot(snapshot)
    state = ensure_state_payload(snapshot) or {}
    outcomes = state.get("outcomes") or []
    intended = 0
    off_target = 0
    neutral = 0
    mass = {"intended": 0.0, "off_target": 0.0, "neutral": 0.0}
    for outcome in outcomes:
        category = _classify_outcome(outcome)
        if category == "intended":
            intended += 1
        elif category == "off_target":
            off_target += 1
        else:
            neutral += 1
        try:
            prob = max(0.0, float(outcome.get("probability", 0.0)))
        except Exception:
            prob = 0.0
        mass[category] += prob
    config = state.get("config")
    if not isinstance(config, Mapping):
        config = {}
    return RunMetrics(
        summary=summary,
        intended_count=intended,
        off_target_count=off_target,
        neutral_count=neutral,
        intended_mass=mass["intended"],
        off_target_mass=mass["off_target"],
        neutral_mass=mass["neutral"],
        dag_stats=_estimate_dag_stats(state.get("dag_from_runtime")),
        top_outcomes=_top_outcomes(outcomes),
        perf=_extract_perf(config),
        config=config,
    )


def format_snapshot_text(metrics: RunMetrics) -> str:
    dag = metrics.dag_stats
    lines = [
        f"Run #{metrics.summary['run_id']} [{metrics.summary['run_kind']}]",
        f"Genome: {metrics.summary['genome']}",
        f"Guide: {metrics.summary['guide']}",
        f"Timestamp: {metrics.summary['timestamp']}",
        f"Draws: {metrics.summary['draws']}",
        f"Outcomes: {metrics.summary['outcomes']} (Intended: {metrics.intended_count}, Off-target: {metrics.off_target_count}, Neutral: {metrics.neutral_count})",
        f"Prob mass intended={metrics.intended_mass:.2f} off-target={metrics.off_target_mass:.2f} neutral={metrics.neutral_mass:.2f}",
        (
            "DAG nodes={nodes} branches={branches} leaves={leaves} depth={depth}".format(
                nodes=dag.get("nodes", 0),
                branches=dag.get("branch_nodes", 0),
                leaves=dag.get("leaves", 0),
                depth=dag.get("depth", 0),
            )
        ),
    ]
    perf_line = _format_perf_line(metrics.perf)
    if perf_line:
        lines.append(perf_line)
    lines.extend(
        [
            "",
            "Config JSON:",
            json.dumps(metrics.config, indent=2),
        ]
    )
    return "\n".join(lines)


def diff_summary_text(lhs: RunMetrics, rhs: RunMetrics) -> str:
    verdict = classify_run_delta(lhs, rhs)
    lines = [
        f"Verdict (Run B vs Run A): {verdict.label} – {verdict.details}",
        f"Δ Intended count: {lhs.intended_count - rhs.intended_count}  Δ Off-target: {lhs.off_target_count - rhs.off_target_count}",
        f"Δ Prob mass intended: {lhs.intended_mass - rhs.intended_mass:.2f}  Δ Off-target: {lhs.off_target_mass - rhs.off_target_mass:.2f}",
    ]
    dag_keys = ("nodes", "branch_nodes", "leaves", "depth")
    dag_deltas = []
    for key in dag_keys:
        dag_deltas.append(
            f"{key}: {lhs.dag_stats.get(key, 0) - rhs.dag_stats.get(key, 0)}"
        )
    lines.append("Δ DAG → " + "  ".join(dag_deltas))
    lines.append(_outcome_diff_text(lhs.top_outcomes, rhs.top_outcomes))
    same_genome = lhs.summary["genome"] == rhs.summary["genome"]
    lines.append(f"Genome match: {'yes' if same_genome else 'no'}")
    perf_delta = _format_perf_delta(lhs.perf, rhs.perf)
    if perf_delta:
        lines.append(perf_delta)
    return "\n".join(lines)


def build_report(snapshot: Snapshot) -> Mapping[str, Any]:
    metrics = compute_run_metrics(snapshot)
    verdict = classify_run_quality(metrics)
    return {
        "summary": metrics.summary,
        "genome": {
            "source": metrics.summary.get("genome_source") or metrics.config.get("genome_source"),
            "uri": metrics.config.get("genome_uri"),
            "hash": metrics.config.get("genome_hash"),
        },
        "run_kind": metrics.summary["run_kind"],
        "label": metrics.summary.get("label"),
        "dag": metrics.dag_stats,
        "outcomes": {
            "total": metrics.intended_count + metrics.off_target_count + metrics.neutral_count,
            "intended": metrics.intended_count,
            "off_target": metrics.off_target_count,
            "neutral": metrics.neutral_count,
        },
        "probability_mass": {
            "intended": metrics.intended_mass,
            "off_target": metrics.off_target_mass,
            "neutral": metrics.neutral_mass,
        },
        "top_outcomes": metrics.top_outcomes,
        "perf": metrics.perf,
        "config": metrics.config,
        "verdict": asdict(verdict),
    }


def render_markdown_report(report: Mapping[str, Any]) -> str:
    summary = report["summary"]
    dag = report["dag"]
    outcomes = report["outcomes"]
    probs = report["probability_mass"]
    verdict = report.get("verdict") or {}
    verdict_label = verdict.get("label")
    verdict_details = verdict.get("details")
    lines = [
        f"# Run {summary['run_id']} ({summary['run_kind']})",
        f"Label: {report.get('label') or '—'}",
        f"Genome: {summary['genome']} (source={report.get('genome', {}).get('source')}, uri={report.get('genome', {}).get('uri')})",
        f"Guide: {summary['guide']}",
        f"Timestamp: {summary['timestamp']}",
        f"Draws: {summary['draws']}",
    ]
    if verdict_label and verdict_details:
        lines.append(f"Verdict: {verdict_label} – {verdict_details}")
    lines.extend(
        [
            "",
            "## Outcome stats",
            f"Total outcomes: {outcomes['total']}",
            f"Intended: {outcomes['intended']}  Off-target: {outcomes['off_target']}  Neutral: {outcomes['neutral']}",
            f"Probability mass intended={probs['intended']:.2f} off-target={probs['off_target']:.2f} neutral={probs['neutral']:.2f}",
            "",
            "## DAG stats",
            (
                "Nodes: {nodes}  Branch nodes: {branches}  Leaves: {leaves}  Depth: {depth}".format(
                    nodes=dag.get("nodes", 0),
                    branches=dag.get("branch_nodes", 0),
                    leaves=dag.get("leaves", 0),
                    depth=dag.get("depth", 0),
                )
            ),
            "",
            "## Top outcomes",
        ]
    )
    for label, prob in report.get("top_outcomes", []):
        lines.append(f"- {label}: {prob:.2f}")
    lines.extend(
        [
            "",
        ]
    )
    perf_line = _format_perf_line(report.get("perf", {}))
    if perf_line:
        lines.append(perf_line)
    lines.extend(
        [
            "## Config",
            "```json",
            json.dumps(report.get("config", {}), indent=2),
            "```",
        ]
    )
    return "\n".join(lines)


def is_intended_outcome(outcome: Mapping[str, Any]) -> bool:
    return _is_intended_outcome(outcome)


def _guide_snippet(state_payload: Mapping[str, Any]) -> str:
    config = state_payload.get("config") or {}
    guide = config.get("guide")
    if isinstance(guide, Mapping):
        label = guide.get("id") or guide.get("name")
        if label:
            return str(label)
        start = guide.get("start")
        end = guide.get("end")
        if start is not None and end is not None:
            return f"{start}-{end}"
    return "—"


def _genome_label(state_payload: Mapping[str, Any]) -> str:
    source = state_payload.get("genome_source") or "inline"
    uri = state_payload.get("genome_uri")
    if source == "file" and uri:
        return Path(str(uri)).name
    if source == "inline":
        return "inline"
    return str(uri or source or "—")


def _extract_perf(config: Mapping[str, Any]) -> dict[str, float]:
    perf: dict[str, float] = {}
    timings = config.get("perf")
    if isinstance(timings, Mapping):
        for key, value in timings.items():
            try:
                perf[str(key)] = float(value)
            except Exception:
                continue
    return perf


def _estimate_dag_stats(payload: object) -> dict[str, int]:
    base = {"nodes": 0, "edges": 0, "leaves": 0, "depth": 0, "branch_nodes": 0, "max_branching": 0}
    if not isinstance(payload, Mapping):
        return base
    nodes = payload.get("nodes")
    edges = payload.get("edges")
    node_count = len(nodes) if isinstance(nodes, list) else int(payload.get("node_count", 0))
    edge_count = len(edges) if isinstance(edges, list) else int(payload.get("edge_count", 0))
    depth = int(payload.get("max_depth", 0))
    leaves = 0
    branch_nodes = 0
    max_branching = 0
    if isinstance(nodes, list) and nodes and isinstance(nodes[0], Mapping):
        for node in nodes:
            children = node.get("children")
            child_count = len(children) if isinstance(children, list) else 0
            if child_count == 0:
                leaves += 1
            elif child_count > 1:
                branch_nodes += 1
            if child_count > max_branching:
                max_branching = child_count
    base.update(
        {
            "nodes": node_count,
            "edges": edge_count,
            "leaves": leaves,
            "depth": depth,
            "branch_nodes": branch_nodes,
            "max_branching": max_branching,
        }
    )
    return base


def _top_outcomes(outcomes: Sequence[Mapping[str, Any]]) -> list[tuple[str, float]]:
    items: list[tuple[str, float]] = []
    for outcome in outcomes:
        label = str(outcome.get("label", "?"))
        try:
            prob = float(outcome.get("probability", 0.0))
        except Exception:
            prob = 0.0
        items.append((label, prob))
    items.sort(key=lambda item: item[1], reverse=True)
    return items[:5]


def _classify_outcome(outcome: Mapping[str, Any]) -> str:
    if _is_intended_outcome(outcome):
        return "intended"
    diff_kind = str((outcome.get("diff") or {}).get("kind", "")).lower()
    label = str(outcome.get("label", "")).lower()
    if "no_cut" in diff_kind or "no edit" in label or diff_kind in {"none", "neutral"}:
        return "neutral"
    return "off_target"


def _outcome_diff_text(a: Sequence[tuple[str, float]], b: Sequence[tuple[str, float]]) -> str:
    left = {label: prob for label, prob in a}
    right = {label: prob for label, prob in b}
    shared = set(left) & set(right)
    only_a = set(left) - shared
    only_b = set(right) - shared
    lines = []
    if shared:
        parts = [f"{label} Δp={left[label] - right[label]:.2f}" for label in shared]
        lines.append("Shared: " + ", ".join(parts))
    if only_a:
        lines.append("Only run A: " + ", ".join(f"{label} ({left[label]:.2f})" for label in only_a))
    if only_b:
        lines.append("Only run B: " + ", ".join(f"{label} ({right[label]:.2f})" for label in only_b))
    return "\n".join(lines) if lines else "Outcome diff: none"


def _is_intended_outcome(outcome: Mapping[str, Any]) -> bool:
    label = str(outcome.get("label", "")).lower()
    tags = {str(tag).lower() for tag in outcome.get("tags", []) if isinstance(tag, str)}
    diff_kind = str((outcome.get("diff") or {}).get("kind", "")).lower()
    return any(token in label or token in tags or token in diff_kind for token in ("intended", "hdr", "scarless"))


def classify_run_quality(metrics: RunMetrics) -> RunVerdict:
    total_mass = metrics.intended_mass + metrics.off_target_mass + metrics.neutral_mass
    intended_frac = metrics.intended_mass / total_mass if total_mass else 0.0
    off_frac = metrics.off_target_mass / total_mass if total_mass else 0.0
    delta_mass = metrics.intended_mass - metrics.off_target_mass
    if intended_frac >= 0.6 and delta_mass > 0.05:
        label = "Intended-leaning"
    elif off_frac >= 0.5 and delta_mass < -0.05:
        label = "Off-target heavy"
    else:
        label = "Mixed outcome"
    details = (
        f"Intended mass {metrics.intended_mass:.2f} vs off-target {metrics.off_target_mass:.2f} "
        f"(Δ={delta_mass:+.2f})"
    )
    return RunVerdict(label=label, details=details, intended_delta=delta_mass, off_target_delta=-delta_mass)


def classify_run_delta(baseline: RunMetrics, contender: RunMetrics) -> RunVerdict:
    intended_delta = contender.intended_mass - baseline.intended_mass
    off_delta = contender.off_target_mass - baseline.off_target_mass
    intended_count_delta = contender.intended_count - baseline.intended_count
    off_count_delta = contender.off_target_count - baseline.off_target_count
    threshold_mass = 0.01
    better_intended = intended_delta > threshold_mass or intended_count_delta > 0
    worse_intended = intended_delta < -threshold_mass or intended_count_delta < 0
    better_off = off_delta < -threshold_mass or off_count_delta < 0
    worse_off = off_delta > threshold_mass or off_count_delta > 0
    if better_intended and better_off:
        label = "Better"
    elif worse_intended and worse_off:
        label = "Worse"
    else:
        label = "Trade-off"
    details = (
        f"Δ intended mass {intended_delta:+.2f} (counts {intended_count_delta:+d}), "
        f"Δ off-target mass {off_delta:+.2f} (counts {off_count_delta:+d})"
    )
    return RunVerdict(label=label, details=details, intended_delta=intended_delta, off_target_delta=off_delta)


def _format_perf_line(perf: Mapping[str, float]) -> str:
    if not perf:
        return ""
    parts: list[str] = []
    sim_ms = perf.get("sim_ms")
    if sim_ms is not None:
        parts.append(f"Sim {sim_ms:.0f} ms")
    viz_ms = perf.get("viz_ms")
    if viz_ms is not None:
        parts.append(f"Viz {viz_ms:.0f} ms")
    if not parts:
        return ""
    return "Perf: " + " | ".join(parts)


def _format_perf_delta(lhs: Mapping[str, float], rhs: Mapping[str, float]) -> str:
    keys = set(lhs) | set(rhs)
    if not keys:
        return ""
    parts: list[str] = []
    for key in sorted(keys):
        parts.append(f"{key}: {(lhs.get(key, 0.0) - rhs.get(key, 0.0)):.0f} ms")
    return "Δ Perf → " + "  ".join(parts)
