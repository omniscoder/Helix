"""Standalone CRISPR DAG micro-verification script."""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, GuideRNA, PAMRule
from helix.crispr.simulator import find_candidate_sites

GUIDE_SEQUENCE = "ACGTACGTACGTACGTACGT"
GENOME_PATH = Path("tests/data/crispr_micro.fna")
MAX_DEPTH = 2
MAX_SITES = 50
MIN_PROB = 1e-4
NO_EDIT_PROB = 0.1
INDEL_WINDOW = 3
PAM_PATTERN = "NGG"
CUT_OFFSET = 3
MAX_MISMATCHES = 3
MISMATCH_WEIGHT = 1.0
PAM_WEIGHT = 2.0
EDGE_EPS = 1e-6
# Probability tolerances reflect tiny numerical drift introduced by DAG deduping.
PROB_TOL = 1e-2
PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class CandidateSite:
    chrom: str
    start: int
    end: int
    strand: int
    sequence: str
    score: float


@dataclass(frozen=True)
class BruteEvent:
    chrom: str
    start: int
    end: int
    replacement: str
    logp_delta: float
    stage: str
    metadata: Dict[str, object]


@dataclass
class BruteNode:
    id: str
    sequences: Dict[str, str]
    log_prob: float
    metadata: Dict[str, object]


@dataclass
class BruteEdge:
    source: str
    target: str
    event: BruteEvent
    rule: str


@dataclass
class BruteDAG:
    root_id: str
    nodes: Dict[str, BruteNode]
    edges: List[BruteEdge]


def read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, List[str]] = {}
    name: str | None = None
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                sequences.setdefault(name, [])
                continue
            if name is None:
                raise ValueError("FASTA content missing header before sequence lines.")
            sequences[name].append(line.strip().upper())
    return {chrom: "".join(chunks) for chrom, chunks in sequences.items()}


def run_helix_dag(tmp_path: Path) -> Dict[str, object]:
    output = tmp_path / "helix_crispr_micro.dag.json"
    env = os.environ.copy()
    src_path = PROJECT_ROOT / "src"
    python_path = env.get("PYTHONPATH")
    env["PYTHONPATH"] = f"{src_path}:{python_path}" if python_path else str(src_path)
    cmd = [
        sys.executable,
        "-m",
        "helix.cli",
        "crispr",
        "dag",
        "--genome",
        str(GENOME_PATH),
        "--guide-sequence",
        GUIDE_SEQUENCE,
        "--max-depth",
        str(MAX_DEPTH),
        "--max-sites",
        str(MAX_SITES),
        "--min-prob",
        str(MIN_PROB),
        "--seed",
        "0",
        "--json",
        str(output),
    ]
    subprocess.run(cmd, cwd=PROJECT_ROOT, env=env, check=True)
    return json.loads(output.read_text())


def sequence_signature(sequences: Mapping[str, str]) -> Tuple[Tuple[str, str], ...]:
    return tuple(sorted((chrom, seq) for chrom, seq in sequences.items()))


def logsumexp(values: Sequence[float]) -> float:
    if not values:
        return float("-inf")
    top = max(values)
    if math.isinf(top):
        return top
    return top + math.log(sum(math.exp(val - top) for val in values))


def normalized_from_logs(log_map: Mapping[Tuple[Tuple[str, str], ...], float]) -> Dict[Tuple[Tuple[str, str], ...], float]:
    log_total = logsumexp(list(log_map.values()))
    if math.isinf(log_total):
        raise AssertionError("Total log weight is infinite.")
    return {key: math.exp(val - log_total) for key, val in log_map.items()}


def collect_terminal_logprobs_from_payload(payload: Mapping[str, object]) -> Dict[Tuple[Tuple[str, str], ...], float]:
    nodes: Mapping[str, Mapping[str, object]] = payload["nodes"]  # type: ignore[assignment]
    edges: Sequence[Mapping[str, object]] = payload["edges"]  # type: ignore[assignment]
    sources = {edge["source"] for edge in edges}
    grouped: Dict[Tuple[Tuple[str, str], ...], List[float]] = {}
    for node_id, node in nodes.items():
        if node_id in sources:
            continue
        signature = sequence_signature(node["sequences"])  # type: ignore[arg-type]
        grouped.setdefault(signature, []).append(float(node["log_prob"]))
    return {sig: logsumexp(vals) for sig, vals in grouped.items()}


def collect_terminal_logprobs_from_brute(dag: BruteDAG) -> Dict[Tuple[Tuple[str, str], ...], float]:
    sources = {edge.source for edge in dag.edges}
    grouped: Dict[Tuple[Tuple[str, str], ...], List[float]] = {}
    for node_id, node in dag.nodes.items():
        if node_id in sources:
            continue
        signature = sequence_signature(node.sequences)
        grouped.setdefault(signature, []).append(node.log_prob)
    return {sig: logsumexp(vals) for sig, vals in grouped.items()}


def assert_structural_invariants(payload: Mapping[str, object]) -> None:
    nodes: Mapping[str, Mapping[str, object]] = payload["nodes"]  # type: ignore[assignment]
    edges: Sequence[Mapping[str, object]] = payload["edges"]  # type: ignore[assignment]
    indegree = {node_id: 0 for node_id in nodes}
    for edge in edges:
        indegree[edge["target"]] += 1  # type: ignore[index]

    queue = [node_id for node_id, deg in indegree.items() if deg == 0]
    visited = 0
    while queue:
        node_id = queue.pop()
        visited += 1
        for edge in edges:
            if edge["source"] != node_id:  # type: ignore[index]
                continue
            indegree[edge["target"]] -= 1  # type: ignore[index]
            if indegree[edge["target"]] == 0:  # type: ignore[index]
                queue.append(edge["target"])  # type: ignore[index]
    assert visited == len(nodes), "Graph contains a cycle."

    for edge in edges:
        parent = nodes[edge["source"]]  # type: ignore[index]
        child = nodes[edge["target"]]  # type: ignore[index]
        parent_depth = parent["metadata"].get("time_step", 0)  # type: ignore[index]
        child_depth = child["metadata"].get("time_step", 0)  # type: ignore[index]
        assert child_depth == parent_depth + 1, "Time-step ordering violated."
        assert child_depth <= MAX_DEPTH, "Node exceeded max depth."


def assert_probability_invariants(payload: Mapping[str, object]) -> None:
    nodes: Mapping[str, Mapping[str, object]] = payload["nodes"]  # type: ignore[assignment]
    edges: Sequence[Mapping[str, object]] = payload["edges"]  # type: ignore[assignment]
    by_source: Dict[str, List[Mapping[str, object]]] = {}
    for edge in edges:
        by_source.setdefault(edge["source"], []).append(edge)  # type: ignore[index]

    for node_id, outgoing in by_source.items():
        parent_log = float(nodes[node_id]["log_prob"])
        weights = [float(nodes[edge["target"]]["log_prob"]) - parent_log for edge in outgoing]
        norm = logsumexp(weights)
        probs = [math.exp(w - norm) for w in weights]
        assert math.isclose(sum(probs), 1.0, rel_tol=1e-8, abs_tol=1e-8)

    terminal_logs = collect_terminal_logprobs_from_payload(payload)
    probs = normalized_from_logs(terminal_logs)
    assert math.isclose(sum(probs.values()), 1.0, rel_tol=1e-8, abs_tol=1e-8)


def score_candidate_sites(genome: Mapping[str, str]) -> List[CandidateSite]:
    cas = CasSystem(
        name="SpCas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern=PAM_PATTERN)],
        cut_offset=CUT_OFFSET,
        max_mismatches=MAX_MISMATCHES,
        weight_mismatch_penalty=MISMATCH_WEIGHT,
        weight_pam_penalty=PAM_WEIGHT,
    )
    guide = GuideRNA(sequence=GUIDE_SEQUENCE)
    genome_view = DigitalGenome(sequences=dict(genome))
    sites = find_candidate_sites(genome_view, cas, guide, max_sites=MAX_SITES)
    return [
        CandidateSite(
            chrom=site.chrom,
            start=site.start,
            end=site.end,
            strand=site.strand,
            sequence=site.sequence,
            score=site.on_target_score or 0.0,
        )
        for site in sites
    ]


def apply_event(sequences: Mapping[str, str], event: BruteEvent) -> Dict[str, str]:
    if event.chrom == "__none__":
        return dict(sequences)
    seq = sequences[event.chrom]
    updated = seq[: event.start] + event.replacement + seq[event.end :]
    new_sequences = dict(sequences)
    new_sequences[event.chrom] = updated
    return new_sequences


def make_node_id(parent_id: str, event: BruteEvent, new_time: int) -> str:
    return "|".join(
        [
            parent_id,
            f"{event.chrom}:{event.start}-{event.end}:{event.replacement}",
            event.stage,
            str(new_time),
        ]
    )


def brute_force_crispr_micro(genome: Mapping[str, str]) -> BruteDAG:
    min_log_prob = math.log(MIN_PROB)
    candidates = score_candidate_sites(genome)
    nodes: Dict[str, BruteNode] = {}
    edges: List[BruteEdge] = []
    root_id = "brute_root"
    nodes[root_id] = BruteNode(
        id=root_id,
        sequences=dict(genome),
        log_prob=0.0,
        metadata={"stage": "root", "time_step": 0},
    )
    frontier = [nodes[root_id]]
    for depth in range(MAX_DEPTH):
        if not frontier:
            break
        future: Dict[str, BruteNode] = {}
        for node in frontier:
            stage = node.metadata.get("stage", "root")
            if stage == "root":
                for event in propose_clean_cut(candidates):
                    emit_event(
                        node,
                        event,
                        "crispr.clean_cut",
                        min_log_prob,
                        nodes,
                        edges,
                        future,
                    )
            if stage == "cut":
                for event in propose_indel(node):
                    emit_event(
                        node,
                        event,
                        "crispr.indel_branch",
                        min_log_prob,
                        nodes,
                        edges,
                        future,
                    )
            no_edit_event = BruteEvent(
                chrom="__none__",
                start=0,
                end=0,
                replacement="",
                logp_delta=-1.0,
                stage=stage,
                metadata={"label": "no_edit_branch", "mechanism": "crispr", "stage": stage},
            )
            emit_event(
                node,
                no_edit_event,
                "crispr.no_edit",
                min_log_prob,
                nodes,
                edges,
                future,
            )
        frontier = list(future.values())
    return BruteDAG(root_id=root_id, nodes=nodes, edges=edges)


def emit_event(
    parent: BruteNode,
    event: BruteEvent,
    rule_name: str,
    min_log_prob: float,
    nodes: Dict[str, BruteNode],
    edges: List[BruteEdge],
    future_map: Dict[str, BruteNode],
) -> None:
    new_log = parent.log_prob + event.logp_delta
    if new_log < min_log_prob:
        return
    new_sequences = apply_event(parent.sequences, event)
    metadata = dict(parent.metadata)
    metadata.update(event.metadata)
    metadata["stage"] = event.stage
    new_time = parent.metadata.get("time_step", 0) + 1
    metadata["time_step"] = new_time
    node_id = make_node_id(parent.id, event, new_time)
    new_node = BruteNode(id=node_id, sequences=new_sequences, log_prob=new_log, metadata=metadata)
    nodes[node_id] = new_node
    edges.append(BruteEdge(source=parent.id, target=node_id, event=event, rule=rule_name))
    future_map[node_id] = new_node


def propose_clean_cut(candidates: Sequence[CandidateSite]) -> List[BruteEvent]:
    events: List[BruteEvent] = []
    if not candidates:
        events.append(
            BruteEvent(
                chrom="__none__",
                start=0,
                end=0,
                replacement="",
                logp_delta=0.0,
                stage="no_target",
                metadata={"label": "no_target", "stage": "no_target"},
            )
        )
        return events
    for site in candidates:
        cut_pos = max(site.start, min(site.end, site.start + CUT_OFFSET))
        metadata = {
            "label": "clean_cut",
            "mechanism": "crispr",
            "stage": "cut",
            "site_chrom": site.chrom,
            "site_start": site.start,
            "site_end": site.end,
            "cut_position": cut_pos,
            "strand": site.strand,
        }
        score = max(site.score, EDGE_EPS)
        events.append(
            BruteEvent(
                chrom=site.chrom,
                start=cut_pos,
                end=cut_pos,
                replacement="",
                logp_delta=math.log(score),
                stage="cut",
                metadata=metadata,
            )
        )
    if NO_EDIT_PROB > 0:
        events.append(
            BruteEvent(
                chrom="__none__",
                start=0,
                end=0,
                replacement="",
                logp_delta=math.log(NO_EDIT_PROB),
                stage="no_edit",
                metadata={"label": "no_edit", "mechanism": "crispr", "stage": "no_edit"},
            )
        )
    return events


def propose_indel(node: BruteNode) -> List[BruteEvent]:
    chrom = node.metadata.get("site_chrom")
    cut_pos = node.metadata.get("cut_position")
    if chrom is None or cut_pos is None:
        return []
    seq = node.sequences[chrom]
    start = max(0, cut_pos - INDEL_WINDOW)
    end = min(len(seq), cut_pos + INDEL_WINDOW)
    original = seq[start:end]
    intended = original[::-1]
    events = [
        BruteEvent(
            chrom=chrom,
            start=start,
            end=end,
            replacement=intended,
            logp_delta=0.0,
            stage="repaired",
            metadata={"label": "intended_edit", "mechanism": "crispr", "stage": "repaired"},
        ),
        BruteEvent(
            chrom=chrom,
            start=cut_pos,
            end=min(cut_pos + 1, len(seq)),
            replacement="",
            logp_delta=-0.5,
            stage="repaired",
            metadata={"label": "indel", "mechanism": "crispr", "stage": "repaired"},
        ),
    ]
    return events


def verify(prob_tol: float = PROB_TOL, *, work_dir: Path | None = None) -> Tuple[float, bool]:
    tmp_dir = work_dir or Path(tempfile.mkdtemp(prefix="helix_micro_verify_"))
    tmp_dir.mkdir(parents=True, exist_ok=True)
    payload = run_helix_dag(tmp_dir)
    genome = read_fasta(GENOME_PATH)
    brute_dag = brute_force_crispr_micro(genome)

    assert_structural_invariants(payload)
    assert_probability_invariants(payload)

    helix_terminal_logs = collect_terminal_logprobs_from_payload(payload)
    brute_terminal_logs = collect_terminal_logprobs_from_brute(brute_dag)
    states_match = set(helix_terminal_logs) == set(brute_terminal_logs)
    if not states_match:
        missing_h = set(brute_terminal_logs) - set(helix_terminal_logs)
        missing_b = set(helix_terminal_logs) - set(brute_terminal_logs)
        raise AssertionError(
            f"Terminal states diverged (only in Helix={len(missing_b)}, only in brute={len(missing_h)})."
        )

    helix_probs = normalized_from_logs(helix_terminal_logs)
    brute_probs = normalized_from_logs(brute_terminal_logs)
    max_delta = 0.0
    for signature, helix_prob in helix_probs.items():
        brute_prob = brute_probs[signature]
        max_delta = max(max_delta, abs(helix_prob - brute_prob))
        if not math.isclose(helix_prob, brute_prob, rel_tol=prob_tol, abs_tol=prob_tol):
            raise AssertionError(
                f"Probability mismatch for signature {signature}: helix={helix_prob}, brute={brute_prob}"
            )
    return max_delta, True


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify CRISPR DAG micro-universe via brute force.")
    parser.add_argument(
        "--prob-tol",
        type=float,
        default=PROB_TOL,
        help=f"Probability tolerance for comparisons (default: {PROB_TOL}).",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        help="Optional directory for intermediate DAG artifacts.",
    )
    args = parser.parse_args(argv)
    try:
        max_delta, _ = verify(prob_tol=args.prob_tol, work_dir=args.work_dir)
    except AssertionError as exc:
        print(f"[FAIL] {exc}")
        return 1
    print(f"[OK] CRISPR micro verification passed (max probability delta={max_delta:.6f}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
