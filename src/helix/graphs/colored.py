"""Colored De Bruijn graphs."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Set

from .dbg import DeBruijnGraph, build_dbg


@dataclass
class ColoredDBG:
    graph: DeBruijnGraph
    samples: List[str]
    presence: Dict[str, Set[str]] = field(default_factory=dict)


def build_colored_dbg(reads_by_sample: Dict[str, Iterable[str]], k: int) -> ColoredDBG:
    all_reads: List[str] = []
    for sample_reads in reads_by_sample.values():
        all_reads.extend(sample_reads)
    base = build_dbg(all_reads, k)
    presence: Dict[str, Set[str]] = {node_id: set() for node_id in base.nodes}
    for sample, reads in reads_by_sample.items():
        for read in reads:
            seq = read.upper()
            for i in range(len(seq) - k + 1):
                node_id = seq[i : i + k - 1]
                if node_id in presence:
                    presence[node_id].add(sample)
    return ColoredDBG(graph=base, samples=list(reads_by_sample.keys()), presence=presence)


def pseudoalign(read: str, colored_dbg: ColoredDBG) -> Set[str]:
    """Return samples containing k-mers from `read` (intersection heuristic)."""
    k = colored_dbg.graph.k
    seq = read.upper()
    candidate_samples: Set[str] | None = None
    for i in range(len(seq) - k + 1):
        node_id = seq[i : i + k - 1]
        node_samples = colored_dbg.presence.get(node_id)
        if not node_samples:
            continue
        if candidate_samples is None:
            candidate_samples = set(node_samples)
        else:
            candidate_samples &= node_samples
        if candidate_samples == set():
            break
    return candidate_samples or set()
