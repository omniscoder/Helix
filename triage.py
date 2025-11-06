"""Utilities for generating combined Helix sequence triage reports."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import bioinformatics
from codon import Orf, find_orfs


@dataclass(frozen=True)
class KmerCluster:
    canonical: str
    count: int
    patterns: Sequence[str]
    positions: Sequence[int]


@dataclass(frozen=True)
class TriageReport:
    sequence: str
    skew: List[int]
    clusters: List[KmerCluster]
    orfs: List[Orf]


def _normalize_sequence(raw: str) -> str:
    return "".join(raw.upper().split())


def compute_triage_report(
    sequence: str,
    *,
    k: int = 5,
    max_diff: int = 1,
    min_orf_length: int = 90,
) -> TriageReport:
    """Return GC skew, k-mer clusters, and ORF data for a sequence."""
    normalized = _normalize_sequence(sequence)
    skew_array = bioinformatics.skew(normalized)
    clusters_dict = bioinformatics.find_kmers_with_differences(normalized, k, max_diff)
    clusters = [
        KmerCluster(
            canonical=name,
            count=info["count"],
            patterns=tuple(info["patterns"]),
            positions=tuple(info["positions"]),
        )
        for name, info in clusters_dict.items()
    ]
    clusters.sort(key=lambda cluster: cluster.count, reverse=True)

    orfs = find_orfs(normalized, min_length=min_orf_length)

    return TriageReport(
        sequence=normalized,
        skew=skew_array.tolist(),
        clusters=clusters,
        orfs=orfs,
    )
