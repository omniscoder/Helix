"""Utilities for generating combined Helix sequence triage reports."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

from . import bioinformatics
from .codon import Orf, find_orfs


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


def compute_triage_report(
    sequence: str,
    *,
    k: int = 5,
    max_diff: int = 1,
    min_orf_length: int = 90,
) -> TriageReport:
    """Return GC skew, k-mer clusters, and ORF data for a sequence."""
    normalized = bioinformatics.normalize_sequence(sequence)
    display_sequence = normalized.replace("T", "U")
    skew_array = bioinformatics.skew(normalized)
    clusters_dict = bioinformatics.find_kmers_with_differences(normalized, k, max_diff)
    clusters = [
        KmerCluster(
            canonical=name.replace("T", "U"),
            count=info["count"],
            patterns=tuple(pattern.replace("T", "U") for pattern in info["patterns"]),
            positions=tuple(info["positions"]),
        )
        for name, info in clusters_dict.items()
    ]
    clusters.sort(key=lambda cluster: cluster.count, reverse=True)

    orfs = find_orfs(normalized, min_length=min_orf_length)

    return TriageReport(
        sequence=display_sequence,
        skew=skew_array.tolist(),
        clusters=clusters,
        orfs=orfs,
    )
