"""Helix DNA utilities: GC metrics, k-mer helpers, and quick-look CLI."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np


DNA_BASES = {"A", "C", "G", "T"}
AMBIGUOUS_BASES = {"N"}
COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# Sample ~500 nt fragment from an E. coli genome used by the demo scripts.
seq = """atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"""


@dataclass(frozen=True)
class GCWindow:
    """Container describing the GC content for a genome window."""

    start: int
    end: int
    gc_fraction: float


def normalize_sequence(raw: str, *, allow_ambiguous: bool = True) -> str:
    """Uppercase, strip FASTA headers, drop whitespace, and convert U->T."""
    if not raw:
        return ""

    allowed = set(DNA_BASES)
    if allow_ambiguous:
        allowed |= AMBIGUOUS_BASES

    cleaned: List[str] = []
    for line in raw.splitlines():
        if line.startswith(">"):
            continue
        for char in line.upper():
            if char in {" ", "\t", "\r"}:
                continue
            if char == "U":
                char = "T"
            if char not in allowed:
                raise ValueError(f"Unexpected base '{char}' in sequence.")
            cleaned.append(char)
    return "".join(cleaned)


def gc_content(sequence: str) -> float:
    """Return the GC fraction for a DNA string."""
    dna = normalize_sequence(sequence)
    return _gc_fraction(dna)


def _gc_fraction(dna: str) -> float:
    if not dna:
        return 0.0
    gc = sum(1 for base in dna if base in {"G", "C"})
    return gc / len(dna)


def windowed_gc_content(sequence: str, window: int, step: int = 1) -> List[GCWindow]:
    """Return GC fractions across sliding windows."""
    if window <= 0:
        raise ValueError("window must be > 0")
    if step <= 0:
        raise ValueError("step must be > 0")

    dna = normalize_sequence(sequence)
    if window > len(dna):
        return []

    windows: List[GCWindow] = []
    for start in range(0, len(dna) - window + 1, step):
        end = start + window
        gc_fraction = _gc_fraction(dna[start:end])
        windows.append(GCWindow(start=start, end=end, gc_fraction=gc_fraction))
    return windows


def _hamming_distance(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def find_kmers(dna: str, filter_size: int, *, min_count: int = 2, include_positions: bool = False):
    """Return recurring k-mers (optionally with their positions)."""
    if filter_size <= 0:
        raise ValueError("filter_size must be > 0")
    dna = normalize_sequence(dna)
    if filter_size > len(dna):
        return {}

    positions_map = defaultdict(list)
    for i in range(0, len(dna) - filter_size + 1):
        positions_map[dna[i : i + filter_size]].append(i)

    results: Dict[str, Dict[str, List[int] | int]] = {}
    for segment, positions in positions_map.items():
        freq = len(positions)
        if freq < min_count:
            continue
        if include_positions:
            results[segment] = {"count": freq, "positions": positions.copy()}
        else:
            results[segment] = {"count": freq}
    return results


def find_kmers_with_differences(dna: str, filter_size: int, max_diff: int):
    """Cluster k-mers allowing for mismatches."""
    if filter_size <= 0:
        raise ValueError("filter_size must be > 0")
    if max_diff < 0:
        raise ValueError("max_diff must be >= 0")

    dna = normalize_sequence(dna)
    if filter_size > len(dna):
        return {}

    kmer_positions = defaultdict(list)
    for i in range(0, len(dna) - filter_size + 1):
        segment = dna[i : i + filter_size]
        kmer_positions[segment].append(i)

    if max_diff == 0:
        return {
            kmer: {
                "count": len(positions),
                "positions": positions.copy(),
                "patterns": [kmer],
            }
            for kmer, positions in kmer_positions.items()
            if len(positions) >= 2
        }

    kmers = list(kmer_positions.keys())
    parent = {kmer: kmer for kmer in kmers}

    def find_parent(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        root_a, root_b = find_parent(a), find_parent(b)
        if root_a == root_b:
            return
        if root_a < root_b:
            parent[root_b] = root_a
        else:
            parent[root_a] = root_b

    for i, kmer_a in enumerate(kmers):
        for kmer_b in kmers[i + 1 :]:
            if _hamming_distance(kmer_a, kmer_b) <= max_diff:
                union(kmer_a, kmer_b)

    clusters = defaultdict(lambda: {"patterns": [], "positions": []})
    for kmer, positions in kmer_positions.items():
        root = find_parent(kmer)
        clusters[root]["patterns"].append(kmer)
        clusters[root]["positions"].extend(positions)

    results = {}
    for root, data in clusters.items():
        total = len(data["positions"])
        if total < 2:
            continue
        canonical = min(data["patterns"])
        results[canonical] = {
            "count": total,
            "positions": sorted(data["positions"]),
            "patterns": sorted(set(data["patterns"])),
        }
    return results


def complement(kmer: str) -> str:
    """Return the complement for a DNA sequence."""
    dna = normalize_sequence(kmer)
    return "".join(COMPLEMENT_MAP[base] for base in dna)


def reverse_complement(kmer: str) -> str:
    """Return the reverse complement for a DNA sequence."""
    return complement(kmer)[::-1]


def compliment(kmer: str) -> str:
    """Backward-compatible alias for complement()."""
    return complement(kmer)


def reverse_compliment(kmer: str) -> str:
    """Backward-compatible alias for reverse_complement()."""
    return reverse_complement(kmer)


def skew(genome: str) -> np.ndarray:
    """Return the cumulative GC skew profile."""
    dna = normalize_sequence(genome)
    skew_data = np.zeros(len(dna) + 1, dtype=int)
    for i, base in enumerate(dna, start=1):
        if base == "G":
            skew_data[i] = skew_data[i - 1] + 1
        elif base == "C":
            skew_data[i] = skew_data[i - 1] - 1
        else:
            skew_data[i] = skew_data[i - 1]
    return skew_data


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Quick-look DNA metrics using Helix helpers.")
    parser.add_argument(
        "--input",
        type=Path,
        help="Optional path to a text/FASTA file; defaults to the embedded sample fragment.",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=5,
        help="k-mer size when summarizing hotspots (default: 5).",
    )
    parser.add_argument(
        "--max-diff",
        type=int,
        default=1,
        help="Maximum mismatches allowed when clustering k-mers (default: 1).",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=200,
        help="Window size for GC content summaries (default: 200). Set to 0 to skip.",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=50,
        help="Step size for GC windows (default: 50).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of k-mer clusters to print (default: 10).",
    )
    parser.add_argument(
        "--plot-skew",
        action="store_true",
        help="Plot the GC skew curve (requires matplotlib).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.input:
        raw = args.input.read_text(encoding="utf-8")
    else:
        raw = seq

    genome = normalize_sequence(raw)
    print(f"Sequence length: {len(genome)} nt")
    print(f"GC content: {gc_content(genome) * 100:.2f}%")

    if args.window > 0 and len(genome) >= args.window:
        windows = windowed_gc_content(genome, args.window, args.step)
        if windows:
            richest = max(windows, key=lambda win: win.gc_fraction)
            poorest = min(windows, key=lambda win: win.gc_fraction)
            print(
                f"GC window extremes ({args.window} nt): "
                f"max={richest.gc_fraction*100:.2f}% [{richest.start}-{richest.end}), "
                f"min={poorest.gc_fraction*100:.2f}% [{poorest.start}-{poorest.end})"
            )
    else:
        print("GC window summary skipped (window larger than sequence or disabled).")

    clusters = find_kmers_with_differences(genome, args.k, args.max_diff)
    sorted_clusters = sorted(clusters.items(), key=lambda item: item[1]["count"], reverse=True)
    print(f"\nTop {min(args.top, len(sorted_clusters))} k-mer clusters (k={args.k}, max_diff={args.max_diff}):")
    for canonical, info in sorted_clusters[: args.top]:
        print(
            f"{canonical}\tcount={info['count']}\tpatterns={','.join(info['patterns'])}"
            f"\tpositions={','.join(map(str, info['positions']))}"
        )
    if not sorted_clusters:
        print("No clusters detected with the current settings.")

    if args.plot_skew:
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover - optional dependency
            print(f"matplotlib unavailable ({exc}); skipping plot.")
        else:  # pragma: no cover - rendering not exercised in tests
            profile = skew(genome)
            plt.plot(range(len(profile)), profile)
            plt.title("GC Skew")
            plt.xlabel("Nucleotide position")
            plt.ylabel("Cumulative skew")
            plt.tight_layout()
            plt.show()


if __name__ == "__main__":
    main()
