"""Count recurring k-mers with Helix's helpers."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from helix import bioinformatics


def clean_sequence(source: Path | None = None) -> str:
    if source is None:
        raw = bioinformatics.seq
    else:
        raw = source.read_text()
    return bioinformatics.normalize_sequence(raw)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Print recurring k-mers from a DNA sequence.")
    parser.add_argument(
        "--input",
        type=Path,
        help="Optional path to a text/FASTA file; defaults to the embedded sample.",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=5,
        help="k-mer length (default: 5).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="How many entries to display (default: 10).",
    )
    parser.add_argument(
        "--max-diff",
        type=int,
        default=0,
        help="Allow up to this many mismatches when grouping k-mers (default: 0).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help="Optional path to write the clusters table as CSV.",
    )
    parser.add_argument(
        "--plot-top",
        type=int,
        default=0,
        help="Plot the counts for the top N clusters (requires matplotlib).",
    )
    return parser.parse_args()


def as_rows(clusters: dict[str, dict]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for canonical, info in clusters.items():
        rows.append(
            {
                "canonical": canonical,
                "count": str(info["count"]),
                "patterns": ",".join(info["patterns"]),
                "positions": ",".join(str(pos) for pos in info["positions"]),
            }
        )
    return rows


def maybe_write_csv(path: Path | None, clusters: dict[str, dict]) -> None:
    if path is None:
        return
    rows = as_rows(clusters)
    if not rows:
        print("No clusters to write.")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["canonical", "count", "patterns", "positions"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {path}")


def maybe_plot(top_n: int, clusters: list[tuple[str, dict]]) -> None:
    if top_n <= 0:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print("matplotlib not available; skipping plot.")
        return

    subset = clusters[:top_n]
    labels = [canonical for canonical, _ in subset]
    counts = [info["count"] for _, info in subset]
    plt.bar(labels, counts)
    plt.title(f"Top {len(subset)} k-mer clusters")
    plt.xlabel("Canonical k-mer")
    plt.ylabel("Occurrences")
    plt.tight_layout()
    plt.show()


def main() -> None:
    args = parse_args()
    genome = clean_sequence(args.input)

    clusters = bioinformatics.find_kmers_with_differences(genome, args.k, args.max_diff)
    print(f"Found {len(clusters)} groups with <= {args.max_diff} mismatches.")
    sorted_clusters = sorted(clusters.items(), key=lambda item: item[1]["count"], reverse=True)
    for canonical, info in sorted_clusters[: args.top]:
        patterns = ",".join(info["patterns"])
        positions = ",".join(str(pos) for pos in info["positions"])
        print(f"{canonical}\tcount={info['count']}\tpatterns=[{patterns}]\tpositions=[{positions}]")

    maybe_write_csv(args.csv, clusters)
    maybe_plot(args.plot_top, sorted_clusters)


if __name__ == "__main__":
    main()
