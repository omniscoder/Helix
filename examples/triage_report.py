"""Generate a combined GC skew + ORF + k-mer hotspot report."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

from helix.codon import detect_frameshifts, frameshifts_to_csv, orfs_to_csv, orfs_to_fasta
from helix.triage import KmerCluster, compute_triage_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Produce a triage plot for a DNA/RNA sequence.")
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Inline sequence (DNA/RNA). If omitted, use --input.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to a text/FASTA file containing a single sequence.",
    )
    parser.add_argument("--k", type=int, default=5, help="k-mer length (default: 5).")
    parser.add_argument(
        "--max-diff",
        type=int,
        default=1,
        help="Maximum mismatches allowed when clustering k-mers (default: 1).",
    )
    parser.add_argument(
        "--min-orf-length",
        type=int,
        default=90,
        help="Minimum ORF length in nucleotides (default: 90).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Top N k-mer clusters to display in the bar plot (default: 10).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("triage_report.png"),
        help="Destination image path (default: triage_report.png).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively in addition to saving.",
    )
    parser.add_argument(
        "--clusters-csv",
        type=Path,
        help="Optional path to write clustered k-mers as CSV.",
    )
    parser.add_argument(
        "--orfs-csv",
        type=Path,
        help="Optional path to write ORF metadata as CSV.",
    )
    parser.add_argument(
        "--orfs-fasta",
        type=Path,
        help="Optional path to write ORF peptides as FASTA.",
    )
    parser.add_argument(
        "--detect-frameshifts",
        action="store_true",
        help="Detect frameshift candidates (overrides when --frameshift-csv is provided).",
    )
    parser.add_argument(
        "--frameshift-csv",
        type=Path,
        help="Optional path to write frameshift candidates as CSV.",
    )
    parser.add_argument(
        "--frameshift-gap",
        type=int,
        default=3,
        help="Maximum nucleotide gap when chaining frameshifts (default: 3).",
    )
    return parser.parse_args()


def load_sequence(args: argparse.Namespace) -> str:
    if args.sequence and args.input:
        raise SystemExit("Provide either an inline sequence or --input, not both.")
    if args.input:
        raw = args.input.read_text(encoding="utf-8")
    elif args.sequence:
        raw = args.sequence
    else:
        raise SystemExit("Provide a sequence (positional) or use --input.")
    return raw


def write_clusters_csv(path: Path, clusters: list[KmerCluster]) -> None:
    rows = [
        {
            "canonical": cluster.canonical,
            "count": cluster.count,
            "patterns": ",".join(cluster.patterns),
            "positions": ",".join(str(pos) for pos in cluster.positions),
        }
        for cluster in clusters
    ]
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def render_plot(report, top_n: int, output: Path, show: bool) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=False)

    axes[0].plot(report.skew)
    axes[0].set_title("GC Skew")
    axes[0].set_xlabel("Position")
    axes[0].set_ylabel("Cumulative skew")

    axes[1].set_title("ORFs")
    axes[1].set_xlabel("Position")
    axes[1].set_ylabel("Strand")
    if report.orfs:
        for orf in report.orfs:
            strand_offset = 0.2 if orf.strand == "+" else -0.2
            axes[1].plot([orf.start, orf.end], [strand_offset, strand_offset], linewidth=4)
        axes[1].set_yticks([-0.2, 0.2])
        axes[1].set_yticklabels(["-", "+"])
    else:
        axes[1].text(0.5, 0.5, "No ORFs found", ha="center", va="center")

    subset = report.clusters[:top_n] if top_n > 0 else report.clusters
    axes[2].set_title(f"Top {len(subset)} k-mer clusters")
    if subset:
        axes[2].bar([c.canonical for c in subset], [c.count for c in subset])
        axes[2].tick_params(axis="x", rotation=45)
    else:
        axes[2].text(0.5, 0.5, "No clusters found", ha="center", va="center")
    axes[2].set_ylabel("Count")

    fig.tight_layout()
    fig.savefig(output)
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    args = parse_args()
    raw_sequence = load_sequence(args)
    report = compute_triage_report(
        raw_sequence,
        k=args.k,
        max_diff=args.max_diff,
        min_orf_length=args.min_orf_length,
    )

    render_plot(report, args.top, args.output, args.show)
    print(f"Triage plot saved to {args.output}")

    if args.clusters_csv:
        write_clusters_csv(args.clusters_csv, report.clusters)
        print(f"Cluster table written to {args.clusters_csv}")

    if args.orfs_csv:
        orf_rows = orfs_to_csv(report.orfs)
        if orf_rows:
            with args.orfs_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=orf_rows[0].keys())
                writer.writeheader()
                writer.writerows(orf_rows)
        else:
            args.orfs_csv.write_text("", encoding="utf-8")
        print(f"ORF table written to {args.orfs_csv}")

    if args.orfs_fasta:
        fasta_text = orfs_to_fasta(report.orfs)
        args.orfs_fasta.write_text(fasta_text, encoding="utf-8")
        print(f"ORF peptides written to {args.orfs_fasta}")

    frameshift_requested = args.detect_frameshifts or args.frameshift_csv
    if frameshift_requested:
        events = detect_frameshifts(
            report.sequence,
            min_orf_length=args.min_orf_length,
            gap_tolerance=args.frameshift_gap,
        )
        if args.frameshift_csv:
            rows = frameshifts_to_csv(events)
            if rows:
                with args.frameshift_csv.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                    writer.writeheader()
                    writer.writerows(rows)
            else:
                args.frameshift_csv.write_text("", encoding="utf-8")
            print(f"Frameshift table written to {args.frameshift_csv}")

        print(f"Frameshift candidates detected: {len(events)}")


if __name__ == "__main__":
    main()
