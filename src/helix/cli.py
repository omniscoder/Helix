"""Unified Helix CLI that wraps DNA, peptide, RNA, protein, and triage helpers."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from . import bioinformatics, cyclospectrum, nussinov_algorithm, triage
from .io import read_fasta
from .string import fm as string_fm
from .string import edit as string_edit
from .seed import minimizers as seed_minimizers
from .seed import syncmers as seed_syncmers
from .seed.extend import SeedMatch, extend_alignment
from .viz import seed as viz_seed
from .graphs import (
    build_dbg as graph_build_dbg,
    clean_dbg as graph_clean_dbg,
    serialize_graph as graph_serialize,
    deserialize_graph as graph_deserialize,
    export_graphml as graph_export_graphml,
    build_colored_dbg,
)

try:
    from . import protein as protein_module

    PROTEIN_AVAILABLE = getattr(protein_module, "BIOPYTHON_AVAILABLE", True)
except ImportError:  # pragma: no cover - protein extras optional
    protein_module = None
    PROTEIN_AVAILABLE = False


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_sequence_arg(sequence: str | None, path: Path | None, *, default: str | None = None) -> str:
    if sequence and path:
        raise ValueError("Provide either an inline sequence or --input path, not both.")
    if path:
        return _read_text(path)
    if sequence:
        return sequence
    if default is not None:
        return default
    raise ValueError("Missing sequence data; provide a positional sequence or use --input.")


def _parse_spectrum(text: str | None) -> List[int]:
    if not text:
        return []
    tokens = text.replace(",", " ").split()
    if not tokens:
        return []
    return [int(token) for token in tokens]


def command_dna(args: argparse.Namespace) -> None:
    raw = _load_sequence_arg(args.sequence, args.input, default=bioinformatics.seq)
    genome = bioinformatics.normalize_sequence(raw)
    print(f"Sequence length: {len(genome)} nt")
    print(f"GC content: {bioinformatics.gc_content(genome) * 100:.2f}%")

    if args.window > 0 and len(genome) >= args.window:
        windows = bioinformatics.windowed_gc_content(genome, args.window, args.step)
        if windows:
            richest = max(windows, key=lambda win: win.gc_fraction)
            poorest = min(windows, key=lambda win: win.gc_fraction)
            print(
                f"GC window extremes ({args.window} nt): "
                f"max={richest.gc_fraction*100:.2f}% [{richest.start}-{richest.end}), "
                f"min={poorest.gc_fraction*100:.2f}% [{poorest.start}-{poorest.end})"
            )
    else:
        print("GC window summary skipped (window disabled or longer than the sequence).")

    clusters = bioinformatics.find_kmers_with_differences(genome, args.k, args.max_diff)
    sorted_clusters = sorted(clusters.items(), key=lambda item: item[1]["count"], reverse=True)
    if not sorted_clusters:
        print("No k-mer clusters detected with the current parameters.")
    else:
        print(f"\nTop {min(args.top, len(sorted_clusters))} clusters (k={args.k}, max_diff={args.max_diff}):")
        for canonical, info in sorted_clusters[: args.top]:
            patterns = ",".join(info["patterns"])
            positions = ",".join(map(str, info["positions"]))
            print(f"{canonical}\tcount={info['count']}\tpatterns=[{patterns}]\tpositions=[{positions}]")


def command_spectrum(args: argparse.Namespace) -> None:
    spectrum = _parse_spectrum(args.spectrum)
    if args.spectrum_file:
        spectrum.extend(_parse_spectrum(_read_text(args.spectrum_file)))

    if args.peptide:
        theoretical = cyclospectrum.theoretical_spectrum(args.peptide, cyclic=not args.linear)
        mode = "cyclic" if not args.linear else "linear"
        print(f"{mode.title()} spectrum for {args.peptide}:")
        print(" ".join(str(mass) for mass in theoretical))
        if spectrum:
            score = cyclospectrum.score_peptide(args.peptide, spectrum, cyclic=not args.linear)
            print(f"Score vs provided spectrum: {score}")

    if spectrum:
        hits = cyclospectrum.leaderboard_cyclopeptide_sequencing(
            spectrum,
            leaderboard_size=args.leaderboard,
        )
        if hits:
            print(f"\nLeaderboard candidates (top {len(hits)}):")
            for peptide, score in hits:
                print(f"{peptide}\tscore={score}")
        else:
            print("No leaderboard candidates matched the spectrum.")
    elif not args.peptide:
        raise SystemExit("Provide at least --peptide or --spectrum/--spectrum-file.")


def command_rna_fold(args: argparse.Namespace) -> None:
    raw = _load_sequence_arg(args.sequence, args.input, default="GGGAAACCC")
    result = nussinov_algorithm.nussinov(
        raw,
        min_loop_length=args.min_loop,
        allow_wobble_pairs=not args.no_wobble,
    )
    print(f"Sequence ({len(result.sequence)} nt): {result.sequence}")
    print(f"Structure score: {result.score()} pairs")
    print(f"Dot-bracket:\n{result.structure}")
    if args.dot_output:
        args.dot_output.write_text(result.structure, encoding="utf-8")
        print(f"Dot-bracket saved to {args.dot_output}")

    if result.pairs:
        print("\nBase pairs (0-indexed):")
        for i, j in result.pairs:
            print(f"{i:>4} - {j:<4} ({result.sequence[i]}-{result.sequence[j]})")


def command_rna_mea(args: argparse.Namespace) -> None:
    from .rna import mccaskill
    from .viz import rna as viz_rna

    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    header, seq = records[0]
    partition_result = mccaskill.partition(seq)
    posterior = partition_result["posterior"]
    mea_result = mccaskill.mea(seq, posterior, gamma=args.gamma)

    output = {
        "sequence_id": header,
        "partition_function": partition_result["partition_function"],
        "free_energy": partition_result["free_energy"],
        "structure": mea_result["structure"],
        "score": mea_result["score"],
        "gamma": args.gamma,
    }
    text = json.dumps(output, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")
    if args.plot:
        viz_rna.plot_dotplot(posterior, Path(args.plot), title=f"{header} dot-plot")


def command_protein(args: argparse.Namespace) -> None:
    if not PROTEIN_AVAILABLE:
        raise SystemExit("Biopython is required for protein helpers. Install it via 'pip install biopython'.")
    raw = _load_sequence_arg(args.sequence, args.input)
    summary = protein_module.summarize_sequence(raw)
    print(f"Length: {summary.length}")
    print(f"Molecular weight: {summary.molecular_weight:.2f} Da")
    print(f"GRAVY: {summary.gravy:.3f}")
    print(f"Aromaticity: {summary.aromaticity:.3f}")
    print(f"Instability index: {summary.instability_index:.2f}")
    print(f"Charge @ pH 7.0: {summary.charge_at_pH7:.2f}")

    if summary.length >= args.window:
        windows = protein_module.hydropathy_profile(
            summary.sequence,
            window=args.window,
            step=args.step,
            scale=args.scale,
        )
        sorted_windows = sorted(windows, key=lambda win: win.score, reverse=True)
        print(f"\nTop {min(args.top, len(sorted_windows))} hydrophobic windows:")
        for window in sorted_windows[: args.top]:
            print(f"{window.start:>4}-{window.end:<4}\tscore={window.score:.3f}")
    else:
        print("Hydropathy profile skipped (sequence shorter than requested window).")


def _orf_to_dict(orf) -> dict:
    return {
        "start": orf.start,
        "end": orf.end,
        "strand": orf.strand,
        "frame": orf.frame,
        "length_nt": orf.length_nt(),
        "length_aa": orf.length_aa(),
        "peptide": orf.peptide,
    }


def _triage_to_dict(report: triage.TriageReport) -> dict:
    return {
        "sequence": report.sequence,
        "skew": report.skew,
        "clusters": [
            {
                "canonical": cluster.canonical,
                "count": cluster.count,
                "patterns": list(cluster.patterns),
                "positions": list(cluster.positions),
            }
            for cluster in report.clusters
        ],
        "orfs": [_orf_to_dict(orf) for orf in report.orfs],
    }


def command_triage(args: argparse.Namespace) -> None:
    raw = _load_sequence_arg(args.sequence, args.input)
    report = triage.compute_triage_report(
        raw,
        k=args.k,
        max_diff=args.max_diff,
        min_orf_length=args.min_orf_length,
    )
    print(f"Sequence length: {len(report.sequence)} nt")
    print(f"Skew span: min={min(report.skew)} max={max(report.skew)}")
    print(f"Detected {len(report.clusters)} k-mer clusters and {len(report.orfs)} ORFs >= {args.min_orf_length} nt.")
    if report.clusters:
        print("\nTop clusters:")
        for cluster in report.clusters[: args.top]:
            patterns = ",".join(cluster.patterns)
            print(f"{cluster.canonical}\tcount={cluster.count}\tpatterns=[{patterns}]")
    if report.orfs:
        print("\nTop ORFs:")
        for orf in report.orfs[: args.top]:
            print(
                f"start={orf.start} end={orf.end} strand={orf.strand} frame={orf.frame} "
                f"length_nt={orf.length_nt()} length_aa={orf.length_aa()}"
            )

    if args.json:
        payload = _triage_to_dict(report)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON report saved to {args.json}")


def command_string_search(args: argparse.Namespace) -> None:
    records = read_fasta(args.input)
    if not records:
        raise SystemExit(f"No sequences found in {args.input}")

    pattern = args.pattern.upper()
    results = []
    for idx, (header, sequence) in enumerate(records):
        label = header or f"seq_{idx}"
        seq = sequence.upper()
        if args.k == 0:
            fm_index = string_fm.build_fm(seq)
            hits = string_fm.search(fm_index, pattern)
            payload = {
                "sequence_id": label,
                "mode": "exact",
                "hits": hits,
            }
        else:
            matches = string_edit.myers_search(pattern, seq, args.k)
            payload = {
                "sequence_id": label,
                "mode": f"myers_k_{args.k}",
                "pattern": pattern,
                "matches": matches,
            }
        results.append(payload)

    output = {
        "meta": {
            "pattern": pattern,
            "k": args.k,
            "sequence_count": len(records),
        },
        "results": results,
    }

    text = json.dumps(output, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_seed_index(args: argparse.Namespace) -> None:
    records = read_fasta(args.input)
    if not records:
        raise SystemExit(f"No sequences found in {args.input}")

    all_results = []
    for header, seq in records:
        if args.method == "minimizer":
            seeds = seed_minimizers(seq, args.k, args.window)
        else:
            seeds = seed_syncmers(seq, args.k, args.sync)
        payload = {
            "sequence_id": header or "seq",
            "length": len(seq),
            "method": args.method,
            "seed_count": len(seeds),
            "seeds": [{"pos": pos, "kmer": kmer, "hash": h} for pos, kmer, h in seeds],
        }
        all_results.append(payload)
        if args.plot:
            if len(records) == 1:
                output = args.plot
            else:
                output = args.plot.with_name(f"{args.plot.stem}_{len(all_results)}{args.plot.suffix}")
            viz_seed.plot_density(seeds, len(seq), output, title=f"{payload['sequence_id']} ({args.method})")

    data = {"meta": {"method": args.method, "k": args.k}, "results": all_results}
    text = json.dumps(data, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_seed_map(args: argparse.Namespace) -> None:
    ref_records = read_fasta(args.ref)
    if not ref_records:
        raise SystemExit(f"No reference sequences in {args.ref}")
    ref_header, ref_seq = ref_records[0]
    ref_seeds = seed_minimizers(ref_seq, args.k, args.window)
    index = {}
    for pos, _, h in ref_seeds:
        index.setdefault(h, []).append(pos)

    read_records = read_fasta(args.reads)
    results = []
    for header, seq in read_records:
        read_seeds = seed_minimizers(seq, args.k, args.window)
        matches = []
        for pos, _, h in read_seeds:
            if h not in index:
                continue
            for ref_pos in index[h][: args.max_matches]:
                seed = SeedMatch(ref_pos=ref_pos, read_pos=pos, length=args.k)
                aln = extend_alignment(seed, ref_seq, seq, band=args.band, xdrop=args.xdrop)
                matches.append({"seed_ref": ref_pos, "seed_read": pos, "alignment": aln})
        results.append(
            {
                "read_id": header or "read",
                "seed_hits": len(matches),
                "alignments": matches,
            }
        )

    payload = {
        "meta": {
            "reference": ref_header,
            "k": args.k,
            "window": args.window,
            "band": args.band,
            "xdrop": args.xdrop,
        },
        "results": results,
    }
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def _load_reads_from_paths(paths: List[Path]) -> List[str]:
    reads: List[str] = []
    for path in paths:
        records = read_fasta(path)
        for _, seq in records:
            reads.append(seq)
    return reads


def command_dbg_build(args: argparse.Namespace) -> None:
    reads = _load_reads_from_paths(args.reads)
    graph = graph_build_dbg(reads, args.k)
    payload = graph_serialize(graph)
    args.graph.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Graph saved to {args.graph} (nodes={len(graph.nodes)})")
    if args.graphml:
        args.graphml.write_text(graph_export_graphml(graph), encoding="utf-8")
        print(f"GraphML saved to {args.graphml}")


def command_dbg_clean(args: argparse.Namespace) -> None:
    graph_json = json.loads(args.graph.read_text(encoding="utf-8"))
    graph = graph_deserialize(graph_json)
    graph_clean_dbg(graph, tips=not args.no_tips, bubbles=not args.no_bubbles, tip_length=args.tip_length)
    payload = graph_serialize(graph)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Cleaned graph written to {args.out}")


def command_dbg_color(args: argparse.Namespace) -> None:
    if args.labels and len(args.labels) != len(args.reads):
        raise SystemExit("Number of labels must match number of read files.")
    labels = args.labels or [path.stem for path in args.reads]
    reads_by_sample = {}
    for label, path in zip(labels, args.reads):
        reads_by_sample[label] = [seq for _, seq in read_fasta(path)]
    colored = build_colored_dbg(reads_by_sample, args.k)
    presence = {node: sorted(samples) for node, samples in colored.presence.items() if samples}
    payload = {
        "k": colored.graph.k,
        "graph": graph_serialize(colored.graph),
        "samples": colored.samples,
        "presence": presence,
    }
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Colored graph written to {args.out}")


def command_workflows(args: argparse.Namespace) -> None:
    from helix_workflows import run_workflow_config

    results = run_workflow_config(
        args.config,
        output_dir=args.output_dir,
        selected=args.name,
    )
    for result in results:
        print(f"Workflow '{result.name}' completed. Logs at {result.output_dir}")
        for step in result.steps:
            print(f"  - {step.command} -> {step.output_path or 'stdout captured'}")


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(f"matplotlib is required for visualization commands ({exc}).")
    return plt


def command_viz_triage(args: argparse.Namespace) -> None:
    plt = _require_matplotlib()
    data = json.loads(args.json.read_text(encoding="utf-8"))
    skew = data.get("skew", [])
    clusters = data.get("clusters", [])
    orfs = data.get("orfs", [])

    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=False)
    axes[0].plot(skew)
    axes[0].set_title("GC Skew")
    axes[0].set_xlabel("Position")
    axes[0].set_ylabel("Cumulative skew")

    subset_orfs = orfs[: args.top]
    if subset_orfs:
        y_pos = range(len(subset_orfs))
        lengths = [entry["length_nt"] for entry in subset_orfs]
        labels = [f"{entry['strand']}:{entry['frame']}" for entry in subset_orfs]
        axes[1].barh(list(y_pos), lengths)
        axes[1].set_yticks(list(y_pos))
        axes[1].set_yticklabels(labels)
        axes[1].set_xlabel("Length (nt)")
        axes[1].set_title("Top ORFs")
    else:
        axes[1].text(0.5, 0.5, "No ORFs", ha="center", va="center")

    subset_clusters = clusters[: args.top]
    if subset_clusters:
        axes[2].bar([c["canonical"] for c in subset_clusters], [c["count"] for c in subset_clusters])
        axes[2].tick_params(axis="x", rotation=45)
        axes[2].set_title("Top k-mer clusters")
        axes[2].set_ylabel("Count")
    else:
        axes[2].text(0.5, 0.5, "No clusters", ha="center", va="center")

    fig.tight_layout()
    fig.savefig(args.output)
    if args.show:  # pragma: no cover - interactive path
        plt.show()
    plt.close(fig)
    print(f"Triage visualization saved to {args.output}")


def command_viz_hydropathy(args: argparse.Namespace) -> None:
    if not PROTEIN_AVAILABLE:
        raise SystemExit("Biopython is required for hydropathy visualization (pip install biopython).")
    plt = _require_matplotlib()
    raw = _load_sequence_arg(args.sequence, args.input)
    windows = protein_module.hydropathy_profile(raw, window=args.window, step=args.step, scale=args.scale)
    if not windows:
        raise SystemExit("Sequence is shorter than the requested window size.")

    xs = [window.start for window in windows]
    ys = [window.score for window in windows]
    plt.figure(figsize=(10, 4))
    plt.plot(xs, ys, marker="o")
    plt.axhline(0, color="black", linewidth=0.5)
    plt.title(f"Hydropathy profile (window={args.window}, scale={args.scale})")
    plt.xlabel("Position")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.savefig(args.output)
    if args.show:  # pragma: no cover - interactive path
        plt.show()
    plt.close()
    print(f"Hydropathy chart saved to {args.output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Helix unified CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    dna = subparsers.add_parser("dna", help="Summarize GC, windows, and k-mer hotspots.")
    dna.add_argument("--sequence", help="Inline DNA string (defaults to the bundled sample).")
    dna.add_argument("--input", type=Path, help="Path to a FASTA/text file.")
    dna.add_argument("--window", type=int, default=200, help="GC content window size (default: 200).")
    dna.add_argument("--step", type=int, default=50, help="GC window stride (default: 50).")
    dna.add_argument("--k", type=int, default=5, help="k-mer size (default: 5).")
    dna.add_argument("--max-diff", type=int, default=1, help="Maximum mismatches when clustering (default: 1).")
    dna.add_argument("--top", type=int, default=10, help="Print this many top clusters (default: 10).")
    dna.set_defaults(func=command_dna)

    spectrum = subparsers.add_parser("spectrum", help="Compute theoretical spectra or run leaderboard sequencing.")
    spectrum.add_argument("--peptide", help="Peptide sequence to analyse.")
    spectrum.add_argument("--linear", action="store_true", help="Use the linear spectrum instead of cyclic.")
    spectrum.add_argument("--spectrum", help="Comma/space-separated experimental masses.")
    spectrum.add_argument("--spectrum-file", type=Path, help="File containing experimental masses.")
    spectrum.add_argument("--leaderboard", type=int, default=5, help="Leaderboard size (default: 5).")
    spectrum.set_defaults(func=command_spectrum)

    rna = subparsers.add_parser("rna", help="RNA folding + ensemble helpers.")
    rna_sub = rna.add_subparsers(dest="rna_command", required=True)

    rna_fold = rna_sub.add_parser("fold", help="Fold RNA/DNA using the Nussinov prototype.")
    rna_fold.add_argument("--sequence", help="Inline RNA/DNA string (defaults to a demo hairpin).")
    rna_fold.add_argument("--input", type=Path, help="Path to a FASTA/text file.")
    rna_fold.add_argument("--min-loop", type=int, default=3, help="Minimum loop length (default: 3).")
    rna_fold.add_argument("--no-wobble", action="store_true", help="Disable wobble base pairs.")
    rna_fold.add_argument("--dot-output", type=Path, help="Optional path to write the dot-bracket string.")
    rna_fold.set_defaults(func=command_rna_fold)

    rna_mea = rna_sub.add_parser("mea", help="Partition function + MEA structure (McCaskill).")
    rna_mea.add_argument("--fasta", type=Path, required=True, help="FASTA file containing a single sequence.")
    rna_mea.add_argument("--gamma", type=float, default=1.0, help="MEA gamma parameter (default: 1.0).")
    rna_mea.add_argument("--json", type=Path, help="Optional output JSON path.")
    rna_mea.add_argument("--plot", type=Path, help="Optional dot-plot path (requires matplotlib).")
    rna_mea.set_defaults(func=command_rna_mea)

    protein = subparsers.add_parser("protein", help="Summarize protein sequences (requires Biopython).")
    protein.add_argument("--sequence", help="Inline amino-acid string.")
    protein.add_argument("--input", type=Path, help="FASTA/text file containing a protein sequence.")
    protein.add_argument("--window", type=int, default=9, help="Hydropathy window size (default: 9).")
    protein.add_argument("--step", type=int, default=1, help="Hydropathy step size (default: 1).")
    protein.add_argument("--scale", default="kd", help="Hydropathy scale key (default: kd).")
    protein.add_argument("--top", type=int, default=5, help="How many windows to print (default: 5).")
    protein.set_defaults(func=command_protein)

    triage_cmd = subparsers.add_parser("triage", help="Run the combined GC/k-mer/ORF triage report.")
    triage_cmd.add_argument("--sequence", help="Inline DNA/RNA sequence.")
    triage_cmd.add_argument("--input", type=Path, help="Path to a FASTA/text file.")
    triage_cmd.add_argument("--k", type=int, default=5, help="k-mer length (default: 5).")
    triage_cmd.add_argument("--max-diff", type=int, default=1, help="Allowed mismatches for k-mer clustering (default: 1).")
    triage_cmd.add_argument("--min-orf-length", type=int, default=90, help="Minimum ORF length in nucleotides (default: 90).")
    triage_cmd.add_argument("--top", type=int, default=5, help="Show this many clusters/ORFs (default: 5).")
    triage_cmd.add_argument("--json", type=Path, help="Optional path to write the entire report as JSON.")
    triage_cmd.set_defaults(func=command_triage)

    workflows_cmd = subparsers.add_parser("workflows", help="Run YAML-defined workflows.")
    workflows_cmd.add_argument("--config", type=Path, required=True, help="Path to a workflow YAML file.")
    workflows_cmd.add_argument(
        "--output-dir",
        type=Path,
        default=Path("workflow_runs"),
        help="Directory for workflow logs/output (default: workflow_runs).",
    )
    workflows_cmd.add_argument("--name", help="Optional workflow name to run.")
    workflows_cmd.set_defaults(func=command_workflows)

    viz_cmd = subparsers.add_parser("viz", help="Visualization helpers.")
    viz_subparsers = viz_cmd.add_subparsers(dest="viz_command", required=True)

    viz_triage = viz_subparsers.add_parser("triage", help="Plot a triage JSON payload.")
    viz_triage.add_argument("--json", type=Path, required=True, help="Path to a triage JSON file.")
    viz_triage.add_argument(
        "--output",
        type=Path,
        default=Path("triage_viz.png"),
        help="Output image path (default: triage_viz.png).",
    )
    viz_triage.add_argument("--top", type=int, default=5, help="Top N clusters/ORFs to visualize (default: 5).")
    viz_triage.add_argument("--show", action="store_true", help="Display interactively.")
    viz_triage.set_defaults(func=command_viz_triage)

    viz_hydro = viz_subparsers.add_parser("hydropathy", help="Plot a hydropathy profile for a protein.")
    viz_hydro.add_argument("--sequence", help="Inline amino-acid string.")
    viz_hydro.add_argument("--input", type=Path, help="FASTA/text file containing a protein sequence.")
    viz_hydro.add_argument("--window", type=int, default=9, help="Window size (default: 9).")
    viz_hydro.add_argument("--step", type=int, default=1, help="Step size (default: 1).")
    viz_hydro.add_argument("--scale", default="kd", help="Hydropathy scale (default: kd).")
    viz_hydro.add_argument(
        "--output",
        type=Path,
        default=Path("hydropathy.png"),
        help="Output image path (default: hydropathy.png).",
    )
    viz_hydro.add_argument("--show", action="store_true", help="Display interactively.")
    viz_hydro.set_defaults(func=command_viz_hydropathy)

    string_cmd = subparsers.add_parser("string", help="String / sequence search helpers.")
    string_sub = string_cmd.add_subparsers(dest="string_command", required=True)
    string_search = string_sub.add_parser("search", help="Exact or <=k edit-distance search.")
    string_search.add_argument("input", type=Path, help="FASTA or raw text file containing sequence(s).")
    string_search.add_argument("--pattern", required=True, help="Pattern to search for.")
    string_search.add_argument("--k", type=int, default=0, help="Maximum edit distance (default: 0).")
    string_search.add_argument("--json", type=Path, help="Optional path to write the JSON output.")
    string_search.set_defaults(func=command_string_search)

    seed_cmd = subparsers.add_parser("seed", help="Seed extraction and mapping helpers.")
    seed_sub = seed_cmd.add_subparsers(dest="seed_command", required=True)

    seed_index = seed_sub.add_parser("index", help="Compute minimizers or syncmers for a sequence.")
    seed_index.add_argument("input", type=Path, help="FASTA file.")
    seed_index.add_argument("--method", choices=["minimizer", "syncmer"], default="minimizer")
    seed_index.add_argument("--k", type=int, default=15, help="k-mer length.")
    seed_index.add_argument("--window", type=int, default=10, help="Window size for minimizers.")
    seed_index.add_argument("--sync", type=int, default=5, help="s-mer size for syncmers.")
    seed_index.add_argument("--json", type=Path, help="Optional JSON output file.")
    seed_index.add_argument("--plot", type=Path, help="Optional density plot path (requires matplotlib).")
    seed_index.set_defaults(func=command_seed_index)

    seed_map = seed_sub.add_parser("map", help="Seed-and-extend mapping (toy).")
    seed_map.add_argument("--ref", type=Path, required=True, help="Reference FASTA.")
    seed_map.add_argument("--reads", type=Path, required=True, help="Reads FASTA.")
    seed_map.add_argument("--k", type=int, default=15, help="k-mer length.")
    seed_map.add_argument("--window", type=int, default=10, help="Window size for minimizers.")
    seed_map.add_argument("--band", type=int, default=64, help="Band size for extension.")
    seed_map.add_argument("--xdrop", type=int, default=10, help="X-drop threshold.")
    seed_map.add_argument("--max-matches", type=int, default=3, help="Cap per-seed matches to avoid blowups.")
    seed_map.add_argument("--json", type=Path, help="Optional output path.")
    seed_map.set_defaults(func=command_seed_map)

    dbg_cmd = subparsers.add_parser("dbg", help="De Bruijn graph helpers.")
    dbg_sub = dbg_cmd.add_subparsers(dest="dbg_command", required=True)

    dbg_build = dbg_sub.add_parser("build", help="Build a DBG from reads.")
    dbg_build.add_argument("--reads", type=Path, nargs="+", required=True, help="FASTA/FASTQ files.")
    dbg_build.add_argument("--k", type=int, required=True, help="k-mer size.")
    dbg_build.add_argument("--graph", type=Path, required=True, help="Output JSON graph path.")
    dbg_build.add_argument("--graphml", type=Path, help="Optional GraphML output path.")
    dbg_build.set_defaults(func=command_dbg_build)

    dbg_clean = dbg_sub.add_parser("clean", help="Clean a DBG JSON (tips/bubbles).")
    dbg_clean.add_argument("--graph", type=Path, required=True, help="Input JSON graph path.")
    dbg_clean.add_argument("--out", type=Path, required=True, help="Output JSON path.")
    dbg_clean.add_argument("--tip-length", type=int, default=2, help="Tip length threshold (default: 2).")
    dbg_clean.add_argument("--no-tips", action="store_true", help="Disable tip removal.")
    dbg_clean.add_argument("--no-bubbles", action="store_true", help="Disable bubble removal.")
    dbg_clean.set_defaults(func=command_dbg_clean)

    dbg_color = dbg_sub.add_parser("color", help="Build a colored DBG from labeled read sets.")
    dbg_color.add_argument("--reads", type=Path, nargs="+", required=True, help="FASTA files per sample.")
    dbg_color.add_argument("--labels", nargs="+", help="Optional sample labels (defaults to filename stems).")
    dbg_color.add_argument("--k", type=int, required=True, help="k-mer size.")
    dbg_color.add_argument("--out", type=Path, required=True, help="Output colored graph JSON.")
    dbg_color.set_defaults(func=command_dbg_color)

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        args.func(args)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":  # pragma: no cover - manual invocation path
    main()
