"""Pythonic helpers that mirror Helix CLI functionality."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence

from . import bioinformatics, cyclospectrum, nussinov_algorithm, triage

try:  # pragma: no cover - optional dependency
    import protein as protein_module

    PROTEIN_AVAILABLE = getattr(protein_module, "BIOPYTHON_AVAILABLE", True)
except ImportError:  # pragma: no cover - handled gracefully
    protein_module = None
    PROTEIN_AVAILABLE = False


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_sequence(sequence: str | None, path: str | Path | None, *, default: str | None = None) -> str:
    if sequence and path:
        raise ValueError("Provide either a sequence or input_path, not both.")
    if path:
        return _read_text(Path(path))
    if sequence:
        return sequence
    if default is not None:
        return default
    raise ValueError("Sequence is required.")


def dna_summary(
    sequence: str | None = None,
    *,
    input_path: str | Path | None = None,
    window: int = 200,
    step: int = 50,
    k: int = 5,
    max_diff: int = 1,
) -> dict:
    raw = _load_sequence(sequence, input_path, default=bioinformatics.seq)
    genome = bioinformatics.normalize_sequence(raw)
    windows = bioinformatics.windowed_gc_content(genome, window, step) if window > 0 else []
    clusters = bioinformatics.find_kmers_with_differences(genome, k, max_diff)
    return {
        "sequence": genome,
        "length": len(genome),
        "gc_content": bioinformatics.gc_content(genome),
        "gc_windows": [window.__dict__ for window in windows],
        "kmer_clusters": clusters,
    }


def triage_report(
    sequence: str | None = None,
    *,
    input_path: str | Path | None = None,
    k: int = 5,
    max_diff: int = 1,
    min_orf_length: int = 90,
) -> dict:
    raw = _load_sequence(sequence, input_path, default=bioinformatics.seq)
    report = triage.compute_triage_report(raw, k=k, max_diff=max_diff, min_orf_length=min_orf_length)
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
        "orfs": [
            {
                "start": orf.start,
                "end": orf.end,
                "strand": orf.strand,
                "frame": orf.frame,
                "length_nt": orf.length_nt(),
                "length_aa": orf.length_aa(),
                "peptide": orf.peptide,
            }
            for orf in report.orfs
        ],
    }


def fold_rna(
    sequence: str,
    *,
    min_loop_length: int = 3,
    allow_wobble_pairs: bool = True,
) -> dict:
    result = nussinov_algorithm.nussinov(
        sequence,
        min_loop_length=min_loop_length,
        allow_wobble_pairs=allow_wobble_pairs,
    )
    return {
        "sequence": result.sequence,
        "score": result.score(),
        "pairs": result.pairs,
        "dot_bracket": result.structure,
    }


def spectrum_leaderboard(
    peptide: str | None = None,
    *,
    experimental_spectrum: Sequence[int] | None = None,
    cyclic: bool = True,
    leaderboard_size: int = 5,
) -> dict:
    theoretical = cyclospectrum.theoretical_spectrum(peptide or "", cyclic=cyclic) if peptide else []
    hits: List[tuple[str, int]] = []
    if experimental_spectrum:
        hits = cyclospectrum.leaderboard_cyclopeptide_sequencing(
            experimental_spectrum,
            leaderboard_size=leaderboard_size,
        )
    return {
        "theoretical_spectrum": theoretical,
        "leaderboard_hits": [{"peptide": pep, "score": score} for pep, score in hits],
    }


def protein_summary(
    sequence: str | None = None,
    *,
    input_path: str | Path | None = None,
    window: int = 9,
    step: int = 1,
    scale: str = "kd",
) -> dict:
    if not PROTEIN_AVAILABLE:
        raise ImportError("Biopython is required for protein summaries (pip install biopython).")
    raw = _load_sequence(sequence, input_path)
    summary = protein_module.summarize_sequence(raw)
    windows = (
        protein_module.hydropathy_profile(summary.sequence, window=window, step=step, scale=scale)
        if summary.length >= window
        else []
    )
    return {
        "sequence": summary.sequence,
        "length": summary.length,
        "molecular_weight": summary.molecular_weight,
        "aromaticity": summary.aromaticity,
        "instability_index": summary.instability_index,
        "gravy": summary.gravy,
        "charge_at_pH7": summary.charge_at_pH7,
        "hydropathy_profile": [window.__dict__ for window in windows],
    }


def run_workflow(config_path: str | Path, *, output_dir: str | Path, name: str | None = None) -> list:
    from helix_workflows import run_workflow_config

    return run_workflow_config(Path(config_path), output_dir=Path(output_dir), selected=name)
