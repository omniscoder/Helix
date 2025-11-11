"""Demonstrate the helix_api helpers in a single script."""
from __future__ import annotations

from helix import datasets
import helix_api


def main() -> None:
    print("=== DNA SUMMARY ===")
    dataset_path = datasets.get_path("dna/plasmid_demo.fna")
    dna = helix_api.dna_summary(input_path=dataset_path, k=4, max_diff=1, window=100, step=25)
    print(f"Length: {dna['length']} nt | GC content: {dna['gc_content']*100:.2f}%")
    print(f"Detected {len(dna['kmer_clusters'])} k-mer clusters")

    print("\n=== TRIAGE REPORT ===")
    triage_report = helix_api.triage_report(input_path=dataset_path, k=4, max_diff=1, min_orf_length=60)
    print(f"Skew entries: {len(triage_report['skew'])}")
    print(f"Top cluster: {triage_report['clusters'][0]['canonical'] if triage_report['clusters'] else 'N/A'}")

    print("\n=== RNA FOLD ===")
    fold = helix_api.fold_rna("GGGAAACCC", min_loop_length=0)
    print(f"Dot-bracket: {fold['dot_bracket']}")

    print("\n=== SPECTRUM LEADERBOARD ===")
    hits = helix_api.spectrum_leaderboard(
        peptide="NQEL",
        experimental_spectrum=[0, 113, 114, 128, 227, 242, 242, 355, 356, 370, 371, 484],
    )
    print(f"Hits: {hits['leaderboard_hits']}")

    if helix_api.PROTEIN_AVAILABLE:
        print("\n=== PROTEIN SUMMARY ===")
        protein_path = datasets.get_path("protein/demo_protein.faa")
        protein = helix_api.protein_summary(input_path=protein_path, window=11)
        print(f"Length: {protein['length']} aa | MW: {protein['molecular_weight']:.2f} Da")
        print(f"Hydropathy points: {len(protein['hydropathy_profile'])}")
    else:
        print("\nBiopython not installed; skipping protein summary.")


if __name__ == "__main__":
    main()
