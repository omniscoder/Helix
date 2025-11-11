import json
from pathlib import Path

import pytest

from helix import datasets
import helix_api


def test_dna_summary_clusters():
    result = helix_api.dna_summary(sequence="ACGTACGT", k=2, max_diff=0, window=4, step=2)
    assert result["length"] == 8
    assert result["kmer_clusters"]


def test_triage_report_structure():
    report = helix_api.triage_report(sequence="AUGGCCUUUUAA", k=3, max_diff=0, min_orf_length=9)
    assert report["sequence"] == "AUGGCCUUUUAA"
    assert report["clusters"]
    assert report["orfs"]


def test_fold_rna_returns_pairs():
    fold = helix_api.fold_rna("GGGAAACCC", min_loop_length=0, allow_wobble_pairs=False)
    assert fold["dot_bracket"].count("(") == 3


def test_spectrum_leaderboard_hits():
    hits = helix_api.spectrum_leaderboard(
        peptide="NQEL",
        experimental_spectrum=[0, 113, 114, 128, 227, 242, 242, 355, 356, 370, 371, 484],
    )
    assert hits["leaderboard_hits"]


def test_protein_summary_optional():
    if not helix_api.PROTEIN_AVAILABLE:
        pytest.skip("Biopython not installed")
    protein_path = datasets.get_path("protein/demo_protein.faa")
    summary = helix_api.protein_summary(input_path=protein_path, window=11)
    assert summary["length"] > 0
    assert summary["hydropathy_profile"]
