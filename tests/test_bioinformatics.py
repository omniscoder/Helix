from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from helix.bioinformatics import (
    complement,
    find_kmers,
    find_kmers_with_differences,
    gc_content,
    normalize_sequence,
    reverse_complement,
    skew,
    windowed_gc_content,
)


def test_find_kmers_with_differences_exact():
    dna = "ACGTACGT"
    result = find_kmers_with_differences(dna, filter_size=3, max_diff=0)
    assert "ACG" in result
    assert result["ACG"]["positions"] == [0, 4]
    assert result["ACG"]["count"] == 2


def test_find_kmers_with_differences_fuzzy():
    dna = "AAAAGAAG"
    result = find_kmers_with_differences(dna, filter_size=4, max_diff=1)
    key = min(result.keys())
    cluster = result[key]
    assert cluster["count"] == 5
    assert cluster["positions"] == [0, 1, 2, 3, 4]
    assert set(cluster["patterns"]) == {"AAAA", "AAAG", "AAGA", "AGAA", "GAAG"}


def test_find_kmers_with_positions():
    dna = "ACGTAC"
    hits = find_kmers(dna, filter_size=2, min_count=2, include_positions=True)
    assert hits["AC"]["count"] == 2
    assert hits["AC"]["positions"] == [0, 4]


def test_normalize_sequence_handles_headers_and_u():
    raw = ">header\naugc\nccnn"
    assert normalize_sequence(raw) == "ATGCCCNN"


def test_gc_metrics_and_windows():
    seq = "GCGTAT"
    assert gc_content(seq) == pytest.approx(3 / 6)
    windows = windowed_gc_content(seq, window=4, step=2)
    assert [round(win.gc_fraction, 2) for win in windows] == [0.75, 0.25]


def test_skew_profile_and_complements():
    seq = "GCGC"
    profile = skew(seq)
    assert profile.tolist() == [0, 1, 0, 1, 0]
    assert complement("ATGC") == "TACG"
    assert reverse_complement("ATGC") == "GCAT"
