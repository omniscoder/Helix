from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bioinformatics import find_kmers_with_differences


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
