from pathlib import Path

from helix.string import fm as string_fm
from helix.string import edit as string_edit


def test_fm_search_basic():
    text = "GATTACA"
    fm_index = string_fm.build_fm(text)
    hits = string_fm.search(fm_index, "TA")
    assert hits["count"] == 1
    assert hits["positions"] == [3]


def test_fm_search_missing():
    fm_index = string_fm.build_fm("ACGTACGT")
    hits = string_fm.search(fm_index, "TTT")
    assert hits["count"] == 0
    assert hits["positions"] == []


def test_myers_distance():
    assert string_edit.myers("GATTACA", "GACTATA") == 2


def test_myers_search_hits():
    matches = string_edit.myers_search("GATT", "GACTTTGATT", max_distance=1)
    assert any(match["start"] == 5 and match["score"] <= 1 for match in matches)
