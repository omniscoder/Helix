from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nussinov_algorithm import fold_to_dot_bracket, nussinov


def test_nussinov_hairpin_pairs_correctly():
    result = nussinov("GGGAAACCC", min_loop_length=0, allow_wobble_pairs=False)
    assert result.structure == "(((...)))"
    assert result.pairs == [(0, 8), (1, 7), (2, 6)]
    assert result.score() == 3


def test_wobble_pairs_toggle():
    seq = "GUG"
    strict = nussinov(seq, min_loop_length=0, allow_wobble_pairs=False)
    assert strict.pairs == []
    wobble = nussinov(seq, min_loop_length=0, allow_wobble_pairs=True)
    assert wobble.pairs == [(1, 2)]


def test_fold_to_dot_bracket_normalizes_sequence():
    dot = fold_to_dot_bracket("ATGCAU", min_loop_length=0)
    assert len(dot) == 6
    assert " " not in dot
