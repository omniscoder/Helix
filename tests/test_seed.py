import random

from helix.seed import minimizers, syncmers, extend_alignment
from helix.seed.extend import SeedMatch


def revcomp(seq: str) -> str:
    rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(rc.get(base, base) for base in reversed(seq))


def test_minimizers_reverse_complement():
    seq = "ACGTGTCAGTACGTAGCTAGCTAGGAT"
    seeds = minimizers(seq, 5, 4)
    seeds_rc = minimizers(revcomp(seq), 5, 4)
    hashes = sorted(h for _, _, h in seeds)
    hashes_rc = sorted(h for _, _, h in seeds_rc)
    assert hashes == hashes_rc


def test_syncmers_density_under_mutation():
    random.seed(0)
    seq = "".join(random.choice("ACGT") for _ in range(200))
    mutated = list(seq)
    for i in range(0, len(mutated), 20):
        mutated[i] = random.choice("ACGT".replace(mutated[i], ""))
    mutated = "".join(mutated)
    seeds_orig = syncmers(seq, 7, 3)
    seeds_mut = syncmers(mutated, 7, 3)
    assert len(seeds_mut) > 0
    density_ratio = len(seeds_mut) / max(1, len(seeds_orig))
    assert 0.5 <= density_ratio <= 1.5


def test_extend_alignment_perfect_match():
    ref = "ACGTACGTACGT"
    read = "CGTACG"
    seed = SeedMatch(ref_pos=1, read_pos=0, length=4)
    result = extend_alignment(seed, ref, read, band=10, xdrop=5)
    assert result["score"] > 0
    assert result["matches"] >= len(read) - 1
