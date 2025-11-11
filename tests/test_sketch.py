from helix.sketch import compute_minhash, mash_distance, compute_hll, union_hll


def test_minhash_consistency_rc():
    seq = "ACGTACGTACGT"
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    sketch1 = compute_minhash(seq, k=3, sketch_size=10)
    sketch2 = compute_minhash(rc, k=3, sketch_size=10)
    assert sketch1.hashes == sketch2.hashes


def test_mash_distance_identical():
    seq = "ACGTACGTACGT"
    sketch = compute_minhash(seq, k=3, sketch_size=10)
    dist = mash_distance(sketch, sketch)
    assert dist == 0.0


def test_hll_union_estimate():
    seq_a = "ACGTACGTACGT"
    seq_b = "ACGTACGTAAAA"
    hll_a = compute_hll(seq_a, k=3, p=8)
    hll_b = compute_hll(seq_b, k=3, p=8)
    union = union_hll(hll_a, hll_b)
    est_a = hll_a.estimate()
    est_b = hll_b.estimate()
    est_union = union.estimate()
    assert est_union >= max(est_a, est_b)
