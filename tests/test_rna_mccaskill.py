import math

from helix.rna import mccaskill


def test_partition_symmetry():
    result = mccaskill.partition("AUGCUA")
    posterior = result["posterior"]
    n = len(posterior)
    for i in range(n):
        for j in range(n):
            assert abs(posterior[i][j] - posterior[j][i]) < 1e-8
    assert result["free_energy"] <= 0


def test_mea_structure_balanced():
    partition_result = mccaskill.partition("GGCCAACC")
    posterior = partition_result["posterior"]
    mea = mccaskill.mea("GGCCAACC", posterior, gamma=1.0)
    structure = mea["structure"]
    assert structure.count("(") == structure.count(")")
    assert len(structure) == len("GGCCAACC")
