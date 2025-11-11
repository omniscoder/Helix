from helix.rna.partition import partition_posteriors, mea_structure


def test_partition_symmetry():
    result = partition_posteriors("AUGCUA")
    posterior = result["P"]
    n = len(posterior)
    for i in range(n):
        for j in range(n):
            assert abs(posterior[i][j] - posterior[j][i]) < 1e-8
    assert result["Q"] > 0


def test_mea_structure_balanced():
    partition_result = partition_posteriors("GGCCAACC")
    posterior = partition_result["P"]
    mea = mea_structure("GGCCAACC", posterior, gamma=1.0)
    structure = mea["dotbracket"]
    assert structure.count("(") == structure.count(")")
    assert len(structure) == len("GGCCAACC")


def test_mfe_hairpin_forms():
    from helix.rna.mfe import mfe_dotbracket

    seq = "GGGAAAUCCC"
    result = mfe_dotbracket(seq)
    db = result["dotbracket"]
    assert db.count("(") == db.count(")")
    assert "(" in db
