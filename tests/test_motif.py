from helix.motif import discover_motifs
from helix.io import read_fasta


def test_motif_solvers_recovers_consensus():
    records = read_fasta("tests/fixtures/motifs/simple.fna")
    sequences = [seq for _, seq in records]
    params = {
        "em": {"iterations": 20},
        "steme": {"iterations": 20, "restarts": 3},
        "online": {"learning_rate": 0.3, "passes": 3},
    }
    for solver, kwargs in params.items():
        result = discover_motifs(sequences, width=3, solver=solver, **kwargs)
        consensus = result.consensus()
        assert len(consensus) == 3
        assert consensus.startswith("AA")
