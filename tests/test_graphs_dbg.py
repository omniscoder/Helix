from pathlib import Path

from helix.graphs import build_dbg, clean_dbg, serialize_graph, deserialize_graph


def read_fasta(path: Path):
    seqs = []
    header = None
    chunks = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                seqs.append("".join(chunks))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header:
        seqs.append("".join(chunks))
    return seqs


def test_build_and_clean_dbg(tmp_path: Path):
    fixture = Path("tests/fixtures/graphs/bubble_reads.fna")
    reads = read_fasta(fixture)
    graph = build_dbg(reads, k=3)
    assert "AA" in graph.nodes
    cleaned = clean_dbg(graph, tips=True, bubbles=True, tip_length=2)
    serialized = serialize_graph(cleaned)
    restored = deserialize_graph(serialized)
    assert restored.k == 3
    assert "AA" in restored.nodes
