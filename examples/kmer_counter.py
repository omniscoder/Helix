"""Count recurring k-mers with Helix's naive helper."""
from collections import Counter

import bioinformatics


def clean_sequence() -> str:
    return "".join(bioinformatics.seq.upper().split())


def main() -> None:
    genome = clean_sequence()
    k = 5
    counts = bioinformatics.find_kmers(genome, k)

    print(f"Top recurring {k}-mers in the sample fragment:")
    for kmer, freq in Counter(counts).most_common(5):
        print(f"{kmer}: {freq}")


if __name__ == "__main__":
    main()
