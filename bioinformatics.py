from collections import Counter, defaultdict
import numpy as np
from ctypes import *
import matplotlib.pyplot as plt
import pandas as pd

#test seq about 500 nucleotides from e-coli genome
seq = """atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"""

#human_dna = pd.read_table("input/dna/human.txt")

base_pairs = {'G': 'C', 'C': 'G', 'T': 'A', 'A': 'T'}

#refactor this to return starting locations in genome as well
def find_kmers(dna, filter_size: int):
    dna = "".join(dna.upper().split())
    dna_len = len(dna)
    print(f'sequencing dna of length {dna_len}')
    counts = Counter(dna[i:i+filter_size] for i in range(0, dna_len - filter_size + 1))
    return {segment: freq for segment, freq in counts.items() if freq >= 2}

def _hamming_distance(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def find_kmers_with_differences(dna, filter_size: int, max_diff: int):
    dna = "".join(dna.upper().split())
    kmer_positions = defaultdict(list)
    for i in range(0, len(dna) - filter_size + 1):
        segment = dna[i:i+filter_size]
        kmer_positions[segment].append(i)

    if max_diff <= 0:
        return {
            kmer: {
                "count": len(positions),
                "positions": positions,
                "patterns": [kmer],
            }
            for kmer, positions in kmer_positions.items()
            if len(positions) >= 2
        }

    kmers = list(kmer_positions.keys())
    parent = {kmer: kmer for kmer in kmers}

    def find_parent(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        root_a, root_b = find_parent(a), find_parent(b)
        if root_a == root_b:
            return
        if root_a < root_b:
            parent[root_b] = root_a
        else:
            parent[root_a] = root_b

    for i, kmer_a in enumerate(kmers):
        for kmer_b in kmers[i+1:]:
            if _hamming_distance(kmer_a, kmer_b) <= max_diff:
                union(kmer_a, kmer_b)

    clusters = defaultdict(lambda: {"patterns": [], "positions": []})
    for kmer, positions in kmer_positions.items():
        root = find_parent(kmer)
        clusters[root]["patterns"].append(kmer)
        clusters[root]["positions"].extend(positions)

    results = {}
    for root, data in clusters.items():
        total = len(data["positions"])
        if total < 2:
            continue
        canonical = min(data["patterns"])
        results[canonical] = {
            "count": total,
            "positions": sorted(data["positions"]),
            "patterns": sorted(data["patterns"]),
        }
    return results

    

def compliment(kmer):
    compliment = ""
    for base in kmer.upper():
        compliment += base_pairs[base]
    return compliment

def reverse_compliment(kmer):
    #since we bind with the reverse 
    reverse = kmer[::-1]
    compliment = ""
    for i in range(0, len(reverse)):
        compliment += base_pairs[reverse[i]]
        
    return compliment


def skew(genome):
    x = np.arange(0, len(genome), 1)
    skewData = np.empty(shape=len(genome) + 1, dtype=int)
    skewData[0] = 0
    print(len(skewData))
    for i in range(1, len(genome) + 1):
        if(genome[i-1] == 'G'):
            skewData[i] = (skewData[i-1] + 1)
        elif(genome[i-1] == 'C'):
            skewData[i] = (skewData[i-1] - 1)
        else:
            skewData[i] = skewData[i-1]
    return skewData
                    



kmer_dict = dict()

def main():
    #most DnaA boxes appear in 3-9mers
    for i in range (3, 9):
        kmer_dict[i] = find_kmers(seq.upper().replace('\n', ''), i)
        
    print(kmer_dict[3])
    skewData = skew(seq.upper())
    plt.plot(skewData)
    plt.xlabel('nucleotide position')
    plt.ylabel('skew')
    plt.show()
    
    
if __name__ == "__main__":
    main()
