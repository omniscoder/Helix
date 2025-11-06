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
    dna_len = len(dna)
    print(f'sequencing dna of length {dna_len}')
    seq_list = []
    kmer_list = dict()
    for i in range(0, (dna_len-filter_size)):
        segment = dna[i:i+filter_size]
        if(seq_list.count(segment) > 0):
            if(kmer_list.get(segment)):
                kmer_list[segment] = kmer_list[segment] + 1
            #once it appears a second time we add it to the kmer list with a count of 2
            else:
                kmer_list[segment] = 2
        #since we havent sequenced this yet we will add it to the list
        else:
            seq_list.append(segment)
            
    #kind of shit until i feel like making a temporary dictionary to put sequences that only appear once
    for key in kmer_list:
        if(kmer_list[key] < 2):
            del kmer_list[key]
            
    return kmer_list

#find kmers with a variability of max_diff SNPs
def find_kmers_with_differences(dna, filter_size: int, max_diff: int):
    dna_len = len(dna)
    print(f'sequencing dna of length {dna_len}')
    seq_list = []
    kmer_list = dict()
    for i in range(0, (dna_len-filter_size)):
        segment = dna[i:i+filter_size]
        if(seq_list.count(segment) > 0):
            if(kmer_list.get(segment)):
                kmer_list[segment] = kmer_list[segment] + 1
            #once it appears a second time we add it to the kmer list with a count of 2
            else:
                kmer_list[segment] = 2
        #since we havent sequenced this yet we will add it to the list
        else:
            seq_list.append(segment)
            
    #kind of shit until i feel like making a temporary dictionary to put sequences that only appear once
    for key in kmer_list:
        if(kmer_list[key] < 2):
            del kmer_list[key]
            
    return kmer_list

    

def compliment(kmer):
    compliment = ""
    for i in range(0, len(kmer)):
        compliment += base_pairs[kmer]
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
