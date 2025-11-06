# Helix Examples

Welcome to the growing scrapbook of Helix weekend projects. Each example favors readability over optimization and is meant to be hacked on. Activate your virtualenv, install the base dependencies from the root README, and then try the recipes below.

## GC Skew Explorer
- File: `examples/gc_skew_demo.py`
- What it does: plots the GC skew curve for the sample E. coli fragment bundled in `bioinformatics.py`.
- Try it: `python examples/gc_skew_demo.py`
- Remix ideas: swap in your own FASTA file by editing the `load_sequence` call, or annotate the plot with candidate origin-of-replication indices.

## K-mer Frequency Sketch
- File: `examples/kmer_counter.py`
- What it does: uses the naive k-mer finder to list recurring 5-mers in the sample genome fragment.
- Try it: `python examples/kmer_counter.py`
- Remix ideas: expose the k-mer length via `argparse` and compare results across several organisms.

## Peptide Mass Lookup
- File: `examples/peptide_mass_lookup.py`
- What it does: sums monoisotopic masses for a peptide sequence using `amino_acids.mass_index_table`.
- Try it: `python examples/peptide_mass_lookup.py GAVLIM`
- Remix ideas: add support for variable modifications or export the results to CSV for a peptide library.

## Share Your Remix
- Drop your notebook, script, or write-up in this directory (or link it in an issue) so other hobbyists can learn from it.
- Format tip: prepend your files with a short description, e.g., `rna_folding_traceback.ipynb`.
- If a prototype grows into something sturdier, open a discussion about porting it into OGN.
