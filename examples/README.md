# Helix Examples

Welcome to the growing scrapbook of Helix weekend projects. Each example favors readability over optimization and is meant to be hacked on. Activate your virtualenv, install the base dependencies from the root README, and then try the recipes below.

## GC Skew Explorer
- File: `examples/gc_skew_demo.py`
- What it does: plots the GC skew curve for the sample E. coli fragment bundled in `bioinformatics.py`.
- Try it: `python examples/gc_skew_demo.py` (add `--input custom.fna` to graph your own fragment)
- Remix ideas: swap in your own FASTA file by editing the `load_sequence` call, or annotate the plot with candidate origin-of-replication indices.

## K-mer Frequency Sketch
- File: `examples/kmer_counter.py`
- What it does: uses the naive k-mer finder to list recurring 5-mers in the sample genome fragment.
- Try it: `python examples/kmer_counter.py` (use `-k 7 --top 20` to explore different settings, add `--max-diff 1` for SNP-tolerant clusters, `--csv clusters.csv` to export, or `--plot-top 10` for a quick bar chart)
- Remix ideas: slice the CSV into pandas for further filtering or layer the plot atop GC skew windows to spot hotspots.

## Peptide Mass Lookup
- File: `examples/peptide_mass_lookup.py`
- What it does: sums monoisotopic masses for peptide sequences using `amino_acids.mass_index_table`.
- Try it: `python examples/peptide_mass_lookup.py GAVLIM NQEQ`
- Remix ideas: add support for variable modifications or export the results to CSV for a peptide library.

## Translate a Sequence
- File: `examples/translate_sequence.py`
- What it does: wraps `codon.translate_rna` to convert DNA/RNA into a protein string.
- Try it: `python examples/translate_sequence.py AUGGCCUUU` or `python examples/translate_sequence.py --input sample.fna`
- Remix ideas: combine with `peptide_mass_lookup` to score translated segments, or integrate with notebooks for ORF scanning.

## ORF Scanner
- File: `examples/find_orfs.py`
- What it does: lists open reading frames (forward + reverse) using `codon.find_orfs` and optional frameshift chaining.
- Try it: `python examples/find_orfs.py --min-length 60 --include-partial --detect-frameshifts --input sample.fna --orf-fasta peptides.faa --orf-csv orfs.csv --frameshift-csv shifts.csv`
- Remix ideas: visualize ORF start/stop coordinates alongside GC skew, cross-check frameshift hits with reference annotations, or batch-export peptides for BLAST searches.

## Triage Report CLI
- File: `examples/triage_report.py`
- What it does: outputs a PNG plot combining GC skew, ORFs, and top k-mer clusters and can emit CSV/FASTA artifacts.
- Try it: `python examples/triage_report.py --input sample.fna --output triage.png --clusters-csv clusters.csv --orfs-fasta peptides.faa`
- Remix ideas: schedule the script to drop a daily report, or tweak the plotting code to annotate known genes.

## Share Your Remix
- Drop your notebook, script, or write-up in this directory (or link it in an issue) so other hobbyists can learn from it.
- Format tip: prepend your files with a short description, e.g., `rna_folding_traceback.ipynb`.
- If a prototype grows into something sturdier, open a discussion about porting it into OGN.
