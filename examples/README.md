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

## Cyclo-spectrum Candidate Finder
- File: `examples/cyclospectrum_demo.py`
- What it does: prints linear/cyclic spectra, scores peptides against an experimental spectrum, and runs the leaderboard-based cyclo-peptide sequencing search.
- Try it: `python examples/cyclospectrum_demo.py --peptide NQEL --spectrum "0,113,114,128,227,242,242,355,356,370,371,484"`
- Remix ideas: benchmark how fast the leaderboard converges for your dataset, or feed the hits into downstream BLAST/proteomics workflows.

## Nussinov Traceback Explorer
- File: `examples/nussinov_trace.py`
- What it does: folds RNA/DNA into dot-bracket notation, prints the base-pair list, and optionally saves the secondary structure string.
- Try it: `python examples/nussinov_trace.py --input hairpin.fasta --min-loop 4`
- Remix ideas: compare wobble vs canonical pairing, overlay the dot-bracket string on a reference alignment, or batch-fold viral genomes for dashboarding.

## Helix API Demo
- File: `examples/helix_api_demo.py`
- What it does: calls `helix_api` to fetch DNA summaries, triage reports, RNA folds, peptide leaderboard hits, and protein metrics without shelling out.
- Try it: `python examples/helix_api_demo.py` (installs optional Biopython if you want the protein section).
- Remix ideas: drop the snippets into notebooks or wrap them inside Streamlit dashboards for quick exploratory reports.

## Share Your Remix
- Drop your notebook, script, or write-up in this directory (or link it in an issue) so other hobbyists can learn from it.
- Format tip: prepend your files with a short description, e.g., `rna_folding_traceback.ipynb`.
- If a prototype grows into something sturdier, open a discussion about porting it into OGN.

## Bonus: Unified CLI
- Command: `helix …`
- What it does: exposes the DNA (GC + k-mers), spectrum leaderboard, RNA folding, protein summary, viz, and triage workflows from a single entry point.
- Try it: `helix dna --sequence ACGTACGT --k 4 --max-diff 1`
- Remix ideas: script repeatable analyses (e.g., `helix triage --input plasmid.fna --json report.json`) or wire the CLI into notebooks via `subprocess.run`.

## Workflow Runner
- File: `workflows/plasmid_screen.yaml`
- What it does: chains the unified CLI to run DNA summaries, triage reports, spectrum scoring, and RNA folding with one config.
- Try it: `helix workflows --config workflows/plasmid_screen.yaml --output-dir workflow_runs`
- Remix ideas: add `viz` steps for hydropathy plots, or point at your own FASTA files to produce per-sample dashboards overnight.

## PCR Amplicon DAG
- Files: `examples/pcr_primers.json`, `examples/pcr_config.json`
- What it does: feeds the new `helix pcr dag` command with a demo genome + primer pair to generate a `helix.pcr.amplicon_dag.v1` artifact. Pair with `helix edit-dag viz/animate` to render PNG/GIF outputs.
- Try it: `python -m helix.cli pcr dag --genome src/helix/datasets/dna/plasmid_demo.fna --primer-config examples/pcr_primers.json --pcr-config examples/pcr_config.json --out pcr_amplicon_dag.json`
- Remix ideas: chain the PCR DAG after a CRISPR or Prime Edit DAG, or tweak the primers/config to explore how efficiency and error rates reshape the amplicon landscape.
