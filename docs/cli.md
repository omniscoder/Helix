# CLI Cheat Sheet

All commands are subcommands of the `helix` console entry point. The most common flags mirror the Python API arguments.

| Command | What it does |
| --- | --- |
| `helix dna --input path.fna` | GC stats + k-mer clusters |
| `helix triage --input path.fna --json triage.json` | Combined GC/k-mer/ORF report (JSON + optional plots) |
| `helix spectrum --peptide NQEL --spectrum "…"` | Theoretical spectra + leaderboard scoring |
| `helix rna mfe --fasta seq.fna --dotbracket out.dbn` | Zuker-style MFE folding |
| `helix rna ensemble --fasta seq.fna --gamma 1.0` | Partition function + MEA/centroid + dot-plot/entropy |
| `helix protein --input protein.faa` | Protein summaries/hydropathy (requires Biopython) |
| `helix viz triage --json triage.json --output triage.png` | Plot triage payloads (requires matplotlib) |
| `helix viz hydropathy --input src/helix/datasets/protein/demo_protein.faa --window 11` | Plot hydropathy profile (Biopython + matplotlib) |
| `helix string search --pattern GATTACA --k 1 seqs.fna` | FM-index exact search (`k=0`) or ≤k Myers hits with JSON output |
| `helix seed index --method minimizer --k 15 --window 10 seq.fna` | Emit minimizers/syncmers + optional density plots |
| `helix seed map --ref ref.fna --reads reads.fna --k 15 --window 10` | Toy seed-and-extend mapping summary |
| `helix dbg build --reads reads.fna --k 31 --graph dbg.json` | Build a DBG + optional GraphML |
| `helix dbg clean --graph dbg.json --out dbg_clean.json` | Remove tips/bubbles via CLI |
| `helix dbg color --reads sample1.fna sample2.fna --k 31` | Produce colored DBG presence JSON |
| `helix motif find --fasta seqs.fna --width 8 --solver steme` | PWM discovery via EM/STEME/online + optional JSON/plot |
| `helix sketch build --method minhash --fasta seq.fna --k 21 --size 1000` | Build a MinHash or HLL sketch |
| `helix sketch compare --method hll --fasta-a a.fna --fasta-b b.fna` | HLL Jaccard/cardinality or Mash distance |
| `helix workflows --config workflows/plasmid_screen.yaml` | Run YAML-defined pipelines |

Tip: add `PYTHONPATH=src` when running from the repo root or install with `pip install -e .`.
