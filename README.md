# Helix

Helix is a hobbyist-first playground for bioinformatics and computational biology. Think of it as a backpack full of lightweight tools, algorithms, and experiments you can remix on evenings, in classrooms, or between lab runs. We embrace rough edges and fast iteration so that ideas can leap from a notebook sketch to a runnable prototype quickly.

Helix complements our production platform OGN rather than competing with it. When a prototype proves its value, you can polish and port it into OGN. Until then, Helix is the sandbox where curiosity rules.

## Why Helix Exists
- Lower the barrier to tinkering: ship batteries-included examples and tiny datasets.
- Showcase approachable implementations of classic algorithms so learners can peek under the hood.
- Encourage sharing and remixing of exploratory workflows without the ceremony of production deployments.
- Offer a bridge to OGN by keeping APIs compatible and providing off-ramps when users need industrial-scale tooling.

## Highlights
- **DNA and motif experiments** (`bioinformatics.py`): quick-and-dirty k-mer counting, SNP-tolerant motif clustering, GC skew plots, and small FASTA helpers.
- **Translation and mass lookups** (`codon.py`, `amino_acids.py`): resilient codon translation, bidirectional ORF scanning, frameshift heuristics, and peptide mass utilities.
- **Peptide spectrum sandbox** (`cyclospectrum.py`): entry point for building leaderboard scoring and spectrum analysis.
- **RNA secondary structure sketches** (`nussinov_algorithm.py`): an annotated Nussinov dynamic-programming prototype.
- **Protein helpers** (`protein.py`): friendly wrappers around Biopython for peeking at sequences and structures.
- **Neural net doodles** (`ann.py`): minimal NumPy-only network for experimenting with small bio datasets.

## Repo Layout
```
.
├── amino_acids.py          # peptide mass lookup table
├── ann.py                  # single-layer neural network example
├── bioinformatics.py       # DNA utilities + GC skew plotting
├── bioinformatics.c        # scratch file (currently placeholder)
├── codon.py                # codon translation helpers
├── cyclospectrum.py        # future peptide spectrum tool
├── input/dna/human.txt     # example labeled DNA sequences
├── nussinov_algorithm.py   # Nussinov algorithm skeleton
└── protein.py              # Biopython-powered protein helpers
```

## Getting Started
### Requirements
- Python 3.10+ (3.11 tested)
- pip or another package manager
- NumPy, pandas, matplotlib, Biopython

Optional extras:
- A matplotlib backend capable of rendering windows (GC skew plots).
- Network access for fetching PDB or AlphaFold structures on the fly.
- pytest (if you want to run the test suite).

### Installation
```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas matplotlib biopython
```

### Run a Script
- **K-mer + skew analysis**
  ```bash
  python bioinformatics.py
  ```
  Edit the constants at the top of `bioinformatics.py` to point at different sequences and tweak window sizes.
  For quick clustering with exports, try `python examples/kmer_counter.py --max-diff 1 --csv clusters.csv --plot-top 10`.

- **Neural net demo**
  ```bash
  python ann.py
  ```
  Prints training progress and final weights for a tiny XOR-style problem.

- **Translate a sequence**
  ```bash
  python examples/translate_sequence.py AUGGCCUUU
  ```
  Add `--no-stop` to continue through stop codons or point to a file with `--input`.

- **Find ORFs**
  ```bash
  python examples/find_orfs.py --min-length 90 --include-partial --detect-frameshifts --input your_sequence.fna --orf-fasta peptides.faa --orf-csv orfs.csv --frameshift-csv shifts.csv
  ```
  Prints coordinates, frames, strands, optional frameshift candidates, and can export FASTA/CSV artifacts.

- **Triage report CLI**
  ```bash
  python examples/triage_report.py --input your_sequence.fna --output triage.png --clusters-csv clusters.csv --orfs-csv orfs.csv
  ```
  Generates a composite plot plus optional CSV/FASTA exports for quick daily snapshots.

- **Notebook triage dashboard**
  Open `notebooks/triage_dashboard.ipynb` to plot GC skew, ORFs, and k-mer hotspots side-by-side for a quick daily scan.

- **Protein sequence peek**
  ```python
  from protein import show_sequence
  show_sequence("1CRN.cif")
  ```
  Requires the target structure file in the working directory (or adjust the loader).

Browse task-specific quickstarts in `examples/README.md`. `input/dna/human.txt` ships a toy labeled dataset for quick experiments with pandas or sklearn.

### Run Tests
```bash
pytest
```
Pytest powers translator and k-mer regression checks; feel free to add more as you create new helpers.

## Weekend Project Ideas
- Plot the GC skew for a bacterial plasmid and compare predicted origins to literature.
- Extend the ORF scanner to sweep reverse complements and test on viral genomes.
- Compare frameshift candidates against known gene models to flag likely sequencing errors.
- Pair the ORF scanner with the GC skew plot to compare predicted origins and coding regions.
- Use the CSV/plot outputs from `examples/kmer_counter.py` to highlight SNP hotspots and share charts with the community.
- Customize `notebooks/triage_dashboard.ipynb` with your own sequences and publish the visuals for lab updates.
- Hook `cyclospectrum.py` into a simple leaderboard scorer and visualize the mass differences.
- Swap the activation function in `ann.py`, log loss curves, and document what changes.
- Build a notebook that fetches a PDB entry, prints its sequence via `protein.py`, and sketches the secondary structure counts.
- Chain `examples/translate_sequence.py` with `peptide_mass_lookup.py` to score translated open reading frames.

Browse ready-to-run snippets in `examples/README.md`, and share your results in `examples/` (add new files freely) or link to gist/notebook URLs in issues so others can remix.

## Design Philosophy
- **Approachable first**: readable code, inline comments when helpful, datasets that fit in memory.
- **Composable**: functions return plain Python data structures so you can plug them into pandas, NumPy, or future OGN pipelines.
- **Biopython-friendly**: we stand on Biopython's shoulders; no wheel reinvention when a stable API exists.
- **Prototype-to-production bridge**: helper scripts should make it easy to migrate successful ideas into OGN when the time comes.

## Roadmap
1. Bundle a CLI command/notebook for combining GC skew, ORFs, and motif clusters into shareable reports.
2. Implement scoring for cyclo-spectrum experiments and publish a walkthrough notebook.
3. Finish the Nussinov traceback to output secondary structure strings and diagrams.
4. Add small CLIs (argparse or Typer) for swapping inputs without editing source files.
5. Draft an `examples/` gallery featuring community notebooks and weekend projects.

## Relationship to OGN
Helix is intentionally lightweight. We do not guarantee production stability, large-scale data orchestration, or SLA-backed support. When your prototype needs robustness, data governance, or integration with lab automation:
1. Package the core logic (functions, notebooks, scripts).
2. Identify equivalent building blocks in OGN or write adapters that call into it.
3. Open an OGN ticket or PR referencing the Helix prototype so we can collaborate on the migration.

This separation keeps Helix nimble while letting OGN remain the home for hardened workflows.

## Contributing
We welcome ideas, experiments, and docs improvements. To keep things playful:
- Open issues with context, references, or notebooks that inspired your idea.
- Tag contributions by complexity (`good-first-experiment`, `deep-dive`, etc.).
- Respect the code of conduct (be kind, give credit, document assumptions).
- If you plan a larger refactor, start a discussion thread so we can pair-program or offer pointers.

Happy hacking!
