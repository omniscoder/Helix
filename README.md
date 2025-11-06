# Helix

Helix is a grab‑bag of bioinformatics and computational biology experiments: quick DNA utilities, spectrum toys, structural biology helpers, and even a tiny neural net. The goal of the project is to give us a playground for prototyping analysis ideas, then iterating toward something more polished together.

## Features
- **Motif discovery sandbox** (`bioinformatics.py`): naive k‑mer counting plus GC skew plotting for origin-of-replication hunting.
- **Transcription & translation helpers** (`codon.py`, `amino_acids.py`): the start of a codon translator and mass lookup table for peptides.
- **Cyclo-spectrum stub** (`cyclospectrum.py`): placeholder for peptide spectrum calculations using the amino-acid mass index.
- **Secondary structure experiments** (`nussinov_algorithm.py`): beginnings of an RNA pairing predictor via the Nussinov dynamic program.
- **3D protein viewer utilities** (`protein.py`): helper wrappers around BioPython’s PDB parser and peptide builders.
- **Tiny ANN demo** (`ann.py`): a minimal NumPy-only neural net training loop for XOR-style datasets.

## Layout
```
.
├── amino_acids.py          # peptide mass lookup table
├── ann.py                  # single-layer neural network example
├── bioinformatics.py       # DNA utilities + GC skew plotting
├── bioinformatics.c        # scratch file (currently placeholder)
├── codon.py                # WIP codon translator
├── cyclospectrum.py        # future peptide spectrum tool
├── input/dna/human.txt     # example labeled DNA sequences
├── nussinov_algorithm.py   # Nussinov algorithm skeleton
└── protein.py              # BioPython-powered protein helpers
```

## Getting Started
### Requirements
- Python 3.10+ (3.11 works fine)
- pip / virtualenv
- NumPy, pandas, matplotlib, Biopython (install via pip)

Optional:
- A matplotlib backend that can display windows (for GC skew plots).
- Network access if you plan to fetch PDB structures dynamically.

### Installation
```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas matplotlib biopython
```

## Usage
- **K-mer + skew analysis**
  ```bash
  python bioinformatics.py
  ```
  Edits inside `bioinformatics.py` determine which genome fragment is analyzed. The script prints k-mer frequency dictionaries and pops up a GC-skew plot.

- **Neural net demo**
  ```bash
  python ann.py
  ```
  Trains a tiny single-layer network on a hard-coded dataset and prints the learned weights.

- **Codon utilities**
  ```python
  from codon import translate_rna
  translate_rna("AUGGCC...")
  ```
  (Function currently WIP; expect rough edges.)

- **Protein helpers**
  ```python
  from protein import show_sequence
  show_sequence("1CRN.cif")
  ```
  Requires the corresponding CIF file to be locally available or fetched beforehand.

The `input/dna/human.txt` file contains sample (sequence, class) pairs you can load via pandas for custom experiments.

## Roadmap / Ideas
1. Finish the codon translator and add unit tests for edge cases.
2. Flesh out `cyclospectrum.py` with scoring and leaderboard generation.
3. Complete the Nussinov implementation with proper scoring + traceback visualization.
4. Add CLIs/argparse to each module so users can swap input files without editing code.
5. Package dependencies in `pyproject.toml` and add CI for lint/tests.

## Contributing
Issues, ideas, and PRs are very welcome. This repo is intentionally informal—feel free to sketch new notebooks, add scripts, or refactor existing modules. Just keep commits focused and document assumptions in code or the README so future contributors can follow along.
