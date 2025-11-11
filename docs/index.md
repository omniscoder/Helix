# Helix Docs (Preview)

Helix is a hobbyist-first bioinformatics playground. Install it in editable mode and try the CLI:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[viz,protein]"
helix --help
```

## Quick Starts
- `helix dna --sequence ACGTACGT --k 4`
- `helix spectrum --peptide NQEL --spectrum "0,113,114,128,227,242,242,355,356,370,371,484"`
- `helix rna --sequence GGGAAACCC`
- `helix viz triage --json report.json --output report.png`

## Python API
```python
from helix import api as hx

summary = hx.dna_summary(sequence="ACGTACGT", k=4)
print(summary["kmer_clusters"].keys())
```

Read the full README for philosophy, datasets, and roadmap while the docs site grows.
