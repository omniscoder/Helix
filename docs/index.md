# Helix Docs (Preview)

Helix is a hobbyist-first bioinformatics playground. Install it in editable mode and try the CLI:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[viz,protein]"
helix --help
```

## Quick Starts
- `helix seed map --ref src/helix/datasets/dna/plasmid_demo.fna --reads src/helix/datasets/dna/plasmid_demo.fna --json map.json`  
  _Plot this result →_ `helix viz alignment-ribbon --input map.json --save map.png` (see [schema](viz.md#alignment-ribbon))
- `helix sketch compare --method hll --fasta-a your_ref.fna --fasta-b your_query.fna --json dist.json`  
  _Plot this result →_ `helix viz distance-heatmap --input dist.json --save dist.png` (see [schema](viz.md#distance-heatmap))
- `helix motif find --fasta your_sequences.fna --width 6 --json motif.json`  
  _Plot this result →_ `helix viz motif-logo --input motif.json --save motif.png` (see [schema](viz.md#motif-logo))
- Need the exact field list? `helix viz schema --kind viz_motif_logo` prints the JSON schema inline.
- Want example PNGs + payloads instantly? `helix demo viz --output demo_viz` renders every visualization and writes the paired `.viz.json` files.
- Grab a sample JSON from [Visualization JSON Schemas](viz.md) (e.g., minimizers or seed-chain) and run `helix viz ... --input sample.json --save sample.png` to experiment with the plotting commands immediately.

## Python API
```python
from helix import api as hx

summary = hx.dna_summary(sequence="ACGTACGT", k=4)
print(summary["kmer_clusters"].keys())
```

Read the full README for philosophy, datasets, and roadmap while the docs site grows.
