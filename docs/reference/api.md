# `helix.api`

These high-level helpers mirror the CLI commands and return JSON-friendly dicts. They are perfect for notebooks or scripts that want schema-ready artifacts without shelling out.

| Function | Description | Returns |
| --- | --- | --- |
| `dna_summary(sequence=None, *, input_path=None, window=200, step=50, k=5, max_diff=1)` | Normalize a DNA sequence, compute GC %, windowed GC, and SNP-tolerant k-mer clusters. | `dict` with `sequence`, `length`, `gc_content`, `gc_windows`, `kmer_clusters`. |
| `triage_report(sequence=None, *, input_path=None, k=5, max_diff=1, min_orf_length=90)` | Run the full triage stack (skew, k-mer clustering, ORF detection). | `dict` containing skew array, clusters, ORFs ready for `helix viz triage`. |
| `fold_rna(sequence=None, *, input_path=None)` | Convenience wrapper around the Zuker-style MFE folding. | `dict` with dot-bracket structure, energy, pairing list. |
| `spectrum_leaderboard(peptide=None, spectrum=None, *, leaderboard_size=5, linear=False)` | Run leaderboard cyclopeptide sequencing. | `dict` listing leaderboard hits, theoretical spectra. |
| `protein_summary(sequence=None, *, input_path=None, window=9, step=1, scale='kd')` | Summarize amino-acid sequences (length, MW, GRAVY, hydropathy windows). | `dict` plus hydropathy windows for plotting. |

## Usage Example

```python
from helix import api as hx

report = hx.triage_report(sequence="AUGGCCUUUUAA", k=3)
print(report["orfs"][0])
```

All outputs are plain dicts/lists, so you can dump them to JSON, feed them to `helix viz ...`, or plug into pandas for further analysis. Each helper accepts either inline `sequence=...` or `input_path=Path(...)` (mutually exclusive). Errors are raised early if inputs are missing or invalid, keeping provenance clean.
