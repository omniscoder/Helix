# `helix.api`

Research workflows often start in notebooks or lightweight scripts where shelling out to the CLI is inconvenient. The `helix.api` module mirrors the CLI surface area but returns plain Python dictionaries and lists, making it straightforward to serialize to JSON, pass through pandas, or feed into downstream visualization components.

The helpers prefer explicit inputs: every function accepts either an inline `sequence="ACGU..."` **or** an `input_path=Path("sample.fasta")` when DNA/protein IO is relevant. Validation happens up front—invalid bases, overlapping arguments, and missing files raise informative `ValueError` or `ImportError` exceptions so that provenance remains audit-friendly.

## Function reference

### `dna_summary(sequence=None, *, input_path=None, window=200, step=50, k=5, max_diff=1)`

Normalize a DNA string, compute GC statistics, and discover k-mer hotspots tolerant of SNP-scale variation.

**Parameters**

| Name | Type | Details |
| --- | --- | --- |
| `sequence` | `str \| None` | Inline DNA input. Cannot be combined with `input_path`. |
| `input_path` | `str \| Path \| None` | File to read (FASTA/plain). Takes precedence over `sequence`. |
| `window` | `int` | Sliding-window size for GC summaries. Set to `0` to skip. |
| `step` | `int` | Advance between windows; smaller steps increase resolution. |
| `k` | `int` | k-mer size used for clustering recurrent motifs. |
| `max_diff` | `int` | Maximum Hamming distance when grouping similar k-mers. |

**Returns**

`dict` with:

- `sequence`: normalized uppercase DNA (U is converted to T, FASTA headers stripped).  
- `length`: integer length of `sequence`.  
- `gc_content`: float fraction (0–1).  
- `gc_windows`: list of `{start, end, gc_fraction}` windows for downstream plotting.  
- `kmer_clusters`: map keyed by canonical motif with `count`, `positions`, `patterns`.  

The GC window calculation and k-mer clustering match the logic behind `helix dna summarize`, ensuring CLI and notebook reports stay interchangeable.

---

### `triage_report(sequence=None, *, input_path=None, k=5, max_diff=1, min_orf_length=90)`

Produce the complete “triage” bundle (GC skew, motif clusters, ORF calls) that backs `helix viz triage`.

**Parameters**

| Name | Type | Details |
| --- | --- | --- |
| `sequence`, `input_path` | see above | Mutually exclusive DNA inputs. |
| `k`, `max_diff` | `int` | Passed through to the underlying k-mer clustering. |
| `min_orf_length` | `int` | Filter ORFs shorter than the threshold (nt). |

**Returns**

`dict` with:

- `sequence`: RNA display sequence (T→U) for visualization.  
- `skew`: list of cumulative skew values (len = sequence+1).  
- `clusters`: each cluster has `canonical`, `count`, `patterns`, `positions`.  
- `orfs`: list of ORF dicts reporting `start`, `end`, `strand`, `frame`, `length_nt`, `length_aa`, and translated `peptide`.  

All payloads are schema-compatible with `helix schema show triage-report`.

---

### `fold_rna(sequence, *, min_loop_length=3, allow_wobble_pairs=True)`

Convenience wrapper around the annotated Nussinov dynamic program for RNA folding.

**Parameters**

| Name | Type | Details |
| --- | --- | --- |
| `sequence` | `str` | RNA/DNA string; U/T are normalized to U internally. |
| `min_loop_length` | `int` | Nussinov “hairpin” constraint (nt separating paired bases). |
| `allow_wobble_pairs` | `bool` | When `True`, GU wobble pairs are permitted in addition to AU/GC. |

**Returns**

`dict` with:

- `sequence`: normalized RNA string used for folding.  
- `score`: optimal base-pair count reported by the DP table.  
- `pairs`: list of `(i, j)` tuples describing paired indices.  
- `dot_bracket`: canonical dot-bracket representation for plotting or downstream structure comparison.  

Because the helper exposes the same knobs as `helix rna fold`, results are reproducible across interfaces.

---

### `spectrum_leaderboard(peptide=None, *, experimental_spectrum=None, cyclic=True, leaderboard_size=5)`

Run the leaderboard cyclopeptide sequencing algorithm and return notebook-friendly results.

**Parameters**

| Name | Type | Details |
| --- | --- | --- |
| `peptide` | `str \| None` | Optional candidate peptide for generating a theoretical spectrum. |
| `experimental_spectrum` | `Sequence[int] \| None` | Observed masses used for leaderboard selection. |
| `cyclic` | `bool` | Whether the theoretical spectrum should be cyclic (`True`) or linear (`False`). |
| `leaderboard_size` | `int` | Maximum number of peptides kept per iteration (ties are preserved). |

**Returns**

`dict` with:

- `theoretical_spectrum`: the computed linear/cyclic spectrum for `peptide` (empty if no peptide supplied).  
- `leaderboard_hits`: list of `{peptide, score}` pairs sorted by score.  

Pass only `peptide` to preview its spectrum, only `experimental_spectrum` to search for best-scoring sequences, or both to compare expectations vs. observed data.

---

### `protein_summary(sequence=None, *, input_path=None, window=9, step=1, scale="kd")`

Summarize amino-acid sequences using Biopython’s `ProtParam` utilities along with Helix hydropathy profiles.

**Parameters**

| Name | Type | Details |
| --- | --- | --- |
| `sequence`, `input_path` | see above | FASTA headers are handled transparently. |
| `window` | `int` | Sliding window for hydropathy averaging. |
| `step` | `int` | Offset between successive hydropathy windows. |
| `scale` | `str` | Hydropathy scale identifier (e.g., `"kd"` for Kyte-Doolittle). |

**Returns**

`dict` with canonical protein metrics: `sequence`, `length`, `molecular_weight`, `aromaticity`, `instability_index`, `gravy`, `charge_at_pH7`, plus `hydropathy_profile` entries (`start`, `end`, `score`). Raises `ImportError` if Biopython is unavailable so that workflows can fail fast with a clear dependency message.

---

### `run_workflow(config_path, *, output_dir, name=None)`

Execute a YAML workflow definition (the same format consumed by `helix workflows run`) from Python. The helper wraps `helix_workflows.run_workflow_config`, returning the list of materialized artifacts. `name` can restrict execution to a single workflow from a manifest, mirroring the CLI `--name` flag.

## Usage patterns

```python
from pathlib import Path
from helix import api as hx

report = hx.triage_report(sequence="AUGGCCUUUUAA", k=3)
gc_bins = hx.dna_summary(input_path=Path("samples/ecoli.fna"), window=500, step=50)
rna = hx.fold_rna("GGGAAACCC", min_loop_length=0)
peptides = hx.spectrum_leaderboard(
    experimental_spectrum=[0, 113, 128, 227, 242, 355, 370, 484],
    leaderboard_size=10,
)
```

Each function returns plain JSON-serializable structures, making it trivial to call `json.dumps(...)`, ship results to `helix viz ...`, or interoperate with scientific Python stacks. Exceptions surface early and with actionable messages—ideal for research-grade reproducibility and provenance.

Need performance baselines? Run `python -m benchmarks.api_benchmarks --repeat 5 --limit 0 --sort mean --out bench/api.json` from the repo root to capture timing data for every helper (or focus on specific functions via `--scenario`). Set `HELIX_BENCH_DNA_FASTA` / `HELIX_BENCH_PROTEIN_FASTA` to swap in larger genomes or proteomes once available, and use `--limit 10000` to mimic CI’s quicker sampling. Each run emits a schema-tagged payload (`bench_result` v1.0) that logs git SHA, BLAS vendor, CPU/threads, RNG seed, and per-case RSS stats so notebooks + CI can do apples-to-apples comparisons. Compare two runs via `scripts/bench_check.py baseline.json current.json --threshold 5` to flag >5% slowdowns automatically, and browse the rolling history at `docs/benchmarks.md`.
