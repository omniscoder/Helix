# CRISPR Off-target Panel Benchmark

This benchmark connects Helix's CRISPR off-target search + scoring engine to assay-like hit counts (e.g., CHANGE-seq / CIRCLE-seq).

- **Spec template:** `templates/crispr_offtarget_panel.helix.yml`
- **Harness:** `benchmarks/crispr_offtarget_panel.py`

## Spec shape

The off-target panel config declares:

- `genome.fasta` – digital DNA reference used for off-target search.
- `defaults` – shared search parameters:
  - `pam` – name of a registered PAM (e.g. `SpCas9-NGG`).
  - `max_mm`, `max_gap` – mismatch / gap limits for enumeration.
- `guides[]` – per-guide entries with:
  - `id`, `sequence`, `pam` override (optional).
  - `search` – optional per-guide search params.
  - `assay` – optional off-target results:
    - `file` – path to a TSV.
    - `columns` – mapping of column names for `chrom`, `start`, `strand`, `read_count`.

## Running the benchmark

Example invocation:

```bash
python -m benchmarks.crispr_offtarget_panel \
  --config templates/crispr_offtarget_panel.helix.yml \
  --out bench-results/crispr_offtarget_panel.json
```

The harness:

- Enumerates off-target candidates using `helix.crispr.score.enumerate_off_targets` over the panel's genome window.
- Scores hits with `helix.crispr.score.score_off_targets` (default weight profile).
- Loads assay hit counts, aggregates by `(strand, start)`, and compares:
  - Per-guide metrics: predicted vs observed site counts and a Pearson correlation between predicted scores and `log1p(read_count)` over overlapping sites.
  - Panel-level summary: number of guides evaluated and mean Pearson correlation.

The JSON output (`crispr_offtarget_panel` schema) is designed for CI drift checks and for supporting methods-text statements about off-target risk calibration.

