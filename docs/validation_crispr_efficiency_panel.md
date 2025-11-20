# CRISPR Efficiency Panel Benchmark

This benchmark ties Helix's CRISPR on-target scoring model to observed edited fractions for a panel of endogenous-like sites.

- **Spec template:** `templates/crispr_efficiency_panel.helix.yml`
- **Harness:** `benchmarks/crispr_efficiency_panel.py`

## Spec shape

The efficiency panel config declares:

- Global `cas` settings (Cas system type, PAM pattern, mismatch weights).
- A `targets` list where each entry defines:
  - `reference_sequence` – digital DNA window for the locus.
  - `guide` – guide `id`, `sequence`, and optional `strand`.
  - `measurements.edited_fraction` – either:
    - `value`: inline edited fraction for quick tests, or
    - `file` + `column`: path to a TSV with an edited fraction column.
  - Optional `chromatin_features` (ATAC scores, histone marks, etc.) passed through to the JSON output.

## Running the benchmark

Example invocation:

```bash
python -m benchmarks.crispr_efficiency_panel \
  --config templates/crispr_efficiency_panel.helix.yml \
  --out bench-results/crispr_efficiency_panel.json
```

The harness:

- Builds a `CasSystem` + `GuideRNA` from the panel's `cas` and per-target `guide`.
- Uses `helix.crispr.simulator.predict_efficiency_for_targets` (the public CRISPR engine batch API) to obtain physics-based on-target scores for each reference window.
- Compares predicted scores to observed edited fractions and reports:
  - Per-target status and absolute error.
  - Panel-level mean absolute error (MAE) and Pearson correlation (when at least two observations are available).

The JSON output (`crispr_efficiency_panel` schema) is suitable for CI drift tracking or methods-text calibration summaries.
