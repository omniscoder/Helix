# CRISPR Plasmid Panel Benchmark

This benchmark ties Helix's CRISPR and prime editing outcome models to a simple synthetic plasmid library abstraction.

- **Spec template:** `templates/lab_crispr_plasmid_bench.helix.yml`
- **Harness:** `benchmarks/crispr_plasmid_panel.py`

## Spec shape

The panel config declares one or more plasmid targets:

- `reference_sequence` or `reference_fasta` – digital DNA window.
- `guides[]` – CRISPR guides with `sequence`, `pam_profile`, and optional `simulation` overrides (`draws`, `priors_profile`, `seed`).
- `prime_pegs[]` – optional prime editing pegRNAs with `spacer`, `pbs`, `rtt`, and optional `simulation` overrides.
- Each guide/peg may carry a `measurements` block pointing at an amplicon counts table (`counts_file`, `label_column`, `count_column`, optional `timepoint_h`).

## Running the benchmark

Example invocation:

```bash
python -m benchmarks.crispr_plasmid_panel \
  --config templates/lab_crispr_plasmid_bench.helix.yml \
  --out bench-results/plasmid_panel.json
```

The harness:

- Runs `helix.crispr.simulate.simulate_cut_repair` for each CRISPR guide.
- Runs `helix.prime.simulator.simulate_prime_edit` for each prime pegRNA.
- Loads counts (when available), normalizes to probability distributions, and reports L1 distance and Jensen–Shannon divergence between simulated and observed outcome spectra.

The JSON output (`crispr_plasmid_panel` schema) records per-target/per-guide metrics plus aggregate means, suitable for CI tracking or methods figures.
