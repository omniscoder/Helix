# Prime Editing Panel Benchmark

This benchmark ties Helix's prime editing outcome model to observed outcome spectra for a panel of targets and pegRNAs.

- **Spec template:** `templates/prime_editing_panel.helix.yml`
- **Harness:** `benchmarks/prime_panel.py`

## Spec shape

The prime panel config declares:

- Global `defaults` for simulation (`draws`, `priors_profile`, `seed`).
- A `targets` list where each entry defines:
  - `reference_sequence` – digital DNA window around the edited site.
  - `peg` – pegRNA fields:
    - `id`, `spacer`, `pbs`, `rtt`
    - optional `simulation` overrides per target.
  - `measurements` – optional outcome counts table:
    - `counts_file` – path to a CSV.
    - `label_column`, `count_column` – per-row outcome label and integer count.
    - optional `timepoint_h` metadata.

## Running the benchmark

Example invocation:

```bash
python3 -m benchmarks.prime_panel \
  --config templates/prime_editing_panel.helix.yml \
  --out bench-results/prime_panel.json
```

The harness:

- Runs `helix.prime.simulator.simulate_prime_edit` on each reference window with the configured pegRNA and prime priors profile.
- Aggregates outcome probabilities by label and, when counts are available, compares simulated vs observed distributions:
  - Per-target metrics: L1 distance and Jensen–Shannon divergence between simulated and observed label distributions.
  - Panel-level summary: mean L1 and JS divergence across all targets.

The JSON output (`prime_editing_panel` schema) is intended for CI drift monitoring and for supporting methods-text claims about prime editing model calibration.

