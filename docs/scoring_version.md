# Scoring Versions

Helix freezes its scoring numerics with version tokens so downstream tooling
can reason about compatibility.

## Naming

- `crispr_scoring_version` → `crispr-vMAJOR.MINOR.PATCH` (currently `1.0.0`).
- `prime_scoring_version`  → `prime-vMAJOR.MINOR.PATCH` (currently `1.0.0`).

The version string must bump whenever the underlying math changes: new physics
models, tweaked logistic parameters, fixture refreshes, etc. Cosmetic refactors
that do not alter numeric results do **not** require a bump.

## Where the version appears

- `helix engine info` JSON.
- `helix engine benchmark --json` payloads (`scoring_versions` node).
- Prime simulation metadata (`meta.prime_scoring_version`).
- Prime physics scores when `--physics-score` is enabled.

## Contracts

- CRISPR clients must treat `crispr_scoring_version` as part of the benchmark
  schema; CI compares versions before diffing MPairs/s.
- Prime clients must treat `prime_scoring_version` and the
  `physics_score` struct as a unit: if the version changes, re-run benchmark and
  regenerate physics baselines.
