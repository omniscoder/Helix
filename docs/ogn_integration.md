# Remote Engine Contracts (OGN-ready)

Helix and OGN share the same JSON contracts so local prototypes and remote
pipelines stay interchangeable.

## Performance benchmark response

Remote CRISPR services must return the exact JSON produced by
`helix engine benchmark --json` (see `docs/engine_architecture.md`). Required
fields:

- `helix_version`, `scoring_versions`, `env`, `seed`, `config`
- `benchmarks.crispr[]` entries with `backend_requested/backend_used/shape/mpairs_per_s`
- `benchmarks.prime[]` entries with `backend_requested/backend_used/workload/predictions_per_s`

## Prime physics response

Remote prime scoring endpoints should echo the `physics_score` structure used by
`helix prime simulate --physics-score`:

```json
{
  "physics_score": {
    "pbs_dG": -8.4,
    "rt_cum_dG": [...],
    "flap_ddG": 0.7,
    "microhomology": 6,
    "mmr_flag": false,
    "nick_distance": 4,
    "P_RT": 0.82,
    "P_flap": 0.61,
    "E_pred": 0.50
  }
}
```

If a remote service adds ML-augmented probabilities, place them adjacent to the
same structure so clients can continue to rely on one schema regardless of
backend.
