# Helix Benchmarks

The `benchmarks` package houses reproducible performance measurements for the public `helix.api` helpers and representative workflows. The canonical entry point is:

```bash
python -m benchmarks.api_benchmarks --repeat 5 --warmup 1 --sort mean --json bench/api.json
```

This reports mean/stdev/min/max timings per scenario and writes a rich JSON payload suitable for CI dashboards or notebook drift visualizations.

## Dataset overrides

By default, the benchmarks rely on the lightweight demo datasets bundled with Helix. To stress-test heavier inputs without committing large files, drop them under `benchmarks/data/` (kept out of version control via `.keep`) or anywhere on your machine and export the following environment variables before running the suite:

| Environment variable | Purpose |
| --- | --- |
| `HELIX_BENCH_DNA_FASTA` | Path to a larger DNA FASTA/plaintext file used for `dna_summary`, `triage_report`, and derived RNA folding samples. |
| `HELIX_BENCH_PROTEIN_FASTA` | Path to a protein FASTA used for `protein_summary` and hydropathy windows. |

When unset, the defaults fall back to `src/helix/datasets/dna/plasmid_demo.fna` and `src/helix/datasets/protein/demo_protein.faa`. Override paths are recorded in the JSON metadata so future dashboards can correlate regressions with input changes.

As we collect research-grade reference genomes and proteomes, drop them under `benchmarks/data/` (ignored from source control by default) and point the env vars at those files to drive stress runs locally, via `workflow_dispatch`, or on self-hosted runners.

## JSON schema

Benchmark runs emit a stable payload shaped as:

```json
{
  "schema": {"kind": "bench_result", "spec_version": "1.0"},
  "run": {
    "timestamp": "...Z",
    "commit": "abcdef",
    "branch": "main",
    "repeat": 3,
    "warmup": 1,
    "dataset": {
      "dna_fasta": "/abs/path/plasmid_demo.fna",
      "protein_fasta": "/abs/path/demo_protein.faa",
      "dna_size_bp": 304,
      "protein_count": 1
    }
  },
  "env": {
    "python": "3.11.9",
    "helix_version": "0.2.0",
    "numpy": "2.1.3",
    "platform": "...",
    "cpu": "AMD EPYC ...",
    "threads": 16,
    "blas": "openblas,...",
    "locale": "C.UTF-8",
    "seed": 1337,
    "omp_threads": "4",
    "mkl_threads": "4",
    "git_commit": "...",
    "git_branch": "..."
  },
  "cases": [
    {
      "name": "helix.fold_rna",
      "params": {"input_length_nt": 400, "min_loop_length": 3, "...": "..."},
      "n": 3,
      "time_s": {"mean": 0.97, "std": 0.01, "min": 0.95, "max": 0.99},
      "rss_mb": {"peak": 1120.0},
      "throughput": {"items_s": 1.02},
      "status": "ok",
      "delta_vs_baseline_pct": 1.2
    }
  ]
}
```

That schema enables apples-to-apples comparisons across machines and datasets.

## Drift checks & CI

- Use the `--limit` flag (or `limit` output from CI) to cap how many nucleotides/amino acids are processed. This keeps CI quick (`--limit 10000`) while still allowing full runs (`--limit 0`) during heavy sweeps.
- `.bench/baseline.json` stores the reference measurements committed to the repo. Update it whenever you intentionally change performance characteristics.
- `scripts/bench_check.py baseline current --threshold 5` compares current runs against the baseline and exits non-zero when the mean runtime regresses by more than the given percentage.
- The GitHub Actions `benchmarks` job pins thread counts (`OMP_NUM_THREADS=4`, `MKL_NUM_THREADS=4`), seeds RNGs, runs the suite (repeat=3 by default, repeat=10 when `BENCH_HEAVY=true`), uploads both `latest.json` and `bench-$SHA.json`, and appends the Markdown summary emitted by the benchmark runner to the workflow summary tab.
- Trigger the heavy datasets path by firing `workflow_dispatch` with `bench_heavy=true` and providing absolute paths to the large FASTA files (they must already exist on the runner or self-hosted machine).

Raw JSON artifacts accumulate under `benchmarks/out/` (ignored from git) so you can build static dashboards or notebooks that plot trends over time.
