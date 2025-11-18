# Headless CLI & Telemetry

Helix snapshots (`.hxs`) capture Studio sessions: the manifest, session JSON, and per-run snapshots. The headless CLI replays those bundles without a GUI so pipelines (Nextflow, Cromwell, Argo, CI) can launch the same runs, compare results, and export reports.

## Snapshot lifecycle

```bash
# Pack a session + run snapshots into a bundle
helix snapshot pack \
  --session session.json \
  --include-run runs/crispr_run.json \
  --include-run runs/prime_run.json \
  --out bundle.hxs

# Inspect a bundle
helix snapshot inspect --bundle bundle.hxs
```

`helix snapshot inspect` prints run IDs, kinds, and status along with a JSON summary when `--json path` is supplied.

## Running simulations

```bash
# Local CRISPR replay (runs actual simulators)
helix run \
  --snapshot bundle.hxs \
  --run-id run_crispr_A \
  --engine local \
  --out out/crispr_local/

# Prime run on OGN (mock client today, real queue later)
helix run \
  --snapshot bundle.hxs \
  --run-id run_prime_B \
  --engine ogn \
  --queue gpu-a100 \
  --out out/prime_ogn/
```

Each run directory contains:

- `snapshot.json` – canonical run snapshot (state + outcomes)
- `metrics.json` – normalized `RunMetrics` payload
- `report.json` + `report.md` – verdict + textual summary
- `run_id.txt` – the manifest-sourced run identifier
- `telemetry.jsonl` – structured telemetry events (one line per event)

## Compare and report

```bash
# Compare two run outputs
helix compare \
  --run out/crispr_local/snapshot.json \
  --run out/prime_ogn/snapshot.json \
  --out compare/prime_vs_crispr.json

# Export a Markdown report for a run
helix report \
  --run out/prime_ogn/snapshot.json \
  --format md \
  --out reports/prime_run.md
```

`helix compare` accepts run snapshot paths or manifest IDs via `--snapshot bundle.hxs`. The result JSON records both summary metrics and the `run_metrics` delta verdict.

## Telemetry schema

Every headless command appends a JSONL line to `telemetry.jsonl` next to its outputs. Fields are flattened for easy ingestion (DuckDB/SQLite/log pipelines):

- Common fields: `ts`, `event` (`run_completed`, `compare_completed`, `report_completed`), `session_id`, `run_id`, `run_kind`.
- Run events: `engine`, `engine_version` (local), `ogn_queue`/`ogn_job_id` (OGN), `ttfv_ms`, `status`, `verdict_label`, `compare_used`, `report_exported`.
- Compare events: `left_run`, `right_run`, `verdict_label`.
- Report events: `format`, `path`.

These logs feed directly into KPI dashboards (TTFV, compare coverage, report exports) without re-parsing nested JSON.

## Params overrides & sweeps

`helix run --params params.json` deep-merges the provided JSON into the run config. Example (Prime sweep):

```json
{
  "run_config": {
    "draws": 64,
    "seed": 42,
    "pam_profile": "SpCas9_NGG",
    "pam_softness": 0.3,
    "priors_profile": "default_indel"
  }
}
```

Python snippet to sweep PBS lengths:

```python
import json, subprocess
pbs_lengths = [8, 10, 12]
for pbs in pbs_lengths:
    params = {"run_config": {"pbs_length": pbs}}
    path = f"params_pbs{pbs}.json"
    with open(path, "w") as fh:
        json.dump(params, fh)
    out = f"out/pbs{pbs}"
    subprocess.run([
        "helix", "run",
        "--snapshot", "prime_bundle.hxs",
        "--run-id", "prime-1",
        "--engine", "local",
        "--params", path,
        "--out", out,
    ], check=True)
```

Each sweep run emits its own snapshot/metrics/report + telemetry line, making downstream analysis (e.g., in DuckDB) trivial.

## Notes

- `--engine ogn` currently uses `MockOGNClient`; the dataclasses (`OGNJobSpec`, `OGNJobStatus`, `OGNJobResult`) are frozen so the real GPU queue can drop in later.
- Local engine currently supports CRISPR, Prime, and PCR experiment kinds; unsupported kinds raise a clear runtime error.
