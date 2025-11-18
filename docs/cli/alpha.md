# CLI Alpha

The Helix CLI α release focuses on reproducible design loops that mirror Studio sessions: pack snapshots, run CRISPR/Prime/PCR engines locally or on OGN, compare results, and emit reports for sharing or automation.

> Status: Behind the feature flag `HELIX_STUDIO_ALPHA=1`. Commands + spec stabilize for `release/v1.1`.

## Core Commands
### Run
```bash
helix run \
  --snapshot session.hxs \
  --kind Prime \
  --params params.json \
  --engine local \
  --out runs/prime/
```
- `--snapshot` accepts `.hxs` bundles (see [Snapshot Spec v1](../snapshot_spec_v1.md)).
- `--kind` maps to runtime plugins (`Prime`, `CRISPR`, `PCR`).
- `--engine` supports `local` or `ogn`; additional engines register via plugin API.
- Outputs include `run_id.txt`, DAG JSON, metrics summaries, provenance logs.

### Compare
```bash
helix compare \
  --run runs/prime/run_id.txt \
  --run runs/crispr/run_id.txt \
  --out compare/prime_vs_crispr.json
```
- Accepts two run IDs or manifest paths.
- Emits verdict payload (`improved|worse|tradeoff|inconclusive`) + metric deltas.
- Optional `--report compare.md` to emit Markdown table.

### Report
```bash
helix report \
  --run $(cat runs/prime/run_id.txt) \
  --format md \
  --out reports/prime_run.md
```
- Formats: `md`, `html`, `json` (PDF later).
- Pulls snapshot metadata + compare verdict (if present) for context blocks.

### Snapshot Pack & Inspect
```bash
helix snapshot pack \
  --session studio_session.json \
  --out session.hxs \
  --include-runs runs/prime runs/crispr

helix snapshot inspect session.hxs
```
- `pack` accepts Studio-exported session JSON and run folders; produces `.hxs` bundles.
- `inspect` validates manifest, prints run kinds, compare verdicts, and report inventory.

## Remote Execution via OGN
```bash
helix run \
  --snapshot session.hxs \
  --kind CRISPR \
  --engine ogn \
  --queue gpu-a100 \
  --out s3://lab-bucket/helix/A123/
```
- Requires `HELIX_OGN_TOKEN`.
- CLI streams job status + perf overlay when `--watch` is set.

## Pipeline Integration
### Nextflow Module
```groovy
process HELIX_PRIME {
  container 'helixstudio/cli:1.0'
  input:
    path snapshot
    path params_json
  output:
    path "out/**"
    path "report.md"
  script:
  """
  helix run \
    --snapshot ${snapshot} \
    --kind Prime \
    --params ${params_json} \
    --engine ogn \
    --out out/
  helix report --run $(cat out/run_id.txt) --format md --out report.md
  """
}
```

### Cromwell / WDL Task
```wdl
task HelixCRISPR {
  input {
    File snapshot
    File params_json
  }
  command <<<'
    helix run --snapshot ~{snapshot} --kind CRISPR --engine ogn --out out/
    helix report --run `cat out/run_id.txt` --format json --out report.json
  >>>
  output {
    File report = "report.json"
    Directory artifacts = "out"
  }
}
```

## Instrumentation & Telemetry
- Each command emits structured logs with `loop_id`, `ttfv_ms`, `compare_used`, `report_exported` for KPI tracking.
- `--trace` flag sends anonymized telemetry (opt-in, documented in `docs/privacy.md`).

## Roadmap to 1.0
- Improve ergonomic defaults (auto-naming runs, implicit `--snapshot` from cwd session).
- Support `helix compare --diff viz` for inline panel exports.
- `helix report --template custom.md` surfaces Template Gallery presets.
- CLI smoke tests gate merges via `helix cli smoke --engine local` in CI.
