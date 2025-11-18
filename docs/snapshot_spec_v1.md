# Snapshot Spec v1

Helix snapshots (`.hxs`) are portable, content-addressed bundles that freeze a design session: inputs, parameters, deterministic seeds, engine build, run metadata, emitted artifacts, and compare/report decisions. Treat them like Terraform state for biology—immutable contracts you can replay locally, headlessly, or on OGN.

## Goals
- **Deterministic re-hydration:** every snapshot includes engine + schema versions, seeds, and parameter blocks so CRISPR/Prime/PCR runs can be reproduced bit-for-bit.
- **Content integrity:** assets and outputs are referenced via SHA-256 and stored alongside a manifest to detect tampering.
- **Schema-versioned contract:** the manifest declares `snapshot_spec`, `schema_version`, and compatible CLI/Studio versions.
- **Multi-run history:** snapshots can record multiple executions (design, rerun, compare) under a single bundle ID.

## Bundle Anatomy
A snapshot is a ZIP file with the following layout:

```
manifest.json                   # top-level JSON manifest (spec here)
assets/                         # optional static inputs named by SHA-256
  <sha256>                      # e.g., reference genomes, primers, pathway graphs
runs/<run_id>/                  # deterministic run folders (Prime, CRISPR, PCR, etc.)
  params.json                   # engine parameters
  provenance.json               # CLI/Studio invocation metadata
  artifacts/<sha>.json          # DAGs, metrics, verdicts, viz specs
reports/<report_id>/report.md   # rendered report exports referenced in manifest
compare/<compare_id>/diff.json  # structured compare/vs verdict outputs
```

Additional files (plots, FASTA, binary attachments) live anywhere beneath `runs/<run_id>/artifacts/` or `reports/` but must be declared in the manifest with `path`, `kind`, and `sha256`.

## Manifest Schema (v1)
Key fields:

```json
{
  "snapshot_spec": "1.0.0",
  "snapshot_id": "hxs_20240215T104455Z_A1B2",
  "created_at": "2024-02-15T10:44:55Z",
  "studio_version": "1.1.0",
  "cli_version": "0.9.0-alpha",
  "seed": 1337,
  "schema_version": "helix.schema/1.1",
  "runs": [
    {
      "id": "run_crispr_A123",
      "kind": "CRISPR",
      "engine": {
        "name": "helix.crispr",
        "version": "2.3.0",
        "build": "sha256:..."
      },
      "params": {
        "$ref": "runs/run_crispr_A123/params.json",
        "sha256": "..."
      },
      "artifacts": [
        {
          "kind": "helix.crispr.edit_dag.v1.1",
          "path": "runs/run_crispr_A123/artifacts/dag.json",
          "sha256": "...",
          "size": 128734,
          "schema": "helix.crispr.edit_dag.v1.1"
        },
        {
          "kind": "metrics.summary",
          "path": "runs/run_crispr_A123/artifacts/metrics.json",
          "sha256": "..."
        }
      ],
      "status": {
        "outcome": "success",
        "duration_ms": 4520
      }
    }
  ],
  "compare": [
    {
      "id": "cmp_A123_B456",
      "left_run": "run_crispr_A123",
      "right_run": "run_prime_B456",
      "verdict": "improved",
      "path": "compare/cmp_A123_B456/diff.json",
      "sha256": "...",
      "thresholds": {
        "edit_efficiency": "+5%",
        "indel_rate": "-2%"
      }
    }
  ],
  "reports": [
    {
      "id": "report_prime_v1",
      "format": "md",
      "source_run": "run_prime_B456",
      "path": "reports/report_prime_v1/report.md",
      "sha256": "..."
    }
  ],
  "assets": [
    {
      "path": "assets/23b47f....fa",
      "sha256": "23b47f...",
      "kind": "fasta",
      "role": "genome"
    }
  ]
}
```

The full JSON Schema lives at `docs/snapshot_spec_v1.schema.json` (to be published alongside the CLI) and includes validation for enums, timestamps, and deterministic hashing metadata.

## Canonicalization & Hashing Rules
1. **SHA-256 of file contents**; for JSON artifacts we canonicalize before hashing (sorted keys, `\n` newline, UTF-8, no trailing spaces).
2. **Manifest canonicalization:** when hashing the manifest itself, exclude the `manifest_sha256` field (computed last) and use RFC 8785 JSON canonical form.
3. **Deterministic ordering:** arrays like `runs`, `compare`, `reports`, and `assets` must sort by `id` (stable) to ensure deterministic bundle diffs.
4. **Compression:** `.hxs` archives use ZIP with `store` or `deflate`, but the manifest records uncompressed sizes.
5. **Content-address map:** everything referenced by SHA must also exist under `assets/<sha>` or `runs/.../artifacts/<sha>`.

## Worked Examples
### CRISPR Snapshot
`examples/snapshots/crispr_peg.hxs` (forthcoming) showcases:
- `run_crispr_A123` referencing genome FASTA + guide library assets.
- Compare entry vs. a Prime run with verdict `tradeoff` (higher efficiency, more indels).
- Markdown report with embedded PNG that links back to the DAG artifact via `viz_spec`.

### Prime Snapshot
`examples/snapshots/prime_batch.hxs` demonstrates multi-peg batch sweeps with:
- Parameter sweep matrix recorded under `runs[].params.grid`.
- Batch-level metadata referencing `batch_id`, `queue`, and `ogn_submission` fields.
- Attached template gallery preset metadata so Studio can rehydrate UI helpers.

## Tooling & Lifecycle
- `helix snapshot pack --session session.json --out myrun.hxs` gathers a Studio session into a bundle (alpha command, behind `HELIX_SNAPSHOT_ALPHA=1`).
- `helix snapshot inspect myrun.hxs` prints manifest summary, run kinds, compare verdicts, and asset integrity.
- `helix run --snapshot myrun.hxs --kind Prime --params overrides.json` replays a stored run locally or forwards to `--engine ogn` with queue selection.
- Snapshots capture CLI/Studio versions; attempting to replay on an incompatible version surfaces a compatibility warning with remediation instructions.

## Versioning Strategy
- Spec increments patch numbers for additive fields, minor for backward-compatible structural changes, major for breaking layout changes.
- Snapshots declare compatible ranges for Studio + CLI + OGN engines via `compatibility.studio = ">=1.0.0 <2.0.0"` (same for CLI/OGN).
- Migrations live in `helix snapshot migrate --from 1.0.0 --to 1.1.0` and operate purely on manifests + asset references (no binary editing).

## Acceptance Checklist
- [ ] Manifest validates against `snapshot_spec_v1.schema.json`.
- [ ] All assets referenced by path exist and have matching SHA-256.
- [ ] Runs reference schema kinds in the Helix registry.
- [ ] Compare/report entries reference existing runs.
- [ ] `manifest_sha256` appended after building the archive.

Snapshot Spec v1 is the baseline for public sharing; all future CLI and Studio exports must conform starting with `release/v1.1`.
