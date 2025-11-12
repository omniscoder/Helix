# Workflows

`helix workflows` runs YAML-defined pipelines:

```yaml
steps:
  - command: dna
    args: { sequence: ACGT..., k: 3 }
    schema:
      kind: viz_minimizers
      output: dna.json
```

- `schema` blocks validate step outputs and stamp `spec_version` + SHA-256.
- `helix workflows --with-schema --as-json` prints provenance tables suitable for dashboards.
- Logs and artifacts live under `--output-dir`, making it easy to archive runs.
