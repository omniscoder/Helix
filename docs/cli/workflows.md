# CLI: Workflows

```bash
helix workflows --config workflows/plasmid.yaml --output-dir runs --with-schema --as-json
```

- `--with-schema` prints a table of `{step, schema_kind, spec_version, sha256}`.
- `--as-json` emits the same data as JSON for dashboards.
- Each `schema: {kind, output}` block in the YAML validates output JSON before the next step.
