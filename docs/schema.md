# Schema Reference

Helix ships a schema registry so every JSON artifact (CLI output, workflow step, visualization payload) can be validated and traced back to the exact `spec_version`.

- **Current spec version:** see `src/helix/schema/spec_manifest.json` (loaded automatically by the CLI and docs). Version `1.0` covers the visualization payloads listed below.
- **Inspect schemas from the CLI:** `helix viz schema` (lists keys) or `helix viz schema --kind viz_alignment_ribbon` (prints the JSON schema for one payload).
- **Export a manifest:** `helix schema manifest --out schemas.json` writes the spec manifest (including every schema, validator availability, and the current `spec_version`). Ideal for publications or workflow lockfiles.
- **Diff manifests:** `helix schema diff --base old_manifest.json [--target new_manifest.json]` summarizes added/removed/modified schemas (use `--format json` for machine-friendly reports).

## Registered Schemas

All schemas live under the `viz_*` namespace today and correspond to the visualization contracts in [Visualization JSON Schemas](viz.md). Use `helix schema manifest` to obtain the machine-readable form; a truncated excerpt:

```bash
$ helix schema manifest | jq '.schemas.viz_alignment_ribbon.schema.properties | keys'
[
  "cigar",
  "metadata",
  "qry_length",
  "qry_start",
  "ref_length",
  "ref_start"
]
```

Whenever we add optional fields, the spec version bumps to `1.x`. Breaking changes (field removal/renames) would trigger a major bump to `2.0`. The manifest + CLI help make these bumps explicit.

## Workflow Provenance

Workflows can declare schema expectations per step:

```yaml
steps:
  - command: ["seed", "map"]
    args: {...}
    schema:
      kind: viz_alignment_ribbon
      output: map.json
```

When you run `helix workflows ... --with-schema`, Helix validates the artifact, stamps `spec_version` + `input_sha256`, and prints a provenance table summarizing `{step, schema kind, version, SHA256, status}`—handy for reviewers and archival logs.

Pair `helix schema manifest` with `helix workflow` outputs to guarantee every JSON artifact in your pipeline references a declared schema.
