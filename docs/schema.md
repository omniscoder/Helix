# Schema Reference

Helix ships a schema registry so every JSON artifact (CLI output, workflow step, visualization payload) can be validated and traced back to the exact `spec_version`.

- **Current spec version:** see `src/helix/schema/spec_manifest.json` (loaded automatically by the CLI and docs). Version `1.0` covers the visualization payloads listed below.
- **Inspect schemas from the CLI:** `helix viz schema` (lists keys) or `helix viz schema --kind viz_alignment_ribbon` (prints the JSON schema for one payload).
- **Export a manifest:** `helix schema manifest --out schemas.json` writes the spec manifest (including every schema, validator availability, and the current `spec_version`). Ideal for publications or workflow lockfiles.
- **Diff manifests:** `helix schema diff --base old_manifest.json [--target new_manifest.json]` summarizes added/removed/modified schemas (use `--format json` for machine-friendly reports).

## Registered Schemas

Most schemas live under the `viz_*` namespace and correspond to the visualization contracts in [Visualization JSON Schemas](viz.md). Use `helix schema manifest` to obtain the machine-readable form; a truncated excerpt:

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

### CRISPR / Prime simulation payloads

The new digital simulators emit their own schema-tagged payloads so downstream tooling can trust the structure:

- `crispr.cut_events` – produced by `helix crispr genome-sim`, contains `{cas, guide, genome summary, params, events[]}` where every event records the candidate `TargetSite`, cut position, score, and serialized guide/Cas metadata.
- `prime.edit_sim` – produced by `helix prime simulate`, bundles `{editor, peg, genome summary, params, outcomes[]}` where each outcome includes the targeted site, edited sequence, logit score, and a short description (`intended_edit`, `indel_loss`, `no_edit`, etc.).

Edit DAG artifacts join the registry as well:

- `helix.crispr.edit_dag.v1.1` – full node/edge graph from `helix crispr dag`, where every node embeds materialized sequences for its genome view and edges record the `EditEvent` + rule metadata.
- `helix.prime.edit_dag.v1.1` – equivalent artifact for prime editing, emitted by `helix prime dag` with RTT-driven branches.

Both schemas participate in the manifest/validator pipeline, so `helix schema manifest` and `helix viz schema --kind ...` will now list them alongside the existing viz payloads.

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
