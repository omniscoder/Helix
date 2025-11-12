# CLI: Schema Tools

- `helix viz --schema` → display contract + sample for a specific renderer.
- `helix schema manifest --out schemas.json` → export registry.
- `helix schema diff --base old.json [--target new.json]` → show changes (JSON/table).
- `helix schema manifest --format json` (via `jq`) helps pin schemas in papers.

Use manifests as lockfiles for workflows or publications so data/figure contracts remain auditable.
