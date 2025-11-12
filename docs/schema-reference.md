# Schema Reference

Use `helix schema manifest --out schemas.json` to export the latest registry. Each entry includes:

- `schema_kind` (e.g., `viz_alignment_ribbon`)
- `spec_version`
- JSON Schema (properties + types)

Keep this file under version control for publications or pipelines. Combine with `helix schema diff` to review changes between releases.
