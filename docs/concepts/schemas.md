# Schemas

- Every CLI command that emits JSON stamps `schema_kind` and `spec_version`.
- `helix schema manifest --out schemas.json` exports the full registry.
- `helix schema diff --base old.json --target new.json` shows additions/removals.
- `helix viz --schema` prints contracts plus sample payloads.

Treat schemas like APIs: bump the version, document diffs, and add workflow validation (`schema: {kind, output}`) whenever a step produces JSON.
