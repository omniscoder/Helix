# Artifacts & Provenance

Veri-Helix records a provenance trail for every artifact:

- Schema manifests define the JSON contract (kind + `spec_version`).
- Viz-spec JSON captures metrics about each plot plus `input_sha256`.
- `<image>.provenance.json` pairs schema kind, viz-spec hash, and image hash.

Together these form a chain of custody you can cite in papers or audits. For implementation details see `helix schema manifest`, `helix schema diff`, and `helix workflows --with-schema`.
