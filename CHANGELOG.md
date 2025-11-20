# Changelog

## [Unreleased]
- Added prime engine scaffolding (batch API, golden fixtures, perf harness) with
  matching scoring-version metadata.
- CRISPR artifacts now record `crispr_scoring_version` and
  `crispr_engine_backend`; Prime outputs expose `prime_scoring_version` and
  `prime_engine_backend`.
- Shared `helix.engine.encoding` helper guarantees consistent DNA encoding for
  CRISPR, Prime, native, and future CUDA backends. Added
  `docs/engine_cuda_plan.md` describing the first CUDA milestone.
- Added a native CUDA backend for CRISPR scoring. Building with
  `-DHELIX_ENGINE_ENABLE_CUDA=ON` enables the GPU kernel exposed via
  `backend="gpu"`; falls back cleanly when CUDA is unavailable.

## [0.2.0] - 2024-12-30
- Added schema manifest exports/diffs plus CLI helpers for schema inspection and workflow provenance.
- Introduced deterministic visualization provenance: viz-spec hashing, input SHA-256 capture, and per-image `*.provenance.json` trailers.
- Added `helix demo viz` assets + schema reference docs, `--schema` flag for every viz subcommand, and JSON workflow provenance.
- Locked extras (`viz`, `protein`, `schema`) and trimmed core dependencies for the 0.2.0 release.
- Bootstrapped CRISPR/prime editing scaffolding (PAM registry, guide discovery CLI + schema, off-target enumeration/scoring, cut/repair simulator, and the first CRISPR viz track).

## [0.1.0] - 2024-10-01
- Initial public release of the Helix bioinformatics toolkit.
