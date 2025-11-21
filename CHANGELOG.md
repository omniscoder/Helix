# Changelog

## [Unreleased]
-

## [0.4.0] - 2025-11-20
- Locked the CRISPR physics architecture around three interchangeable backends
  (`cpu-reference`, `native-cpu`, `gpu`), added CLI/env knobs plus Studio status
  wiring, and stamped every CRISPR artifact with
  `crispr_engine_backend`/`crispr_scoring_version` for provenance.
- Brought the Prime engine up to parity with batch APIs, golden fixtures,
  scoring-version metadata, and perf harnesses so CRISPR/Prime outputs share
  the same guarantees.
- Published `docs/engine_architecture.md` + CUDA plan updates covering backend
  behavior, throughput tables (backend × G × N × L), and the shared encoding
  helpers that let CRISPR, Prime, and native/CUDA kernels stay in sync.
- Added an optional CUDA build of `helix_engine._native`
  (`-DHELIX_ENGINE_ENABLE_CUDA=ON`) with GPU harnesses/tests; falls back to
  CPU reference automatically when CUDA hardware is missing.

## [0.2.0] - 2024-12-30
- Added schema manifest exports/diffs plus CLI helpers for schema inspection and workflow provenance.
- Introduced deterministic visualization provenance: viz-spec hashing, input SHA-256 capture, and per-image `*.provenance.json` trailers.
- Added `helix demo viz` assets + schema reference docs, `--schema` flag for every viz subcommand, and JSON workflow provenance.
- Locked extras (`viz`, `protein`, `schema`) and trimmed core dependencies for the 0.2.0 release.
- Bootstrapped CRISPR/prime editing scaffolding (PAM registry, guide discovery CLI + schema, off-target enumeration/scoring, cut/repair simulator, and the first CRISPR viz track).

## [0.1.0] - 2024-10-01
- Initial public release of the Helix bioinformatics toolkit.
