# Veri‑Helix

Helix is a research‑grade bioinformatics toolbox that “proves what it plots.” Every figure and workflow carries schema validation, versioned specs, and cryptographic provenance so you can trust what you share and reproduce what you see.

- Install: `pip install "veri-helix[viz,schema,protein]"`
- First look: `helix demo viz --all`
- Why it matters → [Artifacts & Provenance](concepts/artifacts.md)

---

## Quickstart

Run a triage snapshot (GC skew + k‑mer hotspots + ORFs) and save a ready‑to‑share plot:

```bash
helix triage --input path/to/sequence.fna --k 5 --min-orf-length 90 --json triage.json
helix viz triage --json triage.json --output triage.png
```

Fold RNA and export dot‑bracket and ensemble plots:

```bash
helix rna mfe --fasta path/to/seq.fna --dotbracket mfe.dbn
helix rna ensemble --fasta path/to/seq.fna --gamma 1.0 \
  --dotplot dotplot.png --arc arc.png --entropy-plot entropy.png
```

Notebook‑friendly API (plain dict/list payloads):

```python
from helix import api as hx
report = hx.triage_report(sequence="AUGGCCUUUUAA", k=3)
fold   = hx.fold_rna("GGGAAACCC", min_loop_length=0)
```

---

## Deep Wiki (Project Map)

Use this “deep wiki” map to jump directly to the concept or tool you need:

- Getting Started: installation, CLI entry points, extras → [Getting Started](getting-started.md)
- Concepts: how Helix treats data and evidence
  - Artifacts & provenance (hashes, versioned specs) → [Artifacts & Provenance](concepts/artifacts.md)
  - Schemas (contracts for JSON/viz payloads) → [Schemas](concepts/schemas.md)
  - Workflows (YAML pipelines + schema hooks) → [Workflows](concepts/workflows.md)
- CLI Guides
  - Visualization commands → [CLI: Visualization](cli/viz.md)
  - Schema tools (show/spec/manifest/diff) → [CLI: Schema Tools](cli/schema.md)
  - Workflows runner and config → [CLI: Workflows](cli/workflows.md)
  - Demos and fixtures → [CLI: Demos](cli/demo.md)
- API Reference (notebook‑first helpers) → [helix.api](reference/api.md)
- Visualization Gallery (spec‑stamped figures)
  - Alignment Ribbon → [Alignment Ribbon](visualization/alignment-ribbon.md)
  - RNA Dot‑Plot → [RNA Dot‑Plot](visualization/rna-dotplot.md)
  - Motif Logo → [Motif Logo](visualization/motif-logo.md)
  - Distance Heatmap → [Distance Heatmap](visualization/distance-heatmap.md)
  - Minimizer Density → [Minimizer Density](visualization/minimizer-density.md)
- Schema Reference (all kinds, versions, and samples) → [Schema Reference](schema-reference.md)
- Performance Dashboard (CI trends for time + memory) → [Benchmarks](benchmarks.md)
- What’s New (highlights across releases) → [What’s New](whats-new.md)
- Contributing (style, tests, ideas) → [Contributing](contributing.md)

---

## Benchmarks & Reproducibility

Helix includes a research‑grade benchmark harness. Every run emits a schema‑tagged JSON that captures git SHA, BLAS vendor, CPU/threads, RNG seed, per‑case timing, and RSS peaks. CI appends results to `docs/data/bench/history.csv`; the live dashboard renders trends without servers or external deps.

- See live trends → [Benchmarks](benchmarks.md)
- Harness CLI → `python -m benchmarks.api_benchmarks --help`

---

## Design Principles

- Evidence‑first: every artifact ships with a contract (schema) and provenance (hashes, spec version, parameters).
- Friendly internals: approachable implementations of classic algorithms, with JSON payloads you can inspect and remix.
- CLI ⇄ API symmetry: what you can do from the shell, you can reproduce in a notebook—and vice versa.
- Research‑grade: optional heavy datasets, performance dashboards, and drift gates in CI for trustworthy iterations.
