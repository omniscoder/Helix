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

### Interactive edit-DAG demos

Two end-to-end recipes ship in `examples/` so you can go from FASTA/configs → DAG → Playground in seconds:

- **CRISPR demo (guide library)**

  ```bash
  helix crispr dag \
    --genome examples/hg19_chr_demo.fa \
    --cas-config examples/cas9.json \
    --guides-file examples/guides.tsv \
    --region chrDemo:1-500 \
    --out-dir out/crispr_dags/
  ```

  Artifacts appear in `out/crispr_dags/` and a prebuilt copy lives at [`docs/data/crispr_demo.edit_dag.json`](docs/data/crispr_demo.edit_dag.json).
  [Open in Playground](playground/?json=docs/data/crispr_demo.edit_dag.json)

- **Prime editing demo (peg sweep)**

  ```bash
  helix prime dag \
    --genome examples/hg19_chr_demo.fa \
    --editor-config examples/prime_editor.json \
    --pegs-file examples/pegs.tsv \
    --region chrDemo:1-500 \
    --out-dir out/prime_dags/
  ```

  Produces `helix.prime.edit_dag.v1.1` artifacts per peg (see `docs/data/prime_demo.edit_dag.json` for a hosted example).
  [Open in Playground](playground/?json=docs/data/prime_demo.edit_dag.json)

- **Experiment configs**

  ```bash
  helix experiment new --type crispr --out experiments/demo_crispr.helix.yml
  helix experiment run --config experiments/demo_crispr.helix.yml --out out/demo_crispr.edit_dag.json
  helix experiment report --config experiments/demo_crispr.helix.yml --out out/demo_crispr.html
  ```

  A single YAML file describes the entire “edit experiment” (genome, Cas/Prime configs, guides/pegs, simulation knobs). Helix can regenerate DAGs, PNGs, and HTML reports from that spec at any time; see `templates/*.helix.yml` for starter files.

### Live realtime Playground

- [Realtime CRISPR Simulator](playground/realtime.html) — explore guides, stream DAG frames, and compare two designs interactively (browser-only, no backend required).
  - Suggested demo: Guide A `ACCCAGGAAACCCGGGTTTT` vs Guide B `TTTACCCAGGAAACCCGGGT` (chrDemo) to see how probability mass shifts live.
  - Branch chart + entropy panel update every frame (`intended`, `indel`, `no-edit`) so you can watch probability flow in real time.
  - Export frames JSONL and hand them to `helix edit-dag generate-dataset --frames-input …` to append real experiments to your ML corpora (see [Edit DAG Frames](edit_dag_frames.md)).
  - Prefer a desktop shell? `pip install veri-helix[gui] && helix gui` launches the PySide6 version with the same streaming frames and Cytoscape canvas.

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
- Edit DAG internals → [Edit DAG Overview](edit_dag_overview.md)
- Schema Reference (all kinds, versions, and samples) → [Schema Reference](schema-reference.md)
- Performance Dashboard (CI trends for time + memory) → [Benchmarks](benchmarks.md)
- What’s New (highlights across releases) → [What’s New](whats-new.md)
- Contributing (style, tests, ideas) → [Contributing](contributing.md)
- Realtime Playground → [Live Simulator](playground/realtime.html)

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
