# Edit DAG Visualization Toolkit

Helix emits deterministic DAG artifacts (`helix.crispr.edit_dag.v*`, `helix.prime.edit_dag.v*`, `helix.pcr.amplicon_dag.v*`). This page documents the premium rendering tools that turn those JSON artifacts into publication-ready figures, GIFs, and shareable reports.

## Rendering options

### CLI PNG render

```bash
helix edit-dag viz --input dag.json --out dag.png \
  --layout timeline \
  --min-prob 0.02 \
  --max-time 8
```

- `--layout`: `timeline` (default) for left→right stages, or `spring`, `kamada-kawai`, etc. for exploratory plots.
- `--min-prob`, `--max-time`: filter low-probability branches or truncate long simulations before plotting.

Output:
- `dag.png` — cinematic lineage figure (stage-colored rings, probability-weighted edges, top-outcome inset).
- `dag.viz.json` — viz-spec with schema metadata.
- `dag.provenance.json` — cryptographic audit trail (schema kind, spec version, checksums, CLI command).

### CLI GIF animation

```bash
helix edit-dag animate --input dag.json --out dag.gif \
  --layout timeline --fps 6 \
  --min-prob 0.02 --max-time 8
```

Produces a frame-per-time-step GIF with glowing highlights for emerging branches, temporal bands, and probability-aware arrow widths.

### Playground

https://omniscoder.github.io/Helix/playground/

- Drag/drop any DAG JSON.
- Use the toolbar buttons to pull the latest CRISPR, Prime, or PCR demos (generated directly from the simulator). You can also append `?demo=crispr|prime|pcr` or `?json=URL` to auto-load a DAG.
- Sliders for probability/time filters.
- Hover tooltips show stage metadata and root→node diffs.
- Export PNG or copy a shareable permalink.
- Rebuild the hosted demos anytime with:

```bash
PYTHONPATH=src python3 scripts/generate_playground_dags.py
```

### Side-by-side compare

```bash
helix edit-dag compare --a dag_before.json --b dag_after.json --out compare.png \
  --label-a "CRISPR" --label-b "PCR" \
  --min-prob 0.02 --max-time 8 \
  --summary diff.json
```

- Unique-to-A branches glow teal, unique-to-B branches glow amber, and shared paths stay neutral.
- Timeline layout + probability filters keep dense DAGs readable.
- Inset panel shows the top terminal nodes from each artifact (with normalized probabilities) so you can quantify shifts at a glance.
- `--summary diff.json` writes a JSON diff (unique/shared nodes, terminal outcomes, top probabilities) for changelog/PR use.

Perfect for baseline vs redesign comparisons, CRISPR → PCR chaining, or “what-if” pegs.

### Dataset generator

```bash
helix edit-dag generate-dataset --mechanism crispr --n 25 --out dag_data.jsonl --topk 5
```

- Emits JSONL where each row bundles the DAG artifact plus top outcomes/metadata.
- Mechanisms: `crispr` (default) or `prime`.
- Deterministic via `--seed`.
- Handy for training downstream models or stocking the Playground with sample DAGs.

## HTML report generator

```bash
python -m helix.cli edit-dag report \
  --input dag.json \
  --out report.html \
  --png dag.png
```

Generates a single HTML artifact with:

1. Summary metrics (entropy, intended vs unintended probability, node/edge counts).
2. Embedded PNG/GIF (or fetches from `--png`).
3. Top outcome table with diffs and probabilities.
4. JSON download links for transparency.

You can host the HTML as-is or embed it in docs dashboards.

## Artifact schema v1.1 highlights

- Deterministic node IDs (content-addressed hash of parents + event + stage + time).
- `seq_hashes` and `diffs` per node for compact, tamper-evident storage.
- Metadata: `schema_version`, `rule_version`, `created_at`.
- Fully compatible with earlier tools; `dag_from_payload` reconstructs views lazily.

## Recommended workflow

1. Run simulation (e.g., `helix crispr dag ... --json dag.json`).
2. Render PNG/GIF (`helix edit-dag viz/animate`).
3. Generate HTML report (`helix edit-dag report`).
4. Share via Playground (host `dag.json` and link `?json=`).

This pipeline yields reproducible, explainable edit twins that teams can inspect, compare, and discuss.
