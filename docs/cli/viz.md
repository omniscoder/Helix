# CLI: Visualization

```bash
helix viz minimizers --input minimizers.json --save mins.png
helix viz minimizers --schema   # print schema + sample payload
```

Flags:
- `--save` writes PNG/SVG/PDF.
- `--save-viz-spec` overrides the `.viz.json` path.
- `--schema` bypasses plotting and prints the JSON contract + sample payload.

All viz commands (minimizers, seed-chain, rna-dotplot, alignment-ribbon, distance-heatmap, motif-logo) mirror this interface and emit `<image>.provenance.json`.

# Edit DAG visualizations

```bash
helix edit-dag viz --input dag.json --out dag.png \
  [--layout timeline|spring|...] [--min-prob 0.05] [--max-time 6]
```
Plots any `helix.*.edit_dag.v*` artifact (CRISPR, Prime, PCR) using the cinematic renderer. `--layout` controls node layering (`timeline` default), while `--min-prob`/`--max-time` filter cluttered graphs before drawing.

### HTML reports

```bash
helix edit-dag report --input dag.json --out report.html --png dag.png
```

Produces a single-page HTML summary (dark theme) with:
- Metrics (node/edge counts, entropy, max probability, top outcomes)
- Embedded PNG (optional) or a placeholder if omitted
- Provenance-friendly JSON links

Great for attaching to PRs or sharing with teammates.

### Compare two DAGs

```bash
helix edit-dag compare --a dag_before.json --b dag_after.json --out compare.png \
  --label-a "baseline" --label-b "design-B" --min-prob 0.02 \
  --summary diff.json
```

- Blue nodes/edges appear only in the first artifact; orange ones only in the second; white are shared.
- Uses the same filtering knobs (`--layout`, `--min-prob`, `--max-time`) as `viz`/`animate`.
- Optional `--summary diff.json` writes a JSON report (unique/shared nodes, terminal outcomes, top probabilities).
- Ideal for CRISPR vs PCR overlays, control vs experimental runs, or before/after optimization.
