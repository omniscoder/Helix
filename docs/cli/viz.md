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
