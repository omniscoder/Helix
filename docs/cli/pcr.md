# CLI: PCR Amplicon DAGs

Helix now ships an in‑silico PCR workflow that produces the same provenance‑stamped DAG artifacts as the CRISPR and Prime editors. Each branch captures primer binding, per-cycle amplification, and optional error events so you can visualize how an amplicon population evolves over time.

This walkthrough uses the bundled plasmid FASTA plus the sample configs in `examples/pcr_*.json`. Everything stays simulation-only—no wet-lab instructions or thermocycler recipes are provided.

---

## Primer + PCR configs

| File | Purpose |
| ---- | ------- |
| `examples/pcr_primers.json` | Defines a primer pair (`forward`/`reverse`) with optional metadata and mismatch tolerances. |
| `examples/pcr_config.json`  | Captures high-level parameters such as cycle count, per-cycle efficiency, error rate, and amplicon length bounds. |

Primer schema (excerpt):

```json
{
  "name": "demo_pcr",
  "forward": {
    "name": "demo_forward",
    "sequence": "ATGCGTACGTTG",
    "max_mismatches": 0
  },
  "reverse": {
    "name": "demo_reverse",
    "sequence": "CGGTCAACGTAC",
    "max_mismatches": 0
  },
  "metadata": {
    "use_case": "demo"
  }
}
```

PCR configuration schema (excerpt):

```json
{
  "cycles": 6,
  "per_cycle_efficiency": 0.85,
  "error_rate": 0.0005,
  "min_amplicon_length": 80,
  "max_amplicon_length": 500,
  "max_amplicons": 32
}
```

---

## Generate an amplicon DAG

```bash
python -m helix.cli pcr dag \
  --genome src/helix/datasets/dna/plasmid_demo.fna \
  --primer-config examples/pcr_primers.json \
  --pcr-config examples/pcr_config.json \
  --out pcr_amplicon_dag.json
```

The output is a standard Helix DAG artifact (`helix.pcr.amplicon_dag.v1`). Nodes store fully materialized genome views, stages (`root`, `binding`, `cycled`, `error`, …), log probabilities, and primer metadata. The artifact feeds directly into:

```bash
python -m helix.cli edit-dag viz --input pcr_amplicon_dag.json --out pcr_dag.png
python -m helix.cli edit-dag animate --input pcr_amplicon_dag.json --out pcr_dag.gif --fps 6
```

---

## Playground link

Want to inspect the DAG interactively? The playground accepts any hosted JSON via the `?json=` parameter. The docs site publishes the demo artifact generated above, so you can open:

```
https://omniscoder.github.io/Helix/playground/?json=assets/viz/pcr_amplicon_dag.json
```

Inside the playground you can:

- Drag to pan/zoom, toggle layouts, or filter by probability/time-step.
- Click nodes to see stage metadata plus root→node sequence diffs.
- Export a PNG snapshot or copy the shareable link with the `?json=` payload.

You can also drag-drop your own `helix.pcr.amplicon_dag.v1` file right into the browser.

---

## Chaining PCR after CRISPR/Prime

Because all edit engines emit compatible DAG artifacts, you can combine workflows end-to-end:

1. Use `helix crispr dag` or `helix prime dag` to model genome edits.
2. Feed the edited genome (or specific branches) into the PCR simulator.
3. Visualize everything via `helix edit-dag viz/animate` or the playground for side-by-side comparisons.

Pair this with `helix workflows` to automate CRISPR → PCR → visualization pipelines in a single YAML config.
