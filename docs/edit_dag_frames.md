## Edit DAG Frame Schema

Helix can stream edit-DAG construction as incremental “frames” instead of a single final artifact. Each frame is a JSON object with the same structure as the full DAG payload, but limited to the nodes/edges added at that time step.

### JSONL framing

Frames are emitted as JSON Lines (`--frames path.jsonl` or `--frames -` for stdout). Each line contains:

```json
{
  "kind": "helix.edit_dag.frame.v1",
  "step": 3,
  "new_nodes": {
    "n_cut_0": { "...": "standard EditNode payload" }
  },
  "new_edges": [
    { "source": "n_root", "target": "n_cut_0", "...": "standard EditEdge payload" }
  ],
  "meta": {
    "mechanism": "crispr",
    "guide_name": "HBB_G1"
  }
}
```

- `step`: monotonically increasing integer (0 = root).
- `new_nodes`: mapping from node_id → full node payload (log_prob, metadata, sequences).
- `new_edges`: list of edges referencing those node IDs.
- `meta`: optional run metadata (guide/peg names, mechanism, experiment name, etc.).

### Reconstructing a DAG from frames

1. Initialize empty dictionaries for nodes/edges.
2. For each frame:
   - Merge `frame.new_nodes` into the node map.
   - Append `frame.new_edges` to the edge list.
3. The final DAG is identical to the standard JSON artifact.

Python example:

```python
import json
from collections import OrderedDict

nodes = OrderedDict()
edges = []

with open("run.frames.jsonl") as f:
    for line in f:
        frame = json.loads(line)
        nodes.update(frame["new_nodes"])
        edges.extend(frame["new_edges"])

dag_payload = {
    "artifact": "helix.crispr.edit_dag.v1.1",
    "version": "1.1",
    "nodes": nodes,
    "edges": edges,
    "root_id": "n_root",
}
```

### Semantic metadata for frontends

Edit DAG frames carry additional metadata that GUI and web clients rely on for
coloring and labeling:

- Node `metadata` commonly includes:
  - `stage`: coarse state label (e.g. `"root"`, `"cut"`, `"repaired"`, `"error"`, `"no_edit"`).
  - `time_step`: integer depth / frame index for the node.
  - `desired`: optional boolean flag for “desired” vs “undesired” outcomes.
  - `edit_class`: coarse edit type derived from the event:
    - `"deletion"`, `"insertion"`, `"substitution"`, `"no_cut"`, or `"other"`.
- Edge `metadata` commonly includes:
  - `mechanism`: human-facing mechanism label inferred from the rule name,
    e.g. `"CRISPR cleavage"`, `"NHEJ"`, `"Prime editing"`, `"Flap resolution"`,
    `"No change"`.
  - `edit_class`: the same coarse edit class as the target node.

These fields are optional, but when present they form the contract that:

- Qt’s realtime DAG view and the web Playground both treat nodes as edit
  outcomes and color them by `edit_class`.
- Edge labels and tooltips prefer `mechanism` over raw rule names when
  available.

The CRISPR DAG runtime populates `edit_class` and `mechanism`
automatically for CRISPR edit rules; other mechanisms can opt into the same
conventions so that all frontends share a consistent visual language.

### CLI quickstart

```
helix crispr dag \
  --genome examples/hg19_chr_demo.fa \
  --cas-config examples/cas9.json \
  --guide-sequence ACGTACGTACGTACGTACGT \
  --frames run.frames.jsonl \
  --out out/final.edit_dag.json
```

The JSONL stream is ideal for:
- Real-time visualization (Cytoscape animations).
- Notebooks that want to inspect intermediate probabilities.
- Streaming to remote clients or logs.

Final note: the frame schema is backward-compatible with the artifact payload; if you can parse `helix.crispr.edit_dag.v1.1`, you already know how to parse each `new_nodes`/`new_edges` entry.

### Frames → dataset rows

Use the dataset generator to include real experiments:

```
helix edit-dag generate-dataset --frames-input run.frames.jsonl --n 0 --out dataset.jsonl
```

This reads `run.frames.jsonl`, reconstructs the final DAG, and appends a dataset record (complete with `top_outcomes` and the full artifact). Combine `--frames-input` with `--n` to mix real experiments and synthetic DAGs in a single JSONL.

Example dataset row emitted from a frame stream:

```json
{
  "id": 12,
  "mechanism": "prime",
  "node_count": 14,
  "edge_count": 13,
  "top_outcomes": [
    {"stage": "repaired", "prob": 0.48, "sequence_hash": "bb51c8b5"},
    {"stage": "error", "prob": 0.37, "sequence_hash": "a14cf12e"}
  ],
  "frame_source": "runs/prime_demo.frames.jsonl",
  "artifact": { "...": "helix.prime.edit_dag.v1.1 payload" }
}
```
