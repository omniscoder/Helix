# Live Runtime

Helix ships an experimental LiveGraph runtime so you can iterate on multi-rate, SCC-condensed HGX models without leaving the CLI. The `helix live run` command ingests an `.hgx` graph, builds solver islands, executes them with the `LiveScheduler`, and writes a COMBINE-style bundle for reproducibility.

### Quick start

```
# 1. Run the synthetic demo graph and stream realtime deltas
helix live run --hgx examples/live_demo.hgx --duration 5 --realtime

# 2. Attach the viewer (while the run above is active)
helix live viz --endpoint tcp://127.0.0.1:8765
# or replay a bundle after the run finishes:
# helix live viz --bundle live_runs/LiveDemo-*

# 3. Inspect the bundle that was written during the run
helix live inspect live_runs/LiveDemo-* --json live_demo_summary.json
```

```bash
helix live run \
  --hgx examples/p53_oscillator.hgx \
  --duration 2.0 \
  --sync-dt 0.5 \
  --default-dt 0.1 \
  --dt p53=0.05 --dt mdm2=0.05 \
  --input dna_damage.energy=0.2 \
  --bundle out/p53_live_run
```

## Arguments

- `--hgx`: path to the HGX/YAML LiveGraph definition (e.g., `examples/live_demo.hgx`; required).
- `--duration`: simulation time horizon.
- `--sync-dt`: operator-split sync cadence between islands.
- `--default-dt`: fallback Δt for any island that does not have an override.
- `--dt node=Δt`: override a node’s Δt (repeat for multiple). The smallest Δt across the SCC sets the island step.
- `--input node.port=value`: inject constant inputs (repeat as needed).
- `--inputs-json path`: load constants or time-series stimuli from JSON (see below).
- `--seed`: deterministic seed for the runtime.
- `--hz`: target UI cadence for the `StateReducer`.
- `--wall-ms`: budget (ms) for the pacing sleep inside the scheduler loop.
- `--bundle`: output directory. Defaults to `live_runs/<model>-<timestamp>` so you can archive runs by hash.
- `--realtime`: emit `helix.live.delta.v1` payloads during the run.
- `--realtime-endpoint`: optional socket/IPC endpoint (e.g., `tcp://127.0.0.1:8765` or `ipc:///tmp/helix.sock`) for external consumers such as `helix live viz`. Defaults to `tcp://127.0.0.1:8765`.

### Time-series inputs

The JSON format maps `node.port` names to either a constant or a list of `{ "t": <time>, "value": <number> }` points. Values are piecewise constant and hold until the next timestamp.

```json
{
  "dna_damage.energy": [
    {"t": 0.0, "value": 0.05},
    {"t": 0.5, "value": 0.3},
    {"t": 1.0, "value": 0.0}
  ],
  "wound_field.signal": 0.1
}
```

Use `--input` for quick overrides; JSON is better for replaying experimental stimuli or scripted perturbations. CLI flags win if both specify the same `node.port`.

The repository ships `examples/live_demo_inputs.json` so you can feed the demo graph a tiny step stimulus:

```
helix live run \
  --hgx examples/live_demo.hgx \
  --inputs-json examples/live_demo_inputs.json \
  --duration 4
```

## Inspect bundles

Use `helix live inspect --bundle <path>` to summarize any stored run. The command prints the model metadata, dt overrides, island breakdown (with up to `--max-islands` rows), and head/tail snapshot summaries so you can sanity-check outputs without opening raw JSON.

```
# Human-readable summary
helix live inspect runs/egfr-demo/ --head 2 --metrics 2

# Automation-friendly JSON
helix live inspect runs/egfr-demo/ --json summary.json
# Filter to specific metrics for CI
helix live inspect runs/egfr-demo/ --json summary.json --metric pERK --metric cell_count
```

Trimmed JSON example (`schema_version` keeps the format stable for CI):

```json
{
  "schema_version": 1,
  "model": {
    "name": "EGFR Monolayer",
    "nodes": 3,
    "edges": 2
  },
  "runtime": {
    "islands": 2,
    "snapshots": 4,
    "wall_time_sec": 0.42
  },
  "nodes": {
    "wound_field": {
      "kind": "field",
      "metrics": {
        "tile": { "min": 0.0, "max": 0.19, "mean": 0.08, "var": 0.009, "std": 0.095 }
      }
    }
  }
}
```

Pair this with `jq`/`rg` when you need advanced slicing, but `inspect` should cover the quick “what happened?” loops whether you prefer text or JSON. Use `--no-metrics` to omit the per-node aggregates when inspecting very large bundles; the structural/runtime metadata stays available. `--metric` lets you keep only the metrics you care about (e.g., `pERK`, `cell_count`) so CI payloads stay tiny.

## Bundle layout

Every run emits a reproducible bundle:

| File | Description |
| --- | --- |
| `model.hgx` | Exact HGX payload that was executed. |
| `snapshots.json` | Time-ordered scheduler snapshots (per sync boundary). |
| `deltas.json` | StateReducer diffs suitable for UI delta-streaming. |
| `events.json` | EventBus stream (timestamps + snapshots). |
| `meta.json` | Seeds, cadence, dt overrides, CLI command, Helix version, timestamp. |

You can hand this bundle around, rerun with the same seed, and plug the JSON directly into dashboards or notebooks.

## Realtime viz + schema

Pass `--realtime` (optionally with `--realtime-endpoint tcp://0.0.0.0:7777`) to broadcast `helix.live.delta.v1` payloads while the scheduler runs. Each payload contains:

- `schema`: the literal string `helix.live.delta.v1`.
- `time` + `runtime`: the latest simulation time and metadata (slice, variant, hz, etc.).
- `snapshot` (first tick only) + `delta`: full node states and incremental diffs thereafter.
- `node_meta`: node kinds + port descriptions so visualizers can discover fields, agents, and observers without bespoke hooks.

Attach the new GUI with `helix live viz --endpoint tcp://127.0.0.1:8765` (or `--bundle runs/egfr_wt` to replay `deltas.json`). The renderer spins up a custom GLFW + moderngl HUD, paints the field heatmap, draws instanced cells with per-agent pERK colors, plots metrics, and exposes bidirectional controls:

- Pause/Resume toggles call into the active scheduler without blocking the sim loop.
- The GRB2 slider updates the SBML `egfr_rules.grb2_scale` parameter in realtime.
- Variant dropdowns dispatch `set_variant` hot swaps (EGFR/GRB2 updates SBML parameters on the fly).

Since the transport is socket-based you can run multiple viz clients (or headless loggers) in parallel. Commands use the same structured schema as the viewer (e.g., `{"kind": "live_control", "type": "pause"}` or `{"kind": "live_control", "type": "set_param", "target": "egfr_rules", "param": "grb2_scale", "value": 1.6}`) so custom dashboards can build on the same channel.

Need a quick feed without running a full model? `helix live demo --hz 30 --duration 120` spins up a synthetic EGF gradient + 200 mock agents that stream over the realtime endpoint, so you can test the renderer/controls end-to-end before wiring in a real slice.

## Plan runs without executing

Use `helix live plan` to dry-run a model/config and inspect the SCC→island breakdown, dt table, and overrides before launching a full simulation. This is perfect for PR reviews or checking `--dt node=...` overrides.

```
helix live plan --hgx examples/live_demo.hgx --dt cells=0.25

Live plan for examples/live_demo.hgx
  model: LiveDemo
  sync_dt: 0.25
  default_dt: 0.1
  dt overrides:
    - cells = 0.25
  graph: 3 nodes / 4 edges
  islands:
    1. island_0 (dt=0.1)
       nodes: egf_field, cells, observer
```

## LiveLab (interactive dev shell)

Need to sketch or tweak a graph without editing files? `helix live dev` drops you into an interactive REPL (LiveLab). You can add/delete nodes, connect ports, set per-node dt, run the scheduler, and export an `.hgx`/bundle — all while the sim stays hot.

```
# start empty
helix live dev

# start from an existing HGX
helix live dev --hgx examples/live_demo.hgx
```

Inside the shell, commands are prefixed with `:` (see `:help`):

```
:add node egf kind=Field diffusion=50.0 decay=0.01 mesh=[128,128,1.0]
:add node cells kind=ABM size=[64,64] max_agents=500
:connect egf.tile cells.field_signal
:plan
:run t=5 realtime
:export hgx builds/my_demo.hgx
:quit
```

The LiveLab session uses the same scheduler, delta stream, and export path as the standard CLI, so anything you design interactively can flow straight into `helix live run`/`viz`/`inspect`.

## Tips

- Use `examples/*.hgx` as starting points (EGFR monolayer, p53 oscillator, toggle switch).
- The node → dt map makes it easy to test stiff vs. slow islands before promoting to GPU kernels.
- `--input` is great for quick stimulus experiments (e.g., step response on `dna_damage.energy`).
- The bundle metadata already captures CLI arguments and Helix version, so you can treat each run as an auditable artifact. Add any additional annotations by writing sidecar JSON into the same directory.
