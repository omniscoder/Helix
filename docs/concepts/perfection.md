# Perfect Helix Definition

This document pins down what “perfect” means for the Helix runtime, how we inspect it, and which slices we ship first. Treat it as the non‑negotiable checklist when proposing new work or deciding that a feature is “done”.

## Non-negotiables

### Scientific fidelity
- Reproduce canonical benchmarks: EGFR adaptation, p53–MDM2 oscillations, toggle-switch bistability, logistic tissue growth, and wound-closure kinetics.
- Parameterize against real assays with posterior predictive checks (PPCs) and held-out validation datasets.

### Numerical guarantees
- Typed ports carry explicit physical units, with strict mass and positivity invariants enforced at boundaries.
- Deterministic seeding plus bitwise-identical trajectories across machines and runs.
- Stable multi-rate stepping across solver islands.

### Performance envelope
- UI holds a steady 60 Hz loop with <10 ms jitter; delta-sync state delivery keeps renderer work bounded.
- GPU kernels run within sub-millisecond budgets using micro-batching and multi-stream execution.

### UX and reproducibility
- Every run emits a COMBINE-style bundle bundling IR, config, seeds, and outputs.
- Hot-swap edits (params, rules, or nodes) reset only the touched state.

## Loops That Keep Us Honest

1. **Model Loop** — refine the IR structure. Graph IR supports cycles; strongly connected components (SCCs) collapse into solver islands. Lift/restrict couplers enforce invariants at the boundaries.
2. **Data Loop** — fit and validate models. Combine gradient-based and sampling calibration, run PPCs, sensitivity analyses (Sobol/Morris), and use active learning to drive down uncertainty.
3. **Performance Loop** — guarantee real-time behavior. Autotune per-island Δt, reuse kernels via CUDA Graphs, stream deltas to the renderer, and track UI frame cadence continuously.

## Architecture Snapshot

- **IR**: general directed graph (not a DAG).
- **Runtime**: compute SCCs, build the condensation DAG, treat each SCC as a solver island.
- **Scheduler**: multi-rate per island with operator-split sync at DAG boundaries, content-addressed caching keyed by `(node hash, inputs, Δt)`, and a hot-swap path that minimizes state resets.
- **Live mode**: StateReducer emits deltas, renderer double-buffers, compute never waits on UI.

### Node palette

| Node | Role |
| --- | --- |
| `ProteinNode` | Patch-level MD/docking, ΔΔG to rate multipliers |
| `RuleNetNode` | ODE/SSA/rule-based hybrid |
| `GRNNode` | Gene regulation (ODE/SSA) |
| `FieldNode` | Reaction-diffusion/PDE tiles |
| `ABMNode` | Agent-based models at cell/agent level |
| `MechNode` | Tissue mechanics |
| `CouplerNode` | Lift/restrict, flux maps, mixing |
| `RewriterNode` | CRISPR/prime/PTM/pathway rewrites |
| `ObserverNode` | Metrics, losses, reports |

## Reference Code

```python
# core/scc.py
from collections import defaultdict

class Graph:
    def __init__(self): self.adj = defaultdict(set)
    def add_edge(self, u, v): self.adj[u].add(v)

def tarjan_scc(g: Graph):
    index, stack, onstack, idx, sccs = {}, [], set(), 0, []
    low = {}
    def strongconnect(v):
        nonlocal idx
        index[v] = low[v] = idx; idx += 1
        stack.append(v); onstack.add(v)
        for w in g.adj[v]:
            if w not in index:
                strongconnect(w); low[v] = min(low[v], low[w])
            elif w in onstack:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop(); onstack.remove(w)
                comp.append(w)
                if w == v: break
            sccs.append(tuple(comp))
    for v in list(g.adj.keys()) | {x for vs in g.adj.values() for x in vs}:
        if v not in index: strongconnect(v)
    return sccs

def condensation_dag(g: Graph):
    comps = tarjan_scc(g)
    comp_index = {v: i for i, c in enumerate(comps) for v in c}
    dag = Graph()
    for u in g.adj:
        for v in g.adj[u]:
            cu, cv = comp_index[u], comp_index[v]
            if cu != cv: dag.add_edge(cu, cv)
    return comps, dag
```

```python
# core/scheduler.py
from time import perf_counter

class Island:
    def __init__(self, name, nodes, dt):
        self.name, self.nodes, self.dt = name, nodes, dt

    def step(self, t, dt, inputs):
        # run ODE/SSA/PDE/ABM kernels; return outputs dict
        return {}

class LiveScheduler:
    def __init__(self, islands, sync_dt, state_reducer):
        self.islands = islands
        self.sync_dt = sync_dt
        self.state_reducer = state_reducer
        self.buffers = {isl.name: {} for isl in islands}

    def run_until(self, T, wall_budget_ms=5):
        t = 0.0
        next_sync = self.sync_dt
        wall_anchor = perf_counter()
        while t < T:
            for isl in self.islands:
                steps = int(self.sync_dt // isl.dt)
                for _ in range(max(1, steps)):
                    out = isl.step(t, isl.dt, self.buffers[isl.name])
                    self.buffers[isl.name] = out
                    t += isl.dt
            snapshot = self.collect_snapshot()
            self.state_reducer.push(snapshot)
            elapsed_ms = (perf_counter() - wall_anchor) * 1e3
            if elapsed_ms < wall_budget_ms:
                pass
            wall_anchor = perf_counter()

    def collect_snapshot(self):
        return {}
```

```python
# live/state_reducer.py
class StateReducer:
    def __init__(self, hz=60):
        self.hz = hz; self._last = None

    def push(self, snapshot):
        delta = diff(snapshot, self._last)
        self._last = snapshot
        self.emit(delta)

    def emit(self, delta): ...
```

## Validation Slices

| Slice | Stack | Targets |
| --- | --- | --- |
| EGFR–ERK with GRB2 edit | RuleNet + RD field + ABM | Dose–response curve, adaptation half-time, wound-closure t50 |
| p53–MDM2 oscillator | RuleNet + GRN | Period/amplitude distribution, noise-induced switching |
| Toggle switch (microfluidic gradient) | GRN + RD | Hysteresis loop, bifurcation map, front pinning |
| Matrix-sensing motility | ABM + Mech + RD | Speed vs stiffness curve, traction stress statistics |

Each slice is an end-to-end LiveGraph model with hot-swappable edits and observers.

### Validation matrix

| Slice | Metrics | Pass criteria |
| --- | --- | --- |
| EGFR | Dose–response, adaptation t½, spatial front speed | Within published CIs, monotone dose–response, adaptation overshoot below threshold |
| p53–MDM2 | Period, amplitude distribution, noise spectrum | Period coefficient of variation in expected band, autocorrelation peak at period, PPC covers held-out traces |
| Toggle switch | Bifurcation diagram, hysteresis area | Bistable region matches theory, measurable hysteresis loop |
| Motility | Speed–stiffness curve, traction stats | Speed rises to an optimum then drops, stress histogram shape matches reference |

## Data, Calibration, and Uncertainty

- Inputs: phospho-proteomics time courses, single-cell traces, imaging fields (AnnData/OME-TIFF), literature-derived parameter priors.
- Fitting: gradient (adjoint ODE where possible) plus stochastic (MCMC/VI) per island; combine outcomes via modular posteriors.
- Uncertainty: propagate distributions through couplers and surface credible bands live.

## Observability

- Per-node dashboards: solver Δt, SSA firing rate, PDE residual norm, agent counts, cache hit rate.
- Assertions: unit/mass/positivity, finite values, CFL checks. Fail fast with labeled errors.

## KPIs

| Category | Metric | Target |
| --- | --- | --- |
| Fit | ELPD / RMSE on held-out time series, posterior coverage | Coverage between 90–95% |
| Dynamics | Canonical time constants, oscillation periods, bifurcation points | Match within confidence bands |
| Invariants | Unit/mass/positivity checks | 100% pass rate across seeds |
| Real-time | UI cadence, frame latency, sim→UI bandwidth | 60 Hz UI, 95th percentile latency <16 ms, <2 MB/s for 100k agents |
| Repro | Bitwise reproducibility | Identical outputs for identical seeds/builds |

## Repo Scaffold

```
helix/
  core/
    graph.py            # IR: nodes/ports/edges (cycles allowed)
    scc.py              # Tarjan/Kosaraju + condensation DAG
    scheduler.py        # Multi-rate DAG-of-islands executor
    invariants.py       # mass/unit/positivity checks
    cache.py            # content-addressed memo
  nodes/
    rulenet.py
    grn.py
    rd.py
    abm.py
    mech.py
    protein.py
    couplers.py
    rewriter.py
    observer.py
  live/
    state_reducer.py
    event_bus.py
    hot_swap.py
  dsl/
    hgx.py
    schema.yaml
  stdlib/
    pathways/egfr_erk.bngl
    grn/p53_mdm2.yaml
    behaviors/epithelial.yaml
    observers/*.yaml
  examples/
    egfr_monolayer.hgx
    p53_oscillator.hgx
    toggle_switch.hgx
  tests/
    test_scc.py
    test_invariants.py
    test_scheduler_islands.py
    test_deltas.py
  cli.py
  README.md
```

## Immediate Build Targets

1. Implement IR + SCC condensation + LiveScheduler stubs above.
2. Ship the EGFR monolayer `.hgx` model with knockout and interface-mutation rewriters.
3. Add observers for adaptation score and wound-closure t50, complete with uncertainty bands.
4. Wire delta-sync capable of handling 100k agents plus a 512×512 RD tile map.
5. Treat the KPI table as the definition of done.
