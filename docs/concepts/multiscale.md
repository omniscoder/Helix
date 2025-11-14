# Multiscale Stack: From Edit to Phenotype

This guide distills a practical pathway for simulating CRISPR edits all the way to tissue-level phenotypes. Treat it as the minimal kit needed for reproducible “edit → signaling → spatial context → cells → structure” studies.

## Layer stack

| Layer | Purpose | Primary tooling |
| --- | --- | --- |
| Rule-based kinetics | Encode post-edit signaling/regulation with combinatorial rigor. | PySB for authoring, BioNetGen + NFsim for ODE/SSA or network-free SSA when site combinatorics explode. |
| Spatial RD | Capture ligand gradients and microenvironment effects. | Smoldyn/ReaDDy 2 for particle-scale BD + reactions; Lattice Microbes (pyLM) for GPU RDME. |
| Cell dynamics | Model proliferation, death, mechanics, and microenvironments. | PhysiCell (agent-based, built-in transport), CompuCell3D (CPM + SBML per-cell), Morpheus (friendly GUI with native ODE/RD coupling). |
| Protein structure/interactions | Evaluate edits that alter binding interfaces/rates. | AlphaFold-Multimer/AF3 for complexes, HADDOCK/ClusPro for docking, OpenMM for short MD refinements. |
| Deterministic SBML engines | Fast gradients and reusable code. | AMICI (CVODES/IDAS) for compiled ODE/DAE; libRoadRunner/Tellurium for lightweight scripting. |
| Standards | Exchangeability and reproducibility. | SBML for models, SED-ML for simulation recipes, COMBINE archives + containers for packaging. |

## Coupling strategy

1. **Edit layer → Pathway model**: Genomic outcomes drive structural changes in rule-based models (remove rules, tweak rates, gate binding sites).
2. **Pathway ↔ Spatial RD**: Operator splitting per global Δt: run RD for `m` substeps (Smoldyn/ReaDDy/Lattice Microbes), then advance pathways using updated concentrations.
3. **Pathway → Cell ABM**: Map pathway observables (pERK, p53, NF-κB, etc.) to phenotypic knobs (division, apoptosis, motility, adhesion).
4. **Cell ABM ↔ RD**: Cells consume/secrete species into the diffusion solver; diffusion fields feed back on cell fate (hypoxia, paracrine cues). PhysiCell and CC3D expose this handshake out of the box.
5. **Protein hooks**: When edits touch interfaces, rebuild/dock complexes (AF-Multimer, HADDOCK/ClusPro) and, if needed, run short OpenMM MD to decide on kinetic parameter shifts (k_on / k_off) feeding the pathway.

## MVP vertical slice

**Objective**: EGFR–RAS–MAPK response after a GRB2 knockout in a 2D monolayer exposed to an EGF stripe.

- **Pathway**: PySB → BioNetGen EGFR→RAS→RAF→MEK→ERK with GRB2 adapter. Provide WT vs KO parameter sets (remove GRB2-mediated bindings). Use ODE mode first, flip to NFsim when multi-site phosphorylation appears.
- **RD**: Start with GPU RDME (Lattice Microbes) to diffuse/consume EGF; pivot to Smoldyn if particle realism near the front matters.
- **Cells**: PhysiCell (throughput) or CC3D (morphology). Let pERK modulate division and migration; low EGF or high p53 increases apoptosis.
- **Protein checkpoint**: If the edit clips an SH2 motif, run AF-Multimer for EGFR tail + GRB2 SH2, sanity-check via HADDOCK/ClusPro, and adjust binding rates accordingly.

## Orchestrator philosophy

Keep module interfaces skinny: one Python driver (“`simhub.py`”) steps the world clock, invokes each solver, and exchanges only the shared variables (local ligand values, handful of pathway observables, phenotype parameters). Each module remains swappable.

```python
class PathwayModel:
    def __init__(self, grb2_present: bool):
        self.state = {"pERK": 0.1}
        self.grb2_present = grb2_present

    def step(self, dt: float, ligand: float):
        gain = (1.0 if self.grb2_present else 0.2) * ligand / (0.5 + ligand)
        self.state["pERK"] += dt * (gain - 0.3 * self.state["pERK"])
        self.state["pERK"] = max(self.state["pERK"], 0.0)

    def phenotype(self):
        p = 0.01 + 0.09 * (self.state["pERK"] / (0.5 + self.state["pERK"]))
        return {"p_divide": p}
```

Swap that toy logic for AMICI/NFsim outputs when ready.

## Calibration, UQ, hygiene

- **Parameter fitting**: pyPESTO + AMICI for gradient-based fitting against phospho-proteomics / live-cell traces (supports replicates & qualitative constraints).
- **Sensitivity**: SALib (Sobol/Morris) via batched orchestrator runs to find the knobs that matter.
- **Provenance**: Store SBML + SED-ML in COMBINE archives, version datasets and containers, and persist seeds/configs.

## Kit list (opinionated)

| Domain | Picks |
| --- | --- |
| Kinetics | PySB, BioNetGen, NFsim, AMICI (for gradients). |
| Spatial | Lattice Microbes (GPU RDME) first; Smoldyn/ReaDDy when particle-level fidelity matters. |
| Cells | PhysiCell for throughput, CompuCell3D when morphology matters, Morpheus for GUI-centric coupling. |
| Protein | AlphaFold-Multimer → HADDOCK/ClusPro → OpenMM (only when kinetics change). |

## Immediate actions

1. Lock the EGFR/GRB2 MVP scope; seed the repo with the orchestrator scaffold, SBML stubs, and synthetic datasets.
2. Stand up AMICI + pyPESTO for fitting and SALib for sensitivity on the combined stack.
3. Pick PhysiCell (speed) or CC3D (shape) as the first ABM and produce an end-to-end run with spatial pERK and growth maps.

Only escalate fidelity (e.g., MD, particle RD) where it changes the story; multi-rate stepping and stable interfaces keep the system tractable.

### EGFR/GRB2 command cheat-sheet

```bash
# Build AMICI backends + calibrate parameters from config.yaml
python -m helixtasks.egfr_grb2_build_amici
python -m helixtasks.egfr_grb2_fit_pypesto

# Run WT vs KO (and GRB2-high) slices and capture bundles
helix live run models/egfr_grb2/config.yaml --variant wt --bundle runs/egfr_wt
helix live run models/egfr_grb2/config.yaml --variant ko --bundle runs/egfr_ko
helix live run models/egfr_grb2/config.yaml --variant hig --bundle runs/egfr_hig
# or use the slice registry shortcut
helix live run --slice egfr_grb2 --variant wt --bundle runs/egfr_wt
helix live run --slice egfr_grb2 --variant ko --bundle runs/egfr_ko
helix live run --slice egfr_grb2 --variant hig --bundle runs/egfr_hig

# Inspect metrics (pERK + cell counts) for CI or debugging
helix live inspect runs/egfr_wt --json --metric pERK --metric agents
helix live inspect runs/egfr_ko --json --metric pERK --metric agents

# Stream deltas for a realtime viewer (logs until a GUI attaches)
helix live run --slice egfr_grb2 --variant wt --realtime
```

The config-driven workflow keeps variants centralized—add new mutants by editing `models/egfr_grb2/config.yaml`, rebuild AMICI models, and the same `helix live run`/`inspect` commands continue to work. Continuous tests compare WT vs KO bundles via the JSON summary (pERK means/max, cell counts) so regressions surface immediately.

### Variant matrix (pERK & growth expectations)

| Variant | GRB2 level | Expected pERK | Growth phenotype |
| --- | --- | --- | --- |
| `ko` | 0× (GRB2 knockout) | Lowest (attenuated) | Slowest proliferation |
| `wt` | 1× (wild-type) | Mid-range | Baseline proliferation |
| `hig` | >1× (overexpression) | Highest (sustained) | Fastest proliferation |

CI enforces these relationships (`ko < wt ≤ hig`) for both pERK mean/max and mean/max cell counts using the LiveGraph bundles.

### Realtime streaming notes

- Passing `--realtime` to `helix live run` spawns a non-blocking event queue fed by the LiveScheduler’s `StateReducer`. Each sync boundary logs message stubs such as `[realtime:EGFR_GRB2/wt] t=12.0 +0 ~2 -0`, which already proves the delta feed is alive.
- The `helix live viz` client attaches to the same queue/socket, renders fields + agents at 60 FPS, and keeps the sim thread non-blocking.
- The custom Helix HUD (GLFW + moderngl) consumes the same delta schema (fields/agents/metrics) without any external UI toolkit.
- `helix live viz --endpoint tcp://127.0.0.1:8765 --hz 60` launches the HUD, decodes `helix.live.delta.v1` payloads, paints the EGF field heatmap, draws instanced cells (positions + pERK colors), and streams pERK/cell-count metrics. Use `--bundle runs/egfr_wt` to replay `bundle/deltas.json` offline.
- The renderer exposes a control dock: pause/resume toggles call back into the running scheduler, the GRB2 slider updates the SBML `egfr_rules.grb2_scale` parameter in realtime, and the variant dropdown fires `set_variant` hot swaps (EGFR/GRB2 maps these to SBML parameter sets so WT/KO/HIG reuse the same UI flow).
- Delta payloads now carry schema + node metadata (`helix.live.delta.v1`): downstream slices (p53, toggle switch, etc.) can rely on `payload["node_meta"][node]["kind"]` to discover fields/agents/observers without bespoke adapters. The transport also supports multi-process sockets (`--realtime-endpoint tcp://0.0.0.0:7777` or `ipc:///tmp/helix.sock`) so local or remote viz clients can subscribe simultaneously.
