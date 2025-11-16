# Helix: Computational Biology Playground + Digital Editor
[![PyPI](https://img.shields.io/pypi/v/veri-helix.svg)](https://pypi.org/project/veri-helix/)
[![Reproducible Viz (spec v1.0)](https://img.shields.io/badge/reproducible%20viz-spec%201.0-6f42c1)](docs/schema.md)
[![Benchmarks](https://img.shields.io/badge/bench-passing-success)](https://omniscoder.github.io/Helix/benchmarks/)
[![Verified by VeriBiota](https://img.shields.io/badge/verified%20by-VeriBiota-10b981?logo=lean&logoColor=white)](docs/veribiota.md)
[![VeriBiota CI](https://github.com/omniscoder/Helix/actions/workflows/veribiota.yml/badge.svg)](https://github.com/omniscoder/Helix/actions/workflows/veribiota.yml)

Helix is a hobbyist-friendly, simulation-only toolkit for computational genomics, CRISPR/Prime editing, PCR amplicon modeling, RNA folding, string algorithms, and reproducible visualizations. It’s the “notebook laboratory” that complements our production platform OGN: when a prototype proves itself here, it can graduate into OGN’s GPU-accelerated pipelines.

Helix aims to make digital biology fun: tiny datasets, clean Python APIs, batteries-included CLI commands, and a full provenance system so every plot, JSON, and DAG can be reproduced down to the SHA-256.

👉 Docs & Playground: https://omniscoder.github.io/Helix/

---

## Why Helix?

- **End-to-end edit DAGs** – CRISPR, Prime, and PCR simulations all produce explicit edit DAGs with per-node genome materialization, log-probabilities, and provenance. The same artifacts drive notebooks, PNGs, and interactive web/desktop viz.
- **Blueprint-level provenance** – every DAG, PNG, and `.viz.json` is tagged with schema kind, spec version, SHA-256, and runtime metadata. You can always answer *“where did this figure come from?”*.
- **Realtime + batch in one toolkit** – run tiny “bench” genomes in a notebook, or stream full edit DAGs into a GUI / Playground via JSONL frames for live inspection.
- **Formal verification hooks** – optional VeriBiota integration re-encodes Helix DAGs in Lean and proves structural + probabilistic invariants. If the badge is green, the math checked out.
- **Teaching-friendly, research-grade** – small, inspectable examples and clear CLI affordances, but under the hood you get FM-indexes, De Bruijn graphs, MinHash/HLL, prime editing physics, and CRISPR scoring engines.

---

## Verified by VeriBiota™

Helix’s CRISPR, Prime, and PCR DAGs carry the **Verified by VeriBiota** badge when every generated artifact satisfies:

- **Structural correctness:** unique root, acyclic topology, monotonic depth, and rule-consistent transitions.
- **Semantic correctness:** every edit event matches the formal CRISPR/Prime/PCR rewrite semantics tracked in Lean.
- **Probability sanity:** outgoing probabilities and terminal probabilities sum to ≈1 with no negative or undefined weights.
- **Reproducibility:** deterministic replays under recorded seeds, provenance-aligned frame streams, and matching sequence hashes.

The contract is documented in [docs/veribiota.md](docs/veribiota.md) and enforced by the `helix veribiota lean-check`, `preflight`, and `export-dags` commands + the `veribiota.yml` GitHub Action badge above. If you see the badge, the math has run.

## What Helix Is

- Think of Helix as:
- A computational wet lab you can run in a notebook.
- A CRISPR/Prime/PCR digital twin simulator using deterministic rules and DAGs.
- A teaching toolkit for genomics algorithms and bioinformatics basics.
- A bridge into OGN, sharing concepts and API style while staying lightweight and experimental.

Everything Helix emits is a **software artifact**: JSON, PNG, YAML, GraphML, viz-specs, or CLI summaries. No wet-lab instructions. No reagent guidance.

## Highlights

### Core genome + sequence tools
- DNA summaries, k-mer counting, GC skew, motif clustering
- FM-index search + approximate search via Myers bit-vector
- Minimizers, syncmers, and seed-extend demos
- ORF detection, frameshift heuristics, translation, protein metrics

### RNA folding + ensembles
- Zuker-style MFE
- McCaskill partition function
- Ensembles, dot-plots, entropy tracks, centroids structures

### CRISPR simulation stack
- Guide discovery with PAM registries
- Off-target scanning (exact + mismatch-tolerant)
- Probabilistic CRISPR cut/repair simulation
- Edit DAGs (helix.crispr.edit_dag.v1.1) with node-level genome snapshots
- Multi-guide, region-specific, FASTA-native support
- Real-time JSONL frames for animation + Playground sync
- Protein-impact annotation via transcript models

### Prime editing simulator
- pegRNA + Prime Editor definitions
- RTT/RTT-branching logic
- Outcome probabilities
- Full Prime Edit DAGs with per-node genome materialization
- Multi-peg batch workflows for large design sets
- PCR Amplicon DAGs

### PCR amplicon DAGs
- Primer binding + per-cycle branching
- Amplicon growth simulation
- Visualization + JSON artifact export

### Graphs & sketching
- Build/clean/color De Bruijin graphs
- MinHash/HLL Sketching for fast genome distance

### Reproducible visualization
- Viz-spec system (structured spec, SHA-verified)
- Every PNG ships a provenance sidebar + sibling .viz.json
- CLI + notebook-friendly

### Workflows & experiments
- YAML-driven experiments: CRISPR, Prime, PCR
- Generates DAGs, PNGs, reports, provenance manifests
- Reproducible from a single .helix.yml

### GUI + realtime viz
- Optional PySide6 desktop shell
- Cytoscape/D3 visualization of live edit DAG streams
- Offline-capable; accepts CLI artifacts and real-time frames

## Quick Start
### Crispr Edit DAG
```
helix crispr dag \
  --genome genome.fna \
  --guide-sequence GGGGTTTAGAGCTATGCT \
  --json crispr_dag.json
```

### Prime Editing Simulator
```
helix prime dag \
  --genome genome.fna \
  --peg-config peg.json \
  --editor-config pe3.json \
  --json prime_dag.json
```

### PCR Amplicon DAG
```
helix pcr dag \
  --genome plasmid.fna \
  --primer-config primers.json \
  --pcr-config pcr_settings.json \
  --out pcr_dag.json
```
### RNA Folding
```
helix rna mfe --fasta hairpin.fna --dotbracket mfe.dbn
helix rna ensemble --fasta hairpin.fna --gamma 1.0 --dotplot plot.png
```

### Sequence Tools
```
helix dna --input seq.fna --window 400 --k 5 --plot-skew
helix seed index seq.fna --method minimizer --k 15 --plot seeds.png
helix string search seqs.fna --pattern GATTACA --k 1
```

## Reproducible Viz (Spec 1.x)
Every plot:
- stamps a provenance footer (Helix version, SHA-256, viz kind)
- emits a .viz.json viz-spec
- stores full command replay in <image>.provenance.json
Use:
```
helix viz --schema
helix schema manifest
helix schema diff old.json new.json
```

## Workflows
Turn a single YAML file into a full simulation:
```
helix experiment new --type crispr --out demo.helix.yml
helix experiment run --config demo.helix.yml --out dag.json
helix experiment viz --config demo.helix.yml --out dag.png
```
Everything—FASTA, guides, Cas config, seed, viz parameters—is captured and reproducible.

## Repo Layout
```
src/helix/
    crispr/         # CRISPR models, guides, simulators, edit DAGs
    prime/          # Prime editing engine + DAGs
    pcr/            # PCR amplicon simulator
    rna/            # MFE + ensemble folding
    bioinformatics/ # DNA utilities, k-mers, motif clustering
    seed/           # minimizers/syncmers + mapping demo
    string/         # FM-index + Myers ED
    graphs/         # DBG tooling
    viz/            # spec-driven visualization system
    workflows/      # YAML runner
    gui/            # optional PySide6 desktop app
benchmarks/
docs/
examples/
tests/
```

## Getting Started
### Requirements
- Python 3.10+ (3.11 tested)
- pip or another package manager
- Optional extras: `matplotlib` for plotting, `biopython` for protein helpers, `pyyaml` for workflow configs (already included in base deps).

### Installation
Stable release from PyPI (installs CLI + package):
```bash
python -m venv .venv
source .venv/bin/activate
pip install "veri-helix[viz,protein,schema]"
```
Need only the core library? Drop the extras (viz/matplotlib, protein/Biopython, schema/pydantic). For local development, clone the repo and run:
```bash
pip install -e ".[dev]"
```
This exposes the `helix` console command and the `helix` Python package (`from helix import bioinformatics`).

### Run a Script
- **K-mer + skew analysis**
  ```bash
  helix dna --input path/to/sequence.fna --window 400 --step 50 --k 5 --plot-skew
  ```
  Change the GC window/step, filter top k-mers, or point at the bundled dataset `src/helix/datasets/dna/plasmid_demo.fna`.
  For quick clustering with exports, try `python examples/kmer_counter.py --max-diff 1 --csv clusters.csv --plot-top 10`.

- **Neural net demo**
  ```bash
  python ann.py
  ```
  Prints training progress and final weights for a tiny XOR-style problem.

- **Translate a sequence**
  ```bash
  python examples/translate_sequence.py AUGGCCUUU
  ```
  Add `--no-stop` to continue through stop codons or point to a file with `--input`.

- **Find ORFs**
  ```bash
  python examples/find_orfs.py --min-length 90 --include-partial --detect-frameshifts --input your_sequence.fna --orf-fasta peptides.faa --orf-csv orfs.csv --frameshift-csv shifts.csv
  ```
  Prints coordinates, frames, strands, optional frameshift candidates, and can export FASTA/CSV artifacts.

- **Cyclo-spectrum playground**
  ```bash
  python examples/cyclospectrum_demo.py --peptide NQEL --spectrum "0,113,114,128,227,242,242,355,356,370,371,484"
  ```
  Print linear/cyclic spectra, score against an experiment, or recover candidate peptides with the leaderboard search.

- **RNA folding trace**
  ```bash
  python examples/nussinov_trace.py --input hairpin.fasta --min-loop 4
  ```
  Outputs the dot-bracket structure, base-pair list, and optional file export using the upgraded Nussinov implementation.

- **Protein summary**
  ```bash
  helix protein --input src/helix/datasets/protein/demo_protein.faa --window 11 --top 8
  ```
  Computes molecular weight, charge, hydropathy windows, and more (requires the `protein` extra / Biopython).

- **Unified Helix CLI**
  ```bash
  helix dna --sequence ACGTACGT --k 4
  helix spectrum --peptide NQEL --spectrum "0,113,114,128,227,242,242,355,356,370,371,484"
  helix rna mfe --fasta src/helix/datasets/dna/plasmid_demo.fna --dotbracket mfe.dbn
  helix rna ensemble --fasta src/helix/datasets/dna/plasmid_demo.fna --gamma 1.0 --dotplot dotplot.png --entropy entropy.png
  ```
  The `helix` entry point wraps the DNA, spectrum, RNA, protein, triage, viz, and workflow helpers so you can run ad-hoc analyses without hunting for scripts.

- **CRISPR guide + off-target scan**
  ```bash
  helix crispr find-guides --fasta target.fna --pam SpCas9-NGG --guide-len 20 --json guides.json
  helix crispr offtargets --fasta genome.fna --guides guides.json --max-mm 3 --json hits.json
  helix crispr score --guides guides.json --hits hits.json --weights weights/cfd-lite.json --json scores.json
  helix crispr simulate --fasta target.fna --guides guides.json --guide-id g1 --draws 1000 --seed 42 --json crispr_sim.json
  helix viz crispr-track --input crispr_sim.json --save crispr_track.png
  ```
  Produces schema-tagged JSON (`crispr.guides`, `crispr.offtargets`, `crispr.sim`) with optional scoring and cut/repair simulations; CLI viz renders a provenance-stamped PNG. Sequences remain masked unless `--emit-sequences` is explicitly passed.

- **CRISPR genome simulation**
  ```bash
  helix crispr genome-sim --genome genome.fna --guide-sequence GGGGTTTAGAGCTATGCT --cas cas9 --json crispr_cut_events.json
  ```
  Loads the genome into a `DigitalGenome`, instantiates a preset (or JSON-defined) `CasSystem`, and calls the in-silico cut simulator so you can inspect potential target sites. Outputs include serialized guides, Cas parameters, and any simulated `CutEvent` entries.

- **CRISPR edit DAG**
  ```bash
  helix crispr dag --genome genome.fna --guide-sequence GGGGTTTAGAGCTATGCT --max-depth 1 --json crispr_edit_dag.json
  ```
  Builds the first-version “digital twin” graph using the new edit runtime. Nodes contain fully materialized genome views, edges capture each clean-cut event, and the JSON artifact (`helix.crispr.edit_dag.v1.1`) can feed notebooks and future viz surfaces.

- **CRISPR edit DAG (FASTA-native, multi-guide)**
  ```bash
  helix crispr dag \
    --genome examples/hg19_chr_demo.fa \
    --cas-config examples/cas9.json \
    --guides-file examples/guides.tsv \
    --region chr7:55000000-55005000 \
    --max-depth 2 \
    --min-prob 1e-4 \
    --max-sites 20 \
    --seed 0 \
    --out-dir out/crispr_dags/
  ```
  - `cas9.json` (Cas/physics config)  
    ```json
    {
      "name": "SpCas9",
      "system_type": "cas9",
      "pam_pattern": "NGG",
      "cut_offset": 3,
      "max_mismatches": 3,
      "weight_mismatch_penalty": 1.0,
      "weight_pam_penalty": 2.0
    }
    ```
  - `guides.tsv` (guide library)  
    ```text
    name	sequence	                region
    G1	  ACGTACGTACGTACGTACGT	    chr7:55000010-55000040
    G2	  TGCATGCATGCATGCATGCA	    chr7:55003000-55003025
    G3	  CCCCCGGGGGAAAAATTTTT	    .
    ```
  The CLI slices the FASTA per guide (falling back to `--region`), fans out over each design, and writes one artifact per guide (e.g., `out/crispr_dags/crispr_001_G1.edit_dag.json`). Artifacts stay compatible with the viz/report/Playground surfaces.

- **Protein-impact annotations**
  ```bash
  helix crispr dag \
    --genome examples/hg19_chr_demo.fa \
    --guide-sequence ACGTACGTACGTACGTACGT \
    --coding-json transcripts/BRCA1_tx.json \
    --coding-transcript BRCA1-201 \
    --out out/crispr_brca1.edit_dag.json
  ```
  Provide a transcript JSON (same schema used by the GUI loader) to annotate SNV outcomes with `protein_impact` metadata (`silent`, `missense`, `nonsense`). Works for both `helix crispr dag` and `helix prime dag`.

- **Prime editing sandbox**
  ```bash
  helix prime simulate --genome genome.fna --peg-config peg.json --editor-config pe3.json --max-outcomes 16 --json prime_edits.json
  ```
  Wraps the new prime-editing models: pegRNA definitions (inline flags or JSON), prime-editor parameters, and the `simulate_prime_edit` entrypoint. Like the CRISPR command, this is purely computational—it emits hypothetical outcomes for downstream notebooks and viz.

- **Prime edit DAG**
  ```bash
  helix prime dag --genome genome.fna --peg-config peg.json --editor-config pe3.json --json prime_edit_dag.json
  ```
  Produces the prime-editing DAG artifact (`helix.prime.edit_dag.v1.1`) so you can inspect RTT-driven branches, log probabilities, and materialized genome snapshots per node.

- **Prime edit DAG (config-driven, multi-peg)**
  ```bash
  helix prime dag \
    --genome examples/hg19_chr_demo.fa \
    --editor-config examples/prime_editor.json \
    --pegs-file examples/pegs.tsv \
    --region chr11:5227000-5227500 \
    --max-depth 3 \
    --min-prob 1e-4 \
    --seed 0 \
    --out-dir out/prime_dags/
  ```
  - `prime_editor.json`
    ```json
    {
      "name": "PE2-like",
      "cas": {
        "name": "SpCas9-H840A",
        "type": "cas9",
        "pam_pattern": "NGG",
        "cut_offset": 3,
        "max_mismatches": 2
      },
      "nick_to_rtt_offset": 0,
      "efficiency_scale": 0.6,
      "mismatch_tolerance": 2,
      "indel_bias": 0.1,
      "metadata": {
        "flap_model": "left>right (demo)"
      }
    }
    ```
  - `pegs.tsv`
    ```text
    name	spacer	                     pbs	      rtt	                        region
    peg1	ACGTACGTACGTACGTACGT	     GCTAGCTA	  TCTGACTCTCTCAGGAGTC	     chr11:5227000-5227100
    peg2	TGCATGCATGCATGCATGCA	     AACCGGTT	  AAGGTTCCGGAACTTG	         chr11:5227200-5227300
    ```
  Each peg row emits a dedicated `helix.prime.edit_dag.v1.1` artifact, mirroring the CRISPR workflow. The Playground buttons (`?demo=prime`) load the same format, so teams can plug their configs directly into visualization/reporting pipelines.

- **Experiment configs (YAML → DAG)**
  ```bash
  helix experiment new --type crispr --out experiments/demo_crispr.helix.yml
  helix experiment run --config experiments/demo_crispr.helix.yml --out out/demo_crispr.edit_dag.json
  helix experiment viz --config experiments/demo_crispr.helix.yml --out out/demo_crispr.png
  ```
  A `*.helix.yml` file captures everything humans care about—FASTA path, optional region, Cas/Prime config, guide or peg design, and simulation knobs. Helix regenerates DAG JSON, PNGs, and HTML reports from that single spec (`helix experiment run/viz/report`). Starter templates live in `templates/`, and `helix experiment new` bootstraps a fresh config with placeholders ready to fill in.

- **Real-time DAG frames (JSONL)**
  ```bash
  helix crispr dag \
    --genome examples/hg19_chr_demo.fa \
    --cas-config examples/cas9.json \
    --guide-sequence ACGTACGTACGTACGTACGT \
    --frames - \
    --out out/crispr_rt.edit_dag.json
  ```
  Streams `helix.edit_dag.frame.v1` JSON lines so you can animate CRISPR edits as they unfold. See [Edit DAG Frames](docs/edit_dag_frames.md) for schema details, or try it live in the [Realtime Playground](docs/playground/realtime.html).

- **CRISPR DAG micro-verification**
  ```bash
  python benchmarks/verify_crispr_micro.py
  pytest tests/test_crispr_dag_micro.py
  ```
  A synthetic `tests/data/crispr_micro.fna` genome (two short chromosomes packed with overlapping NGG PAMs) powers a brute-force verifier that rebuilds the entire probability tree outside of the CRISPR physics engine. The `benchmarks/verify_crispr_micro.py` script emits a pass/fail summary, while `tests/test_crispr_dag_micro.py` runs the same cross-check during CI so we know every DAG leaf and probability mass matches the reference enumeration.
- **Lean/VeriBiota bridge**
  ```bash
  helix crispr dag --genome genome.fna --guide-sequence ACGT... --json out/eg.edit_dag.json
  helix veribiota export \
    --input out/eg.edit_dag.json \
    --out out/eg.lean \
    --dag-name exampleDag \
    --module-name VeriBiota.Bridge
  ```
  Converts any Helix `helix.*edit_dag.v1.*` artifact into a Lean module that defines `exampleDag : EditDAG`, emits the node/edge lists, and inserts a ready-to-run `#eval VeriBiota.check exampleDag` plus a theorem stub (disable with `--skip-theorem` / `--skip-eval`).

  To consolidate several JSON artifacts into one Lean namespace (faster CI, shared proofs):
  ```bash
  helix veribiota export-dags \
    --inputs out/dag1.json out/dag2.json \
    --module-name Helix.CrisprExamples \
    --list-name exampleDags \
    --out veribiota/generated/Helix/CrisprExamples.lean
  ```
  The generated module defines `def dag1 : EditDAG`, `dag2`, bundles them into `def exampleDags : List EditDAG`, and inserts an aggregate theorem stub (`∀ dag ∈ exampleDags, VeriBiota.check dag`) you can turn into a real proof.

- **Lean pipeline glue**
  ```bash
  helix veribiota lean-check --input out/dag1.json --out out/dag1.lean-check.json
  helix veribiota preflight --checks out/*.lean-check.json
  helix veribiota export-dags --inputs out/dag*.json --out veribiota/generated/Helix/CrisprExamples.lean
  ```
  Each DAG JSON gets a companion `.lean-check.json` (hashes, probabilities, metadata). `preflight` validates those summaries (and, optionally, re-hashes the source DAGs) so CI fails fast before Lean boots. Once preflight passes, `export-dags` emits a single Lean module for all DAGs, letting VeriBiota prove shared invariants (`well_formed`, probability sums, hash consistency) in one place.

  When integrating with the external [VeriBiota/VeriBiota](https://github.com/VeriBiota/VeriBiota) repo, point Helix directly at that checkout and let it populate the generated namespace:
  ```bash
  helix veribiota export-suite \
    --inputs out/dag*.json \
    --veribiota-root ../VeriBiota \
    --module-path Biosim/VeriBiota/Helix/MicroSuite.lean \
    --module-name Biosim.VeriBiota.Helix.MicroSuite
  ```
  This writes the Lean file straight into `Biosim/VeriBiota/Helix/MicroSuite.lean`, ready for VeriBiota’s Lake build + proof pipeline.

- **Frames → dataset**
  ```bash
  helix edit-dag generate-dataset --n 0 --frames-input run.frames.jsonl --out dataset.jsonl
  ```
  Converts a JSONL frame stream into a dataset record (mix with random generations by combining `--frames-input` and `--n`). Each row stays human-readable so you can audit what landed in the corpus:
  ```json
  {
    "id": 7,
    "mechanism": "crispr",
    "node_count": 11,
    "edge_count": 10,
    "top_outcomes": [
      {"stage": "repaired", "prob": 0.61, "sequence_hash": "9f4b1cbe"},
      {"stage": "error", "prob": 0.31, "sequence_hash": "72acdc01"}
    ],
    "frame_source": "runs/hbb_demo.frames.jsonl",
    "artifact": { "...": "helix.crispr.edit_dag.v1.1 payload" }
  }
  ```
- **Hero comparison (HBB vs mutant)**
  1. Open the [Realtime Playground](docs/playground/realtime.html).
  2. Guide A: `ACCCAGGAAACCCGGGTTTT`, Guide B: `TTTACCCAGGAAACCCGGGT`, PAM `NGG`.
  3. Hit “Compare” to watch probability mass shift between intended vs indel branches, then export the experiment `.helix.yml` for reproducible CLI runs.
- **Desktop GUI (optional PySide6 extra)**
  ```bash
  pip install veri-helix[gui]
  helix gui
  ```
  Ships a PySide6 desktop shell with a QWebEngineView + Cytoscape canvas. The GUI streams the same JSONL frames as the CLI/Playground, so you can iterate on CRISPR or Prime runs locally (even offline) and export specs later.

- **PCR amplicon DAG**
  ```bash
  python -m helix.cli pcr dag \
    --genome src/helix/datasets/dna/plasmid_demo.fna \
    --primer-config examples/pcr_primers.json \
    --pcr-config examples/pcr_config.json \
    --out pcr_amplicon_dag.json
  python -m helix.cli edit-dag viz --input pcr_amplicon_dag.json --out pcr_dag.png
  ```
  Simulates in-silico amplification (binding → cycles → error branches) and emits `helix.pcr.amplicon_dag.v1`. Drag-drop the JSON into the [Playground](https://omniscoder.github.io/Helix/playground/?json=assets/viz/pcr_amplicon_dag.json) for an interactive tour or animate it via `helix edit-dag animate`.

- **Edit DAG visualization**
  ```bash
  helix edit-dag viz --input examples/crispr_edit_dag.json --out crispr_dag.png
  helix edit-dag viz --input examples/prime_edit_dag.json --out prime_dag.png
  ```
  Renders any DAG artifact to a PNG using the built-in networkx/matplotlib helper. The `examples/` directory ships ready-to-plot JSON fixtures for CRISPR and Prime editing so you can kick the tires immediately.
  For an interactive walkthrough (including sequence diffs), open `docs/notebooks/edit_dag_visual_demo.ipynb`.

- **JSON configs → CLI**
  ```bash
  # Convert the JSON genome to FASTA on the fly
  python - <<'PY' > /tmp/demo_genome.fna
  import json
  cfg = json.load(open("examples/crispr_demo_genome.json"))
  for chrom in cfg["chromosomes"]:
      print(f">{chrom['name']}\\n{chrom['sequence']}")
  PY
  GUIDE=$(jq -r '.sequence' examples/crispr_demo_guide.json)
  helix crispr dag --genome /tmp/demo_genome.fna --guide-sequence "$GUIDE" --json /tmp/demo_dag.json
  ```
  The `examples/crispr_demo_genome.json` + `examples/crispr_demo_guide.json` pair provides a self-contained, copy-pastable config for CLI experiments without needing any external FASTA files.

- **Prime config quickstart**
  ```bash
  python examples/scripts/make_prime_demo_fasta.py --input examples/prime_demo_genome.json --out /tmp/prime_demo.fna
  helix prime dag --genome /tmp/prime_demo.fna \
    --peg-config examples/prime_demo_configs.json \
    --editor-config examples/prime_demo_configs.json \
    --json /tmp/prime_dag.json
  ```
  This uses the bundled `examples/prime_demo_genome.json` plus peg/editor definitions in `examples/prime_demo_configs.json` to build a full prime-edit DAG without any external data.

- **Workflow runner**
  ```bash
  helix workflows --config workflows/plasmid_screen.yaml --output-dir workflow_runs
  ```
  Chains multiple subcommands from YAML, captures per-step logs, and writes artifacts to structured run directories.

- **Visualization helpers**
  ```bash
  helix viz triage --json triage.json --output triage.png
  helix viz hydropathy --input src/helix/datasets/protein/demo_protein.faa --window 11
  ```
  Render plots directly from CLI artifacts (triage JSON, hydropathy windows). Requires matplotlib; hydropathy also needs Biopython.

- **Python API demo**
  ```bash
  python examples/helix_api_demo.py
  ```
  Showcases the `helix_api` module for notebook-friendly access to DNA summaries, triage reports, spectra, RNA folding, and (optionally) protein metrics. For full signatures and payload descriptions, see the [API reference](docs/reference/api.md).

- **Triage report CLI**
  ```bash
  python examples/triage_report.py --input your_sequence.fna --output triage.png --clusters-csv clusters.csv --orfs-csv orfs.csv
  ```
  Generates a composite plot plus optional CSV/FASTA exports for quick daily snapshots.

- **Notebook triage dashboard**
  Open `notebooks/triage_dashboard.ipynb` to plot GC skew, ORFs, and k-mer hotspots side-by-side for a quick daily scan.

- **Protein sequence peek**
  ```python
  from protein import show_sequence
  show_sequence("1CRN.cif")
  ```
  Requires the target structure file in the working directory (or adjust the loader).

Browse task-specific quickstarts in `examples/README.md`. Tiny datasets ship inside the package (see `helix.datasets.available()`), including `dna/human.txt`, `dna/plasmid_demo.fna`, and `protein/demo_protein.faa` for quick experiments with pandas, sklearn, or hydropathy charts.

### Run Tests
```bash
pytest
```
Pytest powers translator and k-mer regression checks; feel free to add more as you create new helpers.

### Benchmarks
```bash
python -m benchmarks.api_benchmarks --repeat 5 --warmup 1 --limit 0 \
  --out bench-results/api.json --summary-md bench-results/api.md
```
The benchmark harness now emits a schema-stamped payload (`bench_result` v1.0) that records commit SHA, dataset provenance, BLAS vendor, CPU/threads, locale, RNG seed, and per-case timing/RSS stats. Use `--scenario dna_summary` to focus on a subset, `--limit 10000` to mimic CI’s faster sweep, and `--summary-md` to capture the Markdown table that CI publishes automatically.

- **CI drift tracking**: the `benchmarks` GitHub Actions job pins `OMP_NUM_THREADS`/`MKL_NUM_THREADS`, seeds RNGs, runs the suite (`repeat=3` by default, `repeat=10` when `bench_heavy=true` via `workflow_dispatch`), uploads `benchmarks/out/bench-$GITHUB_SHA.{json,md}`, and appends the rendered Markdown summary to the workflow summary tab.
- **Regression gate**: `scripts/bench_check.py .bench/baseline.json benchmarks/out/latest.json --threshold 5` enforces a >+5 % slowdown limit and fails the workflow when hit. Update `.bench/baseline.json` whenever you intentionally change performance characteristics.
- **Heavier datasets**: export `HELIX_BENCH_DNA_FASTA` / `HELIX_BENCH_PROTEIN_FASTA` (or pass them as workflow_dispatch inputs) to stress-test larger references. The benchmark JSON records the absolute paths and sizes so dashboards can keep apples-to-apples comparisons. Store private or future references under `benchmarks/data/` and toggle them via `bench_heavy=true` when ready.
- **Dashboard**: CI appends every main-branch run to `docs/data/bench/history.csv`; the published chart lives at [docs/benchmarks.md](docs/benchmarks.md).

## Reproducible Viz & Viz-Spec
- Every `helix viz ...` (and CLI modes that call them) accepts `--save out.png` (PNG/SVG/PDF) and auto-emits a sibling `.viz.json` unless `--save-viz-spec` overrides the path.
- Each plot footer stamps `Helix vX.Y • viz-kind • spec=1.x • key params • timestamp • input_sha256` so shared figures always carry their provenance and the SHA-256 of the original JSON payload.
- The viz-spec JSON captures counts, quantiles, bounds, and the `input_sha256` used for hashing; regressions assert against that structured payload instead of brittle pixel hashes.
- You can feed those viz-specs (plus the original JSON inputs) into docs/notebooks to explain how a figure was produced and which parameters generated it.
- Explore or inspect schemas with `helix viz --schema`, diff manifests with `helix schema diff --base old.json`, export everything via `helix schema manifest --out schemas.json`, or render ready-to-plot payloads via `helix demo viz`.
- Workflows can enforce schemas per step and print provenance tables/JSON with `helix workflows ... --with-schema [--as-json]`.
- Every saved plot writes `<image>.provenance.json` next to the PNG, capturing `{schema_kind, spec_version, input_sha256, viz_spec_sha256, image_sha256, helix_version, command}` for chain-of-custody.
- Full schemas, screenshots, and sample payloads live under [docs/viz.md](docs/viz.md) and the [Schema Reference](docs/schema.md).

## Weekend Project Ideas
- Plot the GC skew for a bacterial plasmid and compare predicted origins to literature.
- Extend the ORF scanner to sweep reverse complements and test on viral genomes.
- Compare frameshift candidates against known gene models to flag likely sequencing errors.
- Pair the ORF scanner with the GC skew plot to compare predicted origins and coding regions.
- Use the CSV/plot outputs from `examples/kmer_counter.py` to highlight SNP hotspots and share charts with the community.
- Customize `notebooks/triage_dashboard.ipynb` with your own sequences and publish the visuals for project updates (digital reports only).
- Hook `cyclospectrum.py` into a simple leaderboard scorer and visualize the mass differences.
- Swap the activation function in `ann.py`, log loss curves, and document what changes.
- Build a notebook that fetches a PDB entry, prints its sequence via `protein.py`, and sketches the secondary structure counts.
- Chain `examples/translate_sequence.py` with `peptide_mass_lookup.py` to score translated open reading frames.

Browse ready-to-run snippets in `examples/README.md`, and share your results in `examples/` (add new files freely) or link to gist/notebook URLs in issues so others can remix.

## Design Philosophy
- **Approachable first**: readable code, inline comments when helpful, datasets that fit in memory.
- **Composable**: functions return plain Python data structures so you can plug them into pandas, NumPy, or future OGN pipelines.
- **Biopython-friendly**: we stand on Biopython's shoulders; no wheel reinvention when a stable API exists.
- **Prototype-to-production bridge**: helper scripts should make it easy to migrate successful ideas into OGN when the time comes.

## Roadmap
- Real-time 2.5D simulation panels (CRISPR + Prime + PCR)
- Unified DAG viewer with edit-diff timelines
- Notebook-to-OGN adapters for migrating prototypes
- More ensemble RNA metrics + stochastic simulators
- Expanded motif discovery solvers + GPU-ready variants
- Community-driven examples gallery

## Contributing
We welcome ideas, experiments, and docs improvements. To keep things playful:
- Open issues with context, references, or notebooks that inspired your idea.
- Tag contributions by complexity (`good-first-experiment`, `deep-dive`, etc.).
- Respect the code of conduct (be kind, give credit, document assumptions).
- If you plan a larger refactor, start a discussion thread so we can pair-program or offer pointers.

Happy hacking!
- **String search**
  ```bash
  helix string search sequences.fna --pattern GATTACA --k 1 --json hits.json
  ```
  Uses the FM-index for exact matches (`k=0`) or Myers bit-vector streaming for ≤k edit-distance hits in FASTA/plaintext inputs.
- **Seed + extend demo**
  ```bash
  helix seed index src/helix/datasets/dna/plasmid_demo.fna --method minimizer --k 15 --window 10 --plot seeds.png
  helix seed map --ref src/helix/datasets/dna/plasmid_demo.fna --reads src/helix/datasets/dna/plasmid_demo.fna --k 15 --window 10 --band 64 --xdrop 10
  ```
  Generates deterministic minimizers (or syncmers) and a simple seed-and-extend JSON summary; `--plot` uses `helix.viz.seed` for density snapshots.

- **DBG toolbox**
  ```bash
  helix dbg build --reads reads1.fna reads2.fna --k 31 --graph dbg.json --graphml dbg.graphml
  helix dbg clean --graph dbg.json --out dbg_clean.json
  helix dbg color --reads sample1.fna sample2.fna --labels case control --k 31 --out colored.json
  ```
  Builds/cleans JSON + GraphML de Bruijn graphs and produces colored DBG presence tables ready for pseudoalignment experiments.

- **Motif discovery**
  ```bash
  helix motif find --fasta promoters.fasta --width 8 --solver steme --iterations 40 --json motif.json --plot pwm.png
  ```
  Runs EM/STEME/online solvers to infer PWMs/log-likelihoods and renders optional probability heatmaps.
- **Sketching (MinHash/HLL)**
  ```bash
  helix sketch build --method minhash --fasta seq.fna --k 21 --size 1000
  helix sketch compare --method hll --fasta-a a.fna --fasta-b b.fna --precision 12
  ```
  Quickly approximate genome distances via Mash-style MinHash or HLL cardinality/Jaccard estimates.
