# Helix
[![PyPI](https://img.shields.io/pypi/v/veri-helix.svg)](https://pypi.org/project/veri-helix/)
[![Reproducible Viz (spec v1.0)](https://img.shields.io/badge/reproducible%20viz-spec%201.0-6f42c1)](docs/schema.md)
[![Benchmarks](https://img.shields.io/badge/bench-passing-success)](https://omniscoder.github.io/Helix/benchmarks/)

Helix is a hobbyist-first playground for **computational** bioinformatics and simulations. Think of it as a backpack full of lightweight tools, algorithms, and in-silico experiments you can remix on evenings, in classrooms, or during notebook sessions. We embrace rough edges and fast iteration so that ideas can leap from a notebook sketch to a reproducible **software** prototype quickly.

Helix complements our production platform OGN rather than competing with it. When a computational prototype proves its value, you can polish and port it into OGN. Until then, Helix is the sandbox where curiosity rules. All outputs are simulations or static analyses; no wet-lab instructions or experimental guidance are provided.

👉 Docs: https://omniscoder.github.io/Helix/

## Why Helix Exists
- Lower the barrier to tinkering: ship batteries-included examples and tiny datasets.
- Showcase approachable implementations of classic algorithms so learners can peek under the hood.
- Encourage sharing and remixing of exploratory workflows without the ceremony of production deployments.
- Offer a bridge to OGN by keeping APIs compatible and providing off-ramps when users need industrial-scale tooling.

## Highlights
- **DNA and motif experiments** (`helix.bioinformatics`): quick-and-dirty k-mer counting, SNP-tolerant motif clustering, GC skew plots, FASTA cleaning, and a CLI for summarizing GC/cluster hotspots.
- **Translation and mass lookups** (`helix.codon`, `helix.amino_acids`): resilient codon translation, bidirectional ORF scanning, frameshift heuristics, and peptide mass utilities.
- **Peptide spectrum sandbox** (`helix.cyclospectrum`): linear + cyclic theoretical spectra, scoring helpers, and a leaderboard CLI for reconstructing peptides.
- **RNA folding + ensembles** (`helix.rna`): Zuker-style MFE plus McCaskill partition/MEA/centroid helpers with dot-plots and entropy tracks.
- **CRISPR/prime scaffolding** (`helix.crispr`, `helix.prime`): PAM registry, guide discovery CLI + schema payloads, and stubs for scoring/simulation/viz.
- **PCR amplification DAGs** (`helix.pcr`): primer pair + PCR config models, DAG construction helpers, and interactive viz outputs that slot into the edit-dag tooling pipeline.
- **Protein helpers** (`helix.protein`): sequence-first summaries (weight, charge, hydropathy windows) with FASTA loading, visualization, and a friendly CLI wrapper.
- **Workflows + API** (`helix.cli`, `helix.workflows`, `helix.api`): YAML-driven automation, visualization hooks, and a pure-Python API for notebooks/scripts (full helper reference lives in [`docs/reference/api.md`](docs/reference/api.md)).
- **Seeding + seed-and-extend** (`helix.seed`): deterministic minimizers/syncmers and banded seed-extend helpers for toy mappers and density visualizations.
- **String/search helpers** (`helix.string`): FM-index construction, exact pattern search, and Myers bit-vector edit-distance for CLI/API explorations.
- **Graphs & DBG tooling** (`helix.graphs`): build/clean De Bruijn graphs, serialize to JSON/GraphML, and prep for colored/pseudoalignment experiments.
- **Motif discovery** (`helix.motif`): EM-based PWM inference (baseline MEME) with CLI/API symmetry and optional PWM plots.
- **Neural net doodles** (`ann.py`): minimal NumPy-only network for experimenting with small bio datasets.
- **Schema + provenance tooling** (`helix.schema`, `helix viz --schema`, `helix schema diff/manifest`, `helix workflows --with-schema`): every JSON artifact and PNG knows its schema kind, spec version, and SHA-256 for audit-ready reproducibility.

## Safety & Scope
- Helix operates entirely in silico: it ingests digital sequences or public reference datasets and emits JSON summaries, plots, or other software artifacts.
- The project never controls lab equipment, prescribes wet-lab steps, or suggests reagent handling; those activities belong in downstream, compliance-reviewed systems.
- Contributions, docs, and examples must stay simulation-only—if an idea drifts toward physical protocols, move it to OGN (or another governed surface) after review.
- Provenance metadata is provided so researchers can audit computational pipelines and reproducible simulations, not to guide lab execution.

## Repo Layout
```
.
├── pyproject.toml              # packaging metadata + extras
├── src/
│   └── helix/
│       ├── __init__.py         # user-facing namespace
│       ├── amino_acids.py
│       ├── bioinformatics.py
│       ├── codon.py
│       ├── cyclospectrum.py
│       ├── datasets/           # bundled FASTA/FAA toy data
│       ├── cli.py              # console entry point (`helix …`)
│       ├── api.py              # notebook-friendly helpers
│       ├── workflows.py        # YAML runner
│       ├── nussinov_algorithm.py
│       ├── protein.py
│       └── triage.py
├── examples/                   # runnable scripts + demos
├── tests/                      # pytest suites
└── README.md
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
  Builds the first-version “digital twin” graph using the new edit runtime. Nodes contain fully materialized genome views, edges capture each clean-cut event, and the JSON artifact (`helix.crispr.edit_dag.v1`) can feed notebooks and future viz surfaces.

- **Prime editing sandbox**
  ```bash
  helix prime simulate --genome genome.fna --peg-config peg.json --editor-config pe3.json --max-outcomes 16 --json prime_edits.json
  ```
  Wraps the new prime-editing models: pegRNA definitions (inline flags or JSON), prime-editor parameters, and the `simulate_prime_edit` entrypoint. Like the CRISPR command, this is purely computational—it emits hypothetical outcomes for downstream notebooks and viz.

- **Prime edit DAG**
  ```bash
  helix prime dag --genome genome.fna --peg-config peg.json --editor-config pe3.json --json prime_edit_dag.json
  ```
  Produces the prime-editing DAG artifact (`helix.prime.edit_dag.v1`) so you can inspect RTT-driven branches, log probabilities, and materialized genome snapshots per node.

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
1. Bundle a CLI command/notebook for combining GC skew, ORFs, and motif clusters into shareable reports.
2. Implement scoring for cyclo-spectrum experiments and publish a walkthrough notebook.
3. Finish the Nussinov traceback to output secondary structure strings and diagrams.
4. Add small CLIs (argparse or Typer) for swapping inputs without editing source files.
5. Draft an `examples/` gallery featuring community notebooks and weekend projects.

## Relationship to OGN
Helix is intentionally lightweight. We do not guarantee production stability, large-scale data orchestration, or SLA-backed support, and we never interface with lab automation hardware. When your prototype needs robustness, data governance, or integration with those systems:
1. Package the core logic (functions, notebooks, scripts).
2. Identify equivalent building blocks in OGN or write adapters that call into it.
3. Open an OGN ticket or PR referencing the Helix prototype so we can collaborate on the migration.

This separation keeps Helix nimble while letting OGN remain the home for hardened workflows.

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
