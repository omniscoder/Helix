# VeriBiota Verification Pipeline

The **Verified by VeriBiota** badge certifies that Helix’s edit DAGs pass a formal verification gauntlet built on Lean + Lake. We keep the standard crisp so the badge remains meaningful.

## Badge Contract

A release (or repo) may display the badge only if every CRISPR/Prime/PCR DAG used in that context satisfies:

1. **Structural correctness** – there is exactly one root, the DAG is acyclic, node depths increase monotonically, and every edge corresponds to a legal transition in the simulator.
2. **Semantic correctness** – each edge’s `EditEvent` matches the formal rewrite rule for the mechanism (CRISPR, Prime, PCR) encoded in Lean. No ad-hoc transitions.
3. **Probability sanity** – for all nodes, outgoing edge probabilities sum to ≈1, and the terminal (leaf) probabilities sum to ≈1 with no negative/`NaN` weights.
4. **Reproducibility** – replays using the recorded RNG seed and provenance reproduce the same DAG; streamed frames equal the static DAG snapshot; sequence hashes line up.

If the goal is “badge everywhere”, every DAG emitted by Helix must pass these constraints.

## CLI Pipeline

1. **Generate DAGs** with the helper scripts in `tests/veribiota/` (they wrap the usual Helix CLIs and drop outputs under `veribiota_work/`):

   ```bash
   python tests/veribiota/gen_crispr_micro.py
   python tests/veribiota/gen_prime_micro.py
   python tests/veribiota/gen_pcr_micro.py
   ```
2. **Emit lean-check metadata** – for each DAG JSON, run:
   ```bash
   helix veribiota lean-check --input dag.json --out dag.lean-check.json
   ```
   This records node/edge counts, log-probabilities, and sequence hashes.
3. **Preflight** – validate the lean-check files (optionally against the source DAGs) before Lean boots:
   ```bash
   helix veribiota preflight --checks out/*.lean-check.json --payloads out/*.json
   ```
   CI fails here if any invariant or hash check fails.
4. **Consolidate DAGs** – once preflight passes, generate one Lean module that defines every DAG plus an aggregate list/theorem:
   ```bash
   helix veribiota export-dags \
     --inputs out/dag*.json \
     --module-name Helix.CrisprExamples \
     --list-name exampleDags \
     --out veribiota/generated/Helix/CrisprExamples.lean
   ```
5. **Lean / Lake proofs** – `lake build` and `lake exe veribiota-check veribiota/generated/Helix/CrisprExamples.lean`. Our simple checker currently type-checks the file; downstream projects can extend it with the full VeriBiota proof stack.

When collaborating with the external [`VeriBiota/VeriBiota`](https://github.com/VeriBiota/VeriBiota) repository, use `helix veribiota export-suite` instead. It writes the Lean file straight into the VeriBiota checkout, ready for its Lake build:
```bash
helix veribiota export-suite \
  --inputs out/dag*.json \
  --veribiota-root ../VeriBiota \
  --module-path Biosim/VeriBiota/Helix/MicroSuite.lean \
  --module-name Biosim.VeriBiota.Helix.MicroSuite \
  --list-name allDags
```
`export-suite` shares the same template as `export-dags` but respects VeriBiota’s namespace conventions (`Biosim.VeriBiota.Helix.*`) and writes directly into the repo so Lake can import `Biosim/VeriBiota/Helix/MicroSuite.lean`.

## GitHub Action

The workflow in `.github/workflows/veribiota.yml` wires these steps together:

1. Install Helix, generate the micro-universe DAG, emit lean-check metadata, and run `preflight`.
2. Checkout the VeriBiota repo at a pinned revision, run `helix veribiota export-suite` to populate `Biosim/VeriBiota/Helix/*.lean`, install Lean + Lake, and run `lake build` + `lake exe veribiota-check`.
3. Surface the action status via a badge (`VeriBiota CI`).

Projects inheriting from Helix can reuse this `veribiota.yml` as a template—just swap the DAG inputs for your suite. Keep the badge honest!
