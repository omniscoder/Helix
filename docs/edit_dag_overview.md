Helix Edit DAG Simulation Engine — Methods Overview
===================================================

1. Introduction
---------------

Helix implements a fully in-silico digital twin of genome editing events. Instead of predicting a single CRISPR or Prime Editing outcome, Helix models editing as a probabilistic branching process over a reference sequence. The process is represented as a Directed Acyclic Graph (DAG), where:

* Nodes represent possible digital genome states.
* Edges represent editing operations (binding, cutting, repair, flap resolution).
* Log-probabilities propagate through edges to define an edit outcome distribution.

This DAG representation provides an interpretable, reproducible, and extensible framework for studying editing behavior computationally. All simulations occur on digital sequences only; no wet-lab procedures, protocols, or experimental parameters are modeled.

2. Inputs
---------

### 2.1 Reference Genome (FASTA)

Helix accepts any FASTA file:

```
--genome reference.fa
```

`DigitalGenome.from_fasta()` loads each FASTA entry as a chromosome or contig. Optionally, users may restrict simulations to a genomic subregion:

```
--region chr7:55000000-55002000
```

This enables efficient local simulations for target-specific analyses.

### 2.2 CRISPR System Configuration (cas-config)

A CRISPR nuclease is defined by a JSON configuration specifying:

* PAM pattern (IUPAC)
* Cut offset relative to PAM or guide
* Mismatch tolerance
* Seed vs. distal mismatch weighting
* Optional bulge allowance

Example (`cas9.json`):

```
{
  "name": "SpCas9",
  "type": "cas9",
  "pam_pattern": "NGG",
  "cut_offset": 3,
  "max_mismatches": 3,
  "physics": {
    "seed_length": 8,
    "seed_mismatch_penalty": 1.5,
    "distal_mismatch_penalty": 0.75,
    "pam_match_bonus": 1.0
  }
}
```

### 2.3 Guide RNAs

Guides can be passed directly:

```
--guide-sequence ACGTACGTACGTACGTACGT
```

or in batch via TSV/JSON:

```
name    sequence
G1      ACGTACGTACGTACGTACGT
G2      TGCATGCATGCATGCATGCA
```

Batch mode enables library-scale simulations using:

```
--guides-file guides.tsv
```

### 2.4 Prime Editing Configurations

Prime Editor configurations include:

* The underlying Cas nuclease
* Nick-to-edit offset
* Efficiency and indel bias parameters
* Flap-resolution logits

Example (`prime_editor.json`):

```
{
  "name": "PE2-like",
  "cas": { "name": "SpCas9-H840A", "pam_pattern": "NGG", "cut_offset": 3 },
  "nick_to_rtt_offset": 0,
  "efficiency_scale": 0.6,
  "indel_bias": 0.1,
  "flap_model": {
    "left_win_logit": 0.0,
    "right_win_logit": -0.3,
    "no_edit_logit": -1.2
  }
}
```

pegRNAs:

```
{
  "name": "peg_demo",
  "spacer": "ACGTACGTACGTACGTACGT",
  "pbs": "GCTAGCTA",
  "rtt": "TCTGACTCTCTCAGGAGTC"
}
```

Batch pegs are supported via TSV, mirroring CRISPR guides.

3. Computational Pipeline
-------------------------

Helix simulations proceed in staged computational steps, each corresponding to a rule in a pluggable rule system.

### 3.1 Target Identification

Candidate on-target and off-target locations are identified using:

* IUPAC-aware PAM recognition
* Seed-weighted mismatch penalties
* Optional bit-parallel ≤k mismatch search
* Optional bulge-aware windows

This step produces a ranked list of candidate sites.

### 3.2 DAG Runtime

The core of Helix is its edit DAG engine, which simulates editing as a state-transition system:

**Nodes**

Each node contains:

* A `DigitalGenomeView` (diff-based virtual genome)
* A `log_prob`
* metadata fields (`stage`, `time_step`, labels)

Nodes are deterministically hashed using a Merkle-like scheme:

```
node_id = SHA256(sorted(parent_ids) || event_hash || stage)[:16]
```

**Edges**

Edges describe:

* Transition rule name (e.g., `crispr.clean_cut`)
* An `EditEvent` (`chrom`, `start`, `end`, `replacement`)
* Metadata passed from physics rules

**Stages**

* CRISPR examples: `root → cut → repaired → terminal states`
* Prime editing: `root → prime_rtt → flap → repaired → terminal`

Each stage corresponds to a deterministic transformation of the genome view.

4. Editing Physics
------------------

Helix decouples physics from graph structure via a pluggable rule system.

### 4.1 CRISPR Physics

Rules include:

* `clean_cut`: Create a cut event at a recognized site
* `indel_branch`: Produce alternative repair outcomes
* `no_edit`: Failure or no-change path

Scoring includes:

* Seed and distal mismatch penalties
* PAM match bonus
* Optional bulge penalties

### 4.2 Prime Editing Physics

Rules include:

* `rtt_clean`: Apply RTT over spacer-matched site
* `flap_resolution`: Split into left/right/no-edit branches
* `indel_side_paths`: Optional micro-indel variants

Log-probabilities derive from the pegRNA parameters and PrimeEditor config.

5. Output Artifacts
-------------------

Each simulation produces a deterministic, content-addressed JSON artifact:

```
{
  "artifact": "helix.crispr.edit_dag.v1.1",
  "schema_version": "1.1",
  "root_id": "n_abc123...",
  "nodes": {
    "n_abc123": {
      "log_prob": 0.0,
      "metadata": { "stage": "root", "time_step": 0 },
      "seq_hashes": { "chrDemo": "deadbeef..." },
      "parent_ids": []
    },
    "n_def456": {
      "log_prob": -0.7,
      "metadata": { "stage": "repaired" },
      "seq_hashes": { "chrDemo": "beadfeed..." },
      "parent_ids": ["n_abc123"]
    }
  },
  "edges": [...]
}
```

Helix also supports:

* Terminal node deduplication (log-sum-exp)
* Sequence hashing for compact artifacts
* Optional sequence reconstruction for visualization

6. Visualization and Analysis
-----------------------------

Helix includes:

* PNG renderer (networkx + matplotlib)
* Interactive web Playground (Cytoscape.js)

Filters for:

* minimum probability
* maximum time-step

Sequence diff viewer and outcome distribution charts.

Metrics:

* entropy of outcome distribution
* effective number of modes
* number of terminal branches
* edit distances

These tools enable rapid comparison of designs and intuitive inspection of editing behavior.

7. Use Cases (In-Silico Only)
-----------------------------

Helix supports computational workflows such as:

* CRISPR/Prime editing behavior exploration
* Comparative guide efficiency evaluation
* pegRNA design searches
* Computational robustness testing
* ML dataset generation from synthetic edits
* Visualization and interpretability research

All performed without any wet-lab protocols or biological manipulation.

8. Extensibility
----------------

Helix is designed for future expansion:

* New nucleases
* New Prime Editing variants
* Digital PCR + sequencing simulators
* Variant caller twins
* Parameter-fitting for calibrating physics models
* GPU acceleration

Rule-based physics and content-addressed DAGs make Helix a flexible platform for in-silico experimentation.

Final Notes
-----------

This document formalizes Helix’s simulation engine as if it were a published-methods section. It positions the CRISPR/Prime DAG engine as:

* Mechanistic
* Deterministic
* Extensible
* In-silico
* Reproducible
* Research-grade
