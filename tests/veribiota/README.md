# VeriBiota Micro-DAG Corpus

These helpers generate the tiny CRISPR/Prime/PCR edit DAGs that feed the VeriBiota proof harness. Run them from the repository root:

```bash
python tests/veribiota/gen_crispr_micro.py
python tests/veribiota/gen_prime_micro.py
python tests/veribiota/gen_pcr_micro.py
```

Each script emits:

```
veribiota_work/<mechanism>_micro.dag.json
veribiota_work/<mechanism>_micro.lean-check.json
```

The DAG JSONs capture the raw Helix edit graphs, and the lean-check summaries encode the structural / probability invariants consumed by the VeriBiota adapter and CI action.
