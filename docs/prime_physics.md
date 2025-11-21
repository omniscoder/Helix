# Prime Physics Score (v0)

Prime physics scoring adds interpretable heuristics on top of prime.edit_sim
outputs. Every simulation can now annotate the metadata with a
`physics_score` object (use `--physics-score` on the CLI) containing:

- **PBS ΔG (`pbs_dG`)** – coarse binding stability of the primer binding site.
  Values near `-10 kcal/mol` indicate strong binding; positive values flag weak
  anchoring.
- **RT path ΔG (`rt_cum_dG`)** – cumulative thermodynamic cost of copying the
  RTT. The final entry shows the total energy change after synthesizing the
  template.
- **Flap ΔΔG (`flap_ddG`) & microhomology** – flap resolution prefers a short,
  favorable ΔΔG and several nucleotides of microhomology. Low microhomology or
  large positive ΔΔG will push probability into reanneal/indel branches.
- **Mismatch/repair hints** – `mmr_flag` toggles when the RTT falls outside the
  editor’s mismatch tolerance, and `nick_distance` reports how far the nick is
  from the RTT terminus.
- **Probabilities and predicted efficiency** – `P_RT` models the probability of
  completing RTT synthesis, `P_flap` is the flap resolution prior, and the final
  `E_pred` headline combines both terms with the editor’s efficiency scale.

## CLI example

```
helix prime simulate \
  --genome /tmp/genome.fna \
  --peg-config peg.json \
  --editor-config editor.json \
  --physics-score \
  --json prime.json
```

Trimmed JSON output from a sample run:

```json
"physics_score": {
  "pbs_dG": 2.50,
  "rt_cum_dG": [-1.00, -0.25, 0.50, 1.25, 1.75, 2.50],
  "flap_ddG": 2.30,
  "microhomology": 1,
  "mmr_flag": true,
  "nick_distance": 5,
  "P_RT": 0.182,
  "P_flap": 0.373,
  "E_pred": 0.054
}
```

Interpretation:

- PBS ΔG is slightly positive → binding is marginal.
- Microhomology is only 1 nt with a positive flap ΔΔG → flap resolution is
  penalty heavy.
- Resulting `E_pred ≈ 0.05` signals this peg is unlikely to perform well.

Tweaking the peg to lengthen the PBS (to 12–14 nt) and increase microhomology
pushes `P_RT`/`P_flap` higher and lifts `E_pred` into the ~0.25–0.35 range. Use
this score to triage candidates before investing in full DAG simulations.
