# Helix Studio Strategy Playbook

## 0. Positioning
Helix Studio is the Visual Studio of bioinformatics: a multi-process, snapshot-driven IDE that designs, simulates, compares, and reports on CRISPR/Prime, PCR, signaling, and cell dynamics—with a headless CLI and GPU backend (OGN) built for real pipelines.

## 1. Strategic Arc (Wedge → Platform → Ecosystem)
- **Wedge (Months 0–6):** Dominate CRISPR/Prime design loops with an end-to-end workflow (design helpers → run → compare → report) and stand up PCR as proof of the multi-process architecture.
- **Platform (Months 6–12):** Integrate OGN for remote GPU runs + batch sweeps, ship a headless CLI with Nextflow/Cromwell modules, and publish a versioned snapshot specification.
- **Ecosystem (Months 9–18):** Deliver a plugin API + scripting surface, a marketplace for metrics/panels/run kinds, and importers for wet-lab outputs (CRISPResso, AMP-seq, BAM/VCF) so real labs can unify wet + dry results.

## 2. North-Star & KPIs
- **North-Star:** Successful design loops per week (create design → run → compare → report).
- **Supporting KPIs:**
  - Time-to-first-verdict ≤ 10 minutes from new project.
  - Snapshot share rate ≥ 30% of sessions export a report.
  - Compare usage ≥ 60% of runs get compared.
  - Headless adoption ≥ 50% of runs via CLI/OGN by Month 9.
  - Engineering reliability: job success ≥ 99%, p95 run metadata ingest < 5 seconds.
- Instrument each metric from day one.

## 3. 12-Month Milestone Plan (Branches)
### Q1 (Months 0–3) – Nail the Wedge
- **v1.0 (stabilize):** CRISPR/Prime IDE spine, snapshots, compare, reports (`release/v1.0`).
- **v1.1:** Design helpers + verdicts + report v1.1 + snapshot migrations (`feature/studio-v1-1-design-helpers` → `release/v1.1`).
- **Go-to-market:** LOIs with 5 design-partner labs, weekly design reviews.

### Q2 (Months 3–6) – Prove Multi-Process
- **v2.0:** PCR engine + metrics + compare + reports (`feature/studio-pcr` → `release/v2.0`).
- **Importers α:** CRISPResso/AMP-seq summary importer to align sims ↔ wet-lab (`feature/studio-importers-alpha`).
- **Docs:** Public snapshot spec v1 + CLI alpha.

### Q3 (Months 6–9) – Become a Platform
- **v5.0 (pulled forward):** OGN integration (remote jobs, perf overlay, batch sweeps, remote snapshot ingest) via `feature/studio-ogn-gpu` → `release/v5.0`.
- **Headless:** CLI 1.0 + Nextflow/Cromwell modules (`feature/studio-cli-headless`).
- **Telemetry:** Privacy-safe metrics + opt-in crash reports + reliability SLOs.

### Q4 (Months 9–12) – Open the Ecosystem
- **v6.0:** Plugin API + scripting console + notebooks (`feature/studio-plugins` → `release/v6.0`).
- **Marketplace β:** Signed plugins, verified authors, revenue share.
- **Case studies:** Prime campaign, PCR assay, OGN batch with figure-quality reports.

### Parallel Track (Optional)
- **v3.0:** Signaling/expression to round out cell-logic story (`feature/studio-signaling`).

## 4. Technical Moats
- Snapshot as Contract: schema-versioned, content-addressed, deterministic rehydration—Terraform for biology.
- Compare/Verdict UX: one-click “Improved/Worse/Tradeoff” with transparent diffs across run kinds.
- OGN Flight Deck: GPU job control with perf HUD + batch sweep grid.
- Plugin safety + determinism: isolated workers, time-boxed calls, reproducible notebooks referencing snapshot IDs.
- Real-data calibration loop: import wet-lab results next to sims, compute calibration error, bake into reports.

## 5. Headless & Pipeline Integration
### CLI (v1 syntax)
```bash
helix run \
  --snapshot session.hxs \
  --kind Prime \
  --params params.json \
  --engine local \
  --out out/

helix compare --run A123 --run B456 --out compare.json

helix report --run B456 --format md --out report.md

helix run \
  --snapshot session.hxs \
  --kind CRISPR \
  --engine ogn --queue gpu-a100 \
  --out s3://my-bucket/helix/A123/
```

### Snapshot Bundle (.hxs)
Zipped JSON manifest + referenced assets (refs, primers, pathway graphs) with a content-address map for reproducibility.

### Nextflow Module
```groovy
process HELIX_PRIME {
  container 'helixstudio/cli:1.0'
  input:
    path snapshot
    path params_json
  output:
    path "out/**"
    path "report.md"
  script:
  """
  helix run \
    --snapshot ${snapshot} \
    --kind Prime \
    --params ${params_json} \
    --engine ogn \
    --out out/
  helix report --run $(cat out/run_id.txt) --format md --out report.md
  """
}
```

### Cromwell/WDL Snippet
```wdl
task HelixCRISPR {
  input {
    File snapshot
    File params_json
  }
  command <<<
    helix run --snapshot ~{snapshot} --kind CRISPR --engine ogn --out out/
    helix report --run `cat out/run_id.txt` --format json --out report.json
  >>>
  output {
    File report = "report.json"
    Directory artifacts = "out"
  }
}
```

## 6. Plugin API Shape
### Manifest (TypeScript)
```ts
export interface HelixPluginManifest {
  name: string; version: string;
  helixApi: { min: "6.0.0"; max: "6.x" };
  capabilities: ("runKind"|"panel"|"metric"|"analysis")[];
  entry: string;
  permissions?: ("fs.read"|"net.fetch"|"gpu.compute")[];
  signature?: string;
}
```

### Activation Example
```ts
export default function activate(ctx: HelixContext) {
  ctx.runs.register("AmpliconQC", {
    validateParams: schemaCheck,
    simulate: async (p, state) => shellOut("amplicon-qc", p, state),
    summarize: metricsToVerdict,
    report: renderReport
  });
  ctx.metrics.register("OnTargetDelta", computeOnTargetDelta);
  ctx.panels.register("AmpliconCloud", AmpliconCloudPanel);
}
```

### Security
Plugins run in isolated workers with structured-clone IPC, resource caps, and signed marketplace packages.

## 7. Market Motion
- **Ideal customers:** Methods teams (CRISPR/Prime/PCR) at biotech/pharma, CROs, advanced academic cores.
- **Design partner program (Q1–Q2):** 5–8 labs, LOIs, shared roadmap, monthly office hours, co-authored case studies.
- **Pricing & packaging (open-core):**
  - Community (free): local engine, CRISPR/Prime, reports, limited importers.
  - Pro ($/seat/mo): PCR, Compare+, advanced reports, local batch sweeps, plugin install.
  - Teams ($/user/mo + org features): shared storage, access control, private plugin registry.
  - Enterprise (site license + support): OGN remote, SSO/SAML, audit logs, cold storage, priority support.
  - Compute: pass-through for OGN GPU hours at transparent rates.
  - Marketplace: 80/20 revenue share for paid plugins.
- **Growth loops:** Shareable reports, “Launch in Helix” badges, template gallery, one-click publish to lab wiki.

## 8. Proof & Credibility
Validation packs, public datasets, golden snapshots, reproducible scripts, calibration curves, preprint on snapshot model, live demos with batch sweeps + OGN perf HUD.

## 9. Enterprise Checklist
Security (content-addressed artifacts, encryption, SSO/SAML, RBAC), auditability (append-only provenance, signed snapshots), compliance (SOC2 path, HIPAA readiness), data governance (project-scoped storage, customer-managed keys).

## 10. Team & Resourcing
Lean, senior team: PM, Tech Lead, 2 engine engineers (CRISPR/Prime/PCR), 1 graphics/viz, 1 infra (OGN/CLI/CI), 1 DevRel/Docs, 1 QA/SDET. Scale in Q3 with plugin platform engineer + importers specialist.

## 11. Risk Register & Mitigations
- **Perception of inaccuracy:** Publish calibration packs, display uncertainty, tunable thresholds, pin engine versions per snapshot.
- **GUI misfit for pipelines:** Snapshot as contract + CLI + Nextflow/Cromwell.
- **Plugin security:** Worker isolation, permission manifest, signed packages, review program.
- **OGN reliability under load:** Backpressure, retries, detailed job states, safe resume for partial runs.
- **Domain fragmentation:** Keep State/Dynamics/Viz pattern; new domains plug into shared history/compare/reports.

## 12. Concrete Artifacts (This Week)
- Snapshot Spec v1 (JSON schema, canonicalization, hashing, worked CRISPR + Prime examples).
- CLI α (helix run/compare/report wired to local engines).
- Design Partner Pack (30-min demo script, 3 sample sessions, logging feedback doc).
- Template Gallery v0 (5 Prime presets, 3 PCR presets).
- OGN Contracts Draft (EngineJobSpec/Status/SnapshotPayload + mock server).
- Website README framing (positioning line, design loop GIF, links to snapshot spec + CLI).

## 13. Branch & Release Hygiene
Long-lived `helix-studio`. Version branches `release/vX.Y`, tags `vX.Y.Z`. Feature branches per epic. CI gates: typecheck, engine tests, golden snapshot diff, report rendering diff, CLI smoke. Every release ships 1 demo session, 1 tutorial, 1 case study outline.

## 14. Opinionated Yes/No
Yes: snapshot as immutable content-addressed bundle; deterministic engines with fixed seeds + pinned versions; compare/verdict as one button with transparent math.
No: adding domains without Metrics/Compare/Report; feature flags that bypass snapshot schema migration; plugins with unbounded FS/network access.

## 15. User Value (Why They Switch)
Design loops end in verdicts; consolidated view across CRISPR/Prime/PCR/signaling/cell outcomes; single artifact (snapshot) that runs locally or on OGN/Nextflow; publication-grade reports; plugin ecosystem where scripts become first-class citizens.

## 16. Immediate Next Steps (Execution Checklist)
- Create repos/docs for Snapshot Spec v1 and CLI α, publish minimal examples.
- Open `feature/studio-v1-1-design-helpers` and land guide/peg suggestions + verdict aggregator behind a flag.
- Start design-partner cadence: weekly sessions, track TTFV/compare usage/report exports.
- Draft OGN client stubs and local fake for job lifecycle testing.
- Write Nextflow and Cromwell examples and check them into the repo.
- Ship v1.1 with magical design helpers + verdicts, then harden headless loops at scale.
