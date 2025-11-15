/* global cytoscape */

import { initBranchChart, updateBranchChart } from "./realtime_branch_chart.js";

const worker = new Worker("js/realtime_worker.js", { type: "module" });
const frames = [];
const comparisonState = {
  worker: new Worker("js/realtime_worker.js", { type: "module" }),
  framesA: [],
  framesB: [],
  pending: null,
};
let cytoscapeInstance;
let comparisonCy;
let branchChart;
let playTimer = null;
let currentIndex = 0;
let rootSeq = null;
let rootChrom = null;

function initCytoscape() {
  cytoscapeInstance = cytoscape({
    container: document.getElementById("rt-cy"),
    style: [
      {
        selector: "node",
        style: {
          "shape": "round-rectangle",
          "background-color": (ele) => {
            const meta = ele.data("raw_meta") || {};
            let editClass = ele.data("edit_class") || meta.edit_class;
            if (!editClass) {
              const stage = ele.data("stage");
              const branch = meta.branch;
              if (stage === "repaired") {
                if (branch === "intended") editClass = "substitution";
                else if (branch === "indel") editClass = "deletion";
              } else if (stage === "no_edit" || stage === "root") {
                editClass = "no_cut";
              }
            }
            if (editClass === "deletion") return "#f97316"; // orange
            if (editClass === "insertion") return "#ec4899"; // magenta
            if (editClass === "substitution") return "#22d3ee"; // turquoise
            if (editClass === "no_cut") return "#64748b"; // gray
            const stage = ele.data("stage");
            if (stage === "prime_rtt") return "#10b981"; // green
            if (stage === "repaired") return "#3b82f6"; // blue
            if (stage === "error") return "#ef4444"; // red
            if (stage && String(stage).startsWith("flap")) return "#f59e0b"; // amber
            return "#4f46e5"; // indigo default
          },
          "width": "mapData(prob, 0, 1, 24, 64)",
          "height": "mapData(prob, 0, 1, 16, 48)",
          "label": "data(label)",
          "color": "#e5e7eb",
          "font-size": 10,
          "text-wrap": "wrap",
          "border-color": (ele) => {
            const d = ele.data("desired");
            if (d === true) return "#10b981"; // green border for desired
            if (d === false) return "#ef4444"; // red border for undesired
            return "#0b1220"; // default
          },
          "border-width": (ele) => {
            const ph = (ele.data("phenotype") || "").toString().toLowerCase();
            if (ph === "viable") return 3;
            if (ph === "lethal") return 1;
            return 2; // unknown
          },
          "border-style": (ele) => {
            const ph = (ele.data("phenotype") || "").toString().toLowerCase();
            return ph ? "solid" : "dashed";
          },
          // subtle halo based on desired flag
          "shadow-blur": (ele) => (typeof ele.data("desired") === "boolean" ? 12 : 0),
          "shadow-color": (ele) => (ele.data("desired") === true ? "#10b981" : ele.data("desired") === false ? "#ef4444" : "#000000"),
          "shadow-opacity": (ele) => (typeof ele.data("desired") === "boolean" ? 0.25 : 0),
          "shadow-offset-x": 0,
          "shadow-offset-y": 0,
        },
      },
      {
        selector: "edge",
        style: {
          "curve-style": "bezier",
          "line-color": (ele) => {
            const rule = ele.data("rule") || "";
            if (rule.includes("prime.rtt_clean")) return "#10b981"; // green
            if (rule.includes("flap")) return "#f59e0b"; // amber
            if (rule.includes("indel") || rule.includes("error")) return "#ef4444"; // red
            if (rule.includes("clean_cut")) return "#8b5cf6"; // violet
            return "#94a3b8";
          },
          "target-arrow-color": "data(line-color)",
          "target-arrow-shape": "triangle",
          "width": "mapData(weight, 0, 1, 1, 4)",
          "opacity": "mapData(weight, 0, 1, 0.25, 1)",
          "label": "data(rule_human)",
          "font-size": 11,
          "color": "#e5e7eb",
          "text-rotation": "autorotate",
          "text-outline-width": 2,
          "text-outline-color": "#0b1220",
          "text-background-color": (ele) => {
            const tag = ele.data("effect_tag");
            if (tag === "frameshift") return "#3b1d1d"; // reddish
            if (tag === "inframe") return "#0d2f2f"; // teal-ish
            if (tag === "snv") return "#1f2937"; // neutral dark
            return "#0b1220";
          },
          "text-background-opacity": 0.55,
          "text-background-padding": 2,
          "text-background-shape": "roundrectangle",
          "text-wrap": "wrap",
          "text-margin-y": -6,
        },
      },
      { selector: ".faded", style: { "opacity": 0.12 } },
      { selector: ".hidden", style: { "display": "none" } },
    ],
    wheelSensitivity: 0.2,
  });
  comparisonCy = cytoscape({
    container: document.getElementById("rt-compare-cy"),
    style: [
      {
        selector: "node",
        style: {
          "background-color": ele => {
            const group = ele.data("group");
            if (group === "A") return "#10b981";
            if (group === "B") return "#f97316";
            return "#fbbf24";
          },
          label: "data(label)",
          color: "#e5e7eb",
          "font-size": 8,
          "text-wrap": "wrap",
        },
      },
      {
        selector: "edge",
        style: {
          "line-color": ele => (ele.data("group") === "A" ? "#059669" : "#ea580c"),
          "target-arrow-color": ele => (ele.data("group") === "A" ? "#059669" : "#ea580c"),
          "target-arrow-shape": "triangle",
          width: 1.2,
        },
      },
    ],
  });
}

function applyFrame(frame) {
  const cy = cytoscapeInstance;
  function fmtNum(n) {
    try { return Number(n).toLocaleString(); } catch { return String(n); }
  }
  function effectLabel(ev) {
    if (!ev) return null;
    const len = Math.max(0, (ev.end || 0) - (ev.start || 0));
    const ins = (ev.replacement || "").length;
    const delta = ins - len;
    if (delta !== 0) return Math.abs(delta) % 3 === 0 ? "in‑frame" : "frameshift";
    if (len > 0 && ins > 0) return len === 1 ? "SNV" : "substitution";
    return null;
  }
  function changeSummary(ev) {
    if (!ev) return "";
    const len = Math.max(0, (ev.end || 0) - (ev.start || 0));
    const ins = (ev.replacement || "").length;
    const pos = `${ev.chrom || "chr"}:${fmtNum((ev.start || 0) + 1)}`;
    const eff = effectLabel(ev);
    if (ins > 0 && len === 0) return `+${ins} bp insertion${eff ? ` (${eff})` : ""} • ${pos}`;
    if (ins === 0 && len > 0) return `−${len} bp deletion${eff ? ` (${eff})` : ""} • ${pos}`;
    if (ins > 0 && len > 0) {
      if (len === 1 && ins === 1) {
        const from = rootSeq ? rootSeq[(ev.start || 0)] : null;
        const to = (ev.replacement || "");
        return `${from && to ? `${from}→${to}` : 'SNV'} • ${pos}`;
      }
      return `substitution (${len}→${ins} bp) • ${pos}`;
    }
    return `edit • ${pos}`;
  }
  function mechanismFromRule(rule) {
    if (!rule) return null;
    if (rule.includes("clean_cut")) return "CRISPR cleavage";
    if (rule.includes("indel")) return "NHEJ";
    if (rule.includes("prime")) return "Prime Editing";
    if (rule.includes("flap")) return "Flap resolution";
    if (rule.includes("no_edit")) return "No Change";
    return null;
  }
  function humanEdgeLabel(rule) {
    if (!rule) return "";
    if (rule.includes("clean_cut")) return "🔪 Cas9 cleavage";
    if (rule.includes("indel")) return "✂️ Repair via NHEJ";
    if (rule.includes("prime.rtt_clean")) return "🧬 Prime RT extension";
    if (rule.includes("flap")) return "🔧 Flap resolution";
    if (rule.includes("no_edit")) return "No incorporation";
    return rule;
  }
  function effectTag(ev) {
    if (!ev) return null;
    const len = Math.max(0, (ev.end || 0) - (ev.start || 0));
    const ins = (ev.replacement || "").length;
    const delta = ins - len;
    if (delta !== 0) return Math.abs(delta) % 3 === 0 ? "inframe" : "frameshift";
    if (len > 0 && ins > 0 && len === 1 && ins === 1) return "snv";
    return null;
  }
  function displayLines(stage, ev, meta, ruleName) {
    const desired = meta?.desired;
    if (stage === "root") return ["Root Population", "Wild‑Type Genome"];
    const mech = mechanismFromRule(ruleName) || (stage === "error" ? "NHEJ" : stage === "prime_rtt" ? "Prime Editing" : stage === "no_edit" ? "No Change" : stage === "cut" ? "CRISPR cleavage" : null);
    if (stage === "cut") return ["Double‑Strand Break", `Cas9 cut • ${changeSummary(ev).replace(/.*@ /,'')}`];
    if (stage === "repaired") return [desired === false ? "NHEJ‑Indel Outcome" : "✨ HDR Result", `${mech || "Repair"} • ${changeSummary(ev).replace('@ ', '• ')}`];
    if (stage === "error") return ["✂️ NHEJ Outcome", `${mech || "NHEJ"} • ${changeSummary(ev).replace('@ ', '• ')}`];
    if (stage === "no_edit") return ["Wild‑Type Allele", "No sequence change"];
    if (stage && String(stage).startsWith("flap")) return ["Flap Resolution", `${mech || "Flap"} • ${changeSummary(ev).replace('@ ', '• ')}`];
    if (stage === "prime_rtt") return ["🧬 Prime Edit", `${mech || "Prime Editing"} • ${changeSummary(ev).replace('@ ', '• ')}`];
    return [String(stage || "State"), `${mech ? mech + " • " : ""}${changeSummary(ev)}`];
  }
  Object.entries(frame.new_nodes || {}).forEach(([nodeId, node]) => {
    if (cy.$id(nodeId).length) return;
    const lines = displayLines(node.metadata?.stage, null, node.metadata || {});
    const label = `${lines[0]}\n${lines[1] || ""}`;
    const sequences = node.sequences || {};
    const prob0 = Math.exp(node.log_prob || 0);
    const meta = node.metadata || {};
    cy.add({
      group: "nodes",
      data: {
        id: nodeId,
        label,
        stage: node.metadata?.stage,
        log_prob: node.log_prob,
        prob: prob0,
        phenotype: meta.phenotype || "",
        time: meta.time_step ?? null,
        sequences,
        on_target: !!meta.on_target,
        desired: typeof meta.desired === "boolean" ? meta.desired : undefined,
        raw_meta: meta,
      },
    });
  });
  const maxProb = Math.max(
    1e-12,
    ...Object.values(frame.new_nodes || {}).map((n) => Math.exp(n.log_prob || 0)),
    ...cy.nodes().map((n) => n.data("prob") || 0)
  );
  Object.entries(frame.new_nodes || {}).forEach(([nodeId, node]) => {
    const prob = Math.exp(node.log_prob || 0);
    const ele = cy.$id(nodeId);
    if (ele && ele.length) {
      ele.data("prob", prob);
      ele.data("phenotype", node.metadata?.phenotype || "");
      ele.data("raw_meta", node.metadata || {});
      return;
    }
    const lines = displayLines(node.metadata?.stage, null, node.metadata || {});
    const label = `${lines[0]}\n${lines[1] || ""}`;
    const sequences = node.sequences || {};
    const chroms = Object.keys(sequences);
    if (frame.step === 0 && node.metadata?.stage === "root" && chroms.length) {
      rootChrom = chroms[0];
      rootSeq = sequences[rootChrom];
    }
    cy.add({
      group: "nodes",
      data: {
        id: nodeId,
        label,
        stage: node.metadata?.stage,
        log_prob: node.log_prob,
        prob,
        phenotype: node.metadata?.phenotype || "",
        time: node.metadata?.time_step ?? null,
        sequences,
        on_target: !!node.metadata?.on_target,
        edit_class: node.metadata?.edit_class,
        desired: typeof node.metadata?.desired === "boolean" ? node.metadata.desired : undefined,
        raw_meta: node.metadata || {},
      },
    });
  });
  (frame.new_edges || []).forEach((edge) => {
    const edgeId = `${edge.source}_${edge.target}_${edge.rule}`;
    if (cy.$id(edgeId).length) return;
    const tgt = cy.$id(edge.target);
    const weight = tgt.length ? tgt.data("prob") || 0 : 0;
    const edgeMeta = edge.metadata || {};
    const mechanism = edgeMeta.mechanism || mechanismFromRule(edge.rule);
    cy.add({
      group: "edges",
      data: {
        id: edgeId,
        source: edge.source,
        target: edge.target,
        rule: edge.rule,
        mechanism,
        rule_human: edgeMeta.human_label || humanEdgeLabel(edge.rule),
        weight,
        event: edge.event || null,
        effect_tag: effectTag(edge.event),
      },
    });
    // attach last event to target node for detail rendering
    if (tgt && tgt.length && edge.event) {
      tgt.data("event", edge.event);
      // update target label now that we know the event
      const tStage = tgt.data("stage");
      const tMeta = cytoscapeInstance.$id(edge.target).data();
      const lines = displayLines(tStage, edge.event, tMeta || {}, edge.rule);
      tgt.data("label", `${lines[0]}\n${lines[1] || ""}`);
    }
  });
  cy.layout({ name: "breadthfirst", directed: true, spacingFactor: 1.05, padding: 20 }).run();
  updateMetrics(frame);
  // apply current filter settings after each frame
  if (typeof applyOutcomeFilters === "function") {
    applyOutcomeFilters();
  }
}

function summarizeBranchProbabilities(nodes) {
  return nodes.reduce(
    (acc, node) => {
      const stage = node.data("stage") || node.data("metadata")?.stage || "unknown";
      const prob = Math.exp(node.data("log_prob") || 0);
      if (stage === "repaired") acc.intended += prob;
      else if (stage === "error") acc.indel += prob;
      else if (stage === "no_edit" || stage === "root") acc.noEdit += prob;
      return acc;
    },
    { intended: 0, indel: 0, noEdit: 0 },
  );
}

function updateMetrics(frame) {
  const metricsBox = document.getElementById("rt-metrics");
  const nodes = cytoscapeInstance.nodes();
  const entropy = nodes
    .toArray()
    .reduce((acc, node) => {
      const prob = Math.exp(node.data("log_prob") || 0);
      return acc - (prob ? prob * Math.log(prob) : 0);
    }, 0);
  const branchTotals = summarizeBranchProbabilities(nodes);
  const total = nodes.length;
  metricsBox.innerHTML = `
    <div>Frame: ${frame.step}</div>
    <div>Nodes: ${total}</div>
    <div>Entropy: ${entropy.toFixed(3)}</div>
    <div>Intended prob: ${branchTotals.intended.toFixed(3)}</div>
    <div>Indel prob: ${branchTotals.indel.toFixed(3)}</div>
    <div>No-edit prob: ${branchTotals.noEdit.toFixed(3)}</div>
  `;
  if (branchChart) updateBranchChart(branchChart, branchTotals.intended, branchTotals.indel, branchTotals.noEdit);
}

function resetPlayback() {
  frames.length = 0;
  currentIndex = 0;
  document.getElementById("rt-slider").value = 0;
  cytoscapeInstance.elements().remove();
  document.getElementById("rt-metrics").textContent = "Load a simulation…";
  if (branchChart) updateBranchChart(branchChart, 0, 0, 0);
}

function playFrames() {
  if (!frames.length) return;
  if (playTimer) {
    clearInterval(playTimer);
    playTimer = null;
    return;
  }
  playTimer = setInterval(() => {
    if (currentIndex >= frames.length) {
      clearInterval(playTimer);
      playTimer = null;
      return;
    }
    applyFrame(frames[currentIndex]);
    document.getElementById("rt-slider").value = currentIndex;
    currentIndex += 1;
  }, Number(document.getElementById("rt-speed").value) || 500);
}

function loadFrameAt(idx) {
  cytoscapeInstance.elements().remove();
  for (let i = 0; i <= idx && i < frames.length; i++) {
    applyFrame(frames[i]);
  }
  currentIndex = idx + 1;
}

function exportExperiment() {
  const config = {
    kind: "helix.crispr.experiment.v1",
    name: "realtime_demo",
    description: "Exported from Helix realtime simulator.",
    genome: {
      fasta: "inline://browser",
      region: null,
    },
    cas: { config: "inline://SpCas9" },
    guide: {
      sequence: document.getElementById("rt-guide").value.trim(),
      name: document.getElementById("rt-guide").value ? "guide" : null,
      pam: document.getElementById("rt-pam").value.trim(),
    },
    simulation: {
      max_depth: 2,
      min_prob: 1e-4,
      max_sites: Number(document.getElementById("rt-max-sites").value),
      seed: 0,
    },
    metadata: {
      genome_inline: document.getElementById("rt-genome").value,
      mismatch_tolerance: Number(document.getElementById("rt-mismatch").value),
      indel_window: Number(document.getElementById("rt-window").value),
    },
  };
  const blob = new Blob([JSON.stringify(config, null, 2)], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "helix_realtime_experiment.json";
  link.click();
  URL.revokeObjectURL(link.href);
}

  function runSimulation() {
    resetPlayback();
  const payload = {
    genome: document.getElementById("rt-genome").value,
    guide: document.getElementById("rt-guide").value,
    pam: document.getElementById("rt-pam").value,
    mismatchTolerance: document.getElementById("rt-mismatch").value,
    maxSites: document.getElementById("rt-max-sites").value,
    window: document.getElementById("rt-window").value,
    cdsStart: document.getElementById("rt-cds-start").value,
    cdsEnd: document.getElementById("rt-cds-end").value,
    cdsStrand: document.getElementById("rt-cds-strand").value,
    exons: document.getElementById("rt-exons").value,
  };
    worker.postMessage({ type: "simulate", payload });
  }

const compareWorker = new Worker("js/realtime_worker.js", { type: "module" });
const compareState = { framesA: [], framesB: [], pending: null };

function runComparison() {
  comparisonState.framesA = [];
  comparisonState.framesB = [];
  comparisonState.pending = "A";
  const basePayload = {
    genome: document.getElementById("rt-genome").value,
    pam: document.getElementById("rt-compare-pam").value,
    mismatchTolerance: document.getElementById("rt-compare-mismatch").value,
    maxSites: document.getElementById("rt-max-sites").value,
    window: document.getElementById("rt-window").value,
  };
  comparisonState.worker.postMessage({
    type: "simulate",
    payload: { ...basePayload, guide: document.getElementById("rt-compare-guide-a").value },
  });
  comparisonState.worker.onmessage = (event) => {
    const { type, frame, error } = event.data || {};
    if (type === "error") {
      document.getElementById("rt-compare-output").textContent = error;
      return;
    }
    if (type === "frame" && comparisonState.pending === "A") {
      comparisonState.framesA.push(frame);
    } else if (type === "frame" && comparisonState.pending === "B") {
      comparisonState.framesB.push(frame);
    } else if (type === "done" && comparisonState.pending === "A") {
      comparisonState.pending = "B";
      comparisonState.worker.postMessage({
        type: "simulate",
        payload: { ...basePayload, guide: document.getElementById("rt-compare-guide-b").value },
      });
    } else if (type === "done" && comparisonState.pending === "B") {
      comparisonState.pending = null;
      renderComparison();
    }
  };
}

function renderComparison() {
  comparisonCy.elements().remove();
  const nodeScoresA = {};
  const nodeScoresB = {};
  comparisonState.framesA.forEach((frame) => {
    Object.entries(frame.new_nodes || {}).forEach(([id, node]) => {
      nodeScoresA[id] = Math.exp(node.log_prob || 0);
    });
  });
  comparisonState.framesB.forEach((frame) => {
    Object.entries(frame.new_nodes || {}).forEach(([id, node]) => {
      nodeScoresB[id] = Math.exp(node.log_prob || 0);
    });
  });
  const union = new Set([...Object.keys(nodeScoresA), ...Object.keys(nodeScoresB)]);
  union.forEach((id) => {
    const aScore = nodeScoresA[id] || 0;
    const bScore = nodeScoresB[id] || 0;
    const diff = aScore - bScore;
    const group = !aScore ? "B" : !bScore ? "A" : "Both";
    comparisonCy.add({
      group: "nodes",
      data: { id, label: `${id}\nΔ${diff.toFixed(3)}`, group, diff },
    });
  });

  document.getElementById("rt-compare-output").innerHTML = `
    <div>Nodes only A: ${
      Object.keys(nodeScoresA).filter((id) => !nodeScoresB[id]).length
    }</div>
    <div>Nodes only B: ${
      Object.keys(nodeScoresB).filter((id) => !nodeScoresA[id]).length
    }</div>
    <div>Shared: ${
      [...union].filter((id) => nodeScoresA[id] && nodeScoresB[id]).length
    }</div>
  `;
  comparisonCy.layout({ name: "cose", animate: false }).run();
}

function exportComparisonYAML() {
  const doc = {
    kind: "helix.compare.v1",
    genome: document.getElementById("rt-genome").value,
    pam: document.getElementById("rt-compare-pam").value,
    mismatch_tolerance: Number(document.getElementById("rt-compare-mismatch").value),
    guide_a: document.getElementById("rt-compare-guide-a").value,
    guide_b: document.getElementById("rt-compare-guide-b").value,
  };
  const blob = new Blob([JSON.stringify(doc, null, 2)], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "helix_realtime_compare.json";
  link.click();
  URL.revokeObjectURL(link.href);
}
worker.addEventListener("message", (event) => {
  const { type, frame, error } = event.data || {};
  if (type === "frame") {
    frames.push(frame);
    document.getElementById("rt-slider").max = frames.length - 1;
  } else if (type === "error") {
    alert(error || "Simulation failed.");
  } else if (type === "done") {
    if (frames.length) {
      loadFrameAt(0);
    }
  }
});

document.addEventListener("DOMContentLoaded", () => {
  initCytoscape();
  branchChart = initBranchChart(document.getElementById("rt-branch-chart"));
  document.getElementById("rt-run").addEventListener("click", runSimulation);
  document.getElementById("rt-export").addEventListener("click", exportExperiment);
  document.getElementById("rt-play").addEventListener("click", playFrames);
  document.getElementById("rt-slider").addEventListener("input", (event) => {
    const idx = Number(event.target.value);
    loadFrameAt(idx);
  });
  document.getElementById("rt-reset").addEventListener("click", resetPlayback);
  document.getElementById("rt-compare-run").addEventListener("click", runComparison);
  document.getElementById("rt-compare-export").addEventListener("click", exportComparisonYAML);

  // Path highlighting: click node to highlight its predecessors path
  cytoscapeInstance.on("tap", "node", (evt) => {
    const sel = evt.target;
    const ancestors = sel.predecessors().add(sel);
    cytoscapeInstance.elements().addClass("faded");
    ancestors.removeClass("faded");
    const details = document.getElementById("rt-details");
    const stage = sel.data("stage");
    const prob = sel.data("prob") || Math.exp(sel.data("log_prob") || 0);
    const time = sel.data("time");
    const ev = sel.data("event");
    const eff = effectLabel(ev);
    details.innerHTML = `
      <div>Stage: ${stage || "?"}</div>
      <div>Time step: ${time ?? "?"}</div>
      <div>Probability: ${prob.toFixed(4)}</div>
      <div>Desired: ${sel.data("desired") === true ? "yes" : sel.data("desired") === false ? "no" : "—"}</div>
      ${eff ? `<div>Effect: ${eff}</div>` : ""}
      ${sel.data("raw_meta")?.protein_impact ? `<div>Protein impact: ${sel.data("raw_meta").protein_impact}</div>` : ""}
    `;
    // sequence strip (simple window around event)
    const seqBox = document.getElementById("rt-seq");
    const sequences = sel.data("sequences") || {};
    const seq = rootChrom ? sequences[rootChrom] : null;
    if (rootSeq && seq && ev && ev.start !== undefined && ev.end !== undefined) {
      const s = Math.max(0, ev.start - 15);
      const e = Math.min(rootSeq.length, (ev.end || ev.start) + 15);
      const refWin = rootSeq.slice(s, e);
      const altWin = seq.slice(s, e);
      const a = refWin.split("");
      const b = altWin.split("");
      const outA = [];
      const outB = [];
      const n = Math.max(a.length, b.length);
      for (let i = 0; i < n; i++) {
        const ca = a[i] || "";
        const cb = b[i] || "";
        if (ca !== cb) {
          outA.push(`[${ca}]`);
          outB.push(`[${cb}]`);
        } else {
          outA.push(ca);
          outB.push(cb);
        }
      }
      seqBox.textContent = `ref ${rootChrom}:${s+1}-${e}\n` + outA.join("") + "\nalt\n" + outB.join("");
    } else {
      seqBox.textContent = "";
    }
  });
  cytoscapeInstance.on("tap", (evt) => {
    if (evt.target === cytoscapeInstance) {
      cytoscapeInstance.elements().removeClass("faded");
    }
  });

  // Hover tooltip for nodes with internal ID and metadata
  const tip = document.getElementById("rt-tooltip");
  function renderTip(ele) {
    const d = ele.data();
    const prob = d.prob ?? Math.exp(d.log_prob || 0);
    const meta = d.raw_meta || {};
    const event = d.event || null;
    const kv = {
      id: d.id,
      stage: d.stage,
      desired: typeof d.desired === "boolean" ? d.desired : undefined,
      on_target: typeof d.on_target === "boolean" ? d.on_target : undefined,
      time_step: d.time,
    };
    const pre = document.createElement("pre");
    pre.style.margin = "6px 0 0";
    pre.style.whiteSpace = "pre-wrap";
    const metaShown = JSON.stringify(meta, null, 2);
    const eventShown = event ? JSON.stringify(event, null, 2) : null;
    tip.innerHTML = `<div style="font-weight:600">${d.label.split("\n")[0]}</div>
      <div style="opacity:.8">p=${prob.toFixed(4)}</div>`;
    if (eventShown) {
      const evDiv = document.createElement("div");
      evDiv.style.marginTop = "6px";
      evDiv.textContent = `event: ${eventShown}`;
      tip.appendChild(evDiv);
    }
    const idDiv = document.createElement("div");
    idDiv.style.marginTop = "6px";
    idDiv.textContent = `id: ${kv.id}`;
    tip.appendChild(idDiv);
    pre.textContent = metaShown;
    tip.appendChild(pre);
  }
  function moveTip(clientX, clientY) {
    const pad = 12;
    tip.style.left = `${clientX + pad}px`;
    tip.style.top = `${clientY + pad}px`;
  }
  cytoscapeInstance.on("mouseover", "node", (evt) => {
    renderTip(evt.target);
    tip.style.display = "block";
    const e = evt.originalEvent || {};
    moveTip(e.clientX || 0, e.clientY || 0);
  });
  cytoscapeInstance.on("mouseout", "node", () => {
    tip.style.display = "none";
  });
  cytoscapeInstance.on("mousemove", "node", (evt) => {
    const e = evt.originalEvent || {};
    moveTip(e.clientX || 0, e.clientY || 0);
  });

  // Outcome filters
  function applyOutcomeFilters() {
    const showIntended = document.getElementById("filter-intended")?.checked ?? true;
    const showError = document.getElementById("filter-error")?.checked ?? true;
    const showNoEdit = document.getElementById("filter-noedit")?.checked ?? true;
    const showOther = document.getElementById("filter-other")?.checked ?? true;
    const showOnTarget = document.getElementById("filter-ontarget")?.checked ?? true;
    const showOffTarget = document.getElementById("filter-offtarget")?.checked ?? true;
    const showDesired = document.getElementById("filter-desired")?.checked ?? true;
    const showUndesired = document.getElementById("filter-undesired")?.checked ?? true;
    const onlyDesired = document.getElementById("filter-only-desired")?.checked ?? false;
    const onlyTerminal = document.getElementById("filter-only-terminal")?.checked ?? false;
    const keepStage = (stage) => {
      if (stage === "repaired") return showIntended;
      if (stage === "error") return showError;
      if (stage === "no_edit") return showNoEdit;
      return showOther; // root, cut, prime_rtt, flap_*
    };
    const nodes = cytoscapeInstance.nodes();
    nodes.forEach((n) => {
      const stage = n.data("stage");
      const desiredFlag = n.data("desired");
      const desiredOk = onlyDesired
        ? desiredFlag === true
        : (desiredFlag === true && showDesired) || (desiredFlag === false && showUndesired) || (typeof desiredFlag === "undefined" && ((stage === "repaired" && showDesired) || (stage === "error" && showUndesired) || (stage !== "repaired" && stage !== "error" && (showDesired || showUndesired))));
      const onTargetFlag = n.data("on_target");
      const onTargetOk = (onTargetFlag === true && showOnTarget) || (onTargetFlag === false && showOffTarget) || (typeof onTargetFlag === "undefined" && (showOnTarget || showOffTarget));
      const isTerminal = n.outgoers('edge').length === 0;
      const keep = (onlyDesired ? true : keepStage(stage)) && desiredOk && onTargetOk && (!onlyTerminal || isTerminal);
      if (keep) n.removeClass("hidden");
      else n.addClass("hidden");
    });
    const edges = cytoscapeInstance.edges();
    edges.forEach((e) => {
      const keep = !e.source().hasClass("hidden") && !e.target().hasClass("hidden");
      if (keep) e.removeClass("hidden");
      else e.addClass("hidden");
    });
  }
  // expose to outer scope for applyFrame hook
  window.applyOutcomeFilters = applyOutcomeFilters;
  ["filter-intended", "filter-error", "filter-noedit", "filter-other", "filter-desired", "filter-undesired", "filter-ontarget", "filter-offtarget", "filter-only-desired", "filter-only-terminal"].forEach((id) => {
    const el = document.getElementById(id);
    if (el) el.addEventListener("change", applyOutcomeFilters);
  });

  // Transcript JSON loader + selector + exporter
  const fileInput = document.getElementById("rt-transcript-file");
  const select = document.getElementById("rt-transcript-select");
  let loadedTranscripts = [];
  function normalizeTranscripts(obj) {
    const list = [];
    if (!obj) return list;
    const arr = Array.isArray(obj.transcripts) ? obj.transcripts : (Array.isArray(obj) ? obj : []);
    for (const t of arr) {
      const name = t.name || t.id || `tx_${list.length + 1}`;
      const strand = (t.strand || "+").toString();
      const cds_start = Number(t.cds_start || t.cdsStart || 0);
      const cds_end = Number(t.cds_end || t.cdsEnd || 0);
      const exons = (t.exons || [])
        .map((e) => ({ start1: Number(e[0] ?? e.start ?? e.start1), end1: Number(e[1] ?? e.end ?? e.end1) }))
        .filter((e) => !Number.isNaN(e.start1) && !Number.isNaN(e.end1));
      if (cds_start > 0 && cds_end > 0 && exons.length) {
        list.push({ id: t.id || name, name, strand, cds_start, cds_end, exons });
      }
    }
    return list;
  }
  function populateTranscriptSelect(list) {
    if (!select) return;
    select.innerHTML = "";
    const opt0 = document.createElement("option");
    opt0.value = "";
    opt0.textContent = "(choose transcript)";
    select.appendChild(opt0);
    list.forEach((t, idx) => {
      const opt = document.createElement("option");
      opt.value = String(idx);
      opt.textContent = t.name;
      select.appendChild(opt);
    });
    select.disabled = list.length === 0;
  }
  function applyTranscriptToFields(t) {
    if (!t) return;
    document.getElementById("rt-cds-start").value = String(t.cds_start);
    document.getElementById("rt-cds-end").value = String(t.cds_end);
    document.getElementById("rt-cds-strand").value = t.strand === "-" ? "-" : "+";
    document.getElementById("rt-exons").value = t.exons.map((e) => `${e.start1}-${e.end1}`).join(",");
  }
  if (fileInput) {
    fileInput.addEventListener("change", (ev) => {
      const f = ev.target.files?.[0];
      if (!f) return;
      const rd = new FileReader();
      rd.onload = () => {
        try {
          const obj = JSON.parse(String(rd.result || "{}"));
          loadedTranscripts = normalizeTranscripts(obj);
          populateTranscriptSelect(loadedTranscripts);
        } catch (e) {
          alert("Failed to parse transcript JSON: " + e.message);
        }
      };
      rd.readAsText(f);
    });
  }
  if (select) {
    select.addEventListener("change", () => {
      const idx = Number(select.value);
      if (!Number.isNaN(idx) && loadedTranscripts[idx]) applyTranscriptToFields(loadedTranscripts[idx]);
    });
  }
  const exportBtn = document.getElementById("rt-transcript-export");
  if (exportBtn) {
    exportBtn.addEventListener("click", () => {
      const cdsStart = Number(document.getElementById("rt-cds-start").value || 0);
      const cdsEnd = Number(document.getElementById("rt-cds-end").value || 0);
      const strand = document.getElementById("rt-cds-strand").value || "+";
      const exonsText = document.getElementById("rt-exons").value || "";
      if (!(cdsStart > 0 && cdsEnd > 0)) {
        alert("Please fill CDS start/end to export.");
        return;
      }
      const exons = (exonsText || "")
        .split(/[,\s]+/)
        .filter(Boolean)
        .map((t) => t.split("-").map((x) => Number(x)))
        .filter((p) => p.length === 2 && !Number.isNaN(p[0]) && !Number.isNaN(p[1]));
      const tx = {
        transcripts: [
          { id: "TX1", name: "TX1", strand, cds_start: cdsStart, cds_end: cdsEnd, exons },
        ],
      };
      const blob = new Blob([JSON.stringify(tx, null, 2)], { type: "application/json" });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = "transcript.json";
      a.click();
      URL.revokeObjectURL(a.href);
    });
  }
});
