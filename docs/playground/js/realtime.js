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

function initCytoscape() {
  cytoscapeInstance = cytoscape({
    container: document.getElementById("rt-cy"),
    style: [
      {
        selector: "node",
        style: {
          "shape": "round-rectangle",
          "background-color": "#3b82f6",
          "width": 40,
          "height": 20,
          "label": "data(label)",
          "color": "#e5e7eb",
          "font-size": 10,
          "text-wrap": "wrap",
        },
      },
      {
        selector: "edge",
        style: {
          "curve-style": "bezier",
          "line-color": "#94a3b8",
          "target-arrow-color": "#94a3b8",
          "target-arrow-shape": "triangle",
          "width": 1.5,
          "label": "data(rule)",
          "font-size": 8,
          "color": "#94a3b8",
          "text-rotation": "autorotate",
        },
      },
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
  Object.entries(frame.new_nodes || {}).forEach(([nodeId, node]) => {
    if (cy.$id(nodeId).length) return;
    const label = `${node.metadata?.stage || "node"}\n${nodeId}`;
    cy.add({
      group: "nodes",
      data: { id: nodeId, label, stage: node.metadata?.stage, log_prob: node.log_prob },
    });
  });
  (frame.new_edges || []).forEach((edge) => {
    const edgeId = `${edge.source}_${edge.target}_${edge.rule}`;
    if (cy.$id(edgeId).length) return;
    cy.add({
      group: "edges",
      data: {
        id: edgeId,
        source: edge.source,
        target: edge.target,
        rule: edge.rule,
      },
    });
  });
  cy.layout({ name: "cose", animate: false, padding: 20 }).run();
  updateMetrics(frame);
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
});
