/* global cytoscape */

const worker = new Worker("js/realtime_worker.js", { type: "module" });
const frames = [];
let cytoscapeInstance;
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

function updateMetrics(frame) {
  const metricsBox = document.getElementById("rt-metrics");
  const nodes = cytoscapeInstance.nodes();
  const entropy = nodes
    .toArray()
    .reduce((acc, node) => {
      const prob = Math.exp(node.data("log_prob") || 0);
      return acc - (prob ? prob * Math.log(prob) : 0);
    }, 0);
  metricsBox.innerHTML = `
    <div>Frame: ${frame.step}</div>
    <div>Nodes: ${nodes.length}</div>
    <div>Entropy: ${entropy.toFixed(3)}</div>
  `;
}

function resetPlayback() {
  frames.length = 0;
  currentIndex = 0;
  document.getElementById("rt-slider").value = 0;
  cytoscapeInstance.elements().remove();
  document.getElementById("rt-metrics").textContent = "Load a simulation…";
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
  document.getElementById("rt-run").addEventListener("click", runSimulation);
  document.getElementById("rt-play").addEventListener("click", playFrames);
  document.getElementById("rt-slider").addEventListener("input", (event) => {
    const idx = Number(event.target.value);
    loadFrameAt(idx);
  });
  document.getElementById("rt-reset").addEventListener("click", resetPlayback);
});
