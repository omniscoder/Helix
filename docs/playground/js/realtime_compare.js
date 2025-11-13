/* global cytoscape */

function initCompareSimulation(container, workerPath) {
  const worker = new Worker(workerPath, { type: "module" });
  const frames = [];
  let cy = cytoscape({
    container,
    style: [
      {
        selector: "node",
        style: {
          shape: "round-rectangle",
          "background-color": "#6366f1",
          label: "data(label)",
          color: "#e2e8f0",
          "font-size": 10,
        },
      },
      {
        selector: "edge",
        style: {
          "line-color": "#94a3b8",
          "target-arrow-color": "#94a3b8",
          "target-arrow-shape": "triangle",
          width: 1.5,
        },
      },
    ],
  });

  worker.addEventListener("message", (event) => {
    const { type, frame, error } = event.data || {};
    if (type === "frame") {
      frames.push(frame);
      applyFrame(cy, frame);
    } else if (type === "error") {
      // eslint-disable-next-line no-alert
      alert(error || "Realtime compare failed.");
    }
  });

  function applyFrame(cyto, frameData) {
    Object.entries(frameData.new_nodes || {}).forEach(([nodeId, node]) => {
      if (cyto.$id(nodeId).length) return;
      cyto.add({
        group: "nodes",
        data: {
          id: nodeId,
          label: `${node.metadata?.stage || "node"}\n${nodeId}`,
        },
      });
    });
    (frameData.new_edges || []).forEach((edge) => {
      const edgeId = `${edge.source}_${edge.target}_${edge.rule}`;
      if (cyto.$id(edgeId).length) return;
      cyto.add({
        group: "edges",
        data: {
          id: edgeId,
          source: edge.source,
          target: edge.target,
        },
      });
    });
    cyto.layout({ name: "cose", animate: false }).run();
  }

  return {
    run(payload) {
      frames.length = 0;
      cy.elements().remove();
      worker.postMessage({ type: "simulate", payload });
    },
  };
}

export { initCompareSimulation };
