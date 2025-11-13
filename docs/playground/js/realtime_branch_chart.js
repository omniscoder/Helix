export function initBranchChart(canvas) {
  return new Chart(canvas, {
    type: "bar",
    data: {
      labels: ["Intended", "Indel", "No Edit"],
      datasets: [
        {
          label: "Probability",
          data: [0, 0, 0],
          backgroundColor: ["#10b981", "#f97316", "#64748b"],
        },
      ],
    },
    options: {
      responsive: true,
      animation: false,
      plugins: { legend: { display: false } },
      scales: {
        y: {
          beginAtZero: true,
          suggestedMax: 1,
          ticks: { color: "#94a3b8" },
          grid: { color: "#1e293b" },
        },
        x: {
          ticks: { color: "#94a3b8" },
          grid: { color: "#1e293b" },
        },
      },
    },
  });
}

export function updateBranchChart(chart, intended, indel, noEdit) {
  if (!chart) return;
  chart.data.datasets[0].data = [intended, indel, noEdit];
  chart.update();
}
