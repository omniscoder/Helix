# Benchmarks

Helix ships a JSON benchmark harness (`python -m benchmarks.api_benchmarks`) and continuously publishes the results to CI. The goals are:

- **Comparable runs:** every payload is tagged with commit SHA, dataset provenance, CPU/threads, BLAS vendor, RNG seed, and per-case RSS stats.
- **Drift detection:** CI rejects PRs when a case slows down by more than 5 % relative to `.bench/baseline.json`.
- **Transparency:** the latest measurements are summarized here, powered by the CSV in `docs/data/bench/history.csv` (populated automatically on `main`).

## Running locally

```bash
python -m benchmarks.api_benchmarks \
  --repeat 5 \
  --warmup 1 \
  --limit 0 \
  --out bench-results/api.json \
  --summary-md bench-results/api.md
```

- `--limit N` keeps only the first `N` nucleotides/aminos (0 = entire dataset). Use this for quick inner-loop runs or to mimic CI’s 10k-sample sweep.
- Override datasets with `HELIX_BENCH_DNA_FASTA=/abs/path/...` and `HELIX_BENCH_PROTEIN_FASTA=/abs/path/...`. The harness records those paths in the JSON payload so dashboards can compare apples-to-apples.
- Pass `--baseline path/to/baseline.json` to compute Δ% vs. a stored run. `scripts/bench_check.py baseline current --threshold 5` is what CI uses to gate regressions.

## Heavy datasets via GitHub Actions

Trigger a manual heavy sweep from the **Actions → CI → Run workflow** button:

1. Set `bench_heavy` to `true` (this bumps repeats to 10 and disables the 10k sampling limit).
2. Optionally provide runner-accessible overrides for `dna_fasta` / `protein_fasta`. On hosted runners you typically leave these blank; on self-hosted boxes you can point at a mounted volume or fetcher script.

Each run publishes:

- `benchmarks/out/bench-<SHA>.json` —  the full schema payload.
- `benchmarks/out/bench-<SHA>.md` — a Markdown table appended to the CI summary.
- `docs/data/bench/history.csv` (main branch only) — an append-only log that powers the chart below.

## Trend (mean seconds)

The gallery below visualizes every `*.mean_s` column recorded in `docs/data/bench/history.csv` and summarizes the latest run.

<style>
.bench-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
  gap: 1rem;
  margin-bottom: 1rem;
}
.bench-card {
  border: 1px solid #e5e7eb;
  border-radius: 0.75rem;
  padding: 0.75rem;
  background: #fff;
  box-shadow: 0 1px 2px rgba(0,0,0,0.05);
}
.bench-card canvas {
  width: 100%;
  height: 200px;
}
.bench-table {
  width: 100%;
  border-collapse: collapse;
  margin-top: 0.75rem;
}
.bench-table th,
.bench-table td {
  border: 1px solid #e5e7eb;
  padding: 0.4rem 0.6rem;
  text-align: left;
  font-size: 0.9rem;
}
.bench-summary {
  margin-top: 1rem;
}
</style>
<script>
(async () => {
  async function loadCSV(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Failed to load ${url}`);
    const text = await resp.text();
    const lines = text.trim().split(/\r?\n/).filter(Boolean);
    if (!lines.length) return [];
    const header = lines[0].split(',');
    return lines.slice(1).map(line => {
      const cols = line.split(',');
      const row = {};
      header.forEach((key, idx) => { row[key] = cols[idx] ?? ''; });
      return row;
    });
  }

  function drawSeries(canvas, label, rows) {
    const ctx = canvas.getContext('2d');
    const values = rows.map(r => parseFloat(r[label])).filter(v => !Number.isNaN(v));
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (!values.length) {
      ctx.fillStyle = '#6b7280';
      ctx.font = 'bold 14px sans-serif';
      ctx.fillText('No data', 10, 20);
      return;
    }
    const min = Math.min(...values);
    const max = Math.max(...values);
    const left = 40;
    const bottom = canvas.height - 20;
    const top = 20;
    const plotWidth = canvas.width - left - 20;
    const plotHeight = bottom - top;

    ctx.strokeStyle = '#d1d5db';
    ctx.beginPath();
    ctx.moveTo(left, top);
    ctx.lineTo(left, bottom);
    ctx.lineTo(canvas.width - 20, bottom);
    ctx.stroke();

    ctx.strokeStyle = '#2563eb';
    ctx.beginPath();
    values.forEach((v, idx) => {
      const x = left + (values.length <= 1 ? 0 : (idx / (values.length - 1)) * plotWidth);
      const y = bottom - ((v - min) / ((max - min) || 1)) * plotHeight;
      if (idx === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();

    ctx.fillStyle = '#111827';
    ctx.font = '12px sans-serif';
    ctx.fillText(`min=${min.toFixed(4)}s max=${max.toFixed(4)}s`, left, top - 4);
  }

  let rows = [];
  try {
    rows = await loadCSV('../data/bench/history.csv');
  } catch (err) {
    console.warn('Unable to load benchmark history:', err);
  }
  const metrics = rows.length ? Object.keys(rows[0]).filter(k => k.endsWith('.mean_s')) : [];
  if (!rows.length || !metrics.length) {
    const note = document.createElement('p');
    note.textContent = 'Benchmark history not available yet.';
    document.currentScript.insertAdjacentElement('afterend', note);
    return;
  }

  const grid = document.createElement('div');
  grid.className = 'bench-grid';
  document.currentScript.insertAdjacentElement('afterend', grid);

  metrics.forEach(label => {
    const card = document.createElement('div');
    card.className = 'bench-card';
    const title = document.createElement('strong');
    title.textContent = label.replace('.mean_s', '');
    const canvas = document.createElement('canvas');
    canvas.width = 320;
    canvas.height = 220;
    card.appendChild(title);
    card.appendChild(canvas);
    grid.appendChild(card);
    drawSeries(canvas, label, rows);
  });

  const latest = rows[rows.length - 1];
  const datasetLabel = latest.dataset ? latest.dataset.split('/').pop() : 'default';
  const summary = document.createElement('div');
  summary.className = 'bench-summary';
  const info = document.createElement('p');
  info.textContent = `Latest run • sha=${latest.sha} • dataset=${datasetLabel} • limit=${latest.limit || 'full'} • repeat=${latest.repeat}`;
  summary.appendChild(info);

  const table = document.createElement('table');
  table.className = 'bench-table';
  const headerRow = document.createElement('tr');
  ['Case', 'Mean (s)', 'Δ vs baseline', 'RSS peak (MB)'].forEach(text => {
    const th = document.createElement('th');
    th.textContent = text;
    headerRow.appendChild(th);
  });
  table.appendChild(headerRow);

  metrics.forEach(label => {
    const mean = parseFloat(latest[label] || 'NaN');
    const delta = latest[label.replace('.mean_s', '.delta_pct')];
    const rss = latest[label.replace('.mean_s', '.rss_mb')];
    const rowEl = document.createElement('tr');
    const cells = [
      label.replace('.mean_s', ''),
      Number.isNaN(mean) ? '—' : mean.toFixed(4),
      delta ? `${parseFloat(delta).toFixed(1)}%` : '—',
      rss ? parseFloat(rss).toFixed(1) : '—',
    ];
    cells.forEach(text => {
      const td = document.createElement('td');
      td.textContent = text;
      rowEl.appendChild(td);
    });
    table.appendChild(rowEl);
  });

  summary.appendChild(table);
  grid.insertAdjacentElement('afterend', summary);
})();
</script>

> CSV source: `docs/data/bench/history.csv`. Commit history contains the raw JSON artifacts under `benchmarks/out/` (uploaded by CI) if you need to recompute metrics offline.
