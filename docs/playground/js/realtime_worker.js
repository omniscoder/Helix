/* eslint-disable no-restricted-globals */

function parseGenome(text) {
  const lines = text.trim().split(/\r?\n/);
  const seq = [];
  for (const line of lines) {
    if (!line.trim()) continue;
    if (line.startsWith(">")) continue;
    seq.push(line.trim().toUpperCase());
  }
  return seq.join("");
}

function pamMatches(seq, idx, pam) {
  if (!pam) return true;
  const chars = pam.toUpperCase();
  for (let i = 0; i < chars.length; i++) {
    const sIdx = idx + i;
    if (sIdx >= seq.length) return false;
    const ch = seq[sIdx];
    const rule = chars[i];
    const table = {
      R: "AG",
      Y: "CT",
      S: "GC",
      W: "AT",
      K: "GT",
      M: "AC",
      H: "ACT",
      B: "CGT",
      V: "ACG",
      D: "AGT",
    };
    if (rule === "N") continue;
    if (table[rule]) {
      if (!table[rule].includes(ch)) return false;
      continue;
    }
    if (rule !== ch) return false;
  }
  return true;
}

function mismatchCount(a, b) {
  let mismatches = 0;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) {
    if (a[i] !== b[i]) mismatches += 1;
  }
  return mismatches + Math.abs(a.length - b.length);
}

function reverseComplement(seq) {
  const map = { A: "T", T: "A", C: "G", G: "C" };
  return seq
    .split("")
    .reverse()
    .map((ch) => map[ch] || ch)
    .join("");
}

function computeFrames(payload) {
  const genome = parseGenome(payload.genome || "");
  const chrom = payload.chrom || "chrDemo";
  const guide = (payload.guide || "").trim().toUpperCase();
  const pam = (payload.pam || "NGG").trim().toUpperCase();
  const mismatchTolerance = Number(payload.mismatchTolerance || 3);
  const maxSites = Number(payload.maxSites || 5);
  const window = Number(payload.window || 5);

  if (!genome || !guide) {
    throw new Error("Genome sequence and guide are required.");
  }

  const frames = [];
  const rootId = "n_root";
  frames.push({
    kind: "helix.edit_dag.frame.v1",
    step: 0,
    new_nodes: {
      [rootId]: {
        log_prob: 0,
        metadata: { stage: "root", time_step: 0 },
        parent_ids: [],
        seq_hashes: { [chrom]: genome },
        sequences: { [chrom]: genome },
      },
    },
    new_edges: [],
    meta: { mechanism: "crispr" },
  });

  const candidates = [];
  for (let pos = 0; pos <= genome.length - guide.length; pos++) {
    const target = genome.slice(pos, pos + guide.length);
    const mismatches = mismatchCount(target, guide);
    if (mismatches > mismatchTolerance) continue;
    const pamStart = pos + guide.length;
    if (pamStart + pam.length > genome.length) continue;
    if (!pamMatches(genome, pamStart, pam)) continue;
    const score = Math.max(1e-6, 1 - mismatches / guide.length);
    candidates.push({ pos, target, score });
  }

  candidates.sort((a, b) => b.score - a.score);
  const selected = candidates.slice(0, maxSites);

  selected.forEach((site, idx) => {
    const cutId = `n_cut_${idx}`;
    frames.push({
      kind: "helix.edit_dag.frame.v1",
      step: frames.length,
      new_nodes: {
        [cutId]: {
          log_prob: Math.log(site.score),
          metadata: { stage: "cut", time_step: 1, position: site.pos },
          parent_ids: [rootId],
          seq_hashes: { [chrom]: genome },
          sequences: { [chrom]: genome },
        },
      },
      new_edges: [
        {
          source: rootId,
          target: cutId,
          rule: "crispr.clean_cut",
          event: {
            chrom,
            start: site.pos,
            end: site.pos,
            replacement: "",
            metadata: {},
          },
          metadata: {},
        },
      ],
      meta: { mechanism: "crispr" },
    });

    const intendedSeq = `${genome.slice(0, site.pos)}${reverseComplement(site.target)}${genome.slice(site.pos + site.target.length)}`;
    const indelSeq = `${genome.slice(0, site.pos)}${genome.slice(Math.min(genome.length, site.pos + window))}`;
    const intendedId = `${cutId}_intended`;
    const indelId = `${cutId}_indel`;
    frames.push({
      kind: "helix.edit_dag.frame.v1",
      step: frames.length,
      new_nodes: {
        [intendedId]: {
          log_prob: Math.log(site.score * 0.7),
          metadata: { stage: "repaired", time_step: 2, branch: "intended" },
          parent_ids: [cutId],
          seq_hashes: { [chrom]: intendedSeq },
          sequences: { [chrom]: intendedSeq },
        },
        [indelId]: {
          log_prob: Math.log(site.score * 0.3 + 1e-6),
          metadata: { stage: "error", time_step: 2, branch: "indel" },
          parent_ids: [cutId],
          seq_hashes: { [chrom]: indelSeq },
          sequences: { [chrom]: indelSeq },
        },
      },
      new_edges: [
        {
          source: cutId,
          target: intendedId,
          rule: "crispr.indel_branch",
          event: {
            chrom,
            start: site.pos,
            end: site.pos + site.target.length,
            replacement: reverseComplement(site.target),
            metadata: { branch: "intended" },
          },
          metadata: {},
        },
        {
          source: cutId,
          target: indelId,
          rule: "crispr.indel_branch",
          event: {
            chrom,
            start: site.pos,
            end: Math.min(genome.length, site.pos + window),
            replacement: "",
            metadata: { branch: "indel" },
          },
          metadata: {},
        },
      ],
      meta: { mechanism: "crispr" },
    });
  });

  return frames;
}

self.addEventListener("message", (event) => {
  const { type, payload } = event.data || {};
  if (type !== "simulate") {
    self.postMessage({ type: "error", error: "Unknown worker request." });
    return;
  }
  try {
    const frames = computeFrames(payload);
    frames.forEach((frame) => self.postMessage({ type: "frame", frame }));
    self.postMessage({ type: "done", frameCount: frames.length });
  } catch (error) {
    self.postMessage({ type: "error", error: error.message || String(error) });
  }
});
