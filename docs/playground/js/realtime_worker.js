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
  // If the payload already contains edit DAG frames, treat them as the
  // canonical source of truth and avoid re-implementing CRISPR logic here.
  if (payload && Array.isArray(payload.frames)) {
    return payload.frames
      .filter((frame) => frame && typeof frame === "object")
      .map((frame, idx) => {
        const copy = { ...frame };
        if (!copy.kind) copy.kind = "helix.edit_dag.frame.v1";
        if (typeof copy.step !== "number") copy.step = idx;
        return copy;
      });
  }

  const genomeInput = payload.genome || "";
  let genome = "";
  if (genomeInput.startsWith(">") || genomeInput.includes("\n")) {
    genome = parseGenome(genomeInput);
  } else {
    genome = genomeInput.trim().toUpperCase();
  }
  const chrom = payload.chrom || "chrDemo";
  const guide = (payload.guide || "").trim().toUpperCase();
  const pam = (payload.pam || "NGG").trim().toUpperCase();
  const mismatchTolerance = Number(payload.mismatchTolerance || 3);
  const maxSites = Number(payload.maxSites || 5);
  const window = Number(payload.window || 5);
  const cdsStart1 = payload.cdsStart ? Number(payload.cdsStart) : null; // 1-based
  const cdsEnd1 = payload.cdsEnd ? Number(payload.cdsEnd) : null; // 1-based
  const cdsStrand = (payload.cdsStrand || "+").trim();
  const exonsText = (payload.exons || "").trim();
  function parseExons(text) {
    if (!text) return null;
    const items = text.split(/[,\s]+/).filter(Boolean);
    const out = [];
    for (const it of items) {
      const m = it.match(/^(\d+)-(\d+)$/);
      if (m) {
        const a = Number(m[1]);
        const b = Number(m[2]);
        if (!Number.isNaN(a) && !Number.isNaN(b)) out.push({ start1: Math.min(a, b), end1: Math.max(a, b) });
      }
    }
    return out.length ? out : null;
  }
  const exons = parseExons(exonsText);

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

  // mark top-scoring as on-target for demo; rest off-target
  const onTargetIdx = selected.length > 0 ? 0 : -1;
  selected.forEach((site, idx) => {
    const cutId = `n_cut_${idx}`;
    frames.push({
      kind: "helix.edit_dag.frame.v1",
      step: frames.length,
      new_nodes: {
        [cutId]: {
          log_prob: Math.log(site.score),
          metadata: { stage: "cut", time_step: 1, position: site.pos, on_target: idx === onTargetIdx, score: site.score },
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
          metadata: { stage: "repaired", time_step: 2, branch: "intended", on_target: idx === onTargetIdx },
          parent_ids: [cutId],
          seq_hashes: { [chrom]: intendedSeq },
          sequences: { [chrom]: intendedSeq },
        },
        [indelId]: {
          log_prob: Math.log(site.score * 0.3 + 1e-6),
          metadata: { stage: "error", time_step: 2, branch: "indel", on_target: idx === onTargetIdx },
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

  // annotate terminal nodes with protein impact if CDS provided and SNV
  if (cdsStart1 != null) {
    const code = {
      TTT:"F",TTC:"F",TTA:"L",TTG:"L",TCT:"S",TCC:"S",TCA:"S",TCG:"S",
      TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",TGT:"C",TGC:"C",TGA:"*",TGG:"W",
      CTT:"L",CTC:"L",CTA:"L",CTG:"L",CCT:"P",CCC:"P",CCA:"P",CCG:"P",
      CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",CGT:"R",CGC:"R",CGA:"R",CGG:"R",
      ATT:"I",ATC:"I",ATA:"I",ATG:"M",ACT:"T",ACC:"T",ACA:"T",ACG:"T",
      AAT:"N",AAC:"N",AAA:"K",AAG:"K",AGT:"S",AGC:"S",AGA:"R",AGG:"R",
      GTT:"V",GTC:"V",GTA:"V",GTG:"V",GCT:"A",GCC:"A",GCA:"A",GCG:"A",
      GAT:"D",GAC:"D",GAA:"E",GAG:"E",GGT:"G",GGC:"G",GGA:"G",GGG:"G"
    };
    // Build coding map if exons and cds start/end provided
    let codingMap = null; // array of genome pos (0-based) in transcript 5'->3' order
    if (exons && cdsEnd1 != null) {
      const codingBlocks = [];
      for (const ex of exons.sort((a,b)=>a.start1-b.start1)) {
        const s = Math.max(ex.start1, cdsStart1);
        const e = Math.min(ex.end1, cdsEnd1);
        if (s <= e) codingBlocks.push({ start1: s, end1: e });
      }
      codingMap = [];
      if (cdsStrand === '+') {
        for (const bl of codingBlocks) {
          for (let p = bl.start1 - 1; p <= bl.end1 - 1; p++) codingMap.push(p);
        }
      } else {
        for (let i = codingBlocks.length - 1; i >= 0; i--) {
          const bl = codingBlocks[i];
          for (let p = bl.end1 - 1; p >= bl.start1 - 1; p--) codingMap.push(p);
        }
      }
    }

    function proteinImpactForEvent(ev, refSeq) {
      if (!ev) return null;
      const len = Math.max(0, (ev.end||0)-(ev.start||0));
      const ins = (ev.replacement||"").length;
      if (!(len===1 && ins===1)) return null; // SNV only
      const pos0 = ev.start; // 0-based genomic position
      let refCodon = null, altCodon = null;
      if (codingMap && codingMap.length >= 3) {
        const tIdx = codingMap.indexOf(pos0);
        if (tIdx < 0) return null; // non-coding
        const frame = tIdx % 3;
        const tCodonStart = tIdx - frame;
        if (tCodonStart + 2 >= codingMap.length) return null;
        function getBaseAtGenome(gp){ return refSeq[gp] || 'N'; }
        function rcBase(b){ const map={A:'T',T:'A',C:'G',G:'C'}; return map[b]||b; }
        const bases = [codingMap[tCodonStart], codingMap[tCodonStart+1], codingMap[tCodonStart+2]].map(getBaseAtGenome);
        if (cdsStrand === '-') {
          // convert to transcript orientation
          refCodon = bases.map(rcBase).reverse().join('');
          const idxInCodon = frame; // since reversed
          const repRC = rcBase((ev.replacement||'N'));
          const arr = refCodon.split('');
          arr[idxInCodon] = repRC;
          altCodon = arr.join('');
        } else {
          refCodon = bases.join('');
          const idxInCodon = frame;
          const arr = refCodon.split('');
          arr[idxInCodon] = (ev.replacement||'N');
          altCodon = arr.join('');
        }
      } else {
        // fallback to simple CDS start method
        if (cdsStrand === '+') {
          const offset = pos0 - (cdsStart1-1);
          if (offset < 0) return null;
          const codonStart = pos0 - (offset % 3);
          refCodon = refSeq.slice(codonStart, codonStart+3);
          if (refCodon.length !== 3) return null;
          altCodon = refCodon.split('').map((c,i)=> (i === (pos0 - codonStart) ? ev.replacement : c)).join('');
        } else {
          const offset = (cdsStart1-1) - pos0;
          if (offset < 0) return null;
          const r = offset % 3;
          const codonStartPlus = pos0 - (2 - r);
          if (codonStartPlus < 0) return null;
          const plusSlice = refSeq.slice(codonStartPlus, codonStartPlus+3);
          if (plusSlice.length !== 3) return null;
          function rc(s){ const map={A:'T',T:'A',C:'G',G:'C'}; return s.split('').reverse().map(ch=>map[ch]||ch).join(''); }
          refCodon = rc(plusSlice);
          const idxInCodon = 2 - r;
          const repRC = rc(ev.replacement);
          const arr = refCodon.split('');
          arr[idxInCodon] = repRC[0];
          altCodon = arr.join('');
        }
      }
      const aa0 = code[refCodon] || '?';
      const aa1 = code[altCodon] || '?';
      if (aa1 === '*') return 'nonsense';
      if (aa0 === aa1) return 'silent';
      return 'missense';
    }
    const last = frames[frames.length-1];
    if (last && last.new_nodes) {
      for (const [nid, node] of Object.entries(last.new_nodes)) {
        const ev = null; // terminal frame may not have event stored
        // try infer event from earlier edge by searching frames
        let foundEv = null;
        for (let i=frames.length-1; i>=0 && !foundEv; i--) {
          const es = frames[i].new_edges||[];
          for (const e of es) if (e.target === nid && e.event) { foundEv = e.event; break; }
        }
        const seq = node.sequences && node.sequences[chrom];
        if (seq && foundEv) {
          const impact = proteinImpactForEvent(foundEv, genome);
          if (impact) {
            node.metadata = node.metadata || {};
            node.metadata.protein_impact = impact;
          }
        }
      }
    }
  }
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
    frames.forEach((frame, idx) => {
      if (!frame || typeof frame !== "object") return;
      const kind = frame.kind || "helix.edit_dag.frame.v1";
      if (kind !== "helix.edit_dag.frame.v1") return;
      const step = typeof frame.step === "number" ? frame.step : idx;
      self.postMessage({ type: "frame", frame: { ...frame, kind, step } });
    });
    self.postMessage({ type: "done", frameCount: frames.length });
  } catch (error) {
    self.postMessage({ type: "error", error: error.message || String(error) });
  }
});
