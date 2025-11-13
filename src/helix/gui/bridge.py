"""Python ↔ JS bridge for the Helix PySide6 GUI."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from PySide6.QtCore import QObject, QThread, Signal, Slot

from helix.crispr.dag_api import build_crispr_edit_dag
from helix.crispr.model import (
    CasSystem,
    CasSystemType,
    DigitalGenome as LegacyGenome,
    GuideRNA,
    PAMRule,
)
from helix.genome.digital import DigitalGenome
from helix.prime.dag_api import build_prime_edit_dag
from helix.prime.model import PegRNA, PrimeEditor
from helix.edit.dag import EditDAGFrame, EditEdge, EditNode


def _parse_inline_fasta(text: str) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    current_name: str | None = None
    buffer: list[str] = []
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("Inline FASTA content was empty.")
    for line in lines:
        if line.startswith(">"):
            if current_name is not None:
                sequences[current_name] = "".join(buffer).upper()
            header = line[1:].strip()
            current_name = header.split()[0] if header else f"sequence_{len(sequences) + 1}"
            buffer = []
        else:
            buffer.append(line)
    if current_name is None:
        # Treat plain sequence without header.
        sequences["sequence_1"] = "".join(lines).upper()
    else:
        sequences[current_name] = "".join(buffer).upper()
    return sequences


def _edit_node_to_payload(node: EditNode) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "log_prob": node.log_prob,
        "metadata": node.metadata,
        "parent_ids": list(node.parents),
        "seq_hashes": node.seq_hashes,
    }
    if node.diffs:
        entry["diffs"] = [
            {
                "chrom": event.chrom,
                "start": event.start,
                "end": event.end,
                "replacement": event.replacement,
                "metadata": event.metadata,
            }
            for event in node.diffs
        ]
    entry["sequences"] = node.genome_view.materialize_all()
    return entry


def _edit_edge_to_payload(edge: EditEdge) -> Dict[str, Any]:
    return {
        "source": edge.source,
        "target": edge.target,
        "rule": edge.rule_name,
        "event": {
            "chrom": edge.event.chrom,
            "start": edge.event.start,
            "end": edge.event.end,
            "replacement": edge.event.replacement,
            "metadata": edge.event.metadata,
        },
        "metadata": edge.metadata,
    }


def _frame_to_payload(frame: EditDAGFrame) -> Dict[str, Any]:
    return {
        "kind": "helix.edit_dag.frame.v1",
        "step": frame.step,
        "new_nodes": {node_id: _edit_node_to_payload(node) for node_id, node in frame.new_nodes.items()},
        "new_edges": [_edit_edge_to_payload(edge) for edge in frame.new_edges],
    }


class SimulationWorker(QThread):
    """Background thread that runs CRISPR/Prime simulations and emits frames."""

    frameProduced = Signal(str)
    finishedSimulation = Signal()
    errorRaised = Signal(str)

    def __init__(self, spec: Dict[str, Any], parent: QObject | None = None):
        super().__init__(parent)
        self.spec = spec

    def _load_core_genome(self) -> DigitalGenome:
        genome_spec = self.spec.get("genome", {})
        fasta_field = genome_spec.get("fasta", "")
        inline_content = genome_spec.get("inline_fasta")
        sequences: Dict[str, str]
        if inline_content:
            sequences = _parse_inline_fasta(inline_content)
            core = DigitalGenome(sequences=sequences)
        elif fasta_field and Path(fasta_field).expanduser().exists():
            core = DigitalGenome.from_fasta(Path(fasta_field).expanduser())
        elif fasta_field.strip().startswith(">"):
            sequences = _parse_inline_fasta(fasta_field)
            core = DigitalGenome(sequences=sequences)
        else:
            raise FileNotFoundError(
                "Genome FASTA was not found. Provide an absolute path or inline FASTA content."
            )
        region = genome_spec.get("region")
        if region:
            core = core.slice_region(region)
        return core

    def _legacy_genome(self) -> LegacyGenome:
        core = self._load_core_genome()
        return LegacyGenome(sequences=dict(core.sequences))

    def _load_cas(self, cfg: Dict[str, Any]) -> CasSystem:
        if isinstance(cfg, (str, Path)):
            data = json.loads(Path(cfg).read_text(encoding="utf-8"))
        else:
            data = dict(cfg)
        return CasSystem(
            name=data.get("name", "SpCas9"),
            system_type=CasSystemType(data.get("system_type", "cas9")),
            pam_rules=[PAMRule(pattern=data.get("pam_pattern", "NGG"))],
            cut_offset=int(data.get("cut_offset", 3)),
            max_mismatches=int(data.get("max_mismatches", 3)),
            weight_mismatch_penalty=float(data.get("weight_mismatch_penalty", 1.0)),
            weight_pam_penalty=float(data.get("weight_pam_penalty", 2.0)),
        )

    def _load_prime_editor(self, cfg: Dict[str, Any]) -> PrimeEditor:
        if isinstance(cfg, (str, Path)):
            data = json.loads(Path(cfg).read_text(encoding="utf-8"))
        else:
            data = dict(cfg)
        cas_cfg = data.get("cas") or {}
        cas = self._load_cas(cas_cfg)
        return PrimeEditor(
            name=data.get("name", "PrimeEditor"),
            cas=cas,
            nick_to_edit_offset=int(data.get("nick_to_edit_offset", data.get("nick_to_rtt_offset", 0))),
            efficiency_scale=float(data.get("efficiency_scale", 0.5)),
            indel_bias=float(data.get("indel_bias", 0.1)),
            mismatch_tolerance=int(data.get("mismatch_tolerance", 2)),
        )

    def _emit_frame(self, frame: EditDAGFrame, mechanism: str, meta: Dict[str, Any]) -> None:
        payload = _frame_to_payload(frame)
        payload["meta"] = {"mechanism": mechanism, **meta}
        self.frameProduced.emit(json.dumps(payload))

    def _parse_coding(self):
        coding = self.spec.get("coding") or {}
        try:
            cds_start = int(coding.get("cds_start")) if coding.get("cds_start") else None
            cds_end = int(coding.get("cds_end")) if coding.get("cds_end") else None
        except Exception:
            cds_start = cds_end = None
        strand = (coding.get("strand") or "+").strip()
        exons_text = (coding.get("exons") or "").strip()
        exons = []
        if exons_text:
            for token in [t for t in exons_text.replace(";", ",").split(",") if t.strip()]:
                m = None
                try:
                    a, b = token.split("-")
                    s = int(a)
                    e = int(b)
                    exons.append((min(s, e), max(s, e)))
                except Exception:
                    continue
        return {
            "cds_start": cds_start,
            "cds_end": cds_end,
            "strand": strand,
            "exons": exons,
        }

    def _build_coding_map(self, chrom: str, ref_seq: str, cds_start: int | None, cds_end: int | None, strand: str, exons: list[tuple[int, int]] | None):
        if not cds_start or not cds_end or cds_end < cds_start:
            return None
        blocks: list[tuple[int, int]] = []
        if exons:
            for s1, e1 in sorted(exons, key=lambda x: x[0]):
                s = max(s1, cds_start)
                e = min(e1, cds_end)
                if s <= e:
                    blocks.append((s, e))
        else:
            blocks.append((cds_start, cds_end))
        cmap: list[int] = []
        if strand == "+":
            for s, e in blocks:
                for p in range(s - 1, e):
                    if 0 <= p < len(ref_seq):
                        cmap.append(p)
        else:
            for s, e in reversed(blocks):
                for p in range(e - 1, s - 1, -1):
                    if 0 <= p < len(ref_seq):
                        cmap.append(p)
        return cmap if cmap else None

    @staticmethod
    def _rc_base(b: str) -> str:
        return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, b)

    def _protein_impact_snv(self, event, ref_seq: str, strand: str, coding_map: list[int] | None) -> str | None:
        # event has .start, .end, .replacement
        if event is None:
            return None
        length = max(0, int(event.end) - int(event.start))
        ins = len(event.replacement or "")
        if not (length == 1 and ins == 1):
            return None
        if not coding_map:
            return None
        pos0 = int(event.start)
        try:
            t_idx = coding_map.index(pos0)
        except ValueError:
            return None
        frame = t_idx % 3
        if t_idx - frame < 0 or t_idx - frame + 2 >= len(coding_map):
            return None
        codon_pos = [coding_map[t_idx - frame + i] for i in range(3)]
        bases = [ref_seq[p] if 0 <= p < len(ref_seq) else "N" for p in codon_pos]
        if strand == "-":
            # convert to transcript orientation: rc and keep order
            ref_codon = "".join(self._rc_base(b) for b in bases)
            idx = frame
            rep = self._rc_base((event.replacement or "N")[0])
            alt = list(ref_codon)
            alt[idx] = rep
            alt_codon = "".join(alt)
        else:
            ref_codon = "".join(bases)
            idx = frame
            alt = list(ref_codon)
            alt[idx] = (event.replacement or "N")[0]
            alt_codon = "".join(alt)
        CODE = {
            "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
            "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
            "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
            "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
            "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
            "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
            "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
            "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
        }
        aa0 = CODE.get(ref_codon)
        aa1 = CODE.get(alt_codon)
        if not aa0 or not aa1:
            return None
        if aa1 == "*":
            return "nonsense"
        if aa0 == aa1:
            return "silent"
        return "missense"

    def run(self) -> None:  # pragma: no cover - GUI thread
        try:
            mode = self.spec.get("mode", "crispr").lower()
            sim_cfg = self.spec.get("simulation", {})
            genome = self._legacy_genome()
            coding_spec = self._parse_coding()
            # build coding maps per chrom
            coding_maps: dict[str, list[int] | None] = {}
            if coding_spec.get("cds_start") and coding_spec.get("cds_end"):
                for chrom, seq in genome.sequences.items():
                    cmap = self._build_coding_map(
                        chrom,
                        seq,
                        coding_spec["cds_start"],
                        coding_spec["cds_end"],
                        coding_spec["strand"],
                        coding_spec.get("exons"),
                    )
                    coding_maps[chrom] = cmap
            if mode == "prime":
                editor_cfg = self.spec.get("editor") or {}
                peg_cfg = self.spec.get("peg") or {}
                for field in ("spacer", "pbs", "rtt"):
                    if field not in peg_cfg:
                        raise ValueError(f"PegRNA configuration missing '{field}'.")
                editor = self._load_prime_editor(editor_cfg)
                peg = PegRNA(
                    spacer=peg_cfg["spacer"].upper(),
                    pbs=peg_cfg["pbs"].upper(),
                    rtt=peg_cfg["rtt"].upper(),
                    name=peg_cfg.get("name") or "PEG",
                    metadata=peg_cfg.get("metadata") or {},
                )

                def consumer(frame: EditDAGFrame, mech: str = "prime", peg_name: str | None = peg.name):
                    # annotate protein impact (SNV) if coding map present
                    for node in frame.new_nodes.values():
                        try:
                            ev = node.diffs[0] if node.diffs else None
                            if ev and coding_maps:
                                cmap = coding_maps.get(ev.chrom)
                                impact = self._protein_impact_snv(ev, genome.sequences.get(ev.chrom, ""), coding_spec.get("strand", "+"), cmap)
                                if impact:
                                    node.metadata["protein_impact"] = impact
                        except Exception:
                            pass
                    self._emit_frame(frame, mech, {"peg_name": peg_name or "PEG"})

                build_prime_edit_dag(
                    genome,
                    editor,
                    peg,
                    rng_seed=int(sim_cfg.get("seed", 0)),
                    max_depth=int(sim_cfg.get("max_depth", 3)),
                    min_prob=float(sim_cfg.get("min_prob", 1e-4)),
                    frame_consumer=consumer,
                )
            else:
                cas_cfg = self.spec.get("cas") or {}
                guide_cfg = self.spec.get("guide") or {}
                if "sequence" not in guide_cfg:
                    raise ValueError("Guide configuration requires a 'sequence' field.")
                cas = self._load_cas(cas_cfg)
                guide = GuideRNA(sequence=guide_cfg["sequence"].upper(), name=guide_cfg.get("name"))
                max_sites = sim_cfg.get("max_sites")
                if max_sites is not None:
                    max_sites = int(max_sites)

                def consumer(frame: EditDAGFrame, mech: str = "crispr", guide_name: str | None = guide.name):
                    for node in frame.new_nodes.values():
                        try:
                            ev = node.diffs[0] if node.diffs else None
                            if ev and coding_maps:
                                cmap = coding_maps.get(ev.chrom)
                                impact = self._protein_impact_snv(ev, genome.sequences.get(ev.chrom, ""), coding_spec.get("strand", "+"), cmap)
                                if impact:
                                    node.metadata["protein_impact"] = impact
                        except Exception:
                            pass
                    self._emit_frame(frame, mech, {"guide_name": guide_name or "GUIDE"})

                build_crispr_edit_dag(
                    genome,
                    cas,
                    guide,
                    rng_seed=int(sim_cfg.get("seed", 0)),
                    max_depth=int(sim_cfg.get("max_depth", 2)),
                    min_prob=float(sim_cfg.get("min_prob", 1e-4)),
                    max_sites=max_sites,
                    frame_consumer=consumer,
                    use_gpu=bool(sim_cfg.get("use_gpu", False)),
                )
        except Exception as exc:  # pragma: no cover - GUI thread
            self.errorRaised.emit(str(exc))
        finally:  # pragma: no cover - GUI thread
            self.finishedSimulation.emit()


class SimulationBridge(QObject):
    """Expose slots/signals to the QWebChannel JS client."""

    frameReceived = Signal(str)
    simFinished = Signal()
    simErrored = Signal(str)

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self._worker: SimulationWorker | None = None

    @Slot(str)
    def runSimulation(self, spec_json: str) -> None:  # noqa: N802
        spec = json.loads(spec_json)
        if self._worker is not None and self._worker.isRunning():
            self.simErrored.emit("Simulation already running. Please wait for it to finish.")
            return
        self._worker = SimulationWorker(spec, parent=self)
        self._worker.frameProduced.connect(self.frameReceived)
        self._worker.finishedSimulation.connect(self.simFinished)
        self._worker.errorRaised.connect(self.simErrored)
        self._worker.finishedSimulation.connect(self._clear_worker)
        self._worker.start()

    @Slot()
    def _clear_worker(self) -> None:
        self._worker = None
