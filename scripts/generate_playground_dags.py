#!/usr/bin/env python3
"""
Generate canonical CRISPR/Prime edit DAG samples for the playground.

This script runs the in-repo simulators so the docs assets always reflect
the latest physics. Outputs land in docs/assets/viz/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from helix.cli import _edit_dag_to_payload  # type: ignore
from helix.crispr.dag_api import build_crispr_edit_dag
from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, GuideRNA, PAMRule
from helix.prime.dag_api import build_prime_edit_dag
from helix.prime.model import PegRNA, PrimeEditor


def _demo_sequence() -> str:
    return "TTTACCCAGGAAACCCGGGTTTTAGGTTT"


def generate_crispr_payload() -> dict:
    genome = DigitalGenome({"chr_demo": _demo_sequence()})
    cas = CasSystem(
        name="SpCas9-demo",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG", description="SpCas9 canonical PAM")],
        cut_offset=3,
        max_mismatches=3,
    )
    guide = GuideRNA(sequence="ACCCAGGAAACCCGGGTTTT")
    dag = build_crispr_edit_dag(genome, cas, guide, rng_seed=7, max_depth=2, max_sites=3)
    metadata = {"mechanism": "crispr", "demo": "playground"}
    return _edit_dag_to_payload(dag, artifact="helix.crispr.edit_dag.v1.1", metadata=metadata)


def generate_prime_payload() -> dict:
    genome = DigitalGenome({"chr_demo": _demo_sequence()})
    cas = CasSystem(
        name="SpCas9-H840A",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
    )
    peg = PegRNA(spacer="ACCCAGGAAACCCGGGTTTT", pbs="GAAAC", rtt="TTTTAA")
    editor = PrimeEditor(
        name="PE2-demo",
        cas=cas,
        nick_to_edit_offset=1,
        efficiency_scale=0.75,
        indel_bias=0.15,
        flap_balance=0.6,
    )
    dag = build_prime_edit_dag(genome, editor, peg, rng_seed=11, max_depth=2)
    metadata = {"mechanism": "prime", "demo": "playground"}
    return _edit_dag_to_payload(dag, artifact="helix.prime.edit_dag.v1.1", metadata=metadata)


def main() -> None:
    viz_dir = ROOT / "docs" / "assets" / "viz"
    data_dir = ROOT / "docs" / "data"
    viz_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    crispr_payload = generate_crispr_payload()
    prime_payload = generate_prime_payload()
    assets = {
        "crispr_playground_dag.json": crispr_payload,
        "prime_playground_dag.json": prime_payload,
    }
    demos = {
        "crispr_demo.edit_dag.json": crispr_payload,
        "prime_demo.edit_dag.json": prime_payload,
    }
    for filename, payload in assets.items():
        path = viz_dir / filename
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"Wrote {path.relative_to(ROOT)}")
    for filename, payload in demos.items():
        path = data_dir / filename
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"Wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
