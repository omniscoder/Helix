"""Microbenchmark for prime editing batch predictions."""
from __future__ import annotations

import argparse
import random
import time

from helix import bioinformatics
from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, PAMRule
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.simulator import PrimeTargetRequest, predict_prime_outcomes_for_targets


def _random_dna(length: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(length))


def _make_cas() -> CasSystem:
    return CasSystem(name="SpCas9", system_type=CasSystemType.CAS9, pam_rules=[PAMRule(pattern="NGG")], cut_offset=3)


def _make_editor() -> PrimeEditor:
    return PrimeEditor(name="PE-demo", cas=_make_cas(), nick_to_edit_offset=1, efficiency_scale=0.8, indel_bias=0.1)


def run_benchmark(targets: int, genome_len: int, spacer_len: int, backend: str) -> None:
    genome = DigitalGenome({"chr_demo": _random_dna(genome_len)})
    editor = _make_editor()
    requests = []
    for idx in range(targets):
        start = random.randint(0, max(0, genome_len - spacer_len))
        spacer = genome.sequences["chr_demo"][start : start + spacer_len]
        peg = PegRNA(spacer=spacer, pbs="GAAAC", rtt="TTTTAA")
        requests.append(PrimeTargetRequest(target_id=f"target_{idx}", genome=genome, peg=peg, editor=editor))
    start = time.perf_counter()
    predict_prime_outcomes_for_targets(requests)
    elapsed = time.perf_counter() - start
    print(f"prime backend={backend} targets={targets} genome_len={genome_len} -> {targets/elapsed:.2f} preds/sec")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prime engine microbenchmark")
    parser.add_argument("--targets", type=int, default=32)
    parser.add_argument("--genome-len", type=int, default=2000)
    parser.add_argument("--spacer-len", type=int, default=20)
    parser.add_argument("--backend", default="cpu-reference")
    args = parser.parse_args()
    run_benchmark(args.targets, args.genome_len, args.spacer_len, args.backend)


if __name__ == "__main__":
    main()
