import json
from pathlib import Path

import numpy as np

from helix.crispr.model import DigitalGenome, CasSystem, CasSystemType, PAMRule
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.simulator import PrimeTargetRequest, predict_prime_outcomes_for_targets


def _cas() -> CasSystem:
    return CasSystem(name="SpCas9", system_type=CasSystemType.CAS9, pam_rules=[PAMRule(pattern="NGG")], cut_offset=3)


def _editor() -> PrimeEditor:
    return PrimeEditor(name="pe-mini", cas=_cas(), nick_to_edit_offset=1, efficiency_scale=0.8)


def test_predict_prime_outcomes_batch():
    genome = DigitalGenome({"chr": "ACGT" * 8})
    peg = PegRNA(spacer="ACGTACGTACGTACG", pbs="GAAAC", rtt="TTTTAA")
    req = PrimeTargetRequest(target_id="t1", genome=genome, peg=peg, editor=_editor())
    preds = predict_prime_outcomes_for_targets([req])
    assert len(preds) == 1
    assert 0.0 <= preds[0].predicted_efficiency <= 1.0


def test_prime_engine_golden_fixtures():
    fixture_path = Path("tests/data/engine_prime/fixtures.json")
    scores_path = Path("tests/data/engine_prime/scores_reference.json")
    fixtures = json.loads(fixture_path.read_text())
    expected = np.array(json.loads(scores_path.read_text()), dtype=float)
    cas = _cas()
    requests = []
    for entry in fixtures["targets"]:
        genome = DigitalGenome(entry["genome"])
        peg = PegRNA(**entry["peg"])
        editor = PrimeEditor(cas=cas, **entry["editor"])
        requests.append(PrimeTargetRequest(target_id=entry["id"], genome=genome, peg=peg, editor=editor))
    preds = predict_prime_outcomes_for_targets(requests)
    values = np.array([pred.predicted_efficiency for pred in preds], dtype=float)
    assert np.allclose(values, expected)
