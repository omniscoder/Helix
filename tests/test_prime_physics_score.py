"""Prime physics scoring heuristics tests."""

from __future__ import annotations

from helix.crispr.model import CasSystem, CasSystemType, PAMRule
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.physics import score_prime_design


def _editor() -> PrimeEditor:
    cas = CasSystem(name="demo-cas9", system_type=CasSystemType.CAS9, pam_rules=[PAMRule(pattern="NGG")], cut_offset=3)
    return PrimeEditor(name="PE-test", cas=cas, nick_to_edit_offset=1, efficiency_scale=0.8, indel_bias=0.1)


def test_pbs_length_window_affects_prediction():
    target = ("C" * 30) + ("ACGT" * 10)
    editor = _editor()
    spacer = target[30:50]
    peg_short = PegRNA(spacer=spacer, pbs="ATATA", rtt=target[28:40])
    peg_mid = PegRNA(spacer=spacer, pbs=target[16:30], rtt=target[28:40])
    peg_long = PegRNA(spacer=spacer, pbs=target[6:30], rtt=target[28:40])
    score_short = score_prime_design(peg_short, target, editor)
    score_mid = score_prime_design(peg_mid, target, editor)
    score_long = score_prime_design(peg_long, target, editor)
    assert score_mid.E_pred > score_short.E_pred


def test_microhomology_influences_flap_probability():
    target = "ACGT" * 12
    editor = _editor()
    spacer = target[0:20]
    peg_good = PegRNA(spacer=spacer, pbs=target[0:12], rtt=target[20:32])
    peg_bad = PegRNA(spacer=spacer, pbs=target[0:12], rtt="T" * 12)
    good_score = score_prime_design(peg_good, target, editor)
    bad_score = score_prime_design(peg_bad, target, editor)
    assert good_score.P_flap > bad_score.P_flap
