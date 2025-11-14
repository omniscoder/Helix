from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))

from helixtasks.egfr_grb2_orchestrator import build_egfr_grb2


def test_egfr_grb2_orchestrator_runs_short_simulation(tmp_path):
    config_path = Path("models/egfr_grb2/config.yaml")
    graph, scheduler = build_egfr_grb2(config_path, variant="wt")
    scheduler.run_until(5.0)
    assert scheduler.snapshots, "Expected snapshots to be recorded"
    last = scheduler.snapshots[-1]
    assert "egfr_rules" in last
    assert last["egfr_rules"]["pERK"] >= 0.0
