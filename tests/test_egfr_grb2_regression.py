import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.append(str(Path(__file__).resolve().parents[1]))

from helix.cli import command_live_run, command_live_inspect


CONFIG = Path("models/egfr_grb2/config.yaml")


def _run_variant_summary(tmp_path, variant: str):
    bundle_dir = tmp_path / f"bundle_{variant}"
    args = SimpleNamespace(
        target=CONFIG,
        hgx=None,
        config=None,
        variant=variant,
        slice="egfr_grb2",
        duration=2.0,
        sync_dt=0.5,
        default_dt=0.1,
        dt_override=[],
        input=[],
        inputs_json=None,
        hz=30.0,
        wall_ms=1.0,
        seed=1234,
        bundle=bundle_dir,
        realtime=False,
    )
    command_live_run(args)

    summary = bundle_dir / "summary.json"
    inspect_args = SimpleNamespace(
        bundle=bundle_dir,
        head=3,
        metrics=2,
        max_islands=5,
        max_nodes=5,
        no_metrics=False,
        metric=["pERK", "agents"],
        json=summary,
    )
    command_live_inspect(inspect_args)
    return json.loads(summary.read_text())


def test_egfr_wt_vs_ko_behavior(tmp_path):
    wt = _run_variant_summary(tmp_path, "wt")
    ko = _run_variant_summary(tmp_path, "ko")
    hig = _run_variant_summary(tmp_path, "hig")

    wt_perk = wt["nodes"]["egfr_rules"]["metrics"]["pERK"]
    ko_perk = ko["nodes"]["egfr_rules"]["metrics"]["pERK"]
    hig_perk = hig["nodes"]["egfr_rules"]["metrics"]["pERK"]

    assert wt_perk["max"] > 0.2
    assert ko_perk["max"] > 0.05
    assert ko_perk["max"] < wt_perk["max"] <= hig_perk["max"]
    assert ko_perk["mean"] < wt_perk["mean"] <= hig_perk["mean"]
    assert hig_perk["max"] >= wt_perk["max"]
    assert hig_perk["mean"] >= wt_perk["mean"]

    wt_cells = wt["nodes"]["epithelium"]["metrics"]["agents"]
    ko_cells = ko["nodes"]["epithelium"]["metrics"]["agents"]
    hig_cells = hig["nodes"]["epithelium"]["metrics"]["agents"]
    assert wt_cells["max"] >= wt_cells["min"]
    assert ko_cells["max"] >= ko_cells["min"]
    assert hig_cells["max"] >= hig_cells["min"]
    assert ko_cells["max"] < wt_cells["max"] <= hig_cells["max"]
    assert ko_cells["mean"] < wt_cells["mean"] <= hig_cells["mean"]
