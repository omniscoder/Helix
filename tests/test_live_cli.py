import json
from pathlib import Path
from types import SimpleNamespace

from helix.cli import command_live_run, command_live_inspect


def test_live_cli_run_bundle(tmp_path):
    hgx = Path("examples/egfr_monolayer.hgx")
    bundle_dir = tmp_path / "bundle"
    args = SimpleNamespace(
        target=None,
        hgx=hgx,
        config=None,
        variant=None,
        slice=None,
        duration=0.5,
        sync_dt=0.25,
        default_dt=0.1,
        dt_override=[],
        input=[],
        inputs_json=None,
        hz=30.0,
        wall_ms=1.0,
        seed=123,
        bundle=bundle_dir,
        realtime=False,
    )
    command_live_run(args)

    model_copy = bundle_dir / "model.hgx"
    snapshots_path = bundle_dir / "snapshots.json"
    meta_path = bundle_dir / "meta.json"

    assert model_copy.exists()
    snapshots = json.loads(snapshots_path.read_text())
    assert len(snapshots) > 0
    meta = json.loads(meta_path.read_text())
    assert meta["model"] == "EGFR Monolayer"
    assert meta["seed"] == 123


def test_live_cli_time_series_inputs(tmp_path):
    hgx = Path("examples/egfr_monolayer.hgx")
    bundle_dir = tmp_path / "bundle_ts"
    inputs_path = tmp_path / "inputs.json"
    inputs_path.write_text(
        json.dumps(
            {
                "egfr_rules.drive": [
                    {"t": 0.0, "value": 0.0},
                    {"t": 0.2, "value": 1.0},
                ]
            }
        ),
        encoding="utf-8",
    )
    args = SimpleNamespace(
        target=None,
        hgx=hgx,
        config=None,
        variant=None,
        slice=None,
        duration=0.6,
        sync_dt=0.2,
        default_dt=0.1,
        dt_override=[],
        input=[],
        inputs_json=inputs_path,
        hz=30.0,
        wall_ms=1.0,
        seed=0,
        bundle=bundle_dir,
        realtime=False,
    )
    command_live_run(args)

    snapshots = json.loads((bundle_dir / "snapshots.json").read_text())
    final_rate = snapshots[-1]["egfr_rules"]["rate"]
    assert final_rate > 0.15


def test_live_cli_inspect(tmp_path, capsys):
    hgx = Path("examples/egfr_monolayer.hgx")
    bundle_dir = tmp_path / "bundle_inspect"
    args = SimpleNamespace(
        target=None,
        hgx=hgx,
        config=None,
        variant=None,
        slice=None,
        duration=0.4,
        sync_dt=0.2,
        default_dt=0.1,
        dt_override=[],
        input=[],
        inputs_json=None,
        hz=30.0,
        wall_ms=1.0,
        seed=42,
        bundle=bundle_dir,
        realtime=False,
    )
    command_live_run(args)

    inspect_args = SimpleNamespace(
        bundle=bundle_dir,
        head=2,
        metrics=2,
        max_islands=2,
        max_nodes=2,
        metric=[],
        no_metrics=False,
        json=None,
    )
    command_live_inspect(inspect_args)
    captured = capsys.readouterr()
    assert "Model: EGFR Monolayer" in captured.out
    assert "snapshots" in captured.out

    json_summary = tmp_path / "inspect.json"
    inspect_args.json = json_summary
    command_live_inspect(inspect_args)
    payload = json.loads(json_summary.read_text())
    assert payload["schema_version"] == 1
    runtime = payload["runtime"]
    assert isinstance(runtime["islands"], int)
    assert runtime["snapshots"] >= 1
    nodes = payload["nodes"]
    assert "wound_field" in nodes
    tile_metrics = nodes["wound_field"]["metrics"].get("tile")
    assert tile_metrics["min"] <= tile_metrics["mean"] <= tile_metrics["max"]
    assert "var" in tile_metrics
    assert "std" in tile_metrics

    inspect_args.no_metrics = True
    inspect_args.json = tmp_path / "inspect_no_metrics.json"
    command_live_inspect(inspect_args)
    payload_no_metrics = json.loads(inspect_args.json.read_text())
    assert payload_no_metrics["nodes"]["wound_field"]["metrics"] == {}


def test_live_run_with_config(tmp_path):
    args = SimpleNamespace(
        target=None,
        hgx=None,
        config=None,
        slice="egfr_grb2",
        variant="wt",
        duration=0.5,
        sync_dt=0.25,
        default_dt=0.1,
        dt_override=[],
        input=[],
        inputs_json=None,
        hz=30.0,
        wall_ms=1.0,
        seed=1,
        bundle=tmp_path / "egfr_bundle",
        realtime=False,
    )
    command_live_run(args)
    meta = json.loads((tmp_path / "egfr_bundle" / "meta.json").read_text())
    assert meta["variant"] == "wt"
    assert meta["target_type"] == "config"
