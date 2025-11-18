from __future__ import annotations

import json
import shutil
from argparse import Namespace
from pathlib import Path

from helix.cli import command_headless_run, command_headless_compare, command_headless_report
from helix.snapshot_bundle import pack_snapshot_bundle

CRISPR_GENOME = "ACGT" * 16
CRISPR_GUIDE = {
    "id": "guide-1",
    "sequence": "GGGGTTTAGAGCTATGCT",
    "start": 4,
    "end": 23,
    "strand": "+",
    "gc_content": 0.5,
}
CRISPR_RUN_CONFIG = {
    "draws": 32,
    "seed": 17,
    "pam_profile": "SpCas9_NGG",
    "pam_softness": 0.2,
    "priors_profile": "default_indel",
}
PCR_TEMPLATE = "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT"
PCR_CONFIG = {
    "fwd_primer": "AAACCCGG",
    "rev_primer": "CCCGGGTT",
    "cycles": 10,
    "polymerase_name": "TestPol",
    "base_error_rate": 1e-5,
    "indel_fraction": 0.1,
    "efficiency": 0.8,
    "initial_copies": 1,
}
PRIME_GENOME = "ACGTACGTACGTNNNN" * 2
PRIME_PEG = {
    "name": "peg-1",
    "spacer": "GGGGTTTAGAGCTATGCTA",
    "pbs": "ACGTAC",
    "rtt": "TTTTAAAACC",
}
PRIME_RUN_CONFIG = {
    "draws": 16,
    "seed": 33,
    "pam_profile": "SpCas9_NGG",
    "pam_softness": 0.2,
    "priors_profile": "default_indel",
}
PRIME_EDITOR = {
    "name": "PrimeEditor",
    "nick_to_edit_offset": 3,
    "efficiency_scale": 1.0,
    "indel_bias": 0.1,
    "mismatch_tolerance": 3,
    "flap_balance": 0.5,
    "reanneal_bias": 0.1,
    "metadata": {},
    "cas": {"name": "SimPrime", "pam_rules": [{"pattern": "NGG"}], "system_type": "CAS9", "cut_offset": 3},
}


def _write_json(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _make_bundle(tmp_path: Path) -> Path:
    session = {"version": 2, "state": {"run_kind": "CRISPR"}}
    run_payload = {
        "state": {
            "run_kind": "CRISPR",
            "run_id": 1,
            "viz_dirty": False,
            "config": {
                "guide": CRISPR_GUIDE,
                "run_config": dict(CRISPR_RUN_CONFIG),
                "sim_type": "CRISPR",
            },
            "genome_source": "inline",
            "genome": {"chrGui": CRISPR_GENOME},
            "outcomes": [],
        }
    }
    session_path = _write_json(tmp_path / "session.json", session)
    run_path = _write_json(tmp_path / "run.json", run_payload)
    bundle_path = tmp_path / "bundle.hxs"
    pack_snapshot_bundle(session_path=session_path, run_paths=[run_path], out_path=bundle_path)
    return bundle_path


def _make_prime_bundle(tmp_path: Path) -> Path:
    session = {"version": 2, "state": {"run_kind": "PRIME"}}
    run_payload = {
        "state": {
            "run_kind": "PRIME",
            "run_id": "prime-1",
            "viz_dirty": False,
            "peg": PRIME_PEG,
            "editor": PRIME_EDITOR,
            "config": {
                "run_config": dict(PRIME_RUN_CONFIG),
                "sim_type": "PRIME",
            },
            "genome_source": "inline",
            "genome": {"chrPrime": PRIME_GENOME},
            "outcomes": [],
        }
    }
    session_path = _write_json(tmp_path / "prime_session.json", session)
    run_path = _write_json(tmp_path / "prime_run.json", run_payload)
    bundle_path = tmp_path / "prime_bundle.hxs"
    pack_snapshot_bundle(session_path=session_path, run_paths=[run_path], out_path=bundle_path)
    return bundle_path


def _make_pcr_bundle(tmp_path: Path) -> Path:
    session = {"version": 2, "state": {"run_kind": "PCR"}}
    run_payload = {
        "state": {
            "run_kind": "PCR",
            "run_id": "pcr-1",
            "viz_dirty": False,
            "pcr_config": {
                "template_seq": PCR_TEMPLATE,
                **PCR_CONFIG,
            },
            "config": {
                "run_config": dict(PCR_CONFIG),
                "sim_type": "PCR",
            },
            "genome_source": "inline",
            "outcomes": [],
        }
    }
    session_path = _write_json(tmp_path / "pcr_session.json", session)
    run_path = _write_json(tmp_path / "pcr_run.json", run_payload)
    bundle_path = tmp_path / "pcr_bundle.hxs"
    pack_snapshot_bundle(session_path=session_path, run_paths=[run_path], out_path=bundle_path)
    return bundle_path


def test_headless_run_local_emits_artifacts_and_telemetry(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    out_dir = tmp_path / "out_local"
    args = Namespace(
        snapshot=bundle,
        run_id="1",
        kind=None,
        out=out_dir,
        params=None,
        engine="local",
        queue="gpu-a100",
    )
    command_headless_run(args)
    assert (out_dir / "snapshot.json").exists()
    telemetry_path = out_dir / "telemetry.jsonl"
    record = json.loads(telemetry_path.read_text().strip())
    assert record["event"] == "run_completed"
    assert record["engine"] == "local"


def test_local_run_produces_metrics(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    out_dir = tmp_path / "metrics_local"
    args = Namespace(
        snapshot=bundle,
        run_id="1",
        kind=None,
        out=out_dir,
        params=None,
        engine="local",
        queue="gpu-a100",
    )
    command_headless_run(args)
    metrics = json.loads((out_dir / "metrics.json").read_text())
    assert metrics["counts"]["intended"] + metrics["counts"]["off_target"] >= 0


def test_prime_run_executes_locally(tmp_path: Path) -> None:
    bundle = _make_prime_bundle(tmp_path)
    out_dir = tmp_path / "prime_local"
    args = Namespace(
        snapshot=bundle,
        run_id="prime-1",
        kind="Prime",
        out=out_dir,
        params=None,
        engine="local",
        queue="gpu-a100",
    )
    command_headless_run(args)
    metrics = json.loads((out_dir / "metrics.json").read_text())
    assert metrics["counts"]["neutral"] >= 0


def test_pcr_run_executes_locally(tmp_path: Path) -> None:
    bundle = _make_pcr_bundle(tmp_path)
    out_dir = tmp_path / "pcr_local"
    args = Namespace(
        snapshot=bundle,
        run_id="pcr-1",
        kind="PCR",
        out=out_dir,
        params=None,
        engine="local",
        queue="gpu-a100",
    )
    command_headless_run(args)
    metrics = json.loads((out_dir / "metrics.json").read_text())
    pcr_metrics = metrics.get("pcr")
    assert pcr_metrics is not None
    assert pcr_metrics.get("final_mass_ng", 0) > 0


def test_headless_run_ogn_submits_and_logs(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    out_dir = tmp_path / "out_ogn"
    args = Namespace(
        snapshot=bundle,
        run_id="1",
        kind=None,
        out=out_dir,
        params=None,
        engine="ogn",
        queue="gpu-h100",
    )
    command_headless_run(args)
    telemetry_path = out_dir / "telemetry.jsonl"
    record = json.loads(telemetry_path.read_text().strip())
    assert record["event"] == "run_completed"
    assert record["engine"] == "ogn"
    assert record["ogn_queue"] == "gpu-h100"
    shutil.rmtree(Path("mock_ogn_runs"), ignore_errors=True)


def test_headless_compare_and_report_emit_outputs(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    local_dir = tmp_path / "local"
    remote_dir = tmp_path / "remote"
    command_headless_run(
        Namespace(
            snapshot=bundle,
            run_id="1",
            kind=None,
            out=local_dir,
            params=None,
            engine="local",
            queue="gpu-a100",
        )
    )
    command_headless_run(
        Namespace(
            snapshot=bundle,
            run_id="1",
            kind=None,
            out=remote_dir,
            params=None,
            engine="local",
            queue="gpu-a100",
        )
    )
    local_snapshot = str(local_dir / "snapshot.json")
    remote_snapshot = str(remote_dir / "snapshot.json")
    compare_out = tmp_path / "compare" / "compare.json"
    compare_out.parent.mkdir(parents=True, exist_ok=True)
    compare_args = Namespace(runs=[local_snapshot, remote_snapshot], snapshot=None, out=compare_out)
    command_headless_compare(compare_args)
    assert compare_out.exists()
    telemetry_path = compare_out.parent / "telemetry.jsonl"
    assert telemetry_path.exists()
    compare_record = json.loads(telemetry_path.read_text().strip())
    assert compare_record["event"] == "compare_completed"

    report_out = tmp_path / "reports" / "report.md"
    report_out.parent.mkdir(parents=True, exist_ok=True)
    report_args = Namespace(run=Path(local_snapshot), format="md", out=report_out)
    command_headless_report(report_args)
    assert report_out.exists()
    report_log = (report_out.parent / "telemetry.jsonl").read_text().strip()
    assert "report_completed" in report_log
