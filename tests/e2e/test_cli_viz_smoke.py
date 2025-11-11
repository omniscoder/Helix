import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


pytestmark = pytest.mark.skipif(
    shutil.which(sys.executable) is None,
    reason="Python executable not found for CLI tests",
)


PAYLOADS = {
    "minimizers": {
        "data": {"sequence_length": 1000, "minimizers": [5, [10, "AAA", 1], {"pos": 20}]},
        "extra": ["--bins", "10"],
    },
    "seed-chain": {
        "data": {
            "ref_length": 500,
            "qry_length": 480,
            "chains": [
                [{"ref_start": 10, "ref_end": 40, "qry_start": 5, "qry_end": 35}],
                [{"ref_start": 60, "ref_end": 90, "qry_start": 55, "qry_end": 85}],
            ],
        },
        "extra": [],
    },
    "rna-dotplot": {
        "data": {"posterior": [[0.0, 0.6], [0.6, 0.0]]},
        "extra": ["--vmin", "0", "--vmax", "1"],
    },
    "alignment-ribbon": {
        "data": {
            "ref_length": 300,
            "qry_length": 290,
            "ref_start": 20,
            "qry_start": 15,
            "cigar": "20M2I10M3D15M",
            "metadata": {"name": "read_001"},
        },
        "extra": [],
    },
    "distance-heatmap": {
        "data": {
            "labels": ["A", "B", "C"],
            "matrix": [[0.0, 0.05, 0.1], [0.05, 0.0, 0.08], [0.1, 0.08, 0.0]],
        },
        "extra": [],
    },
    "motif-logo": {
        "data": {
            "alphabet": ["A", "C", "G", "T"],
            "pwm": [
                [0.25, 0.25, 0.25, 0.25],
                [0.05, 0.05, 0.85, 0.05],
                [0.6, 0.1, 0.1, 0.2],
            ],
            "background": [0.25, 0.25, 0.25, 0.25],
            "consensus": "AGT",
        },
        "extra": [],
    },
}


def test_cli_viz_smoke(tmp_path):
    env = os.environ.copy()
    project_root = Path(__file__).resolve().parents[2]
    src_path = project_root / "src"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(src_path)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
    )

    for command, payload in PAYLOADS.items():
        input_path = tmp_path / f"{command}.json"
        input_path.write_text(json.dumps(payload["data"]), encoding="utf-8")
        output_img = tmp_path / f"{command}.png"
        args = [
            sys.executable,
            "-m",
            "helix.cli",
            "viz",
            command,
            "--input",
            str(input_path),
            "--save",
            str(output_img),
        ] + payload["extra"]
        proc = subprocess.run(args, capture_output=True, text=True, env=env, check=False)
        assert proc.returncode == 0, proc.stderr
        assert output_img.exists(), f"{command} should produce an image"
        viz_json = tmp_path / f"{command}.viz.json"
        assert viz_json.exists(), f"{command} should emit a viz-spec JSON"
        spec = json.loads(viz_json.read_text(encoding="utf-8"))
        assert spec["kind"], "viz-spec must include a kind"
        assert "input_sha256" in spec.get("meta", {})
        prov = tmp_path / f"{command}.provenance.json"
        assert prov.exists(), f"{command} should emit a provenance JSON"
        provenance = json.loads(prov.read_text(encoding="utf-8"))
        assert provenance["schema_kind"]
        assert provenance["image_sha256"]
