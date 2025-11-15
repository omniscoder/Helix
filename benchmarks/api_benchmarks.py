"""Targeted benchmarks for the helix.api surface area."""
from __future__ import annotations

import argparse
import json
import locale
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Callable, Dict, List, Mapping, Sequence

import numpy as np
import yaml

try:
    import psutil
except ImportError:  # pragma: no cover - optional dependency
    psutil = None


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = REPO_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from helix import api as hx  # noqa: E402
from helix import datasets  # noqa: E402

SCHEMA_KIND = "bench_result"
SCHEMA_VERSION = "1.0"
DEFAULT_SEED = 1337


# Base datasets ---------------------------------------------------------------

def _load_sequence(default_relative: str, env_var: str) -> tuple[str, Path, str | None]:
    override = os.environ.get(env_var)
    if override:
        path = Path(override).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"{env_var} points to missing file: {path}")
        return path.read_text(encoding="utf-8"), path, override
    path = datasets.get_path(default_relative)
    return path.read_text(encoding="utf-8"), path, None


def _load_file_path(default_relative: str, env_var: str) -> tuple[Path, str | None]:
    override = os.environ.get(env_var)
    if override:
        path = Path(override).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"{env_var} points to missing file: {path}")
        return path, override
    return datasets.get_path(default_relative), None


def _count_fasta_residues(path: Path) -> int:
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def _wrap_fasta(seq: str, width: int = 60) -> str:
    if not seq:
        return ""
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


BASE_DNA_SEQUENCE, DNA_SOURCE_PATH, DNA_ENV = _load_sequence("dna/plasmid_demo.fna", "HELIX_BENCH_DNA_FASTA")
BASE_PROTEIN_PATH, PROTEIN_ENV = _load_file_path("protein/demo_protein.faa", "HELIX_BENCH_PROTEIN_FASTA")
BASE_PROTEIN_RESIDUES = _count_fasta_residues(BASE_PROTEIN_PATH)

ACTIVE_DNA_SEQUENCE = BASE_DNA_SEQUENCE
ACTIVE_DNA_FILE_PATH = DNA_SOURCE_PATH
ACTIVE_RNA_SEQUENCE = ACTIVE_DNA_SEQUENCE.replace("T", "U")
ACTIVE_PROTEIN_PATH = BASE_PROTEIN_PATH
ACTIVE_PROTEIN_RESIDUES = BASE_PROTEIN_RESIDUES
ACTIVE_LIMIT = 0
ACTIVE_DATASET_INFO: Dict[str, Any] = {
    "dna_fasta": str(DNA_SOURCE_PATH),
    "protein_fasta": str(BASE_PROTEIN_PATH),
    "dna_size_bp": len(ACTIVE_DNA_SEQUENCE),
    "protein_residues": BASE_PROTEIN_RESIDUES,
    "limit": None,
}


# Scenario plumbing ----------------------------------------------------------

@dataclass
class Scenario:
    name: str
    runner: Callable[[], Any]
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class BenchmarkResult:
    name: str
    samples: List[float]
    metadata: Dict[str, Any]
    rss_peak_mb: float | None
    status: str = "ok"
    error: str | None = None

    @property
    def mean_s(self) -> float:
        return statistics.fmean(self.samples) if self.samples else 0.0

    @property
    def stdev_s(self) -> float:
        if len(self.samples) < 2:
            return 0.0
        return statistics.stdev(self.samples)

    @property
    def min_s(self) -> float:
        return min(self.samples, default=0.0)

    @property
    def max_s(self) -> float:
        return max(self.samples, default=0.0)


def _run_workflow() -> None:
    config_path = REPO_ROOT / "workflows" / "plasmid_screen.yaml"
    bench_config_dir = REPO_ROOT / ".bench"
    bench_config_dir.mkdir(parents=True, exist_ok=True)
    config_copy = bench_config_dir / "workflow_bench.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    for wf in data.get("workflows", []):
        for step in wf.get("steps", []):
            args = step.get("args", {})
            if ACTIVE_DNA_FILE_PATH:
                if "input" in args:
                    args["input"] = str(ACTIVE_DNA_FILE_PATH)
                if "fasta" in args:
                    args["fasta"] = str(ACTIVE_DNA_FILE_PATH)
    config_copy.write_text(yaml.safe_dump(data), encoding="utf-8")
    with TemporaryDirectory(prefix="helix_bench_workflow_") as tmp:
        tmpdir = Path(tmp)
        hx.run_workflow(config_copy, output_dir=tmpdir, name="plasmid_screen")


def _dna_summary_runner() -> None:
    hx.dna_summary(sequence=ACTIVE_DNA_SEQUENCE, window=200, step=40, k=6, max_diff=1)


def _triage_runner() -> None:
    hx.triage_report(sequence=ACTIVE_DNA_SEQUENCE, k=5, max_diff=1, min_orf_length=60)


def _fold_runner() -> None:
    hx.fold_rna(ACTIVE_RNA_SEQUENCE, min_loop_length=3, allow_wobble_pairs=True)


def _spectrum_runner() -> None:
    hx.spectrum_leaderboard(
        peptide="NQEL",
        experimental_spectrum=[0, 113, 114, 128, 227, 242, 242, 355, 356, 370, 371, 484],
        cyclic=True,
        leaderboard_size=10,
    )


def _protein_runner() -> None:
    hx.protein_summary(input_path=ACTIVE_PROTEIN_PATH, window=11, step=1, scale="kd")


DEF_SCENARIOS: List[Scenario] = [
    Scenario("helix.dna_summary", _dna_summary_runner),
    Scenario("helix.triage_report", _triage_runner),
    Scenario("helix.fold_rna", _fold_runner),
    Scenario("helix.spectrum_leaderboard", _spectrum_runner),
    Scenario("helix.protein_summary", _protein_runner),
    Scenario("workflow.plasmid_screen", _run_workflow, metadata={"config": "workflows/plasmid_screen.yaml", "steps": 4}),
]


# Dataset limiting -----------------------------------------------------------

def _write_limited_sequence(seq: str, dest: Path) -> None:
    dest.write_text(f">helix_bench\n{_wrap_fasta(seq)}\n", encoding="utf-8")


def _prepare_protein_file(limit: int, tmp_dir: Path | None) -> tuple[Path, int]:
    if limit <= 0 or not tmp_dir:
        return BASE_PROTEIN_PATH, BASE_PROTEIN_RESIDUES
    out_path = tmp_dir / "protein_limit.faa"
    remaining = limit
    written = 0
    with BASE_PROTEIN_PATH.open("r", encoding="utf-8") as src, out_path.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                if remaining <= 0:
                    break
                dst.write(line)
                continue
            if remaining <= 0:
                continue
            chunk = line.strip()
            if not chunk:
                continue
            take = chunk[:remaining]
            dst.write(take + "\n")
            written += len(take)
            remaining -= len(take)
            if remaining <= 0:
                break
    if written == 0:
        # Fall back to base dataset if we somehow failed to copy residues.
        return BASE_PROTEIN_PATH, BASE_PROTEIN_RESIDUES
    return out_path, written


def _configure_datasets(limit: int, tmp_dir: Path | None) -> None:
    global ACTIVE_DNA_SEQUENCE, ACTIVE_RNA_SEQUENCE, ACTIVE_DNA_FILE_PATH
    global ACTIVE_PROTEIN_PATH, ACTIVE_PROTEIN_RESIDUES, ACTIVE_LIMIT, ACTIVE_DATASET_INFO

    seq = BASE_DNA_SEQUENCE if limit <= 0 else BASE_DNA_SEQUENCE[:limit]
    ACTIVE_DNA_SEQUENCE = seq
    ACTIVE_RNA_SEQUENCE = seq.replace("T", "U")
    if limit > 0 and tmp_dir is not None:
        dna_path = tmp_dir / "dna_limit.fna"
        _write_limited_sequence(seq, dna_path)
    else:
        dna_path = DNA_SOURCE_PATH
    ACTIVE_DNA_FILE_PATH = dna_path

    protein_path, protein_residues = _prepare_protein_file(limit, tmp_dir)
    ACTIVE_PROTEIN_PATH = protein_path
    ACTIVE_PROTEIN_RESIDUES = protein_residues
    ACTIVE_LIMIT = limit
    ACTIVE_DATASET_INFO = {
        "dna_fasta": str(dna_path),
        "protein_fasta": str(protein_path),
        "dna_size_bp": len(seq),
        "protein_residues": protein_residues,
        "limit": limit if limit > 0 else None,
    }
    _refresh_metadata()


def _refresh_metadata() -> None:
    for scenario in DEF_SCENARIOS:
        if scenario.name == "helix.dna_summary":
            scenario.metadata.update({
                "length_nt": len(ACTIVE_DNA_SEQUENCE),
                "window": 200,
                "step": 40,
                "k": 6,
                "max_diff": 1,
            })
        elif scenario.name == "helix.triage_report":
            scenario.metadata.update({
                "length_nt": len(ACTIVE_DNA_SEQUENCE),
                "k": 5,
                "max_diff": 1,
                "min_orf_length": 60,
            })
        elif scenario.name == "helix.fold_rna":
            scenario.metadata.update({
                "input_length_nt": len(ACTIVE_RNA_SEQUENCE),
                "min_loop_length": 3,
                "allow_wobble_pairs": True,
            })
        elif scenario.name == "helix.spectrum_leaderboard":
            scenario.metadata.update({
                "peptide": "NQEL",
                "experimental_spectrum_len": 11,
                "leaderboard_size": 10,
            })
        elif scenario.name == "helix.protein_summary":
            scenario.metadata.update({
                "window": 11,
                "step": 1,
                "scale": "kd",
                "residues": ACTIVE_PROTEIN_RESIDUES,
            })


_refresh_metadata()


# CLI ------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=5, help="Measured iterations per scenario (default: 5).")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup iterations discarded before timing.")
    parser.add_argument("--scenario", action="append", help="Optional scenario names to run (defaults to all).")
    parser.add_argument("--json", type=Path, help="Deprecated alias for --out.")
    parser.add_argument(
        "--out",
        dest="out",
        type=Path,
        default=Path("bench-results/api_benchmarks.json"),
        help="Output JSON path.",
    )
    parser.add_argument("--summary", dest="summary", type=Path, help="Optional path for Markdown summary.")
    parser.add_argument(
        "--summary-md",
        dest="summary_md",
        type=Path,
        help="Alias for --summary (Markdown).",
    )
    parser.add_argument("--sort", choices=("name", "mean"), default="name", help="Sort results by name or mean time.")
    parser.add_argument("--baseline", type=Path, help="Compare against a baseline JSON to compute deltas.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="RNG seed (default: 1337).")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Maximum nucleotides/amino acids to keep (0 = full dataset).",
    )
    return parser.parse_args()


def _select_scenarios(filter_names: Sequence[str] | None) -> List[Scenario]:
    scenarios = DEF_SCENARIOS
    if filter_names:
        names = set(filter_names)
        scenarios = [scenario for scenario in DEF_SCENARIOS if scenario.name in names]
        if not scenarios:
            raise SystemExit("No matching scenarios found.")
    return scenarios


def _current_rss_mb() -> float | None:
    if not psutil:
        return None
    return psutil.Process().memory_info().rss / (1024 * 1024)


def run_scenario(scenario: Scenario, repeat: int, warmup: int) -> BenchmarkResult:
    rss_peak = _current_rss_mb()
    samples: List[float] = []
    try:
        for _ in range(max(warmup, 0)):
            scenario.runner()
        for _ in range(repeat):
            start = time.perf_counter()
            scenario.runner()
            elapsed = time.perf_counter() - start
            samples.append(elapsed)
            current_rss = _current_rss_mb()
            if current_rss is not None and (rss_peak is None or current_rss > rss_peak):
                rss_peak = current_rss
        return BenchmarkResult(name=scenario.name, samples=samples, metadata=scenario.metadata.copy(), rss_peak_mb=rss_peak)
    except Exception as exc:  # pragma: no cover - surfaced via CLI
        return BenchmarkResult(
            name=scenario.name,
            samples=samples,
            metadata=scenario.metadata.copy(),
            rss_peak_mb=rss_peak,
            status="error",
            error=str(exc),
        )


def _git_info() -> Dict[str, str | None]:
    def _cmd(args: List[str]) -> str | None:
        try:
            return subprocess.check_output(args, cwd=REPO_ROOT).decode().strip()
        except Exception:
            return None

    return {
        "commit": _cmd(["git", "rev-parse", "HEAD"]),
        "branch": _cmd(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
    }


def _collect_env(seed: int, git: Mapping[str, str | None]) -> dict:
    try:
        helix_version = metadata.version("veri-helix")
    except metadata.PackageNotFoundError:
        helix_version = "unknown"

    def _pkg_version(pkg: str) -> str | None:
        try:
            return metadata.version(pkg)
        except metadata.PackageNotFoundError:
            return None

    blas_info = None
    try:
        info = np.__config__.get_info("blas_opt_info")
        libs = info.get("libraries") if isinstance(info, dict) else None
        blas_info = ",".join(libs) if libs else None
    except Exception:  # pragma: no cover
        blas_info = None

    locale_name = "".join(filter(None, locale.getlocale())) or locale.getpreferredencoding(False)

    return {
        "python": platform.python_version(),
        "helix_version": helix_version,
        "numpy": _pkg_version("numpy"),
        "pydantic": _pkg_version("pydantic"),
        "matplotlib": _pkg_version("matplotlib"),
        "platform": platform.platform(),
        "cpu": platform.processor() or platform.machine(),
        "threads": os.cpu_count(),
        "blas": blas_info,
        "locale": locale_name,
        "seed": seed,
        "omp_threads": os.environ.get("OMP_NUM_THREADS"),
        "mkl_threads": os.environ.get("MKL_NUM_THREADS"),
        "git_commit": git.get("commit"),
        "git_branch": git.get("branch"),
    }


def _dataset_metadata() -> dict:
    info = dict(ACTIVE_DATASET_INFO)
    info.setdefault("dna_fasta", str(DNA_SOURCE_PATH))
    info.setdefault("protein_fasta", str(BASE_PROTEIN_PATH))
    return info


def _baseline_map(path: Path | None) -> Dict[str, float]:
    if not path or not path.exists():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        case["name"]: case.get("time_s", {}).get("mean", 0.0)
        for case in data.get("cases", [])
        if case.get("status") == "ok"
    }


def _build_report(
    results: Sequence[BenchmarkResult],
    repeat: int,
    warmup: int,
    seed: int,
    baseline: Mapping[str, float],
) -> dict:
    git = _git_info()
    env = _collect_env(seed, git)
    timestamp = datetime.now(timezone.utc).isoformat()
    cases = []
    for result in results:
        throughput = 1.0 / result.mean_s if result.mean_s else None
        baseline_mean = baseline.get(result.name)
        delta_pct = None
        if baseline_mean and baseline_mean > 0:
            delta_pct = (result.mean_s - baseline_mean) / baseline_mean * 100.0
        cases.append(
            {
                "name": result.name,
                "params": result.metadata,
                "n": len(result.samples),
                "time_s": {
                    "mean": result.mean_s,
                    "std": result.stdev_s,
                    "min": result.min_s,
                    "max": result.max_s,
                },
                "rss_mb": {"peak": result.rss_peak_mb},
                "throughput": {"items_s": throughput} if throughput is not None else None,
                "status": result.status,
                "error": result.error,
                "delta_vs_baseline_pct": delta_pct,
            }
        )
    return {
        "schema": {"kind": SCHEMA_KIND, "spec_version": SCHEMA_VERSION},
        "run": {
            "timestamp": timestamp,
            "commit": git.get("commit"),
            "branch": git.get("branch"),
            "repeat": repeat,
            "warmup": warmup,
            "dataset": _dataset_metadata(),
        },
        "env": env,
        "cases": cases,
    }


def _render_markdown(report: dict, baseline: Mapping[str, float]) -> str:
    env = report.get("env", {})
    run = report.get("run", {})
    dataset = run.get("dataset", {})
    commit = env.get("git_commit") or "unknown"
    dataset_label = Path(dataset.get("dna_fasta", "(default)" )).name
    header = f"### Bench Summary @ {commit}\n"
    header += f"Dataset: {dataset_label} • Repeats: {run.get('repeat')} • Warmup: {run.get('warmup')} • Limit: {dataset.get('limit') or 'full'}\n\n"
    table = ["| Case | mean (s) | std (s) | Δ vs baseline | status |", "| --- | --- | --- | --- | --- |"]
    for case in report.get("cases", []):
        name = case["name"]
        mean = case["time_s"]["mean"]
        std = case["time_s"]["std"]
        delta = case.get("delta_vs_baseline_pct")
        if delta is None and name in baseline and baseline[name] > 0:
            delta = (mean - baseline[name]) / baseline[name] * 100.0
        delta_str = f"{delta:+.1f}%" if delta is not None else "N/A"
        table.append(f"| {name} | {mean:.4f} | {std:.4f} | {delta_str} | {case.get('status')} |")
    return header + "\n".join(table) + "\n"


def main() -> None:
    args = parse_args()
    json_path = args.json or args.out
    summary_path = args.summary or args.summary_md
    limit = max(args.limit, 0)
    random.seed(args.seed)
    try:
        np.random.seed(args.seed)
    except Exception:  # pragma: no cover
        pass

    with TemporaryDirectory(prefix="helix_bench_data_") as data_tmp:
        tmp_dir = Path(data_tmp)
        _configure_datasets(limit, tmp_dir)
        scenarios = _select_scenarios(args.scenario)
        results: List[BenchmarkResult] = []
        for scenario in scenarios:
            results.append(run_scenario(scenario, repeat=args.repeat, warmup=args.warmup))

        if args.sort == "mean":
            results.sort(key=lambda r: r.mean_s, reverse=True)
        else:
            results.sort(key=lambda r: r.name)

        baseline = _baseline_map(args.baseline)
        report = _build_report(results, args.repeat, args.warmup, args.seed, baseline)

        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

        summary = _render_markdown(report, baseline)
        print(summary)
        if summary_path:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(summary, encoding="utf-8")
        else:
            table_lines = [
                f"{res.name}: mean={res.mean_s:.4f}s stdev={res.stdev_s:.4f}s status={res.status}"
                for res in results
            ]
            print("\n".join(table_lines))

        if any(res.status != "ok" for res in results):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
