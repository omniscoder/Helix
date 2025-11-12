#!/usr/bin/env python3
"""Append benchmark metrics to docs/data/bench/history.csv."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List

BASE_COLUMNS = ["sha", "timestamp", "python", "cpu", "repeat", "limit", "bench_heavy", "dataset"]


def _read_rows(path: Path) -> tuple[List[Dict[str, str]], List[str]]:
    if not path.exists():
        return [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return rows, reader.fieldnames or []


def _write_rows(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bench_json", type=Path)
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--sha", required=True)
    parser.add_argument("--limit", default="0")
    parser.add_argument("--bench-heavy", dest="bench_heavy", default="false")
    args = parser.parse_args()

    payload = json.loads(args.bench_json.read_text(encoding="utf-8"))
    case_details = payload.get("cases", [])
    run = payload.get("run", {})
    env = payload.get("env", {})
    dataset = run.get("dataset", {})

    row: Dict[str, str] = {
        "sha": args.sha,
        "timestamp": run.get("timestamp", ""),
        "python": env.get("python", ""),
        "cpu": env.get("cpu", ""),
        "repeat": str(run.get("repeat", "")),
        "limit": str(run.get("dataset", {}).get("limit") or args.limit),
        "bench_heavy": args.bench_heavy,
        "dataset": dataset.get("dna_fasta", ""),
    }
    for case in case_details:
        name = case["name"]
        time_info = case.get("time_s", {}) or {}
        mean = time_info.get("mean")
        delta = case.get("delta_vs_baseline_pct")
        rss = (case.get("rss_mb") or {}).get("peak")
        if mean is not None:
            row[f"{name}.mean_s"] = f"{mean:.6f}"
        if delta is not None:
            row[f"{name}.delta_pct"] = f"{delta:.2f}"
        if rss is not None:
            row[f"{name}.rss_mb"] = f"{rss:.2f}"

    existing_rows, existing_fields = _read_rows(args.csv_path)
    all_fields: List[str] = []
    for column in BASE_COLUMNS:
        if column not in all_fields:
            all_fields.append(column)
    for column in existing_fields:
        if column not in all_fields:
            all_fields.append(column)
    new_columns = [col for col in row.keys() if col not in all_fields]
    new_columns.sort()
    for column in new_columns:
        all_fields.append(column)

    # Ensure every row has every field.
    for existing in existing_rows:
        for column in all_fields:
            existing.setdefault(column, "")
    for column in all_fields:
        row.setdefault(column, "")
    existing_rows.append(row)

    _write_rows(args.csv_path, all_fields, existing_rows)


if __name__ == "__main__":
    main()
