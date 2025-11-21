#!/usr/bin/env python3
"""Compare current benchmark results against a baseline."""

from __future__ import annotations

import json
import sys
from typing import Dict, Tuple, Any

THRESHOLD = 0.85  # allow up to 15% slowdown


def load(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _index(entries, key_fields):
    indexed = {}
    for entry in entries:
        key = tuple(entry.get(field) for field in key_fields)
        indexed[key] = entry
    return indexed


def _compare(current: Dict[str, Any], baseline: Dict[str, Any], *, field: str) -> list[Tuple[Tuple[str, ...], float, float, float]]:
    failures = []
    for key, base_entry in baseline.items():
        base_value = float(base_entry.get(field, 0.0))
        if base_value <= 0:
            continue
        cur_entry = current.get(key)
        if not cur_entry:
            continue
        cur_value = float(cur_entry.get(field, 0.0))
        ratio = cur_value / base_value if base_value else 1.0
        if ratio < THRESHOLD:
            failures.append((key, base_value, cur_value, ratio))
    return failures


def main(current_path: str, baseline_path: str) -> int:
    current = load(current_path)
    baseline = load(baseline_path)

    cur_crispr = _index(current["benchmarks"]["crispr"], ["backend_used", "shape"])
    base_crispr = _index(baseline["benchmarks"]["crispr"], ["backend_used", "shape"])

    cur_prime = _index(current["benchmarks"]["prime"], ["backend_used", "workload"])
    base_prime = _index(baseline["benchmarks"]["prime"], ["backend_used", "workload"])

    failures = []
    failures.extend(_compare(cur_crispr, base_crispr, field="mpairs_per_s"))
    failures.extend(_compare(cur_prime, base_prime, field="predictions_per_s"))

    if failures:
        print("Benchmark regressions detected:")
        for key, base_value, cur_value, ratio in failures:
            print(
                f"  {key}: baseline={base_value:.3f}, current={cur_value:.3f}, ratio={ratio:.3f}"
            )
        return 1

    print("No benchmark regressions beyond threshold.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: check_benchmark_regression.py CURRENT BASELINE")
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
