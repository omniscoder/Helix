#!/usr/bin/env python3
"""Compare benchmark JSON outputs and flag regressions.

Adds guardrails for micro-bench noise: only consider cases whose baseline
runtime exceeds a minimum, and whose absolute slowdown exceeds a minimum,
with optional include/exclude filters.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Iterable


def _load_cases(path: Path) -> Dict[str, float]:
    data = json.loads(path.read_text())
    cases = {}
    for case in data.get("cases", []):
        if case.get("status") != "ok":
            continue
        mean = case.get("time_s", {}).get("mean")
        if mean is not None:
            cases[case["name"]] = mean
    return cases


def _filter(names: Iterable[str], include: Iterable[str] | None, exclude: Iterable[str] | None) -> set[str]:
    selected = set(names)
    if include:
        include_set = set(include)
        selected = {n for n in selected if n in include_set}
    if exclude:
        exclude_set = set(exclude)
        selected = {n for n in selected if n not in exclude_set}
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path, help="Path to baseline bench JSON")
    parser.add_argument("current", type=Path, help="Path to current bench JSON")
    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help="Percent slowdown tolerated before failing (default: 5.0).",
    )
    parser.add_argument(
        "--min-abs-slowdown",
        type=float,
        default=0.002,
        help="Only flag if absolute slowdown exceeds this many seconds (default: 0.002s).",
    )
    parser.add_argument(
        "--min-runtime",
        type=float,
        default=0.01,
        help="Only consider cases with baseline mean >= this many seconds (default: 0.01s).",
    )
    parser.add_argument("--include", action="append", help="Restrict checks to specific case names (repeatable).")
    parser.add_argument("--exclude", action="append", help="Skip specific case names (repeatable).")
    args = parser.parse_args()

    if not args.baseline.exists():
        print(f"Baseline '{args.baseline}' not found; skipping drift check.")
        return
    if not args.current.exists():
        raise SystemExit(f"Current benchmark file '{args.current}' not found.")

    baseline = _load_cases(args.baseline)
    current = _load_cases(args.current)

    if not baseline:
        print("Baseline has no comparable cases; skipping drift check.")
        return

    failures = []
    skipped = []
    candidates = _filter(sorted(set(baseline) & set(current)), args.include, args.exclude)
    for name in candidates:
        old = baseline[name]
        new = current[name]
        if old <= 0:
            continue
        abs_slowdown = new - old
        change = (abs_slowdown) / old * 100.0
        # Skip micro-bench noise unless both runtime and absolute delta exceed thresholds
        if old < args.min_runtime or abs_slowdown < args.min_abs_slowdown:
            skipped.append((name, old, new, change))
            continue
        if change > args.threshold:
            failures.append((name, old, new, change))

    if failures:
        print(f"Benchmark regressions (>{args.threshold:.1f}% slowdown):")
        for name, old, new, change in failures:
            print(f" - {name}: {old:.4f}s -> {new:.4f}s ({change:+.1f}%)")
        sys.exit(1)
    if skipped:
        print(
            f"Skipped {len(skipped)} case(s) due to min-runtime ({args.min_runtime}s) "
            f"or min-abs-slowdown ({args.min_abs_slowdown}s) guardrails."
        )
    print(f"No significant regressions (>{args.threshold:.1f}%).")


if __name__ == "__main__":
    main()
