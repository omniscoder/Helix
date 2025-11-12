#!/usr/bin/env python3
"""Compare benchmark JSON outputs and flag regressions."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict


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
    for name in sorted(set(baseline) & set(current)):
        old = baseline[name]
        new = current[name]
        if old <= 0:
            continue
        change = (new - old) / old * 100.0
        if change > args.threshold:
            failures.append((name, old, new, change))

    if failures:
        print(f"Benchmark regressions (>{args.threshold:.1f}% slowdown):")
        for name, old, new, change in failures:
            print(f" - {name}: {old:.4f}s -> {new:.4f}s ({change:+.1f}%)")
        sys.exit(1)

    print(f"No significant regressions (>{args.threshold:.1f}%).")


if __name__ == "__main__":
    main()
