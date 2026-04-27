#!/usr/bin/env python3
"""Compare DQMC scalar output files against documented ED reference values."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable


def read_series(path: Path) -> list[float]:
    if not path.exists():
        raise FileNotFoundError(f"missing DQMC output file: {path}")
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        values.append(float(stripped.split()[0]))
    if not values:
        raise ValueError(f"no numeric values found in {path}")
    return values


def last_values(run_dir: Path, files: Iterable[str]) -> list[float]:
    return [read_series(run_dir / name)[-1] for name in files]


def compute_actual(run_dir: Path, operation: str, files: list[str], lq: int) -> float:
    values = last_values(run_dir, files)
    if operation == "last":
        if len(values) != 1:
            raise ValueError("operation 'last' requires exactly one file")
        return values[0]
    if operation == "sum_last":
        return sum(values)
    if operation == "last_times_lq":
        if len(values) != 1:
            raise ValueError("operation 'last_times_lq' requires exactly one file")
        return values[0] * float(lq)
    raise ValueError(f"unsupported comparison operation: {operation}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="JSON reference file")
    parser.add_argument("--run-dir", required=True, type=Path, help="Directory containing DQMC scalar output files")
    args = parser.parse_args()

    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    params = reference["parameters"]
    lq = int(params["Lx"]) * int(params["Ly"])

    failures: list[str] = []
    print(f"Reference: {args.reference}")
    print(f"DQMC run directory: {args.run_dir}")
    print("")

    for name, spec in reference["observables"].items():
        dqmc = spec["dqmc"]
        actual = compute_actual(args.run_dir, dqmc["operation"], dqmc["files"], lq)
        expected = float(spec["value"])
        atol = float(spec.get("atol", 0.0))
        rtol = float(spec.get("rtol", 0.0))
        passed = math.isclose(actual, expected, abs_tol=atol, rel_tol=rtol)
        status = "PASS" if passed else "FAIL"
        print(
            f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
            f"abs_diff={abs(actual - expected):.3g} atol={atol:.3g} rtol={rtol:.3g}"
        )
        if not passed:
            failures.append(name)

    if failures:
        print("")
        print("Failed observables: " + ", ".join(failures))
        return 1

    print("")
    print("All benchmark comparisons passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
