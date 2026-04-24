#!/usr/bin/env python3
"""Regression checks for the matched native-mesh/OpenVDB baseline."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def safe_ratio(num: float, den: float) -> float:
    if den <= 1.0e-14:
        return math.inf if num > 1.0e-14 else 0.0
    return num / den


def require(condition: bool, name: str, details: str, results: list[dict[str, str]]) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"{status} {name}: {details}")
    results.append({"check": name, "status": status, "details": details})
    if not condition:
        raise AssertionError(f"{name}: {details}")


def check(args: argparse.Namespace) -> list[dict[str, str]]:
    project_root: Path = args.project_root.resolve()
    comparison_path = args.comparison or project_root / "out" / "milestone_27" / "matched_mesh_comparison.csv"
    rows = read_csv(comparison_path)
    if not rows:
        raise AssertionError(f"{comparison_path} has no rows")

    results: list[dict[str, str]] = []
    scenario_names = [row.get("scenario", "") for row in rows]

    require(
        len(rows) >= args.min_scenarios,
        "matched_mesh_scenario_count",
        f"{len(rows)} scenarios present ({', '.join(scenario_names)}) >= {args.min_scenarios}",
        results,
    )

    min_count_ratio = min(
        safe_ratio(f(row, "native_mean_abs_contact_count_change"), f(row, "field_mean_abs_patch_count_change"))
        for row in rows
    )
    require(
        min_count_ratio > args.count_churn_ratio_threshold,
        "native_contact_count_churn_larger",
        f"min native/field mean count-churn ratio={min_count_ratio:.8f} "
        f"> {args.count_churn_ratio_threshold:.8f}",
        results,
    )

    min_normal_ratio = min(f(row, "native_contact_normal_to_field_normal_jump_ratio") for row in rows)
    require(
        min_normal_ratio > args.normal_jump_ratio_threshold,
        "field_normal_jump_lower_than_native_contact_normal",
        f"min native contact-normal jump / field normal-jump ratio={min_normal_ratio:.8f} "
        f"> {args.normal_jump_ratio_threshold:.8f}",
        results,
    )

    max_coulomb = max(f(row, "field_max_tangential_force_ratio") for row in rows)
    require(
        max_coulomb <= 1.0 + args.tolerance,
        "field_coulomb_bound",
        f"max field Coulomb ratio={max_coulomb:.8f} <= {1.0 + args.tolerance:.8f}",
        results,
    )

    max_pos_error = max(f(row, "native_max_position_error") for row in rows)
    max_vel_error = max(f(row, "native_max_velocity_error") for row in rows)
    max_acc_error = max(f(row, "native_max_acceleration_error") for row in rows)
    require(
        max_pos_error <= args.max_position_error,
        "native_position_tracking_controllable",
        f"max position error={max_pos_error:.8e} <= {args.max_position_error:.8e}",
        results,
    )
    require(
        max_vel_error <= args.max_velocity_error,
        "native_velocity_tracking_controllable",
        f"max velocity error={max_vel_error:.8e} <= {args.max_velocity_error:.8e}",
        results,
    )
    require(
        max_acc_error <= args.max_acceleration_error,
        "native_acceleration_tracking_controllable",
        f"max acceleration error={max_acc_error:.8e} <= {args.max_acceleration_error:.8e}",
        results,
    )

    return results


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--min-scenarios", type=int, default=2)
    parser.add_argument("--count-churn-ratio-threshold", type=float, default=10.0)
    parser.add_argument("--normal-jump-ratio-threshold", type=float, default=1.5)
    parser.add_argument("--max-position-error", type=float, default=3.5e-3)
    parser.add_argument("--max-velocity-error", type=float, default=1.0e-3)
    parser.add_argument("--max-acceleration-error", type=float, default=1.0e-1)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=default_root / "out" / "milestone_27" / "matched_mesh_regression_summary.json",
    )
    args = parser.parse_args()

    try:
        results = check(args)
    except (AssertionError, FileNotFoundError) as exc:
        print(f"MATCHED MESH REGRESSION FAILED: {exc}", file=sys.stderr)
        return 1

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.summary_output}")
    print("MATCHED MESH REGRESSION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
