#!/usr/bin/env python3
"""Regression checks for Milestone 28 free-dynamics field contact response."""

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


def require(condition: bool, name: str, details: str, results: list[dict[str, str]]) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"{status} {name}: {details}")
    results.append({"check": name, "status": status, "details": details})
    if not condition:
        raise AssertionError(f"{name}: {details}")


def finite(value: float) -> bool:
    return math.isfinite(value)


def check(args: argparse.Namespace) -> list[dict[str, str]]:
    project_root: Path = args.project_root.resolve()
    comparison_path = args.comparison or project_root / "out" / "milestone_28" / "free_dynamics_comparison.csv"
    summary_path = args.summary or project_root / "out" / "milestone_28" / "free_dynamics_summary.csv"
    comparison = read_csv(comparison_path)
    summary = read_csv(summary_path)
    if not comparison:
        raise AssertionError(f"{comparison_path} has no rows")
    if not summary:
        raise AssertionError(f"{summary_path} has no rows")

    results: list[dict[str, str]] = []
    scenario_names = [row.get("scenario", "") for row in comparison]

    require(
        len(comparison) >= args.min_scenarios,
        "free_dynamics_scenario_count",
        f"{len(comparison)} scenarios present ({', '.join(scenario_names)}) >= {args.min_scenarios}",
        results,
    )

    max_energy_budget = max(
        max(f(row, "field_energy_budget_ratio"), f(row, "native_energy_budget_ratio")) for row in comparison
    )
    require(
        finite(max_energy_budget) and max_energy_budget <= args.max_energy_budget_ratio,
        "no_nonphysical_energy_explosion",
        f"max field/native energy budget ratio={max_energy_budget:.8f} <= {args.max_energy_budget_ratio:.8f}",
        results,
    )

    max_field_coulomb = max(f(row, "field_max_coulomb_ratio") for row in comparison)
    max_native_coulomb = max(f(row, "native_max_coulomb_ratio") for row in comparison)
    require(
        max_field_coulomb <= 1.0 + args.tolerance and max_native_coulomb <= 1.0 + args.tolerance,
        "coulomb_bound",
        f"field max={max_field_coulomb:.8f}, native max={max_native_coulomb:.8f} <= {1.0 + args.tolerance:.8f}",
        results,
    )

    max_field_pen = max(f(row, "field_max_penetration") for row in comparison)
    max_any_pen = max(max(f(row, "field_max_penetration"), f(row, "native_max_penetration")) for row in comparison)
    require(
        max_field_pen <= args.max_field_penetration,
        "field_penetration_bounded",
        f"max field penetration={max_field_pen:.8f} <= {args.max_field_penetration:.8f}",
        results,
    )
    require(
        max_any_pen <= args.max_any_penetration,
        "all_penetration_bounded",
        f"max field/native penetration={max_any_pen:.8f} <= {args.max_any_penetration:.8f}",
        results,
    )

    max_smoothness_ratio = max(f(row, "field_to_native_acceleration_jump_ratio") for row in comparison)
    require(
        max_smoothness_ratio <= args.max_acceleration_jump_ratio,
        "field_trajectory_smoothness_not_worse",
        f"max field/native acceleration-jump ratio={max_smoothness_ratio:.8f} "
        f"<= {args.max_acceleration_jump_ratio:.8f}",
        results,
    )

    regular_rows = [row for row in comparison if not row.get("scenario", "").startswith("free_chrono_original_")]
    original_mesh_rows = [row for row in comparison if row.get("scenario", "").startswith("free_chrono_original_")]

    max_final_delta = max(f(row, "trajectory_final_position_delta") for row in regular_rows or comparison)
    require(
        max_final_delta <= args.max_final_position_delta,
        "trajectory_response_comparable_for_controlled_meshes",
        f"max regular-scene final field/native position delta={max_final_delta:.8f} "
        f"<= {args.max_final_position_delta:.8f}",
        results,
    )

    if original_mesh_rows:
        max_original_delta = max(f(row, "trajectory_final_position_delta") for row in original_mesh_rows)
        require(
            max_original_delta <= args.max_original_mesh_final_position_delta,
            "original_mesh_trajectory_divergence_bounded",
            f"max original-mesh final field/native position delta={max_original_delta:.8f} "
            f"<= {args.max_original_mesh_final_position_delta:.8f}",
            results,
        )

        max_original_penetration_ratio = max(
            f(row, "field_max_penetration") / max(1.0e-14, f(row, "native_max_penetration"))
            for row in original_mesh_rows
        )
        require(
            max_original_penetration_ratio <= args.max_original_mesh_penetration_ratio,
            "original_mesh_field_penetration_lower",
            f"max original-mesh field/native penetration ratio={max_original_penetration_ratio:.8f} "
            f"<= {args.max_original_mesh_penetration_ratio:.8f}",
            results,
        )

        max_original_accel_ratio = max(f(row, "field_to_native_acceleration_jump_ratio") for row in original_mesh_rows)
        require(
            max_original_accel_ratio <= args.max_original_mesh_acceleration_jump_ratio,
            "original_mesh_field_acceleration_jump_lower",
            f"max original-mesh field/native acceleration-jump ratio={max_original_accel_ratio:.8f} "
            f"<= {args.max_original_mesh_acceleration_jump_ratio:.8f}",
            results,
        )

        max_original_count_change_ratio = max(
            f(row, "field_mean_abs_count_change") / max(1.0e-14, f(row, "native_mean_abs_count_change"))
            for row in original_mesh_rows
        )
        require(
            max_original_count_change_ratio <= args.max_original_mesh_count_change_ratio,
            "original_mesh_field_contact_churn_lower",
            f"max original-mesh field/native count-change ratio={max_original_count_change_ratio:.8f} "
            f"<= {args.max_original_mesh_count_change_ratio:.8f}",
            results,
        )

    max_field_contact_work = max(
        f(row, "contact_work")
        for row in summary
        if row.get("variant") == "field_primitive_accumulator"
    )
    require(
        max_field_contact_work <= args.max_positive_net_contact_work,
        "field_contact_net_work_noninjecting",
        f"max field net contact work={max_field_contact_work:.8f} <= {args.max_positive_net_contact_work:.8f}",
        results,
    )

    min_field_active_ratio = min(
        f(row, "active_frames") / max(1.0, f(row, "frames"))
        for row in summary
        if row.get("variant") == "field_primitive_accumulator"
    )
    require(
        min_field_active_ratio >= args.min_field_active_ratio,
        "field_contact_remains_active",
        f"min field active-frame ratio={min_field_active_ratio:.8f} >= {args.min_field_active_ratio:.8f}",
        results,
    )

    max_stick_slip_ratio = max(
        f(row, "field_stick_slip_switches") / max(1.0, f(row, "native_stick_slip_switches")) for row in comparison
    )
    require(
        max_stick_slip_ratio <= args.max_stick_slip_switch_ratio,
        "field_stick_slip_sequence_less_chattery",
        f"max field/native stick-slip switch ratio={max_stick_slip_ratio:.8f} "
        f"<= {args.max_stick_slip_switch_ratio:.8f}",
        results,
    )

    return results


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--min-scenarios", type=int, default=2)
    parser.add_argument("--max-energy-budget-ratio", type=float, default=0.75)
    parser.add_argument("--max-field-penetration", type=float, default=0.020)
    parser.add_argument("--max-any-penetration", type=float, default=0.070)
    parser.add_argument("--max-acceleration-jump-ratio", type=float, default=0.15)
    parser.add_argument("--max-final-position-delta", type=float, default=0.080)
    parser.add_argument("--max-original-mesh-final-position-delta", type=float, default=0.120)
    parser.add_argument("--max-original-mesh-penetration-ratio", type=float, default=0.25)
    parser.add_argument("--max-original-mesh-acceleration-jump-ratio", type=float, default=0.08)
    parser.add_argument("--max-original-mesh-count-change-ratio", type=float, default=0.01)
    parser.add_argument("--max-positive-net-contact-work", type=float, default=0.02)
    parser.add_argument("--min-field-active-ratio", type=float, default=0.95)
    parser.add_argument("--max-stick-slip-switch-ratio", type=float, default=0.10)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=default_root / "out" / "milestone_28" / "free_dynamics_regression_summary.json",
    )
    args = parser.parse_args()

    try:
        results = check(args)
    except (AssertionError, FileNotFoundError) as exc:
        print(f"FREE DYNAMICS REGRESSION FAILED: {exc}", file=sys.stderr)
        return 1

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.summary_output}")
    print("FREE DYNAMICS REGRESSION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
