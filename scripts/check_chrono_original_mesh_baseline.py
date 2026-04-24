#!/usr/bin/env python3
"""Regression checks for the Chrono-original mesh/OpenVDB baseline."""

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
    comparison_path = (
        args.comparison or project_root / "out" / "milestone_31" / "chrono_original_mesh_comparison.csv"
    )
    geometry_path = args.geometry or project_root / "out" / "milestone_31" / "chrono_original_mesh_geometry.csv"

    comparison = read_csv(comparison_path)
    geometry = read_csv(geometry_path)
    if not comparison:
        raise AssertionError(f"{comparison_path} has no rows")
    if not geometry:
        raise AssertionError(f"{geometry_path} has no rows")

    results: list[dict[str, str]] = []
    scenario_names = [row.get("scenario", "") for row in comparison]

    require(
        len(comparison) >= args.min_scenarios and len(geometry) >= args.min_scenarios,
        "chrono_original_mesh_scenario_count",
        f"{len(comparison)} comparison rows and {len(geometry)} geometry rows >= {args.min_scenarios}",
        results,
    )

    missing_assets = []
    for row in geometry:
        target = Path(row.get("target_obj", ""))
        source_demo = row.get("source_demo", "")
        if not target.exists() or "demo_MBS_" not in source_demo:
            missing_assets.append(f"{row.get('scenario', '')}: target={target}, source={source_demo}")
    require(
        not missing_assets,
        "chrono_original_assets_recorded",
        "all target OBJ assets exist and reference original Chrono MBS demos",
        results,
    )

    min_faces = min(f(row, "target_faces") for row in geometry)
    min_voxels = min(f(row, "active_voxels") for row in geometry)
    require(
        min_faces >= args.min_target_faces,
        "target_mesh_resolution_nontrivial",
        f"min target faces={min_faces:.0f} >= {args.min_target_faces:.0f}",
        results,
    )
    require(
        min_voxels >= args.min_active_voxels,
        "openvdb_sdf_backend_active",
        f"min OpenVDB active voxels={min_voxels:.0f} >= {args.min_active_voxels:.0f}",
        results,
    )

    min_field_active = min(f(row, "field_active_frames") for row in comparison)
    min_native_active = min(f(row, "native_active_frames") for row in comparison)
    require(
        min_field_active >= args.min_field_active_frames and min_native_active >= args.min_native_active_frames,
        "both_methods_make_contact",
        f"min field/native active frames={min_field_active:.0f}/{min_native_active:.0f}",
        results,
    )

    min_count_ratio = min(f(row, "native_to_field_count_churn_ratio") for row in comparison)
    min_normal_ratio = min(f(row, "native_contact_normal_to_field_normal_jump_ratio") for row in comparison)
    min_force_ratio = min(f(row, "native_to_field_force_oscillation_ratio") for row in comparison)
    min_torque_ratio = min(f(row, "native_to_field_torque_oscillation_ratio") for row in comparison)
    require(
        math.isfinite(min_count_ratio) and min_count_ratio >= args.min_count_churn_ratio,
        "native_contact_count_churn_larger",
        f"min native/field count churn ratio={min_count_ratio:.8f} >= {args.min_count_churn_ratio:.8f}",
        results,
    )
    require(
        math.isfinite(min_normal_ratio) and min_normal_ratio >= args.min_normal_jump_ratio,
        "field_normal_jump_lower",
        f"min native/field normal jump ratio={min_normal_ratio:.8f} >= {args.min_normal_jump_ratio:.8f}",
        results,
    )
    require(
        math.isfinite(min_force_ratio) and min_force_ratio >= args.min_force_oscillation_ratio,
        "field_force_oscillation_lower",
        f"min native/field force oscillation ratio={min_force_ratio:.8f} >= {args.min_force_oscillation_ratio:.8f}",
        results,
    )
    require(
        math.isfinite(min_torque_ratio) and min_torque_ratio >= args.min_torque_oscillation_ratio,
        "field_torque_oscillation_lower",
        f"min native/field torque oscillation ratio={min_torque_ratio:.8f} >= {args.min_torque_oscillation_ratio:.8f}",
        results,
    )

    max_coulomb = max(f(row, "field_max_tangential_force_ratio") for row in comparison)
    require(
        max_coulomb <= 1.0 + args.tolerance,
        "field_coulomb_bound",
        f"max field Coulomb ratio={max_coulomb:.8f} <= {1.0 + args.tolerance:.8f}",
        results,
    )

    max_pos_error = max(f(row, "native_max_position_error") for row in comparison)
    max_vel_error = max(f(row, "native_max_velocity_error") for row in comparison)
    max_acc_error = max(f(row, "native_max_acceleration_error") for row in comparison)
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

    print("scenarios: " + ", ".join(scenario_names))
    print(
        "derived min native/field count-churn check="
        f"{min(safe_ratio(f(row, 'native_mean_abs_contact_count_change'), f(row, 'field_mean_abs_patch_count_change')) for row in comparison):.8f}"
    )
    return results


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--min-scenarios", type=int, default=2)
    parser.add_argument("--min-target-faces", type=float, default=60.0)
    parser.add_argument("--min-active-voxels", type=float, default=6.0e5)
    parser.add_argument("--min-field-active-frames", type=float, default=60.0)
    parser.add_argument("--min-native-active-frames", type=float, default=35.0)
    parser.add_argument("--min-count-churn-ratio", type=float, default=40.0)
    parser.add_argument("--min-normal-jump-ratio", type=float, default=3.0)
    parser.add_argument("--min-force-oscillation-ratio", type=float, default=5.0)
    parser.add_argument("--min-torque-oscillation-ratio", type=float, default=3.0)
    parser.add_argument("--max-position-error", type=float, default=2.0e-3)
    parser.add_argument("--max-velocity-error", type=float, default=1.0e-4)
    parser.add_argument("--max-acceleration-error", type=float, default=8.0e-1)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=default_root / "out" / "milestone_31" / "chrono_original_mesh_regression_summary.json",
    )
    args = parser.parse_args()

    try:
        results = check(args)
    except (AssertionError, FileNotFoundError) as exc:
        print(f"CHRONO ORIGINAL MESH REGRESSION FAILED: {exc}", file=sys.stderr)
        return 1

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.summary_output}")
    print("CHRONO ORIGINAL MESH REGRESSION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
