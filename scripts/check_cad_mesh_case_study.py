#!/usr/bin/env python3
"""Regression checks for Milestone 29 real CAD/mesh case study."""

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
    comparison_path = args.comparison or project_root / "out" / "milestone_29" / "cad_mesh_case_study_comparison.csv"
    geometry_path = args.geometry or project_root / "out" / "milestone_29" / "cad_mesh_case_study_geometry.csv"
    comparison = read_csv(comparison_path)
    geometry = read_csv(geometry_path)
    if not comparison:
        raise AssertionError(f"{comparison_path} has no rows")
    if not geometry:
        raise AssertionError(f"{geometry_path} has no rows")

    results: list[dict[str, str]] = []
    row = comparison[0]
    geom = geometry[0]
    scenario = row.get("scenario", "")

    target = Path(row.get("target_obj", ""))
    moving = Path(row.get("moving_obj", ""))
    require(
        target.exists() and moving.exists() and "data" in target.parts and "data" in moving.parts,
        "real_mesh_assets_exist",
        f"target={target}, moving={moving}",
        results,
    )

    require(
        f(row, "target_faces") >= args.min_target_faces and f(row, "moving_faces") >= args.min_moving_faces,
        "mesh_resolution_nontrivial",
        f"{scenario}: target_faces={f(row, 'target_faces'):.0f}, moving_faces={f(row, 'moving_faces'):.0f}",
        results,
    )

    active_voxels = f(row, "openvdb_active_voxels")
    require(
        active_voxels >= args.min_openvdb_active_voxels,
        "openvdb_sdf_backend_active",
        f"active voxels={active_voxels:.0f} >= {args.min_openvdb_active_voxels:.0f}",
        results,
    )

    field_active = f(row, "field_active_frames")
    native_active = f(row, "native_active_frames")
    require(
        field_active >= args.min_field_active_frames and native_active >= args.min_native_active_frames,
        "both_methods_contact",
        f"field_active={field_active:.0f}, native_active={native_active:.0f}",
        results,
    )

    max_field_pen = f(row, "field_max_penetration")
    max_native_pen = f(row, "native_max_penetration")
    require(
        max_field_pen >= args.min_field_penetration and max_field_pen <= args.max_field_penetration,
        "field_penetration_controlled",
        f"field max penetration={max_field_pen:.8f} in [{args.min_field_penetration:.8f}, "
        f"{args.max_field_penetration:.8f}]",
        results,
    )
    require(
        max_native_pen <= args.max_native_penetration,
        "native_penetration_controlled",
        f"native max penetration={max_native_pen:.8f} <= {args.max_native_penetration:.8f}",
        results,
    )

    max_coulomb = max(f(row, "field_max_coulomb_ratio"), f(row, "native_max_coulomb_ratio"))
    require(
        max_coulomb <= 1.0 + args.tolerance,
        "coulomb_bound",
        f"max field/native Coulomb ratio={max_coulomb:.8f} <= {1.0 + args.tolerance:.8f}",
        results,
    )

    count_ratio = f(row, "native_to_field_count_churn_ratio")
    normal_ratio = f(row, "native_to_field_normal_jump_ratio")
    force_ratio = f(row, "native_to_field_force_oscillation_ratio")
    torque_ratio = f(row, "native_to_field_torque_oscillation_ratio")
    require(
        finite(count_ratio) and count_ratio >= args.min_count_churn_ratio,
        "native_contact_count_churn_larger",
        f"native/field count churn ratio={count_ratio:.8f} >= {args.min_count_churn_ratio:.8f}",
        results,
    )
    require(
        finite(normal_ratio) and normal_ratio >= args.min_normal_jump_ratio,
        "field_normal_jump_lower",
        f"native/field normal jump ratio={normal_ratio:.8f} >= {args.min_normal_jump_ratio:.8f}",
        results,
    )
    require(
        finite(force_ratio) and force_ratio >= args.min_force_oscillation_ratio,
        "field_force_oscillation_lower",
        f"native/field force oscillation ratio={force_ratio:.8f} >= {args.min_force_oscillation_ratio:.8f}",
        results,
    )
    require(
        finite(torque_ratio) and torque_ratio >= args.min_torque_oscillation_ratio,
        "field_torque_oscillation_lower",
        f"native/field torque oscillation ratio={torque_ratio:.8f} >= {args.min_torque_oscillation_ratio:.8f}",
        results,
    )

    moving_scale = f(geom, "moving_scale")
    require(
        moving_scale > 0.0,
        "moving_mesh_scale_recorded",
        f"moving scale={moving_scale:.8f}",
        results,
    )

    return results


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--min-target-faces", type=float, default=60.0)
    parser.add_argument("--min-moving-faces", type=float, default=1000.0)
    parser.add_argument("--min-openvdb-active-voxels", type=float, default=100000.0)
    parser.add_argument("--min-field-active-frames", type=float, default=60.0)
    parser.add_argument("--min-native-active-frames", type=float, default=8.0)
    parser.add_argument("--min-field-penetration", type=float, default=5.0e-4)
    parser.add_argument("--max-field-penetration", type=float, default=8.0e-3)
    parser.add_argument("--max-native-penetration", type=float, default=8.0e-3)
    parser.add_argument("--min-count-churn-ratio", type=float, default=3.0)
    parser.add_argument("--min-normal-jump-ratio", type=float, default=1.5)
    parser.add_argument("--min-force-oscillation-ratio", type=float, default=3.0)
    parser.add_argument("--min-torque-oscillation-ratio", type=float, default=2.5)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=default_root / "out" / "milestone_29" / "cad_mesh_case_study_regression_summary.json",
    )
    args = parser.parse_args()

    try:
        results = check(args)
    except (AssertionError, FileNotFoundError) as exc:
        print(f"CAD MESH CASE STUDY REGRESSION FAILED: {exc}", file=sys.stderr)
        return 1

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.summary_output}")
    print("CAD MESH CASE STUDY REGRESSION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
