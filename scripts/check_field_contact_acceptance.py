#!/usr/bin/env python3
"""Acceptance-level regression checks for field-contact paper conclusions."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
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


def find_row(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    for row in rows:
        if all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise AssertionError(f"missing row matching {criteria}")


def run_command(cmd: list[str], cwd: Path) -> None:
    print("RUN " + " ".join(str(part) for part in cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def find_demo(project_root: Path, name: str) -> Path:
    suffix = ".exe" if sys.platform.startswith("win") else ""
    candidates = [
        project_root / "build" / "bin" / "Release" / f"{name}{suffix}",
        project_root / "build" / "bin" / "Debug" / f"{name}{suffix}",
        project_root / "build" / "bin" / f"{name}{suffix}",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"cannot locate {name}; checked: {', '.join(str(c) for c in candidates)}")


def require(condition: bool, name: str, details: str, results: list[dict[str, str]]) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"{status} {name}: {details}")
    results.append({"check": name, "status": status, "details": details})
    if not condition:
        raise AssertionError(f"{name}: {details}")


def max_coulomb_ratio(rows: list[dict[str, str]]) -> float:
    return max((f(row, "max_tangential_force_ratio") for row in rows), default=0.0)


def check_acceptance(args: argparse.Namespace) -> list[dict[str, str]]:
    project_root: Path = args.project_root.resolve()

    if not args.check_only:
        feature_demo = args.feature_demo or find_demo(project_root, "demo_CH_field_contact_feature_switching")
        multislip_demo = args.multislip_demo or find_demo(project_root, "demo_CH_field_contact_multislip_interlock")
        native_demo = args.native_demo or find_demo(project_root, "demo_CH_chrono_native_contact_baseline")

        run_command([str(feature_demo)], project_root)
        run_command([str(multislip_demo)], project_root)
        run_command([str(native_demo)], project_root)
        run_command(
            [
                sys.executable,
                str(args.compare_script),
                "--native",
                str(project_root / "out" / "milestone_23" / "chrono_native_baseline_frames.csv"),
                "--field",
                str(project_root / "out" / "milestone_21" / "field_contact_multislip_interlock_frames.csv"),
                "--output",
                str(project_root / "out" / "milestone_23"),
            ],
            project_root,
        )

    m20 = read_csv(project_root / "out" / "milestone_20" / "field_contact_feature_switching_summary.csv")
    m21_summary = read_csv(project_root / "out" / "milestone_21" / "field_contact_multislip_interlock_summary.csv")
    m21_compare = read_csv(project_root / "out" / "milestone_21" / "field_contact_multislip_interlock_comparison.csv")
    native_compare = read_csv(project_root / "out" / "milestone_23" / "chrono_native_field_comparison.csv")

    results: list[dict[str, str]] = []
    tol = args.tolerance

    direct = find_row(m20, scenario="narrow_groove_entrance", variant="direct_inherit")
    direct_ratio = f(direct, "max_inherited_energy_ratio")
    require(
        direct_ratio > args.direct_inherit_energy_threshold,
        "direct_inherit_energy_amplifies",
        f"narrow_groove_entrance direct_inherit max_inherited_energy_ratio={direct_ratio:.8f} "
        f"> {args.direct_inherit_energy_threshold:.8f}",
        results,
    )

    minimal = find_row(m20, scenario="narrow_groove_entrance", variant="minimal_rotation_gate_split_merge")
    minimal_ratio = f(minimal, "max_inherited_energy_ratio")
    require(
        minimal_ratio <= 1.0 + tol,
        "minimal_rotation_gate_split_merge_non_amplifying",
        f"narrow_groove_entrance minimal_rotation_gate_split_merge max_inherited_energy_ratio={minimal_ratio:.8f} "
        f"<= {1.0 + tol:.8f}",
        results,
    )

    normal_bin_ratios = [
        f(row, "rms_normal_jump_ratio")
        for row in m21_compare
        if row.get("baseline_variant") == "normal_bin_fragments"
    ]
    min_normal_bin_ratio = min(normal_bin_ratios) if normal_bin_ratios else 0.0
    require(
        normal_bin_ratios and min_normal_bin_ratio > args.normal_bin_normal_jump_threshold,
        "normal_bin_fragments_normal_jump_regression",
        f"min rms_normal_jump_ratio over normal_bin_fragments={min_normal_bin_ratio:.8f} "
        f"> {args.normal_bin_normal_jump_threshold:.8f}",
        results,
    )

    native_oscillation_candidates = []
    for row in native_compare:
        native_oscillation_candidates.append(("force", row["scenario"], f(row, "native_to_field_force_oscillation_ratio")))
        native_oscillation_candidates.append(("torque", row["scenario"], f(row, "native_to_field_torque_oscillation_ratio")))
    best_kind, best_scenario, best_native_osc = max(native_oscillation_candidates, key=lambda item: item[2])
    require(
        best_native_osc > args.native_oscillation_threshold,
        "chrono_native_force_or_torque_oscillation_regression",
        f"max native_to_field force/torque oscillation ratio={best_native_osc:.8f} "
        f"({best_kind}, {best_scenario}) > {args.native_oscillation_threshold:.8f}",
        results,
    )

    m20_coulomb = max_coulomb_ratio(m20)
    require(
        m20_coulomb <= 1.0 + tol,
        "milestone_20_all_variants_coulomb_bound",
        f"max_tangential_force_ratio={m20_coulomb:.8f} <= {1.0 + tol:.8f}",
        results,
    )

    m21_coulomb = max_coulomb_ratio(m21_summary)
    require(
        m21_coulomb <= 1.0 + tol,
        "milestone_21_all_variants_coulomb_bound",
        f"max_tangential_force_ratio={m21_coulomb:.8f} <= {1.0 + tol:.8f}",
        results,
    )

    m24_summary = project_root / "out" / "milestone_24" / "mesh_sdf_nonconvex_summary.csv"
    if m24_summary.exists():
        m24_coulomb = max_coulomb_ratio(read_csv(m24_summary))
        require(
            m24_coulomb <= 1.0 + tol,
            "milestone_24_all_variants_coulomb_bound",
            f"max_tangential_force_ratio={m24_coulomb:.8f} <= {1.0 + tol:.8f}",
            results,
        )

    return results


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--feature-demo", type=Path)
    parser.add_argument("--multislip-demo", type=Path)
    parser.add_argument("--native-demo", type=Path)
    parser.add_argument("--compare-script", type=Path, default=default_root / "scripts" / "compare_chrono_native_field.py")
    parser.add_argument("--check-only", action="store_true", help="do not run demos; only validate existing CSV outputs")
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument("--direct-inherit-energy-threshold", type=float, default=1.5)
    parser.add_argument("--normal-bin-normal-jump-threshold", type=float, default=1.5)
    parser.add_argument("--native-oscillation-threshold", type=float, default=1.03)
    parser.add_argument("--summary-output", type=Path, default=default_root / "out" / "milestone_26" / "acceptance_regression_summary.json")
    args = parser.parse_args()

    try:
        results = check_acceptance(args)
    except (AssertionError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ACCEPTANCE REGRESSION FAILED: {exc}", file=sys.stderr)
        return 1

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.summary_output}")
    print("ACCEPTANCE REGRESSION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
