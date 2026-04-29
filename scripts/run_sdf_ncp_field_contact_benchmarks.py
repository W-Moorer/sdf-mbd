#!/usr/bin/env python3
"""Run and summarize SDF-NCP benchmarks based on existing field_contact assets."""

from __future__ import annotations

import csv
import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "build" / "bin" / "Release" / "demo_CH_sdf_ncp_benchmarks.exe"
DEMO_OPENVDB = ROOT / "build" / "bin" / "Release" / "demo_CH_sdf_ncp_benchmarks_openvdb.exe"
RESULT_DIR = ROOT / "results" / "sdf_ncp_benchmarks"
SUMMARY_PATH = RESULT_DIR / "benchmark_summary.csv"
COMPARISON_SUMMARY_PATH = RESULT_DIR / "benchmark_comparison_summary.csv"

DEFAULT_CASES = [
    "headon_spheres",
    "headon_spheres_mass_ratio",
    "cam",
    "eccentric_roller",
    "onset_stress",
    "simple_gear",
    "rev_joint_clearance",
]
REFERENCE_CASES = [
    "cam_recurdyn_solid_contact",
    "rev_joint_clearance_ggeomcontact_calibration",
    "rev_joint_clearance_ggeomcontact_hht",
    "rev_joint_clearance_ggeomcontact_euler_substep",
    "rev_joint_clearance_ggeomcontact_hht_substep",
]
MANUAL_CASES = [
    "rev_joint_clearance_hht_substep",
]
SIMPLE_GEAR_DT_SWEEP_CASES = [
    "simple_gear_dt_001",
    "simple_gear_dt_0005",
    "simple_gear_dt_0001",
]
OPENVDB_CASES = {
    "cam",
    "eccentric_roller",
    "onset_stress",
    "simple_gear",
    "rev_joint_clearance",
    *REFERENCE_CASES,
    *MANUAL_CASES,
    *SIMPLE_GEAR_DT_SWEEP_CASES,
}


def run_case(case_name: str) -> None:
    demo = DEMO_OPENVDB if case_name in OPENVDB_CASES else DEMO
    if not demo.exists():
        target = "demo_CH_sdf_ncp_benchmarks_openvdb" if case_name in OPENVDB_CASES else "demo_CH_sdf_ncp_benchmarks"
        print(
            f"Missing benchmark executable: {demo}\n"
            "Build it first with:\n"
            f"  cmake --build build --config Release --target {target}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    print(f"[sdf-ncp-benchmarks] Running {case_name}", flush=True)
    completed = subprocess.run([str(demo), case_name], cwd=ROOT)
    if completed.returncode != 0:
        print(f"Benchmark failed: {case_name}, return code {completed.returncode}", file=sys.stderr)
        raise SystemExit(completed.returncode)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        rows = []
        for row in csv.DictReader(f):
            extras = row.pop(None, None)
            if extras:
                if "notes" in row:
                    row["notes"] = ",".join([row.get("notes", ""), *extras]).strip(",")
                else:
                    row["extra"] = ",".join(extras)
            rows.append(row)
        return rows


def write_merged(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge_outputs(cases: list[str]) -> None:
    summary_rows: list[dict[str, str]] = []
    comparison_rows: list[dict[str, str]] = []
    for case_name in cases:
        case_dir = RESULT_DIR / case_name
        summary_rows.extend(read_csv(case_dir / "summary.csv"))
        comparison_rows.extend(read_csv(case_dir / "comparison.csv"))
    write_merged(SUMMARY_PATH, summary_rows)
    write_merged(COMPARISON_SUMMARY_PATH, comparison_rows)
    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote {COMPARISON_SUMMARY_PATH}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        nargs="*",
        default=None,
        help="Benchmark cases to run. Defaults to the stable SDF-NCP benchmark set.",
    )
    parser.add_argument(
        "--include-manual",
        action="store_true",
        help="Include manually classified cases when such cases exist. Currently all implemented cases are stable.",
    )
    parser.add_argument(
        "--include-reference",
        action="store_true",
        help="Include long-form RecurDyn reference-alignment benchmarks such as cam_recurdyn_solid_contact.",
    )
    parser.add_argument(
        "--include-simple-gear-dt-sweep",
        action="store_true",
        help="Include simple_gear dt stability cases: 0.001, 0.0005, and 0.0001 seconds.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cases = list(args.cases) if args.cases else list(DEFAULT_CASES)
    if args.include_manual:
        for case_name in MANUAL_CASES:
            if case_name not in cases:
                cases.append(case_name)
    if args.include_reference:
        for case_name in REFERENCE_CASES:
            if case_name not in cases:
                cases.append(case_name)
    if args.include_simple_gear_dt_sweep:
        for case_name in SIMPLE_GEAR_DT_SWEEP_CASES:
            if case_name not in cases:
                cases.append(case_name)

    for case_name in cases:
        run_case(case_name)
    merge_outputs(cases)


if __name__ == "__main__":
    main()
