#!/usr/bin/env python3
"""Run the simple_gear SDF-NCP time-step stability sweep and check acceptance metrics."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "bin" / "Release" / "demo_CH_sdf_ncp_benchmarks_openvdb.exe"
RESULT_ROOT = ROOT / "results" / "sdf_ncp_benchmarks"
SWEEP_CASES = [
    "simple_gear_dt_001",
    "simple_gear_dt_0005",
    "simple_gear_dt_0001",
]
SWEEP_SUMMARY = RESULT_ROOT / "simple_gear_dt_sweep_summary.csv"


def read_one_csv(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"CSV is empty: {path}")
    return rows[0]


def read_float(row: dict[str, str], key: str, default: float = math.nan) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def compute_mae_from_timeseries(path: Path) -> float:
    if not path.exists():
        return math.nan
    errors: list[float] = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("omega_abs_error", "") != "":
                errors.append(float(row["omega_abs_error"]))
    return sum(errors) / len(errors) if errors else math.nan


def run_case(case_name: str) -> None:
    if not EXE.exists():
        print(
            f"Missing executable: {EXE}\n"
            "Build first with:\n"
            "  cmake --build build --config Release --target demo_CH_sdf_ncp_benchmarks_openvdb",
            file=sys.stderr,
        )
        raise SystemExit(1)
    print(f"[simple-gear-dt-sweep] running {case_name}", flush=True)
    completed = subprocess.run([str(EXE), case_name], cwd=ROOT)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def collect_case(case_name: str) -> dict[str, str]:
    case_dir = RESULT_ROOT / case_name
    summary = read_one_csv(case_dir / "summary.csv")
    analytic = read_one_csv(case_dir / "gear22_analytic_comparison_summary.csv")
    mae = read_float(analytic, "mae")
    if not math.isfinite(mae):
        mae = compute_mae_from_timeseries(case_dir / "gear22_analytic_comparison_timeseries.csv")

    return {
        "case_name": case_name,
        "dt": summary["dt"],
        "num_steps": summary["num_steps"],
        "success_rate": summary["success_rate"],
        "max_penetration": summary["max_penetration"],
        "mean_iterations": summary["mean_iterations"],
        "rmse": analytic["rmse"],
        "mae": f"{mae:.17g}",
        "max_abs_error": analytic["max_abs_error"],
        "final_value": analytic["final_value"],
        "final_abs_error": analytic["final_abs_error"],
    }


def write_rows(rows: list[dict[str, str]]) -> None:
    SWEEP_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SWEEP_SUMMARY.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {SWEEP_SUMMARY}")


def validate(rows: list[dict[str, str]], max_penetration: float, rmse_factor: float, mae_factor: float) -> None:
    baseline_path = RESULT_ROOT / "simple_gear" / "gear22_analytic_comparison_summary.csv"
    if baseline_path.exists():
        baseline = read_one_csv(baseline_path)
        baseline_rmse = read_float(baseline, "rmse")
        baseline_mae = read_float(baseline, "mae")
        if not math.isfinite(baseline_mae):
            baseline_mae = compute_mae_from_timeseries(
                RESULT_ROOT / "simple_gear" / "gear22_analytic_comparison_timeseries.csv"
            )
    else:
        baseline_rmse = read_float(rows[0], "rmse")
        baseline_mae = read_float(rows[0], "mae")

    failures: list[str] = []
    for row in rows:
        case = row["case_name"]
        success_rate = read_float(row, "success_rate")
        penetration = read_float(row, "max_penetration")
        rmse = read_float(row, "rmse")
        mae = read_float(row, "mae")

        if success_rate != 1.0:
            failures.append(f"{case}: success_rate={success_rate}, expected 1.0")
        if not math.isfinite(penetration) or penetration > max_penetration:
            failures.append(f"{case}: max_penetration={penetration}, limit {max_penetration}")
        if not math.isfinite(rmse) or not math.isfinite(mae):
            failures.append(f"{case}: non-finite MAE/RMSE")

        if case == "simple_gear_dt_001":
            if rmse > baseline_rmse * (1.0 + 1.0e-10):
                failures.append(f"{case}: RMSE {rmse} worse than baseline {baseline_rmse}")
            if math.isfinite(baseline_mae) and mae > baseline_mae * (1.0 + 1.0e-10):
                failures.append(f"{case}: MAE {mae} worse than baseline {baseline_mae}")
        else:
            if rmse > max(rmse_factor * baseline_rmse, 0.25):
                failures.append(f"{case}: RMSE {rmse} exceeds non-explosion limit")
            if math.isfinite(baseline_mae) and mae > max(mae_factor * baseline_mae, 0.15):
                failures.append(f"{case}: MAE {mae} exceeds non-explosion limit")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        raise SystemExit(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-run", action="store_true", help="Only read existing CSV files.")
    parser.add_argument("--max-penetration", type=float, default=1.0e-5)
    parser.add_argument("--rmse-factor", type=float, default=2.0)
    parser.add_argument("--mae-factor", type=float, default=2.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.no_run:
        for case_name in SWEEP_CASES:
            run_case(case_name)

    rows = [collect_case(case_name) for case_name in SWEEP_CASES]
    write_rows(rows)
    validate(rows, args.max_penetration, args.rmse_factor, args.mae_factor)
    print("simple_gear SDF-NCP dt sweep acceptance passed.")


if __name__ == "__main__":
    main()
