#!/usr/bin/env python3
"""Build a compact CSV summary table for the first SDF-NCP paper figures."""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PENALTY_SUMMARY = ROOT / "results" / "sdf_ncp" / "penalty_sensitivity" / "summary.csv"
EPS_SUMMARY = ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity" / "summary.csv"
RIGIDBODY_CSV = ROOT / "results" / "sdf_ncp_cpp" / "rigidbody2d_rollout.csv"
OUT_DIR = ROOT / "results" / "sdf_ncp_paper1" / "tables"
OUT_PATH = OUT_DIR / "method_summary.csv"

FIELDNAMES = [
    "experiment",
    "method",
    "parameter",
    "max_penetration",
    "max_complementarity_error",
    "mean_iterations",
    "success_rate",
]


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input is missing: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def add_penalty_rows(rows: list[dict[str, str]], mean_iterations_by_eps: dict[float, str]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for row in rows:
        mean_iterations = "0.0"
        if row["method"] == "sdf_ncp":
            mean_iterations = mean_iterations_by_eps.get(float(row["parameter"]), "")
        output.append(
            {
                "experiment": "penalty_sensitivity",
                "method": row["method"],
                "parameter": row["parameter"],
                "max_penetration": row["max_penetration"],
                "max_complementarity_error": row["max_complementarity_error"],
                "mean_iterations": mean_iterations,
                "success_rate": row["success_rate"],
            }
        )
    return output


def add_epsilon_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for row in rows:
        output.append(
            {
                "experiment": "epsilon_sensitivity",
                "method": "sdf_ncp",
                "parameter": row["eps"],
                "max_penetration": row["max_penetration"],
                "max_complementarity_error": row["max_complementarity_error"],
                "mean_iterations": row["mean_iterations"],
                "success_rate": row["success_rate"],
            }
        )
    return output


def add_rigidbody_row(rows: list[dict[str, str]]) -> dict[str, str]:
    if not rows:
        raise ValueError(f"No rows found in {RIGIDBODY_CSV}")
    max_penetration = max(float(row["max_penetration"]) for row in rows)
    max_complementarity_error = max(float(row["max_complementarity_error"]) for row in rows)
    mean_iterations = statistics.fmean(float(row["iterations"]) for row in rows)
    success_rate = statistics.fmean(float(row["success"]) for row in rows)
    return {
        "experiment": "rigidbody2d_rollout",
        "method": "sdf_ncp_cpp_rigidbody2d",
        "parameter": "eps=1e-6",
        "max_penetration": f"{max_penetration:.17g}",
        "max_complementarity_error": f"{max_complementarity_error:.17g}",
        "mean_iterations": f"{mean_iterations:.17g}",
        "success_rate": f"{success_rate:.17g}",
    }


def build_summary() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    eps_rows = read_csv(EPS_SUMMARY)
    mean_iterations_by_eps = {float(row["eps"]): row["mean_iterations"] for row in eps_rows}
    rows.extend(add_penalty_rows(read_csv(PENALTY_SUMMARY), mean_iterations_by_eps))
    rows.extend(add_epsilon_rows(eps_rows))
    rows.append(add_rigidbody_row(read_csv(RIGIDBODY_CSV)))
    return rows


def main() -> None:
    try:
        rows = build_summary()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
