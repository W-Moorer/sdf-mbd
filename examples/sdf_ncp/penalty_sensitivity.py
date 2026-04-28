#!/usr/bin/env python3
"""Penalty stiffness sensitivity compared with SDF-NCP smoothing choices."""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sdf_mbd.sdf import PlaneSDF
from sdf_mbd.solvers import simulate_ncp_point_mass, simulate_penalty_point_mass


MASS = 1.0
GRAVITY = 9.81
Q0 = [0.0, 1.0]
V0 = [0.0, 0.0]
DT = 1.0e-3
T_END = 1.0

RESULT_DIR = ROOT / "results" / "sdf_ncp" / "penalty_sensitivity"
SUMMARY_PATH = RESULT_DIR / "summary.csv"

FIELDNAMES = [
    "method",
    "parameter",
    "max_penetration",
    "mean_penetration",
    "max_force",
    "mean_complementarity_error",
    "max_complementarity_error",
    "mean_residual_norm",
    "success_rate",
    "runtime_seconds",
]


def summarize(method: str, parameter: float, trajectory, runtime_seconds: float) -> dict[str, float | str]:
    diagnostics = [diag for _, _, diag in trajectory]
    count = max(1, len(diagnostics))
    return {
        "method": method,
        "parameter": parameter,
        "max_penetration": max(diag.penetration for diag in diagnostics),
        "mean_penetration": sum(diag.penetration for diag in diagnostics) / count,
        "max_force": max(abs(diag.lambda_n) for diag in diagnostics),
        "mean_complementarity_error": sum(diag.complementarity_error for diag in diagnostics) / count,
        "max_complementarity_error": max(diag.complementarity_error for diag in diagnostics),
        "mean_residual_norm": sum(diag.residual_norm for diag in diagnostics) / count,
        "success_rate": sum(1.0 for diag in diagnostics if diag.solver_success) / count,
        "runtime_seconds": runtime_seconds,
    }


def run_sensitivity() -> list[dict[str, float | str]]:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    rows: list[dict[str, float | str]] = []

    for k_n in (1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6):
        start = time.perf_counter()
        trajectory = simulate_penalty_point_mass(ground, Q0, V0, MASS, GRAVITY, DT, T_END, k_n)
        rows.append(summarize("penalty", k_n, trajectory, time.perf_counter() - start))

    for eps in (1.0e-2, 1.0e-4, 1.0e-6):
        start = time.perf_counter()
        trajectory = simulate_ncp_point_mass(ground, Q0, V0, MASS, GRAVITY, DT, T_END, eps)
        rows.append(summarize("sdf_ncp", eps, trajectory, time.perf_counter() - start))

    return rows


def write_summary(rows: list[dict[str, float | str]]) -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with SUMMARY_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows: list[dict[str, float | str]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for method, marker in (("penalty", "o"), ("sdf_ncp", "s")):
        series = sorted((r for r in rows if r["method"] == method), key=lambda r: float(r["parameter"]))
        ax.plot(
            [float(r["parameter"]) for r in series],
            [float(r["max_penetration"]) for r in series],
            marker=marker,
            linewidth=1.8,
            label=method,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("parameter (k_n for penalty, eps for SDF-NCP)")
    ax.set_ylabel("max penetration (m)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "max_penetration_vs_parameter.png", dpi=180)
    plt.close(fig)


def main() -> None:
    rows = run_sensitivity()
    write_summary(rows)
    make_plot(rows)
    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote {RESULT_DIR / 'max_penetration_vs_parameter.png'}")


if __name__ == "__main__":
    main()

