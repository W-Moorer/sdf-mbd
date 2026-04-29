#!/usr/bin/env python3
"""Smooth Fischer-Burmeister epsilon sensitivity for point-mass SDF-NCP contact."""

from __future__ import annotations

import csv
import math
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
from sdf_mbd.solvers import simulate_ncp_point_mass


MASS = 1.0
GRAVITY = 9.81
Q0 = [0.0, 1.0]
V0 = [0.0, 0.0]
DT = 1.0e-3
T_END = 1.0
EPS_LIST = [1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7]

RESULT_DIR = ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity"
SUMMARY_PATH = RESULT_DIR / "summary.csv"

FIELDNAMES = [
    "eps",
    "max_penetration",
    "mean_penetration",
    "max_complementarity_error",
    "mean_complementarity_error",
    "max_ncp_residual",
    "mean_ncp_residual",
    "mean_iterations",
    "success_rate",
    "max_lambda",
    "runtime_seconds",
]


def summarize_eps(eps: float) -> dict[str, float]:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    start = time.perf_counter()
    try:
        trajectory = simulate_ncp_point_mass(ground, Q0, V0, MASS, GRAVITY, DT, T_END, eps)
    except Exception:
        return {
            "eps": eps,
            "max_penetration": math.nan,
            "mean_penetration": math.nan,
            "max_complementarity_error": math.nan,
            "mean_complementarity_error": math.nan,
            "max_ncp_residual": math.nan,
            "mean_ncp_residual": math.nan,
            "mean_iterations": math.nan,
            "success_rate": 0.0,
            "max_lambda": math.nan,
            "runtime_seconds": time.perf_counter() - start,
        }

    diagnostics = [diag for _, _, diag in trajectory]
    count = max(1, len(diagnostics))
    return {
        "eps": eps,
        "max_penetration": max(diag.penetration for diag in diagnostics),
        "mean_penetration": sum(diag.penetration for diag in diagnostics) / count,
        "max_complementarity_error": max(diag.complementarity_error for diag in diagnostics),
        "mean_complementarity_error": sum(diag.complementarity_error for diag in diagnostics) / count,
        "max_ncp_residual": max(diag.ncp_residual for diag in diagnostics),
        "mean_ncp_residual": sum(diag.ncp_residual for diag in diagnostics) / count,
        "mean_iterations": sum(diag.solver_iterations for diag in diagnostics) / count,
        "success_rate": sum(1.0 for diag in diagnostics if diag.solver_success) / count,
        "max_lambda": max(abs(diag.lambda_n) for diag in diagnostics),
        "runtime_seconds": time.perf_counter() - start,
    }


def write_summary(rows: list[dict[str, float]]) -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with SUMMARY_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def finite_series(rows: list[dict[str, float]], key: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for row in sorted(rows, key=lambda r: r["eps"]):
        value = row[key]
        if math.isfinite(value):
            xs.append(row["eps"])
            ys.append(value)
    return xs, ys


def plot_single(rows: list[dict[str, float]], key: str, ylabel: str, filename: str) -> None:
    xs, ys = finite_series(rows, key)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(xs, ys, marker="o", linewidth=1.8)
    ax.set_xscale("log")
    ax.set_xlabel("smooth FB eps")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / filename, dpi=180)
    plt.close(fig)


def plot_complementarity(rows: list[dict[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for key, label, marker in (
        ("max_complementarity_error", "max complementarity error", "o"),
        ("mean_complementarity_error", "mean complementarity error", "s"),
    ):
        xs, ys = finite_series(rows, key)
        ax.plot(xs, ys, marker=marker, linewidth=1.8, label=label)
    ax.set_xscale("log")
    ax.set_xlabel("smooth FB eps")
    ax.set_ylabel("complementarity error")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "complementarity_vs_eps.png", dpi=180)
    plt.close(fig)


def make_plots(rows: list[dict[str, float]]) -> None:
    plot_single(rows, "max_penetration", "max penetration (m)", "max_penetration_vs_eps.png")
    plot_complementarity(rows)
    plot_single(rows, "mean_iterations", "mean nonlinear solver iterations", "iterations_vs_eps.png")


def main() -> None:
    rows = [summarize_eps(eps) for eps in EPS_LIST]
    write_summary(rows)
    make_plots(rows)
    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote plots under {RESULT_DIR}")


if __name__ == "__main__":
    main()
