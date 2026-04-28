#!/usr/bin/env python3
"""Point mass falling onto a ground-plane SDF.

Run from the repository root:

    python examples/sdf_ncp/point_mass_plane.py
"""

from __future__ import annotations

import csv
import sys
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

RESULT_DIR = ROOT / "results" / "sdf_ncp" / "point_mass_plane"
CSV_PATH = RESULT_DIR / "point_mass_plane.csv"

FIELDNAMES = [
    "time",
    "method",
    "parameter",
    "x",
    "y",
    "vx",
    "vy",
    "gap",
    "penetration",
    "lambda_n",
    "complementarity_error",
    "ncp_residual",
    "solver_success",
    "residual_norm",
]


def row_dict(time: float, method: str, parameter: float, state, diagnostics) -> dict[str, float | str | bool]:
    return {
        "time": time,
        "method": method,
        "parameter": parameter,
        "x": state.q[0],
        "y": state.q[1],
        "vx": state.v[0],
        "vy": state.v[1],
        "gap": diagnostics.gap,
        "penetration": diagnostics.penetration,
        "lambda_n": diagnostics.lambda_n,
        "complementarity_error": diagnostics.complementarity_error,
        "ncp_residual": diagnostics.ncp_residual,
        "solver_success": diagnostics.solver_success,
        "residual_norm": diagnostics.residual_norm,
    }


def run_all_methods() -> list[dict[str, float | str | bool]]:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    rows: list[dict[str, float | str | bool]] = []

    for k_n in (1.0e3, 1.0e4, 1.0e5):
        trajectory = simulate_penalty_point_mass(ground, Q0, V0, MASS, GRAVITY, DT, T_END, k_n)
        rows.extend(row_dict(t, "penalty", k_n, state, diag) for t, state, diag in trajectory)

    for eps in (1.0e-4, 1.0e-6):
        trajectory = simulate_ncp_point_mass(ground, Q0, V0, MASS, GRAVITY, DT, T_END, eps)
        rows.extend(row_dict(t, "sdf_ncp", eps, state, diag) for t, state, diag in trajectory)

    return rows


def write_csv(rows: list[dict[str, float | str | bool]]) -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def group_rows(rows: list[dict[str, float | str | bool]]) -> dict[str, list[dict[str, float | str | bool]]]:
    grouped: dict[str, list[dict[str, float | str | bool]]] = {}
    for row in rows:
        label = f"{row['method']} {float(row['parameter']):.0e}"
        grouped.setdefault(label, []).append(row)
    return grouped


def plot_quantity(
    rows: list[dict[str, float | str | bool]],
    key: str,
    ylabel: str,
    filename: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for label, series in group_rows(rows).items():
        series = sorted(series, key=lambda r: float(r["time"]))
        ax.plot([float(r["time"]) for r in series], [float(r[key]) for r in series], label=label, linewidth=1.6)
    ax.set_xlabel("time (s)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / filename, dpi=180)
    plt.close(fig)


def make_plots(rows: list[dict[str, float | str | bool]]) -> None:
    plot_quantity(rows, "y", "height y (m)", "height_vs_time.png")
    plot_quantity(rows, "gap", "gap phi(q) (m)", "gap_vs_time.png")
    plot_quantity(rows, "penetration", "penetration max(0, -g) (m)", "penetration_vs_time.png")
    plot_quantity(rows, "lambda_n", "normal contact force lambda_n (N)", "contact_force_vs_time.png")
    plot_quantity(rows, "complementarity_error", "complementarity error", "complementarity_error_vs_time.png")


def main() -> None:
    rows = run_all_methods()
    write_csv(rows)
    make_plots(rows)
    print(f"Wrote {CSV_PATH}")
    print(f"Wrote plots under {RESULT_DIR}")


if __name__ == "__main__":
    main()

