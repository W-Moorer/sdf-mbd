#!/usr/bin/env python3
"""Plot the C++ planar rigid-body SDF-NCP rollout CSV."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = ROOT / "results" / "sdf_ncp_cpp" / "rigidbody2d_rollout.csv"
FIGURE_DIR = ROOT / "results" / "sdf_ncp_cpp" / "figures"


def read_rows() -> list[dict[str, float]]:
    if not CSV_PATH.exists():
        raise FileNotFoundError(
            f"{CSV_PATH} does not exist. Run "
            r"build\bin\Release\demo_CH_sdf_ncp_regression.exe rigidbody2d_export first."
        )
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [{key: float(value) for key, value in row.items()} for row in reader]


def series(rows: list[dict[str, float]], key: str) -> list[float]:
    return [row[key] for row in rows]


def plot_pose(rows: list[dict[str, float]]) -> Path:
    path = FIGURE_DIR / "rigidbody2d_pose_vs_time.png"
    time = series(rows, "time")
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 7.0), sharex=True)
    for ax, key, ylabel in zip(axes, ("x", "y", "theta"), ("x (m)", "y (m)", "theta (rad)")):
        ax.plot(time, series(rows, key), linewidth=1.6)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("time (s)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_gaps(rows: list[dict[str, float]]) -> Path:
    path = FIGURE_DIR / "rigidbody2d_gaps_vs_time.png"
    time = series(rows, "time")
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(time, series(rows, "gap0"), label="gap0", linewidth=1.6)
    ax.plot(time, series(rows, "gap1"), label="gap1", linewidth=1.6)
    ax.axhline(0.0, color="#111827", linewidth=0.9)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("SDF gap (m)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_lambdas(rows: list[dict[str, float]]) -> Path:
    path = FIGURE_DIR / "rigidbody2d_lambdas_vs_time.png"
    time = series(rows, "time")
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(time, series(rows, "lambda0"), label="lambda0", linewidth=1.6)
    ax.plot(time, series(rows, "lambda1"), label="lambda1", linewidth=1.6)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("normal multiplier (N)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_complementarity(rows: list[dict[str, float]]) -> Path:
    path = FIGURE_DIR / "rigidbody2d_complementarity_vs_time.png"
    time = series(rows, "time")
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(time, series(rows, "max_complementarity_error"), label="max complementarity error", linewidth=1.6)
    ax.plot(time, series(rows, "ncp_residual_norm"), label="NCP residual norm", linewidth=1.6)
    ax.plot(time, series(rows, "residual_norm"), label="total residual norm", linewidth=1.2, alpha=0.8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("residual / error")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    try:
        rows = read_rows()
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        plot_pose(rows),
        plot_gaps(rows),
        plot_lambdas(rows),
        plot_complementarity(rows),
    ]
    for path in outputs:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
