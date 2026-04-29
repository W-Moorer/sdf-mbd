#!/usr/bin/env python3
"""Plot SDF-NCP field-contact benchmark CSV outputs."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "sdf_ncp_benchmarks"
FIGURE_DIR = RESULT_DIR / "figures"
CASES = [
    "headon_spheres",
    "headon_spheres_mass_ratio",
    "cam",
    "cam_recurdyn_solid_contact",
    "eccentric_roller",
    "onset_stress",
    "simple_gear",
]


def read_rows(case_name: str) -> list[dict[str, str]]:
    path = RESULT_DIR / case_name / "trajectory.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run python scripts/run_sdf_ncp_field_contact_benchmarks.py first."
        )
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def unique_contact_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[tuple[str, str]] = set()
    out: list[dict[str, str]] = []
    for row in rows:
        key = (row["time"], row["contact_id"])
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def by_body(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["body_id"], []).append(row)
    for body_rows in grouped.values():
        body_rows.sort(key=lambda r: f(r, "time"))
    return grouped


def angle_about_z(row: dict[str, str]) -> float:
    return 2.0 * math.atan2(f(row, "q3"), f(row, "q0"))


def angle_about_x(row: dict[str, str]) -> float:
    return 2.0 * math.atan2(f(row, "q1"), f(row, "q0"))


def plot_contact_quantity(case_name: str, rows: list[dict[str, str]], key: str, ylabel: str, filename: str) -> Path:
    path = FIGURE_DIR / case_name / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    contact_rows = unique_contact_rows(rows)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot([f(r, "time") for r in contact_rows], [f(r, key) for r in contact_rows], linewidth=1.6)
    ax.set_xlabel("time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(case_name)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_total_contact_force(case_name: str, rows: list[dict[str, str]], filename: str) -> Path:
    path = FIGURE_DIR / case_name / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    force_key = "lambda_force" if rows and "lambda_force" in rows[0] else "lambda_n"
    by_time: dict[float, float] = {}
    for row in rows:
        if row.get("contact_id", "-1") == "-1":
            continue
        time = f(row, "time")
        by_time[time] = by_time.get(time, 0.0) + f(row, force_key)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    times = sorted(by_time)
    ax.plot(times, [by_time[t] for t in times], linewidth=1.6)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("total normal force contribution")
    ax.set_title(case_name)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_headon_extra(case_name: str, rows: list[dict[str, str]]) -> list[Path]:
    case_dir = FIGURE_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    grouped = by_body(rows)

    outputs: list[Path] = []
    if "sphere_a" in grouped and "sphere_b" in grouped:
        a_rows = grouped["sphere_a"]
        b_rows = grouped["sphere_b"]
        path = case_dir / "relative_distance_vs_time.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(
            [f(r, "time") for r in a_rows],
            [abs(f(b, "px") - f(a, "px")) for a, b in zip(a_rows, b_rows)],
            linewidth=1.6,
        )
        ax.set_xlabel("time (s)")
        ax.set_ylabel("relative center distance (m)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

        path = case_dir / "velocity_vs_time.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot([f(r, "time") for r in a_rows], [f(r, "vx") for r in a_rows], label="sphere_a vx", linewidth=1.6)
        ax.plot([f(r, "time") for r in b_rows], [f(r, "vx") for r in b_rows], label="sphere_b vx", linewidth=1.6)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("x velocity (m/s)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)
    return outputs


def plot_cam_extra(case_name: str, rows: list[dict[str, str]]) -> list[Path]:
    case_dir = FIGURE_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    grouped = by_body(rows)
    outputs: list[Path] = []

    if "cam" in grouped:
        cam_rows = grouped["cam"]
        path = case_dir / "cam_angle_or_position_vs_time.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot([f(r, "time") for r in cam_rows], [angle_about_z(r) for r in cam_rows], linewidth=1.6)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("cam angle about z (rad)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    if "follower" in grouped:
        follower_rows = grouped["follower"]
        path = case_dir / "follower_response_vs_time.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot([f(r, "time") for r in follower_rows], [f(r, "py") for r in follower_rows], linewidth=1.6)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("follower y (m)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    outputs.append(plot_total_contact_force(case_name, rows, "contact_force_vs_time.png"))
    return outputs


def plot_case(case_name: str) -> list[Path]:
    rows = read_rows(case_name)
    outputs = [
        plot_contact_quantity(case_name, rows, "gap", "SDF gap (m)", "gap_vs_time.png"),
        plot_contact_quantity(case_name, rows, "lambda_n", "normal multiplier", "lambda_vs_time.png"),
        plot_contact_quantity(case_name, rows, "penetration", "penetration (m)", "penetration_vs_time.png"),
        plot_contact_quantity(
            case_name,
            rows,
            "complementarity_error",
            "complementarity error",
            "complementarity_vs_time.png",
        ),
    ]
    if case_name.startswith("headon_spheres"):
        outputs.extend(plot_headon_extra(case_name, rows))
    if case_name.startswith("cam"):
        outputs.extend(plot_cam_extra(case_name, rows))
    if case_name in {"eccentric_roller", "onset_stress"}:
        outputs.extend(plot_cam_extra(case_name, rows))
    if case_name == "simple_gear":
        grouped = by_body(rows)
        case_dir = FIGURE_DIR / case_name
        case_dir.mkdir(parents=True, exist_ok=True)
        gear_rows = grouped.get("gear22", [])
        outputs.append(plot_total_contact_force(case_name, rows, "contact_force_vs_time.png"))
        if gear_rows:
            path = case_dir / "gear_angle_vs_time.png"
            fig, ax = plt.subplots(figsize=(8.0, 4.8))
            ax.plot([f(r, "time") for r in gear_rows], [angle_about_x(r) for r in gear_rows], linewidth=1.6)
            ax.set_xlabel("time (s)")
            ax.set_ylabel("gear22 angle about x (rad)")
            ax.set_title(case_name)
            ax.grid(True, alpha=0.25)
            fig.tight_layout()
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)

            path = case_dir / "angular_velocity_vs_time.png"
            fig, ax = plt.subplots(figsize=(8.0, 4.8))
            ax.plot([f(r, "time") for r in gear_rows], [f(r, "wx") for r in gear_rows], linewidth=1.6)
            ax.set_xlabel("time (s)")
            ax.set_ylabel("gear22 omega_x (rad/s)")
            ax.set_title(case_name)
            ax.grid(True, alpha=0.25)
            fig.tight_layout()
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)
    return outputs


def main() -> None:
    try:
        outputs: list[Path] = []
        for case_name in CASES:
            outputs.extend(plot_case(case_name))
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)

    for path in outputs:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
