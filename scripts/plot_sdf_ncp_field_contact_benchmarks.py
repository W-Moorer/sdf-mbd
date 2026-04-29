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
    "rev_joint_clearance",
    "rev_joint_clearance_ggeomcontact_calibration",
    "rev_joint_clearance_ggeomcontact_hht",
    "rev_joint_clearance_ggeomcontact_euler_substep",
    "rev_joint_clearance_ggeomcontact_hht_substep",
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


def read_optional_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        print(f"Skipping missing optional CSV: {path}", file=sys.stderr)
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def interpolate_series(times: list[float], source_times: list[float], source_values: list[float]) -> list[float]:
    if not source_times:
        return [math.nan for _ in times]
    values: list[float] = []
    j = 0
    last = len(source_times) - 1
    for t in times:
        if t <= source_times[0]:
            values.append(source_values[0])
            continue
        if t >= source_times[last]:
            values.append(source_values[last])
            continue
        while j + 1 < last and source_times[j + 1] < t:
            j += 1
        t0 = source_times[j]
        t1 = source_times[j + 1]
        v0 = source_values[j]
        v1 = source_values[j + 1]
        alpha = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
        values.append(v0 + alpha * (v1 - v0))
    return values


def plot_rev_joint_clearance_extra(case_name: str, rows: list[dict[str, str]]) -> list[Path]:
    case_dir = FIGURE_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = [plot_total_contact_force(case_name, rows, "contact_force_vs_time.png")]

    reference_rows = read_optional_csv(RESULT_DIR / case_name / "recurdyn_reference_comparison_timeseries.csv")
    if not reference_rows:
        return outputs
    ideal_rows = read_optional_csv(ROOT / "assets" / "rev_joint_clearance" / "data" / "body2_ideal.csv")

    def rf(row: dict[str, str], key: str) -> float:
        return float(row[key])

    times = [rf(r, "time") for r in reference_rows]
    ideal_time_key = "X:Pos_TX-Body2-rev_clearance_joint_ideal(m)"
    ideal_y_key = "Y:Pos_TY-Body2-rev_clearance_joint_ideal(m)"
    ideal_z_key = "Y:Pos_TZ-Body2-rev_clearance_joint_ideal(m)"
    ideal_times: list[float] = []
    ideal_y: list[float] = []
    ideal_z: list[float] = []
    if ideal_rows and all(key in ideal_rows[0] for key in (ideal_time_key, ideal_y_key, ideal_z_key)):
        ideal_times = [rf(r, ideal_time_key) for r in ideal_rows]
        ideal_y = [rf(r, ideal_y_key) for r in ideal_rows]
        ideal_z = [rf(r, ideal_z_key) for r in ideal_rows]

    result_label = "GGEOM calibration" if "ggeomcontact" in case_name else "SDF-NCP"

    path = case_dir / "body2_position_vs_recurdyn.png"
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
    axes[0].plot(times, [rf(r, "body2_y_ref") for r in reference_rows], label="RecurDyn clearance", linewidth=1.4)
    axes[0].plot(times, [rf(r, "body2_y_sdf_ncp") for r in reference_rows], label=result_label, linewidth=1.4)
    axes[1].plot(times, [rf(r, "body2_z_ref") for r in reference_rows], label="RecurDyn clearance", linewidth=1.4)
    axes[1].plot(times, [rf(r, "body2_z_sdf_ncp") for r in reference_rows], label=result_label, linewidth=1.4)
    if ideal_times:
        axes[0].plot(ideal_times, ideal_y, label="Ideal hinge", linewidth=1.2, linestyle="--")
        axes[1].plot(ideal_times, ideal_z, label="Ideal hinge", linewidth=1.2, linestyle="--")
    axes[0].set_ylabel("body2 y (m)")
    axes[1].set_ylabel("body2 z (m)")
    axes[1].set_xlabel("time (s)")
    axes[0].set_title(case_name)
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    if ideal_times:
        ideal_y_on_ref = interpolate_series(times, ideal_times, ideal_y)
        ideal_z_on_ref = interpolate_series(times, ideal_times, ideal_z)
        path = case_dir / "body2_error_vs_ideal_hinge.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(
            times,
            [rf(r, "body2_y_sdf_ncp") - y for r, y in zip(reference_rows, ideal_y_on_ref)],
            label=f"{result_label} body2 y - ideal",
            linewidth=1.4,
        )
        ax.plot(
            times,
            [rf(r, "body2_z_sdf_ncp") - z for r, z in zip(reference_rows, ideal_z_on_ref)],
            label=f"{result_label} body2 z - ideal",
            linewidth=1.2,
        )
        ax.plot(
            times,
            [rf(r, "body2_y_ref") - y for r, y in zip(reference_rows, ideal_y_on_ref)],
            label="RecurDyn clearance body2 y - ideal",
            linewidth=1.0,
            linestyle="--",
        )
        ax.plot(
            times,
            [rf(r, "body2_z_ref") - z for r, z in zip(reference_rows, ideal_z_on_ref)],
            label="RecurDyn clearance body2 z - ideal",
            linewidth=1.0,
            linestyle="--",
        )
        ax.set_xlabel("time (s)")
        ax.set_ylabel("position error vs ideal hinge (m)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    path = case_dir / "body3_position_vs_recurdyn.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(times, [rf(r, "body3_y_ref") for r in reference_rows], label="RecurDyn body3 y", linewidth=1.4)
    ax.plot(times, [rf(r, "body3_y_sdf_ncp") for r in reference_rows], label=f"{result_label} body3 y", linewidth=1.4)
    ax.plot(times, [rf(r, "body3_z_ref") for r in reference_rows], label="RecurDyn body3 z", linewidth=1.0)
    ax.plot(times, [rf(r, "body3_z_sdf_ncp") for r in reference_rows], label=f"{result_label} body3 z", linewidth=1.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("position (m)")
    ax.set_title(case_name)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    path = case_dir / "position_error_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(times, [rf(r, "body2_y_error") for r in reference_rows], label="body2 y error", linewidth=1.4)
    ax.plot(times, [rf(r, "body2_z_error") for r in reference_rows], label="body2 z error", linewidth=1.2)
    ax.plot(times, [rf(r, "body3_y_error") for r in reference_rows], label="body3 y error", linewidth=1.4)
    ax.plot(times, [rf(r, "body3_z_error") for r in reference_rows], label="body3 z error", linewidth=1.2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("SDF-NCP - RecurDyn position (m)")
    ax.set_title(case_name)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    path = case_dir / "body3_xy_orbit_vs_recurdyn.png"
    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    ax.plot(
        [rf(r, "body3_x_ref") for r in reference_rows],
        [rf(r, "body3_y_ref") for r in reference_rows],
        label="RecurDyn body3",
        linewidth=1.4,
    )
    ax.plot(
        [rf(r, "body3_x_sdf_ncp") for r in reference_rows],
        [rf(r, "body3_y_sdf_ncp") for r in reference_rows],
        label=f"{result_label} body3",
        linewidth=1.4,
    )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(case_name)
    ax.grid(True, alpha=0.25)
    ax.axis("equal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    kinematic_rows = read_optional_csv(RESULT_DIR / case_name / "recurdyn_chrono_kinematic_diagnostics.csv")
    if kinematic_rows:
        ktimes = [rf(r, "time") for r in kinematic_rows]

        path = case_dir / "relative_vector_angle_error_vs_time.png"
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(
            ktimes,
            [rf(r, "rel_vector_angle_error_deg") for r in kinematic_rows],
            label="Body2-Body3 relative-vector angle error",
            linewidth=1.4,
        )
        ax.set_xlabel("time (s)")
        ax.set_ylabel("angle error (deg)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

        path = case_dir / "contact_wrench_vs_time.png"
        fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
        axes[0].plot(ktimes, [rf(r, "contact_force_norm") for r in kinematic_rows], label="|contact force|")
        axes[1].plot(ktimes, [rf(r, "contact_torque_norm") for r in kinematic_rows], label="|contact torque|")
        axes[0].set_ylabel("force (N)")
        axes[1].set_ylabel("torque (N m)")
        axes[1].set_xlabel("time (s)")
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.legend()
        axes[0].set_title(case_name)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

        path = case_dir / "action_marker_contact_center_yz.png"
        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        ax.plot(
            [rf(r, "action_marker_sim_y") for r in kinematic_rows],
            [rf(r, "action_marker_sim_z") for r in kinematic_rows],
            label=f"{result_label} ActionMarker",
            linewidth=1.2,
        )
        ax.plot(
            [rf(r, "action_marker_ref_proxy_y") for r in kinematic_rows],
            [rf(r, "action_marker_ref_proxy_z") for r in kinematic_rows],
            label="RecurDyn ActionMarker proxy",
            linewidth=1.0,
            linestyle="--",
        )
        ax.plot(
            [rf(r, "contact_center_y") for r in kinematic_rows if rf(r, "contact_force_norm") > 0.0],
            [rf(r, "contact_center_z") for r in kinematic_rows if rf(r, "contact_force_norm") > 0.0],
            label="contact center",
            linewidth=1.0,
            alpha=0.8,
        )
        ax.set_xlabel("y (m)")
        ax.set_ylabel("z (m)")
        ax.set_title(case_name)
        ax.grid(True, alpha=0.25)
        ax.axis("equal")
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

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
    if case_name.startswith("rev_joint_clearance"):
        outputs.extend(plot_rev_joint_clearance_extra(case_name, rows))
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
