"""Compare ideal revolute joint reaction wrench against contact patch wrenches.

This is a diagnostic script for the rev_joint_clearance benchmark.  The ideal
Chrono revolute joint provides the constraint wrench needed by the same RMD
front-end without clearance contact.  The SDF-NCP and GGEOMCONTACT calibration
runs provide the actual integrated contact patch wrench.  Comparing the two
isolates contact-patch force/torque delivery from front-end kinematics.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = ROOT / "results" / "sdf_ncp_benchmarks"
IDEAL_CASE = "rev_joint_clearance_ideal_revolute"
CONTACT_CASES = {
    "sdf_ncp": "rev_joint_clearance",
    "ggeomcontact_calibration": "rev_joint_clearance_ggeomcontact_calibration",
    "descriptor_patch_delassus": "rev_joint_clearance_descriptor_patch_delassus",
    "descriptor_patch_delassus_laplacian": "rev_joint_clearance_descriptor_patch_delassus_laplacian",
    "descriptor_patch_velocity_delassus": "rev_joint_clearance_descriptor_patch_velocity_delassus",
    "descriptor_patch_velocity_delassus_laplacian": "rev_joint_clearance_descriptor_patch_velocity_delassus_laplacian",
    "local_patch_delassus_wrench_closure": "rev_joint_clearance_local_patch_wrench_closure",
    "bidirectional_local_patch_wrench_closure": "rev_joint_clearance_bidirectional_local_patch_wrench_closure",
}
OUT_DIR = RESULT_ROOT / "rev_joint_clearance_wrench_comparison"
FIG_DIR = RESULT_ROOT / "figures" / "rev_joint_clearance_wrench_comparison"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def norm3(x: float, y: float, z: float) -> float:
    return math.sqrt(x * x + y * y + z * z)


def dot3(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cosine(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    an = norm3(*a)
    bn = norm3(*b)
    if an < 1.0e-12 or bn < 1.0e-12:
        return float("nan")
    return max(-1.0, min(1.0, dot3(a, b) / (an * bn)))


def interp_scalar(xs: list[float], ys: list[float], x: float) -> float:
    if not xs:
        return float("nan")
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo = 0
    hi = len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - xs[lo]) / max(xs[hi] - xs[lo], 1.0e-30)
    return ys[lo] * (1.0 - t) + ys[hi] * t


def rms(values: list[float]) -> float:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return float("nan")
    return math.sqrt(sum(v * v for v in finite) / len(finite))


def max_abs(values: list[float]) -> float:
    finite = [abs(v) for v in values if math.isfinite(v)]
    return max(finite) if finite else float("nan")


def mean(values: list[float]) -> float:
    finite = [v for v in values if math.isfinite(v)]
    return sum(finite) / len(finite) if finite else float("nan")


def extract_ideal() -> list[dict[str, float]]:
    path = RESULT_ROOT / IDEAL_CASE / "joint_reaction_wrench.csv"
    try:
        rows = read_rows(path)
    except FileNotFoundError as exc:
        raise SystemExit(
            f"缺少 {path}\n"
            "请先运行：build\\bin\\Release\\demo_CH_sdf_ncp_benchmarks_openvdb.exe "
            "rev_joint_clearance_ideal_revolute"
        ) from exc

    out = []
    for row in rows:
        if row.get("body") != "body3":
            continue
        fx, fy, fz = f(row, "force_x"), f(row, "force_y"), f(row, "force_z")
        tx = f(row, "torque_about_body_ref_x")
        ty = f(row, "torque_about_body_ref_y")
        tz = f(row, "torque_about_body_ref_z")
        out.append(
            {
                "time": f(row, "time"),
                "force_x": fx,
                "force_y": fy,
                "force_z": fz,
                "force_norm": norm3(fx, fy, fz),
                "torque_x": tx,
                "torque_y": ty,
                "torque_z": tz,
                "torque_norm": norm3(tx, ty, tz),
            }
        )
    if not out:
        raise SystemExit(f"{path} 中没有 body3 reaction wrench 行。")
    return out


def extract_contact(case_name: str) -> list[dict[str, float]]:
    path = RESULT_ROOT / case_name / "recurdyn_chrono_kinematic_diagnostics.csv"
    try:
        rows = read_rows(path)
    except FileNotFoundError as exc:
        raise SystemExit(
            f"缺少 {path}\n"
            f"请先运行：build\\bin\\Release\\demo_CH_sdf_ncp_benchmarks_openvdb.exe {case_name}"
        ) from exc
    out = []
    for row in rows:
        fx, fy, fz = f(row, "contact_force_x"), f(row, "contact_force_y"), f(row, "contact_force_z")
        tx = f(row, "contact_torque_world_x")
        ty = f(row, "contact_torque_world_y")
        tz = f(row, "contact_torque_world_z")
        out.append(
            {
                "time": f(row, "time"),
                "force_x": fx,
                "force_y": fy,
                "force_z": fz,
                "force_norm": norm3(fx, fy, fz),
                "torque_x": tx,
                "torque_y": ty,
                "torque_z": tz,
                "torque_norm": norm3(tx, ty, tz),
            }
        )
    if not out:
        raise SystemExit(f"{path} 没有可用 contact wrench 行。")
    return out


def interpolate_contact(contact_rows: list[dict[str, float]], time: float) -> dict[str, float]:
    xs = [r["time"] for r in contact_rows]
    keys = ("force_x", "force_y", "force_z", "force_norm", "torque_x", "torque_y", "torque_z", "torque_norm")
    return {key: interp_scalar(xs, [r[key] for r in contact_rows], time) for key in keys}


def write_comparison(ideal_rows: list[dict[str, float]], contacts: dict[str, list[dict[str, float]]]) -> tuple[Path, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    comparison_path = OUT_DIR / "ideal_vs_contact_wrench.csv"
    summary_path = OUT_DIR / "wrench_summary.csv"

    fieldnames = [
        "time",
        "method",
        "ideal_force_x",
        "ideal_force_y",
        "ideal_force_z",
        "ideal_force_norm",
        "contact_force_x",
        "contact_force_y",
        "contact_force_z",
        "contact_force_norm",
        "force_error_x",
        "force_error_y",
        "force_error_z",
        "force_error_norm",
        "force_cosine",
        "ideal_torque_x",
        "ideal_torque_y",
        "ideal_torque_z",
        "ideal_torque_norm",
        "contact_torque_x",
        "contact_torque_y",
        "contact_torque_z",
        "contact_torque_norm",
        "torque_error_x",
        "torque_error_y",
        "torque_error_z",
        "torque_error_norm",
        "torque_cosine",
    ]
    summary_rows = []
    with comparison_path.open("w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=fieldnames)
        writer.writeheader()
        for method, contact_rows in contacts.items():
            force_errors = []
            torque_errors = []
            force_norm_errors = []
            torque_norm_errors = []
            force_cosines = []
            torque_cosines = []
            contact_force_norms = []
            contact_torque_norms = []
            for ideal in ideal_rows:
                c = interpolate_contact(contact_rows, ideal["time"])
                ferr = (
                    c["force_x"] - ideal["force_x"],
                    c["force_y"] - ideal["force_y"],
                    c["force_z"] - ideal["force_z"],
                )
                terr = (
                    c["torque_x"] - ideal["torque_x"],
                    c["torque_y"] - ideal["torque_y"],
                    c["torque_z"] - ideal["torque_z"],
                )
                fnerr = norm3(*ferr)
                tnerr = norm3(*terr)
                fcos = cosine(
                    (ideal["force_x"], ideal["force_y"], ideal["force_z"]),
                    (c["force_x"], c["force_y"], c["force_z"]),
                )
                tcos = cosine(
                    (ideal["torque_x"], ideal["torque_y"], ideal["torque_z"]),
                    (c["torque_x"], c["torque_y"], c["torque_z"]),
                )
                force_errors.extend(ferr)
                torque_errors.extend(terr)
                force_norm_errors.append(fnerr)
                torque_norm_errors.append(tnerr)
                force_cosines.append(fcos)
                torque_cosines.append(tcos)
                contact_force_norms.append(c["force_norm"])
                contact_torque_norms.append(c["torque_norm"])
                writer.writerow(
                    {
                        "time": ideal["time"],
                        "method": method,
                        "ideal_force_x": ideal["force_x"],
                        "ideal_force_y": ideal["force_y"],
                        "ideal_force_z": ideal["force_z"],
                        "ideal_force_norm": ideal["force_norm"],
                        "contact_force_x": c["force_x"],
                        "contact_force_y": c["force_y"],
                        "contact_force_z": c["force_z"],
                        "contact_force_norm": c["force_norm"],
                        "force_error_x": ferr[0],
                        "force_error_y": ferr[1],
                        "force_error_z": ferr[2],
                        "force_error_norm": fnerr,
                        "force_cosine": fcos,
                        "ideal_torque_x": ideal["torque_x"],
                        "ideal_torque_y": ideal["torque_y"],
                        "ideal_torque_z": ideal["torque_z"],
                        "ideal_torque_norm": ideal["torque_norm"],
                        "contact_torque_x": c["torque_x"],
                        "contact_torque_y": c["torque_y"],
                        "contact_torque_z": c["torque_z"],
                        "contact_torque_norm": c["torque_norm"],
                        "torque_error_x": terr[0],
                        "torque_error_y": terr[1],
                        "torque_error_z": terr[2],
                        "torque_error_norm": tnerr,
                        "torque_cosine": tcos,
                    }
                )
            summary_rows.append(
                {
                    "method": method,
                    "num_samples": len(ideal_rows),
                    "ideal_force_norm_max": max(r["force_norm"] for r in ideal_rows),
                    "contact_force_norm_max": max(contact_force_norms),
                    "force_error_norm_rmse": rms(force_norm_errors),
                    "force_error_norm_max": max_abs(force_norm_errors),
                    "force_component_rmse": rms(force_errors),
                    "force_cosine_mean": mean(force_cosines),
                    "ideal_torque_norm_max": max(r["torque_norm"] for r in ideal_rows),
                    "contact_torque_norm_max": max(contact_torque_norms),
                    "torque_error_norm_rmse": rms(torque_norm_errors),
                    "torque_error_norm_max": max_abs(torque_norm_errors),
                    "torque_component_rmse": rms(torque_errors),
                    "torque_cosine_mean": mean(torque_cosines),
                }
            )

    with summary_path.open("w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    return comparison_path, summary_path


def plot_outputs(ideal_rows: list[dict[str, float]], contacts: dict[str, list[dict[str, float]]]) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    times = [r["time"] for r in ideal_rows]

    def contact_series(method: str, key: str) -> list[float]:
        return [interpolate_contact(contacts[method], t)[key] for t in times]

    path = FIG_DIR / "force_norm_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(times, [r["force_norm"] for r in ideal_rows], label="ideal revolute reaction", linewidth=1.5)
    for method in contacts:
        ax.plot(times, contact_series(method, "force_norm"), label=method, linewidth=1.2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("force norm (N)")
    ax.set_title("rev_joint_clearance wrench comparison")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    path = FIG_DIR / "torque_norm_about_body3_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(times, [r["torque_norm"] for r in ideal_rows], label="ideal revolute reaction", linewidth=1.5)
    for method in contacts:
        ax.plot(times, contact_series(method, "torque_norm"), label=method, linewidth=1.2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("torque norm about Body3 ref (N m)")
    ax.set_title("rev_joint_clearance wrench comparison")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    for kind, keys, ylabel, filename in [
        ("force", ("force_x", "force_y", "force_z"), "force (N)", "force_components_vs_time.png"),
        ("torque", ("torque_x", "torque_y", "torque_z"), "torque about Body3 ref (N m)", "torque_components_vs_time.png"),
    ]:
        path = FIG_DIR / filename
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 7.2), sharex=True)
        labels = ("x", "y", "z")
        for ax, key, label in zip(axes, keys, labels):
            ax.plot(times, [r[key] for r in ideal_rows], label=f"ideal {kind}_{label}", linewidth=1.4)
            for method in contacts:
                ax.plot(times, contact_series(method, key), label=f"{method} {kind}_{label}", linewidth=1.0)
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=8)
        axes[-1].set_xlabel("time (s)")
        axes[0].set_title("rev_joint_clearance wrench components")
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    path = FIG_DIR / "wrench_error_norm_vs_time.png"
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
    for method in contacts:
        force_err = []
        torque_err = []
        for ideal in ideal_rows:
            c = interpolate_contact(contacts[method], ideal["time"])
            force_err.append(
                norm3(
                    c["force_x"] - ideal["force_x"],
                    c["force_y"] - ideal["force_y"],
                    c["force_z"] - ideal["force_z"],
                )
            )
            torque_err.append(
                norm3(
                    c["torque_x"] - ideal["torque_x"],
                    c["torque_y"] - ideal["torque_y"],
                    c["torque_z"] - ideal["torque_z"],
                )
            )
        axes[0].plot(times, force_err, label=method)
        axes[1].plot(times, torque_err, label=method)
    axes[0].set_ylabel("force error norm (N)")
    axes[1].set_ylabel("torque error norm (N m)")
    axes[1].set_xlabel("time (s)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[0].set_title("contact wrench minus ideal revolute reaction")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    return outputs


def extract_patch_global(case_name: str) -> list[dict[str, float]]:
    path = RESULT_ROOT / case_name / "contact_patch_moment_diagnostics.csv"
    if not path.exists():
        return []
    rows = []
    for row in read_rows(path):
        if row.get("patch_id") != "-1":
            continue
        rows.append(
            {
                "time": f(row, "time"),
                "active_samples": f(row, "active_samples"),
                "raw_force_norm": f(row, "raw_force_norm"),
                "weighted_force_norm": f(row, "weighted_force_norm"),
                "effective_force_norm": f(row, "effective_force_norm"),
                "raw_torque_norm": f(row, "raw_torque_norm"),
                "weighted_torque_norm": f(row, "weighted_torque_norm"),
                "effective_torque_norm": f(row, "effective_torque_norm"),
                "resultant_line_offset_m": f(row, "resultant_line_offset_m"),
                "closure_ratio_max_arm": f(row, "closure_ratio_max_arm"),
            }
        )
    return rows


def plot_patch_moment_diagnostics() -> list[Path]:
    outputs: list[Path] = []
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    patch_rows = {method: extract_patch_global(case_name) for method, case_name in CONTACT_CASES.items()}
    patch_rows = {method: rows for method, rows in patch_rows.items() if rows}
    if not patch_rows:
        return outputs

    path = FIG_DIR / "patch_resultant_line_offset_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for method, rows in patch_rows.items():
        ax.plot([r["time"] for r in rows], [r["resultant_line_offset_m"] for r in rows], label=method)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("|T| / |F| about Body3 ref (m)")
    ax.set_title("rev_joint_clearance patch wrench line offset")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    path = FIG_DIR / "patch_closure_ratio_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for method, rows in patch_rows.items():
        ax.plot([r["time"] for r in rows], [r["closure_ratio_max_arm"] for r in rows], label=method)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("|T| / (|F| max |r|)")
    ax.set_title("rev_joint_clearance patch torque closure ratio")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    if "sdf_ncp" in patch_rows:
        rows = patch_rows["sdf_ncp"]
        times = [r["time"] for r in rows]
        path = FIG_DIR / "sdf_ncp_raw_vs_weighted_patch_wrench.png"
        fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
        axes[0].plot(times, [r["raw_force_norm"] for r in rows], label="raw descriptor multiplier")
        axes[0].plot(times, [r["weighted_force_norm"] for r in rows], label="area-weighted diagnostic")
        axes[0].plot(times, [r["effective_force_norm"] for r in rows], label="effective residual force", linestyle="--")
        axes[1].plot(times, [r["raw_torque_norm"] for r in rows], label="raw descriptor multiplier")
        axes[1].plot(times, [r["weighted_torque_norm"] for r in rows], label="area-weighted diagnostic")
        axes[1].plot(times, [r["effective_torque_norm"] for r in rows], label="effective residual torque", linestyle="--")
        axes[0].set_ylabel("force norm (N)")
        axes[1].set_ylabel("torque norm (N m)")
        axes[1].set_xlabel("time (s)")
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.legend()
        axes[0].set_title("SDF-NCP descriptor raw vs area-weighted patch wrench")
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    path = FIG_DIR / "active_samples_vs_time.png"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for method, rows in patch_rows.items():
        ax.plot([r["time"] for r in rows], [r["active_samples"] for r in rows], label=method)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("active samples")
    ax.set_title("rev_joint_clearance active contact samples")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    outputs.append(path)

    return outputs


def main() -> int:
    ideal_rows = extract_ideal()
    contacts = {method: extract_contact(case_name) for method, case_name in CONTACT_CASES.items()}
    comparison_path, summary_path = write_comparison(ideal_rows, contacts)
    figures = plot_outputs(ideal_rows, contacts)
    figures.extend(plot_patch_moment_diagnostics())
    print(f"Wrote {comparison_path}")
    print(f"Wrote {summary_path}")
    for path in figures:
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
