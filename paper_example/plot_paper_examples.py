import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".", help="Repository root")
    return parser.parse_args()


def load_csv(path):
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def load_json(path):
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def rows_for(rows, case_name):
    return [row for row in rows if row["case"] == case_name]


def col(rows, name):
    return [float(row[name]) for row in rows]


def elastic_post(v1, v2, m1, m2):
    post_1 = ((m1 - m2) / (m1 + m2)) * v1 + (2.0 * m2 / (m1 + m2)) * v2
    post_2 = (2.0 * m1 / (m1 + m2)) * v1 + ((m2 - m1) / (m1 + m2)) * v2
    return post_1, post_2


def headon_reference_velocities(model, times):
    radius = float(model["sphere_radius"])
    ax0 = float(model["sphere_a_init_pos"][0])
    bx0 = float(model["sphere_b_init_pos"][0])
    av0 = float(model["sphere_a_init_vel"][0])
    bv0 = float(model["sphere_b_init_vel"][0])
    rho_a = float(model["sphere_a_density"])
    rho_b = float(model["sphere_b_density"])
    volume = 4.0 * math.pi * radius**3 / 3.0
    av_post, bv_post = elastic_post(av0, bv0, rho_a * volume, rho_b * volume)
    impact_time = (bx0 - ax0 - 2.0 * radius) / (av0 - bv0)

    ref_av = []
    ref_bv = []
    for time in times:
        if time <= impact_time:
            ref_av.append(av0)
            ref_bv.append(bv0)
        else:
            ref_av.append(av_post)
            ref_bv.append(bv_post)
    return ref_av, ref_bv, impact_time


def save_figure(fig, output_base):
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_cam_case(rows, case_name, output_dir):
    times = col(rows, "time")
    backend = col(rows, "backend_y")
    reference = col(rows, "reference_y")
    error = col(rows, "y_error")
    min_phi = col(rows, "min_phi")

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.0), sharex=True)
    axes[0].plot(times, backend, label="backend dynamic", linewidth=1.8)
    axes[0].plot(times, reference, "--", label="reference", linewidth=1.5)
    axes[0].set_ylabel("follower y (m)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(times, error, color="#b7410e", linewidth=1.4)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("y error (m)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, min_phi, color="#2f6f4e", linewidth=1.4)
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("min phi (m)")
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(case_name)
    fig.tight_layout()
    save_figure(fig, output_dir / f"{case_name}_curves")


def plot_headon_case(rows, model, case_name, output_dir):
    times = col(rows, "time")
    ax = col(rows, "sphere_a_x")
    bx = col(rows, "sphere_b_x")
    av = col(rows, "sphere_a_vx")
    bv = col(rows, "sphere_b_vx")
    ref_ax = col(rows, "reference_sphere_a_x")
    ref_bx = col(rows, "reference_sphere_b_x")
    ref_av, ref_bv, impact_time = headon_reference_velocities(model, times)

    err_a = [a - r for a, r in zip(ax, ref_ax)]
    err_b = [b - r for b, r in zip(bx, ref_bx)]

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.5), sharex=True)
    axes[0].plot(times, ax, label="A backend", linewidth=1.8)
    axes[0].plot(times, bx, label="B backend", linewidth=1.8)
    axes[0].plot(times, ref_ax, "--", label="A reference", linewidth=1.4)
    axes[0].plot(times, ref_bx, "--", label="B reference", linewidth=1.4)
    axes[0].axvline(impact_time, color="#666666", linestyle=":", linewidth=1.2)
    axes[0].set_ylabel("x (m)")
    axes[0].legend(ncol=2)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(times, av, label="A backend", linewidth=1.8)
    axes[1].plot(times, bv, label="B backend", linewidth=1.8)
    axes[1].plot(times, ref_av, "--", label="A reference", linewidth=1.4)
    axes[1].plot(times, ref_bv, "--", label="B reference", linewidth=1.4)
    axes[1].axvline(impact_time, color="#666666", linestyle=":", linewidth=1.2)
    axes[1].set_ylabel("vx (m/s)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, err_a, label="A x error", linewidth=1.4)
    axes[2].plot(times, err_b, label="B x error", linewidth=1.4)
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("x error (m)")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(case_name)
    fig.tight_layout()
    save_figure(fig, output_dir / f"{case_name}_curves")


def plot_simple_gear_case(rows, case_name, output_dir):
    times = col(rows, "time")
    omega = col(rows, "backend_y")
    reference = col(rows, "reference_y")
    error = col(rows, "y_error")
    patch_count = col(rows, "patch_count")

    fig, axes = plt.subplots(4, 1, figsize=(9.0, 9.5), sharex=False)
    axes[0].plot(times, omega, label="backend dynamic", linewidth=1.6)
    axes[0].plot(times, reference, "--", label="reference", linewidth=1.4)
    axes[0].set_ylabel(r"GEAR22 $\omega_x$ (rad/s)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    zoom_rows = [i for i, time in enumerate(times) if time <= 0.05]
    axes[1].plot([times[i] for i in zoom_rows], [omega[i] for i in zoom_rows],
                 label="backend dynamic", linewidth=1.6)
    axes[1].plot([times[i] for i in zoom_rows], [reference[i] for i in zoom_rows],
                 "--", label="reference", linewidth=1.4)
    axes[1].set_ylabel(r"0-0.05 s $\omega_x$")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, error, color="#b7410e", linewidth=1.2)
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_ylabel("error (rad/s)")
    axes[2].grid(True, alpha=0.3)

    axes[3].step(times, patch_count, where="post", color="#2f6f4e", linewidth=1.2)
    axes[3].set_xlabel("time (s)")
    axes[3].set_ylabel("patch count")
    axes[3].grid(True, alpha=0.3)

    fig.suptitle(case_name)
    fig.tight_layout()
    save_figure(fig, output_dir / f"{case_name}_curves")


def main():
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    paper_dir = project_root / "paper_example"
    out_dir = project_root / "out" / "paper_example_dynamic_benchmarks"
    figure_dir = paper_dir / "figures"

    frames = load_csv(out_dir / "sparse_sdf_frames.csv")

    for case_name in ("eccentric_roller",):
        plot_cam_case(rows_for(frames, case_name), case_name, figure_dir)

    for case_name in ("headon_spheres", "headon_spheres_mass_ratio"):
        model = load_json(paper_dir / "cases" / case_name / f"{case_name}_model.json")
        plot_headon_case(rows_for(frames, case_name), model, case_name, figure_dir)

    plot_simple_gear_case(rows_for(frames, "simple_gear"), "simple_gear", figure_dir)

    print(f"Wrote figures to {figure_dir}")


if __name__ == "__main__":
    main()
