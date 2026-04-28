import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "SimSun", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.9,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.unicode_minus": False,
    }
)

LINE_STACK = [
    {"lw": 1.45, "ls": "-", "color": "#4D4D4D", "zorder": 1},
    {"lw": 1.32, "ls": (0, (10, 4)), "color": "#E64B35", "zorder": 2},
    {"lw": 1.22, "ls": (0, (10, 3, 2, 3)), "color": "#2CA02C", "zorder": 3},
    {"lw": 1.12, "ls": (0, (6, 2)), "color": "#7A3DB8", "zorder": 4},
    {"lw": 1.05, "ls": (0, (6, 2, 1.5, 2)), "color": "#1F4EFA", "zorder": 5},
    {"lw": 0.95, "ls": (0, (2.5, 1.8)), "color": "#00CFEA", "zorder": 6},
]
CHINESE_FONT = FontProperties(family="SimSun")


def apply_axes_style(ax):
    ax.grid(True, alpha=0.16, linewidth=0.6)
    ax.tick_params(direction="in", width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)


def add_legend(ax, **kwargs):
    defaults = {
        "frameon": True,
        "facecolor": "white",
        "edgecolor": "0.55",
        "framealpha": 1.0,
        "borderpad": 0.35,
        "handlelength": 2.8,
    }
    defaults.update(kwargs)
    return ax.legend(**defaults)


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


def optional_col(rows, name):
    values = []
    for row in rows:
        value = row.get(name, "")
        if value == "":
            return None
        values.append(float(value))
    return values


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


def save_figure(fig, output_base, tight=True):
    output_base.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(output_base.with_suffix(".pdf"), **save_kwargs)
    plt.close(fig)


def apply_compact_stack_layout(fig):
    fig.subplots_adjust(left=0.205, right=0.985, top=0.972, bottom=0.118, hspace=0.155)


def plot_cam_case(rows, case_name, output_dir):
    times = col(rows, "time")
    backend = col(rows, "backend_y")
    reference = col(rows, "reference_y")
    error = col(rows, "y_error")
    min_phi = col(rows, "min_phi")

    compact_figsize = (3.25, 4.45)
    fig, axes = plt.subplots(3, 1, figsize=compact_figsize, sharex=True)
    axes[0].plot(times, backend, label="backend dynamic", **LINE_STACK[4])
    axes[0].plot(times, reference, label="reference", **LINE_STACK[5])
    axes[0].set_ylabel("follower y (m)")
    add_legend(axes[0], loc="best", fontsize=6.5, handlelength=2.2)
    apply_axes_style(axes[0])

    axes[1].plot(times, error, color="#E64B35", linewidth=1.1, linestyle=(0, (10, 4)))
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("y error (m)")
    apply_axes_style(axes[1])

    axes[2].plot(times, min_phi, color="#2CA02C", linewidth=1.1, linestyle=(0, (10, 3, 2, 3)))
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("min phi (m)")
    apply_axes_style(axes[2])

    apply_compact_stack_layout(fig)
    save_figure(fig, output_dir / f"{case_name}_curves", tight=False)


def plot_headon_case(rows, model, case_name, output_dir, compact=False):
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

    figsize = (3.25, 4.45) if compact else (6.8, 4.8)
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    axes[0].plot(times, ax, label="A backend", **LINE_STACK[4])
    axes[0].plot(times, bx, label="B backend", **LINE_STACK[1])
    axes[0].plot(times, ref_ax, label="A reference", **LINE_STACK[5])
    axes[0].plot(times, ref_bx, label="B reference", **LINE_STACK[3])
    axes[0].axvline(impact_time, color="#666666", linestyle=":", linewidth=1.2)
    axes[0].set_ylabel("x (m)")
    add_legend(axes[0], ncol=1 if compact else 2, loc="best", fontsize=6.2 if compact else 7.5, handlelength=2.1)
    apply_axes_style(axes[0])

    axes[1].plot(times, av, label="A backend", **LINE_STACK[4])
    axes[1].plot(times, bv, label="B backend", **LINE_STACK[1])
    axes[1].plot(times, ref_av, label="A reference", **LINE_STACK[5])
    axes[1].plot(times, ref_bv, label="B reference", **LINE_STACK[3])
    axes[1].axvline(impact_time, color="#666666", linestyle=":", linewidth=1.2)
    axes[1].set_ylabel("vx (m/s)")
    apply_axes_style(axes[1])

    axes[2].plot(times, err_a, label="A x error", color="#1F4EFA", linewidth=1.05, linestyle=(0, (6, 2, 1.5, 2)))
    axes[2].plot(times, err_b, label="B x error", color="#E64B35", linewidth=1.15, linestyle=(0, (10, 4)))
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("x error (m)")
    add_legend(axes[2], loc="best", fontsize=6.2 if compact else 7.5, handlelength=2.1)
    apply_axes_style(axes[2])

    if compact:
        apply_compact_stack_layout(fig)
        save_figure(fig, output_dir / f"{case_name}_curves", tight=False)
    else:
        fig.tight_layout(pad=0.25)
        save_figure(fig, output_dir / f"{case_name}_curves")


def plot_simple_gear_case(rows, case_name, output_dir):
    times = col(rows, "time")
    omega = col(rows, "backend_y")
    analytic = col(rows, "reference_y")
    commercial_baseline = optional_col(rows, "recur" + "dyn_y")
    error = [w - r for w, r in zip(omega, analytic)]
    commercial_error = None
    if commercial_baseline is not None:
        commercial_error = [r - a for r, a in zip(commercial_baseline, analytic)]
    patch_count = col(rows, "patch_count")

    fig, axes = plt.subplots(4, 1, figsize=(9.0, 9.5), sharex=False)
    axes[0].plot(times, omega, label="backend dynamic", **LINE_STACK[4])
    axes[0].plot(times, analytic, label="analytic -1 rad/s", **LINE_STACK[5])
    if commercial_baseline is not None:
        axes[0].plot(times, commercial_baseline, label="commercial software", **LINE_STACK[1])
    axes[0].set_ylabel("follower gear omega_x (rad/s)")
    add_legend(axes[0], loc="best")
    apply_axes_style(axes[0])

    zoom_rows = [i for i, time in enumerate(times) if time <= 0.05]
    axes[1].plot([times[i] for i in zoom_rows], [omega[i] for i in zoom_rows],
                 label="backend dynamic", **LINE_STACK[4])
    axes[1].plot([times[i] for i in zoom_rows], [analytic[i] for i in zoom_rows],
                 label="analytic -1 rad/s", **LINE_STACK[5])
    if commercial_baseline is not None:
        axes[1].plot([times[i] for i in zoom_rows], [commercial_baseline[i] for i in zoom_rows],
                     label="commercial software", **LINE_STACK[1])
    axes[1].set_ylabel(r"0-0.05 s $\omega_x$")
    add_legend(axes[1], loc="best")
    apply_axes_style(axes[1])

    axes[2].plot(times, error, color="#1F4EFA", linewidth=1.05, linestyle=(0, (6, 2, 1.5, 2)), label="backend - analytic")
    if commercial_error is not None:
        axes[2].plot(times, commercial_error, color="#E64B35", linestyle=(0, (10, 4)), linewidth=1.15,
                     label="commercial software - analytic")
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_ylabel("error vs analytic (rad/s)")
    add_legend(axes[2], loc="best")
    apply_axes_style(axes[2])

    axes[3].step(times, patch_count, where="post", color="#2CA02C", linewidth=1.0)
    axes[3].set_xlabel("time (s)")
    axes[3].set_ylabel("patch count")
    apply_axes_style(axes[3])

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

    plot_headon_case(
        rows_for(frames, "headon_spheres"),
        load_json(paper_dir / "cases" / "headon_spheres" / "headon_spheres_model.json"),
        "headon_spheres",
        figure_dir,
        compact=True,
    )

    for case_name in ("headon_spheres_mass_ratio",):
        model = load_json(paper_dir / "cases" / case_name / f"{case_name}_model.json")
        plot_headon_case(rows_for(frames, case_name), model, case_name, figure_dir, compact=True)

    plot_simple_gear_case(rows_for(frames, "simple_gear"), "simple_gear", figure_dir)

    print(f"Wrote figures to {figure_dir}")


if __name__ == "__main__":
    main()
