"""Plot JCND validation and ablation CSVs as paper-ready PDFs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def first_existing(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        if key in row and row[key] != "":
            return key
    raise KeyError(keys[0])


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_headon(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.35))

    physics = next(r for r in rows if r["suite"] == "physics")
    metrics = [
        ("e", as_float(physics, "restitution_error")),
        ("energy", as_float(physics, "energy_rel_error")),
        ("vel.", as_float(physics, "velocity_l2_error")),
    ]
    axes[0].bar([m[0] for m in metrics], [m[1] for m in metrics], color="#4c78a8")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("error")
    axes[0].set_title("impact physics")
    axes[0].grid(True, axis="y", alpha=0.25)

    sdf_rows = [r for r in rows if r["suite"] == "sdf_resolution"]
    sdf_rows.sort(key=lambda r: as_float(r, "voxel_size"), reverse=True)
    axes[1].plot(
        [as_float(r, "voxel_size") * 1e6 for r in sdf_rows],
        [as_float(r, "velocity_l2_error") for r in sdf_rows],
        marker="o",
        lw=1.4,
        color="#2f6f4e",
        label="voxel sweep",
    )
    dt_rows = [r for r in rows if r["suite"] == "time_step"]
    dt_rows.sort(key=lambda r: as_float(r, "contact_substep"), reverse=True)
    ax2 = axes[1].twiny()
    ax2.plot(
        [as_float(r, "contact_substep") * 1e6 for r in dt_rows],
        [as_float(r, "velocity_l2_error") for r in dt_rows],
        marker="s",
        lw=1.2,
        color="#b7410e",
        label="step sweep",
    )
    axes[1].set_xlabel("voxel size (um)")
    ax2.set_xlabel("substep (microsecond)")
    axes[1].set_ylabel("velocity L2 error")
    axes[1].set_title("resolution / step")
    axes[1].grid(True, alpha=0.25)

    ablation = [r for r in rows if r["suite"] == "bidirectional_ablation"]
    labels = [r["mode"].replace("_", "\n") for r in ablation]
    values = [as_float(r, "velocity_l2_error") for r in ablation]
    axes[2].bar(labels, values, color=["#4c78a8", "#b7410e", "#6f4c9b"])
    axes[2].set_ylabel("velocity L2 error")
    axes[2].set_title("one-way ablation")
    axes[2].grid(True, axis="y", alpha=0.25)

    fig.tight_layout(w_pad=1.0)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_gear(summary_rows: list[dict[str, str]], frame_rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.5))
    axes = axes.ravel()

    modes = [r["mode"] for r in summary_rows]
    colors = {
        "bidirectional_distributed": "#4c78a8",
        "one_way_gear22_surface": "#b7410e",
        "one_way_gear21_surface": "#6f4c9b",
        "symmetric_single_point": "#2f6f4e",
        "sdf_max_penetration_point": "#8c564b",
    }

    for mode in modes:
        rows = [r for r in frame_rows if r["mode"] == mode]
        rows.sort(key=lambda r: as_float(r, "time"))
        axes[0].plot(
            [as_float(r, "time") for r in rows],
            [as_float(r, first_existing(r, "analytic_error", "error")) for r in rows],
            lw=1.1,
            label=mode.replace("_", " "),
            color=colors.get(mode),
        )
    axes[0].axhline(0.0, lw=0.8, color="0.2")
    axes[0].set_title("GEAR22 RX error vs analytic")
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("rad/s")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False, ncol=1)

    x = range(len(summary_rows))
    axes[1].bar(x, [as_float(r, "rms_error") for r in summary_rows], color=[colors.get(r["mode"]) for r in summary_rows])
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels([r["mode"].replace("_", "\n") for r in summary_rows], rotation=0)
    axes[1].set_ylabel("RMS error")
    axes[1].set_title("short-window RMS")
    axes[1].grid(True, axis="y", alpha=0.25)

    axes[2].bar(
        x,
        [as_float(r, "max_omega_jump") for r in summary_rows],
        color=[colors.get(r["mode"]) for r in summary_rows],
    )
    axes[2].set_xticks(list(x))
    axes[2].set_xticklabels([r["mode"].replace("_", "\n") for r in summary_rows], rotation=0)
    axes[2].set_ylabel("max jump")
    axes[2].set_title("omega continuity")
    axes[2].grid(True, axis="y", alpha=0.25)

    bidir = [r for r in frame_rows if r["mode"] == "bidirectional_distributed"]
    bidir.sort(key=lambda r: as_float(r, "time"))
    axes[3].step(
        [as_float(r, "time") for r in bidir],
        [as_float(r, "patch_count") for r in bidir],
        where="post",
        lw=1.1,
        color="#4c78a8",
    )
    axes[3].set_xlabel("time (s)")
    axes[3].set_ylabel("patch count")
    axes[3].set_title("bidirectional patches")
    axes[3].grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_gear_sensitivity(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.75))
    axes = axes.ravel()
    groups = [
        ("voxel", "voxel size"),
        ("time_step", "time step"),
        ("pressure", "normal pressure"),
        ("damping", "normal damping"),
        ("bpen", "BPEN"),
    ]
    baseline = next((r for r in rows if r["group"] == "baseline"), None)
    baseline_rms = as_float(baseline, "rms_error") if baseline else None

    for ax, (group, title) in zip(axes, groups):
        group_rows = [r for r in rows if r["group"] == group]
        labels = [r["label"] for r in group_rows]
        rms = [as_float(r, "rms_error") for r in group_rows]
        jumps = [as_float(r, "max_omega_jump") for r in group_rows]
        ax.bar(labels, rms, color="#4c78a8", label="RMS")
        ax.plot(labels, jumps, color="#b7410e", marker="o", lw=1.0, label="max jump")
        if baseline_rms is not None:
            ax.axhline(baseline_rms, color="0.25", lw=0.8, ls=":")
        ax.set_title(title)
        ax.set_ylabel("rad/s")
        ax.tick_params(axis="x", rotation=20)
        ax.grid(True, axis="y", alpha=0.25)
        if group == "voxel":
            ax.legend(frameon=False)

    corr_rows = [r for r in rows if r["group"] != "baseline"]
    labels = [r["group"] + "\n" + r["label"] for r in corr_rows]
    corr = [abs(as_float(r, "patch_jump_omega_jump_correlation")) for r in corr_rows]
    axes[5].bar(range(len(labels)), corr, color="#2f6f4e")
    axes[5].set_xticks(range(len(labels)))
    axes[5].set_xticklabels(labels, rotation=90)
    axes[5].set_ylabel("|corr|")
    axes[5].set_title("patch jump vs omega jump")
    axes[5].grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_spatial_convergence(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    rows = sorted(rows, key=lambda r: as_float(r, "voxel_size"), reverse=True)
    x = [as_float(r, "voxel_size") * 1e6 for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.55))
    axes[0].loglog(x, [as_float(r, "rms_distance_error") for r in rows], marker="o", lw=1.3, label="distance RMS")
    axes[0].loglog(x, [as_float(r, "max_distance_error") for r in rows], marker="s", lw=1.1, label="distance max")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("voxel size (um)")
    axes[0].set_ylabel("distance error (m)")
    axes[0].set_title("SDF distance")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].loglog(x, [as_float(r, "rms_normal_angle") for r in rows], marker="o", lw=1.3, label="normal RMS")
    axes[1].loglog(x, [as_float(r, "force_rel_error") for r in rows], marker="s", lw=1.1, label="force rel.")
    axes[1].loglog(x, [as_float(r, "torque_rel_error") for r in rows], marker="^", lw=1.1, label="torque rel.")
    axes[1].invert_xaxis()
    axes[1].set_xlabel("voxel size (um)")
    axes[1].set_ylabel("error")
    axes[1].set_title("normal / wrench")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_patch_wrench_convergence(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    sampling = sorted(
        [r for r in rows if r["suite"] == "surface_sampling"],
        key=lambda r: as_float(r, "sample_count"),
    )
    voxel = sorted(
        [r for r in rows if r["suite"] == "sdf_voxel"],
        key=lambda r: as_float(r, "voxel_size"),
        reverse=True,
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.65))
    axes[0].loglog(
        [as_float(r, "sample_count") for r in sampling],
        [as_float(r, "exact_quad_force_rel_error") for r in sampling],
        marker="o",
        lw=1.2,
        label="force",
    )
    axes[0].loglog(
        [as_float(r, "sample_count") for r in sampling],
        [as_float(r, "exact_quad_torque_rel_error") for r in sampling],
        marker="s",
        lw=1.2,
        label="torque",
    )
    axes[0].set_xlabel("surface samples")
    axes[0].set_ylabel("relative error")
    axes[0].set_title("surface quadrature")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)

    x_voxel = [as_float(r, "voxel_size") * 1e6 for r in voxel]
    axes[1].loglog(x_voxel, [as_float(r, "sdf_force_rel_error") for r in voxel], marker="o", lw=1.2, label="SDF force")
    axes[1].loglog(x_voxel, [as_float(r, "sdf_torque_rel_error") for r in voxel], marker="s", lw=1.2, label="SDF torque")
    axes[1].loglog(x_voxel, [as_float(r, "center_error") for r in voxel], marker="^", lw=1.2, label="center (m)")
    axes[1].invert_xaxis()
    axes[1].set_xlabel("voxel size (um)")
    axes[1].set_ylabel("error")
    axes[1].set_title("SDF query in patch")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_gear_physics(summary_rows: list[dict[str, str]], frame_rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    bidir = [r for r in frame_rows if r["mode"] == "bidirectional_distributed"]
    bidir.sort(key=lambda r: as_float(r, "time"))
    time = [as_float(r, "time") for r in bidir]
    torque = [as_float(r, "torque") for r in bidir]
    force = [as_float(r, "force_norm") for r in bidir]
    penetration = [as_float(r, "max_effective_penetration") * 1e6 for r in bidir]
    power = [as_float(r, "contact_power") for r in bidir]
    work = [0.0]
    for idx in range(1, len(time)):
        work.append(work[-1] + 0.5 * (power[idx] + power[idx - 1]) * (time[idx] - time[idx - 1]))

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.35))
    axes = axes.ravel()
    axes[0].plot(time, torque, color="#4c78a8", lw=1.1)
    axes[0].set_title("contact torque")
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("N m")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(time, force, color="#2f6f4e", lw=1.1)
    axes[1].set_title("normal force norm")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("N")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(time, penetration, color="#b7410e", lw=1.1)
    axes[2].set_title("effective penetration")
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("um")
    axes[2].grid(True, alpha=0.25)

    axes[3].plot(time, work, color="#6f4c9b", lw=1.1)
    axes[3].set_title("contact work on GEAR22")
    axes[3].set_xlabel("time (s)")
    axes[3].set_ylabel("J")
    axes[3].grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_profile(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    row = rows[0]
    labels = ["SDF query", "patch/history/wrench", "pair assembly"]
    values = [
        as_float(row, "mean_query_ms"),
        as_float(row, "mean_patch_history_wrench_ms"),
        as_float(row, "mean_pair_assembly_ms"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.55))
    axes[0].bar(labels, values, color=["#4c78a8", "#2f6f4e", "#b7410e"])
    axes[0].set_ylabel("mean time (ms)")
    axes[0].set_title("per-evaluation profile")
    axes[0].tick_params(axis="x", rotation=15)
    axes[0].grid(True, axis="y", alpha=0.25)

    axes[1].bar(["active samples", "patches"], [as_float(row, "mean_active_samples"), as_float(row, "mean_patch_count")],
                color=["#4c78a8", "#2f6f4e"])
    axes[1].set_title("contact work size")
    axes[1].grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def mean_by_baseline(rows: list[dict[str, str]], key: str) -> tuple[list[str], list[float]]:
    values: dict[str, list[float]] = {}
    for row in rows:
        values.setdefault(row["baseline_variant"], []).append(as_float(row, key))
    labels = list(values)
    means = [sum(values[label]) / len(values[label]) for label in labels]
    return labels, means


def plot_nonconvex_baselines(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55))
    specs = [
        ("mean_patch_count_change_ratio", "patch churn ratio"),
        ("rms_normal_jump_ratio", "normal jump ratio"),
        ("torque_oscillation_ratio", "torque osc. ratio"),
    ]
    for ax, (key, title) in zip(axes, specs):
        labels, values = mean_by_baseline(rows, key)
        ax.bar([label.replace("_", "\n") for label in labels], values, color="#4c78a8")
        ax.axhline(1.0, color="0.2", lw=0.8)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_chrono_native(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55))
    labels = [r["scenario"].replace("chrono_original_", "").replace("_", "\n") for r in rows]
    specs = [
        ("native_to_field_count_churn_ratio", "contact churn"),
        ("native_contact_normal_to_field_normal_jump_ratio", "normal jump"),
        ("native_to_field_torque_oscillation_ratio", "torque oscillation"),
    ]
    for ax, (key, title) in zip(axes, specs):
        ax.bar(labels, [as_float(r, key) for r in rows], color="#b7410e")
        ax.axhline(1.0, color="0.2", lw=0.8)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_free_friction_dynamics(rows: list[dict[str, str]], output: Path) -> None:
    configure_matplotlib()
    labels = [r["scenario"].replace("free_", "").replace("_", "\n") for r in rows]
    count_ratios = [
        as_float(r, "native_mean_abs_count_change") / max(as_float(r, "field_mean_abs_count_change"), 1e-12)
        for r in rows
    ]
    accel_ratios = [1.0 / max(as_float(r, "field_to_native_acceleration_jump_ratio"), 1e-12) for r in rows]
    switch_ratios = [
        as_float(r, "native_stick_slip_switches") / max(as_float(r, "field_stick_slip_switches"), 1e-12)
        for r in rows
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55))
    specs = [
        (count_ratios, "count-change ratio"),
        (accel_ratios, "accel-jump ratio"),
        (switch_ratios, "stick-slip switch ratio"),
    ]
    for ax, (values, title) in zip(axes, specs):
        ax.bar(labels, values, color="#4c78a8")
        ax.axhline(1.0, color="0.2", lw=0.8)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    input_dir = root / "out" / "jcnd_validation"
    output_dir = root / "paper" / "jcnd" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_headon(read_csv(input_dir / "headon_validation.csv"), output_dir / "jcnd_headon_validation.pdf")
    plot_spatial_convergence(read_csv(input_dir / "sdf_spatial_convergence.csv"), output_dir / "jcnd_sdf_spatial_convergence.pdf")
    plot_patch_wrench_convergence(
        read_csv(input_dir / "patch_wrench_convergence.csv"),
        output_dir / "jcnd_patch_wrench_convergence.pdf",
    )
    plot_profile(read_csv(input_dir / "profile_summary.csv"), output_dir / "jcnd_profile_summary.pdf")
    gear_summary = read_csv(input_dir / "gear_ablation_summary.csv")
    gear_frames = read_csv(input_dir / "gear_ablation_frames.csv")
    plot_gear(gear_summary, gear_frames, output_dir / "jcnd_gear_ablation.pdf")
    plot_gear_physics(gear_summary, gear_frames, output_dir / "jcnd_gear_physical_diagnostics.pdf")
    plot_gear_sensitivity(
        read_csv(input_dir / "gear_sensitivity_summary.csv"),
        output_dir / "jcnd_gear_sensitivity.pdf",
    )
    plot_nonconvex_baselines(
        read_csv(root / "out" / "milestone_24" / "mesh_sdf_nonconvex_comparison.csv"),
        output_dir / "jcnd_nonconvex_sdf_baselines.pdf",
    )
    plot_chrono_native(
        read_csv(root / "out" / "milestone_31" / "chrono_original_mesh_comparison.csv"),
        output_dir / "jcnd_chrono_native_baseline.pdf",
    )
    plot_free_friction_dynamics(
        read_csv(root / "out" / "milestone_28" / "free_dynamics_comparison.csv"),
        output_dir / "jcnd_free_friction_dynamics.pdf",
    )


if __name__ == "__main__":
    main()
