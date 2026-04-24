#!/usr/bin/env python3
"""Plot field-contact friction validation cases as paper-ready PDFs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def set_times_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.figsize": (6.6, 2.6),
        }
    )


def plot_inclined(input_dir: Path, output_dir: Path) -> None:
    rows = read_csv(input_dir / "inclined_plane.csv")
    stick = [r for r in rows if r["case"] == "stick_15deg"]
    slip = [r for r in rows if r["case"] == "slip_30deg"]

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5))

    axes[0].plot([f(r, "time") for r in stick], [1e3 * f(r, "s") for r in stick], lw=1.2)
    axes[0].axhline(0.0, color="0.3", lw=0.7, ls="--")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Down-slope displacement (mm)")
    axes[0].set_title("15 deg, tan(theta) < mu")

    axes[1].plot([f(r, "time") for r in slip], [f(r, "acceleration") for r in slip], lw=1.0, label="Computed")
    axes[1].plot(
        [f(r, "time") for r in slip],
        [f(r, "expected_acceleration") for r in slip],
        lw=0.9,
        ls="--",
        label="Analytic",
    )
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Down-slope acceleration (m/s^2)")
    axes[1].set_title("30 deg, tan(theta) > mu")
    axes[1].legend(frameon=False)

    fig.tight_layout(pad=0.35, w_pad=1.0)
    fig.savefig(output_dir / "friction_inclined_plane.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_spring_pull(input_dir: Path, output_dir: Path) -> None:
    rows = read_csv(input_dir / "spring_pull.csv")
    t = [f(r, "time") for r in rows]
    drive = [f(r, "drive_force") for r in rows]
    contact = [-f(r, "contact_x") for r in rows]
    limit = [0.4 * f(r, "normal_force") for r in rows]
    slip_times = [f(r, "time") for r in rows if r["state"] == "slip"]
    first_slip = slip_times[0] if slip_times else None

    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    ax.plot(t, drive, lw=1.1, label="Drive force")
    ax.plot(t, contact, lw=1.1, label="Contact resistance")
    ax.plot(t, limit, lw=0.9, ls="--", color="0.2", label="mu N")
    if first_slip is not None:
        ax.axvline(first_slip, color="0.3", lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Force (N)")
    ax.legend(frameon=False)
    fig.tight_layout(pad=0.35)
    fig.savefig(output_dir / "friction_spring_pull.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_transport(input_dir: Path, output_dir: Path) -> None:
    rows = read_csv(input_dir / "objective_transport.csv")
    angle = [f(r, "angle_deg") for r in rows]
    minimal_norm_error = [f(r, "minimal_norm_error") for r in rows]
    untransported_tangent = [f(r, "untransported_tangent_error") for r in rows]
    projected_energy = [f(r, "projected_naive_energy_ratio") for r in rows]

    fig, ax1 = plt.subplots(figsize=(4.2, 2.6))
    ax1.semilogy(angle, minimal_norm_error, lw=1.1, label="Minimal-rotation norm error")
    ax1.set_xlabel("Normal rotation angle (deg)")
    ax1.set_ylabel("Relative error")
    ax1.set_ylim(1e-17, 2.0)

    ax2 = ax1.twinx()
    ax2.plot(angle, untransported_tangent, lw=1.0, ls="--", color="0.25", label="Untransported normal component")
    ax2.plot(angle, projected_energy, lw=1.0, ls="-.", color="0.55", label="Projected-naive energy ratio")
    ax2.set_ylabel("Ratio")
    ax2.set_ylim(0.0, 1.05)

    lines = ax1.get_lines() + ax2.get_lines()
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, frameon=False, loc="lower right")
    fig.tight_layout(pad=0.35)
    fig.savefig(output_dir / "friction_objective_transport.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="out/friction_validation")
    parser.add_argument("--output-dir", default="paper/figures")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    set_times_style()
    plot_inclined(input_dir, output_dir)
    plot_spring_pull(input_dir, output_dir)
    plot_transport(input_dir, output_dir)


if __name__ == "__main__":
    main()
