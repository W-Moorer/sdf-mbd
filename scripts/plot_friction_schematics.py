#!/usr/bin/env python3
"""Draw vector schematics for the tangential friction validation cases."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.transforms import Affine2D


BLUE = "#1F4EFA"
RED = "#E64B35"
GREEN = "#2CA02C"
DARK = "#3A3A3A"
GRAY = "#7A7A7A"
LIGHT = "#F3F5F7"


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def arrow(ax, start, end, color=DARK, lw=0.85, ms=8, ls="-", zorder=4) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "mutation_scale": ms,
            "lw": lw,
            "color": color,
            "linestyle": ls,
            "shrinkA": 0,
            "shrinkB": 0,
        },
        zorder=zorder,
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        0.02,
        0.97,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
    )


def setup_axis(ax) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.75)
        spine.set_color("0.25")


def draw_incline(ax) -> None:
    setup_axis(ax)
    panel_label(ax, "(a) Inclined plane")

    theta = math.radians(24)
    p0 = (0.12, 0.25)
    p1 = (0.92, 0.25 + 0.80 * math.tan(theta))
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=DARK, lw=0.9)
    ax.plot([p0[0], 0.92], [p0[1], p0[1]], color="0.75", lw=0.55)
    ax.add_patch(
        patches.Arc(
            p0,
            0.24,
            0.24,
            theta1=0,
            theta2=math.degrees(theta),
            color=GRAY,
            lw=0.65,
        )
    )
    ax.text(0.25, 0.29, r"$\theta$", fontsize=8)

    center = (0.50, 0.25 + (0.50 - 0.12) * math.tan(theta) + 0.08)
    block_w, block_h = 0.22, 0.12
    rect = patches.Rectangle(
        (center[0] - block_w / 2, center[1] - block_h / 2),
        block_w,
        block_h,
        facecolor=LIGHT,
        edgecolor=DARK,
        lw=0.85,
    )
    rect.set_transform(Affine2D().rotate_around(*center, theta) + ax.transData)
    ax.add_patch(rect)

    t = (math.cos(theta), math.sin(theta))
    n = (-math.sin(theta), math.cos(theta))
    c = center
    arrow(ax, c, (c[0] + 0.22 * n[0], c[1] + 0.22 * n[1]), BLUE)
    ax.text(
        c[0] + 0.12 * n[0] + 0.10 * t[0],
        c[1] + 0.12 * n[1] + 0.10 * t[1],
        r"$N$",
        color=BLUE,
        fontsize=8,
        ha="center",
    )
    arrow(ax, c, (c[0], c[1] - 0.24), DARK)
    ax.text(c[0] + 0.025, c[1] - 0.22, r"$mg$", fontsize=8)
    arrow(ax, (c[0] + 0.04 * t[0], c[1] + 0.04 * t[1]), (c[0] - 0.22 * t[0], c[1] - 0.22 * t[1]), RED)
    ax.text(c[0] - 0.25 * t[0], c[1] - 0.25 * t[1], r"$f_t$", color=RED, fontsize=8, ha="center")
    arrow(ax, (c[0] + 0.05 * t[0], c[1] + 0.05 * t[1]), (c[0] + 0.28 * t[0], c[1] + 0.28 * t[1]), GREEN, lw=0.8)
    ax.text(c[0] + 0.31 * t[0], c[1] + 0.31 * t[1], r"$a$", color=GREEN, fontsize=8)

    ax.text(0.08, 0.82, r"stick: $\tan\theta < \mu$", fontsize=8, color=DARK)
    ax.text(0.08, 0.73, r"slip: $a=g(\sin\theta-\mu\cos\theta)$", fontsize=8, color=DARK)
    ax.text(0.60, 0.18, r"$f_t \leq \mu N$", fontsize=8, color=RED)


def draw_spring(ax) -> None:
    setup_axis(ax)
    panel_label(ax, "(b) Spring pull")

    y0 = 0.28
    ax.plot([0.08, 0.94], [y0, y0], color=DARK, lw=0.85)
    ax.plot([0.08, 0.08], [y0, 0.72], color=DARK, lw=0.85)

    block = patches.Rectangle((0.52, y0), 0.18, 0.16, facecolor=LIGHT, edgecolor=DARK, lw=0.85)
    ax.add_patch(block)

    xs = [0.08, 0.13]
    ys = [0.42, 0.42]
    ncoil = 8
    x_start, x_end = 0.13, 0.52
    for i in range(2 * ncoil + 1):
        x = x_start + (x_end - x_start) * i / (2 * ncoil)
        y = 0.42 + (0.035 if i % 2 else -0.035)
        xs.append(x)
        ys.append(y)
    xs.append(0.52)
    ys.append(0.42)
    ax.plot(xs, ys, color=DARK, lw=0.85)
    ax.text(0.27, 0.50, r"$k$", fontsize=8)

    arrow(ax, (0.70, 0.39), (0.89, 0.39), BLUE)
    ax.text(0.79, 0.44, r"$x_d=Vt$", fontsize=8, color=BLUE, ha="center")
    arrow(ax, (0.58, 0.28), (0.42, 0.28), RED)
    ax.text(0.43, 0.20, r"$f_t$", fontsize=8, color=RED)
    arrow(ax, (0.66, 0.28), (0.85, 0.28), GREEN)
    ax.text(0.82, 0.20, r"$v$", fontsize=8, color=GREEN)

    ax.plot([0.15, 0.85], [0.77, 0.77], color="0.65", lw=0.65)
    ax.plot([0.48, 0.48], [0.735, 0.805], color=RED, lw=0.85)
    ax.text(0.16, 0.83, r"stick: $kx_d < \mu N$", fontsize=8)
    ax.text(0.50, 0.83, r"slip: $kx_d = \mu N$", fontsize=8)
    ax.text(0.43, 0.67, r"$t_s=\mu N/(kV)$", fontsize=8, color=RED, ha="center")


def draw_transport(ax) -> None:
    setup_axis(ax)
    panel_label(ax, "(c) Objective transport")

    origin = (0.42, 0.36)
    n0_angle = math.radians(78)
    n1_angle = math.radians(42)
    n0 = (math.cos(n0_angle), math.sin(n0_angle))
    n1 = (math.cos(n1_angle), math.sin(n1_angle))
    t0 = (-n0[1], n0[0])
    t1 = (-n1[1], n1[0])

    ax.add_patch(patches.Circle(origin, 0.055, facecolor=LIGHT, edgecolor=DARK, lw=0.75))
    for vec, color, label, off in [
        (n0, GRAY, r"$n_0$", (0.01, 0.02)),
        (n1, BLUE, r"$n_1$", (0.02, -0.02)),
    ]:
        end = (origin[0] + 0.33 * vec[0], origin[1] + 0.33 * vec[1])
        arrow(ax, origin, end, color=color)
        ax.text(end[0] + off[0], end[1] + off[1], label, fontsize=8, color=color)

    for vec, color in [(t0, "0.70"), (t1, "0.55")]:
        ax.plot(
            [origin[0] - 0.25 * vec[0], origin[0] + 0.25 * vec[0]],
            [origin[1] - 0.25 * vec[1], origin[1] + 0.25 * vec[1]],
            color=color,
            lw=0.65,
            ls="--",
        )

    hist0_start = (origin[0] - 0.02 * n0[0], origin[1] - 0.02 * n0[1])
    hist0_end = (hist0_start[0] + 0.22 * t0[0], hist0_start[1] + 0.22 * t0[1])
    hist1_start = (origin[0] - 0.02 * n1[0], origin[1] - 0.02 * n1[1])
    hist1_end = (hist1_start[0] + 0.22 * t1[0], hist1_start[1] + 0.22 * t1[1])
    arrow(ax, hist0_start, hist0_end, RED, lw=0.85)
    arrow(ax, hist1_start, hist1_end, GREEN, lw=0.85)
    ax.text(hist0_end[0] - 0.03, hist0_end[1] + 0.03, r"$\xi_t^0$", fontsize=8, color=RED)
    ax.text(hist1_end[0] + 0.01, hist1_end[1] + 0.02, r"$\xi_t^1$", fontsize=8, color=GREEN)

    arc = patches.Arc(origin, 0.43, 0.43, theta1=42, theta2=78, color=BLUE, lw=0.75)
    ax.add_patch(arc)
    ax.text(0.61, 0.70, r"$\Delta R$", fontsize=8, color=BLUE)
    ax.plot(
        [hist0_end[0], origin[0] + 0.22 * t1[0]],
        [hist0_end[1], origin[1] + 0.22 * t1[1]],
        color=RED,
        lw=0.75,
        ls=(0, (3, 2)),
    )
    ax.text(0.08, 0.83, "minimal-rotation transport", fontsize=8, color=GREEN)
    ax.text(0.08, 0.74, "naive projection changes history", fontsize=8, color=RED)


def make_figure(output_dir: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.35))
    draw_incline(axes[0])
    draw_spring(axes[1])
    draw_transport(axes[2])
    fig.tight_layout(pad=0.25, w_pad=0.55)
    fig.savefig(output_dir / "friction_validation_schematics.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="paper/jcnd/figures")
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    make_figure(output_dir)


if __name__ == "__main__":
    main()
