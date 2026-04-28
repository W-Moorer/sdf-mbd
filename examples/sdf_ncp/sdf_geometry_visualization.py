#!/usr/bin/env python3
"""Visualize analytical SDF contours and normals for method figures."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sdf_mbd.contact import sdf_gap, sdf_normal
from sdf_mbd.sdf import CircleSDF, PlaneSDF


RESULT_DIR = ROOT / "results" / "sdf_ncp" / "geometry"
FIGURE_PATH = RESULT_DIR / "sdf_contours_normals.png"


def plot_sdf(ax, sdf, title: str, example_point: np.ndarray) -> None:
    x = np.linspace(-1.2, 1.2, 161)
    y = np.linspace(-0.5, 1.4, 161)
    xx, yy = np.meshgrid(x, y)
    points = np.stack([xx, yy], axis=-1)
    values = sdf.phi(points)

    contour = ax.contour(xx, yy, values, levels=np.linspace(-0.4, 0.8, 7), colors="#64748b", linewidths=0.8)
    ax.clabel(contour, inline=True, fontsize=7, fmt="%.2f")
    ax.contour(xx, yy, values, levels=[0.0], colors="#dc2626", linewidths=2.0)

    sample_points = []
    for px in np.linspace(-0.9, 0.9, 7):
        for py in np.linspace(-0.2, 1.1, 6):
            p = np.array([px, py])
            if abs(sdf_gap(sdf, p)) < 0.18:
                sample_points.append(p)
    if sample_points:
        pts = np.array(sample_points)
        normals = np.array([sdf_normal(sdf, p) for p in pts])
        ax.quiver(pts[:, 0], pts[:, 1], normals[:, 0], normals[:, 1], color="#0f766e", angles="xy", scale=18)

    gap = sdf_gap(sdf, example_point)
    normal = sdf_normal(sdf, example_point)
    ax.scatter([example_point[0]], [example_point[1]], color="#111827", s=32, zorder=5)
    ax.arrow(
        example_point[0],
        example_point[1],
        0.22 * normal[0],
        0.22 * normal[1],
        width=0.006,
        color="#111827",
        length_includes_head=True,
    )
    ax.text(example_point[0] + 0.04, example_point[1] + 0.04, f"phi={gap:.3f}", fontsize=9)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-0.5, 1.4)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.18)


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    plane = PlaneSDF([0.0, 1.0], 0.0)
    circle = CircleSDF([0.0, 0.55], 0.35)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    plot_sdf(axes[0], plane, "Ground plane SDF", np.array([0.35, 0.24]))
    plot_sdf(axes[1], circle, "Circle SDF", np.array([0.42, 0.78]))
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=200)
    plt.close(fig)
    print(f"Wrote {FIGURE_PATH}")


if __name__ == "__main__":
    main()

