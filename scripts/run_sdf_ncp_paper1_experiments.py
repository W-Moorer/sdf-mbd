#!/usr/bin/env python3
"""Run the first-paper SDF-NCP reproducibility experiments."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CPP_DEMO = ROOT / "build" / "bin" / "Release" / "demo_CH_sdf_ncp_regression.exe"

GENERATED_FILES = [
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "point_mass_plane.csv",
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "height_vs_time.png",
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "gap_vs_time.png",
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "penetration_vs_time.png",
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "contact_force_vs_time.png",
    ROOT / "results" / "sdf_ncp" / "point_mass_plane" / "complementarity_error_vs_time.png",
    ROOT / "results" / "sdf_ncp" / "penalty_sensitivity" / "summary.csv",
    ROOT / "results" / "sdf_ncp" / "penalty_sensitivity" / "max_penetration_vs_parameter.png",
    ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity" / "summary.csv",
    ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity" / "max_penetration_vs_eps.png",
    ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity" / "complementarity_vs_eps.png",
    ROOT / "results" / "sdf_ncp" / "epsilon_sensitivity" / "iterations_vs_eps.png",
    ROOT / "results" / "sdf_ncp" / "geometry" / "sdf_contours_normals.png",
    ROOT / "results" / "sdf_ncp_cpp" / "rigidbody2d_rollout.csv",
    ROOT / "results" / "sdf_ncp_cpp" / "figures" / "rigidbody2d_pose_vs_time.png",
    ROOT / "results" / "sdf_ncp_cpp" / "figures" / "rigidbody2d_gaps_vs_time.png",
    ROOT / "results" / "sdf_ncp_cpp" / "figures" / "rigidbody2d_lambdas_vs_time.png",
    ROOT / "results" / "sdf_ncp_cpp" / "figures" / "rigidbody2d_complementarity_vs_time.png",
    ROOT / "results" / "sdf_ncp_paper1" / "tables" / "method_summary.csv",
]


def run_step(name: str, command: list[str]) -> None:
    print(f"[sdf-ncp-paper1] Running {name}", flush=True)
    completed = subprocess.run(command, cwd=ROOT)
    if completed.returncode != 0:
        print(
            f"[sdf-ncp-paper1] Step failed: {name}, return code {completed.returncode}",
            file=sys.stderr,
        )
        raise SystemExit(completed.returncode)


def main() -> None:
    if not CPP_DEMO.exists():
        print(
            f"C++ demo not found: {CPP_DEMO}\n"
            "Build it first with:\n"
            r"  cmake --build build --config Release --target demo_CH_sdf_ncp_regression",
            file=sys.stderr,
        )
        raise SystemExit(1)

    python = sys.executable
    steps = [
        ("point_mass_plane", [python, "examples/sdf_ncp/point_mass_plane.py"]),
        ("penalty_sensitivity", [python, "examples/sdf_ncp/penalty_sensitivity.py"]),
        ("epsilon_sensitivity", [python, "examples/sdf_ncp/epsilon_sensitivity.py"]),
        ("sdf_geometry_visualization", [python, "examples/sdf_ncp/sdf_geometry_visualization.py"]),
        ("rigidbody2d_export", [str(CPP_DEMO), "rigidbody2d_export"]),
        ("plot_cpp_rigidbody2d_rollout", [python, "examples/sdf_ncp/plot_cpp_rigidbody2d_rollout.py"]),
        ("make_sdf_ncp_paper1_summary", [python, "scripts/make_sdf_ncp_paper1_summary.py"]),
    ]

    for name, command in steps:
        run_step(name, command)

    print("[sdf-ncp-paper1] Generated files:")
    for path in GENERATED_FILES:
        if path.exists():
            print(f"  {path}")
        else:
            print(f"  MISSING: {path}")


if __name__ == "__main__":
    main()
