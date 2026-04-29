#!/usr/bin/env python3
"""Compare SDF-NCP benchmark trajectories against RecurDyn reference CSV data.

The comparison code is intentionally separated from the C++ solver.  Case-specific
configuration is limited to reference column selection and output body selection;
the interpolation, metrics, and plotting are shared.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "sdf_ncp_benchmarks"


@dataclass(frozen=True)
class QuantityConfig:
    name: str
    reference_column: int
    trajectory_column: str
    ylabel: str


@dataclass(frozen=True)
class CompareConfig:
    case_name: str
    result_case: str
    reference_path: Path
    trajectory_path: Path
    body_id: str
    time_column: int
    quantities: tuple[QuantityConfig, ...]


def cam_config(result_case: str) -> CompareConfig:
    return CompareConfig(
        case_name="cam",
        result_case=result_case,
        reference_path=ROOT / "assets" / "cam" / "data" / "cam_data.csv",
        trajectory_path=RESULT_DIR / result_case / "trajectory.csv",
        body_id="follower",
        time_column=0,
        quantities=(
            QuantityConfig("follower_y", 10, "py", "follower y (m)"),
            QuantityConfig("follower_vy", 13, "vy", "follower y velocity (m/s)"),
        ),
    )


def read_reference(config: CompareConfig) -> list[dict[str, float]]:
    if not config.reference_path.exists():
        raise FileNotFoundError(f"Missing RecurDyn reference CSV: {config.reference_path}")
    out: list[dict[str, float]] = []
    with config.reference_path.open(newline="", encoding="utf-8-sig", errors="replace") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if not row:
                continue
            values: dict[str, float] = {"time": float(row[config.time_column])}
            for quantity in config.quantities:
                values[quantity.name] = float(row[quantity.reference_column])
            out.append(values)
    out.sort(key=lambda row: row["time"])
    return out


def read_trajectory(config: CompareConfig) -> list[dict[str, float]]:
    if not config.trajectory_path.exists():
        raise FileNotFoundError(
            f"Missing SDF-NCP trajectory CSV: {config.trajectory_path}. "
            "Run the matching C++ benchmark first."
        )
    rows: list[dict[str, float]] = []
    seen_times: set[float] = set()
    with config.trajectory_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("body_id") != config.body_id:
                continue
            time = float(row["time"])
            # Each contact sample can duplicate the same body state. Keep one state
            # per time for trajectory comparison.
            if time in seen_times:
                continue
            seen_times.add(time)
            values = {"time": time}
            for quantity in config.quantities:
                values[quantity.name] = float(row[quantity.trajectory_column])
            rows.append(values)
    rows.sort(key=lambda row: row["time"])
    return rows


def interpolate(rows: list[dict[str, float]], key: str, time: float) -> float:
    if not rows:
        raise ValueError("Cannot interpolate an empty trajectory.")
    if time <= rows[0]["time"]:
        return rows[0][key]
    if time >= rows[-1]["time"]:
        return rows[-1][key]
    lo = 0
    hi = len(rows) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rows[mid]["time"] <= time:
            lo = mid
        else:
            hi = mid
    t0 = rows[lo]["time"]
    t1 = rows[hi]["time"]
    if abs(t1 - t0) < 1e-15:
        return rows[lo][key]
    alpha = (time - t0) / (t1 - t0)
    return (1.0 - alpha) * rows[lo][key] + alpha * rows[hi][key]


def compare(config: CompareConfig) -> tuple[list[dict[str, float]], list[dict[str, str]]]:
    reference = read_reference(config)
    trajectory = read_trajectory(config)
    if not reference:
        raise ValueError("RecurDyn reference CSV is empty.")
    if not trajectory:
        raise ValueError("SDF-NCP trajectory CSV has no rows for the configured body.")

    t_min = max(reference[0]["time"], trajectory[0]["time"])
    t_max = min(reference[-1]["time"], trajectory[-1]["time"])
    aligned: list[dict[str, float]] = []
    for ref in reference:
        time = ref["time"]
        if time < t_min - 1e-14 or time > t_max + 1e-14:
            continue
        row: dict[str, float] = {"time": time}
        for quantity in config.quantities:
            ref_value = ref[quantity.name]
            sdf_value = interpolate(trajectory, quantity.name, time)
            row[f"{quantity.name}_recurdyn"] = ref_value
            row[f"{quantity.name}_sdf_ncp"] = sdf_value
            row[f"{quantity.name}_error"] = sdf_value - ref_value
        aligned.append(row)

    summary: list[dict[str, str]] = []
    for quantity in config.quantities:
        errors = [row[f"{quantity.name}_error"] for row in aligned]
        refs = [row[f"{quantity.name}_recurdyn"] for row in aligned]
        sdfs = [row[f"{quantity.name}_sdf_ncp"] for row in aligned]
        rmse = math.sqrt(sum(e * e for e in errors) / max(1, len(errors)))
        max_abs = max((abs(e) for e in errors), default=0.0)
        final_ref = refs[-1] if refs else float("nan")
        final_sdf = sdfs[-1] if sdfs else float("nan")
        final_abs = abs(final_sdf - final_ref)
        scale = max(max((abs(v) for v in refs), default=0.0), 1e-14)
        summary.append(
            {
                "case_name": config.case_name,
                "result_case": config.result_case,
                "quantity": quantity.name,
                "time_start": f"{t_min:.17g}",
                "time_end": f"{t_max:.17g}",
                "num_samples": str(len(aligned)),
                "rmse": f"{rmse:.17g}",
                "max_abs_error": f"{max_abs:.17g}",
                "final_recurdyn": f"{final_ref:.17g}",
                "final_sdf_ncp": f"{final_sdf:.17g}",
                "final_abs_error": f"{final_abs:.17g}",
                "relative_rmse_reference_scale": f"{rmse / scale:.17g}",
            }
        )
    return aligned, summary


def write_csvs(out_dir: Path, aligned: list[dict[str, float]], summary: list[dict[str, str]]) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    timeseries_path = out_dir / "recurdyn_reference_comparison_timeseries.csv"
    summary_path = out_dir / "recurdyn_reference_comparison_summary.csv"

    if aligned:
        with timeseries_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(aligned[0].keys()))
            writer.writeheader()
            writer.writerows(aligned)
    else:
        timeseries_path.write_text("time\n", encoding="utf-8")

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "case_name",
            "result_case",
            "quantity",
            "time_start",
            "time_end",
            "num_samples",
            "rmse",
            "max_abs_error",
            "final_recurdyn",
            "final_sdf_ncp",
            "final_abs_error",
            "relative_rmse_reference_scale",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary)
    return timeseries_path, summary_path


def plot_comparison(config: CompareConfig, out_dir: Path, aligned: list[dict[str, float]]) -> list[Path]:
    fig_dir = out_dir / "figures" / "reference_comparison"
    fig_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    if not aligned:
        return outputs
    times = [row["time"] for row in aligned]

    for quantity in config.quantities:
        path = fig_dir / f"{quantity.name}_overlay.png"
        fig, ax = plt.subplots(figsize=(9.0, 5.0))
        ax.plot(
            times,
            [row[f"{quantity.name}_recurdyn"] for row in aligned],
            color="#111111",
            linewidth=2.2,
            label="RecurDyn reference",
        )
        ax.plot(
            times,
            [row[f"{quantity.name}_sdf_ncp"] for row in aligned],
            color="#d95f02",
            linewidth=1.8,
            linestyle="--",
            label="SDF-NCP benchmark",
        )
        ax.set_xlabel("time (s)")
        ax.set_ylabel(quantity.ylabel)
        ax.set_title(f"{config.case_name}: {quantity.name}")
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        outputs.append(path)

        path = fig_dir / f"{quantity.name}_error.png"
        fig, ax = plt.subplots(figsize=(9.0, 4.6))
        ax.axhline(0.0, color="#111111", linewidth=1.0)
        ax.plot(
            times,
            [row[f"{quantity.name}_error"] for row in aligned],
            color="#1b9e77",
            linewidth=1.8,
            label="SDF-NCP - RecurDyn",
        )
        ax.set_xlabel("time (s)")
        ax.set_ylabel(f"{quantity.ylabel} error")
        ax.set_title(f"{config.case_name}: {quantity.name} error")
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        outputs.append(path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="cam", choices=["cam"])
    parser.add_argument(
        "--result-case",
        default="cam_recurdyn_compare",
        help="Subdirectory under results/sdf_ncp_benchmarks containing trajectory.csv.",
    )
    args = parser.parse_args()

    if args.case == "cam":
        config = cam_config(args.result_case)
    else:
        raise ValueError(args.case)

    out_dir = RESULT_DIR / args.result_case
    aligned, summary = compare(config)
    csv_paths = write_csvs(out_dir, aligned, summary)
    figure_paths = plot_comparison(config, out_dir, aligned)

    for path in csv_paths + tuple(figure_paths):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
