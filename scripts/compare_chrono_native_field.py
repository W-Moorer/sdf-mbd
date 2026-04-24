#!/usr/bin/env python3
"""Compare Chrono native SMC baseline against field primitive CSV outputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import generate_field_contact_figures as figures


figures.COLORS["chrono_native_smc"] = "#b91c1c"
figures.COLORS["native_command"] = "#6b7280"
figures.COLORS["native_actual"] = "#b91c1c"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def i(row: dict[str, str], key: str, default: int = 0) -> int:
    try:
        return int(float(row.get(key, default)))
    except (TypeError, ValueError):
        return default


def v(row: dict[str, str], prefix: str) -> tuple[float, float, float]:
    return (f(row, f"{prefix}_x"), f(row, f"{prefix}_y"), f(row, f"{prefix}_z"))


def norm(vec: tuple[float, float, float]) -> float:
    return math.sqrt(sum(x * x for x in vec))


def diff(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def group_by_scenario(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        out[row["scenario"]].append(row)
    for scenario_rows in out.values():
        scenario_rows.sort(key=lambda row: i(row, "frame"))
    return out


def field_rows(frames: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    return group_by_scenario([row for row in frames if row.get("variant") == "field_primitive"])


def vector_jump(rows: list[dict[str, str]], prefix: str) -> list[float]:
    jumps = []
    previous = None
    for row in rows:
        current = v(row, prefix)
        jumps.append(0.0 if previous is None else norm(diff(current, previous)))
        previous = current
    return jumps


def mean_abs_count_change(rows: list[dict[str, str]], key: str) -> float:
    if len(rows) < 2:
        return 0.0
    total = 0.0
    for a, b in zip(rows[:-1], rows[1:]):
        total += abs(i(b, key) - i(a, key))
    return total / float(len(rows) - 1)


def oscillation_index(rows: list[dict[str, str]], prefix: str) -> tuple[float, float, float]:
    jumps = vector_jump(rows, prefix)
    rms = math.sqrt(sum(x * x for x in jumps) / float(len(jumps))) if jumps else 0.0
    max_norm = max((norm(v(row, prefix)) for row in rows), default=0.0)
    index = rms / max_norm if max_norm > 1.0e-14 else 0.0
    return rms, max_norm, index


def finite_diff_acceleration(rows: list[dict[str, str]], vel_prefix: str) -> list[tuple[float, tuple[float, float, float]]]:
    out = []
    for idx, row in enumerate(rows):
        time = f(row, "time")
        if idx == 0 or idx + 1 >= len(rows):
            out.append((time, (0.0, 0.0, 0.0)))
            continue
        dt = f(rows[idx + 1], "time") - f(rows[idx - 1], "time")
        if abs(dt) < 1.0e-14:
            out.append((time, (0.0, 0.0, 0.0)))
            continue
        vp = v(rows[idx + 1], vel_prefix)
        vm = v(rows[idx - 1], vel_prefix)
        out.append((time, ((vp[0] - vm[0]) / dt, (vp[1] - vm[1]) / dt, (vp[2] - vm[2]) / dt)))
    return out


def points_norm(rows: list[dict[str, str]], prefix: str) -> list[tuple[float, float]]:
    return [(f(row, "time"), norm(v(row, prefix))) for row in rows]


def points_scalar(rows: list[dict[str, str]], key: str) -> list[tuple[float, float]]:
    return [(f(row, "time"), f(row, key)) for row in rows]


def make_comparison(native: dict[str, list[dict[str, str]]],
                    field: dict[str, list[dict[str, str]]],
                    output: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for scenario, native_rows in native.items():
        field_scenario = field.get(scenario, [])
        if not field_scenario:
            continue

        field_force_rms, field_force_max, field_force_osc = oscillation_index(field_scenario, "total_force")
        field_torque_rms, field_torque_max, field_torque_osc = oscillation_index(field_scenario, "total_torque")
        native_force_rms, native_force_max, native_force_osc = oscillation_index(native_rows, "contact_force")
        native_torque_rms, native_torque_max, native_torque_osc = oscillation_index(native_rows, "contact_torque")

        field_count_change = mean_abs_count_change(field_scenario, "patch_count")
        native_count_change = mean_abs_count_change(native_rows, "contact_count")

        row = {
            "scenario": scenario,
            "field_force_oscillation_index": field_force_osc,
            "native_force_oscillation_index": native_force_osc,
            "native_to_field_force_oscillation_ratio": native_force_osc / field_force_osc if field_force_osc > 1.0e-14 else 0.0,
            "field_rms_force_jump": field_force_rms,
            "native_rms_force_jump": native_force_rms,
            "native_to_field_rms_force_jump_ratio": native_force_rms / field_force_rms if field_force_rms > 1.0e-14 else 0.0,
            "field_torque_oscillation_index": field_torque_osc,
            "native_torque_oscillation_index": native_torque_osc,
            "native_to_field_torque_oscillation_ratio": native_torque_osc / field_torque_osc if field_torque_osc > 1.0e-14 else 0.0,
            "field_rms_torque_jump": field_torque_rms,
            "native_rms_torque_jump": native_torque_rms,
            "native_to_field_rms_torque_jump_ratio": native_torque_rms / field_torque_rms if field_torque_rms > 1.0e-14 else 0.0,
            "field_max_force_norm": field_force_max,
            "native_max_force_norm": native_force_max,
            "field_max_torque_norm": field_torque_max,
            "native_max_torque_norm": native_torque_max,
            "field_mean_patch_count_change": field_count_change,
            "native_mean_contact_count_change": native_count_change,
            "native_to_field_count_change_ratio": native_count_change / field_count_change if field_count_change > 1.0e-14 else 0.0,
            "field_max_patch_count": max((i(row, "patch_count") for row in field_scenario), default=0),
            "native_max_contact_count": max((i(row, "contact_count") for row in native_rows), default=0),
            "native_active_frames": sum(1 for row in native_rows if i(row, "contact_count") > 0),
            "field_active_frames": sum(1 for row in field_scenario if i(row, "patch_count") > 0),
            "native_max_position_error": max((f(row, "pos_error_norm") for row in native_rows), default=0.0),
            "native_rms_position_error": math.sqrt(sum(f(row, "pos_error_norm") ** 2 for row in native_rows) / float(len(native_rows))),
            "native_max_velocity_error": max((f(row, "vel_error_norm") for row in native_rows), default=0.0),
            "native_rms_velocity_error": math.sqrt(sum(f(row, "vel_error_norm") ** 2 for row in native_rows) / float(len(native_rows))),
            "native_max_acceleration_error": max((f(row, "acc_error_norm") for row in native_rows), default=0.0),
            "native_rms_acceleration_error": math.sqrt(sum(f(row, "acc_error_norm") ** 2 for row in native_rows) / float(len(native_rows))),
        }
        rows.append(row)

    keys = list(rows[0].keys()) if rows else []
    with output.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def generate_figures(native: dict[str, list[dict[str, str]]],
                     field: dict[str, list[dict[str, str]]],
                     output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for scenario, native_rows in native.items():
        field_rows_for_scenario = field.get(scenario, [])
        if not field_rows_for_scenario:
            continue

        base = figures.sanitize_name(scenario)
        series_force = [
            ("field_primitive", [(f(row, "time"), jump) for row, jump in zip(field_rows_for_scenario, vector_jump(field_rows_for_scenario, "total_force"))]),
            ("chrono_native_smc", [(f(row, "time"), jump) for row, jump in zip(native_rows, vector_jump(native_rows, "contact_force"))]),
        ]
        path = output_dir / f"{base}_force_jump_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Force Jump", "||Delta F||", series_force, y_min=0.0)
        written.append(path)

        series_torque = [
            ("field_primitive", [(f(row, "time"), jump) for row, jump in zip(field_rows_for_scenario, vector_jump(field_rows_for_scenario, "total_torque"))]),
            ("chrono_native_smc", [(f(row, "time"), jump) for row, jump in zip(native_rows, vector_jump(native_rows, "contact_torque"))]),
        ]
        path = output_dir / f"{base}_torque_jump_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Torque Jump", "||Delta tau||", series_torque, y_min=0.0)
        written.append(path)

        series_count = [
            ("field_primitive", points_scalar(field_rows_for_scenario, "patch_count")),
            ("chrono_native_smc", points_scalar(native_rows, "contact_count")),
        ]
        path = output_dir / f"{base}_count_timeline_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Contact Count Timeline", "count", series_count, y_min=0.0)
        written.append(path)

        series_pos = [
            ("field_primitive", points_norm(field_rows_for_scenario, "center")),
            ("native_command", points_norm(native_rows, "cmd_pos")),
            ("native_actual", points_norm(native_rows, "actual_pos")),
        ]
        path = output_dir / f"{base}_displacement_norm_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Displacement Norm", "||x||", series_pos, y_min=0.0)
        written.append(path)

        series_vel = [
            ("field_primitive", points_norm(field_rows_for_scenario, "velocity")),
            ("native_command", points_norm(native_rows, "cmd_vel")),
            ("native_actual", points_norm(native_rows, "actual_vel")),
        ]
        path = output_dir / f"{base}_velocity_norm_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Velocity Norm", "||v||", series_vel, y_min=0.0)
        written.append(path)

        field_acc = [(time, norm(acc)) for time, acc in finite_diff_acceleration(field_rows_for_scenario, "velocity")]
        native_cmd_acc = points_norm(native_rows, "cmd_acc")
        native_actual_acc = points_norm(native_rows, "actual_acc")
        series_acc = [
            ("field_primitive", field_acc),
            ("native_command", native_cmd_acc),
            ("native_actual", native_actual_acc),
        ]
        path = output_dir / f"{base}_acceleration_norm_field_vs_native.svg"
        figures.draw_line_chart(path, f"{scenario}: Acceleration Norm", "||a||", series_acc, y_min=0.0)
        written.append(path)

    figures.write_index(output_dir, written)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=Path("out/milestone_23/chrono_native_baseline_frames.csv"))
    parser.add_argument("--field", type=Path, default=Path("out/milestone_21/field_contact_multislip_interlock_frames.csv"))
    parser.add_argument("--output", type=Path, default=Path("out/milestone_23"))
    args = parser.parse_args()

    native = group_by_scenario(read_csv(args.native))
    field = field_rows(read_csv(args.field))
    args.output.mkdir(parents=True, exist_ok=True)

    comparison_path = args.output / "chrono_native_field_comparison.csv"
    rows = make_comparison(native, field, comparison_path)
    figure_paths = generate_figures(native, field, args.output / "figures")

    print(comparison_path)
    print(f"Generated {len(figure_paths)} SVG figures")
    for row in rows:
        print(
            f"{row['scenario']}: force_ratio={row['native_to_field_force_oscillation_ratio']:.6g} "
            f"torque_ratio={row['native_to_field_torque_oscillation_ratio']:.6g} "
            f"count_ratio={row['native_to_field_count_change_ratio']:.6g} "
            f"max_pos_err={row['native_max_position_error']:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
