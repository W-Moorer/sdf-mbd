#!/usr/bin/env python3
"""Generate SVG figures for field-contact primitive benchmark CSV files.

The script intentionally uses only the Python standard library so paper figures
can be regenerated on a clean developer machine after running the C++ demos.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
from collections import defaultdict
from pathlib import Path


VARIANT_ORDER = [
    "field_primitive",
    "convex_decomposition_proxy",
    "chrono_traditional_proxy",
    "normal_bin_fragments",
    "raw_sample_contacts",
]

COLORS = {
    "field_primitive": "#0f766e",
    "convex_decomposition_proxy": "#7c3aed",
    "chrono_traditional_proxy": "#dc2626",
    "normal_bin_fragments": "#ea580c",
    "raw_sample_contacts": "#2563eb",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def safe_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def safe_int(row: dict[str, str], key: str, default: int = 0) -> int:
    try:
        return int(float(row.get(key, default)))
    except (TypeError, ValueError):
        return default


def sanitize_name(value: str) -> str:
    out = []
    for ch in value.lower():
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "figure"


def grouped_frames(rows: list[dict[str, str]]) -> dict[str, dict[str, list[dict[str, str]]]]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        grouped[row["scenario"]][row["variant"]].append(row)

    for variants in grouped.values():
        for variant_rows in variants.values():
            variant_rows.sort(key=lambda r: safe_int(r, "frame"))
    return grouped


def ordered_variants(variants: dict[str, list[dict[str, str]]]) -> list[str]:
    known = [v for v in VARIANT_ORDER if v in variants]
    extra = sorted(v for v in variants if v not in known)
    return known + extra


def vector_jump_points(rows: list[dict[str, str]], prefix: str) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    previous: tuple[float, float, float] | None = None
    for row in rows:
        current = (
            safe_float(row, f"total_{prefix}_x"),
            safe_float(row, f"total_{prefix}_y"),
            safe_float(row, f"total_{prefix}_z"),
        )
        if previous is None:
            jump = 0.0
        else:
            jump = math.sqrt(sum((current[i] - previous[i]) ** 2 for i in range(3)))
        points.append((safe_float(row, "time"), jump))
        previous = current
    return points


def scalar_points(rows: list[dict[str, str]], key: str) -> list[tuple[float, float]]:
    return [(safe_float(row, "time"), safe_float(row, key)) for row in rows]


def event_points(rows: list[dict[str, str]]) -> list[tuple[float, float]]:
    points = []
    for row in rows:
        count = (
            safe_int(row, "newborn_count")
            + safe_int(row, "merge_count")
            + safe_int(row, "split_count")
            + safe_int(row, "death_count")
        )
        points.append((safe_float(row, "time"), float(count)))
    return points


def format_tick(value: float) -> str:
    if abs(value) >= 1000.0:
        return f"{value:.0f}"
    if abs(value) >= 10.0:
        return f"{value:.1f}"
    if abs(value) >= 1.0:
        return f"{value:.2f}"
    return f"{value:.3f}"


def draw_line_chart(
    path: Path,
    title: str,
    y_label: str,
    series: list[tuple[str, list[tuple[float, float]]]],
    y_min: float | None = None,
    y_max: float | None = None,
    reference_lines: list[tuple[float, str]] | None = None,
) -> None:
    width = 1120
    height = 620
    left = 88
    right = 300
    top = 56
    bottom = 76
    plot_w = width - left - right
    plot_h = height - top - bottom
    reference_lines = reference_lines or []

    all_points = [pt for _, pts in series for pt in pts]
    if not all_points:
        return

    xmin = min(x for x, _ in all_points)
    xmax = max(x for x, _ in all_points)
    ymin = min(y for _, y in all_points) if y_min is None else y_min
    ymax = max(y for _, y in all_points) if y_max is None else y_max
    for ref, _ in reference_lines:
        ymin = min(ymin, ref)
        ymax = max(ymax, ref)
    if abs(xmax - xmin) < 1.0e-12:
        xmax = xmin + 1.0
    if abs(ymax - ymin) < 1.0e-12:
        pad = max(1.0, abs(ymax) * 0.1)
        ymin -= pad
        ymax += pad
    else:
        pad = 0.05 * (ymax - ymin)
        ymin = ymin if y_min is not None else ymin - pad
        ymax = ymax if y_max is not None else ymax + pad

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sy(y: float) -> float:
        return top + (ymax - y) / (ymax - ymin) * plot_h

    parts: list[str] = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append(f'<text x="{left}" y="30" font-family="Arial" font-size="22" font-weight="700" fill="#111827">{html.escape(title)}</text>')
    parts.append(f'<text x="{left + plot_w / 2}" y="{height - 24}" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">time (s)</text>')
    parts.append(f'<text transform="translate(24 {top + plot_h / 2}) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">{html.escape(y_label)}</text>')

    for i in range(6):
        tx = xmin + (xmax - xmin) * i / 5.0
        x = sx(tx)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#e5e7eb" stroke-width="1"/>')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_h + 24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#4b5563">{format_tick(tx)}</text>')

    for i in range(6):
        ty = ymin + (ymax - ymin) * i / 5.0
        y = sy(ty)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#e5e7eb" stroke-width="1"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#4b5563">{format_tick(ty)}</text>')

    parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#111827" stroke-width="1.2"/>')

    for ref, label in reference_lines:
        y = sy(ref)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#6b7280" stroke-dasharray="6 4" stroke-width="1.4"/>')
        parts.append(f'<text x="{left + plot_w + 8}" y="{y + 4:.2f}" font-family="Arial" font-size="12" fill="#6b7280">{html.escape(label)}</text>')

    for label, pts in series:
        if not pts:
            continue
        color = COLORS.get(label, "#111827")
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in pts)
        parts.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/>')

    legend_x = left + plot_w + 32
    legend_y = top + 10
    for i, (label, _) in enumerate(series):
        y = legend_y + i * 28
        color = COLORS.get(label, "#111827")
        parts.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 24}" y2="{y}" stroke="{color}" stroke-width="3" stroke-linecap="round"/>')
        parts.append(f'<text x="{legend_x + 34}" y="{y + 5}" font-family="Arial" font-size="13" fill="#111827">{html.escape(label)}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def draw_grouped_bar_chart(
    path: Path,
    title: str,
    y_label: str,
    rows: list[dict[str, str]],
    metric: str,
    y_min: float = 0.0,
) -> None:
    width = 1120
    height = 620
    left = 92
    right = 300
    top = 56
    bottom = 128
    plot_w = width - left - right
    plot_h = height - top - bottom

    scenarios = list(dict.fromkeys(row["scenario"] for row in rows))
    variants = [v for v in VARIANT_ORDER if any(row["variant"] == v for row in rows)]
    values = {(row["scenario"], row["variant"]): safe_float(row, metric) for row in rows}
    ymax = max([values.get((s, v), 0.0) for s in scenarios for v in variants] + [1.0])
    ymax *= 1.12

    def sy(y: float) -> float:
        return top + (ymax - y) / (ymax - y_min) * plot_h

    parts: list[str] = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append(f'<text x="{left}" y="30" font-family="Arial" font-size="22" font-weight="700" fill="#111827">{html.escape(title)}</text>')
    parts.append(f'<text transform="translate(24 {top + plot_h / 2}) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">{html.escape(y_label)}</text>')

    for i in range(6):
        ty = y_min + (ymax - y_min) * i / 5.0
        y = sy(ty)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#e5e7eb" stroke-width="1"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#4b5563">{format_tick(ty)}</text>')

    group_w = plot_w / max(1, len(scenarios))
    bar_gap = 3
    bar_w = max(5.0, (group_w * 0.72 - bar_gap * (len(variants) - 1)) / max(1, len(variants)))
    for si, scenario in enumerate(scenarios):
        group_x = left + si * group_w + group_w * 0.14
        for vi, variant in enumerate(variants):
            value = values.get((scenario, variant), 0.0)
            x = group_x + vi * (bar_w + bar_gap)
            y = sy(value)
            color = COLORS.get(variant, "#111827")
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{top + plot_h - y:.2f}" fill="{color}" opacity="0.88"/>')
        label_x = left + si * group_w + group_w / 2.0
        parts.append(f'<text transform="translate({label_x:.2f} {top + plot_h + 18}) rotate(20)" text-anchor="start" font-family="Arial" font-size="12" fill="#374151">{html.escape(scenario)}</text>')

    parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#111827" stroke-width="1.2"/>')

    legend_x = left + plot_w + 32
    legend_y = top + 10
    for i, variant in enumerate(variants):
        y = legend_y + i * 28
        color = COLORS.get(variant, "#111827")
        parts.append(f'<rect x="{legend_x}" y="{y - 10}" width="22" height="14" fill="{color}" opacity="0.88"/>')
        parts.append(f'<text x="{legend_x + 34}" y="{y + 2}" font-family="Arial" font-size="13" fill="#111827">{html.escape(variant)}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def generate_time_series_figures(frames: list[dict[str, str]], output_dir: Path) -> list[Path]:
    written: list[Path] = []
    groups = grouped_frames(frames)
    charts = [
        ("force_jump", "Force Jump", "||Delta F||", lambda rows: vector_jump_points(rows, "force"), []),
        ("torque_jump", "Torque Jump", "||Delta tau||", lambda rows: vector_jump_points(rows, "torque"), []),
        ("patch_count", "Patch Count", "patches", lambda rows: scalar_points(rows, "patch_count"), []),
        ("event_timeline", "Topology Event Timeline", "events per frame", event_points, []),
        ("coulomb_ratio", "Coulomb Ratio", "max ||Ft||/(mu Fn)", lambda rows: scalar_points(rows, "max_tangential_force_ratio"), [(1.0, "limit")]),
        ("energy_ratio", "Inherited Energy Ratio", "max inherited energy ratio", lambda rows: scalar_points(rows, "max_inherited_energy_ratio"), [(1.0, "non-amplification")]),
    ]

    for scenario, variants in groups.items():
        for suffix, label, ylabel, point_builder, refs in charts:
            series = []
            for variant in ordered_variants(variants):
                series.append((variant, point_builder(variants[variant])))
            filename = output_dir / f"{sanitize_name(scenario)}_{suffix}.svg"
            draw_line_chart(filename, f"{scenario}: {label}", ylabel, series, y_min=0.0, reference_lines=refs)
            written.append(filename)
    return written


def generate_summary_figures(summary: list[dict[str, str]], output_dir: Path) -> list[Path]:
    written: list[Path] = []
    specs = [
        ("summary_force_oscillation.svg", "Force Oscillation Index", "RMS force jump / max force", "force_oscillation_index"),
        ("summary_torque_oscillation.svg", "Torque Oscillation Index", "RMS torque jump / max torque", "torque_oscillation_index"),
        ("summary_patch_churn.svg", "Mean Patch Count Change", "mean |Delta patch count|", "mean_abs_patch_count_change"),
        ("summary_normal_jump.svg", "RMS Normal Jump", "radians", "rms_normal_jump_angle"),
        ("summary_energy_ratio.svg", "Max Inherited Energy Ratio", "energy ratio", "max_inherited_energy_ratio"),
    ]
    for filename, title, ylabel, metric in specs:
        path = output_dir / filename
        draw_grouped_bar_chart(path, title, ylabel, summary, metric)
        written.append(path)
    return written


def write_index(output_dir: Path, figures: list[Path]) -> Path:
    index = output_dir / "index.html"
    lines = [
        "<!doctype html>",
        "<html>",
        "<head><meta charset=\"utf-8\"><title>Field Contact Figures</title>",
        "<style>body{font-family:Arial,sans-serif;margin:28px;color:#111827} img{max-width:100%;border:1px solid #e5e7eb;margin:12px 0 28px} h2{margin-top:28px}</style>",
        "</head><body>",
        "<h1>Field Contact Benchmark Figures</h1>",
    ]
    for figure in figures:
        rel = figure.name
        title = figure.stem.replace("_", " ")
        lines.append(f"<h2>{html.escape(title)}</h2>")
        lines.append(f"<img src=\"{html.escape(rel)}\" alt=\"{html.escape(title)}\">")
    lines.append("</body></html>")
    index.write_text("\n".join(lines), encoding="utf-8")
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("out/milestone_21"), help="directory containing benchmark CSV files")
    parser.add_argument("--output", type=Path, default=Path("out/milestone_22/figures"), help="directory for generated SVG figures")
    parser.add_argument(
        "--prefix",
        default="field_contact_multislip_interlock",
        help="CSV file prefix, for example field_contact_multislip_interlock or mesh_sdf_nonconvex",
    )
    args = parser.parse_args()

    frames_path = args.input / f"{args.prefix}_frames.csv"
    summary_path = args.input / f"{args.prefix}_summary.csv"
    if not frames_path.exists():
        raise FileNotFoundError(frames_path)
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)

    args.output.mkdir(parents=True, exist_ok=True)
    frames = read_csv(frames_path)
    summary = read_csv(summary_path)

    figures = []
    figures.extend(generate_time_series_figures(frames, args.output))
    figures.extend(generate_summary_figures(summary, args.output))
    index = write_index(args.output, figures)

    print(f"Generated {len(figures)} SVG figures")
    print(index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
