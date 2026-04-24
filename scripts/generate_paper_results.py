#!/usr/bin/env python3
"""Generate paper-ready result tables and figure snippets for Milestone 25."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


SCENE_LABELS = {
    "folded_seam": "折线拼接面",
    "concave_polyline_corner": "非凸折线角",
    "narrow_groove_entrance": "窄槽入口",
    "guide_rail_sliding": "导轨/槽道滑移",
    "nested_interlock": "嵌套互锁接触",
    "multi_patch_rolling_sliding": "多 patch 滚滑",
}

VARIANT_LABELS = {
    "current_field_primitive": "当前 field primitive",
    "minimal_rotation_gate_split_merge": "最小旋转+gate+split/merge",
    "no_history": "无历史",
    "direct_inherit": "直接继承历史",
    "projection_transport": "简单投影 transport",
    "field_primitive": "field primitive",
    "raw_sample_contacts": "逐采样点接触",
    "normal_bin_fragments": "法向分箱碎片",
    "convex_decomposition_proxy": "凸分解代理",
    "chrono_traditional_proxy": "Chrono 点流形代理",
}

VARIANT_ORDER = [
    "field_primitive",
    "convex_decomposition_proxy",
    "chrono_traditional_proxy",
    "normal_bin_fragments",
    "raw_sample_contacts",
]

COLORS = {
    "field_primitive": (0.06, 0.46, 0.43),
    "convex_decomposition_proxy": (0.49, 0.23, 0.93),
    "chrono_traditional_proxy": (0.86, 0.15, 0.15),
    "normal_bin_fragments": (0.92, 0.35, 0.04),
    "raw_sample_contacts": (0.15, 0.39, 0.92),
}

REPRESENTATIVE_FIGURES = [
    (
        "guide_rail_sliding",
        "patch_count",
        "guide_rail_sliding_patch_count.pdf",
        "fig:guide_patch_count",
        "导轨/槽道滑移场景中的 patch count 时间序列。baseline 包括逐采样点接触、法向分箱碎片、凸分解代理和 Chrono 风格点流形代理；指标为每一时间步的接触原语数量。结果显示，field primitive 将多点接触组织成少量连续 patch，而逐点和碎片化 baseline 产生明显更大的接触集合波动。",
    ),
    (
        "guide_rail_sliding",
        "force_jump",
        "guide_rail_sliding_force_jump.pdf",
        "fig:guide_force_jump",
        "导轨/槽道滑移场景中的合力跳变量时间序列。baseline 与图 \\ref{fig:guide_patch_count} 相同；指标为相邻时间步总接触力差的范数。该图用于观察接触集合切换向力响应的传播，结果表明 field primitive 在保持紧凑接触集合的同时避免了 Chrono 风格点流形代理中的大幅力跳变。",
    ),
    (
        "nested_interlock",
        "torque_jump",
        "nested_interlock_torque_jump.pdf",
        "fig:nested_torque_jump",
        "嵌套互锁场景中的力矩跳变量时间序列。baseline 包括逐点接触、法向碎片、凸分解代理和 Chrono 风格点流形代理；指标为相邻时间步总接触力矩差的范数。结果显示，互锁几何中的局部特征切换会放大点流形代理的力矩波动，而 field primitive 保持较稳定的 patch 级力矩响应。",
    ),
    (
        "multi_patch_rolling_sliding",
        "event_timeline",
        "multi_patch_rolling_sliding_event_timeline.pdf",
        "fig:multipatch_event_timeline",
        "多 patch 滚滑场景中的事件时间线。baseline 与其他非凸实验一致；指标为每步 newborn、death、split 和 merge 事件总数。结果显示，逐点接触和法向碎片化 baseline 持续产生大量拓扑事件，而 field primitive 的事件主要集中在真实接触区域出现、消失或合并的时间段。",
    ),
]


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


def fmt(value: float, digits: int = 3) -> str:
    if abs(value) >= 100.0:
        return f"{value:.1f}"
    if abs(value) >= 10.0:
        return f"{value:.2f}"
    return f"{value:.{digits}f}"


def tex_escape(value: str) -> str:
    return (
        value.replace("\\", "\\textbackslash{}")
        .replace("_", "\\_")
        .replace("%", "\\%")
        .replace("&", "\\&")
        .replace("#", "\\#")
    )


def scene_label(name: str) -> str:
    return SCENE_LABELS.get(name, name.replace("_", " "))


def variant_label(name: str) -> str:
    return VARIANT_LABELS.get(name, name.replace("_", " "))


def begin_resize_table(caption: str, label: str, columns: str, header: str) -> list[str]:
    return [
        "\\begin{table}[H]",
        "\\centering",
        "\\small",
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        "\\resizebox{\\textwidth}{!}{%",
        f"\\begin{{tabular}}{{{columns}}}",
        "\\toprule",
        header,
        "\\midrule",
    ]


def end_resize_table() -> list[str]:
    return [
        "\\bottomrule",
        "\\end{tabular}%",
        "}",
        "\\end{table}",
        "",
    ]


def write_m20_ablation_table(rows: list[dict[str, str]]) -> str:
    chosen = {
        "minimal_rotation_gate_split_merge",
        "no_history",
        "direct_inherit",
        "projection_transport",
    }
    lines = begin_resize_table(
        "特征切换敏感场景的切向历史消融结果。表中 $E_{\\mathrm{inh}}$ 为继承弹性能相对来源能量上界的最大比值；Coulomb 比值不超过 1 表示切向力满足摩擦圆约束。",
        "tab:feature_switching_ablation",
        "llrrrrrr",
        "场景 & 方法 & $\\overline{N}_{p}$ & $N_{p}^{\\max}$ & merge & split & $\\rho_T^{\\max}$ & $E_{\\mathrm{inh}}^{\\max}$ \\\\",
    )
    for row in rows:
        variant = row["variant"]
        if variant not in chosen:
            continue
        scene = scene_label(row["scenario"])
        merge_count = safe_int(row, "total_merge_patches")
        split_count = safe_int(row, "total_split_patches")
        lines.append(
            f"{tex_escape(scene)} & {tex_escape(variant_label(variant))} & "
            f"{fmt(safe_float(row, 'mean_patch_count'))} & {safe_int(row, 'max_patch_count')} & "
            f"{merge_count} & {split_count} & {fmt(safe_float(row, 'max_tangential_force_ratio'))} & "
            f"{fmt(safe_float(row, 'max_inherited_energy_ratio'))} \\\\"
        )
    lines.extend(end_resize_table())
    return "\n".join(lines)


def write_m21_field_summary(rows: list[dict[str, str]]) -> str:
    lines = begin_resize_table(
        "非凸多点接触场景中 field primitive 的总体指标。$\\Delta N_p$ 表示相邻时间步 primitive 数量变化的平均绝对值，$J_F$ 和 $J_\\tau$ 分别为合力和力矩振荡指数。",
        "tab:field_summary",
        "lrrrrrrrr",
        "场景 & active & $\\overline{N}_{p}$ & $N_{p}^{\\max}$ & $\\Delta N_p$ & RMS 法向跳变 & $J_F$ & $J_\\tau$ & Coulomb \\\\",
    )
    for row in rows:
        if row["variant"] != "field_primitive":
            continue
        lines.append(
            f"{tex_escape(scene_label(row['scenario']))} & {safe_int(row, 'active_frames')} & "
            f"{fmt(safe_float(row, 'mean_patch_count'))} & {safe_int(row, 'max_patch_count')} & "
            f"{fmt(safe_float(row, 'mean_abs_patch_count_change'))} & "
            f"{fmt(safe_float(row, 'rms_normal_jump_angle'))} & "
            f"{fmt(safe_float(row, 'force_oscillation_index'))} & "
            f"{fmt(safe_float(row, 'torque_oscillation_index'))} & "
            f"{fmt(safe_float(row, 'max_tangential_force_ratio'))} \\\\"
        )
    lines.extend(end_resize_table())
    return "\n".join(lines)


def write_m21_ratio_table(rows: list[dict[str, str]]) -> str:
    preferred = {"raw_sample_contacts", "normal_bin_fragments", "convex_decomposition_proxy", "chrono_traditional_proxy"}
    lines = begin_resize_table(
        "非凸多点接触场景的 baseline ratio。所有比值均为 baseline 除以 field primitive；大于 1 表示 baseline 在该指标上更不稳定。",
        "tab:baseline_ratios",
        "llrrrrrr",
        "场景 & baseline & $\\Delta N_p$ ratio & $N_p^{\\max}$ ratio & RMS 法向 ratio & $J_F$ ratio & $J_\\tau$ ratio & 事件 ratio \\\\",
    )
    for row in rows:
        baseline = row["baseline_variant"]
        if baseline not in preferred:
            continue
        lines.append(
            f"{tex_escape(scene_label(row['scenario']))} & {tex_escape(variant_label(baseline))} & "
            f"{fmt(safe_float(row, 'mean_patch_count_change_ratio'))} & "
            f"{fmt(safe_float(row, 'max_patch_count_change_ratio'))} & "
            f"{fmt(safe_float(row, 'rms_normal_jump_ratio'))} & "
            f"{fmt(safe_float(row, 'force_oscillation_ratio'))} & "
            f"{fmt(safe_float(row, 'torque_oscillation_ratio'))} & "
            f"{fmt(safe_float(row, 'topology_event_ratio'))} \\\\"
        )
    lines.extend(end_resize_table())
    return "\n".join(lines)


def group_frame_rows(rows: list[dict[str, str]]) -> dict[str, dict[str, list[dict[str, str]]]]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        grouped[row["scenario"]][row["variant"]].append(row)
    for variants in grouped.values():
        for values in variants.values():
            values.sort(key=lambda r: safe_int(r, "frame"))
    return grouped


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


def scalar_points(rows: list[dict[str, str]], key: str) -> list[tuple[float, float]]:
    return [(safe_float(row, "time"), safe_float(row, key)) for row in rows]


def points_for_metric(rows: list[dict[str, str]], metric: str) -> list[tuple[float, float]]:
    if metric == "patch_count":
        return scalar_points(rows, "patch_count")
    if metric == "force_jump":
        return vector_jump_points(rows, "force")
    if metric == "torque_jump":
        return vector_jump_points(rows, "torque")
    if metric == "event_timeline":
        return event_points(rows)
    raise ValueError(metric)


def metric_title(metric: str) -> tuple[str, str]:
    if metric == "patch_count":
        return "Patch Count", "patches"
    if metric == "force_jump":
        return "Force Jump", "||Delta F||"
    if metric == "torque_jump":
        return "Torque Jump", "||Delta tau||"
    if metric == "event_timeline":
        return "Topology Events", "events / step"
    return metric, metric


def draw_pdf_chart(path: Path, title: str, y_label: str, series: list[tuple[str, list[tuple[float, float]]]]) -> None:
    try:
        from reportlab.pdfgen import canvas
    except Exception as exc:  # pragma: no cover - environment fallback
        raise RuntimeError("reportlab is required for PDF figure generation") from exc

    path.parent.mkdir(parents=True, exist_ok=True)
    width = 640.0
    height = 360.0
    left = 62.0
    right = 176.0
    top = 38.0
    bottom = 52.0
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_points = [point for _, points in series for point in points]
    xmin = min(x for x, _ in all_points)
    xmax = max(x for x, _ in all_points)
    ymin = min(y for _, y in all_points)
    ymax = max(y for _, y in all_points)
    if abs(xmax - xmin) < 1.0e-12:
        xmax = xmin + 1.0
    if abs(ymax - ymin) < 1.0e-12:
        pad = max(1.0, abs(ymax) * 0.1)
        ymin -= pad
        ymax += pad
    else:
        pad = 0.06 * (ymax - ymin)
        ymin -= pad
        ymax += pad

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sy(y: float) -> float:
        return bottom + (y - ymin) / (ymax - ymin) * plot_h

    c = canvas.Canvas(str(path), pagesize=(width, height))
    c.setFont("Helvetica-Bold", 14)
    c.drawString(left, height - 24, title)
    c.setFont("Helvetica", 9)
    c.drawCentredString(left + plot_w / 2.0, 18, "time (s)")
    c.saveState()
    c.translate(18, bottom + plot_h / 2.0)
    c.rotate(90)
    c.drawCentredString(0, 0, y_label)
    c.restoreState()

    c.setStrokeColorRGB(0.88, 0.90, 0.93)
    c.setLineWidth(0.5)
    for i in range(6):
        x = left + plot_w * i / 5.0
        c.line(x, bottom, x, bottom + plot_h)
        tick = xmin + (xmax - xmin) * i / 5.0
        c.setFillColorRGB(0.28, 0.32, 0.38)
        c.drawCentredString(x, bottom - 16, fmt(tick, 2))
    for i in range(6):
        y = bottom + plot_h * i / 5.0
        c.setStrokeColorRGB(0.88, 0.90, 0.93)
        c.line(left, y, left + plot_w, y)
        tick = ymin + (ymax - ymin) * i / 5.0
        c.setFillColorRGB(0.28, 0.32, 0.38)
        c.drawRightString(left - 8, y - 3, fmt(tick, 2))

    c.setStrokeColorRGB(0.07, 0.09, 0.15)
    c.setLineWidth(0.8)
    c.rect(left, bottom, plot_w, plot_h, fill=0)

    for variant, points in series:
        r, g, b = COLORS.get(variant, (0.07, 0.09, 0.15))
        c.setStrokeColorRGB(r, g, b)
        c.setLineWidth(1.4)
        if len(points) >= 2:
            path_obj = c.beginPath()
            path_obj.moveTo(sx(points[0][0]), sy(points[0][1]))
            for x, y in points[1:]:
                path_obj.lineTo(sx(x), sy(y))
            c.drawPath(path_obj, stroke=1, fill=0)

    legend_x = left + plot_w + 22.0
    legend_y = bottom + plot_h - 12.0
    c.setFont("Helvetica", 8.5)
    for i, (variant, _) in enumerate(series):
        y = legend_y - 17.0 * i
        r, g, b = COLORS.get(variant, (0.07, 0.09, 0.15))
        c.setStrokeColorRGB(r, g, b)
        c.setLineWidth(2.2)
        c.line(legend_x, y, legend_x + 18.0, y)
        c.setFillColorRGB(0.07, 0.09, 0.15)
        c.drawString(legend_x + 24.0, y - 3.0, variant)

    c.showPage()
    c.save()


def write_figures(frame_rows: list[dict[str, str]], fig_output: Path) -> str:
    grouped = group_frame_rows(frame_rows)
    lines = [
        "% Auto-generated by scripts/generate_paper_results.py",
        "% Representative figures generated from out/milestone_21 frame CSV.",
        "",
    ]
    for scenario, metric, filename, label, caption in REPRESENTATIVE_FIGURES:
        variants = grouped[scenario]
        series = []
        for variant in VARIANT_ORDER:
            if variant in variants:
                series.append((variant, points_for_metric(variants[variant], metric)))
        title, y_label = metric_title(metric)
        draw_pdf_chart(fig_output / filename, f"{scene_label(scenario)} - {title}", y_label, series)
        tex_path = Path("figures") / "milestone_25" / filename
        repo_tex_path = Path("paper") / tex_path
        lines.extend(
            [
                "\\begin{figure}[H]",
                "\\centering",
                f"\\IfFileExists{{{repo_tex_path.as_posix()}}}"
                f"{{\\includegraphics[width=0.96\\textwidth]{{{repo_tex_path.as_posix()}}}}}"
                f"{{\\includegraphics[width=0.96\\textwidth]{{{tex_path.as_posix()}}}}}",
                f"\\caption{{{caption}}}",
                f"\\label{{{label}}}",
                "\\end{figure}",
                "",
            ]
        )
    return "\n".join(lines)


def write_caption_index(m22_figures: Path, output: Path) -> None:
    lines = ["# Milestone 25 Representative Figure Captions", ""]
    for scenario, metric, filename, label, caption in REPRESENTATIVE_FIGURES:
        source_svg = m22_figures / filename.replace(".pdf", ".svg")
        lines.append(f"## `{filename}`")
        lines.append("")
        lines.append(f"- Label: `{label}`")
        lines.append(f"- Data: `out/milestone_21/field_contact_multislip_interlock_frames.csv`")
        lines.append(f"- Matching Milestone 22 SVG: `{source_svg.as_posix()}`")
        lines.append(f"- Scene: `{scenario}`")
        lines.append(f"- Metric: `{metric}`")
        lines.append(f"- Caption: {caption}")
        lines.append("")
    output.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m20", type=Path, default=Path("out/milestone_20"))
    parser.add_argument("--m21", type=Path, default=Path("out/milestone_21"))
    parser.add_argument("--m22-figures", type=Path, default=Path("out/milestone_22/figures"))
    parser.add_argument("--output", type=Path, default=Path("paper/generated/milestone_25"))
    parser.add_argument("--fig-output", type=Path, default=Path("paper/figures/milestone_25"))
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    args.fig_output.mkdir(parents=True, exist_ok=True)

    m20_summary = read_csv(args.m20 / "field_contact_feature_switching_summary.csv")
    m21_summary = read_csv(args.m21 / "field_contact_multislip_interlock_summary.csv")
    m21_comparison = read_csv(args.m21 / "field_contact_multislip_interlock_comparison.csv")
    m21_frames = read_csv(args.m21 / "field_contact_multislip_interlock_frames.csv")

    tables = [
        "% Auto-generated by scripts/generate_paper_results.py",
        "",
        write_m20_ablation_table(m20_summary),
        write_m21_field_summary(m21_summary),
        write_m21_ratio_table(m21_comparison),
    ]
    (args.output / "results_tables.tex").write_text("\n\n".join(tables), encoding="utf-8")

    figures_tex = write_figures(m21_frames, args.fig_output)
    (args.output / "representative_figures.tex").write_text(figures_tex, encoding="utf-8")
    write_caption_index(args.m22_figures, args.output / "figure_captions.md")

    print(args.output / "results_tables.tex")
    print(args.output / "representative_figures.tex")
    print(args.output / "figure_captions.md")
    print(args.fig_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
