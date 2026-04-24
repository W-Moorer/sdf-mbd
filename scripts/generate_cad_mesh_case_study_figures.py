#!/usr/bin/env python3
"""Generate SVG paper figures for Milestone 29 CAD/mesh case study."""

from __future__ import annotations

import argparse
import csv
import html
import math
from collections import defaultdict
from pathlib import Path


COLORS = {
    "field_primitive_openvdb": "#0f766e",
    "native_mesh_baseline": "#dc2626",
    "trajectory": "#374151",
}

LABELS = {
    "field_primitive_openvdb": "field primitive",
    "native_mesh_baseline": "Chrono native",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def grouped(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        out[row.get("variant", "")].append(row)
    for variant_rows in out.values():
        variant_rows.sort(key=lambda row: f(row, "frame"))
    return out


def line_points(rows: list[dict[str, str]], x_key: str, y_key: str) -> list[tuple[float, float]]:
    return [(f(row, x_key), f(row, y_key)) for row in rows]


def tick(value: float) -> str:
    if abs(value) >= 100.0:
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
) -> None:
    width = 1120
    height = 620
    left = 96
    right = 260
    top = 58
    bottom = 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    points = [pt for _, pts in series for pt in pts]
    if not points:
        return

    xmin = min(x for x, _ in points)
    xmax = max(x for x, _ in points)
    ymin = min(y for _, y in points) if y_min is None else y_min
    ymax = max(y for _, y in points)
    if abs(xmax - xmin) < 1.0e-12:
        xmax = xmin + 1.0
    if abs(ymax - ymin) < 1.0e-12:
        ymax = ymin + 1.0
    pad = 0.06 * (ymax - ymin)
    ymax += pad
    if y_min is None:
        ymin -= pad

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sy(y: float) -> float:
        return top + (ymax - y) / (ymax - ymin) * plot_h

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left}" y="32" font-family="Arial" font-size="22" font-weight="700" fill="#111827">{html.escape(title)}</text>',
        f'<text x="{left + plot_w / 2}" y="{height - 24}" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">time (s)</text>',
        f'<text transform="translate(26 {top + plot_h / 2}) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">{html.escape(y_label)}</text>',
    ]

    for i in range(6):
        tx = xmin + (xmax - xmin) * i / 5.0
        x = sx(tx)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_h + 24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#4b5563">{tick(tx)}</text>')
    for i in range(6):
        ty = ymin + (ymax - ymin) * i / 5.0
        y = sy(ty)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#4b5563">{tick(ty)}</text>')
    parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#111827" stroke-width="1.2"/>')

    for label, pts in series:
        if not pts:
            continue
        color = COLORS.get(label, "#111827")
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in pts)
        parts.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"/>')

    legend_x = left + plot_w + 34
    legend_y = top + 18
    for i, (label, _) in enumerate(series):
        y = legend_y + i * 30
        color = COLORS.get(label, "#111827")
        parts.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 26}" y2="{y}" stroke="{color}" stroke-width="3.2" stroke-linecap="round"/>')
        parts.append(f'<text x="{legend_x + 36}" y="{y + 5}" font-family="Arial" font-size="13" fill="#111827">{html.escape(LABELS.get(label, label))}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def draw_trajectory(path: Path, rows: list[dict[str, str]]) -> None:
    field_rows = [row for row in rows if row.get("variant") == "field_primitive_openvdb"]
    if not field_rows:
        return
    points = [(f(row, "center_x"), f(row, "center_z"), f(row, "patch_count")) for row in field_rows]
    width = 760
    height = 620
    left = 84
    right = 48
    top = 62
    bottom = 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    xmin = min(x for x, _, _ in points)
    xmax = max(x for x, _, _ in points)
    zmin = min(z for _, z, _ in points)
    zmax = max(z for _, z, _ in points)
    if abs(xmax - xmin) < 1.0e-12:
        xmax = xmin + 1.0
    if abs(zmax - zmin) < 1.0e-12:
        zmax = zmin + 1.0
    xpad = 0.06 * (xmax - xmin)
    zpad = 0.10 * (zmax - zmin)
    xmin -= xpad
    xmax += xpad
    zmin -= zpad
    zmax += zpad

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sz(z: float) -> float:
        return top + (zmax - z) / (zmax - zmin) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left}" y="32" font-family="Arial" font-size="22" font-weight="700" fill="#111827">CAD mesh path and active patches</text>',
        f'<text x="{left + plot_w / 2}" y="{height - 24}" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">x position (m)</text>',
        f'<text transform="translate(26 {top + plot_h / 2}) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="14" fill="#374151">z position (m)</text>',
    ]
    for i in range(6):
        tx = xmin + (xmax - xmin) * i / 5.0
        x = sx(tx)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_h + 24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#4b5563">{tick(tx)}</text>')
        tz = zmin + (zmax - zmin) * i / 5.0
        y = sz(tz)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#4b5563">{tick(tz)}</text>')
    parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#111827" stroke-width="1.2"/>')
    path_poly = " ".join(f"{sx(x):.2f},{sz(z):.2f}" for x, z, _ in points)
    parts.append(f'<polyline points="{path_poly}" fill="none" stroke="{COLORS["trajectory"]}" stroke-width="1.8"/>')
    for x, z, patches in points[::6]:
        if patches <= 0:
            continue
        r = min(7.0, 2.2 + 1.4 * math.sqrt(patches))
        parts.append(f'<circle cx="{sx(x):.2f}" cy="{sz(z):.2f}" r="{r:.2f}" fill="#0f766e" fill-opacity="0.42" stroke="#0f766e" stroke-width="0.7"/>')
    parts.append('<text x="500" y="86" font-family="Arial" font-size="13" fill="#0f766e">dot size: field patch count</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_manifest(path: Path, captions: list[tuple[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["figure", "caption"])
        writer.writerows(captions)


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--frames", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    root = args.project_root.resolve()
    frames_path = args.frames or root / "out" / "milestone_29" / "cad_mesh_case_study_frames.csv"
    out_dir = args.output_dir or root / "out" / "milestone_29" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_csv(frames_path)
    by_variant = grouped(rows)
    series_counts = [
        ("field_primitive_openvdb", line_points(by_variant.get("field_primitive_openvdb", []), "time", "patch_count")),
        ("native_mesh_baseline", line_points(by_variant.get("native_mesh_baseline", []), "time", "contact_count")),
    ]
    series_force = [
        ("field_primitive_openvdb", line_points(by_variant.get("field_primitive_openvdb", []), "time", "total_force_norm")),
        ("native_mesh_baseline", line_points(by_variant.get("native_mesh_baseline", []), "time", "total_force_norm")),
    ]
    series_torque = [
        ("field_primitive_openvdb", line_points(by_variant.get("field_primitive_openvdb", []), "time", "total_torque_norm")),
        ("native_mesh_baseline", line_points(by_variant.get("native_mesh_baseline", []), "time", "total_torque_norm")),
    ]

    counts_path = out_dir / "cad_mesh_case_study_counts.svg"
    force_path = out_dir / "cad_mesh_case_study_force_norm.svg"
    torque_path = out_dir / "cad_mesh_case_study_torque_norm.svg"
    trajectory_path = out_dir / "cad_mesh_case_study_trajectory.svg"
    draw_line_chart(counts_path, "CAD mesh contact-set timeline", "patch/contact count", series_counts, y_min=0.0)
    draw_line_chart(force_path, "CAD mesh contact force", "force norm (N)", series_force, y_min=0.0)
    draw_line_chart(torque_path, "CAD mesh contact torque", "torque norm (N m)", series_torque, y_min=0.0)
    draw_trajectory(trajectory_path, rows)

    captions = [
        (
            counts_path.name,
            "Track-shoe CAD mesh case study: field primitives keep a low-churn patch timeline while Chrono native contact appears as sparse contact events on the same commanded path.",
        ),
        (
            force_path.name,
            "Contact force norm for the same target mesh and moving probe; the field primitive force changes more smoothly than the native contact baseline.",
        ),
        (
            torque_path.name,
            "Contact torque norm for the real mesh case; lower field oscillation supports the paper claim that primitive history and transport reduce moment chatter.",
        ),
        (
            trajectory_path.name,
            "Top-view commanded path over the real mesh target; marker size indicates active field patch count during the slide.",
        ),
    ]
    manifest_path = out_dir / "cad_mesh_case_study_figure_manifest.csv"
    write_manifest(manifest_path, captions)
    print(f"wrote {counts_path}")
    print(f"wrote {force_path}")
    print(f"wrote {torque_path}")
    print(f"wrote {trajectory_path}")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
