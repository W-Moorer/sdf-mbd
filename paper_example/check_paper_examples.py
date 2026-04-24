import argparse
import csv
import json
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".", help="Repository root")
    return parser.parse_args()


def load_manifest(paper_dir):
    with (paper_dir / "manifest.json").open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_summary(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing summary file: {path}")
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def load_frames(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing frame file: {path}")
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def require_files(project_root, manifest):
    missing = []
    paper_dir = project_root / "paper_example"
    for case in manifest["cases"]:
        paths = [case["model"]] + case["models"]
        for rel in paths:
            path = paper_dir / rel
            if not path.exists():
                missing.append(str(path))
    if missing:
        raise FileNotFoundError("Missing locked paper-example inputs:\n" + "\n".join(missing))


def validate_summary(manifest, rows):
    row_by_key = {(row["case"], row["quantity"]): row for row in rows}
    failures = []
    for case in manifest["cases"]:
        for check in case["checks"]:
            key = (case["name"], check["quantity"])
            row = row_by_key.get(key)
            if row is None:
                failures.append(f"Missing summary row for {key[0]}:{key[1]}")
                continue
            max_error = float(row["max_abs_error"])
            tolerance = float(check["tolerance"])
            row_tolerance = float(row["tolerance"])
            passed = row["passed"].strip().lower() == "true"
            if abs(row_tolerance - tolerance) > max(1.0e-12, tolerance * 1.0e-9):
                failures.append(f"{key[0]}:{key[1]} tolerance changed: {row_tolerance} != {tolerance}")
            if max_error > tolerance or not passed:
                failures.append(f"{key[0]}:{key[1]} failed: max_error={max_error} tolerance={tolerance}")
    if failures:
        raise AssertionError("\n".join(failures))


def validate_frames(manifest, rows):
    if not rows:
        raise AssertionError("Frame CSV is empty")
    expected_cases = {case["name"] for case in manifest["cases"]}
    observed_cases = {row["case"] for row in rows}
    missing_cases = expected_cases - observed_cases
    if missing_cases:
        raise AssertionError("Missing frame data for cases: " + ", ".join(sorted(missing_cases)))
    bad_backends = sorted({row["backend"] for row in rows if row["backend"] != "sparse_sdf_contact_force_dynamics"})
    if bad_backends:
        raise AssertionError("Unexpected backend labels: " + ", ".join(bad_backends))

    active_counts = {case: 0 for case in expected_cases}
    for row in rows:
        if int(row["patch_count"]) > 0:
            active_counts[row["case"]] += 1
    inactive = [case for case, count in active_counts.items() if count == 0]
    if inactive:
        raise AssertionError("No contact patches were reported for cases: " + ", ".join(sorted(inactive)))


def validate_outputs(project_root, manifest):
    output_dir = project_root / manifest["output_dir"]
    missing = []
    for rel in manifest.get("outputs", []):
        path = output_dir / rel
        if not path.exists():
            missing.append(str(path))
    if missing:
        raise FileNotFoundError("Missing paper-example outputs:\n" + "\n".join(missing))


def main():
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    paper_dir = project_root / "paper_example"
    manifest = load_manifest(paper_dir)
    require_files(project_root, manifest)

    output_dir = project_root / manifest["output_dir"]
    summary = load_summary(output_dir / "comparison_summary.csv")
    frames = load_frames(output_dir / "sparse_sdf_frames.csv")
    validate_summary(manifest, summary)
    validate_frames(manifest, frames)
    validate_outputs(project_root, manifest)

    print("paper_example regression passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"paper_example regression failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
