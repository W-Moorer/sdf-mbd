# Milestone 25: Paper Tables, Figures, and Result Section

This milestone turns the Milestone 20/21/22 outputs into paper-ready material.

## Generator

Script:

```text
scripts/generate_paper_results.py
```

Run:

```powershell
python scripts\generate_paper_results.py
```

Inputs:

```text
out/milestone_20/field_contact_feature_switching_summary.csv
out/milestone_21/field_contact_multislip_interlock_summary.csv
out/milestone_21/field_contact_multislip_interlock_comparison.csv
out/milestone_21/field_contact_multislip_interlock_frames.csv
out/milestone_22/figures/*.svg
```

Outputs:

```text
paper/generated/milestone_25/results_tables.tex
paper/generated/milestone_25/representative_figures.tex
paper/generated/milestone_25/figure_captions.md
paper/figures/milestone_25/*.pdf
```

The PDF figures are generated directly from the frame CSV using ReportLab, so
the TeX build does not depend on SVG support, Inkscape, or shell escape.

## Paper Integration

Updated TeX:

```text
paper/field_based_contact_primitives_cn_draft_v25_tangential_revised.tex
```

Changes:

- The abstract now states experimental results instead of expected results.
- The implementation section now describes the current OpenVDB/Chrono runtime.
- The old `实验设计` and `预期结果与讨论` sections were replaced with
  `实验设置与结果`.
- The result section inputs generated tables and representative figures.
- The conclusion now states the verified result scope.

## Generated Tables

The generated table snippet contains:

- `tab:feature_switching_ablation`: Milestone 20 tangential-history ablation.
- `tab:field_summary`: Milestone 21 field-primitive summary table.
- `tab:baseline_ratios`: Milestone 21/22 baseline ratio table.

Key paper points:

- In `narrow_groove_entrance`, direct history inheritance reaches
  `max_inherited_energy_ratio=2.000`, while minimal-rotation + gate +
  split/merge remains at `1.000`.
- Field primitive keeps `max_tangential_force_ratio=1.000` in the reported
  experiments.
- Baselines generally show much larger patch-count churn and topology-event
  ratios than `field_primitive`.

## Generated Figures

Representative figures:

- `guide_rail_sliding_patch_count.pdf`
- `guide_rail_sliding_force_jump.pdf`
- `nested_interlock_torque_jump.pdf`
- `multi_patch_rolling_sliding_event_timeline.pdf`

Each figure has a fixed caption covering:

- scenario;
- baseline set;
- metric meaning;
- conclusion.

## Verification

The paper was compiled with:

```powershell
xelatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/field_based_contact_primitives_cn_draft_v25_tangential_revised.tex
```

The final compile succeeds.  Remaining messages are only underfull-box warnings
from long bibliography URLs.
