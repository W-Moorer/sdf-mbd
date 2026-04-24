# Milestone 29: Real CAD/Mesh Case Study

## Scope

This milestone moves the field-contact primitive pipeline from generated analytic
height fields and code-generated meshes to repository mesh assets:

- Target mesh: `data/models/bulldozer/shoe_collision.obj`
- Moving mesh: `data/models/semicapsule.obj`
- SDF backend: OpenVDB level set built directly from the target triangle mesh
- Native baseline: Chrono SMC static triangle mesh target plus a dynamic convex
  hull generated from the same moving mesh vertices

The target mesh is shared by the field SDF and Chrono native collision. The
moving mesh is convex and scaled to `0.22`, so the native convex-hull
representation remains a close geometric counterpart while using Chrono's
supported convex-concave collision path.

## Outputs

- `out/milestone_29/cad_mesh_case_study_frames.csv`
- `out/milestone_29/cad_mesh_case_study_comparison.csv`
- `out/milestone_29/cad_mesh_case_study_geometry.csv`
- `out/milestone_29/cad_mesh_case_study_regression_summary.json`
- `out/milestone_29/figures/cad_mesh_case_study_counts.svg`
- `out/milestone_29/figures/cad_mesh_case_study_force_norm.svg`
- `out/milestone_29/figures/cad_mesh_case_study_torque_norm.svg`
- `out/milestone_29/figures/cad_mesh_case_study_trajectory.svg`
- `out/milestone_29/figures/cad_mesh_case_study_figure_manifest.csv`

## Current Results

Scenario: `cad_track_shoe_semicapsule_slide`

| Metric | Field primitive | Chrono native | Native / field |
|---|---:|---:|---:|
| active frames | 88 | 14 | 0.159 |
| mean patch/contact count | 0.125535 | 0.019971 | 0.159 |
| mean count churn | 0.004286 | 0.025714 | 6.000 |
| RMS normal jump angle | 0.042067 | 0.091808 | 2.182 |
| force oscillation index | 0.011430 | 0.068697 | 6.010 |
| torque oscillation index | 0.011324 | 0.055464 | 4.898 |
| max penetration | 0.002009 | 0.001949 | - |
| max Coulomb ratio | 1.000000 | 0.99999997 | - |

Geometry metadata:

| Quantity | Value |
|---|---:|
| target vertices / faces | 40 / 72 |
| moving vertices / faces | 3767 / 7530 |
| moving scale | 0.22 |
| OpenVDB active voxels | 1,288,920 |
| voxel size | 0.0015 m |

## Interpretation

This case is no longer an analytic height-field validation. The target contact
field is generated from a real repository triangle mesh, then consumed by the
same primitive tracker, history transport, gate, and split/merge machinery used
in the earlier milestones.

The result supports the TeX claim in a narrow but concrete setting: when a
moving probe slides across a non-smooth mesh-derived target, field primitives
produce a contact description with lower topology churn, lower normal jump, and
lower force/torque oscillation than the native point-contact baseline, while
respecting the Coulomb bound and keeping penetration controlled.

The important regression checks are:

- real mesh assets exist under `data/models`
- OpenVDB SDF contains a nontrivial active narrow band
- both field and native methods produce contact
- field and native Coulomb ratios stay at or below `1 + tol`
- native contact-count churn is at least `3x` field patch churn
- native normal jump is at least `1.5x` field normal jump
- native force oscillation is at least `3x` field force oscillation
- native torque oscillation is at least `2.5x` field torque oscillation

## Limitations

Chrono native dynamic triangle mesh against static triangle mesh did not produce
stable contacts for this case. The native baseline therefore uses a dynamic
convex hull built from the moving mesh vertices. This is acceptable for the
current probe because the moving mesh is convex-like, but future case studies
with a genuinely nonconvex moving body should use either native convex
decomposition or multiple native convex hulls.
