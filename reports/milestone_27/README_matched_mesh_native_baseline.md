# Milestone 27: matched native mesh baseline

## Goal

This milestone removes the main geometry-mismatch risk in the previous native
baseline. Each scenario now builds one triangle mesh and feeds it to both:

- Chrono native SMC collision through `ChCollisionShapeTriangleMesh`.
- Field-contact primitive evaluation through an OpenVDB level-set SDF generated
  from the same triangle mesh.

## Scenarios

| Scenario | Target mesh | Moving body | Vertices | Faces | OpenVDB active voxels |
|---|---:|---:|---:|---:|---:|
| `matched_u_channel_sphere_slide` | U-channel / guide rail | sphere, r=0.12 | 396 | 772 | 721588 |
| `matched_staggered_pad_track` | staggered pad track | sphere, r=0.12 | 546 | 1000 | 309739 |

The native mesh collision radius is zero in both cases, so both backends share
the same explicit target geometry.

## Output Files

- `out/milestone_27/matched_mesh_field_frames.csv`
- `out/milestone_27/matched_mesh_native_frames.csv`
- `out/milestone_27/matched_mesh_comparison.csv`
- `out/milestone_27/matched_mesh_geometry.csv`
- `out/milestone_27/matched_mesh_regression_summary.json`

The frame CSVs include displacement, velocity, acceleration, force, torque, and
contact/patch count. Native frames additionally record actual-vs-commanded
position, velocity, and acceleration errors.

## Key Results

| Scenario | Field mean patch churn | Native mean contact churn | Churn ratio | Field RMS normal jump | Native RMS contact-normal jump | Normal/jump ratio |
|---|---:|---:|---:|---:|---:|---:|
| `matched_u_channel_sphere_slide` | 0.032000 | 2.472000 | 77.25 | 0.073788 | 0.592200 | 8.03 |
| `matched_staggered_pad_track` | 0.025926 | 0.842593 | 32.50 | 0.054315 | 0.105188 | 1.94 |

Tracking remains controlled under native SMC:

- Max position error: `2.28823e-3`.
- Max velocity error: `8.36e-6`.
- Max acceleration error: `3.629392e-2`.
- Field Coulomb ratio: `1.0`.

## Regression Criteria

`scripts/check_matched_mesh_baseline.py` now checks:

- At least two matched-mesh scenarios are present.
- Native contact-count churn is more than `10x` field patch-count churn.
- Native contact-normal jump is more than `1.5x` field normal jump.
- Field Coulomb ratio stays within `1 + 1e-6`.
- Native trajectory errors remain bounded:
  - position <= `3.5e-3`
  - velocity <= `1e-3`
  - acceleration <= `1e-1`

CTest targets:

- `demo_CH_field_contact_matched_mesh_native_baseline`
- `demo_CH_field_contact_matched_mesh_regression`

Current result: both tests pass.

## Interpretation

The comparison is now geometry-matched at the triangle mesh level. The result
supports the paper claim that field primitives reduce discrete contact-set
churn and smooth the response-normal evolution, without relying on an analytic
height field or a native box/cylinder surrogate. The remaining distinction is
backend representation: Chrono native uses mesh collision/contact generation,
while the field method uses the OpenVDB level-set generated from the same mesh.
