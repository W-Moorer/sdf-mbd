# Milestone 20: Feature-Switching Sensitive Scenes

This milestone adds a deterministic analytic demo for the first non-smooth and
non-convex feature-switching experiments from the TeX design.

## Demo

Target:

`src/demos/core/demo_CH_field_contact_feature_switching.cpp`

Build and run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_feature_switching
.\build\bin\Release\demo_CH_field_contact_feature_switching.exe
```

Outputs:

```text
out/milestone_20/field_contact_feature_switching_frames.csv
out/milestone_20/field_contact_feature_switching_patches.csv
out/milestone_20/field_contact_feature_switching_summary.csv
```

## Scenes

The SDF backend is an analytic height field, so the experiment isolates contact
primitive topology and tangential history handling rather than OpenVDB meshing
or voxelization.

- `folded_seam`: a sphere slides across a folded polyline seam with a sharp
  normal switch.
- `concave_polyline_corner`: a sphere crosses a local concave polyline notch.
- `narrow_groove_entrance`: a sphere traverses a wide-groove to narrow-neck to
  wide-groove transition, producing real primitive merge and split events.

## Compared Variants

- `current_field_primitive`: the current runtime tracker.
- `no_history`: topology is still tracked, but tangential elastic history is
  discarded.
- `direct_inherit`: a deliberately naive baseline that directly sums inherited
  source histories without source weights or transport.
- `projection_transport`: weighted source inheritance with simple tangent-plane
  projection.
- `minimal_rotation_gate_split_merge`: weighted source inheritance with minimal
  rotation transport, gate, and split/merge aggregation.

## Current Run

Key summary values from the current run:

```text
folded_seam,current_field_primitive:
  active_frames=601
  primitive_persistence_ratio=0.99833611
  max_normal_jump_angle=0.01859210
  force_oscillation_index=0.00372532
  max_tangential_force_ratio=1.00000000

concave_polyline_corner,current_field_primitive:
  active_frames=570
  primitive_persistence_ratio=0.99649123
  max_normal_jump_angle=0.00241181
  force_oscillation_index=0.00395893
  max_tangential_force_ratio=1.00000000

narrow_groove_entrance,current_field_primitive:
  active_frames=428
  max_patch_count=3
  total_merge_patches=43
  total_split_patches=3
  primitive_persistence_ratio=0.99458484
  max_normal_jump_angle=0.35593444
  force_oscillation_index=0.00504809
  max_tangential_force_ratio=1.00000000
  max_inherited_energy_ratio=1.00000000
```

The narrow-groove scene is the primary Milestone 20 discriminator. It drives the
primitive count through `0/1/2/3`, triggers real merge and split classification,
and exposes the naive direct-inheritance failure:

```text
narrow_groove_entrance,direct_inherit:
  max_inherited_energy_ratio=1.99618868

narrow_groove_entrance,projection_transport:
  max_inherited_energy_ratio=1.00000000

narrow_groove_entrance,minimal_rotation_gate_split_merge:
  max_inherited_energy_ratio=1.00000000
```

This fixes a concrete CI-checkable interpretation: direct source summation can
nearly double inherited elastic energy at merge, while the weighted projection
and minimal-rotation split/merge paths preserve the non-amplification bound.

## Interpretation

The runtime field primitive and the local `minimal_rotation_gate_split_merge`
ablation agree on these scenarios, which is expected because both use the same
source weighting, minimal-rotation transport, and Coulomb return mapping.

`projection_transport` is intentionally close to minimal rotation in these
mostly translational height-field scenes; Milestone 17 remains the sharper
objectivity discriminator for projection versus minimal rotation. Milestone 20
instead establishes the first feature-switching benchmark where non-convex
topology changes, split/merge classification, Coulomb feasibility, and history
non-amplification are all visible in one repeatable run.
