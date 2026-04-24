# Milestone 24: Mesh-to-OpenVDB SDF Non-Convex Experiments

This milestone moves the field-contact primitive experiments from analytic
height fields to a real mesh SDF backend.

## What Was Added

New benchmark:

```text
src/demos/core/demo_CH_field_contact_mesh_sdf_nonconvex_openvdb.cpp
```

The benchmark builds explicit triangle meshes, converts them to OpenVDB
level-set SDFs with `meshToLevelSet`, and reuses the existing
`FieldContactPrimitiveTracker` runtime.

Scenes:

- `u_channel_sphere_slide`
- `toothed_rail_sphere_slide`
- `key_slot_insertion`
- `staggered_pad_track`

Variants:

- `field_primitive`
- `raw_sample_contacts`
- `normal_bin_fragments`
- `convex_decomposition_proxy`
- `chrono_traditional_proxy`

The figure generator now accepts a CSV prefix:

```text
scripts/generate_field_contact_figures.py --prefix mesh_sdf_nonconvex
```

## Commands

```powershell
cmake --build build --config Release --target demo_CH_field_contact_mesh_sdf_nonconvex_openvdb
.\build\bin\Release\demo_CH_field_contact_mesh_sdf_nonconvex_openvdb.exe

python scripts\generate_field_contact_figures.py `
  --input out\milestone_24 `
  --output out\milestone_24\figures `
  --prefix mesh_sdf_nonconvex
```

Outputs:

```text
out/milestone_24/mesh_sdf_nonconvex_frames.csv
out/milestone_24/mesh_sdf_nonconvex_patches.csv
out/milestone_24/mesh_sdf_nonconvex_summary.csv
out/milestone_24/mesh_sdf_nonconvex_comparison.csv
out/milestone_24/mesh_sdf_nonconvex_geometry.csv
out/milestone_24/figures/index.html
```

The current run generated 29 SVG figures.

## Geometry Backend

The static target is no longer an analytic height function.  The current
triangle mesh and OpenVDB grid sizes are:

```text
u_channel_sphere_slide:
  target_vertices=396
  target_faces=772
  active_voxels=721588

toothed_rail_sphere_slide:
  target_vertices=320
  target_faces=584
  active_voxels=391247

key_slot_insertion:
  target_vertices=386
  target_faces=744
  active_voxels=470417

staggered_pad_track:
  target_vertices=546
  target_faces=1000
  active_voxels=309739
```

## Current Results

All variants preserve the Coulomb bound:

```text
max_tangential_force_ratio=1.00000000
```

Topology churn is strongly reduced by `field_primitive`.  Ratios below are
baseline divided by `field_primitive`; values above 1 mean the baseline is worse.

```text
u_channel_sphere_slide:
  raw_sample_contacts mean_patch_count_change_ratio=44.20000000
  normal_bin_fragments mean_patch_count_change_ratio=22.80000000
  chrono_traditional_proxy rms_normal_jump_ratio=1.71848414

toothed_rail_sphere_slide:
  raw_sample_contacts mean_patch_count_change_ratio=15.29166667
  normal_bin_fragments rms_normal_jump_ratio=1.53268748
  chrono_traditional_proxy force_oscillation_ratio=1.78887089

key_slot_insertion:
  raw_sample_contacts mean_patch_count_change_ratio=95.00000000
  normal_bin_fragments mean_patch_count_change_ratio=131.00000000
  chrono_traditional_proxy topology_event_ratio=759.60000000

staggered_pad_track:
  raw_sample_contacts mean_patch_count_change_ratio=19.16666667
  normal_bin_fragments rms_normal_jump_ratio=1.72204848
  chrono_traditional_proxy torque_oscillation_ratio=2.29319977
```

Force and torque oscillation are scene-dependent.  Some fragmented proxy
variants have lower force oscillation because their point/contact selection
under-integrates or smooths away parts of the distributed response.  The robust
claim from this milestone is therefore:

```text
field primitives keep contact collections compact and persistent on real
mesh-derived SDF geometry, while preserving Coulomb feasibility and inherited
history non-amplification.
```

## Meaning For The Paper

This milestone closes the main gap left by the analytic height-field demos:
the primitive runtime is now exercised against an OpenVDB SDF generated from
explicit non-convex triangle meshes.

It supports the TeX claim that the method is a field/SDF primitive formulation,
not a construction that only works for hand-written analytic height functions.

The next paper-strengthening step is to pair this backend with a matched native
collision baseline using the same target triangle mesh, so force/torque claims
can be made without the geometry mismatch caveat from Milestone 23.
