# Milestone 15: Real Split/Merge Primitive Demo

This milestone adds a real split/merge test for the field-based contact
primitive implementation.

## Executable

`demo_CH_field_contact_split_merge_openvdb`

## Scene

- Target field: OpenVDB narrow-band SDF built from two overlapping triangulated
  sphere meshes.
- Moving surface: triangulated plane surface graph with vertex area weights and
  edge adjacency.
- Motion: the plane height follows a down-up trajectory. At high position the
  active contact set has two disconnected components; at lower position the two
  components merge; when the plane rises again, the merged component splits.

This produces a true runtime topology sequence:

`2 patches -> 1 patch -> 2 patches`

The event is caused by active-set connectivity on the surface graph, not by
manual event injection.

## Outputs

- `out/milestone_15/field_contact_split_merge_frames.csv`
- `out/milestone_15/field_contact_split_merge_patches.csv`
- `out/milestone_15/field_contact_split_merge_summary.csv`

## Current Results

- frames: 800
- plane vertices: 5376
- plane faces: 10450
- target vertices: 12324
- target faces: 24640
- target SDF active voxels: 340649
- frames with one patch: 367
- frames with two patches: 433
- frames with merge: 1
- frames with split: 1
- total newborn primitives: 2
- max patch count: 2
- max source count: 2
- max previous reuse: 2
- max tangential force ratio: 1.0

Event frames:

- Frame 258: merge. One current primitive inherits from two previous sources:
  `0:0.4994;1:0.4945`.
- Frame 625: split. Two current primitives inherit from the previous merged
  primitive with weight about `0.4938` each.

Energy checks:

- Merge inherited-energy ratio: about `0.9779`.
- Split inherited-energy ratio: about `0.4938` for each child.

These ratios are relative to the weighted source-history energy bound. Values
below one verify that overlap inheritance does not amplify stored tangential
history through the split/merge transition.

## Interpretation

This is the first demo where the TeX paper's persistent primitive history model
is exercised by an actual topology change. The implementation now demonstrates:

- active-set connectivity can naturally create split/merge primitive events;
- merge can aggregate two previous tangential histories;
- split can let multiple children inherit from one previous primitive while
  assigning unique current persistent IDs;
- inherited history weights remain bounded;
- tangential force remains inside the Coulomb disk.

## Remaining Work

- Add external mesh asset loading instead of generated meshes.
- Add bidirectional split/merge de-duplication across both bodies.
- Convert this kinematic demo into a Chrono body-pair force accumulation demo.
- Add automated regression tests for objectivity, Coulomb feasibility,
  non-amplification, and split/merge event classification.
