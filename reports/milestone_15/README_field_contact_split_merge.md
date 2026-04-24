# Milestone 15: Real Split/Merge Primitive Demo

This milestone adds a real split/merge test for the field-based contact
primitive implementation.

## Executable

`demo_CH_field_contact_split_merge_openvdb`

## Scene

- Target field: OpenVDB narrow-band SDF built from two overlapping triangulated
  sphere meshes.
- Moving body: a Chrono `ChBody` carrying a triangulated plane surface graph
  with vertex area weights and edge adjacency.
- Force coupling: primitive normal/tangential force and torque are accumulated
  into the moving body with `AccumulateForce` and `AccumulateTorque`.
- Motion: a low-frequency PD guide force drives the body through a down-up
  trajectory, while the actual body position and velocity are integrated by
  Chrono and perturbed by the primitive contact force.

This produces a true runtime topology sequence:

`2 patches -> 1 patch -> 2 patches`

The event is caused by active-set connectivity on the surface graph, not by
manual event injection.

## Outputs

- `out/milestone_15/field_contact_split_merge_frames.csv`
- `out/milestone_15/field_contact_split_merge_patches.csv`
- `out/milestone_15/field_contact_split_merge_summary.csv`

## Current Results

- frames: 1600
- plane vertices: 5376
- plane faces: 10450
- target vertices: 12324
- target faces: 24640
- target SDF active voxels: 340649
- frames with one patch: 716
- frames with two patches: 884
- frames with merge: 1
- frames with split: 1
- total newborn primitives: 2
- max patch count: 2
- max source count: 2
- max previous reuse: 2
- max tangential force ratio: 1.0
- max tracking error: about 0.007 m
- max contact force: about 5.62e4 N
- max guide force: about 5.64e4 N

Event frames:

- Frame 524: merge. One current primitive inherits from two previous sources:
  `1:0.4988;0:0.4951`.
- Frame 1240: split. Two current primitives inherit from the previous merged
  primitive with weights about `0.5144` and `0.4782`.

Energy checks:

- Merge inherited-energy ratio: about `0.9937`.
- Split inherited-energy ratios: about `0.5144` and `0.4782`.

These ratios are relative to the weighted source-history energy bound. Values
below one verify that overlap inheritance does not amplify stored tangential
history through the split/merge transition.

## Interpretation

This is the first demo where the TeX paper's persistent primitive history model
is exercised by an actual topology change inside a Chrono force-accumulation
loop. The implementation now demonstrates:

- active-set connectivity can naturally create split/merge primitive events;
- merge can aggregate two previous tangential histories;
- split can let multiple children inherit from one previous primitive while
  assigning unique current persistent IDs;
- inherited history weights remain bounded;
- tangential force remains inside the Coulomb disk.
- primitive contact force and torque enter the rigid-body dynamics through the
  same accumulator path used by the earlier SDF contact demos.

## Regression Checks

`demo_CH_field_contact_regression` fixes the key primitive invariants as
CI-checkable pass/fail metrics without requiring OpenVDB or GoogleTest:

- Coulomb feasibility: the tangential update must project force into the
  Coulomb disk and keep it tangent to the contact normal.
- History non-amplification: aggregated inherited tangential energy must not
  exceed the weighted source-history energy bound.
- Merge classification: two previous primitives merging into one current
  primitive must be detected with bounded source weights.
- Split classification: two current children inheriting from one previous
  primitive must be detected without total inherited source weight exceeding
  one.

Build and run directly:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
.\build\bin\Release\demo_CH_field_contact_regression.exe
```

CTest entry when `CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS=ON`:

```powershell
ctest --test-dir build -C Release -R demo_CH_field_contact_regression --output-on-failure
```

A mirrored GoogleTest file also exists at
`src/tests/unit_tests/collision/utest_COLL_field_contact_primitives.cpp`, but
this workspace currently disables unit-test targets because the googletest
submodule is unavailable. The standalone regression executable is therefore the
robust CI gate for the current checkout.

## Remaining Work

- Add external mesh asset loading instead of generated meshes.
- Add bidirectional split/merge de-duplication across both bodies.
- Replace the guide force with a controlled actuator or an unconstrained
  benchmark once the contact law is ready for full dynamics validation.
- Add automated objectivity regression once the transport/frame benchmark is
  promoted from demo diagnostics to a stable numeric fixture.
