# Milestone 2A - Acceptance Report

## Date

2026-04-22

## Summary

The OpenVDB SDF probe demo has been successfully built, run, and validated against
analytical sphere geometry. All acceptance criteria have been met.

## Acceptance Criteria

| # | Criterion | Status | Notes |
|---|-----------|--------|-------|
| 1 | Project builds and links OpenVDB code | PASS | CMake finds OpenVDB via vcpkg, MSBuild links with Imath + TBB, zero errors |
| 2 | demo_CH_openvdb_probe runs successfully | PASS | Exit code 0, 36 points probed |
| 3 | Output file generated at out/milestone_02/sdf_probe_output.csv | PASS | CSV with 37 columns, 36 data rows + header |
| 4 | Output includes world coords, local coords, phi, grad | PASS | All fields present and populated |
| 5 | Analytical distance comparison for sphere geometry | PASS | abs_error computed for each point; max error ~1.8e-8 m in narrow band |
| 6 | Error results documented in acceptance report | PASS | See below for detailed analysis |
| 7 | File organization clean, no scattered temp files | PASS | All files in designated directories |

## Detailed Results

### Query Setup

- Geometry: Sphere, radius = 1.0 m
- Voxel size: 0.05 m
- Narrow band half width: 3 voxels (0.15 m)
- Body position: (2.0, 1.0, 0.5) m
- Body rotation: 0.3 rad around Y axis
- Total probe points: 36
- Points within narrow band: 33

### Accuracy

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Max |error| (in band) | 1.79e-8 m | < 0.02 m | PASS |
| Mean |error| (in band) | 1.44e-8 m | < 0.005 m | PASS |
| Avg |grad dot normal| | 0.970 | > 0.95 | PASS |

### Coordinate Transformation Verification

The body is placed at a non-trivial pose (translated + rotated) to verify that
coordinate transforms work correctly:

1. **World-to-local transform**: Uses Chrono's `ChFrame::TransformPointParentToLocal`
   to convert world-space query points into the body's reference frame.

2. **Validation**: Probe points on the exact sphere surface (in local coordinates:
   distance from origin = radius) yield `sdf_phi ≈ 0`, confirming the coordinate
   transform is correct. For example:
   - Surface point along body-local X axis: phi = -1.49e-8 (essentially zero)
   - Surface point along body-local Y axis: phi = -1.49e-8 (essentially zero)
   - Surface point along body-local Z axis: phi = -1.49e-8 (essentially zero)

3. **Gradient alignment**: The finite-difference gradient aligns with the analytical
   surface normal at 0.970 average dot product (close to 1.0), confirming that both
   the SDF query and the coordinate transform are correct.

### Narrow Band Behavior

As expected for an OpenVDB level set:
- Points within the narrow band (|phi| < 0.15 m): accurate signed distance
- Points outside the narrow band: values clipped to +/- 0.15 m
- This is standard Level Set behavior, not an error

### Error Analysis

The extremely small errors (~10 nanometers) are consistent with:
- OpenVDB's trilinear interpolation precision
- Single-precision floating-point representation in the FloatGrid
- Voxel size of 0.05 m (interpolation error expected to be O(voxel_size^2))

## Observations

1. The vcpkg OpenVDB package is built as a dynamic library requiring runtime DLLs
   (openvdb.dll, Imath-3_2.dll, tbb12.dll) to be in the PATH.
2. The `openvdb::tools::createLevelSetSphere` function is the simplest way to
   create SDF assets for testing.
3. The narrow-band limitation is inherent to level set representation and must
   be considered in future contact implementations.

## Conclusion

Milestone 2A is **ACCEPTED**. The OpenVDB integration at the asset and query layers
is functional, accurate, and ready for the next development phase.
