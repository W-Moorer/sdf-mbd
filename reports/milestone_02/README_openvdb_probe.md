# Milestone 2A - OpenVDB SDF Probe Demo

## Overview

This milestone verifies the OpenVDB integration with Chrono at the asset and query layers.
It does NOT implement contact forces or modify the contact pipeline.

The demo validates:
1. OpenVDB builds and links successfully in the project
2. A Chrono body can be bound to an SDF asset
3. World -> body local -> SDF local coordinate transformation is correct
4. Phi (signed distance) and gradient queries work correctly
5. Results match analytical sphere distance within numerical precision

## Build

```powershell
cmake -B build -S . -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target demo_CH_openvdb_probe -j 8
```

Dependencies:
- OpenVDB (via vcpkg)
- TBB (via vcpkg)
- Imath (via vcpkg)

## Run

```powershell
# Default output to out/milestone_02/sdf_probe_output.csv
build\bin\Release\demo_CH_openvdb_probe.exe

# Custom output path
build\bin\Release\demo_CH_openvdb_probe.exe custom/path.csv
```

## Output

CSV file at `out/milestone_02/sdf_probe_output.csv` containing:
- point_id: query point index
- world_x/y/z: world-space coordinates
- local_x/y/z: body-local (SDF) coordinates
- sdf_phi: OpenVDB signed distance query
- analytic_phi: analytical sphere signed distance
- abs_error: |sdf_phi - analytic_phi|
- in_band: whether the point is within the narrow band
- grad_x/y/z: finite-difference gradient
- analytic_nx/ny/nz: analytical surface normal
- grad_dot_normal: dot product between gradient and analytical normal

## Results Summary

- Points probed: 36
- Points in band: 33
- Max |error| (in band): ~1.8e-8 m
- Mean |error| (in band): ~1.4e-8 m
- Gradient-normal alignment: 0.97

## Architecture

### SDFAsset
A lightweight wrapper around `openvdb::FloatGrid` that tracks:
- The OpenVDB grid pointer
- Voxel size and narrow band half width
- Asset name

### BodySDFBinding
Binds a Chrono `ChBody` to an `SDFAsset`. The SDF is defined in the body's
local coordinate frame. Coordinate transformation:
- `WorldToSDFLocal`: uses `ChFrame::TransformPointParentToLocal` to transform
  world coordinates into the body's reference frame, which is the SDF's domain.

### Query Methods
- `ProbePhi`: trilinear interpolation in OpenVDB index space via `GridSampler`
- `ProbeGradient`: central finite difference using world-space probes
