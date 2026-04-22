# Milestone 2B - One-way Pointwise SDF Penalty Contact

## Overview

This milestone extends the verified OpenVDB query layer (Milestone 2A) into a minimal
one-way pointwise penalty contact loop. A dynamic Chrono body carries surface sample
points that are queried against a static SDF target each timestep.

## Build

```powershell
cmake -B build -S . -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target demo_CH_sdf_point_contact -j 8
```

## Run

```powershell
build\bin\Release\demo_CH_sdf_point_contact.exe
```

## Output

CSV file at `out/milestone_03/sdf_point_contact_output.csv`:
- time, position (x/y/z), velocity (x/y/z)
- active_sample_count, min_phi, mean_active_phi
- total_force (x/y/z), total_torque (x/y/z)

## Scenario

- **Static SDF target**: Large sphere (radius=100m) with surface at y=0, acting as a flat ground
- **Dynamic body**: Small sphere (radius=0.2m, mass=33.5kg) dropped from 2m height
- **Sample points**: 40 points (5x8 theta-phi) on the dynamic sphere surface
- **Penalty method**: f = (k * penetration + c * damping_vn) * n for active points (phi < 0.3m)

## Parameters

| Parameter | Value |
|-----------|-------|
| Penalty stiffness | 500,000 N/m |
| Penalty damping | 500 N*s/m |
| Activation delta | 0.3 m |
| SDF voxel size | 0.05 m |
| SDF band half width | 20 voxels (1.0 m) |
| Time step | 0.0005 s |
| Total simulation time | 2.0 s |

## Architecture

The demo reuses the SDFAsset and query infrastructure from Milestone 2A:
- `SDFAsset`: wraps `openvdb::FloatGrid` (SDF storage)
- `SDFContactTarget`: wraps SDFAsset with world frame transform + probe methods
- `SurfaceSample`: local position + area weight on dynamic body
- Simulation loop: transform samples → probe SDF → compute penalty forces → apply to body via `ChBody::AccumulateForce`/`AccumulateTorque`
