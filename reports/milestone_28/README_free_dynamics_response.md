# Milestone 28: free dynamics response

## Goal

This milestone moves beyond controlled kinematic-path evaluation. The moving
sphere has real mass and inertia, and its trajectory is integrated by Chrono.
The two compared variants use the same external drive:

- `native_mesh_smc`: Chrono native SMC contact against the matched triangle mesh.
- `field_primitive_accumulator`: native collision disabled; field primitive
  force and torque are accumulated on the moving `ChBody` through
  `AccumulateForce` and `AccumulateTorque`.

The external drive is force based, not a prescribed path:

- x direction: velocity servo force.
- y direction: normal preload plus damping.
- z direction: weak guide force for rail/pad traversal.

## Scenarios

| Scenario | Target mesh | Moving body | Steps | dt |
|---|---|---:|---:|---:|
| `free_u_channel_sphere_slide` | U-channel / guide rail | sphere, r=0.12 | 2400 | 5e-4 |
| `free_staggered_pad_track` | staggered pad track | sphere, r=0.12 | 2500 | 5e-4 |

## Output Files

- `out/milestone_28/free_dynamics_frames.csv`
- `out/milestone_28/free_dynamics_summary.csv`
- `out/milestone_28/free_dynamics_comparison.csv`
- `out/milestone_28/free_dynamics_geometry.csv`
- `out/milestone_28/free_dynamics_regression_summary.json`

The frame CSV records position, velocity, acceleration, external force, contact
force/torque, kinetic energy, external/contact work, penetration envelope,
Coulomb ratio, and stick/slip counts.

## Key Results

| Scenario | Field max K | Native max K | Field max penetration | Native max penetration | Field/native accel-jump ratio |
|---|---:|---:|---:|---:|---:|
| `free_u_channel_sphere_slide` | 0.492741 | 0.651347 | 0.017531 | 0.062407 | 0.0164 |
| `free_staggered_pad_track` | 0.781244 | 0.856195 | 0.015563 | 0.029137 | 0.0444 |

Additional checks:

- Field Coulomb ratio: `1.0`.
- Native Coulomb ratio: `0.99999997`.
- Field net contact work is negative in both scenarios:
  - U-channel: `-0.738224`
  - staggered pad: `-1.637083`
- Final trajectory deltas against native:
  - U-channel: `0.038391`
  - staggered pad: `0.051705`
- Field stick/slip switches are much less chattery:
  - U-channel: `9` vs native `465`
  - staggered pad: `27` vs native `658`

## Regression Criteria

`scripts/check_free_dynamics_response.py` checks:

- At least two free-dynamics scenarios are present.
- No nonphysical energy explosion:
  - max field/native energy budget ratio <= `0.75`
- Coulomb bound for both variants:
  - max ratio <= `1 + 1e-6`
- Penetration bounded:
  - field max penetration <= `0.020`
  - field/native max penetration <= `0.070`
- Field trajectory smoothness is not worse than native:
  - field/native RMS acceleration-jump ratio <= `0.15`
- Final trajectory remains comparable:
  - final position delta <= `0.080`
- Field contact net work is non-injecting:
  - net contact work <= `0.02`
- Field contact remains active:
  - active-frame ratio >= `0.95`
- Field stick-slip sequence is less chattery:
  - field/native switch ratio <= `0.10`

CTest targets:

- `demo_CH_field_contact_free_dynamics_response`
- `demo_CH_field_contact_free_dynamics_regression`

Current result: both tests pass.

## Interpretation

The field primitive model is now exercised as a force-producing contact law in a
Chrono dynamics loop. The moving body is not prescribed along the target path:
contact forces, torque, external drive, and inertia jointly determine the
trajectory. Under the matched mesh/SDF geometry, the field accumulator response
keeps energy bounded, respects the Coulomb cap, bounds penetration, and produces
substantially smoother acceleration and stick/slip timelines than the native
mesh-contact baseline.
