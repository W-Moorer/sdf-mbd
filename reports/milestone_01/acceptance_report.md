# Milestone 01 - Acceptance Report

## Date

2026-04-22

## Summary

The Chrono baseline contact demo has been successfully built, run, and validated.
All acceptance criteria have been met.

## Acceptance Criteria

| # | Criterion | Status | Notes |
|---|-----------|--------|-------|
| 1 | Project builds successfully | PASS | CMake config + MSBuild compile, zero errors |
| 2 | Simulation program runs to completion | PASS | Exit code 0, 4001 steps completed |
| 3 | Baseline scenario is reproducible | PASS | Two runs produce identical CSV output (verified by hash comparison) |
| 4 | Position and contact data are output | PASS | CSV includes: time, pos(x/y/z), vel(x/y/z), contact_count, contact_force(x/y/z), contact_torque(x/y/z) |
| 5 | No divergence, crash, or empty run | PASS | Stable physics, contact force matches gravity at equilibrium (5136.5 N) |
| 6 | Next-step OpenVDB insertion points documented | PASS | See "Next Steps" below |

## Test Results

### Run 1
- Exit code: 0
- Steps: 4001
- Output rows: 201
- Final Y position: 0.4901 m
- Final contact force Y: 5136.51 N

### Run 2 (reproducibility check)
- Exit code: 0
- Steps: 4001
- Output rows: 201
- Final Y position: 0.4901 m
- Final contact force Y: 5136.51 N

### Hash comparison
- Files are IDENTICAL

## Next Steps (OpenVDB Integration)

The recommended insertion point for SDF-based contact is the narrowphase collision
detection layer:

1. Add a new collision shape type (`ChCollisionShapeSDF`) in `src/chrono/collision/`
2. Add a new narrowphase algorithm enum value in `ChNarrowphase::Algorithm`
3. Implement SDF query logic using OpenVDB to replace penetration depth and normal
   computation when both contacting bodies use SDF shapes
4. Integrate with existing SMC contact force computation pipeline

## Conclusion

Milestone 01 is **ACCEPTED**. The repository is in a clean, stable state ready for
the next development phase.
