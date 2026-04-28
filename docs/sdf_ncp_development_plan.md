# SDF-NCP Development Plan

## Repository audit

### Current repository structure

This repository is a Project Chrono based C++ tree with local field-contact
extensions and paper reproduction scripts.

- `src/chrono/`: Chrono core source tree, including collision, physics,
  solver, timestepper, geometry, FEA, vehicle, sensor, and utility modules.
- `src/chrono/collision/`: Collision shapes, Bullet integration, and the
  field-contact prototype headers.
- `src/demos/core/`: Core Chrono demos, including field contact, SDF point
  contact, OpenVDB SDF contact, and paper-style experiments.
- `src/tests/unit_tests/`: C++ and Python unit tests registered through the
  Chrono CMake test structure.
- `scripts/`: Python regression, comparison, and plotting utilities for field
  contact and paper figures.
- `paper_example/`: OpenVDB/SDF paper example inputs, C++ runner, check script,
  plotting script, and generated PDF figures.
- `assets/`: Geometry, mesh, and reference data for existing field-contact
  examples.
- `paper/`, `reports/`, `out/`: Paper drafts and generated experiment outputs.

There is no existing Python package named `sdf_mbd`, no `pyproject.toml`, and
no root-level `pytest` configuration at the time of this audit.

### Existing SDF-related modules

- `src/chrono/collision/ChFieldContactPrimitives.h` contains SDF/field-query
  oriented data structures and utilities. The key compatibility structure is
  `chrono::fieldcontact::FieldSampleQuery` with `phi`, `grad`, `world_pos`, and
  `world_vel`.
- `src/chrono/collision/ChFieldContactRuntime.h` turns per-sample field queries
  into persistent field-contact patch responses.
- `src/demos/core/demo_CH_sdf_point_contact.cpp` includes an analytical plane
  SDF path and a penalty-style point/surface force calculation.
- `src/demos/core/*openvdb*.cpp` contains OpenVDB SDF demos and mesh-based
  sparse SDF examples.
- Bullet also contains `cbtMiniSDF` and `cbtSdfCollisionShape`, but those are
  third-party collision-shape support code and are not the best extension point
  for this first NCP prototype.

### Existing multibody dynamics modules

- Chrono rigid-body dynamics are under `src/chrono/physics`, with systems such
  as `ChSystemSMC`, bodies such as `ChBody` and `ChBodyEasySphere`, and existing
  solvers under `src/chrono/solver`.
- Existing C++ demos assemble simple multibody simulations in `src/demos/core`
  and more specialized demos in other `src/demos/*` folders.
- The current field-contact demos use Chrono time stepping and apply penalty
  contact forces through body force accumulators.

### Existing contact or collision modules

- Classical collision shape infrastructure is in `src/chrono/collision`.
- Field-contact primitive and runtime logic is in
  `ChFieldContactPrimitives.h` and `ChFieldContactRuntime.h`.
- Existing field-contact CTest registrations live in `src/demos/core/CMakeLists.txt`.
- Existing field-contact unit tests are in
  `src/tests/unit_tests/collision/utest_COLL_field_contact_primitives.cpp`.

### Existing examples or scripts that can be extended

- `src/demos/core/demo_CH_sdf_point_contact.cpp` is the closest C++ reference
  for analytical plane SDF contact and penalty-force behavior.
- `scripts/generate_field_contact_figures.py` and related `scripts/check_*.py`
  files show the repository's pattern for reproducible CSV-driven figures.
- `paper_example/plot_paper_examples.py` and `paper_example/check_paper_examples.py`
  show paper reproduction patterns.

### Planned SDF-NCP addition points

The first SDF-NCP paper prototype will be added as a small Python research layer:

- `sdf_mbd/sdf/primitives.py`: analytical SDF primitives and compatibility
  helpers for `phi(x)` and `grad(x)`.
- `sdf_mbd/contact/ncp.py`: smoothed Fischer-Burmeister residual, gradient, and
  complementarity diagnostics.
- `sdf_mbd/contact/sdf_contact.py`: SDF gap, normal, point contact Jacobian, and
  finite-difference gradient checks.
- `sdf_mbd/contact/penalty.py`: penalty-contact baseline from SDF gap and
  normal.
- `sdf_mbd/solvers/point_mass.py`: implicit Euler SDF-NCP stepper for a 2D point
  mass contacting an SDF plane, plus a penalty baseline stepper.
- `examples/sdf_ncp/`: reproducible scripts for point-mass contact, penalty
  sensitivity, and SDF geometry visualization.
- `tests/`: focused pytest tests for the Python prototype.

This avoids destabilizing the existing Chrono C++ tree while providing a
reproducible, testable NCP foundation that can later be ported into Chrono C++
contact assembly.

### Files to avoid touching unless necessary

- Existing Chrono core source under `src/chrono/physics`, `src/chrono/solver`,
  `src/chrono/timestepper`, and existing collision-shape implementations.
- Existing OpenVDB demo files under `src/demos/core/*openvdb*.cpp`.
- Generated outputs under `build/`, `out/`, `paper/generated/`, and existing
  paper figures unless an experiment explicitly writes new SDF-NCP results.
- Existing `paper_example/` benchmark inputs and generated figures.
- Third-party code under `src/chrono_thirdparty/`, `openvdb-master/`, and Bullet
  internals.

### Risks and assumptions

- The repository does not currently expose a lightweight Python package, so the
  prototype will add one at the repository root rather than modifying Chrono's
  C++ build.
- SciPy may or may not be installed. The point-mass NCP solver should prefer
  `scipy.optimize.root` when available and provide a small Newton fallback for
  the one-contact plane case.
- Matplotlib may be unavailable in minimal environments. Example scripts should
  report a clear import error if plotting dependencies are missing.
- Existing C++ tests depend on the local build and optional OpenVDB/TBB/Imath
  configuration; this prototype should not silently change those tests.
- The first implementation is frictionless and limited to analytical SDFs and a
  2D point mass. Rigid-body, patch, friction, and deformable-body integration are
  future work.
