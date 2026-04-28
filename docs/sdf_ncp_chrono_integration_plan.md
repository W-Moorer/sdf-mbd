# SDF-NCP Chrono Integration Plan

## 1. Current C++ Contact-Related Layout

The repository is a Project Chrono based C++ tree. The contact path relevant to
this integration is:

- `src/chrono/collision/`
  - Collision shapes, collision models, collision systems, and contact reports.
  - `ChCollisionInfo` carries contact normal and signed distance, where negative
    distance means penetration.
  - `ChFieldContactPrimitives.h` and `ChFieldContactRuntime.h` implement the
    current field-contact prototype.
- `src/chrono/physics/`
  - `ChSystem`, `ChSystemSMC`, `ChSystemNSC`
  - `ChContactContainer`, `ChContactContainerSMC`, `ChContactContainerNSC`
  - `ChContactSMC`, `ChContactNSC`
  - bodies, force accumulators, constraints, descriptors, and time stepping.
- `src/demos/core/`
  - Existing field-contact and SDF demos.
  - OpenVDB demos are guarded by `find_package(OpenVDB)`.
  - Non-OpenVDB field-contact regression demos are registered under the
    `field_contact` CTest label when `CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS`
    is enabled.
- `src/tests/unit_tests/collision/`
  - Existing GoogleTest coverage for field-contact primitive invariants.

The active local build has `CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS=ON` and
`BUILD_TESTING=OFF`, so the currently runnable `field_contact` CTest suite is
demo/regression based rather than the full GoogleTest tree.

## 2. Current `field_contact` Test Coverage

The current `ctest -L field_contact` run covers the registered demo/regression
tests in `src/demos/core/CMakeLists.txt`, including:

- field-contact primitive invariants,
- friction validation,
- feature switching and acceptance comparisons,
- matched mesh baselines,
- free dynamics response,
- CAD mesh case study,
- Chrono-original mesh baseline,
- OpenVDB paper examples and figure regeneration.

The C++ GoogleTest file
`src/tests/unit_tests/collision/utest_COLL_field_contact_primitives.cpp` covers
primitive-level invariants such as Coulomb feasibility, history aggregation,
merge/split classification, tangent transport, bidirectional pair conservation,
and force accumulation. In the current local build these tests are not in CTest
because `BUILD_TESTING=OFF`.

## 3. How Existing SDF Contact Enters the System

Existing SDF/field-contact code is penalty/force based:

- `ChFieldContactPrimitives.h` defines `FieldSampleQuery` with `phi`, `grad`,
  `world_pos`, and `world_vel`.
- Active samples are grouped into connected `PrimitivePatch` objects.
- `ComputeNormalContactIntegral` computes normal pressure as
  `stiffness * penetration + damping * closing_speed`, integrates force over
  sample area, and adds tangential friction through the runtime layer.
- Demos apply the resulting wrench through Chrono force accumulators or compare
  it to native Chrono contact paths.
- `ChContactSMC` is also penalty-force based.
- `ChContactNSC` is complementarity/constraint based in Chrono's native
  constraint solver, but it is built around collision-generated contact
  constraints and Coulomb tangential constraints, not the smooth
  Fischer-Burmeister residual proposed here.

Therefore, the current field-contact SDF path does not yet solve the
configuration-level SDF Signorini condition

```text
0 <= phi(x(q)) perpendicular lambda_n >= 0
```

with a smooth NCP residual.

## 4. Best First SDF-NCP Integration Point

The least intrusive first step is an independent C++ numerical module under
`src/chrono/collision/`, next to the field-contact headers:

```text
src/chrono/collision/ChSdfNcpContact.h
```

This keeps the first SDF-NCP implementation close to the current SDF/contact
code but avoids changing `ChSystem`, `ChContactContainerSMC`,
`ChContactContainerNSC`, or the Chrono solver interfaces.

The first CTest closure should be a standalone core demo/regression executable:

```text
src/demos/core/demo_CH_sdf_ncp_regression.cpp
```

registered with labels:

```text
collision;field_contact;sdf_ncp
```

This matches the current build, where `field_contact` regression tests are
already demo based.

## 5. Files That Can Be Reused

- `src/chrono/collision/ChFieldContactPrimitives.h`
  - SDF sign convention, safe normalization patterns, and field-query naming.
- `src/chrono/collision/ChFieldContactRuntime.h`
  - Runtime diagnostics and future patch-level extension patterns.
- `src/demos/core/CMakeLists.txt`
  - Existing demo and CTest registration style.
- `src/demos/core/demo_CH_field_contact_regression.cpp`
  - Standalone regression executable style with local `Check` helpers.
- Python reference prototype:
  - `sdf_mbd/contact/ncp.py`
  - `sdf_mbd/solvers/point_mass.py`

Note: the Python prototype requested in the original planning text as
`sdf_mbd/solvers/smooth_fb_solver.py` is implemented in this repository as
`sdf_mbd/solvers/point_mass.py`, with the NCP scalar functions in
`sdf_mbd/contact/ncp.py`.

## 6. Files To Avoid Modifying In This Phase

Do not modify unless a later integration stage requires it:

- `src/chrono/physics/ChSystem*.cpp/.h`
- `src/chrono/physics/ChContactContainer*.cpp/.h`
- `src/chrono/physics/ChContactSMC.cpp/.h`
- `src/chrono/physics/ChContactNSC.cpp/.h`
- Bullet internals under `src/chrono/collision/bullet/`
- OpenVDB demo implementations under `src/demos/core/*openvdb*.cpp`
- generated outputs under `build/`, `out/`, and existing `paper_example/figures/`

## 7. First-Stage C++ Minimal Implementation

The first C++ implementation should contain:

1. Smooth Fischer-Burmeister utility:

```text
Phi_eps(g, lambda) = sqrt(g*g + lambda*lambda + eps*eps) - g - lambda
```

with derivatives:

```text
dPhi_dg      = g / sqrt(g*g + lambda*lambda + eps*eps) - 1
dPhi_dlambda = lambda / sqrt(g*g + lambda*lambda + eps*eps) - 1
```

and diagnostic:

```text
max(0, -g) + max(0, -lambda) + abs(g * lambda)
```

2. A local point-mass/plane residual matching the Python prototype:

```text
q_next = q + dt * v_next
R_v = M(v_next - v) - dt(Q + J^T lambda)
R_lambda = Phi_eps(phi(q_next), lambda)
```

For a 2D point mass over the plane `y = 0`:

```text
R0 = m * (vx_next - vx)
R1 = m * (vy_next - vy) - dt * (-m*g0 + lambda_n)
R2 = Phi_eps(y + dt*vy_next, lambda_n)
```

3. A finite-difference Jacobian and local Newton solve, used only by the
   regression test. This must not alter Chrono's global solver path.

4. CTest entries:

```text
test_sdf_ncp_fb
test_sdf_ncp_point_mass_plane_step
test_sdf_ncp_point_mass_plane_rollout
```

with labels:

```text
collision;field_contact;sdf_ncp
```

## 8. Python/C++ Formula Alignment

Python reference:

- `sdf_mbd/contact/ncp.py`
- `sdf_mbd/solvers/point_mass.py`

C++ first-stage target:

- `src/chrono/collision/ChSdfNcpContact.h`
- `src/demos/core/demo_CH_sdf_ncp_regression.cpp`

Both implement the same one-contact frictionless equations:

```text
q_next = q + dt v_next
R_v = M(v_next - v) - dt(Q + J^T lambda)
R_lambda = Phi_eps(phi(q_next), lambda)
```

The C++ version is intentionally a minimal numerical kernel and regression
closure. It is not yet a replacement for `ChContactSMC`, `ChContactNSC`, or the
field-contact runtime.

## 9. 通用 SDF-NCP contact assembly 层设计

### Why `PointMassPlaneResidual` must be generalized

The current `PointMassPlaneResidual` is useful as a first regression closure, but
it hard-codes three assumptions:

- the SDF is the plane `phi(q) = y`,
- the contact normal and Jacobian are always `[0, 1]`,
- the state is a 2D point mass.

Those assumptions prevent reuse for later rigid-body, mesh-SDF, patch, and
field-contact assembly experiments. The next layer should accept a generic
surface query `phi(x)` and `grad phi(x)` plus a point kinematics map `x(q)`,
then assemble the same NCP equation without knowing the specific surface type.

### Proposed abstractions

The first generic layer remains intentionally small:

- SDF surface/provider:
  - `ChSdfSurface`
  - `Phi(x)` returns signed distance.
  - `Grad(x)` returns the SDF gradient.
  - First concrete implementation: `ChPlaneSdfSurface`.
  - Optional analytical helper: `ChSphereSdfSurface`.

- Contact kinematics:
  - `ChSdfPointKinematicsState`
  - `ChSdfPointKinematics`
  - First implementation is point translation only, with `x(q) = q` and
    `dx/dq = I`.
  - Rigid-body orientation Jacobians are explicitly deferred.

- Contact constraint:
  - `ChSdfNcpContactConstraint`
  - Stores gap, unit normal, point-position Jacobian row, normal multiplier,
    smoothed FB residual, penetration, and complementarity error.
  - For the point-translation case, `jacobian_position = normal`.

- Contact assembly/problem:
  - `ChSdfNcpPointMassState`
  - `ChSdfNcpPointMassSettings`
  - `ChSdfNcpPointMassResidual`
  - `ComputeSdfNcpPointMassResidual(surface, state, settings, z)`
  - `SolveSdfNcpPointMassStep(surface, state, settings)`

This layer is a small assembly kernel, not a new Chrono solver path.

### Why this phase still avoids `ChContactContainer`

Chrono's production contact path requires descriptor bookkeeping, active
constraint offsets, warm-start reaction transfer, solver coupling, and
time-stepper interactions. Connecting a nonlinear smooth FB residual directly
to those systems is a larger design decision. This phase avoids changing:

- `ChSystem`,
- `ChContactContainer`,
- `ChContactSMC`,
- `ChContactNSC`,
- global solver descriptor assembly.

Keeping the generic SDF-NCP layer independent lets the equations and tests
stabilize before choosing how to integrate with Chrono's production contact
container architecture.

### Later integration path

There are two plausible next integration paths:

- NSC-like nonlinear constraint path:
  - expose SDF gap and `J_phi(q)` as a normal unilateral constraint,
  - add smooth FB residual evaluation around the normal multiplier,
  - connect to descriptor-level unknowns and warm-started reactions.

- Field-contact NCP container:
  - reuse `FieldSampleQuery` and primitive/patch extraction,
  - build one or more SDF-NCP constraints from patch centers or sample points,
  - assemble the NCP residuals in a dedicated field-contact container.

The current generic point-mass layer is meant to be the smallest common kernel
for both paths.

### Validation case for this phase

The validation remains a single frictionless point mass against an SDF plane,
but now uses:

```text
gap = surface.Phi(q_next)
normal = normalize(surface.Grad(q_next))
J = normal^T
R_v = m(v_next - v) - dt(Q + normal * lambda)
R_lambda = Phi_eps(gap, lambda)
```

The concrete regression uses `ChPlaneSdfSurface([0, 1, 0], 0)` and 3D point
translation. Existing 2D `PointMassPlaneResidual` tests remain as backwards
compatibility checks.

## 10. Risks

- Smooth FB with `eps > 0` approximates exact complementarity; open contact and
  closed contact may have tiny positive products of order `eps^2`.
- Newton convergence depends on the initial multiplier guess and tolerance.
- The first residual is plane-specific and 2D point-mass specific.
- Direct insertion into Chrono's contact containers would require careful
  descriptor, variable, and constraint bookkeeping; this phase avoids that risk.
- Existing `field_contact` tests regenerate paper example figures; those
  generated artifacts should not be committed as part of SDF-NCP integration.
- A later full Chrono integration must decide whether SDF-NCP belongs in a new
  contact container, an NSC-like nonlinear constraint path, or a dedicated
  field-contact NCP assembly layer.

