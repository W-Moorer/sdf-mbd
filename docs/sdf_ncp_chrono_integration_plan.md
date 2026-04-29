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

## 11. 多点 SDF-NCP contact set assembly

### Why move from one constraint to a contact set

The single-point constraint proves the residual and local Newton solve for one
gap and one normal multiplier. A useful contact layer must also support several
simultaneous unilateral constraints before it can become a field-contact NCP
container. Even a point mass can hit a floor and a side wall at the same time,
and later patch/field contact will naturally produce multiple candidate contact
points.

The next abstraction is therefore a small contact set: a list of independent
frictionless SDF normal constraints that share the same point-mass velocity
unknowns and are solved in one local NCP system.

### Multi-contact residual

For the current point-mass-only version, all contact points use the same point
kinematics:

```text
x_i(q) = q
```

For contact `i`,

```text
g_i = surface_i.Phi(q_next)
n_i = normalize(surface_i.Grad(q_next))
```

The unknown vector is

```text
z = [vx_next, vy_next, vz_next, lambda_1, ..., lambda_N]
```

with

```text
q_next = q + dt * v_next
```

The assembled residual is

```text
R_v = m(v_next - v) - dt(Q + sum_i n_i * lambda_i)
R_lambda_i = Phi_eps(g_i, lambda_i)
```

For `N = 0`, the residual degenerates to the free point-mass implicit Euler
velocity equation:

```text
R_v = m(v_next - v) - dt Q
```

### Why this still stays outside Chrono's main solver

The contact set solve is still a local research kernel. It does not yet manage
Chrono descriptor offsets, active-set persistence, solver warm starting, or
interaction with `ChContactContainer`. Keeping it independent allows the NCP
assembly to be tested and compared against the Python prototype without
destabilizing existing SMC/NSC contact code.

### Role as a future field-contact NCP abstraction

The contact set is a natural precursor to a future field-contact NCP container:

- SDF/field sampling creates candidate points or patch representatives.
- Each candidate becomes one `ChSdfNcpContactPoint`.
- The contact set evaluates all gaps, normals, multipliers, and NCP residuals.
- A later container can map the same residual rows into Chrono's solver
  descriptor or a dedicated nonlinear complementarity solve.

### Validation case for this phase

The regression checks use two orthogonal planes:

```text
y = 0, normal [0, 1, 0]
x = 0, normal [1, 0, 0]
```

The tests verify direct contact-set evaluation, summed contact force, agreement
between the `N=1` contact set and the existing single-contact residual, and a
two-plane rollout with bounded penetration and nonnegative normal multipliers.

### Current limitations

- Friction is still absent.
- All contacts share point-mass translation kinematics, `x_i(q) = q`.
- The local Newton solve uses finite-difference Jacobians.
- The contact set is intended for small `N`; it is not optimized for large
  contact clouds.
- There is still no Chrono global descriptor/container integration.

## 12. 2D 刚体局部接触点 SDF-NCP residual

### Why rigid-body local points follow the point-mass contact set

The point-mass contact set verifies multi-contact complementarity assembly, but
its kinematics are still `x_i(q) = q`. The next minimum mechanics step is a
rigid body with local contact points. This verifies the contact Jacobian
definition used in the paper:

```text
J_phi(q) = grad phi(x_i(q))^T * dx_i(q)/dq
```

This phase keeps the solve local and small, while adding the rotational
contribution required before any later Chrono rigid-body descriptor mapping.

### 2D rigid-body state

The local prototype uses planar generalized coordinates and velocities:

```text
q = [x, y, theta]
v = [vx, vy, omega]
```

The mass matrix is

```text
M = diag(m, m, Iz)
```

and the external generalized force is

```text
Q = [0, -m*g0, 0]
```

### Local point kinematics

For local body point `r_i = [rx, ry]`, the world point is

```text
x_i(q) = p + R(theta) r_i
```

or explicitly,

```text
x_world = x + cos(theta) * rx - sin(theta) * ry
y_world = y + sin(theta) * rx + cos(theta) * ry
```

The point Jacobian with respect to `[x, y, theta]` is

```text
dx_i/dq =
    [1, 0, -sin(theta)*rx - cos(theta)*ry]
    [0, 1,  cos(theta)*rx - sin(theta)*ry]
```

### SDF gap, normal, and contact Jacobian

For contact point `i`,

```text
g_i(q) = surface_i.Phi(x_i(q))
n_i = normalize(surface_i.Grad(x_i(q)))
J_i = n_i^T dx_i/dq
```

The regression tests compare `J_i` against a finite-difference derivative of
`surface_i.Phi(x_i(q))`.

### 2D rigid-body residual

The unknown vector for `N` normal contacts is

```text
z = [vx_next, vy_next, omega_next, lambda_1, ..., lambda_N]
```

with implicit Euler update

```text
q_next = q + dt * v_next
```

The residual is

```text
R_v = M(v_next - v) - dt(Q + sum_i J_i^T * lambda_i)
R_lambda_i = Phi_eps(g_i, lambda_i)
```

This is the rigid-body analogue of the point-mass contact set, with contact
forces entering through generalized force rows `[Fx, Fy, torque_z]`.

### Current limitations

- The model is planar 2D.
- Contact is frictionless and normal-only.
- The nonlinear solve is still a local Newton iteration.
- The residual Jacobian used by Newton is finite-difference based.
- The implementation is independent from Chrono's global contact
  container/descriptor path.

## 13. 第一篇论文数值验证与复现实验计划

### Current SDF-NCP method layers

The repository now contains the independent method layers needed for the first
traditional SDF-NCP paper:

- Python prototype for SDF gap, normal, point-mass penalty contact, and
  smooth-Fischer-Burmeister SDF-NCP contact.
- C++ single-contact residual for a point mass on a plane.
- C++ multi-contact point-mass contact set assembly.
- C++ planar rigid-body local-contact residual that verifies
  `J_phi(q) = grad phi(x_i(q))^T * dx_i(q)/dq`.

These layers remain outside Chrono's production contact containers and global
solver descriptor. That separation keeps the first-paper experiments
reproducible and focused on the numerical formulation.

### Questions the first paper should answer

The numerical experiments should answer the following bounded questions without
overclaiming:

- Are the SDF gap and normal evaluations correct for simple analytical
  surfaces?
- Does the SDF contact Jacobian match finite differences of the SDF gap?
- Does the SDF-NCP formulation reduce penalty-stiffness sensitivity compared
  with a normal penalty baseline?
- How does the smooth Fischer-Burmeister parameter `eps` affect penetration,
  complementarity error, residuals, and solver iterations?
- Does a small planar rigid body with two local SDF contact points remain
  stable in a local rollout?

### CSV outputs and figures

The reproducibility suite should generate:

- `results/sdf_ncp/point_mass_plane/point_mass_plane.csv` and time-history
  figures for height, gap, penetration, contact force, and complementarity
  error.
- `results/sdf_ncp/penalty_sensitivity/summary.csv` and
  `max_penetration_vs_parameter.png`.
- `results/sdf_ncp/epsilon_sensitivity/summary.csv` plus plots for maximum
  penetration, complementarity error, and solver iterations versus `eps`.
- `results/sdf_ncp/geometry/sdf_contours_normals.png` for the method section.
- `results/sdf_ncp_cpp/rigidbody2d_rollout.csv` from the C++ local rigid-body
  residual.
- `results/sdf_ncp_cpp/figures/rigidbody2d_pose_vs_time.png`,
  `rigidbody2d_gaps_vs_time.png`, `rigidbody2d_lambdas_vs_time.png`, and
  `rigidbody2d_complementarity_vs_time.png`.
- `results/sdf_ncp_paper1/tables/method_summary.csv`, a compact table joining
  the penalty, epsilon, and C++ planar rigid-body statistics.

### Reproduction entry point

The intended paper reproduction command is:

```text
python scripts/run_sdf_ncp_paper1_experiments.py
```

It runs the Python point-mass experiments, the SDF geometry visualization, the
C++ planar rigid-body CSV export, the C++ rollout plotting step, and the summary
table generation. The script assumes `demo_CH_sdf_ncp_regression` has already
been built in `build/bin/Release`.

### Still out of scope

This paper-one reproducibility phase still does not include:

- 3D rigid-body contact.
- Friction cone complementarity.
- Flexible-body contact.
- Chrono descriptor or production contact-container integration.
- AI, neural networks, or machine-learning dependencies.

## 14. SDF 几何前端与接触后端分离

当前仓库中的 SDF 接触需要维护两条后端路线：

1. `field_contact` / pressure-field 后端。
2. `sdf_ncp` / nonlinear-complementarity 后端。

二者可以复用同一个几何前端，例如 mesh-to-sparse-SDF、surface samples、
surface graph、SDF query 和法向梯度；但后端的力学装配必须分开。

### 共享几何前端

共享层位于：

```text
src/chrono/collision/ChSdfContactGeometry.h
```

该层只包含后端无关数据：

```text
ChSdfContactSurfaceSample
ChSdfContactSurfaceGraph
ChSdfTriangleFace
ChSdfContactSampleQuery
MakeTriangleMeshSurfaceGraph
MakeSphereSurfaceGraph
```

它只表示：

```text
sample local position
sample area
sample adjacency
phi(x)
grad phi(x)
world sample position
world sample velocity
```

它不包含 pressure law、patch force、history transfer、Fischer-Burmeister
residual、lambda 或 Newton solve。

### pressure-field 后端

pressure-field 后端仍位于：

```text
src/chrono/collision/ChFieldContactPrimitives.h
src/chrono/collision/ChFieldContactRuntime.h
```

`ChFieldContactPrimitives.h` 现在通过别名消费共享几何类型：

```text
SurfaceGraph      = sdfcontact::ChSdfContactSurfaceGraph
TriangleFace      = sdfcontact::ChSdfTriangleFace
FieldSampleQuery  = sdfcontact::ChSdfContactSampleQuery
```

pressure-field 后端独有内容包括：

```text
PrimitivePatch
PatchExtractionSettings
NormalContactSettings
ComputeNormalContactIntegral
TangentialContactSettings
FieldContactPrimitiveTracker
history inheritance / split-merge classification
```

这些内容不应被 SDF-NCP 后端依赖。

### SDF-NCP 后端

SDF-NCP 后端仍位于：

```text
src/chrono/collision/ChSdfNcpContact.h
```

该后端独有内容包括：

```text
SmoothFischerBurmeister
ComplementarityError
ChSdfNcpContactConstraint
ChSdfNcpContactSet
point-mass residuals
2D rigid-body residuals
local Newton solves
```

SDF-NCP 可以直接消费共享几何 query：

```text
EvaluateSdfNcpContactQuery(ChSdfContactSampleQuery, lambda, eps)
```

这条入口只把共享几何中的 `phi` 和 `grad phi` 转换为 NCP gap、normal、
FB residual 和 complementarity diagnostics；它不调用 pressure-field
patch extraction 或 pressure force integration。

### 依赖方向

允许的依赖方向是：

```text
mesh / OpenVDB / sparse SDF frontend
        -> ChSdfContactGeometry
        -> ChFieldContactPrimitives / ChFieldContactRuntime

mesh / OpenVDB / sparse SDF frontend
        -> ChSdfContactGeometry
        -> ChSdfNcpContact
```

不允许的依赖方向是：

```text
ChSdfNcpContact -> ChFieldContactRuntime
ChSdfNcpContact -> pressure-field force law
ChFieldContactRuntime -> SDF-NCP lambda / FB residual
```

这样可以保证两个研究分支共享几何前端，但后端力学装配和验证路径保持独立。


## 15. 通用 SDF-NCP 后端与 RecurDyn/资产映射分离

本阶段将 SDF-NCP benchmark 的后端进一步收敛为统一的 generalized residual 层，避免为 `cam`、`eccentric_roller`、`onset_stress`、`simple_gear` 分别维护不同的 NCP 求解器。

### 通用后端输入

通用后端位于：

```text
src/chrono/collision/ChSdfNcpContact.h
```

新增核心结构：

```text
ChSdfNcpGeneralizedContact
ChSdfNcpGeneralizedProblem
ChSdfNcpGeneralizedResidual
ChSdfNcpGeneralizedDiagnostics
SolveSdfNcpGeneralizedProblem
```

它只消费：

```text
current_velocity
mass_diagonal
external_force
dt
eps
contact_count
evaluate_contacts(v_next) -> { gap_i, J_i, weight_i, lambda_i metadata }
```

其中 `evaluate_contacts` 是前端回调。回调可以来自解析 SDF、OBJ/OpenVDB SDF、RecurDyn marker 映射或未来 SparseSDF 查询，但通用 NCP 后端不读取 OBJ、RMD、OpenVDB 或 field-contact pressure patch 数据结构。

### 通用 residual

```text
z = [v_next, lambda_1, ..., lambda_n]

R_v[j] =
    M_j (v_next[j] - v[j])
    - dt (Q_j + sum_i J_i[j] * lambda_i * weight_i)

R_lambda_i = Phi_eps(g_i, lambda_i)

Phi_eps(g, lambda) = sqrt(g^2 + lambda^2 + eps^2) - g - lambda
```

`weight_i` 的语义统一为 active contact patch 的 quadrature weight。当前 OpenVDB benchmark 对进入 patch 的 surface samples 使用面积归一化权重，使 `sum_i weight_i = 1`。CSV 中：

```text
lambda_n       = sample normal NCP multiplier
contact_weight = active-patch quadrature weight
lambda_force   = lambda_n * contact_weight
```

这避免 cam 使用一种力尺度、gear 使用另一种力尺度。后续如果需要真实面积积分压力，应作为统一的后端建模选项引入，而不是按 case 特判。

### 前端映射职责

RecurDyn/资产映射层仍然是 case-specific，但它不是 NCP 后端。它负责：

```text
OBJ / OpenVDB / SparseSDF loading
RMD marker, CM, part rotation parsing
driven coordinate and free coordinate selection
local/world coordinate transform
active patch sample selection
J_i = grad phi(x_i(q))^T dx_i(q)/dq
```

例如：

- cam-like case 映射为 1 个 follower y-velocity DOF。
- simple_gear 映射为 GEAR22 RX velocity DOF，当前使用 GEAR21 驱动、GEAR22 自由。
- gear contact patch 已做双向候选：`GEAR22 -> GEAR21 SDF` 和 `GEAR21 -> GEAR22 SDF`。

### 与 RecurDyn 曲线的一致性状态

当前实现使用完整 OBJ/OpenVDB 几何和自由动力学局部 NCP residual，但还不能声称与 RecurDyn 参考曲线完全对应。主要差异包括：

1. 当前 SDF-NCP 是无摩擦法向 NCP；RecurDyn/field-contact reference 可能包含不同的接触律、阻尼、摩擦、恢复系数或 pressure-field 分布。
2. 当前 benchmark 只映射了必要的驱动/自由自由度，没有接入 RecurDyn/Chrono 的完整多体约束系统。
3. 当前 active patch 与 RecurDyn 的接触候选、接触区域和接触力分布并非逐点等价。
4. 当前局部 Newton 使用有限差分 Jacobian，尚未使用完整解析 Jacobian 或全局 descriptor。

因此，本阶段结果应解释为：通用 SDF-NCP 后端、完整几何 SDF 查询、多点 patch 候选和局部自由动力学已跑通；要达到“与 RecurDyn 参考建模一模一样”，下一步必须逐项对齐 RMD 约束、驱动、初始状态、接触律和输出坐标定义。
## 2026-04-28 更新：SDF-NCP 进入 Chrono 多体约束求解路径

针对 cam benchmark，当前实现已经从局部 reduced ODE / generalized residual 推进到 Chrono 多体约束描述符路径。新增的独立后端为：

```text
src/chrono/physics/ChSdfNcpConstraintContact.h
```

该后端以 `ChPhysicsItem` 形式存在，不修改 `ChSystem`、`ChContactContainerNSC`、`ChContactContainerSMC` 或 Chrono 主 contact container。它通过固定容量的 `ChConstraintTwoBodies` unilateral constraints 将 SDF-NCP 法向接触样本注入 Chrono descriptor，由 Chrono 的多体约束求解器与关节、驱动和刚体变量一起求解。

### 接入方式

SDF / OpenVDB / OBJ 前端仍负责生成接触样本：

```text
body_a, body_b, point_abs, normal_abs, gap, weight, contact_id
```

descriptor 后端负责构造无摩擦法向约束 Jacobian：

```text
g_dot = n^T(v_A + omega_A x r_A - v_B - omega_B x r_B)

Cq_A = [ n^T,  R_A^T (r_A x n)^T ]
Cq_B = [-n^T, -R_B^T (r_B x n)^T ]
```

其中角速度行按 Chrono body-local angular velocity 变量写入 `ChConstraintTwoBodies`。求解后的约束乘子通过 `ConstraintsFetch_react` 取回，并用于计算：

```text
lambda_force = lambda_n * contact_weight
penetration = max(0, -gap)
Phi_eps(g, lambda_force)
complementarity_error(g, lambda_force)
```

### cam 当前建模状态

`demo_CH_sdf_ncp_benchmarks_openvdb.exe cam` 当前使用：

```text
ChSystemNSC
ChBodyAuxRef cam body
ChBodyAuxRef follower body
ChLinkMotorRotationSpeed for RevJoint1.RMotion
ChLinkMateGeneric for TraJoint1 follower translation
ChSdfNcpConstraintContactSet for SDF-NCP normal contact
full cam OBJ/OpenVDB SDF geometry
```

也就是说，cam 已进入 Chrono 多体约束求解路径；旧的局部 generalized/reduced 入口保留为 `cam_reduced`，仅用于对照和回归。

### 仍未完成的对齐项

当前仍不能声称与 RecurDyn 曲线完全对应，原因不是 SDF-NCP 仍停留在局部 reduced 求解，而是以下建模项尚未逐项完全等价：

1. `RevJoint1` 与 `TraJoint1` 的 marker frame / REULER 方向需要继续按 RMD 精确映射到 Chrono joint frames。
2. RecurDyn `SolidContact1` 中的 penalty/damping/contact-law 参数与当前 frictionless hard NCP 模型不同。
3. 输出坐标需要确认是否与 reference CSV 的 marker/body 坐标定义完全一致。
4. 当前 descriptor 接入仍是独立 SDF-NCP physics item，不是 Chrono 主 `ChContactContainerNSC/SMC` 路径。

因此本阶段结论应写为：SDF-NCP 已具备进入 Chrono 多体约束求解路径的独立后端；下一步重点是逐项对齐 RecurDyn 前端建模和接触律，而不是继续扩展 case-specific reduced residual。

## 2026-04-29 更新：SDF-NCP 后端数值规划，减少有限差分依赖

当前 SDF-NCP 后端已经覆盖 full OBJ/OpenVDB 几何、多点 patch、Chrono 多体约束路径、simple gear 的双向 patch 和基准化输出。下一阶段的核心不应继续堆叠 case-specific 修正，而应把后端数值结构向成熟 NSC/NCP 路径收敛。

### 优先级 1：保留通用后端，不做算例特化

后端统一消费以下数据：

```text
gap_i
normal_i
J_i
weight_i
sample_id / patch_id
body and marker mapping
```

所有 cam、eccentric_roller、simple_gear、onset_stress 和 sphere benchmark 都应走同一套 residual、scaling、warm start、patch history 和 diagnostics。算例差异只能出现在前端建模映射层，例如 OBJ/OpenVDB 查询、RMD marker/joint/driver 解析、候选 sample 生成和输出坐标定义。

### 优先级 2：从 force-level lambda 过渡到 impulse / velocity-level NCP

成熟 NSC 方法通常在速度层处理不可穿透：

```text
0 <= w_i + beta * g_i / dt + cfm * p_i  ⟂  p_i >= 0
```

其中 `p_i` 是接触冲量或局部接触强度，`w_i = J_i v_next` 是法向相对速度。该形式比直接在位置 gap 与 force lambda 上做 smooth-FB 更接近传统 NSC，也能减少小时间步下 `dt` 缩放导致的病态。

本轮已新增一个最小 AD 原型：

```text
ComputeSdfNcpGeneralizedImpulseMixedAd
ChSdfNcpImpulseMixedSettings
ChSdfNcpAdScalar
```

该原型暂时不替换现有 benchmark 默认求解器，只验证 impulse/mixed complementarity residual 与 Jacobian 可以通过前向自动微分直接得到。

### 优先级 3：优先使用自动微分或解析导数，不再把有限差分作为主路径

有限差分适合早期验证，但不适合作为 full geometry benchmark 的长期 Newton Jacobian：

1. active patch 数量变化时，有限差分会把几何切换误差混入导数；
2. OpenVDB 查询噪声会被差分步长放大；
3. 多齿、多 patch 接触时，差分成本随自由度和接触数快速增长；
4. 小时间步下 residual scaling 更敏感，差分 Jacobian 更容易失真。

规划路径如下：

```text
阶段 A：后端代数 AD
    对 mass block、J^T p、scaled Fischer-Burmeister、Baumgarte/CFM、pressure scaling 做 AD。

阶段 B：运动学解析导数
    对 rigid body marker/local point 到 world point 的映射提供解析 dx/dq 和 dJ/dq。

阶段 C：OpenVDB/SparseSDF 几何导数
    从插值 SDF 中提供 phi、grad phi、Hessian 或至少一致的 grad/J derivative。

阶段 D：半光滑 active-set Jacobian
    对 patch 进入/退出使用 hysteresis 和固定 active set Newton；active set 变化只在外层更新。
```

这意味着未来 benchmark 默认路径应尽量使用 AD/解析 Jacobian；有限差分只保留为 debug checker 或 fallback。

### 优先级 4：scaled FB 与自动尺度

simple gear 当前误差受接触尺度、patch 分区和 NCP 条件数影响。后端应统一提供自动尺度：

```text
velocity_scale ~= characteristic_gap / dt
impulse_scale  ~= mass_effective * velocity_scale / contact_weight_scale
gap_scale      ~= SDF voxel size or active-band scale
```

smooth-FB 应作用在无量纲变量上，避免不同 benchmark 的 lambda、pressure、impulse 物理量混用。

### 优先级 5：与 RecurDyn 对齐的边界

如果曲线仍有差距，首先检查前端等价性：

1. RMD joint、marker frame、driver 方向和自由度；
2. 接触律是否为 hard NCP、penalty、damping 或 GGEOMCONTACT；
3. 输出坐标是否是同一 marker/body/方向；
4. OBJ/OpenVDB 离散化精度、active band 和 patch 采样；
5. 时间积分器、步长和求解器容差。

只有在这些项逐项对齐后，剩余误差才应归因于 SDF-NCP 后端理论或数值算法。

## 2026-04-29 更新：simple_gear 时间步稳定性阶段 1-7 落地

本轮围绕 simple_gear 的小时间步发散问题，按后端通用化原则完成以下七个阶段：

### 阶段 1：基线冻结

保留 `dt = 0.001 s` 下已有 hard SDF-NCP force-level 后端作为当前基线，不让新算法在默认步长下变差。当前基线为：

```text
case: simple_gear_dt_001
dt: 0.001
success_rate: 1
max_penetration: 8.2869746620417573e-09
MAE: 0.065929736198861666
RMSE: 0.10246898597734395
final omega22: -1.0511661888263963
```

### 阶段 2：通用 impulse / velocity-level NCP 后端

新增通用后端：

```text
SolveSdfNcpGeneralizedImpulseMixedProblem
ComputeSdfNcpGeneralizedImpulseMixedAd
ChSdfNcpImpulseMixedSettings
```

该后端使用速度层 mixed complementarity：

```text
R_v = M(v_next - v_current) - dt Q - sum_i J_i^T p_i weight_i

w_i = J_i v_next + b_i + beta g_i / dt + cfm p_i

0 <= w_i ⟂ p_i >= 0
```

其中 `b_i` 是 prescribed/driver 自由度造成的 affine normal velocity offset，例如 simple_gear 中主动轮 GEAR21 的驱动角速度。

### 阶段 3：自动微分 Jacobian

mixed 后端的 residual/Jacobian 使用轻量前向自动微分 `ChSdfNcpAdScalar` 计算，不再对该代数块使用黑箱有限差分。OpenVDB/SDF 查询仍作为前端 oracle 提供 `gap`、`J` 和 active samples；OpenVDB 几何 Hessian / active-set 半光滑导数仍是后续工作。

### 阶段 4：driver gap-rate offset

simple_gear 前端为每个 contact sample 提供：

```text
normal_velocity_offset = -J_driver omega21
```

后端只消费该 offset，不知道具体是齿轮、cam 还是其他机构。这避免在 velocity-level NCP 中漏掉 prescribed motion。

### 阶段 5：dt-scaled impulse weight

旧 force-level 后端使用 `force_weight`。mixed 后端使用冲量未知量，因此统一转换为：

```text
impulse_weight = dt * force_weight
```

由于原自动 force scale 随 `1/dt` 变化，该转换使小时间步下的 impulse scale 近似保持有界，避免 `dt -> 0` 时接触块病态放大。

### 阶段 6：dt-adaptive backend selection

默认规则：

```text
dt >= 0.001  : 保留当前 force-level hard SDF-NCP 基线
dt <  0.001  : 启用 impulse / velocity-level mixed SDF-NCP + AD Jacobian
```

也可通过环境变量覆盖：

```text
SDF_NCP_SIMPLE_GEAR_IMPULSE_MIXED
SDF_NCP_SIMPLE_GEAR_IMPULSE_BETA
SDF_NCP_SIMPLE_GEAR_IMPULSE_CFM
SDF_NCP_SIMPLE_GEAR_VELOCITY_SCALE
SDF_NCP_SIMPLE_GEAR_IMPULSE_SCALE
```

当前 mixed 默认参数：

```text
beta = 0.3
cfm = 1.0e-3
velocity_scale = 1.0
impulse_scale = 1.0
```

### 阶段 7：验收脚本和结果

新增验收脚本：

```text
scripts/validate_sdf_ncp_simple_gear_dt_sweep.py
```

输出：

```text
results/sdf_ncp_benchmarks/simple_gear_dt_sweep_summary.csv
```

当前验收结果：

| case | dt | success_rate | max_penetration | MAE | RMSE | final omega22 |
|---|---:|---:|---:|---:|---:|---:|
| simple_gear_dt_001 | 0.001 | 1 | 8.2869746620417573e-09 | 0.065929736198861666 | 0.10246898597734395 | -1.0511661888263963 |
| simple_gear_dt_0005 | 0.0005 | 1 | 1.0006613138102693e-07 | 0.063336256860321918 | 0.13101143283368741 | -0.96575938464339262 |
| simple_gear_dt_0001 | 0.0001 | 1 | 1.0003219585996703e-07 | 0.055549608448871655 | 0.11795684023219731 | -0.96901542949800457 |

结论：`dt=0.001` 未比当前基线变差；`dt=0.0005` 和 `dt=0.0001` 不再发散；三组 `success_rate = 1`；穿透量有界；MAE/RMSE 没有随时间步减小而爆炸。
