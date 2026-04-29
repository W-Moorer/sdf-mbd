# SDF-NCP 与 field_contact 资产复用和对比计划

## 1. 目标

本阶段目标是基于仓库已有 `field_contact` 资产建立 SDF-NCP benchmark，而不是孤立 toy example。

当前实现遵守后端分离原则：

1. OpenVDB/OBJ/SparseSDF 前端加载和查询逻辑可以复用。
2. pressure-field 接触后端和 SDF-NCP 后端分离。
3. 共享前端输出 `ChSdfContactSampleQuery`、surface graph、triangle mesh 数据。
4. SDF-NCP 后端只消费 `phi(x)`、`grad phi(x)`、候选接触点和 Jacobian，不调用 pressure law、patch pressure integration 或 Chrono 主 contact descriptor。

## 2. 资产审计表

| asset | field_contact 使用位置 | 物理问题 | 几何类型 | SDF/field 表示 | SDF-NCP 当前状态 | 备注 |
| --- | --- | --- | --- | --- | --- | --- |
| `assets/headon_spheres` | `paper_example/paper_examples_openvdb.cpp::RunHeadonCase` | 两球迎面碰撞 | 球体 | 解析球 SDF 等价完整球体几何 | 已实现 | 3D 无摩擦自由动力学，单法向 NCP |
| `assets/headon_spheres_mass_ratio` | `paper_example/paper_examples_openvdb.cpp::RunHeadonCase` | 不同质量比两球碰撞 | 球体 | 解析球 SDF 等价完整球体几何 | 已实现 | 复用质量比和初始速度参数 |
| `assets/cam` | `assets/cam/cam_model.json`、`simple_cam.rmd`、`cam_data.csv` | cam-follower | `cam_body1.obj`、`cam_body2.obj` | OpenVDB full mesh SDF | 已实现 | cam speed-driven，follower y 自由动力学，active-band connected patch |
| `assets/eccentric_roller` | `paper_example/paper_examples_openvdb.cpp::RunCamCase` | 偏心 cam 与 roller follower | `eccentric_disk_cam.obj`、`roller_follower.obj` | OpenVDB full mesh SDF | 已实现 | roller surface active-band connected patch，follower y 自由动力学 |
| `assets/onset_stress` | `paper_example/paper_examples_openvdb.cpp::RunCamCase(check_onset=true)` | contact onset / stress onset | `onset_cam.obj`、`roller_follower.obj` | OpenVDB full mesh SDF | 已实现 | active-band connected patch，记录 onset time |
| `assets/simple_gear` | `demo_CH_simple_gear_openvdb.cpp`、`paper_example/paper_examples_openvdb.cpp::RunSimpleGearCase` | 双齿轮啮合 | `gear_21.obj`、`gear_22.obj`、RMD 参数 | OpenVDB full tooth mesh SDF | 已实现 | GEAR21 速度驱动，GEAR22 RX 自由动力学，齿面 active-band connected patch |
| `assets/rev_joint_clearance` | `demo_CH_rev_joint_clearance_openvdb.cpp` | 转动副间隙 | bore/pin OBJ mesh | OpenVDB/SparseSDF | 未实现 | 需要 3D 刚体姿态和关节约束组合，后续处理 |

## 3. 统一前端和后端边界

新增共享前端数据位于：

```text
src/chrono/collision/ChSdfContactGeometry.h
src/chrono/collision/ChSdfOpenVDB.h
```

共享前端负责：

```text
OBJ mesh -> triangle mesh
triangle mesh -> surface graph
triangle mesh -> OpenVDB level set
sample point -> ChSdfContactSampleQuery(phi, grad, world_pos, world_vel)
```

SDF-NCP 后端负责：

```text
g_i(q) = phi(x_i(q))
n_i = normalize(grad phi(x_i(q)))
J_i(q) = n_i^T dx_i(q)/dq
R_v = M(v_next - v) - dt(Q + sum_i J_i^T lambda_i)
R_lambda_i = Phi_eps(g_i, lambda_i)
```

pressure-field 后端仍保留自己的 pressure law、patch extraction 和 force integration，不与 SDF-NCP multiplier/residual 混合。

## 4. 当前 benchmark 入口

解析球体 benchmark：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe all
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres_mass_ratio
```

`headon_spheres` 和 `headon_spheres_mass_ratio` 已收录为 analytic-reference benchmark。当前 frictionless SDF-NCP 不包含 restitution，因此正式解析对比使用 no-restitution common normal velocity；`comparison.csv` 中保留 elastic reference 仅作为非目标对照。

OpenVDB full geometry benchmark：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe all
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe eccentric_roller
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe onset_stress
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe simple_gear
```

正式收录的 analytic-reference benchmark：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres_mass_ratio
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe eccentric_roller
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe onset_stress
```

其中 `onset_stress` 额外输出：

```text
results/sdf_ncp_benchmarks/onset_stress/analytic_comparison.csv
```

统一 Python runner：

```text
python scripts\run_sdf_ncp_field_contact_benchmarks.py
```

## 5. 各 case 建模方式

### cam

使用完整 OBJ：

```text
assets/cam/models/cam_body1.obj
assets/cam/models/cam_body2.obj
```

建模方式：

1. `cam_body1.obj` 生成 OpenVDB SDF。
2. `cam_body2.obj` 生成 follower surface graph。
3. 每步扫描完整 follower surface，找到 `min_phi`，再取 `phi <= min_phi + activation_band` 的连通样本作为 active patch。
4. follower 沿 y 方向自由动力学，cam 按资产中的转速驱动。
5. active patch 内全部样本进入 SDF-NCP residual，不做 top-N 压缩；如果当前窄激活带内只有一个样本，则该时间步自然退化为单样本 patch。

### eccentric_roller

使用完整 OBJ：

```text
assets/eccentric_roller/models/eccentric_disk_cam.obj
assets/eccentric_roller/models/roller_follower.obj
```

建模方式：

1. eccentric cam OBJ 生成 OpenVDB SDF。
2. roller follower OBJ 生成 surface graph。
3. 每步先扫描完整 roller surface，找到 `min_phi`，再取 `phi <= min_phi + activation_band` 的所有连通样本作为 active patch。
4. follower 沿 y 方向自由动力学，cam 按资产 motor speed 驱动。
5. 多点 SDF-NCP residual 使用 active patch 内全部样本的 `sum_i n_{i,y} lambda_i`。

### onset_stress

使用完整 OBJ：

```text
assets/onset_stress/models/onset_cam.obj
assets/onset_stress/models/roller_follower.obj
```

建模方式与 `eccentric_roller` 相同，但记录 min-gap 从正到非正的 onset time。JSON 中的 `phase` 只用于解析 envelope comparison；SDF 网格姿态按 field-contact paper example 从 `theta=0` 开始。

该算例属于 analytic-reference benchmark，不需要 RecurDyn/RMD 前端映射。一致性检查内容为：

1. `onset_stress_model.json` 中的解析参数、相位、转速和 `target_onset_time`。
2. 完整 OBJ/OpenVDB SDF 与解析 envelope 的 follower y 对比。
3. min-gap 过零得到的 observed onset time 与 `target_onset_time` 的误差。
4. `gap`、`lambda_n`、`ncp_residual` 和 `complementarity_error` 是否有界。

### simple_gear

使用完整 OBJ 和 RMD：

```text
assets/simple_gear/simple gear.rmd
assets/simple_gear/gear_21.obj
assets/simple_gear/gear_22.obj
assets/simple_gear/data/Gear22.csv
```

建模方式：

1. 复用 RMD marker、CM、part rotation 解析逻辑，将 OBJ 转换到 body-local 坐标。
2. GEAR21 full tooth mesh 生成 OpenVDB SDF。
3. GEAR22 full tooth mesh 生成 surface graph。
4. 每步扫描完整 GEAR22 齿面，找到 `min_phi`，再取 `phi <= min_phi + activation_band` 的连通齿面样本作为 active patch 查询 GEAR21 SDF。
5. GEAR21 RX 方向速度驱动，GEAR22 RX 为自由动力学未知量。
6. patch 内每个样本的广义 Jacobian 为：

```text
J_i = (r_i^world cross n_i)_x
```

7. 当前第一版使用单向 `GEAR22 surface -> GEAR21 SDF` 查询；后续可升级为双向一致化。

## 6. 输出规范

每个 case 输出到：

```text
results/sdf_ncp_benchmarks/<case_name>/
```

每个目录包含：

```text
trajectory.csv
summary.csv
comparison.csv
```

统一 trajectory 字段：

```text
time,body_id,px,py,pz,q0,q1,q2,q3,vx,vy,vz,wx,wy,wz,
contact_id,gap,lambda_n,penetration,ncp_residual,
complementarity_error,residual_norm,iterations,success
```

统一 summary 字段：

```text
case_name,method,asset,dt,t_end,num_steps,num_contacts,
max_penetration,mean_penetration,max_lambda_n,
max_ncp_residual,mean_ncp_residual,
max_complementarity_error,mean_complementarity_error,
mean_iterations,success_rate,runtime_seconds,notes
```

## 7. CTest 规范

当前注册的 benchmark tests：

```text
sdf_ncp_benchmark_headon_spheres
sdf_ncp_benchmark_headon_spheres_mass_ratio
sdf_ncp_benchmark_cam_full_geometry
sdf_ncp_benchmark_eccentric_roller_full_geometry
sdf_ncp_benchmark_onset_stress_full_geometry
sdf_ncp_benchmark_simple_gear_full_geometry
```

统一标签：

```text
collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

## 8. 当前限制

1. 仍然无摩擦。
2. `cam`、`eccentric_roller`、`onset_stress`、`simple_gear` 都使用 active-band connected patch；如果当前几何只有一个样本落入窄激活带，则该时间步退化为单样本 patch。
3. `simple_gear` 当前是单向 tooth mesh SDF query，不是 field-contact 的双向 pressure patch 等价模型。
4. 当前局部 Newton 使用有限差分 Jacobian。
5. 尚未接入 Chrono 主 contact container/descriptor。
6. `rev_joint_clearance` 尚未纳入 SDF-NCP benchmark。

## 9. 2026-04-28 修订：通用 NCP 后端与完整几何 benchmark 约定

本修订明确：SDF-NCP benchmark 不再把不同资产写成不同 NCP 后端。后端统一为 `ChSdfNcpGeneralizedProblem`，只接收 `gap_i`、`J_i`、`weight_i`、`lambda_i` 和质量/外力/时间步数据。OBJ、OpenVDB、SparseSDF、RecurDyn marker/RMD 解析、active patch 采样都属于前端映射层。

### 后端统一公式

```text
R_v[j] = M_j (v_next[j] - v[j]) - dt (Q_j + sum_i J_i[j] lambda_i weight_i)
R_lambda_i = Phi_eps(g_i, lambda_i)
Phi_eps(g, lambda) = sqrt(g^2 + lambda^2 + eps^2) - g - lambda
```

当前 OpenVDB full-geometry benchmark 对 active patch surface samples 使用统一面积归一化 quadrature weight：

```text
sum_i weight_i = 1
lambda_force_i = lambda_i * weight_i
```

这样 `cam`、`eccentric_roller`、`onset_stress` 和 `simple_gear` 不再混用不同的 patch 力尺度。

### 当前完整几何覆盖

已纳入 SDF-NCP benchmark 体系的资产：

```text
assets/cam
assets/eccentric_roller
assets/headon_spheres
assets/headon_spheres_mass_ratio
assets/onset_stress
assets/simple_gear
```

`simple_gear` 当前使用完整 `gear_21.obj` / `gear_22.obj` OpenVDB SDF，并采用双向 patch 候选：

```text
GEAR22 surface -> GEAR21 SDF
GEAR21 surface -> GEAR22 SDF
```

`rev_joint_clearance` 仍未实现 SDF-NCP benchmark，原因是它需要将间隙关节约束与 3D 刚体自由度一起映射到同一 generalized residual；这应作为下一批完整 RecurDyn 映射任务处理，而不是再写 reduced toy model。

### 参考类型与一致性要求

benchmark 必须先按参考来源分类，再决定一致性检查内容：

1. RecurDyn/RMD reference benchmark：资产本身包含 RMD、RecurDyn 输出 CSV 或等价 RecurDyn 前端定义。此类算例需要逐项对齐 RecurDyn 前端。
2. Analytic-reference benchmark：资产本身由解析几何、解析 envelope 或程序化参数定义，没有 RecurDyn 前端映射。此类算例应对齐解析解，不应强行补 RMD/marker/joint/motion。
3. Full-mesh numerical benchmark：资产有完整 OBJ/OpenVDB 几何和稳定 rollout，但 reference 可能是 field-contact 输出或暂未解析完的外部数值曲线。此类算例先用于验证 full-mesh SDF query、多点 patch 和通用 NCP 后端。

对于 RecurDyn/RMD reference benchmark，后续不能只看 SDF 几何是否完整，还必须逐项对齐：

1. RMD 中的 body CM、surface marker、part rotation、joint、motion。
2. 驱动自由度和自由动力学自由度。
3. 初始位置、初始速度、质量、惯量、重力和时间步。
4. 接触律，包括是否存在摩擦、阻尼、恢复系数或 pressure-field 分布。
5. 输出坐标的定义和 reference CSV 的列含义。

对于 analytic-reference benchmark，例如 `eccentric_roller`、`onset_stress`、`headon_spheres` 和 `headon_spheres_mass_ratio`，一致性要求是解析公式、解析 envelope、质量/速度参数和输出量定义一致，不需要也不应该补造 RecurDyn 前端映射。

当前结果只说明完整几何 SDF 查询、多点 patch 和通用 NCP 后端已跑通；只有被明确归类为 RecurDyn/RMD reference 的算例，才需要进一步声明是否与 RecurDyn 曲线一致。
## 2026-04-28 更新：descriptor 后端接入后的 benchmark 约定

SDF-NCP benchmark 后端现在分为两层：

1. `ChSdfNcpGeneralizedProblem`：用于局部 generalized residual、参数扫描和非 Chrono descriptor 的数值验证。
2. `ChSdfNcpConstraintContactSet`：用于将 SDF-NCP 法向接触作为 unilateral constraints 注入 Chrono descriptor，与 Chrono 刚体、关节和驱动一起求解。

`cam` benchmark 当前应优先使用第二种路径。命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam
```

该命令运行完整 cam OBJ/OpenVDB 几何、Chrono `ChSystemNSC`、cam rotation motor、follower translational joint 和 SDF-NCP descriptor contact constraints。旧 reduced/generalized cam 入口保留为：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_reduced
```

后续所有需要与 RecurDyn 对齐的 cam 结果，都应以 `cam` descriptor path 为主，而不是 `cam_reduced`。

### 对 field_contact 分支的后端分离说明

pressure-field / sdf-field contact 与 sdf-ncp 现在应保持后端分离：

- 前端可复用 OBJ、OpenVDB、SparseSDF、surface graph、RMD marker 解析和 sample query。
- pressure-field 分支继续进入既有 field_contact force/pressure 计算路径。
- sdf-ncp 分支进入 NCP residual 或 descriptor unilateral constraint 路径。
- 不应为不同 asset 编写不同 NCP 后端；case-specific 代码只应存在于前端映射层。

当前 `cam` 已验证 sdf-ncp 后端可以进入 Chrono descriptor，但仍未替代或修改 Chrono 主 contact container。

## 2026-04-28 更新：RecurDyn 前端映射与 SDF-NCP 后端边界

当前边界定义如下：

- RecurDyn/RMD 前端映射层负责解析 PART、MARKER、JOINT、MOTION、CSURFACE、SOLIDCONTACT、GRAVITY 和 reference CSV 输出定义。
- SDF/OpenVDB 前端负责从 OBJ/CSURFACE reference marker 构造 full mesh SDF query 和 active patch samples。
- SDF-NCP 后端只接收 body pair、point、normal、gap、weight，并注入 descriptor unilateral constraints。
- benchmark 对比层负责把 SDF-NCP trajectory 与 RecurDyn reference CSV 对齐，输出 `comparison.csv` 和 `rmd_mapping.csv`。

这保证 cam、gear、eccentric roller 等 asset 不需要各自实现不同 NCP 后端。case-specific 内容只允许存在于“资产实体名选择”和“reference 输出列定义”中；动力学约束、SDF query、patch assembly 和 NCP descriptor contact 必须走通用代码路径。

## 2026-04-29 更新：AABB BVH 粗检测规划

本轮加入的 AABB BVH 只服务于 SDF 前端候选筛选，不改变 SDF-NCP 后端 residual、Jacobian、FB complementarity 或 pressure/intensity 物理量。目标是减少每步对无关 surface samples 的 OpenVDB 查询次数，同时保持 active patch 的几何结果不变。

### 分层查询路径

```text
OBJ / mesh vertices
  -> ChSdfContactSurfaceGraph
  -> ChSdfContactSampleBvh(AABB over local surface samples)
  -> broad-phase candidate sample ids
  -> phi-only OpenVDB narrow query for active-band decision
  -> full OpenVDB phi + grad (+ Hessian where needed) only for active samples
  -> generic SDF-NCP backend
```

### 设计原则

1. BVH 放在 `ChSdfContactGeometry.h`，只依赖 sample position/area，不依赖 NCP 或 pressure-field。
2. `ChOpenVdbSdfGrid` 保存 mesh local AABB，并提供 `QueryLocalPhiOnly` / `SamplePhi` 给 active-band 粗筛使用。
3. `cam`、`eccentric_roller`、`onset_stress` 使用 follower sample BVH 与 cam SDF local bounds 做粗检测。
4. `simple_gear` 使用 GEAR21 和 GEAR22 各自的 sample BVH，双向查询分别筛选：

```text
GEAR22 surface BVH vs GEAR21 SDF bounds
GEAR21 surface BVH vs GEAR22 SDF bounds
```

5. 粗检测 AABB 使用 `max_activation_band + hysteresis + voxel padding` 扩张，避免漏掉 active band 内的 sample。
6. 若 BVH 候选无法满足最小 patch 样本数，代码回退到全 sample phi-only 扫描，优先保证精度不损失。
7. BVH 只减少候选搜索成本；真实接触力、力矩、NCP unknown、warm start 和 patch 权重仍由通用 SDF-NCP assembly 计算。

### 验收标准

1. `simple_gear_dt_001` 的 MAE/RMSE 不比引入 BVH 前变差。
2. `simple_gear_dt_0005` 和 `simple_gear_dt_0001` 不发散，`success_rate = 1`。
3. `max_penetration`、`ncp_residual`、`complementarity_error` 保持有界。
4. `ctest -L sdf_ncp_benchmark`、`ctest -L sdf_ncp`、`ctest -L field_contact` 和 `pytest` 保持通过。
