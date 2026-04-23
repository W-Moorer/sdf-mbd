# 里程碑 9 验收报告：Patch-level Tangential State / Stick-Slip 最小模型

## 1. 从 Tangential-Damping 到 Patch Stick-Slip 的升级

### 复用的内容
- [demo_CH_sdf_patch_constitutive_stickslip_openvdb.cpp](file:///e:/workspace/Multi-body%20Dynamics%20Solver/Multi-body%20Dynamics%20Solver/src/demos/core/demo_CH_sdf_patch_constitutive_stickslip_openvdb.cpp) 完整复用了 milestone 7/8 的以下框架：
  - Patch grouping（BFS 连通分量）
  - Patch descriptors（center, normal, area, mean_phi, effective_penetration）
  - Patch persistence / persistent ID 匹配（`PersistentPatchTrack`）
  - Patch normal constitutive law（EMA 平滑的 penetration_eff）
  - Case A / Case B 场景定义与五模式对照框架

### 新增的内容
- `StickSlipState` 枚举（Stick / Slip）
- `PersistentPatchTrack` 新增字段：
  - `tangential_displacement`（ChVector3d）：跨帧切向位移累积状态
  - `stick_slip_state`：当前 stick/slip 状态
  - `stick_steps` / `slip_steps`：各状态累计步数
  - `stick_to_slip_transitions`：状态转换计数
- `ComputePatchStickSlipForce()` 函数：完整实现 trial force + Coulomb clamp

### 进入 law 的 tangential state
- **tangential_displacement (xi_t)**：这是唯一跨帧保存的切向历史量
  - 更新规则：`xi_t^{k+1} = xi_t^k + vt * dt`
  - 参与力计算：`F_t_trial = -k_t * xi_t - c_t * vt`
  - Slip 时投影：`xi_t = F_t / (-k_t)`

## 2. 新增/修改的文件

| 文件路径 | 修改目的 | 修改摘要 |
|---------|---------|---------|
| `src/demos/core/demo_CH_sdf_patch_constitutive_stickslip_openvdb.cpp` | 新建里程碑 9 主 demo | 实现完整 stick-slip 最小模型 |
| `src/demos/core/CMakeLists.txt` | 构建配置 | 添加 `demo_CH_sdf_patch_constitutive_stickslip_openvdb` 到 `DEMOS_WITH_OPENVDB` |
| `out/milestone_09/sdf_patch_stickslip_case_A.csv` | Case A 五模式输出 | 包含 y_error, torque, force_std, stick/slip 统计 |
| `out/milestone_09/sdf_patch_stickslip_case_B.csv` | Case B 五模式输出 | 包含 y_error, torque, force_std, stick/slip 统计 |
| `out/milestone_09/sdf_patch_stickslip_compare.csv` | 跨场景对比 | Case A/B 各模式关键指标汇总 |
| `out/milestone_09/sdf_patch_stickslip_tracks.csv` | Patch track 统计 | 每个 patch 的 stick/slip 状态历史 |
| `reports/milestone_09/README_patch_stickslip.md` | README 文档 | 构建/运行方法、stick-slip law 定义、结果分析 |
| `reports/milestone_09/acceptance_report.md` | 验收报告 | 本文件 |

## 3. 新 Stick-Slip Law 定义与实现

### F_n（法向力）定义
```
penetration_eff_smooth = 0.3 * current_penetration + 0.7 * prev_track.effective_penetration
F_n = (k * penetration_eff_smooth + c * max(-vn_eff, 0)) * n
```
- k = 1e5 N/m（法向刚度）
- c = 500 Ns/m（法向阻尼）
- vn_eff = v_patch · n

### F_t_trial（试探切向力）定义
```
vt = v_patch - (v_patch · n) * n              （切向速度分量）
xi_t^{k+1} = xi_t^k + vt * dt                 （切向位移累积）
F_t_trial = -k_t * xi_t - c_t * vt
```
- k_t = 1e4 N/m（切向刚度，法向刚度的 1/10）
- c_t = 10 Ns/m（切向阻尼，法向阻尼的 1/50）

### Stick 条件
```
if |F_t_trial| <= mu * |F_n|:
    state = Stick
    F_t = F_t_trial
```
- mu = 0.3（Coulomb 摩擦系数）

### Slip 条件
```
if |F_t_trial| > mu * |F_n|:
    state = Slip
    F_t = -mu * |F_n| * normalize(vt)         （Coulomb 限幅）
    xi_t = F_t / (-k_t)                       （投影到等效状态）
```

### State 更新规则
- 每次接触步进：
  1. 计算 v_patch（patch center 速度）
  2. 分解 vn, vt
  3. 更新 xi_t
  4. 计算 F_t_trial
  5. 判断 stick/slip，计算 F_t
  6. 更新统计量（stick_steps, slip_steps, transitions）

### Torque 定义
```
T_patch = (patch_center - COM) × F_patch
F_patch = F_n + F_t
```

## 4. Case A 对照结果（旋转接触，ω_z = 5 rad/s）

### Milestone 8 Baseline vs 新模型

| 指标 | PatchConstit (MS6) | TangentialDamping (MS8) | StickSlip (MS9) |
|------|-------------------|------------------------|-----------------|
| final_y | 1.1496 m | 1.1165 m | 0.6459 m |
| y_error | 0.0471 m | 0.0802 m | **0.5508 m** ❌ |
| force_std_dev | 795.1 N | 790.6 N | 791.3 N |
| avg_torque_z | -2.91 Nm | -3.95 Nm | **+2.30 Nm** ❌ |
| avg_tangential_force | N/A | 19.59 N | 39.88 N |
| tangential_force_ratio | N/A | 4.33 | 0.094 |
| stick_steps | - | - | 360,564 |
| slip_steps | - | - | 563,976 |
| stick→slip transitions | - | - | 7,140 |
| avg_tangential_displacement | - | - | 0.008 m |
| stable | YES | NO | NO |

### Stick/Slip 行为分析
- Stick 比例：39%（360,564 / 924,540）
- Slip 比例：61%（563,976 / 924,540）
- 高频状态转换：平均每 130 步发生一次 stick→slip 转换
- 切向位移量级：0.008 m，对于接触尺度偏大

### 是否改善
**否，严重恶化。**
- y_error 是 tangential-damping 的 **6.9 倍**（0.080 → 0.551 m）
- torque_z 符号反转（-3.95 → +2.30 Nm），表明切向力方向异常
- 切向力幅值翻倍（19.6 → 39.9 N），但方向错误导致力学响应恶化
- 高频 stick-slip 转换引入数值振荡

## 5. Case B 对照结果（多 patch，盒体旋转）

### Baseline vs 新模型

| 指标 | PatchConstit (MS6) | TangentialDamping (MS8) | StickSlip (MS9) |
|------|-------------------|------------------------|-----------------|
| final_y | 1.1252 m | 1.1058 m | 1.1762 m |
| y_error | 0.0222 m | 0.0416 m | **0.0288 m** ✅ |
| force_std_dev | 486.6 N | 486.0 N | 482.7 N |
| avg_torque_z | +0.034 Nm | -0.153 Nm | **-0.024 Nm** ✅ |
| avg_tangential_force | N/A | 17.27 N | 16.89 N |
| tangential_force_ratio | N/A | 0.047 | 0.033 |
| stick_steps | - | - | 7,886,621 |
| slip_steps | - | - | 14,459,843 |
| stick→slip transitions | - | - | 44,768 |
| avg_tangential_displacement | - | - | 0.001 m |
| stable | YES | YES | YES |

### 稳定性
**保持稳定。**
- y_error 介于 PatchConstit 和 TangentialDamping 之间
- force_std_dev 最低（482.7 N）
- torque_z 幅值最小（0.024 Nm）
- 切向力合理（16.9 N）

### 多 Patch 行为
- 5 个 patch 分散切向力，避免单点集中
- 切向位移量级小（0.001 m），试探力在合理范围
- 65% 时间处于 slip 状态，但状态转换相对平缓

### Stick/Slip 统计
- Stick 比例：35%
- Slip 比例：65%
- 平均每 500 步发生一次转换（远低于 Case A 的 130 步）

### 是否保持
**是，Case B 保持稳定且指标合理。**

## 6. 适用边界总结

### 新模型在哪些场景下更好
1. **多 patch 场景**：多 patch 平均效应使切向力分散，避免单点集中放大模型缺陷
2. **低速接触**：切向速度小，xi_t 累积缓慢，试探力不会过大
3. **平坦接触面**：盒体几何使切向速度分布更均匀

### 哪些场景下仍有问题
1. **高速旋转接触**：Case A 中 ω_z = 5 rad/s 导致 xi_t 快速累积，频繁触发 slip
2. **单 patch 场景**：所有切向力集中在单一 patch 上，模型缺陷被放大
3. **高频 stick-slip 转换**：状态投影（xi_t = F_t / (-k_t)）过于简化，引入数值振荡

### 下一步最值得做什么
1. **引入速率依赖的切向状态管理**：xi_t 更新应考虑速率限幅或衰减机制
2. **改进 slip 状态的投影方法**：使用平滑过渡而非瞬时重置，避免切向力突变
3. **降低参数敏感性**：k_t 可能需要根据接触尺度自适应调整
4. **探索更复杂的摩擦理论**：如 rate-and-state friction 或 smooth Coulomb friction

## 7. 最终文件结构

```
src/demos/core/
└── demo_CH_sdf_patch_constitutive_stickslip_openvdb.cpp  (新建)

out/milestone_09/
├── sdf_patch_constitutive_tangential_case_A.csv           (已有，milestone 8)
├── sdf_patch_constitutive_tangential_case_B.csv           (已有，milestone 8)
├── sdf_patch_constitutive_tangential_compare.csv          (已有，milestone 8)
├── sdf_patch_stickslip_case_A.csv                         (新建，milestone 9)
├── sdf_patch_stickslip_case_B.csv                         (新建，milestone 9)
├── sdf_patch_stickslip_compare.csv                        (新建，milestone 9)
└── sdf_patch_stickslip_tracks.csv                         (新建，milestone 9)

reports/milestone_09/
├── README_patch_stickslip.md                              (新建)
└── acceptance_report.md                                   (新建，本文件)
```

## 8. 构建方法

```bash
cd "e:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build"
cmake --build . --config Release --target demo_CH_sdf_patch_constitutive_stickslip_openvdb
```

## 9. 运行方法

```bash
Set-Location "e:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_patch_constitutive_stickslip_openvdb.exe
```

## 10. 验收结论

### 是否通过里程碑 9
**工程层面：通过。**
- ✅ 工程成功构建
- ✅ 新 demo 成功运行
- ✅ 输出文件正确生成在 `out/milestone_09/`
- ✅ 新 stick-slip law 的数学定义清楚
- ✅ persistent state 中 tangential_displacement 真正进入 law
- ✅ Case A 有明确否定结论（y_error 0.080→0.551 m，恶化 6.9 倍）
- ✅ Case B 结果被检查且保持稳定（y_error 0.029 m，stable=YES）
- ✅ 文件归档正确，无散落临时文件

**科学层面：部分通过。**
- ❌ Case A 未改善，反而严重恶化（不符合 prompt 期望的"明确改善"）
- ✅ 但 prompt 明确要求"如果 stick-slip 没有改善，也必须如实报告"
- ✅ Case B 保持稳定，符合"不会破坏已有收益"的要求

### 是否适合进入更复杂 Patch Friction / Constitutive Model
**不适合直接升级，需要先解决以下问题：**

1. **切向位移累积机制不稳定**：高速接触下 xi_t 增长过快，需要速率限幅或衰减
2. **Slip 状态投影过于简化**：xi_t = F_t / (-k_t) 假设在连续 slip 物理中不成立
3. **参数敏感性高**：k_t=1e4 N/m 可能不适合所有接触尺度

**推荐路径：**
1. 先实现速率依赖的 xi_t 更新（如 `xi_t += clamp(vt * dt, -max_dx, max_dx)`）
2. 再实现 smooth stick-slip transition（如使用连续过渡函数替代硬阈值）
3. 最后再考虑引入更复杂的摩擦理论（如 rate-and-state）

**结论：最小 stick-slip 模型验证了 tangential state 的可行性，但在高速旋转场景下失效。需要更精细的状态管理才能进入生产使用。**
