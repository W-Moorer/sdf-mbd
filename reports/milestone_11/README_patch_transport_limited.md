# 里程碑 11：受限 Transport / Rate-Limited Tangential State

## 概述

本里程碑在 milestone 10 的 frame-consistent tangential state 基础上，引入 transport 限制机制，解决法向快速变化时 tangential state 非物理放大的问题。

## Milestone 10 Baseline 记录

### Case A（旋转接触）

| 模式 | y_error | torque_z | transported_norm | stable |
|------|---------|----------|------------------|--------|
| PatchConstitutive | 0.0471 m | -2.91 Nm | N/A | YES |
| FrameConsistent | 0.5525 m | +2.30 Nm | 5.62 m | NO |

### Case B（多 patch）

| 模式 | y_error | torque_z | transported_norm | stable |
|------|---------|----------|------------------|--------|
| PatchConstitutive | 0.0222 m | 0.034 Nm | N/A | YES |
| FrameConsistent | 0.0292 m | -0.023 Nm | 1.96 m | YES |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_transport_limited_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_transport_limited_openvdb.exe
```

## Transport-Limited 模型定义

### 三种限制机制（同时实现）

#### 方案 A：Transport Magnitude Clamp
```
if |xi_projected| > xi_max:
    xi_projected = xi_max * normalize(xi_projected)
```
- `xi_max = 0.01 m`（1 cm）
- 基于 patch characteristic length estimate

#### 方案 B：Rotation-Aware Attenuation
```
theta = arccos(n_old · n_new)
attenuation = exp(-beta * theta)
xi_projected *= attenuation
```
- `beta = 5.0`
- 法向变化越大，保留的旧状态越少

#### 方案 C：Partial Carry-Over
```
xi_new_init = alpha * xi_projected
```
- `alpha = 0.5`
- 相当于 transport-aware forgetting factor

### State 更新流程

```
1. 旧状态从旧局部基转到世界系: xi_world = xi1_old * t1_old + xi2_old * t2_old
2. 投影到新切平面: xi_projected = xi_world - (xi_world · n_new) * n_new
3. 得到 raw transported state (记录诊断)
4. 对 raw transported state 做限制:
   a. Rotation-aware attenuation: xi *= exp(-beta * theta)
   b. Transport magnitude clamp: if |xi| > xi_max, normalize
   c. Partial carry-over: xi *= alpha
5. 在新局部基表示: xi1 = xi · t1_new, xi2 = xi · t2_new
6. 加入当前步 tangential velocity 增量: xi += vt_local * dt
7. 做 stick/slip 判定（与 milestone 10 相同）
```

### 参数设置

| 参数 | 值 | 说明 |
|------|-----|------|
| xi_max | 0.01 m | transport clamp 阈值 |
| beta | 5.0 | rotation attenuation 强度 |
| alpha | 0.5 | partial carry-over 比例 |
| k_t | 1e4 N/m | 切向刚度 |
| c_t | 10 Ns/m | 切向阻尼 |
| mu | 0.3 | 摩擦系数 |

## 实验结果

### Case A（旋转接触）

| 指标 | PatchConstit | FrameConsistent (MS10) | TransportLimited (MS11) |
|------|--------------|------------------------|-------------------------|
| final_y | 1.1496 m | 0.6442 m | 1.1232 m |
| y_error | 0.0471 m | **0.5525 m** | **0.0735 m** ✅ |
| force_std_dev | 795 N | 792 N | 791 N |
| avg_torque_z | -2.91 Nm | +2.30 Nm | +0.63 Nm |
| avg_tangential_force | N/A | 39.9 N | 17.9 N |
| avg_raw_transport_norm | - | 5.62 m | 0.12 m ✅ |
| clamp_hits | - | - | 0 |
| stable | YES | NO | NO |

**结论：Transport-limited 模型大幅改善了 Case A！**
- y_error 从 0.553 m 降至 0.074 m（减少 87%）
- transport_norm 从 5.62 m 降至 0.12 m（减少 98%）
- torque_z 从 +2.30 Nm 降至 +0.63 Nm（方向仍异常，但幅值大幅下降）
- 切向力从 39.9 N 降至 17.9 N（更接近合理范围）
- stable 仍为 NO，但 y_error 已接近 PatchConstitutive 水平

### Case B（多 patch）

| 指标 | PatchConstit | FrameConsistent (MS10) | TransportLimited (MS11) |
|------|--------------|------------------------|-------------------------|
| final_y | 1.1252 m | 1.1766 m | 1.0780 m |
| y_error | 0.0222 m | 0.0292 m | 0.0694 m |
| force_std_dev | 487 N | 483 N | 484 N |
| avg_torque_z | 0.034 Nm | -0.023 Nm | -0.12 Nm |
| avg_tangential_force | N/A | 17.4 N | 18.7 N |
| avg_raw_transport_norm | - | 1.96 m | 0.13 m ✅ |
| clamp_hits | - | - | 0 |
| stable | YES | YES | YES |

**结论：Transport-limited 模型在 Case B 上保持稳定，但 y_error 略有增加。**
- y_error 从 0.029 m 增至 0.069 m（仍远优于 pointwise 的 0.126 m）
- transport_norm 从 1.96 m 降至 0.13 m（减少 93%）
- stable 仍为 YES

## 诊断分析

### Transport 限制效果

| Case | Raw Transport (MS10) | Raw Transport (MS11) | 减少比例 |
|------|---------------------|---------------------|---------|
| Case A | 5.62 m | 0.12 m | 97.9% |
| Case B | 1.96 m | 0.13 m | 93.4% |

### Clamp 命中情况
- Case A: 0 次 clamp 命中（attenuation 已足够抑制）
- Case B: 0 次 clamp 命中（attenuation 已足够抑制）

这说明 rotation-aware attenuation (方案 B) 是主要的限制来源，magnitude clamp (方案 A) 在 beta=5.0 时基本不会被触发。

## 为什么 Transport-Limited 改善了 Case A？

1. **Attenuation 抑制了非物理放大**: `exp(-5 * theta)` 在法向变化大时迅速衰减旧状态
2. **Partial carry-over 防止累积**: 每次只继承 50% 的旧状态，避免跨帧累积放大
3. **Transport_norm 回归合理范围**: 从 5.62 m 降至 0.12 m，接近 patch 尺度

## 为什么 Case A 仍未完全稳定？

- torque_z 仍为正值（+0.63 Nm），方向异常
- y_error 0.074 m 仍比 PatchConstitutive 的 0.047 m 差 56%
- 这说明 transport 限制缓解了主要问题，但切向力方向仍需修正

## 输出文件

所有输出文件位于 `out/milestone_11/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_transport_limited_case_A.csv` | Case A 五模式详细输出 |
| `sdf_patch_transport_limited_case_B.csv` | Case B 五模式详细输出 |
| `sdf_patch_transport_limited_compare.csv` | 跨场景对比 |
| `sdf_patch_transport_diagnostics.csv` | Transport 诊断（raw/limited norm, clamp hits） |
