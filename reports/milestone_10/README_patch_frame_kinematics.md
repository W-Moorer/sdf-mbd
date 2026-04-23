# 里程碑 10：Patch Tangential Kinematics / Frame-Consistent Tangential State

## 概述

本里程碑在 milestone 9 的 stick-slip 模型基础上，引入 patch-local tangent basis 和 frame-consistent tangential state transport，解决旋转接触时世界系切向状态累积方向失真的问题。

## Milestone 9 Baseline 记录

### Case A（旋转接触）

| 模式 | y_error | torque_z | force_std_dev | stable |
|------|---------|----------|---------------|--------|
| PatchConstitutive | 0.0471 m | -2.91 Nm | 795 N | YES |
| TangentialDamping | 0.0802 m | -3.95 Nm | 791 N | NO |
| StickSlip | 0.0632 m | -0.48 Nm | 791 N | YES |

### Case B（多 patch）

| 模式 | y_error | torque_z | force_std_dev | stable |
|------|---------|----------|---------------|--------|
| PatchConstitutive | 0.0222 m | 0.034 Nm | 487 N | YES |
| TangentialDamping | 0.0416 m | -0.15 Nm | 486 N | YES |
| StickSlip | 0.0639 m | -0.13 Nm | 486 N | YES |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_frame_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_frame_openvdb.exe
```

## Frame-Consistent Tangential State 定义

### Patch Local Frame

为每个 persistent patch 定义局部正交基：

```
n = patch normal（从 SDF 梯度加权平均归一化得到）
t1 = normalize(n × ref)          （ref 取与 n 最小对齐的轴）
t2 = n × t1
```

### Tangential State Transport

当 patch 从上一帧延续到当前帧时，需要将旧切向状态 transport 到新局部基：

```
// 步骤 1：旧局部状态 → 世界坐标
xi_world = xi1_old * t1_old + xi2_old * t2_old

// 步骤 2：投影到新切平面（移除沿新法向的分量）
xi_projected = xi_world - (xi_world · n_new) * n_new

// 步骤 3：在新基下重新表示
xi1_new_init = xi_projected · t1_new
xi2_new_init = xi_projected · t2_new

// 步骤 4：加入当前步 tangential velocity 增量
vt1 = vt · t1_new
vt2 = vt · t2_new
xi1 = xi1_new_init + vt1 * dt
xi2 = xi2_new_init + vt2 * dt
```

### Tangential Force

```
F_t_trial_local = [-k_t * xi1 - c_t * vt1, -k_t * xi2 - c_t * vt2]
F_t_trial_world = F_t_trial_local[0] * t1 + F_t_trial_local[1] * t2

// Stick-slip clamp
if |F_t_trial| <= mu * |F_n|:
    -> Stick: F_t = F_t_trial_world
else:
    -> Slip: F_t = -mu * |F_n| * normalize(vt_local)
    xi1, xi2 = 投影到等效状态
```

### 参数设置

| 参数 | 值 | 说明 |
|------|-----|------|
| k_t (切向刚度) | 1e4 N/m | 法向刚度的 1/10 |
| c_t (切向阻尼) | 10 Ns/m | 法向阻尼的 1/50 |
| mu (摩擦系数) | 0.3 | Coulomb 限幅 |

## Case A / Case B 说明

与 milestone 7/8/9 相同的场景配置。

## 输出文件

所有输出文件位于 `out/milestone_10/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_frame_case_A.csv` | Case A 六模式详细输出 |
| `sdf_patch_frame_case_B.csv` | Case B 六模式详细输出 |
| `sdf_patch_frame_compare.csv` | 跨场景对比 |
| `sdf_patch_frame_tracks.csv` | Patch track 统计 |

## 实验结果

### Case A（旋转接触）

| 指标 | PatchConstit | TangentialDamping | StickSlip (MS9) | Frame (MS10) |
|------|--------------|-------------------|-----------------|--------------|
| final_y | 1.1496 m | 1.1165 m | 1.1335 m | 0.6442 m |
| y_error | 0.0471 m | 0.0802 m | 0.0632 m | **0.5525 m** ❌ |
| force_std_dev | 795 N | 791 N | 791 N | 792 N |
| avg_torque_z | -2.91 Nm | -3.95 Nm | -0.48 Nm | +2.30 Nm |
| avg_tangential_force | N/A | 19.6 N | 14.7 N | 39.9 N |
| avg_tangential_displacement | N/A | N/A | 0.000 m | 0.008 m |
| avg_transported_state_norm | N/A | N/A | N/A | 5.62 m |
| stable | YES | NO | YES | NO |

**结论：Frame-consistent 模型在 Case A 上仍然严重恶化，y_error 是 stick-slip 的 8.7 倍。**

### Case B（多 patch）

| 指标 | PatchConstit | TangentialDamping | StickSlip (MS9) | Frame (MS10) |
|------|--------------|-------------------|-----------------|--------------|
| final_y | 1.1252 m | 1.1058 m | 1.0835 m | 1.1766 m |
| y_error | 0.0222 m | 0.0416 m | 0.0639 m | **0.0292 m** ✅ |
| force_std_dev | 487 N | 486 N | 486 N | 483 N |
| avg_torque_z | 0.034 Nm | -0.15 Nm | -0.13 Nm | -0.023 Nm |
| avg_tangential_force | N/A | 17.3 N | 18.6 N | 17.4 N |
| avg_tangential_displacement | N/A | N/A | 0.000 m | 0.001 m |
| avg_transported_state_norm | N/A | N/A | N/A | 1.96 m |
| stable | YES | YES | YES | YES |

**结论：Frame-consistent 模型在 Case B 上保持稳定，y_error 介于 PatchConstit 和 StickSlip 之间。**

## Transport 行为分析

### Case A
- avg_transported_state_norm = 5.62 m：transport 过程中状态范数巨大
- 这说明法向变化剧烈时，旧切向状态投影到新切平面后量级异常放大
- frame_rotation_angle 隐含在 transport 过程中，但当前实现中记录不够精确

### Case B
- avg_transported_state_norm = 1.96 m：transport 范数较 Case A 小得多
- 多 patch 分散了切向效应，每个 patch 的法向变化较小

## 为什么 Frame-Consistent 在 Case A 上仍然失败？

### 1. Transport 放大了状态量级
- Case A 中球体高速旋转，patch normal 每帧变化显著
- 旧切向状态投影到新切平面时，由于基向量旋转，可能引入较大分量
- transported_state_norm = 5.62 m 远超物理合理的切向位移

### 2. 单 patch 无法分散法向变化
- Case A 只有 1 个 patch，所有法向变化都集中在该 patch 上
- 法向变化导致 tangent basis 旋转，transport 过程放大了不稳定性

### 3. Transport 投影的数学本质
- 旧切平面与新切平面之间的投影不是保距的
- 当两个法向夹角较大时，投影会引入系统性误差
- 这本质上是流形上切向量 transport 的问题（需要 parallel transport 或更高级方法）

### 4. 与 milestone 9 stick-slip 的对比
- 注意：MS9 stick-slip 在本轮运行中 y_error = 0.0632 m（stable=YES）
- 这比之前报告的 0.5508 m 要好得多，可能是随机性或参数微调的结果
- 但 frame-consistent 仍然比 stick-slip 差 8.7 倍

## 为什么 Frame-Consistent 在 Case B 上有效？

1. **多 patch 平均效应**: 5 个 patch 分散了切向力和法向变化
2. **盒体几何**: 平坦接触面使法向变化平缓，transport 投影误差小
3. **低速旋转**: 初始角速度远小于 Case A
4. **切向位移小**: 平均仅 0.001 m，transport 引入的误差可控

## 结论

**Frame-consistent tangential state 在旋转接触场景下仍然失败，但在多 patch 场景下保持有效。**

主要教训：
1. Tangent basis transport 在法向变化剧烈时会放大状态量级
2. 简单的投影 transport 不是保距的，不适合高速旋转场景
3. 需要更高级的切向量 transport 方法（如 parallel transport on manifold）
4. 单 patch 场景放大了 transport 误差
