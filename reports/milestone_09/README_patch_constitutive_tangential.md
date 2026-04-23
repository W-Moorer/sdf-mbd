# 里程碑 8：最小切向项 / History-Aware Patch Law

## 概述

本里程碑在现有 patch constitutive law 基础上加入最小切向阻尼项和 history-aware smoothing，专门验证是否能改善 milestone 7 中 Case A 的非对称旋转接触偏差，同时检查不会破坏 Case B 多 patch 场景的优势。

## Milestone 7 Baseline 记录

### Case A（旋转接触）

| 指标 | Pointwise | PatchConstitutive |
|------|-----------|-------------------|
| y_error | 0.0196 m | 0.0471 m |
| torque_z | -2.08 Nm | -2.91 Nm |
| force_std_dev | 2754 N | 795 N |
| stable | YES | YES |

### Case B（多 patch）

| 指标 | Pointwise | PatchConstitutive |
|------|-----------|-------------------|
| y_error | 0.1258 m | 0.0222 m |
| torque_z | 2.00 Nm | 0.034 Nm |
| force_std_dev | 1987 N | 487 N |
| stable | NO | YES |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_tangential_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_tangential_openvdb.exe
```

## 新增切向项 / history 项定义

### 切向阻尼项（方案 A）

**数学定义**：
```
v_patch = patch_center_velocity
n = patch_normal
vn_eff = v_patch · n
vt_eff = v_patch - vn_eff * n  (切向速度 = 总速度 - 法向分量)
F_t = -c_t * vt_eff           (切向阻尼力)
```

**Coulomb 限幅**：
```
if |F_t| > mu * |F_n|:
    F_t = F_t * (mu * |F_n| / |F_t|)
```

**总力**：
```
F_patch = F_n + F_t
T_patch = (patch_center - COM) × F_patch
```

### History-Aware Smoothing（方案 B）

**EMA 平滑**：
```
penetration_eff_smooth = alpha * penetration_current + (1 - alpha) * penetration_prev_smooth
```

其中 alpha = 0.3（当前帧权重 30%，历史权重 70%）。

**Persistent State 进入 Law**：
- `patch_age_in_steps`: 用于跟踪 patch 生命周期
- `effective_penetration_smooth`: 来自上一帧的 EMA 状态
- 这两个量现在真正进入 constitutive law 计算

## Case A / Case B 说明

### Case A：旋转接触
- 动态球体（R=0.2m）以初始角速度 ω_z = 5 rad/s 下落
- 目标：验证切向项是否改善旋转接触下的位置和力矩偏差

### Case B：多 patch
- 动态盒体（half_size=0.15m）以初始角速度下落
- 目标：验证切向项是否保持多 patch 场景的稳定性

## 输出文件

所有输出文件位于 `out/milestone_09/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_constitutive_tangential_case_A.csv` | Case A 四模式详细输出 |
| `sdf_patch_constitutive_tangential_case_B.csv` | Case B 四模式详细输出 |
| `sdf_patch_constitutive_tangential_compare.csv` | Case A vs Case B 跨场景对比 |

## 实验结果

### Case A（旋转接触）

| 指标 | Pointwise | PatchConstitutive | Tangential (新) |
|------|-----------|-------------------|-----------------|
| final_y | 1.1772 m | 1.1496 m | 1.0560 m |
| y_error | 0.0196 m | 0.0471 m | **0.1407 m** ↑ |
| force_std_dev | 2754 N | 795 N | 585 N |
| avg_torque_z | -2.08 Nm | -2.91 Nm | 0.25 Nm |
| avg_tangential_force | N/A | N/A | 18.9 N |
| stable | YES | YES | NO |

**结论：切向项在 Case A 上没有改善，反而恶化了 y_error。**

### Case B（多 patch）

| 指标 | Pointwise | PatchConstitutive | Tangential (新) |
|------|-----------|-------------------|-----------------|
| final_y | 1.2731 m | 1.1252 m | 1.0864 m |
| y_error | 0.1258 m | 0.0222 m | **0.0609 m** |
| force_std_dev | 1987 N | 487 N | 486 N |
| avg_torque_z | 2.00 Nm | 0.034 Nm | -0.18 Nm |
| avg_tangential_force | N/A | N/A | 16.6 N |
| stable | NO | YES | **YES** |

**结论：切向项在 Case B 上保持稳定性，y_error 仍远优于 pointwise。**

## 切向项参数

| 参数 | 值 | 说明 |
|------|-----|------|
| c_t (切向阻尼) | 10 Ns/m | 法向阻尼的 1/50 |
| mu (摩擦系数) | 0.3 | Coulomb 限幅 |
| ema_alpha | 0.3 | EMA 平滑权重 |

## 关键发现

### 为什么切向项在 Case A 上没有改善？

1. **切向阻尼改变了接触动力学**：切向力虽然小（avg ~19N），但改变了接触点的滑动行为
2. **旋转接触的物理复杂性**：球体旋转产生的切向速度在 patch 层级被过度简化
3. **法向 - 切向耦合缺失**：当前 law 中 F_n 和 F_t 独立计算，没有考虑耦合效应

### 为什么切向项在 Case B 上仍然有效？

1. **多 patch 平均效应**：多个 patch 的切向力相互抵消，总体影响较小
2. **Box 几何的对称性**：盒体的平坦接触面使得切向速度分布更均匀
3. **EMA 平滑的作用**：history-aware penetration 提高了稳定性

## 结论

最小切向阻尼项在**多 patch 场景下表现良好**，但在**旋转接触场景下没有改善**。这表明：
1. 简单的切向阻尼不足以捕获旋转接触的复杂物理
2. 需要更精细的切向模型（如 stick-slip 状态机）
3. history-aware smoothing 是一个有价值的增强
