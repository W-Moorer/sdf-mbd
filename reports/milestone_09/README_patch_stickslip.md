# 里程碑 9：Patch-level Tangential State / Stick-Slip 最小模型

## 概述

本里程碑在 milestone 8 的切向阻尼项基础上，引入 patch-level tangential displacement state 和最小 stick-slip 规则，验证是否能改善 Case A 旋转接触问题，同时保持 Case B 多 patch 场景的稳定性。

## Milestone 8 Baseline 记录

### Case A（旋转接触）

| 模式 | y_error | torque_z | force_std_dev | stable |
|------|---------|----------|---------------|--------|
| PatchConstitutive | 0.0471 m | -2.91 Nm | 795 N | YES |
| TangentialDamping | 0.0802 m | 0.25 Nm | 791 N | NO |

### Case B（多 patch）

| 模式 | y_error | torque_z | force_std_dev | stable |
|------|---------|----------|---------------|--------|
| PatchConstitutive | 0.0222 m | 0.034 Nm | 487 N | YES |
| TangentialDamping | 0.0416 m | -0.15 Nm | 486 N | YES |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_stickslip_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_stickslip_openvdb.exe
```

## Stick-Slip Law 定义

### 法向部分（保持不变）
```
penetration_eff_smooth = 0.3 * current + 0.7 * previous (EMA)
F_n = (k * penetration_eff_smooth + c * max(-vn_eff, 0)) * n
```

### 切向状态部分（新增）
```
xi_t^{k+1} = xi_t^k + vt * dt           (切向位移累积)
F_t_trial = -k_t * xi_t - c_t * vt      (试探切向力)
```

### Stick-Slip 规则
```
if |F_t_trial| <= mu * |F_n|:
    -> Stick: F_t = F_t_trial
    xi_t 保持累积
else:
    -> Slip: F_t = -mu * |F_n| * normalize(vt)
    xi_t = F_t / (-k_t)                 (投影到等效状态)
```

### 总力
```
F_patch = F_n + F_t
T_patch = (patch_center - COM) × F_patch
```

### 参数设置
| 参数 | 值 | 说明 |
|------|-----|------|
| k_t (切向刚度) | 1e4 N/m | 法向刚度的 1/10 |
| c_t (切向阻尼) | 10 Ns/m | 法向阻尼的 1/50 |
| mu (摩擦系数) | 0.3 | Coulomb 限幅 |

## Case A / Case B 说明

与 milestone 7/8 相同的场景配置。

## 输出文件

所有输出文件位于 `out/milestone_09/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_stickslip_case_A.csv` | Case A 五模式详细输出 |
| `sdf_patch_stickslip_case_B.csv` | Case B 五模式详细输出 |
| `sdf_patch_stickslip_compare.csv` | 跨场景对比 |
| `sdf_patch_stickslip_tracks.csv` | Patch track 统计 |

## 实验结果

### Case A（旋转接触）

| 指标 | PatchConstit | TangentialDamping | StickSlip (新) |
|------|--------------|-------------------|----------------|
| final_y | 1.1496 m | 1.1165 m | 0.6459 m |
| y_error | 0.0471 m | 0.0802 m | **0.5508 m** ❌ |
| force_std_dev | 795 N | 791 N | 791 N |
| avg_torque_z | -2.91 Nm | -3.95 Nm | 2.30 Nm |
| avg_tangential_force | N/A | 19.6 N | 39.9 N |
| stick_steps | - | - | 360,564 |
| slip_steps | - | - | 563,976 |
| stable | YES | NO | NO |

**结论：Stick-slip 模型在 Case A 上严重恶化，y_error 是 tangential-damping 的 6.9 倍。**

### Case B（多 patch）

| 指标 | PatchConstit | TangentialDamping | StickSlip (新) |
|------|--------------|-------------------|----------------|
| final_y | 1.1252 m | 1.1058 m | 1.1762 m |
| y_error | 0.0222 m | 0.0416 m | **0.0288 m** ✅ |
| force_std_dev | 487 N | 486 N | 483 N |
| avg_torque_z | 0.034 Nm | -0.153 Nm | -0.024 Nm |
| avg_tangential_force | N/A | 17.3 N | 16.9 N |
| stick_steps | - | - | 7,886,621 |
| slip_steps | - | - | 14,459,843 |
| stable | YES | YES | YES |

**结论：Stick-slip 模型在 Case B 上保持稳定性，y_error 介于 PatchConstit 和 TangentialDamping 之间。**

## Stick-Slip 行为分析

### Case A
- Stick 比例: 39% (360,564 / 924,540)
- Slip 比例: 61% (563,976 / 924,540)
- Stick→Slip 转换: 7,140 次
- 平均切向位移: 0.008 m

### Case B
- Stick 比例: 35% (7,886,621 / 22,346,464)
- Slip 比例: 65% (14,459,843 / 22,346,464)
- Stick→Slip 转换: 44,768 次
- 平均切向位移: 0.001 m

## 为什么 Stick-Slip 在 Case A 上失败？

### 1. 切向位移累积过快
- Case A 中球体高速旋转（ω_z = 5 rad/s），切向速度大
- 切向位移 xi_t 快速累积到 0.008m
- 导致试探力 F_t_trial 过大，频繁进入 slip 状态

### 2. Slip 状态下的力投影问题
- Slip 时 xi_t = F_t / (-k_t) 的投影假设过于简化
- 实际物理中 slip 是连续过程，不是瞬时重置
- 投影导致切向力方向突变，引入数值振荡

### 3. 参数敏感性
- k_t = 1e4 N/m 可能过大，导致切向力响应过激
- mu = 0.3 的摩擦限幅在高速旋转下频繁触发
- 缺少切向力的平滑过渡机制

### 4. 单 patch 局限
- Case A 只有 1 个 patch，切向效应无法分散
- 所有切向力集中在单一 patch 上，放大了模型缺陷

## 为什么 Stick-Slip 在 Case B 上有效？

1. **多 patch 平均效应**: 5 个 patch 分散了切向力
2. **盒体几何**: 平坦接触面使切向速度分布更均匀
3. **低速旋转**: 初始角速度 (0.5, 0, 0.3) rad/s 远小于 Case A
4. **切向位移小**: 平均仅 0.001m，试探力在合理范围

## 结论

**最小 stick-slip 模型在旋转接触场景下失败，但在多 patch 场景下保持有效。**

主要教训：
1. 切向位移累积机制在高速接触下不稳定
2. Slip 状态的力投影过于简化
3. 需要更精细的切向状态管理（如速率依赖、平滑过渡）
4. 单 patch 场景放大了模型缺陷
