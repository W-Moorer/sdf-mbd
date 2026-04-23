# 里程碑 12：Transport-Limited 模型消融实验 / 最优简化版筛选

## 概述

本里程碑对 milestone 11 中的三种 transport 限制机制进行严格的消融实验，分离各机制对 Case A / Case B 的贡献，找到"最小但有效"的 transport-limited 版本。

## Milestone 10/11 Baseline 记录

### Case A（旋转接触）

| 里程碑 | 模式 | y_error | torque_z | transported_norm | stable |
|--------|------|---------|----------|------------------|--------|
| MS10 | FrameConsistent (B0) | 0.553 m | +2.30 Nm | 5.62 m | NO |
| MS11 | Full (A+B+C) | 0.074 m | +0.63 Nm | 0.12 m | NO |

### Case B（多 patch）

| 里程碑 | 模式 | y_error | torque_z | transported_norm | stable |
|--------|------|---------|----------|------------------|--------|
| MS10 | FrameConsistent (B0) | 0.029 m | -0.023 Nm | 1.96 m | YES |
| MS11 | Full (A+B+C) | 0.069 m | -0.12 Nm | 0.13 m | YES |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_transport_ablation_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_transport_ablation_openvdb.exe
```

## 消融开关说明

三种机制可独立启用/禁用：

| 开关 | 机制 | 描述 |
|------|------|------|
| A | Magnitude Clamp | `if |xi| > 0.01m, clamp to 0.01m` |
| B | Rotation-Aware Attenuation | `attenuation = exp(-5 * theta)` |
| C | Partial Carry-Over | `xi *= 0.5`（每次只继承 50%） |

## 消融配置

| 配置 | A (Clamp) | B (Attenuation) | C (Carry-Over) | 描述 |
|------|-----------|-----------------|----------------|------|
| B0 | ✗ | ✗ | ✗ | Frame-consistent baseline (无限制) |
| A-only | ✓ | ✗ | ✗ | 仅 clamp |
| B-only | ✗ | ✓ | ✗ | 仅 attenuation |
| C-only | ✗ | ✗ | ✓ | 仅 carry-over |
| A+B | ✓ | ✓ | ✗ | clamp + attenuation |
| A+C | ✓ | ✗ | ✓ | clamp + carry-over |
| B+C | ✗ | ✓ | ✓ | attenuation + carry-over |
| Full | ✓ | ✓ | ✓ | 三种全开 |

## 实验结果

### Case A（旋转接触）消融结果

| 配置 | y_error | torque_z | trans_norm | clamp_hits | stable |
|------|---------|----------|------------|------------|--------|
| B0 (baseline) | 0.5916 m | +2.15 Nm | 5.55 m | 0 | NO |
| A-only | 0.3795 m | +1.45 Nm | 3.56 m | 419107 | NO |
| B-only | 0.1121 m | -1.39 Nm | 7.96 m | 0 | NO |
| **C-only** | **0.0750 m** | +0.53 Nm | 0.12 m | 0 | NO |
| A+B | 0.3273 m | +0.19 Nm | 3.61 m | 399152 | NO |
| A+C | 0.0750 m | +0.53 Nm | 0.12 m | 0 | NO |
| **B+C** | **0.0743 m** | +0.62 Nm | 0.12 m | 0 | NO |
| Full (A+B+C) | 0.0743 m | +0.62 Nm | 0.12 m | 0 | NO |

### Case B（多 patch）消融结果

| 配置 | y_error | torque_z | trans_norm | clamp_hits | stable |
|------|---------|----------|------------|------------|--------|
| B0 (baseline) | 0.0458 m | -0.020 Nm | 2.66 m | 0 | YES |
| **A-only** | **0.0311 m** | -0.025 Nm | 2.01 m | 83 | YES |
| B-only | 0.0345 m | -0.014 Nm | 2.34 m | 0 | YES |
| C-only | 0.0740 m | -0.121 Nm | 0.13 m | 0 | YES |
| **A+B** | **0.0320 m** | -0.026 Nm | 2.08 m | 45 | YES |
| A+C | 0.0740 m | -0.121 Nm | 0.13 m | 0 | YES |
| B+C | 0.0726 m | -0.120 Nm | 0.13 m | 0 | YES |
| Full (A+B+C) | 0.0726 m | -0.120 Nm | 0.13 m | 0 | YES |

## 各机制贡献分析

### 机制 A：Magnitude Clamp
- **Case A**: 改善 y_error 36% (0.592→0.380 m)，transport_norm 减少 36%
- **Case B**: 改善 y_error 32% (0.046→0.031 m)，是 Case B 最佳单机制
- **特征**: clamp 命中次数高（Case A: 419k, Case B: 83），说明经常触发
- **副作用**: 单独使用不足以达到最优

### 机制 B：Rotation-Aware Attenuation
- **Case A**: 改善 y_error 81% (0.592→0.112 m)，但 transport_norm **增加 43%** (5.55→7.96 m)
- **Case B**: 改善 y_error 25% (0.046→0.035 m)
- **特征**: 不触发 clamp，仅做软衰减
- **关键发现**: 单独使用时 transport_norm 反而增大，说明衰减不够强或方向问题

### 机制 C：Partial Carry-Over
- **Case A**: 改善 y_error 87% (0.592→0.075 m)，transport_norm **减少 98%** (5.55→0.12 m)
- **Case B**: 恶化 y_error 62% (0.046→0.074 m)
- **特征**: 对 transport_norm 抑制最有效
- **关键发现**: 是 Case A 改善的主导因素，但对 Case B 有副作用

## 关键发现

### 1. Carry-Over (C) 是 Case A 的主导改善因素
- C-only 在 Case A 上 y_error = 0.075 m，与 Full 几乎相同 (0.074 m)
- transport_norm 从 5.55 m 降至 0.12 m（减少 98%）
- A 和 B 对 C 的贡献几乎无额外增益（A+C ≈ C-only，B+C ≈ C-only）

### 2. Clamp (A) 是 Case B 的主导改善因素
- A-only 在 Case B 上 y_error = 0.031 m，是 Case B 最佳配置
- B-only 也改善 Case B (0.035 m)，但不如 A
- C 单独使用反而恶化 Case B (0.074 m)

### 3. Full 模型不是最优
- Case A: Full ≈ B+C ≈ C-only ≈ A+C (都在 0.074-0.075 m)
- Case B: Full ≈ B+C ≈ C-only ≈ A+C (都在 0.073-0.074 m)，远不如 A-only (0.031 m)
- 说明添加更多机制并不总是更好

### 4. Case A vs Case B 的根本矛盾
- **Case A 需要强 carry-over**（α=0.5）来抑制 transport 放大
- **Case B 不需要强 carry-over**，因为法向变化小，carry-over 反而过度衰减了有效历史
- **Clamp 对 Case B 有正面效果**，但对 Case A 只能提供中等改善

## 输出文件

所有输出文件位于 `out/milestone_12/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_transport_ablation_case_A.csv` | Case A 八配置详细输出 |
| `sdf_patch_transport_ablation_case_B.csv` | Case B 八配置详细输出 |
| `sdf_patch_transport_ablation_summary.csv` | 跨场景八配置汇总 |
| `sdf_patch_transport_ablation_compare.csv` | 指标对比矩阵 |
