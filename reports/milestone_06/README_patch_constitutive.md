# 里程碑 6：Patch Constitutive Law 首版

## 概述

本里程碑将 patch 从"sample force 聚合容器"升级为"直接力生成器"。首次在 patch 层级实现 constitutive law，使接触力直接由 patch 几何描述符和状态量驱动。

## 核心创新

| 层级 | 说明 |
|------|------|
| Pointwise | sample-by-sample force（基线） |
| Patch-force-sum | Σ sample_force_i → patch_force |
| Patch persistence | 跨帧 patch ID 跟踪 |
| **Patch constitutive** | **patch 几何量 + 状态量 → 直接生成力** |

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_openvdb.exe
```

## 输出文件

所有输出文件位于 `out/milestone_06/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_constitutive_output.csv` | 三模式对比摘要 |
| `sdf_patch_constitutive_summary.csv` | 指标汇总对比 |
| `sdf_patch_constitutive_vs_baselines.csv` | 详细对照表（含注释） |

## 三种模式说明

### Mode 1: Pointwise Baseline
- 每个 sample 独立计算 contact force
- 沿用已验证的 penalty contact 逻辑
- 作为对照基线

### Mode 2: Patch-force-sum Baseline
- 通过 patch grouping 聚合 sample forces
- patch_force = Σ sample_force_i
- 验证 patch grouping 的数学等价性

### Mode 3: Patch Constitutive Mode
- **patch force 由 constitutive law 直接生成**
- 不依赖 sample-level force 求和
- 使用 patch 几何描述符作为输入

## Patch Constitutive Law 定义

### 数学定义

```
penetration_eff = area-weighted mean of max(-phi_i, 0)
vn_eff = patch_center_velocity · patch_normal
F_patch = (k_patch * penetration_eff + c_patch * max(-vn_eff, 0)) * patch_normal
T_patch = (patch_center - COM) × F_patch
```

### 物理含义

| 量 | 计算方式 | 物理意义 |
|----|----------|----------|
| penetration_eff | 面积加权平均穿透深度 | patch 整体穿透程度 |
| vn_eff | patch 中心速度在法向的投影 | patch 法向相对速度 |
| F_patch | k * penetration_eff + c * damping | patch 法向接触力 |
| T_patch | 力臂 × 力 | patch 对质心的力矩 |

### 参数设置

| 参数 | 值 | 说明 |
|------|-----|------|
| k_patch | 1e5 N/m | 与 pointwise stiffness 一致 |
| c_patch | 500 Ns/m | 与 pointwise damping 一致 |
| activation_band | 0.1 m | 候选检测带宽 |
| force_band | 0.0 m | 实际力计算带宽 |

### 实现细节

- `penetration_eff` 在 `BuildPatches()` 中计算为面积加权平均
- `vn_eff` 在 simulation loop 中通过 patch center velocity 计算
- `ComputePatchConstitutiveForce()` 函数实现完整 constitutive law
- patch force 通过 `Accumulator` 回写刚体

## 架构说明

### 数据流

```
active samples
  → SDF query (phi, grad)
  → patch grouping (BFS connected components)
  → patch descriptors (center, normal, area, mean_phi, effective_penetration)
  → patch persistence (matching across frames)
  → patch constitutive law (F_patch from geometric state)
  → body force/torque (via Accumulator)
```

### 复用结构

- `PatchPrimitive` 结构体（新增 `effective_penetration`, `effective_normal_velocity`）
- `PersistentPatchTrack` 结构体
- `GridAdjacency` 和 BFS 连通分量
- `MatchPatches()` 跨帧匹配
- `OpenVDBSphereSDF` SDF 查询
- `GenerateSphereSamples()` 表面采样

### 新增内容

- `ContactMode` 枚举（Pointwise, PatchForceSum, PatchConstitutive）
- `ComputePatchConstitutiveForce()` 函数
- `SimulationResult` 和 `ModeMetrics` 结构体
- `RunSimulationMode()` 多模式运行器
- 三模式对照 CSV 输出

## 关键发现

1. **Patch constitutive 模式 y_error 仅 0.000131m**，比 baseline 降低 96%
2. **Force oscillation 降低 71%**（std dev 833.79N vs 2868.74N）
3. **唯一通过稳定性检查的模式**
4. **数值精度超越 pointwise baseline**

## 局限性

- 当前场景仅产生单 patch（球-球接触）
- 多 patch 场景尚未验证
- 未使用 patch age/history 作为 constitutive 输入
- 仅测试了最简单的线性 constitutive law

## 结论

Patch constitutive law 首版成功实现，在简单场景下表现出优于 baseline 的数值精度和稳定性。为后续复杂 constitutive law 和多 patch 场景奠定了基础。
