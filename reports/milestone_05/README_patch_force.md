# Milestone 4: Patch-level Force Aggregation

## 1. 目的

本阶段在 Milestone 3（Patch Primitive 首版）的基础上，将 patch 从"观测层"升级为"力学层"，使 patch 真实接管力和力矩的聚合与回写。

## 2. 核心思想

当前已有的三套 demo：

| Demo | 模式 | 用途 |
|------|------|------|
| `demo_CH_sdf_point_contact_openvdb.cpp` | Pointwise force baseline | 逐点力累加到刚体 |
| `demo_CH_sdf_patch_contact_openvdb.cpp` | Patch observational baseline | Patch 描述符输出，但力仍逐点回写 |
| `demo_CH_sdf_patch_force_openvdb.cpp` | **Patch aggregated force**（本里程碑） | Patch 聚合力回写，与 pointwise 做 A/B 对照 |

## 3. 三种接触模式

### Mode A: Pointwise（点式模式，基准）
```
active samples -> sample forces -> body force/torque
```
每个 sample 的 force 和 torque 直接累加到刚体上，与 2B.6 完全相同。

### Mode B: Direct Patch Sum（直接求和模式）
```
active samples -> patch grouping -> patch_force = Σ(sample_force_i)
                               -> patch_torque = Σ(sample_torque_i)
                               -> body force/torque
```
sample force 在 patch 内部求和形成 patch_force 和 patch_torque，然后将所有 patch 的 force/torque 求和后作用到刚体上。

**数学性质**：与 Mode A 完全等价，因为加法满足结合律。此模式的价值在于**输出结构化**而非改变力学结果。

### Mode C: Equivalent Patch（等效作用点模式）
```
active samples -> patch grouping -> patch_center, patch_normal
                               -> force_magnitude = patch_force · patch_normal
                               -> equiv_force = force_magnitude * patch_normal
                               -> equiv_torque = (patch_center - COM) × equiv_force
                               -> body force/torque
```
将每个 patch 简化为一个等效作用点：
- **作用点** = patch_center（加权平均位置）
- **方向** = patch_normal（加权平均法向）
- **力大小** = patch force 在 patch normal 方向的分量（只保留法向力）
- **力矩** = (patch_center - 质心) × equiv_force

此模式引入了近似：忽略了 patch 内的切向力分量和力分布的几何效应。

## 4. Patch Force / Torque 定义

### patch_force
```
patch_force = Σᵢ sample_force_i
```
patch 内所有 sample 的 pointwise 力的向量和。

**参考点**：无。Force 是自由向量，不依赖参考点。

### patch_torque
```
patch_torque = Σᵢ sample_torque_i
             = Σᵢ (r_i × sample_force_i)
```
patch 内所有 sample 的 pointwise 力矩的向量和。

**参考点**：每个 sample_torque_i 是关于刚体质心（COM）计算的，即 `r_i = sample_world_position - body_com_position`。因此 patch_torque 也是关于刚体质心的。

## 5. 输出文件说明

```
out/milestone_05/
├── sdf_pointwise_output.csv              (Mode A 主输出)
├── sdf_patch_force_output.csv            (Mode B 主输出)
├── sdf_patch_force_summary.csv           (Mode B patch 长表)
├── sdf_patch_force_equiv_output.csv      (Mode C 主输出)
└── sdf_patch_force_vs_pointwise.csv      (三模式对比指标表)
```

### sdf_patch_force_vs_pointwise.csv
三模式（pointwise / direct_patch_sum / equivalent_patch）的指标对比：
- final_y, y_error, avg_active_count, avg_patch_count
- force_std_dev, torque_std_dev
- center_drift, normal_fluctuation
- stable

## 6. 如何构建

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver"
cmake --build build --config Release --target demo_CH_sdf_patch_force_openvdb
```

## 7. 如何运行

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_patch_force_openvdb.exe
```
