# Milestone 3: Patch Primitive 首版实现

## 1. 目的

本阶段在 Milestone 2B.6 的 OpenVDB 点式接触基础上，引入第一版 patch primitive，使接触从"离散 active sample 力累加"升级为"patch 级组织与 patch 级输出"。

## 2. 核心思想

不是替换 pointwise 接触，而是在其上增加一层：

```
active sample set -> spatial/graph grouping -> patch primitive -> patch descriptors -> patch-level statistics
```

第一版 patch primitive 的目标是：
- patch 级表示
- patch 级统计
- patch 级时序分析
- 对照 pointwise 和 patch 的稳定性

## 3. 场景定义

### 静态目标
- **OpenVDB sphere level set**: 中心在原点，半径 1.0m，voxel size 0.05m，half width 3 voxels（band width ±0.15m）

### 动态刚体
- **Sphere**: 半径 0.2m，密度 1000 kg/m³，质量 33.51 kg
- **初始位置**: y = 3.0m 自由落体

### 理论接触位置
- 动态球心理论平衡 Y = 1.0 + 0.2 = 1.2m

### 接触参数（沿用 2B.6 最佳）
- stiffness = 1e5 N/m, damping = 500 N*s/m, force_band = 0.0m, activation_band = 0.1m

## 4. Patch 数据结构

```cpp
struct PatchPrimitive {
    int patch_id;                      // 唯一标识
    std::vector<int> sample_ids;       // 包含的 active sample 索引
    int active_sample_count;           // active sample 数量
    ChVector3d center;                 // 加权平均位置
    ChVector3d normal;                 // 加权平均法向（归一化）
    double total_area;                 // sample weight 之和
    double min_phi;                    // 最小 phi 值
    double mean_phi;                   // 加权平均 phi
    double max_penetration;            // 最大穿透量
    ChVector3d force;                  // patch 内 pointwise 力之和
    ChVector3d torque;                 // patch 内 pointwise 力矩之和
};
```

## 5. Patch Grouping 策略

### Adjacency 定义
基于球面参数采样网格的邻接关系：
- 采样网格为 (n_theta=8) × (n_phi=16)，共 128 个 sample
- 相邻定义：sample (i,j) 与 (i±1,j) 和 (i,j±1) 相邻
- Phi 方向采用 wrap-around（首尾相连）
- Theta 方向不 wrap-around（两极不相连）

### Grouping 算法
1. 找出所有 active sample 的索引集合
2. 基于邻接图寻找连通分量（BFS）
3. 每个连通分量成为一个 patch

### 示例结果
在球-球接触场景下，接触区域通常是连续的，因此大多数情况下只形成 **1 个 patch**（所有 active sample 连通）。

## 6. 如何构建

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver"
cmake --build build --config Release --target demo_CH_sdf_patch_contact_openvdb
```

## 7. 如何运行

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_patch_contact_openvdb.exe
```

## 8. 输出目录

```
out/milestone_04/
├── sdf_patch_contact_output.csv              (主仿真输出：时间序列)
├── sdf_patch_summary.csv                     (Patch 摘要：长表格式，每行一个 patch)
├── sdf_patch_vs_pointwise.csv                (Pointwise 对照输出)
└── sdf_patch_vs_pointwise_summary.csv        (Patch vs Pointwise 对比指标)

reports/milestone_04/
├── README_patch_contact.md   (本报告)
└── acceptance_report.md      (验收报告)
```

## 9. 输出文件说明

### sdf_patch_contact_output.csv
每行一个时间步输出：
- time, pos_x/y/z, vel_x/y/z
- expected_y, y_error, normalized_y_error
- active_sample_count, patch_count, total_patch_area
- total_force_x/y/z, force_error

### sdf_patch_summary.csv
长表格式，每行一个 patch：
- time, patch_id
- center_x/y/z, normal_x/y/z
- area, min_phi, mean_phi, max_penetration
- active_sample_count
- force_x/y/z, force_magnitude

当无 patch 时，patch_id = -1，表示无接触。

### sdf_patch_vs_pointwise.csv
Pointwise 版本的时间序列输出（与 patch 对比用）。

### sdf_patch_vs_pointwise_summary.csv
对比指标汇总表：
- final_y, y_error, avg_active_count, avg_patch_count
- avg_max_patch_samples, force_std_dev
- center_drift_total, normal_fluctuation_total
- stable 标志

## 10. 设计决策

### 为什么新建独立 demo 而非复用 pointwise demo
Patch primitive 需要在 pointwise 基础上增加 patch grouping、patch 描述符计算和时序追踪等大量新逻辑。复用会导致代码混乱。新 demo 保留 pointwise 力回写逻辑，仅添加 patch 观测层，保持核心力学模型不变。

### 为什么 patch_id 全局递增
为支持 patch 持久性追踪（persistence tracking），每个新生成的 patch 分配唯一 ID。未来可基于 ID 追踪同一 patch 的生命周期。

### 为什么第一版不替换 pointwise 力学
patch primitive 第一版的目标是验证表示和组织能力，不急于改变力学模型。保持 pointwise 力回写确保动力学行为与 2B.6 一致，patch 层仅作为观测和分析工具。
