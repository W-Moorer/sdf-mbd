# Milestone 3 验收报告

## 验收日期
2026-04-23

## 验收人员
AI Agent

---

## 1. 验收项目检查清单

| # | 验收项目 | 状态 | 证据 |
|---|---------|------|------|
| 1 | 工程能成功构建 | PASS | `cmake --build build --config Release --target demo_CH_sdf_patch_contact_openvdb` 编译成功 |
| 2 | demo 能成功运行 | PASS | `.\demo_CH_sdf_patch_contact_openvdb.exe` 正常退出，输出 "PASS" |
| 3 | 输出文件正确生成在 `out/milestone_04/` | PASS | 4 个 CSV 文件已生成 |
| 4 | 成功输出 patch-level 数据 | PASS | `sdf_patch_summary.csv` 含 48 行 patch 记录 |
| 5 | 成功输出 pointwise vs patch 对照 | PASS | `sdf_patch_vs_pointwise_summary.csv` |
| 6 | patch 数量、center、normal 时序指标可读 | PASS | 见 §2 和 §3 |
| 7 | 仿真不崩溃、不明显发散 | PASS | stable=YES, y_error=0.0003m |
| 8 | 文件归档正确 | PASS | 无散落临时文件 |

---

## 2. Patch 定义与实现

### 2.1 Patch 数据结构

```cpp
struct PatchPrimitive {
    int patch_id;                      // 唯一标识（全局递增）
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

### 2.2 Patch Adjacency

基于球面参数采样网格 `(n_theta=8) × (n_phi=16)` 的邻接关系：

| 方向 | 邻接规则 |
|------|---------|
| Theta | (i,j) 与 (i±1,j) 相邻，不 wrap-around |
| Phi   | (i,j) 与 (i,j±1) 相邻，wrap-around（首尾相连） |

### 2.3 Grouping 算法

1. 收集所有 active sample 的索引（phi < activation_band）
2. 在邻接图上执行 BFS 找连通分量
3. 每个连通分量成为一个 patch

### 2.4 Patch 时序指标

| 指标 | 含义 | 值 |
|------|------|-----|
| avg_patch_count | 平均 patch 数量 | 0.24（大部分时间步无接触，约 25% 时间步有 1 个 patch） |
| avg_max_patch_samples | 最大 patch 的平均 sample 数 | 6.69 |
| center_drift_total | 最大 patch center 的累计漂移 | 1.62m |
| normal_fluctuation_total | 最大 patch normal 的累计波动 | 0.0 |

**关键观察：**
- 球-球接触是单点接触，始终只有 1 个 patch
- Normal 完全稳定在 (0,1,0)，因为对称球面接触的法向始终垂直向上
- Center drift 主要来源于弹跳过程中的位置变化，而非接触区域的不稳定

### 2.5 Patch 几何量示例

从 `sdf_patch_summary.csv` 中提取的典型 patch 数据：

| 时间 | patch_id | center_y | normal_y | area | active_samples | force |
|------|---------|----------|----------|------|----------------|-------|
| 0.60s | 20 | 1.054 | 1.000 | 0.250 | 32 | 0 N |
| 0.61s | 40 | 1.024 | 1.000 | 0.375 | 48 | 38778 N |
| 1.32s | 181 | 1.035 | 1.000 | 0.250 | 32 | 0 N |
| 1.33s | 201 | 1.035 | 1.000 | 0.375 | 48 | 31550 N |
| 1.75s | 438 | 1.032 | 1.000 | 0.375 | 48 | 13629 N |
| 2.00s | 938 | 1.020 | 1.000 | 0.250 | 32 | 0 N |

**分析：**
- Patch area 在 0.125~0.375 之间变化，对应 16~48 个 active sample
- Force 仅在 phi < 0（真正穿透）时产生，此时 active_sample_count = 48
- Normal 始终为 (0,1,0)，与球-球接触的几何对称性一致

---

## 3. Pointwise vs Patch 对照结果

### 3.1 核心指标对比

| 指标 | Patch 版本 | Pointwise 版本 | 差异 |
|------|-----------|----------------|------|
| final_y | 1.200327 m | 1.200327 m | 0 |
| y_error | 0.000327 m | 0.000327 m | 0 |
| avg_active_count | 6.69 | 6.69 | 0 |
| avg_patch_count | 0.24 | N/A | - |
| force_std_dev | 3639.35 N | 3639.35 N | 0 |
| stable | YES | YES | - |

**说明：** patch 版本和 pointwise 版本结果完全一致，这是因为第一版 patch primitive 仅作为观测层，不改变力回写逻辑。这是正确的设计——patch 层验证通过后，下一阶段的 patch-level force aggregation 才有可信的基准。

### 3.2 力时序对比

两版本的力时序完全相同（force_std_dev = 3639.35 N），说明：
- Patch 层的引入没有改变任何动力学行为
- 力的波动来自于 bounce 过程中的瞬态响应（damping-dominated 冲击），而非 patch 组织
- 在稳定平衡阶段（t > 1.5s），力趋于 0（因为 force_band = 0.0，只有穿透才施力）

### 3.3 哪个更稳定、哪个更可解释

**稳定性**：两者相同（因为力学模型未变）

**可解释性**：Patch 版本显著更优
- Pointwise：128 个独立 sample 点的力和 phi 值
- Patch：1 个 patch，包含 center、normal、area、force 等聚合描述符
- Patch 描述符直接回答"接触区域在哪里、多大、朝哪、力多大"

---

## 4. 当前局限

### 4.1 Patch 分组策略简单
- 当前使用基于网格索引的连通分量，对于球-球接触（单接触区域）效果良好
- 对于多接触区域（如边缘接触、多面接触），连通分量可能将物理上不相关的样本连在一起

### 4.2 未实现 patch persistence 追踪
- 当前 patch_id 全局递增，每个时间步新建的 patch 都有新 ID
- 未实现跨时间步的 patch 匹配（如基于 center 最近邻或 sample 重叠率）
- 因此无法追踪同一 patch 的生命周期和演化

### 4.3 未实现 patch-level force aggregation
- 当前仍使用 pointwise force 回写到刚体
- 未来的 patch force model 可以：
  - 基于 patch 几何量计算等效 Hertz 接触力
  - 平滑 patch 间的力分布
  - 减少力震荡

### 4.4 力震荡问题仍然存在
- force_std_dev = 3639 N，约为重力 (329 N) 的 11 倍
- 这主要是高刚度 + bounce 过程中的阻尼冲击导致的
- Patch-level force model 可以帮助平滑，但不是 patch primitive 层的问题

---

## 5. 最终文件结构

```
src/demos/core/
├── demo_CH_sdf_point_contact.cpp            (2B.5 解析平面参考)
├── demo_CH_sdf_point_contact_openvdb.cpp    (2B.6 OpenVDB pointwise)
├── demo_CH_sdf_patch_contact_openvdb.cpp    (3.0 OpenVDB patch primitive)
└── CMakeLists.txt                            (更新)

out/milestone_04/
├── sdf_patch_contact_output.csv             (主仿真输出)
├── sdf_patch_summary.csv                    (Patch 摘要)
├── sdf_patch_vs_pointwise.csv               (Pointwise 对照)
└── sdf_patch_vs_pointwise_summary.csv       (对比指标)

reports/milestone_04/
├── README_patch_contact.md
└── acceptance_report.md                      (本报告)
```

---

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

---

## 8. 验收结论

### 里程碑 3 状态：**PASSED** ✅

所有 8 项验收标准均已满足。Patch primitive 首版实现成功验证了：

1. **Active sample 能组织成 patch primitive** ✅（BFS 连通分量）
2. **Patch 级几何量可计算** ✅（center, normal, area, phi stats, force/torque）
3. **Patch 时序指标可读** ✅（patch count, center drift, normal fluctuation）
4. **Patch 与 pointwise 对照一致** ✅（力学模型不变，结果完全相同）
5. **仿真稳定，y_error < 0.001m** ✅

### 是否适合进入 patch-level force model

**条件性可以进入。**

**支撑理由：**
1. Patch 数据结构稳定，所有描述符计算正确
2. Grouping 策略在单接触区域场景下工作良好
3. Patch 时序行为可解释（center drift 来自弹跳，normal 始终稳定）

**需要注意：**
1. Patch persistence 追踪尚未实现，patch-level force model 需要跨时间步的 patch 状态
2. 力震荡问题需要在 patch force model 阶段解决（建议引入 patch-level 平滑）
3. 当前的连通分量策略在多接触区域场景下可能需要改进（如距离阈值聚类）

**建议：可以开始设计 patch-level force model，但需同步实现 patch persistence 追踪。**
