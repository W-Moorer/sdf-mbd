# Milestone 5: Patch Persistence / Patch Tracking

## 1. 目的

本阶段在 Milestone 4（Patch-level Force Aggregation）的基础上，将 patch 从"每帧重新计算的局部集合"升级为"跨时间步可识别、可匹配、可统计的稳定对象"。

## 2. 核心思想

已有的三层关系：

| 层级 | 描述 | Demo |
|------|------|------|
| 1. Pointwise | Active sample 级力 | `demo_CH_sdf_point_contact_openvdb.cpp` |
| 2. Patch-force | Patch 聚合级力 | `demo_CH_sdf_patch_force_openvdb.cpp` |
| 3. **Patch-tracking**（本里程碑） | 跨时间步状态级 | `demo_CH_sdf_patch_persistence_openvdb.cpp` |

Patch-tracking 在 patch-force 基础上增加：
- **Patch matching**: 相邻时间步间的 patch 对应关系
- **Persistent ID**: 跨帧的唯一标识
- **Birth / Death 检测**: patch 的出现与消失
- **Track 时序记录**: 每个 patch 的完整生命周期

## 3. Matching 规则

### 3.1 匹配算法
**基于 patch_center 最近邻的贪心匹配**：

1. 按 force magnitude 从大到小排序当前帧的 patches
2. 对每个当前 patch，在上一帧的 persistent tracks 中寻找最优匹配
3. 匹配条件：
   - **距离阈值**: center 距离 < 0.5 × dyn_sphere_radius（0.1m）
   - **法向相似度**: cos(angle) > 0.9
4. 每个上一帧 patch 最多匹配一个当前 patch（贪心一一匹配）
5. 无匹配的新 patch → **birth**，分配新 persistent ID
6. 未匹配的上一帧 patch → **death**

### 3.2 为什么使用贪心而非 Hungarian
- 当前场景最多只有 1 个 patch（球-球接触是单点接触）
- 贪心算法在单 patch 场景下与最优匹配等价
- 即使多 patch 场景，贪心按 force 排序后优先匹配大 patch 也是合理的

### 3.3 Persistent ID 分配规则
- 新 patch: 分配 `next_persistent_id++`
- 匹配成功的 patch: 继承上一帧的 persistent ID
- 消失后重新出现: 分配新的 persistent ID（不尝试重识别）

## 4. 输出文件说明

```
out/milestone_06/
├── sdf_patch_tracking_output.csv    (帧级追踪: patch_count, birth_count, dead_count, 刚体状态)
├── sdf_patch_tracks.csv             (Patch track 长表: 每个 patch 的生命周期记录)
└── sdf_patch_tracking_summary.csv   (追踪指标汇总: matched_ratio, lifetime, drift, fluctuation)
```

### sdf_patch_tracking_output.csv
每行一个输出时间步：
- time, patch_count, active_sample_count
- matched_count, born_count, dead_count
- pos_x/y/z, vel_y, total_force_y, force_error

### sdf_patch_tracks.csv（长表格式）
每行一个 patch 在一个时间步的状态：
- time, persistent_patch_id, frame_patch_id, status（born/alive/dead）
- center_x/y/z, normal_x/y/z
- area, mean_phi, min_phi, force_magnitude
- age_in_steps

**death 记录**: 当一个 patch 消失时，额外写入一行 status=dead，记录其最后已知状态。

### sdf_patch_tracking_summary.csv
汇总指标：
| 指标 | 说明 |
|------|------|
| total_born | 整个仿真中 patch 出生总数 |
| total_dead | 整个仿真中 patch 死亡总数 |
| total_frames | 总仿真帧数 |
| frames_with_patches | 有 patch 存在的帧数 |
| avg_patch_count | 平均每帧 patch 数 |
| max_patch_count | 最大 patch 数 |
| avg_matched_ratio | 平均匹配率 = matched / max(prev_count, curr_count) |
| max_lifetime_steps | 最长 patch 生命周期（步数） |
| avg_lifetime_steps | 平均 patch 生命周期（步数） |
| max_center_drift | 最大累计 center 漂移 |
| max_normal_fluctuation | 最大累计 normal 波动 |

## 5. 如何构建

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver"
cmake --build build --config Release --target demo_CH_sdf_patch_persistence_openvdb
```

## 6. 如何运行

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_patch_persistence_openvdb.exe
```
