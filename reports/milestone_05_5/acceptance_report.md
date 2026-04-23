# Milestone 5 验收报告

## 验收日期
2026-04-23

## 验收人员
AI Agent

---

## 1. 验收项目检查清单

| # | 验收项目 | 状态 | 证据 |
|---|---------|------|------|
| 1 | 工程能成功构建 | PASS | 编译成功，无警告 |
| 2 | demo 能成功运行 | PASS | 退出码 0，输出 "Overall: PASS" |
| 3 | 输出文件正确生成在 `out/milestone_06/` | PASS | 3 个 CSV 文件已生成 |
| 4 | 成功输出 patch tracking 数据 | PASS | `sdf_patch_tracks.csv` 含完整生命周期记录 |
| 5 | persistent patch ID 跨帧稳定出现 | PASS | avg_matched_ratio = 0.9947 |
| 6 | 正确识别 birth / death 行为 | PASS | 3 births, 2 deaths 正确检测 |
| 7 | tracking summary 指标可读 | PASS | 13 项指标全部有值 |
| 8 | 文件归档正确 | PASS | 无散落临时文件 |

---

## 2. Patch Tracking 定义与实现

### 2.1 Persistent ID

| 规则 | 实现 |
|------|------|
| 新 patch | `next_persistent_id++`（从 0 开始递增） |
| 匹配成功 | 继承上一帧 persistent ID |
| 消失后重现 | 分配新 ID（不尝试重识别） |

### 2.2 Matching 规则

| 条件 | 阈值 | 说明 |
|------|------|------|
| Center 距离 | < 0.1m（0.5 × dyn_sphere_radius） | 接触区域不应突变 |
| 法向相似度 | cos(angle) > 0.9 | 接触面方向不应突变 |
| 匹配策略 | 贪心（按 force 降序） | 单 patch 场景下与最优等价 |
| 一一匹配 | 每个 prev patch 最多匹配一个 current | 避免多对一 |

### 2.3 Birth / Death 规则

| 事件 | 条件 | 处理 |
|------|------|------|
| **Birth** | 当前 patch 无匹配的上一帧 patch | 分配新 persistent ID，status = born |
| **Alive** | 当前 patch 匹配到上一帧 patch | 继承 ID，status = alive |
| **Death** | 上一帧 patch 无匹配的当前 patch | 写入 death 记录，status = dead |

### 2.4 Track 输出格式

**长表格式**（time, persistent_patch_id, ...），每个 patch 在每个时间步一条记录。
- Born/Alive: 正常写入
- Dead: 额外写入一行（记录最后状态）

---

## 3. Tracking 结果摘要

### 3.1 核心指标

| 指标 | 值 | 说明 |
|------|-----|------|
| total_born | 3 | 3 次接触建立 |
| total_dead | 2 | 2 次接触消失 |
| total_frames | 4001 | 仿真总帧数 |
| frames_with_patches | 940 | 23.5% 帧有接触 |
| avg_patch_count | 0.2349 | 平均不到 1 个 patch（单接触区域） |
| max_patch_count | 1 | 最多同时只有 1 个 patch |
| **avg_matched_ratio** | **0.9947** | 匹配率极高，tracking 非常稳定 |
| max_lifetime_steps | 621 | 最长 track 存活 621 步 |
| avg_lifetime_steps | 207 | 平均存活 207 步 |
| max_center_drift | 0.3374 m | 累计 center 漂移 |
| max_normal_fluctuation | 0.0 | **Normal 完全稳定** |

### 3.2 Patch 生命周期记录

| Track ID | Birth Time | Death Time | Lifetime (steps) | 说明 |
|----------|-----------|------------|-----------------|------|
| **0** | t=0.590s | t=0.6495s | 118 | 第一次接触，弹起 |
| **1** | t=1.290s | t=1.3885s | 198 | 第二次接触，弹起 |
| **2** | t=1.690s | (still alive at t=2.0s) | 621+ | 第三次接触，持续接触 |

### 3.3 Birth / Death 行为验证

**Birth 正确性**：
- Track 0 在 t=0.59s birth，此时动态球首次接触静态球（自由落体约 0.6s 到达接触面），正确
- Track 1 在 t=1.29s birth，此时球从第一次弹起后回落到接触面，正确
- Track 2 在 t=1.69s birth，第三次接触，正确

**Death 正确性**：
- Track 0 在 t=0.6495s death，此时球从接触面弹起，active sample 消失，正确
- Track 1 在 t=1.3885s death，第二次弹起，正确

**Continuous contact**：
- Track 2 从 t=1.69s 持续到仿真结束（t=2.0s），说明球最终稳定在接触状态，符合 2B.6 的平衡结论

### 3.4 是否稳定

**是。** avg_matched_ratio = 0.9947 说明几乎每一帧的 patch 都能正确匹配到上一帧。在单 patch 场景下这是预期的——只要 contact 持续存在，patch center 变化很小（< 0.1m 阈值），normal 始终为 (0,1,0)。

---

## 4. Pointwise / Patch-force / Patch-tracking 层级关系

```
Level 1: Pointwise (active sample 级)
         └─ 每个 sample 点有独立 force/torque
         └─ 直接累加到刚体

Level 2: Patch-force (patch 聚合级)
         └─ sample → patch grouping → patch_force = Σ(sample_force)
         └─ 输出结构化的 patch 描述符

Level 3: Patch-tracking (跨时间步状态级)  ← 本里程碑
         └─ 相邻帧 patch matching → persistent ID
         └─ birth / death 检测
         └─ track 时序记录（center drift, normal fluctuation, lifetime）
```

这三个层级是**递进关系**，每一层在前一层的基础上增加新的抽象：
- Level 1 是物理基础
- Level 2 是空间组织
- Level 3 是时间组织

Persistence 不是另一个不相关模块，而是 patch-force 层在时间维度上的自然扩展。

---

## 5. 当前局限

### 5.1 场景过于简单
- 当前只有单接触区域（球-球接触），始终最多 1 个 patch
- Matching 算法在单 patch 下几乎是 trivial 的
- 需要多接触区域场景（如 edge contact）验证 matching 的 robustness

### 5.2 未实现 patch 重识别
- 消失后重新出现的 patch 分配新 ID，不尝试重识别
- 在某些场景中，同一个接触区域可能在间隙后重新出现
- 但当前场景中这不会造成问题（每个接触区域确实是不同的事件）

### 5.3 未实现 split / merge 检测
- 当前 patch grouping 使用连通分量，如果接触区域从连通变为不连通，会创建新 patch
- 但不会显式检测 "split"（一个 patch 分为两个）或 "merge"（两个 patch 合并为一个）
- 可以通过分析 parent/child persistent ID 关系来实现，但不是当前优先任务

### 5.4 Center drift 阈值可能需要场景自适应
- 当前阈值 = 0.5 × dyn_sphere_radius = 0.1m
- 对于更复杂的几何或更大的接触区域，可能需要自适应阈值

---

## 6. 最终文件结构

```
src/demos/core/
├── demo_CH_sdf_point_contact.cpp                    (2B.5 解析平面)
├── demo_CH_sdf_point_contact_openvdb.cpp            (2B.6 OpenVDB pointwise)
├── demo_CH_sdf_patch_contact_openvdb.cpp            (3.0 OpenVDB patch observational)
├── demo_CH_sdf_patch_force_openvdb.cpp              (4.0 OpenVDB patch force)
├── demo_CH_sdf_patch_persistence_openvdb.cpp        (5.0 OpenVDB patch tracking，新增)
└── CMakeLists.txt                                    (更新)

out/milestone_06/
├── sdf_patch_tracking_output.csv                    (帧级追踪)
├── sdf_patch_tracks.csv                             (Patch track 长表)
└── sdf_patch_tracking_summary.csv                   (追踪指标汇总)

reports/milestone_06/
├── README_patch_tracking.md                          (说明文档)
└── acceptance_report.md                              (本报告)
```

---

## 7. 如何构建

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver"
cmake --build build --config Release --target demo_CH_sdf_patch_persistence_openvdb
```

## 8. 如何运行

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_patch_persistence_openvdb.exe
```

---

## 9. 验收结论

### 里程碑 5 状态：**PASSED** ✅

所有 8 项验收标准均已满足。Patch persistence / tracking 首版实现成功验证了：

1. **Patch 可以跨时间步匹配** ✅（avg_matched_ratio = 0.9947）
2. **Persistent ID 稳定延续** ✅（3 tracks，每个都有连续的 ID）
3. **Birth / Death 正确检测** ✅（3 births at correct times, 2 deaths at lift-off）
4. **Track 时序指标可读** ✅（center drift, normal fluctuation, lifetime 全部有值）
5. **Tracking 机制足够稳定** ✅（normal 完全稳定，center drift 可接受）

### 是否适合进入 patch constitutive law

**可以进入 ✅**

**支撑理由：**
1. ✅ Patch tracking 已验证稳定性，可以作为 patch constitutive law 的状态承载体
2. ✅ Persistent ID 确保同一 patch 的历史状态可以累积（如接触深度历史、力历史）
3. ✅ Birth / death 检测确保 constitutive law 可以在接触建立时初始化、在消失时清理

**建议的 constitutive law 阶段优先事项：**
1. **基于 patch area + mean_depth 的等效 Hertz 力模型**：用 patch descriptor 替代 pointwise force 聚合
2. **利用 track 历史**：如接触深度变化率、接触面积变化趋势
3. **平滑机制**：利用 persistent track 的时序连续性来平滑 patch force
