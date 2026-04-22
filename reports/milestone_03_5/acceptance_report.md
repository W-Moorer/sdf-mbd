# Milestone 2B.5 验收报告

## 验收日期
2026-04-22

## 验收人员
AI Agent

---

## 1. 验收项目检查清单

| # | 验收项目 | 状态 | 证据 |
|---|---------|------|------|
| 1 | 工程能成功构建 | PASS | `cmake --build build --config Release --target demo_CH_sdf_point_contact` 编译成功 |
| 2 | 现有 demo 能成功运行 | PASS | `.\demo_CH_sdf_point_contact.exe` 正常退出，输出 "Overall: PASS" |
| 3 | 生成参数扫描输出文件 | PASS | `out/milestone_03_5/sdf_point_contact_sweep.csv` (45 行数据) |
| 4 | 理论平衡位置已定义 | PASS | expected_y = 0.200000 m |
| 5 | 最终位置误差已报告 | PASS | y_error = 0.004969 m |
| 6 | 推荐参数窗口已报告 | PASS | stiffness=[5e4, 2e5], damping=[300, 600], force_band=[0.0, 0.01] |
| 7 | 几何误差较 2B 明显下降 | PASS | 2B: 0.117m → 2B.5: 0.005m (下降 96%) |
| 8 | 文件归档正确 | PASS | 见下方文件结构 |

---

## 2. 文件结构验证

### 源代码
- [x] `src/demos/core/demo_CH_sdf_point_contact.cpp` (主 demo)
- [x] `src/demos/core/CMakeLists.txt` (更新)

### 输出文件
- [x] `out/milestone_03_5/sdf_point_contact_sweep.csv`
- [x] `out/milestone_03_5/sdf_point_contact_best_run.csv`

### 报告文件
- [x] `reports/milestone_03_5/README_geometry_consistency.md`
- [x] `reports/milestone_03_5/acceptance_report.md` (本报告)

### 临时文件清理
- [x] 无散落临时文件
- [x] CSV 输出均在 `out/milestone_03_5/` 目录下
- [x] 报告均在 `reports/milestone_03_5/` 目录下

---

## 3. 关键指标对比

| 指标 | Milestone 2B | Milestone 2B.5 | 改善幅度 |
|------|-------------|----------------|---------|
| 最终 Y 位置 | 0.317 m | 0.195 m | - |
| 理论 Y 位置 | 0.200 m | 0.200 m | - |
| Y 误差 | 0.117 m | 0.005 m | **95.7% 下降** |
| 归一化 Y 误差 | 0.585 | 0.025 | **95.7% 下降** |
| 接触力 (最终) | 328.7 N | 波动 (均值 ≈ mg) | - |
| 力的来源 | damping 主导 | stiffness 主导 | **本质改善** |
| SDF 类型 | OpenVDB 大球 | 解析平面 | **本质改善** |
| 参数可重复性 | 单点 | 窗口 | **本质改善** |

---

## 4. 根因分析总结

### 主要问题（2B 中的位置误差）
**activation delta / force_band 过大导致提前托住**

在 Milestone 2B 中，force_band = 0.3m 意味着当球的 sample 点距离表面还有 0.3m 时就开始计算力。虽然此时 `penetration = 0`，但 damping 项 `damping * max(-vn, 0)` 仍然产生向上的力。当这个阻尼力与重力平衡时，球就"悬停"在了空中，从未真正接触表面。

### 解决方案（2B.5）
1. **解析平面 SDF**：消除大球曲率带来的几何不确定性
2. **两带分离**：activation_band (0.1m) 仅用于候选检测，force_band (0.005m) 用于实际施力
3. **force_band 收窄**：从 0.3m 缩小到 0.005m，确保只有真正接近表面的点才施力
4. **参数扫描**：45 个配置覆盖不同参数组合，找到稳健窗口

---

## 5. 推荐参数窗口

### 最佳点
```
stiffness       = 1e5 N/m
damping         = 500 N*s/m
force_band      = 0.005 m
activation_band = 0.1 m (固定)
sample_count    = 128
```

### 稳健窗口
```
stiffness       ∈ [5e4, 2e5] N/m
damping         ∈ [300, 600] N*s/m
force_band      ∈ [0.0, 0.01] m
activation_band = 0.1 m (建议固定)
```

### 窗口内保证
- y_error < 0.02m
- 不会提前托住 (y < 0.25)
- 不会深度穿透 (y > 0.15)
- 仿真稳定 (不爆炸)

---

## 6. 进入 patch primitive 的适用性评估

### 当前优势
1. **几何可信**：y_error < 2.5% 特征长度
2. **物理正确**：stiffness 主导的接触力，仅穿透后施力
3. **参数稳健**：有明确的推荐参数窗口

### 需要注意
1. **力震荡**：平衡点附近接触力有较大波动（100~800N），需要进一步阻尼或滤波
2. **稀疏接触**：最佳配置下仅 ~2 个 force_sample，可能导致力分布不连续
3. **场景单一**：仅验证了球-平面接触，未验证边缘/多面接触

### 结论
**可以继续进入 patch primitive 阶段**，但建议在 patch 阶段同时解决力震荡和接触稀疏问题。

---

## 7. 验收结论

### 里程碑 2B.5 状态：**PASSED** ✅

所有 8 项验收标准均已满足。点式接触模型已从"能跑"提升到"更可信"，为后续 patch primitive 阶段奠定了可靠基础。
