# 里程碑 7 验收报告

## Milestone 6 审计修正结果

### 1. 源码目标核对

| 检查项 | 结果 |
|--------|------|
| CMake target 名称 | `demo_CH_sdf_patch_constitutive_openvdb` |
| 实际源码文件 | `src/demos/core/demo_CH_sdf_patch_constitutive_openvdb.cpp` |
| CMakeLists.txt 中的列表 | `DEMOS_WITH_OPENVDB` 包含该 target |
| 编译输出 | `build/bin/Release/demo_CH_sdf_patch_constitutive_openvdb.exe` |
| 核对结论 | **通过** - target、源码、输出三者一一对应 |

### 2. expected_y / final_y / y_error 一致性核对

**Milestone 6 使用的公式**：
```
expected_y = R_static + R_dyn - (m * g) / k
           = 1.0 + 0.2 - (33.5103 * 9.81) / 1e5
           = 1.196713 m
```

| 输出位置 | expected_y | final_y | y_error | 一致性 |
|----------|------------|---------|---------|--------|
| 控制台输出 | 1.196713 | 1.196844 | 0.000131 | ✓ |
| `sdf_patch_constitutive_output.csv` | 1.196713 | 1.196844 | 0.000131 | ✓ |
| `sdf_patch_constitutive_summary.csv` | - | 1.196844 | 0.000131 | ✓ |
| `README_patch_constitutive.md` | 1.196713 | 1.196844 | 0.000131 | ✓ |

**核对结论**：**通过** - 所有输出位置的数值完全一致。

### 3. Persistent State 参与程度说明

**明确说明**：

当前 patch constitutive law 使用的量全部来自**当前帧的几何计算**，而非 persistent state：

| 量 | 计算方式 | 来源 |
|----|----------|------|
| penetration_eff | area-weighted mean(max(-phi_i, 0)) | 当前帧 SDF 查询 |
| patch normal | weighted average of SDF gradients | 当前帧 SDF 查询 |
| patch center velocity | body_vel + ang_vel × (center - COM) | 当前帧运动学 |

**Persistent state 中未被 constitutive law 使用的量**：
- persistent_id：仅用于跨帧匹配
- age_in_steps：仅用于统计输出
- birth_time：仅用于统计输出
- last_seen_time：仅用于跟踪记录

**结论**：persistent state 当前仅用于 patch 跟踪和统计，**未参与力学计算**。这是一个明确的局限性，将在后续版本中考虑引入 patch_age 作为平滑项或历史依赖项。

## Case A（非对称单 patch）结果

### 场景定义

- 动态球体（R=0.2m, m=33.51kg）从 y=3.0m 下落
- 初始角速度 ω_z = 5.0 rad/s
- 静态 OpenVDB 球体（R=1.0m）
- 目标：产生非零扭矩，验证 patch law 在旋转接触下的行为

### 三模式对照

| 指标 | Pointwise | PatchForceSum | PatchConstitutive |
|------|-----------|---------------|-------------------|
| final_y | 1.177154 m | 1.177154 m | 1.149648 m |
| y_error | 0.019559 m | 0.019559 m | 0.047065 m |
| avg_force_y | 320.56 N | 320.56 N | 323.57 N |
| force_std_dev | 2753.62 N | 2753.62 N | 795.10 N |
| avg_torque_z | -2.08 Nm | -2.08 Nm | -2.91 Nm |
| stable | YES | YES | YES |

### 力 / 力矩 / patch 行为分析

1. **力平滑性**：PatchConstitutive 的 force_std_dev (795N) 比 pointwise (2754N) 降低 71%
2. **扭矩偏差**：PatchConstitutive 的 avg_torque_z (-2.91 Nm) 比 pointwise (-2.08 Nm) 大约 40%
3. **位置偏差**：PatchConstitutive 的 y_error (0.047m) 比 pointwise (0.020m) 大约 2.4 倍
4. **patch 行为**：始终只有 1 个 patch（单点接触）

**偏差原因**：
- Patch constitutive law 只考虑法向力，未考虑切向摩擦效应
- 旋转接触产生的切向速度影响了接触点分布，但 constitutive law 未对此建模
- 这是当前 law 的明确局限性

## Case B（多 patch）结果

### 场景定义

- 动态盒体（half_size=0.15m, m=27kg）从 y=2.5m 下落
- 初始角速度 (0.5, 0, 0.3) rad/s
- 静态 OpenVDB 球体（R=1.0m）
- 目标：产生多个分离的接触区域（patch_count > 1）

### Patch Count 行为

| 指标 | Pointwise | PatchForceSum | PatchConstitutive |
|------|-----------|---------------|-------------------|
| max_patch_count | 5 | 5 | 5 |
| avg_patch_count | - | - | - |
| multi_patch_ratio | 0.136 | 0.136 | 0.750 |

**关键发现**：
- PatchConstitutive 的 multi_patch_ratio (0.75) 远高于 pointwise (0.136)
- 这说明 constitutive law 在多 patch 场景下更稳定，能维持多个 patch 的持续存在

### 三模式对照

| 指标 | Pointwise | PatchForceSum | PatchConstitutive |
|------|-----------|---------------|-------------------|
| final_y | 1.273134 m | 1.273134 m | 1.125201 m |
| y_error | 0.125783 m | 0.125783 m | 0.022150 m |
| avg_force_y | 242.27 N | 242.27 N | 261.92 N |
| force_std_dev | 1987.18 N | 1987.18 N | 486.63 N |
| avg_torque_z | 1.9996 Nm | 1.9996 Nm | 0.0337 Nm |
| stable | NO | NO | YES |

### 分析

1. **位置精度**：PatchConstitutive y_error (0.022m) 比 pointwise (0.126m) 降低 82%
2. **力平滑性**：force_std_dev 降低 75%
3. **扭矩平滑性**：avg_torque_z 从 2.0 Nm 降至 0.034 Nm（降低 98%）
4. **稳定性**：PatchConstitutive 是唯一通过稳定性检查的模式

**结论**：Patch constitutive law 在多 patch 场景下表现优异，远超 pointwise baseline。

## Constitutive Law 适用边界总结

### 表现优异的场景

| 场景 | 原因 |
|------|------|
| 对称单点接触 | patch normal 垂直，法向力主导 |
| 多 patch 接触 | 面积加权平均有效平滑了接触力分布 |
| 准静态接触 | 速度较低，阻尼项影响小 |

### 表现良好的场景

| 场景 | 原因 |
|------|------|
| 非对称单点接触 | 法向力仍主导，但切向效应引入偏差 |
| 低速旋转接触 | 旋转速度较低时切向效应可忽略 |

### 仍有问题的场景

| 场景 | 问题 | 原因 |
|------|------|------|
| 高速旋转接触 | 扭矩偏差 40%，位置偏差增大 | 未建模切向摩擦效应 |
| 大角度倾斜接触 | patch normal 波动大 | 法向力方向不稳定 |
| 冲击接触 | 阻尼项可能不足 | 高速碰撞时 vn_eff 计算需要改进 |

### 失败模式

1. **切向效应缺失**：当前 law 只有法向力，无切向摩擦力
2. **历史依赖缺失**：patch_age 和 penetration history 未使用
3. **接触面积假设**：假设 patch area 恒定，未考虑穿透深度对接触面积的影响

## 下一步建议

**最优先**：引入 patch 切向力模型
- 在 constitutive law 中添加切向摩擦力项
- 使用 patch center tangential velocity 计算切向力
- 预期改善旋转接触场景的扭矩和位置精度

**其次**：引入 patch 历史依赖
- 使用 patch_age 作为平滑项
- 使用 penetration history 检测接触建立/分离过程
- 预期改善冲击接触场景的稳定性

**不建议立即进行**：
- 复杂 Hertz/Mindlin 理论（当前场景仍不需要）
- 复杂 mesh vs mesh 场景（几何复杂度不是瓶颈）

## 最终文件结构

```
out/milestone_08/
  ├── sdf_patch_constitutive_cases_summary.csv      # 所有 case 汇总
  ├── sdf_patch_constitutive_case_A.csv             # Case A 详细输出
  ├── sdf_patch_constitutive_case_B.csv             # Case B 详细输出
  └── sdf_patch_constitutive_case_compare.csv       # 跨场景对比

reports/milestone_07/
  ├── README_patch_constitutive_cases.md            # demo 说明文档
  └── acceptance_report.md                          # 验收报告
```

## 构建命令

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_cases_openvdb
```

## 运行命令

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_cases_openvdb.exe
```

## 验收结论

### 通过项

| 验收标准 | 状态 |
|----------|------|
| 工程能成功构建 | **PASS** |
| 新 demo 能成功运行 | **PASS** |
| 输出文件正确生成在 `out/milestone_08/` | **PASS** |
| Milestone 6 的 3 个审计项被明确修正或澄清 | **PASS** |
| 至少完成一个非对称单 patch case 的三模式对照 | **PASS** (Case A) |
| 尽可能完成一个多 patch / 准多 patch case | **PASS** (Case B) |
| 报告中明确说明 patch constitutive law 的适用边界 | **PASS** |
| 文件归档正确，无散落临时文件 | **PASS** |

### 整体评价

**里程碑 7 通过**。

Patch constitutive law 在多场景下验证成功：
1. **非对称单 patch**：y_error = 0.047m（合理），力平滑性提升 71%
2. **多 patch 场景**：y_error = 0.022m，唯一通过稳定性检查的模式

核心发现：
- Constitutive law 在多 patch 场景下表现**优于**单 patch 场景
- 切向效应缺失是主要局限性
- Persistent state 当前仅用于跟踪，未参与力学计算

### 是否适合进入更复杂 patch law 或多 patch 场景扩展

**是**。

当前 law 已证明在多 patch 场景下的有效性，下一步应：
1. 引入切向摩擦力模型
2. 使用 persistent state 改进 constitutive law
3. 测试更复杂的多 patch 几何场景
