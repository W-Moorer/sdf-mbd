# 里程碑 7：Patch Constitutive Law 多场景验证 + 审计修正

## 概述

本里程碑将 patch constitutive law 从单一对称球-球接触场景扩展到非对称接触和多 patch 场景，同时完成 milestone 6 的 3 项审计修正。

## Milestone 6 审计修正

### 1. 源码目标核对

| 项目 | 值 |
|------|-----|
| CMake target | `demo_CH_sdf_patch_constitutive_openvdb` |
| 源码文件 | `src/demos/core/demo_CH_sdf_patch_constitutive_openvdb.cpp` |
| CMakeLists.txt | `DEMOS_WITH_OPENVDB` 列表包含该 target |
| 状态 | **已验证 - target 和源码文件一一对应** |

### 2. 数值一致性核对

| 项目 | 公式/值 |
|------|--------|
| expected_y | R_static + R_dyn - (m*g)/k = 1.0 + 0.2 - (33.51*9.81)/1e5 = **1.196713 m** |
| final_y | 模拟结束时动态球体的 y 坐标 |
| y_error | \|final_y - expected_y\| |
| 状态 | **已验证 - summary/CSV/report 三者公式一致** |

### 3. Persistent State 参与程度核对

| 量类型 | 来源 | 是否进入 constitutive law |
|--------|------|--------------------------|
| persistent_id | tracking | 否（仅用于跨帧匹配） |
| birth_time | tracking | 否 |
| age_in_steps | tracking | **否** |
| effective_penetration | 当前帧计算 | **是**（核心输入） |
| patch normal | 当前帧计算 | **是**（力的方向） |
| patch center velocity | 当前帧计算 | **是**（阻尼项输入） |

**结论**：当前 constitutive law 仅使用当前帧的几何描述符，未使用 persistent state（age/history）。tracking 仅用于 patch 跨帧 ID 匹配和统计。

## 构建方法

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_cases_openvdb
```

## 运行方法

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_cases_openvdb.exe
```

## Case 集合说明

### Case A：非对称单 patch（旋转接触）

| 参数 | 值 |
|------|-----|
| 场景 | 动态球体以初始角速度 ω_z = 5 rad/s 下落到静态 OpenVDB 球体 |
| 目的 | patch normal 不再完全垂直，patch torque 非零 |
| 动态体 | 球体 R=0.2m，质量 33.51kg |
| 初始位置 | (0, 3.0, 0) |
| 初始角速度 | (0, 0, 5.0) rad/s |

### Case B：多 patch（盒体旋转接触）

| 参数 | 值 |
|------|-----|
| 场景 | 动态盒体以初始角速度下落到静态 OpenVDB 球体 |
| 目的 | 产生多个分离的接触区域（patch_count > 1） |
| 动态体 | 盒体 half_size=0.15m，质量 27kg |
| 初始位置 | (0, 2.5, 0) |
| 初始角速度 | (0.5, 0, 0.3) rad/s |

## 三种模式说明

与 milestone 6 一致：
1. **Pointwise**: sample-by-sample force
2. **PatchForceSum**: patch_force = Σ sample_force_i
3. **PatchConstitutive**: patch force 由 constitutive law 直接生成

## 输出文件

所有输出文件位于 `out/milestone_08/`：

| 文件名 | 内容 |
|--------|------|
| `sdf_patch_constitutive_cases_summary.csv` | 所有 case 和模式的汇总对比 |
| `sdf_patch_constitutive_case_A.csv` | Case A 三模式详细输出 |
| `sdf_patch_constitutive_case_B.csv` | Case B 三模式详细输出 |
| `sdf_patch_constitutive_case_compare.csv` | Case A vs Case B 跨场景对比 |

## 关键发现

### Case A（非对称单 patch）
- Pointwise 和 PatchForceSum 完全一致（y_error = 0.0196m）
- PatchConstitutive y_error = 0.047m（合理但略大）
- 所有模式都产生了非零扭矩（约 -2 ~ -3 Nm）
- Force std dev: constitutive (795N) < pointwise (2754N)，降低 71%

### Case B（多 patch）
- PatchConstitutive y_error = 0.022m，远优于 pointwise (0.126m)
- PatchConstitutive 是唯一通过稳定性检查的模式
- Max patch count = 5，multi_patch_ratio = 0.75（constitutive 模式）
- PatchConstitutive 的平均扭矩显著低于 pointwise（0.034 vs 2.0 Nm）

## Patch Constitutive Law 适用边界

| 场景类型 | 表现 | 说明 |
|----------|------|------|
| 对称单 patch | 优异 | y_error = 0.0001m，最稳定 |
| 非对称单 patch | 良好 | y_error = 0.047m，力更平滑但位置偏差稍大 |
| 多 patch | 优异 | y_error = 0.022m，唯一稳定模式 |
| 强旋转接触 | 局限 | 未考虑切向效应，扭矩计算有偏差 |

## 结论

Patch constitutive law 在多 patch 场景下表现优异，在非对称单 patch 场景下表现良好。当前 law 的局限性在于未考虑切向效应，导致旋转接触场景下位置偏差稍大。建议在下一步中引入 patch 切向力模型。
