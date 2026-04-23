# 里程碑 6 验收报告

## Patch Constitutive Law 数学定义

### 核心公式

```
penetration_eff = area-weighted mean of max(-phi_i, 0)
vn_eff = patch_center_velocity · patch_normal
F_patch = (k_patch * penetration_eff + c_patch * max(-vn_eff, 0)) * patch_normal
T_patch = (patch_center - COM) × F_patch
```

### 详细说明

- **penetration_eff**: 对 patch 内所有 active samples 的穿透深度 `max(-phi_i, 0)` 进行面积加权平均
- **vn_eff**: patch 中心点在 patch 法向上的投影速度（相对速度法向分量）
- **F_patch**: 法向力大小由刚度项和阻尼项组成，方向沿 patch_normal
- **T_patch**: 力矩由 patch 中心到质心的力臂与 patch 力叉乘得到

### 与 pointwise 的区别

| 特征 | Pointwise | Patch Constitutive |
|------|-----------|-------------------|
| 力计算单元 | 单个 sample | 整个 patch |
| 输入量 | sample phi, sample velocity | patch mean penetration, patch center velocity |
| 力方向 | sample grad 方向 | patch normal 方向 |
| 聚合方式 | 直接求和 | 单一法向力 |

## 参数设置

| 参数 | 值 | 单位 | 说明 |
|------|-----|------|------|
| 静态球半径 | 1.0 | m | OpenVDB level set sphere |
| 动态球半径 | 0.2 | m | Chrono ChBodyEasySphere |
| 动态球质量 | 33.5103 | kg | 密度 1000 kg/m³ |
| 下落高度 | 3.0 | m | 初始 y 位置 |
| 网格体素 | 0.05 | m | OpenVDB voxel size |
| 采样分辨率 | 8×16 | - | n_theta × n_phi |
| 刚度 k_patch | 1e5 | N/m | 与 pointwise 一致 |
| 阻尼 c_patch | 500 | Ns/m | 与 pointwise 一致 |
| 激活带宽 | 0.1 | m | candidate detection |
| 力计算带宽 | 0.0 | m | actual penalty |
| 时间步长 | 5e-4 | s | simulation step |
| 总模拟时间 | 2.0 | s | total duration |

## 与 Baselines 的对照

### 关键指标对比

| 指标 | Pointwise | PatchForceSum | PatchConstitutive | 期望值 |
|------|-----------|---------------|-------------------|--------|
| final_y | 1.200327 m | 1.200327 m | 1.196844 m | 1.196713 m |
| y_error | 0.003615 m | 0.003615 m | 0.000131 m | 0 m |
| avg_force_y | 309.25 N | 309.25 N | 325.81 N | 328.73 N (mg) |
| force_std_dev | 2868.74 N | 2868.74 N | 833.79 N | - |
| stable | NO | NO | YES | YES |

### 分析

#### 位置精度
- **Patch Constitutive y_error = 0.000131m**，相比 baseline 降低 **96.4%**
- 这是所有里程碑中最好的位置精度
- Pointwise 和 PatchForceSum 完全一致（数学等价性验证通过）

#### 力精度
- **Patch Constitutive avg_force = 325.81N**，与 mg=328.73N 相差仅 0.9%
- Baseline avg_force = 309.25N，与 mg 相差 5.9%
- Patch constitutive 的力平衡更精确

#### 稳定性
- **Patch Constitutive 是唯一通过稳定性检查的模式**
- Force std dev 降低 71%（833.79N vs 2868.74N）
- 说明 patch-level 聚合有效平滑了数值振荡

#### 为什么 PatchForceSum 与 Pointwise 完全一致？
- 因为 PatchForceSum 模式下 patch_force = Σ sample_force_i
- 根据力的叠加原理，两者数学等价
- 这验证了 patch grouping 逻辑的正确性

#### 为什么 Patch Constitutive 更优？
1. **单一法向方向**: 消除了 sample grad 方向不一致引起的横向振荡
2. **面积加权平均**: 穿透深度聚合减少了离散采样噪声
3. **中心速度投影**: 使用 patch 整体运动而非单个 sample 速度，更稳定

## 当前局限

1. **单 patch 场景**: 球-球接触始终只产生 1 个 patch，多 patch 交互尚未测试
2. **线性 law**: 当前仅实现最简单的线性刚度+阻尼，未包含非线性效应
3. **无历史项**: 未使用 patch age、历史穿透等时间相关信息
4. **未使用 persistence 状态**: 虽然实现了匹配，但 constitutive law 未利用 persistent track 信息
5. **Hertz 理论缺失**: 未引入接触面积与穿透深度的非线性关系

## 验收结论

### 通过项

| 验收标准 | 状态 |
|----------|------|
| 工程能成功构建 | PASS |
| demo 能成功运行 | PASS |
| 输出文件正确生成 | PASS |
| 成功输出 patch constitutive mode 数据 | PASS |
| 完成 constitutive vs baselines 对照 | PASS |
| Constitutive mode 不崩溃、不发散 | PASS |
| 报告清楚说明 patch law 定义与结果 | PASS |
| 文件归档正确，无散落临时文件 | PASS |

### 整体评价

**里程碑 6 通过**。

Patch constitutive law 首版成功实现了从 patch 几何量到接触力的直接映射。在简单球-球接触场景下：
- 位置精度超越 baseline 96%
- 力振荡降低 71%
- 是唯一通过稳定性检查的模式

这表明 patch-level constitutive law 是一个正确的方向。

### 是否适合进入更复杂 patch law 或多 patch 场景

**是，但有条件**：

**适合进入多 patch 场景的理由**：
1. 单 patch 场景已充分验证 constitutive law 的正确性
2. 多 patch 场景能测试 patch grouping 和 force 分配逻辑
3. 当前框架已支持多 patch（BFS 连通分量）

**建议先完成的工作**：
1. 在 constitutive law 中引入 patch persistence 状态（age、history）
2. 测试简单多 patch 场景（如方盒落在粗糙表面上）
3. 验证 patch 间力分配是否合理

**不建议立即进行的工作**：
1. 复杂非线性 constitutive law（Hertz、Mindlin）
2. 复杂 mesh vs mesh 场景
3. 粘滑摩擦模型
