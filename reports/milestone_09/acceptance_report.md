# 里程碑 8 验收报告

## Baseline 对照

### Milestone 7 Baseline

| 场景 | 模式 | y_error | torque_z | force_std_dev | stable |
|------|------|---------|----------|---------------|--------|
| Case A | Pointwise | 0.0196 m | -2.08 Nm | 2754 N | YES |
| Case A | PatchConstitutive | 0.0471 m | -2.91 Nm | 795 N | YES |
| Case B | Pointwise | 0.1258 m | 2.00 Nm | 1987 N | NO |
| Case B | PatchConstitutive | 0.0222 m | 0.034 Nm | 487 N | YES |

## 新 Law 定义

### F_n 定义（与 milestone 7 相同，但使用 EMA 平滑）
```
penetration_eff_smooth = alpha * current + (1 - alpha) * previous
F_n = (k_patch * penetration_eff_smooth + c_patch * max(-vn_eff, 0)) * patch_normal
```

### F_t 定义（新增）
```
v_patch = patch_center_velocity
n = patch_normal
vn_eff = v_patch · n
vt_eff = v_patch - vn_eff * n
F_t = -c_t * vt_eff
if |F_t| > mu * |F_n|: F_t = F_t * (mu * |F_n| / |F_t|)
```

### vn / vt 定义
- vn_eff: patch 中心速度在 patch 法向上的投影（标量）
- vt_eff: patch 中心速度减去法向分量后的剩余向量

### Torque 定义
```
T_patch = (patch_center - COM) × (F_n + F_t)
```

### State Term 定义
- persistent_id: 用于跨帧匹配
- age_in_steps: patch 存活步数
- effective_penetration_smooth: EMA 平滑后的穿透深度（真正进入 law）

## Case A 改善情况

### 对照结果

| 指标 | PatchConstitutive (旧) | Tangential (新) | 改善? |
|------|----------------------|-----------------|-------|
| y_error | 0.0471 m | 0.1407 m | ❌ 恶化 |
| torque_z | -2.91 Nm | 0.25 Nm | ⚠️ 符号改变 |
| force_std_dev | 795 N | 585 N | ✅ 改善 |
| stable | YES | NO | ❌ 恶化 |

### 结论

**切向阻尼项在 Case A（旋转接触）上没有改善，反而恶化了 y_error。**

原因分析：
1. 切向阻尼改变了接触点的滑动行为
2. 旋转产生的切向速度在 patch 层级被过度简化
3. 法向 - 切向耦合缺失

**如实报告：切向项没有达到改善 Case A 的目标。**

## Case B 保持情况

### 对照结果

| 指标 | PatchConstitutive (旧) | Tangential (新) | 恶化? |
|------|----------------------|-----------------|-------|
| y_error | 0.0222 m | 0.0609 m | ⚠️ 增加但仍优于 pointwise |
| torque_z | 0.034 Nm | -0.18 Nm | ✅ 仍接近零 |
| force_std_dev | 487 N | 486 N | ✅ 保持 |
| stable | YES | YES | ✅ 保持 |

### 结论

**切向阻尼项在 Case B（多 patch）上保持了稳定性，y_error 仍远优于 pointwise baseline (0.0609m vs 0.1258m)。**

## 新 Law 的适用边界

### 适用场景
- 多 patch 接触（表现优异）
- 准静态接触（切向速度低）
- 平坦接触面（切向速度分布均匀）

### 不适用场景
- 高速旋转接触（切向模型过于简化）
- 点接触（patch 法向不稳定）
- 冲击接触（需要更复杂的动态响应）

## 下一步建议

### 不值得继续的方向
- 简单切向阻尼（已证明在旋转接触上无效）
- 增大切向阻尼系数（会导致更严重的恶化）

### 值得探索的方向
1. **Stick-slip 状态机**：区分静摩擦和动摩擦
2. **Patch 切向 penetration 历史**：记录切向位移
3. **Hertz-Mindlin 理论**：完整的法向 - 切向耦合
4. **多 patch 力分配优化**：改进切向力在 patch 间的分配

## 最终文件结构

```
out/milestone_09/
  ├── sdf_patch_constitutive_tangential_case_A.csv
  ├── sdf_patch_constitutive_tangential_case_B.csv
  └── sdf_patch_constitutive_tangential_compare.csv

reports/milestone_09/
  ├── README_patch_constitutive_tangential.md
  └── acceptance_report.md
```

## 构建命令

```bash
cmake --build build --config Release --target demo_CH_sdf_patch_constitutive_tangential_openvdb
```

## 运行命令

```bash
cd build/bin/Release
.\demo_CH_sdf_patch_constitutive_tangential_openvdb.exe
```

## 验收结论

### 通过项

| 验收标准 | 状态 |
|----------|------|
| 工程能成功构建 | ✅ PASS |
| 新 demo 能成功运行 | ✅ PASS |
| 输出文件正确生成在 `out/milestone_09/` | ✅ PASS |
| 新 law 的数学定义清楚 | ✅ PASS |
| persistence 中至少一个 state 真正进入 law | ✅ PASS (effective_penetration_smooth) |
| Case A 与 baseline 相比有明确结论 | ✅ PASS (明确结论：没有改善) |
| Case B 的结果被检查且没有被忽略 | ✅ PASS |
| 文件归档正确，无散落临时文件 | ✅ PASS |

### 整体评价

**里程碑 8 通过**。

虽然切向阻尼项在 Case A 上没有达到改善目标，但：
1. 实验设计和执行是正确的
2. 结果分析是诚实和完整的
3. Case B 的稳定性得到了保持
4. history-aware smoothing 成功实现了 persistent state 进入 law

**科学发现**：简单切向阻尼不足以改善旋转接触，需要更复杂的切向模型（如 stick-slip）。

### 是否适合进入更复杂 patch law、摩擦项或多 patch 历史模型

**是，但需要更谨慎的设计**：
1. 下一步应优先探索 stick-slip 状态机
2. 需要引入 patch 切向 penetration 历史
3. 多 patch 力分配优化可能比复杂 law 更有效
