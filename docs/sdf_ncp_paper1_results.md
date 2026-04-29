# SDF-NCP 第一篇论文数值结果说明

本文档整理当前仓库中已经生成的第一篇传统 SDF-NCP 论文数值结果。目标是为论文写作提供图表编号、表格来源、结果解读和复现实验命令，而不是继续扩展算法功能。

## 1. 论文目标

第一篇论文的目标是建立一个非 AI、传统计算力学意义下的 Signed Distance Field + Nonlinear Complementarity Problem 接触动力学基础方法。当前工作强调 SDF 提供几何 gap 和 normal，NCP 提供 Signorini 单边接触约束的等式化残差。

核心公式为：

```text
g_i(q) = phi(x_i(q))
n_i = grad phi(x_i(q))
J_i(q) = grad phi(x_i(q))^T dx_i(q)/dq

0 <= g_i(q) ⟂ lambda_i >= 0

Phi_eps(g_i, lambda_i) =
    sqrt(g_i^2 + lambda_i^2 + eps^2) - g_i - lambda_i
```

论文应明确说明：本阶段没有使用神经网络、机器学习或 AI 修正项；SDF-NCP 是一个经典接触动力学离散与求解框架。

## 2. 当前实验覆盖范围

当前数值结果覆盖以下内容：

1. SDF 几何可视化：SDF 等值线、零水平集和梯度法向。
2. 点质量-平面接触：penalty 与 SDF-NCP 在同一落体问题上的时间响应。
3. Penalty stiffness sensitivity：不同 `k_n` 下的最大穿透和互补误差。
4. Smooth FB epsilon sensitivity：不同 `eps` 下的穿透、互补误差、NCP residual 和迭代统计。
5. C++ 2D 刚体双底点 rollout：包含姿态、角速度、双接触点 gap、lambda 和互补误差。
6. C++/Python 回归测试：覆盖 SDF、NCP、点质量、多点 contact set 和 2D 刚体局部点 residual。

当前不覆盖：

1. 3D 刚体接触。
2. 摩擦或摩擦锥互补。
3. 柔性体接触。
4. Chrono 主 solver descriptor 或 production contact container 接入。
5. AI/ML、神经网络、参数学习或数据驱动 SDF 修正。

## 3. 推荐论文图表编号

### Fig. 1: SDF 几何层示意

来源：

```text
results/sdf_ncp/geometry/sdf_contours_normals.png
```

用于说明：SDF 零水平集、等值线以及由 `grad phi` 给出的法向方向。

论文中要表达的结论：SDF 将接触几何转化为连续标量场与梯度场，使 gap 与 normal 可以由同一个几何函数获得。

### Fig. 2: 点质量落地接触响应

来源：

```text
results/sdf_ncp/point_mass_plane/height_vs_time.png
results/sdf_ncp/point_mass_plane/gap_vs_time.png
results/sdf_ncp/point_mass_plane/contact_force_vs_time.png
results/sdf_ncp/point_mass_plane/complementarity_error_vs_time.png
```

可选补充图：

```text
results/sdf_ncp/point_mass_plane/penetration_vs_time.png
```

用于说明：SDF-NCP 可处理基本接触激活、闭合接触和接触力响应过程。

论文中要表达的结论：对平面 SDF，gap 直接等于点质量高度；NCP residual 将非穿透与法向乘子互补关系加入时间离散求解。

### Fig. 3: penalty stiffness sensitivity

来源：

```text
results/sdf_ncp/penalty_sensitivity/max_penetration_vs_parameter.png
```

用于说明：penalty 方法对 `k_n` 的选择敏感。

论文中要表达的结论：在当前点质量落地算例中，penalty 最大穿透随接触刚度变化明显；这支持将 penalty 作为基线，而不是作为最终非穿透约束。

### Fig. 4: smooth FB epsilon sensitivity

来源：

```text
results/sdf_ncp/epsilon_sensitivity/max_penetration_vs_eps.png
results/sdf_ncp/epsilon_sensitivity/complementarity_vs_eps.png
results/sdf_ncp/epsilon_sensitivity/iterations_vs_eps.png
```

用于说明：`eps` 控制互补约束逼近与数值求解难度之间的平衡。

论文中要表达的结论：较小 `eps` 通常对应更严格的平滑 FB 逼近，但数值代价应通过 `success_rate` 和平均迭代数一起报告，避免只看误差曲线。

### Fig. 5: 2D rigid body rollout

来源：

```text
results/sdf_ncp_cpp/figures/rigidbody2d_pose_vs_time.png
results/sdf_ncp_cpp/figures/rigidbody2d_gaps_vs_time.png
results/sdf_ncp_cpp/figures/rigidbody2d_lambdas_vs_time.png
results/sdf_ncp_cpp/figures/rigidbody2d_complementarity_vs_time.png
```

用于说明：SDF-NCP 框架可扩展到包含姿态、局部点运动学和广义接触力矩的 2D 刚体局部点接触。

论文中要表达的结论：2D 刚体算例验证了

```text
J_phi(q) = grad phi(x_i(q))^T dx_i(q)/dq
```

并展示多点 assembly 可以同时处理两个底部局部接触点。

## 4. 推荐论文表格

### Table 1: penalty vs NCP 参数敏感性摘要

来源：

```text
results/sdf_ncp_paper1/tables/method_summary.csv
```

列含义：

| 列名 | 含义 |
| --- | --- |
| `experiment` | 数据来源实验 |
| `method` | penalty、SDF-NCP 或 C++ 2D rigid body SDF-NCP |
| `parameter` | `k_n`、`eps` 或固定实验参数 |
| `max_penetration` | 该实验配置下的最大穿透 |
| `max_complementarity_error` | 最大互补误差 |
| `mean_iterations` | 平均非线性求解迭代数；penalty 为 0 |
| `success_rate` | 步进或求解成功比例 |

用途：汇总 penalty stiffness 和 NCP epsilon 敏感性，并把 C++ 2D 刚体 rollout 结果放在同一张表中。

### Table 2: C++/Python 回归测试覆盖

建议论文中列出：

| 测试命令 | 当前结果 | 用途 |
| --- | --- | --- |
| `pytest` | 11/11 passed | Python SDF、NCP、penalty 和点质量原型 |
| `ctest --test-dir build -C Release -L sdf_ncp --output-on-failure` | 16/16 passed | C++ SDF-NCP residual、contact set 和 2D 刚体局部点 |
| `ctest --test-dir build -C Release -L field_contact --output-on-failure` | 30/30 passed | 确认新增 SDF-NCP 层未破坏已有 field_contact 回归 |

用途：支撑方法实现的可复现性，并说明新增研究代码没有破坏已有 field-contact 测试。

## 5. 每个实验回答的论文问题

| 实验 | 文件/脚本 | 结果路径 | 回答的问题 | 论文中建议位置 |
| --- | --- | --- | --- | --- |
| SDF geometry visualization | `examples/sdf_ncp/sdf_geometry_visualization.py` | `results/sdf_ncp/geometry/sdf_contours_normals.png` | SDF gap、零水平集和 normal 是否直观一致？ | SDF-based Contact Geometry |
| Point-mass plane response | `examples/sdf_ncp/point_mass_plane.py` | `results/sdf_ncp/point_mass_plane/` | SDF-NCP 是否能处理基本落地接触响应？ | Numerical Experiments |
| Penalty sensitivity | `examples/sdf_ncp/penalty_sensitivity.py` | `results/sdf_ncp/penalty_sensitivity/summary.csv` | penalty 基线对 `k_n` 是否敏感？ | Numerical Experiments |
| Epsilon sensitivity | `examples/sdf_ncp/epsilon_sensitivity.py` | `results/sdf_ncp/epsilon_sensitivity/summary.csv` | `eps` 对互补误差、穿透和迭代数有何影响？ | Numerical Experiments / Discussion |
| C++ rigidbody2d export | `build/bin/Release/demo_CH_sdf_ncp_regression.exe rigidbody2d_export` | `results/sdf_ncp_cpp/rigidbody2d_rollout.csv` | C++ 2D 刚体局部点 residual 是否能稳定 rollout？ | Numerical Experiments |
| C++ rigidbody2d plotting | `examples/sdf_ncp/plot_cpp_rigidbody2d_rollout.py` | `results/sdf_ncp_cpp/figures/` | 2D 刚体姿态、gap、lambda 和 residual 如何随时间变化？ | Numerical Experiments |
| Paper1 reproduction | `scripts/run_sdf_ncp_paper1_experiments.py` | `results/sdf_ncp_paper1/tables/method_summary.csv` | 是否可以一键复现实验数据和论文图表？ | Reproducibility paragraph / Appendix |

## 6. 关键数值摘要

以下统计来自：

```text
results/sdf_ncp_paper1/tables/method_summary.csv
```

当前 summary 表包含 15 条记录，其中：

- `penalty_sensitivity`: 8 条；
- `epsilon_sensitivity`: 6 条；
- `rigidbody2d_rollout`: 1 条。

Penalty stiffness sensitivity 中，penalty 方法的 `max_penetration` 范围为：

```text
min = 0.004443929999999979
max = 0.5517827999146252
```

Epsilon sensitivity 中，SDF-NCP 的 `success_rate` 范围为：

```text
min = 0.9920079920079921
max = 1.0
```

Epsilon sensitivity 中，平均迭代数范围为：

```text
min = 11.606393606393606
max = 15.177822177822177
```

Epsilon sensitivity 中，最大互补误差范围为：

```text
min = 5.44229736422673e-10
max = 5.000190615342035e-05
```

所有 summary 记录的 `success_rate` 范围为：

```text
min = 0.9920079920079921
max = 1.0
```

因此，并非所有记录的 `success_rate` 都等于 1.0。论文中应客观报告成功率范围，不应写成所有参数均完全收敛。

C++ 2D rigid-body rollout 的 summary 记录为：

```text
max_penetration = 1.8270662760500045e-12
max_complementarity_error = 1.2867972801741118e-09
mean_iterations = 1.2547452547452547
success_rate = 1.0
```

该 CSV 含 1001 个时间步记录。上述数值仅对当前参数、当前局部 Newton 设置和当前双底点平面接触算例成立。

## 7. 结果解读建议

### SDF 几何层

SDF 函数 `phi` 直接提供接触 gap，`grad phi` 提供接触 normal。对于解析平面和圆形 SDF，图中的零水平集和 normal 场可用于说明 SDF 将几何接触查询转化为标量场和梯度场求值。论文中可以强调该做法减少了对显式网格最近点投影的依赖，但不要声称当前实现已经解决所有复杂几何搜索问题。

### NCP 层

Signorini 条件

```text
0 <= g_i(q) ⟂ lambda_i >= 0
```

通过平滑 Fischer-Burmeister 函数写成等式 residual。与 penalty 方法相比，NCP 不是通过增大 `k_n` 来近似非穿透，而是把 gap 与法向乘子的互补关系直接放入非线性方程组。论文中应把 penalty 作为 stiffness-sensitive baseline，而不是把当前 NCP 结果描述成全局最优或无条件稳定。

### Epsilon sensitivity

`eps` 越小，平滑 FB residual 越接近非平滑互补条件；但较小 `eps` 可能改变 Newton 求解难度。当前结果应同时报告 `max_complementarity_error`、`success_rate` 和 `mean_iterations`。从 summary 表看，当前 epsilon sweep 的成功率范围为 `0.9920079920079921` 到 `1.0`，说明仍应保留 solver diagnostic，而不是只画误差曲线。

### 2D 刚体局部点

2D 刚体算例验证了局部点映射、接触 Jacobian 和广义接触力：

```text
x_i(q) = p + R(theta) r_i
J_i(q) = grad phi(x_i(q))^T dx_i(q)/dq
Q_contact = sum_i J_i^T lambda_i
```

双底点 rollout 说明当前 SDF-NCP 装配层可以处理多个局部接触点，并产生包含平动和转动力矩的广义接触响应。论文中应明确该算例仍是 2D、无摩擦、小规模局部 Newton，不是 Chrono production solver 集成。

## 8. 第一篇论文建议结构

1. **Introduction**
   - 说明 mesh contact 的接触对切换、normal 不连续、contact search 和 penalty 参数敏感问题。
   - 引出 SDF-NCP 的目标：用 SDF 定义 gap/normal，用 NCP 定义非穿透互补约束。
   - 可引用 Fig. 1 作为方法直观图。

2. **Related Work**
   - 讨论 multibody contact dynamics、Signorini 条件、NCP/LCP 方法、penalty contact 和 SDF/implicit surface contact。
   - 明确本文不是 AI/ML SDF 工作。

3. **SDF-based Contact Geometry**
   - 放入 `g_i(q) = phi(x_i(q))`、`n_i = grad phi(x_i(q))`。
   - 放入 Fig. 1。
   - 说明局部接触点到世界点映射。

4. **SDF-NCP Contact Formulation**
   - 放入 `J_i(q) = grad phi(x_i(q))^T dx_i(q)/dq`。
   - 放入 Signorini 条件与平滑 FB residual。
   - 说明点质量、多点 contact set、2D 刚体局部点的统一结构。

5. **Time Discretization and Local Newton Solver**
   - 放入 implicit Euler residual：

```text
R_v = M(v_next - v) - dt(Q + sum_i J_i^T lambda_i)
R_lambda_i = Phi_eps(g_i, lambda_i)
```

   - 说明当前使用有限差分 residual Jacobian 和局部 Newton。

6. **Numerical Experiments**
   - Fig. 2：点质量落地响应。
   - Fig. 3：penalty stiffness sensitivity。
   - Fig. 4：epsilon sensitivity。
   - Fig. 5：2D rigid body rollout。
   - Table 1：summary table。
   - Table 2：测试覆盖。

7. **Discussion**
   - 讨论 penalty 与 NCP 的差异。
   - 讨论 `eps` 与 solver cost 的 tradeoff。
   - 讨论当前 2D 刚体结果与未来 3D/摩擦/descriptor 集成之间的距离。

8. **Conclusion**
   - 总结 SDF gap/normal/Jacobian + NCP residual 的可复现原型。
   - 明确未来工作：3D 刚体、摩擦、柔性体、Chrono global solver 集成。

## 9. 当前结果不足与后续补充建议

当前不足：

1. 还没有 3D 刚体接触。
2. 还没有摩擦或摩擦锥互补。
3. 还没有柔性体接触。
4. 还没有接入 Chrono 全局 descriptor。
5. 还没有真实工程复杂几何 SDF 的系统实验。
6. 还没有解析 Hertz 接触对比。

第一篇论文必须补齐或保留的内容：

- 方法描述：SDF gap、normal、contact Jacobian、Signorini 条件和平滑 FB residual。
- Penalty vs NCP 对比。
- Epsilon sensitivity。
- 2D 刚体 Jacobian 验证和双底点 rollout。
- 可复现测试命令和生成脚本。

可以放到未来工作：

- 3D 刚体接触。
- 摩擦。
- 柔性体。
- AI-assisted SDF correction 或参数识别。
- Chrono 全局 descriptor / contact container 接入。
- 工程复杂几何和 Hertz 解析对比。

## 10. 复现实验命令

C++ demo 需要先 build：

```powershell
cmake --build build --config Release --target demo_CH_sdf_ncp_regression
```

生成第一篇论文当前所有数值结果和图表：

```powershell
python scripts\run_sdf_ncp_paper1_experiments.py
```

运行 Python 回归测试：

```powershell
pytest
```

运行 SDF-NCP C++ 回归测试：

```powershell
ctest --test-dir build -C Release -L sdf_ncp --output-on-failure
```

运行已有 field_contact 回归测试：

```powershell
ctest --test-dir build -C Release -L field_contact --output-on-failure
```

单独重新生成 summary table：

```powershell
python scripts\make_sdf_ncp_paper1_summary.py
```

如果本文档引用的 CSV 或图表缺失，优先运行：

```powershell
python scripts\run_sdf_ncp_paper1_experiments.py
```

再重新检查 `results/sdf_ncp/`、`results/sdf_ncp_cpp/` 和 `results/sdf_ncp_paper1/`。
