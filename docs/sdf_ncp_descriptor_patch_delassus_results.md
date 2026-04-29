# SDF-NCP descriptor patch Delassus block solve 结果记录

本记录说明本轮将 `local patch Delassus pressure solve` 从显式 force backend 推进到 Chrono
descriptor/block constraint 后端的实现和结果。

## 1. 实现目标

目标不是继续通过 `IntLoadResidual_F` 显式施加 patch pressure force，而是让 patch pressure、body
velocity、joint constraints 和 driver constraints 在同一个 Chrono descriptor 速度层同步迭代。

新的 descriptor block 形式为：

```text
W = J M^{-1} J^T + CFM + alpha L

0 <= p_i  perpendicular  r_i + sum_j W_ij (p_j - p_j_old) >= 0
```

其中：

1. `p_i` 是 sample pressure/local intensity；
2. `J_i` 已包含真实面积尺度；
3. 真实接触力为 `force_i = p_i area_i n_i`；
4. 力矩为 `torque_i = (x_i - x_body) cross force_i`；
5. `W = J M^{-1} J^T` 在 Chrono PSOR group projection 内由 descriptor variables 扰动得到；
6. `L` 是可选 patch graph Laplacian pressure regularization。

## 2. 修改位置

核心后端：

```text
src/chrono/physics/ChSdfNcpConstraintContact.h
```

新增/扩展内容：

1. `ChSdfNcpPatchProjectedConstraint` 的 patch block solve 可加入 Laplacian compliance；
2. `ChSdfNcpConstraintContactSet::SetPatchBlockLaplacianCompliance(...)`；
3. PSOR group projection 中同一 patch 的 sample rows 解局部 Delassus block，并一次性更新 descriptor state。

rev_joint_clearance 入口：

```text
src/demos/core/demo_CH_sdf_ncp_benchmarks_openvdb.cpp
```

新增命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_delassus
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_delassus_laplacian
```

## 3. 输出文件

每个 case 输出：

```text
results/sdf_ncp_benchmarks/<case>/trajectory.csv
results/sdf_ncp_benchmarks/<case>/summary.csv
results/sdf_ncp_benchmarks/<case>/comparison.csv
results/sdf_ncp_benchmarks/<case>/contact_patch_moment_summary.csv
results/sdf_ncp_benchmarks/<case>/wrench_alignment.csv
results/sdf_ncp_benchmarks/<case>/wrench_alignment_summary.csv
```

统一 wrench 对比脚本：

```text
python scripts\compare_rev_joint_clearance_wrench.py
```

输出：

```text
results/sdf_ncp_benchmarks/rev_joint_clearance_wrench_comparison/wrench_summary.csv
results/sdf_ncp_benchmarks/figures/rev_joint_clearance_wrench_comparison/*.png
```

## 4. 三秒 rev_joint_clearance 对比摘要

| case | max penetration | force RMSE vs ideal wrench | torque RMSE vs ideal wrench | Body2 y RMSE | Body2 z RMSE | open positive rows | success |
|---|---:|---:|---:|---:|---:|---:|---:|
| default descriptor SDF-NCP | `1.092120073735714e-04` | `1.4069907230224984e+05` | `5.485076845841062e+04` | `0.664689965769521` | `0.4154834931757474` | `715` | `1` |
| descriptor patch Delassus block | `1.3728335034102201e-04` | `1.1563282324972148e+05` | `4.875358138845546e+04` | `0.6699946583657322` | `0.43234136533413503` | `695` | `1` |
| descriptor patch Delassus block + Laplacian | `1.3514823513105512e-04` | `1.1390754824133523e+05` | `4.90253659845047e+04` | `0.6674019363056807` | `0.43154870547073626` | `591` | `1` |
| local Delassus force backend | `2.457255031913519e-04` | `1.1476310345308966e+05` | `1.9787467864230468e+04` | `0.5916130397757652` | `0.43994188983870625` | `57` | `1` |

## 5. 结论

1. descriptor-level patch Delassus block 已经进入 Chrono 全局速度约束求解路径，结构上满足本阶段目标。
2. 相比默认 descriptor SDF-NCP，它降低了 ideal wrench force RMSE，并小幅降低 torque RMSE。
3. Laplacian regularization 进一步降低 open positive rows 和 force RMSE，但没有明显改善 torque RMSE。
4. local Delassus force backend 的 torque RMSE 仍更好，但它不是最终结构，因为它没有与 joint/driver constraints 在 descriptor 中同步迭代。
5. 下一步不应回退到显式 force backend，而应继续改 descriptor residual/RHS 的尺度一致性、active-set stabilization 和 patch pressure scaling。

## 6. 验证命令

```text
cmake --build build --config Release --target demo_CH_sdf_ncp_benchmarks_openvdb -- /m:4
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_delassus
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_delassus_laplacian
python scripts\compare_rev_joint_clearance_wrench.py
ctest --test-dir build -C Release -L sdf_ncp_benchmark --output-on-failure
ctest --test-dir build -C Release -L sdf_ncp --output-on-failure
ctest --test-dir build -C Release -L field_contact -E "paper_examples" --output-on-failure
pytest
```
