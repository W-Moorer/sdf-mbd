# SDF-NCP patch pressure aggregation 试验记录

本文记录 2026-04-29 对以下方案的实现和验证结果：

```text
使用成组 patch 级耦合约束或局部 patch pressure solve，而不是继续调单点参数。
```

## 已实现的试验模式

通用 SDF-NCP descriptor 后端新增了一个可选开关：

```text
ChSdfNcpConstraintContactSet::SetPatchPressureAggregation(true)
```

该模式不是针对某个算例的特化逻辑。它按 `patch_id` 将 active samples 分组，并把每个 connected active patch 聚合为一个 pressure/intensity 未知量。装配的广义接触力和力矩使用真实面积积分：

```text
F_patch = sum_i p_patch * area_i * n_i
T_patch = sum_i (x_i - x_body) cross (p_patch * area_i * n_i)
```

默认 `rev_joint_clearance` benchmark 仍保留逐 sample SDF-NCP 后端。patch aggregation 试验通过单独命令运行：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_patch_pressure
```

## 结果对比

与当前默认逐 sample pressure SDF-NCP 结果对比：

| metric | default sample pressure | patch pressure aggregation |
|---|---:|---:|
| active contact rows | `37273` | `3707` |
| open positive-force rows | `720` | `3` |
| num_contacts | `42` | `3` |
| max_penetration | `1.1383846867829561e-04` | `9.9698733538389206e-04` |
| max_effective_force_norm | `2.2890600498895184e+06` | `2.7665270590949859e+07` |
| max_effective_torque_norm | `9.8166663485015277e+05` | `7.8273600418514833e+06` |
| body2_y_rmse_vs_recurdyn | `0.65458890886616761` | `0.79270601530225326` |
| body2_z_rmse_vs_recurdyn | `0.41057562111363477` | `0.49087907197087843` |
| success_rate | `1` | `1` |

## 结论

简单的“一 patch 一个常量 pressure 未知量”在数值上可运行，`success_rate = 1`，但目前不是有效的精度改进方案。它确实显著减少了 open-gap 正接触力诊断行数，但同时过度压缩了 patch 内压力分布自由度，导致峰值接触力、峰值力矩、最大穿透和 Body2 相对 RecurDyn 的误差都变大。

因此该模式目前只保留为显式试验命令，不作为默认 `rev_joint_clearance` benchmark。

下一步如果继续沿 patch pressure solve 方向推进，不应继续只调单个标量参数，而应实现真正的局部 patch pressure solve：

1. patch 内保留 sample 级非负 pressure 未知量；
2. 使用 patch graph Laplacian 或面积质量矩阵耦合压力；
3. 保持真实面积积分的 force/torque；
4. 将 active-set 半光滑项纳入局部 patch Jacobian；
5. 在新局部 solve 明确改善指标前，继续以当前逐 sample SDF-NCP 结果作为基线。

## 局部 patch pressure solve 试验

随后实现了第二个实验后端：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_pressure
```

它与“一 patch 一个常量 pressure”不同，patch 内保留 sample 级 pressure unknown。每个 patch 内求解：

```text
0 <= p_i  perpendicular  g_i + c p_i + alpha (L p)_i >= 0
```

其中：

1. `p_i` 是 sample 级 pressure/local intensity；
2. `L` 是按 sample 面积加权、按行归一化的 patch graph Laplacian；
3. `c` 是局部 pressure compliance；
4. `alpha` 是 Laplacian pressure coupling compliance；
5. 非负与互补条件通过 smooth Fischer-Burmeister residual 求解；
6. 局部 Newton Jacobian 包含 `dPhi/dg * (c I + alpha L) + dPhi/dp`，因此 patch 内压力耦合进入局部半光滑 Jacobian；
7. 实际接触力仍使用真实面积积分 `force_i = p_i * area_i * n_i`。

实现位置：

```text
src/demos/core/demo_CH_sdf_ncp_benchmarks_openvdb.cpp
class ChSdfLocalPatchPressureContactForceSet
```

本轮测试了三个局部柔度尺度：

```text
rev_joint_clearance_local_patch_pressure          scale = 0.1
rev_joint_clearance_local_patch_pressure_balanced scale = 0.03
rev_joint_clearance_local_patch_pressure_stiff    scale = 0.01
```

与默认逐 sample descriptor 和上一个 patch aggregation 试验对比：

| case | contacts | max penetration | max force | max torque | closure ratio | Body2 y RMSE | Body2 z RMSE | relative angle RMSE | open positive-force rows | success |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| default sample descriptor | `42` | `1.1383846867829561e-04` | `2.2890600498895184e+06` | `9.8166663485015277e+05` | `0.95471260395959223` | `0.65458890886616761` | `0.41057562111363477` | `0.29545893287934782` | `720` | `1` |
| local patch pressure, scale 0.1 | `52` | `7.0523843169212341e-04` | `2.0789954435421305e+06` | `3.6085254474174359e+05` | `0.89477433735915135` | `0.6174563239033376` | `0.30036937444860412` | `0.26338260546805053` | `1984` | `1` |
| local patch pressure, scale 0.03 | `46` | `3.1717470847070217e-04` | `4.3961551596769076e+06` | `6.6854877228281822e+05` | `0.89489686326193862` | `0.74942497298675237` | `0.50849402109770425` | `0.34642991350802516` | `1919` | `1` |
| local patch pressure, scale 0.01 | `52` | `1.8800237739924341e-04` | `8.416215276364008e+06` | `2.3205898369636917e+06` | `0.9307449661996644` | `0.70507083684540339` | `0.36183253160570766` | `0.30399783087893717` | `2054` | `1` |
| one-pressure patch aggregation | `3` | `9.9698733538389206e-04` | `2.7665270590949859e+07` | `7.8273600418514833e+06` | `0.95191730379979389` | `0.79270601530225326` | `0.49087907197087843` | `0.35735784003422344` | `3` | `1` |

结论：

1. 局部 patch pressure solve 比“一 patch 一个常量 pressure”明显更合理，因为它保留了 sample 级压力分布自由度。
2. soft 版本在 Body2 y/z RMSE、相对向量相位误差和 patch torque closure 上优于默认逐 sample descriptor。
3. 代价是最大穿透和 open positive-force rows 变大，说明当前局部 solve 仍是显式 force backend，没有进入 Chrono descriptor 的隐式速度层。
4. 更硬的 pressure compliance 能压低穿透，但会增大力峰值并恶化曲线对齐。
5. 因此该方向“可行且有信号”，但当前实现还不能替代默认 benchmark。下一步应把局部 patch pressure solve 接成 descriptor/velocity-level coupling，或者在局部 solve 中加入 body effective mass/Delassus 近似，而不仅是几何 gap + pressure compliance。

## velocity-level Delassus 局部 patch solve 试验

进一步实现了 velocity-level effective-mass 近似：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_delassus_soft
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_delassus
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_delassus_stiff
```

该模式仍保留 sample 级 pressure unknown，但局部 residual 不再只是几何 gap compliance，而是 velocity-level 形式：

```text
0 <= p_i  perpendicular  v_n,i + beta min(g_i, 0) / dt
                         + sum_j W_ij area_j dt p_j
                         + c p_i
                         + alpha (L p)_i >= 0
```

其中：

```text
W = J M^{-1} J^T
```

`W_ij` 使用当前 Chrono body 的质量、局部惯量和 contact point moment arm 计算。力仍按面积积分施加：

```text
force_i = p_i area_i n_i
```

对比结果：

| case | contacts | max penetration | max force | max torque | closure ratio | Body2 y RMSE | Body2 z RMSE | relative angle RMSE | open positive-force rows | success |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| default sample descriptor | `42` | `1.1383846867829561e-04` | `2.2890600498895184e+06` | `9.8166663485015277e+05` | `0.95471260395959223` | `0.65458890886616761` | `0.41057562111363477` | `0.29545893287934782` | `720` | `1` |
| geometric local patch pressure, scale 0.1 | `52` | `7.0523843169212341e-04` | `2.0789954435421305e+06` | `3.6085254474174359e+05` | `0.89477433735915135` | `0.6174563239033376` | `0.30036937444860412` | `0.26338260546805053` | `1984` | `1` |
| Delassus local patch, scale 0.1 | `57` | `3.4915638389065862e-04` | `2.8537524935763376e+06` | `2.3133315258857101e+05` | `0.63243545853398475` | `0.6942807109248087` | `0.48515661269718402` | `0.32353107812075071` | `19` | `1` |
| Delassus local patch, scale 0.03 | `68` | `2.4572550319135189e-04` | `3.276637288523925e+06` | `4.2828262114734715e+05` | `0.87457607245495417` | `0.59161303977576518` | `0.43994188983870625` | `0.28123516161841067` | `57` | `1` |
| Delassus local patch, scale 0.01 | `47` | `2.0960810070391744e-04` | `3.2802429124773457e+06` | `6.1189928800567565e+05` | `0.89503349269627575` | `0.68645270293969296` | `0.50089366024869253` | `0.32438029363355053` | `212` | `1` |

结论：

1. Delassus/effective-mass 近似显著减少 open-gap 正压力诊断行数，说明 velocity-level residual 比纯几何 pressure compliance 更接近 NCP 接触律。
2. Delassus soft 版本的 patch torque closure 最好，`closure_ratio` 从默认 `0.9547` 降到 `0.6324`，但 Body2 曲线误差变差。
3. Delassus scale 0.03 版本在 Body2 y RMSE 和 relative angle RMSE 上优于默认 descriptor，并且 open positive-force rows 从 `720` 降到 `57`，但 Body2 z RMSE 和最大穿透仍差于默认。
4. 该结果说明：`W=J M^{-1} J^T` 是正确方向，但当前作为显式 force backend 仍不是完整 descriptor-level NCP。若要成为默认后端，需要把该 patch block 作为真正的 velocity-level constraint block 接入 descriptor，或者实现自定义 block solver，使 `W`、active-set 和 pressure coupling 在全局速度求解中同步迭代。

## descriptor-level patch block LCP 试验

随后把 patch 耦合推进到 Chrono descriptor/PSOR 速度层。新增的核心机制是：

```text
ChConstraint::ProjectGroupAndIncrementState(...)
ChSdfNcpConstraintContactSet::SetPatchPressureProjection(...)
```

实现方式不是外部显式力，也不是每步后处理压力，而是在 PSOR 的 unilateral projection 阶段对同一 patch 内的 active SDF samples 组装一个局部 Delassus block：

```text
W_ij = J_i M^{-1} J_j^T + cfm_i delta_ij

0 <= p_i  perpendicular  r_i + sum_j W_ij (p_j - p_j_old) >= 0
```

其中 `p_i` 仍是 sample pressure/local intensity，真实接触力仍为：

```text
force_i = p_i * area_i * n_i
torque_i = (x_i - x_body) cross force_i
```

为了保持 Chrono 全局速度状态一致，block solve 完成后会对 patch 内所有行调用 `IncrementState(delta_p_i)`。非 leader 行在同一 PSOR sweep 中会丢弃自己的 scalar tentative update，以避免同一 patch 被重复标量投影。

新增命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_projection
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_projection_strong
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_descriptor_patch_lcp
```

当前 3 秒 `rev_joint_clearance` 结果如下：

| case | max penetration | max torque | closure ratio | Body2 y RMSE | Body2 z RMSE | relative angle RMSE | open positive-force rows | success |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| default sample descriptor | `1.1383846867829561e-04` | `9.8166663485015277e+05` | `0.95471260395959223` | `0.65458890886616761` | `0.41057562111363477` | `0.29545893287934782` | `720` | `1` |
| descriptor patch block, relax 0.15 | `1.3281201245263219e-04` | `7.0830641117351246e+05` | `0.89507318112156575` | `0.66284644036658136` | `0.4216407467787377` | `0.30036274714546496` | `856` | `1` |
| descriptor patch block, relax 0.35 | `9.9448719993233681e-05` | `1.2803858147722422e+06` | `0.95287399478602408` | `0.66480974128545112` | `0.42412509929788417` | `0.30151039512504085` | `886` | `1` |
| descriptor patch block, relax 1.0 | `1.4315225416794419e-04` | `1.2233724414053506e+06` | `0.89509504856871158` | `0.66742758193319152` | `0.43244420438516362` | `0.30403926017997834` | `801` | `1` |
| local Delassus force backend, scale 0.03 | `2.4572550319135189e-04` | `4.2828262114734715e+05` | `0.87457607245495417` | `0.59161303977576518` | `0.43994188983870625` | `0.28123516161841067` | `57` | `1` |

结论：

1. descriptor-level patch block LCP 已经接入全局速度求解路径，算法位置比局部 force backend 正确。
2. 但当前结果没有改善与 RecurDyn/ideal revolute 的整体轨迹对齐。0.15 和 1.0 版本降低了 patch torque closure 指标，但 Body2 位置和相位误差略差于默认 descriptor。
3. 这说明当前主要误差不再只是“patch 内压力没有耦合”。更可能的主因是 contact manifold/active-set 的物理等价性不足：sample patch 的接触点、法向分布、开闭切换和 RecurDyn GGEOMCONTACT/理想铰链所需要的等效约束 wrench 仍不一致。
4. 因此 descriptor patch block LCP 保留为实验后端，不替代默认 `rev_joint_clearance` benchmark。下一步如果继续追求该算例对齐，应优先做 contact manifold 生成和 persistent patch tracking，而不是继续调单个 pressure relaxation 参数。

## 通用 contact manifold 修复

随后实现了通用 contact manifold 管理层，目标不是针对 `rev_joint_clearance` 特化，而是把 SDF 查询前端和 SDF-NCP descriptor 后端之间缺失的状态连续性补齐。

新增核心组件：

```text
ChSdfContactManifoldManager
```

实现位置：

```text
src/chrono/physics/ChSdfNcpConstraintContact.h
```

它负责：

1. 以 `contact_id` 维护 persistent sample state；
2. 把每步 connected-component 产生的 transient patch id 映射到 persistent patch id；
3. 使用 `gap_on/gap_off`、gap rate、预测 gap 和 release hysteresis 管理 active set；
4. 对 opening open-gap sample 清理 warm-start multiplier；
5. 保持 sample 面积语义一致，descriptor 后端继续使用

```text
force_i = pressure_i * area_i * normal_i
torque_i = (x_i - x_body) cross force_i
```

同时修复 descriptor warm start：

```text
旧实现：slot i 继承 slot i 的 lambda
新实现：contact_id 相同的 sample 继承上一帧 lambda
```

这避免了 active sample 顺序变化时把压力错误继承到另一个几何点。

`rev_joint_clearance` 新增输出：

```text
results/sdf_ncp_benchmarks/rev_joint_clearance/contact_manifold_diagnostics.csv
```

字段包括：

```text
time
candidate_count
active_count
patch_count
reused_patch_count
new_patch_count
released_sample_count
active_area
open_positive_lambda_rows
open_positive_lambda_force_sum
open_positive_lambda_force_max
```

本轮还新增了一个对照命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_manifold_no_prediction
```

该模式关闭 predictive activation，用于判断 open-gap positive multiplier 是否主要来自预测接触行。

3 秒 `rev_joint_clearance` 对比：

| case | max penetration | Body2 y RMSE | Body2 z RMSE | relative angle RMSE | closure ratio | open positive rows | success |
|---|---:|---:|---:|---:|---:|---:|---:|
| previous default descriptor baseline | `1.1383846867829561e-04` | `0.65458890886616761` | `0.41057562111363477` | `0.29545893287934782` | `0.95471260395959223` | `720` | `1` |
| persistent manifold, prediction on | `1.092120073735714e-04` | `0.664689965769521` | `0.4154834931757474` | `0.29977397861504324` | `0.93732100234833027` | `715` | `1` |
| persistent manifold, prediction off | `1.2926622002851218e-04` | `0.69709593911206746` | `0.51314637663604201` | `0.33045418318086572` | `0.8950721171101188` | `0` | `1` |

结论：

1. persistent manifold 修复没有破坏稳定性，`success_rate = 1`。
2. 默认 prediction-on 模式略微改善最大穿透和 patch torque closure，但 Body2 轨迹误差略大于上一版 baseline。
3. prediction-off 模式完全消除 open positive-force rows，但明显恶化 Body2 轨迹和相位误差。
4. 因此目前不能简单把“open positive rows = 0”作为唯一目标。对该间隙转动副，预测激活行虽然在严格 Signorini 诊断上不好看，但它在当前一阶 descriptor 装配里承担了接触响应的提前稳定作用。
5. 当前最大误差仍主要来自 contact manifold 与 RecurDyn/理想铰链所需等效 contact wrench 的不一致，而不是单纯来自 multiplier warm start 或 patch id 跳变。

本轮保留 prediction-on 作为默认 benchmark，因为它不使核心回归变差太多，并且 closure 有小幅改善；prediction-off 保留为诊断命令，不作为默认结果。

同一轮还把 simple_gear 的 generalized SDF-NCP 路径从“component 最小 sample id 作为 patch id”改为基于 overlap 的 persistent patch id。该路径本来已经有双向 sample persistent id、active sample memory 和面积积分，本轮补齐的是 patch id 连续性：

```text
GearPatchMemory::patch_sample_ids
GearPatchMemory::next_patch_id
```

`simple_gear` 一秒结果仍通过 benchmark：

```text
max_penetration = 8.2869746620417573e-09
success_rate = 1
gear22 omega analytic MAE = 0.065929736198861666
gear22 omega analytic RMSE = 0.10246898597734395
```

这说明 persistent patch id 改动没有损坏齿轮 benchmark，但齿轮相对解析解的剩余误差仍主要来自接触律/active-set/几何导数近似，而不是 patch id 跳变。

## rev_joint_clearance 双向 full mesh manifold 与 wrench 对齐

按后续排查顺序，新增了 `wrench_alignment.csv`，用于把 SDF-NCP patch 合力/合力矩与 Chrono ideal revolute 的 reaction wrench 对齐比较。ideal wrench 来源于：

```text
results/sdf_ncp_benchmarks/rev_joint_clearance_ideal_revolute/joint_reaction_wrench.csv
```

每个 `rev_joint_clearance` 变体现在输出：

```text
results/sdf_ncp_benchmarks/<case>/wrench_alignment.csv
```

字段包括：

```text
sdf_force
sdf_torque_about_body3_ref
ideal_force
ideal_torque_about_body3_ref
force_error_norm
torque_error_norm
force_cosine
torque_cosine
```

同时实现了双向 full mesh SDF manifold：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_bidirectional
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_bidirectional_no_prediction
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_bidirectional_patch_lcp
```

双向模式包含两类 sample：

```text
Body3 pin samples -> Body1 bore SDF
Body1 bore samples -> Body3 pin SDF
```

二者进入同一个 persistent manifold manager，但通过不同 `contact_id` 偏移保持 sample identity 不冲突。

3 秒结果摘要：

| case | max penetration | Body2 y RMSE | Body2 z RMSE | relative angle RMSE | force error RMSE | torque error RMSE | open positive rows | wall time | success |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| single direction manifold | `1.092120073735714e-04` | `0.664689965769521` | `0.4154834931757474` | `0.29977397861504324` | `1.406990723022499e+05` | `5.48507684584106e+04` | `715` | `21.06 s` | `1` |
| bidirectional manifold | `1.0937870683846995e-04` | `0.65374403756521127` | `0.40891517079854883` | `0.29484582137633564` | `1.2207962491264021e+05` | `5.793933350598337e+04` | `864` | `183.61 s` | `1` |
| bidirectional + descriptor patch LCP | `8.752616122364998e-05` | `0.68599106173205249` | `0.44244312543461084` | `0.31217254127366317` | `1.1315549717010533e+05` | `8.041379920757323e+04` | `908` | `194.12 s` | `1` |
| bidirectional, prediction off | `1.4453020412474871e-04` | `0.73838719738149328` | `0.54579726161995712` | `0.35078150822064108` | `1.7838966314884508e+05` | `6.526575125308152e+04` | `0` | `193.06 s` | `1` |

结论：

1. 双向 full mesh manifold 是有效方向：Body2 y/z RMSE、相对向量相位 RMSE 和 ideal wrench force error 都比单向默认略好。
2. 代价很高：墙钟从约 `21 s` 增加到约 `184 s`，说明 bore surface 反向采样需要更强的局部候选筛选/BVH narrow band。
3. 双向没有消除 open positive rows，反而从 `715` 增加到 `864`。这说明反向 manifold 扩大了接触候选，但当前 NCP active-set/预测行仍会在正 gap 区产生正压力。
4. 关闭 prediction 可以把 open positive rows 降为 `0`，但轨迹和 wrench error 明显变差；因此不能把 open-gap 清零作为唯一准则。
5. 在双向 manifold 上重新打开 descriptor patch LCP 后，最大穿透和 force error 继续下降，但 Body2 轨迹、相位和 torque error 变差。因此 patch LCP 仍不能作为默认后端。
6. 目前最可靠的下一步不是继续调 patch LCP，而是改进双向 manifold 的 contact candidate 分布和压力作用点：需要让反向 bore samples 只覆盖真实局部接触环带，并减少远离真实接触区域的预测 active rows。
