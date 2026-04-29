# rev_joint_clearance patch wrench closure 结果记录

日期：2026-04-30

## 目的

本轮目标是让 SDF-NCP 后端在完整几何接触 patch 上自行求出更合理的压力分布，而不是继续调单点刚度或使用参考曲线作为求解目标。

实现仍保持非 AI、无摩擦、无 Chrono 主 contact container 改动。

## 实现变化

1. `rev_joint_clearance` 的 active sample 选择从 `min_phi + band` 改为支持绝对 active band。
   这样接触 patch 不再长期退化成过少 sample。

2. `ChSdfLocalPatchPressureContactForceSet` 的 patch wrench closure 从局部 patch 质心改为关于接触 body 参考点闭合。
   这是因为 Body3 动力学实际接收的是关于 Body3 参考点的 wrench。

3. patch pressure closure 增加非负二次投影：

   ```text
   minimize ||F(p) - F0||^2 + w_tau ||tau_ref(p)||^2 + r ||p - p0||^2
   subject to p_i >= 0
   ```

   其中：
   - `p_i` 是 sample pressure；
   - `F0` 是原 SDF-NCP/Delassus pressure solve 得到的 patch 合力；
   - `tau_ref` 是关于接触 body 参考点的合力矩；
   - 该步骤不使用 RecurDyn 或 ideal wrench 作为求解目标。

## 新增命令

```bat
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque100
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque250
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque500
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque1000
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance_local_patch_wrench_closure_torque2000
```

## 当前结果摘要

与 `rev_joint_clearance_local_patch_delassus` 相比，`torque1000` 是当前最好的非 oracle 结果：

| case | Body2 y RMSE | Body2 z RMSE | force RMSE | torque RMSE | max penetration | success |
|---|---:|---:|---:|---:|---:|---:|
| local patch Delassus | 0.681902 | 0.514181 | 97581.99 | 5926.53 | 8.13e-05 | 1 |
| wrench closure torque250 | 0.426940 | 0.345989 | 104305.25 | 3235.75 | 4.65e-04 | 1 |
| wrench closure torque500 | 0.417620 | 0.343288 | 104067.79 | 2787.36 | 4.49e-04 | 1 |
| wrench closure torque1000 | 0.398955 | 0.333527 | 105776.34 | 2578.76 | 3.67e-04 | 1 |
| wrench closure torque2000 | 0.567611 | 0.356549 | 133025.58 | 2630.00 | 4.95e-04 | 1 |

重新归一化 wrench residual 后，`torque1000` 的复跑结果为：

| case | Body2 y RMSE | Body2 z RMSE | force RMSE | torque RMSE | max penetration | success |
|---|---:|---:|---:|---:|---:|---:|
| wrench closure torque1000 | 0.409301 | 0.336346 | 108230.79 | 2765.86 | 5.74e-04 | 1 |

## Augmented residual 尝试

本轮还尝试把 wrench closure 从后处理推进到 patch 主残差中：

```text
R(p) = [
  Phi_eps(c_i(p), p_i)
  sqrt(w_F) * (F(p) - F0)
  sqrt(w_T) * Tau_body_ref(p)
]
```

其中 `F0` 仍来自同一 patch 的非 oracle FB/Delassus 解，不使用 RecurDyn 或 ideal wrench。

结果：

| case | Body2 y RMSE | Body2 z RMSE | force RMSE | torque RMSE | max penetration | success |
|---|---:|---:|---:|---:|---:|---:|
| augmented wrench closure torque500 | 0.758367 | 0.505663 | 114717.63 | 19205.04 | 1.61e-03 | 1 |
| augmented wrench closure torque1000 | 0.651628 | 0.477326 | 96386.94 | 9905.06 | 5.33e-04 | 1 |

结论：当前 augmented 主残差版不如后处理式 `torque1000`。它降低了一部分 force RMSE，但 torque RMSE 和 Body2 轨迹均变差，说明 FB complementarity residual 与 wrench residual 的尺度/目标耦合还没有处理好。

## 结论

1. 扩大 active patch 后，逐帧冻结扫描证明接触几何 patch 具备更好的 wrench 可达性。
2. 关于 Body3 参考点的 pressure redistribution 明显改善了 Body2 运动趋势。
3. 当前最优仍是后处理式 `torque1000`；augmented residual 版本目前无效。
4. 当前最优 `torque1000` 仍未完全对齐 RecurDyn/ideal 曲线，说明还缺少更完整的全局 velocity-level descriptor/block constraint 耦合。
5. OpenVDB voxel 从 `0.002` 提升到 `0.001` 没有改善逐帧 oracle wrench 误差，因此当前主误差不是 SDF 分辨率。

## 下一步

将上述非负 pressure redistribution 从 local force backend 推进到 descriptor/block constraint 后端，并让 patch pressure 与 bilateral joint/driver constraints 在同一速度层同步迭代。
