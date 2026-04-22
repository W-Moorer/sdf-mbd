# Milestone 2B.6: OpenVDB SDF Re-integration

## 1. 目的

本阶段将 2B.5 验证过的清洗后点式接触逻辑（activation band / force band 分离、几何一致性指标、参数扫描框架）重新接回 OpenVDB SDF 主线，验证几何一致性改进在 OpenVDB 场景下依然成立。

## 2. 场景定义

### 静态目标
- **OpenVDB sphere level set**: 中心在原点，半径 1.0m，voxel size 0.05m，half width 3 voxels（band width ±0.15m）

### 动态刚体
- **Sphere**: 半径 0.2m，密度 1000 kg/m³，质量 33.51 kg
- **初始位置**: y = 3.0m 自由落体

### 理论接触位置
- 动态球心理论平衡 Y = 1.0 + 0.2 = 1.2m

### 对照基线
- **解析球面 SDF**: 相同几何参数但无离散化误差的解析球面

## 3. 关键改进（沿用 2B.5）

1. **两带分离**: activation_band = 0.1m（候选检测），force_band ∈ {-0.01, 0.0, 0.005}m（实际施力）
2. **几何指标**: expected_y, y_error, normalized_y_error, min_phi, mean_active_phi
3. **参数扫描**: stiffness ∈ {1e4, 5e4, 1e5, 5e5}，damping ∈ {100, 200, 500}，force_band ∈ {-0.01, 0.0, 0.005}
4. **Reference vs OpenVDB 对比**: 在相同最佳参数下直接比较

## 4. 如何构建

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver"
cmake --build build --config Release --target demo_CH_sdf_point_contact_openvdb
```

## 5. 如何运行

```bash
cd "E:\workspace\Multi-body Dynamics Solver\Multi-body Dynamics Solver\build\bin\Release"
.\demo_CH_sdf_point_contact_openvdb.exe
```

## 6. 输出目录

```
out/milestone_03_6/
├── sdf_point_contact_openvdb_sweep.csv      (36 个 OpenVDB 配置的扫描结果)
├── sdf_point_contact_openvdb_best_run.csv   (OpenVDB 最佳参数完整仿真)
└── sdf_point_contact_reference_vs_openvdb.csv (解析参考 vs OpenVDB 对比表)

reports/milestone_03_6/
├── README_openvdb_reintegration.md   (本报告)
└── acceptance_report.md               (验收报告)
```

## 7. 为什么新建独立 demo 而非复用原 demo

`demo_CH_sdf_point_contact_openvdb.cpp` 是新建的，因为：
1. OpenVDB 需要 `openvdb::initialize()`、grid 创建、voxel 参数等初始化步骤，与解析平面版本完全不同
2. OpenVDB 需要额外的 `#include`（`openvdb/openvdb.h` 等）
3. 保持 2B.5 的解析平面 demo 干净，作为永久 reference baseline
4. 两者共享核心仿真逻辑（`RunSimulation` 函数和接触链路），只是 SDF 源不同

## 8. 执行顺序

1. 创建 OpenVDB sphere level set 作为静态接触目标
2. 在单一参数（2B.5 最佳参数）下验证接触逻辑
3. 对 OpenVDB 版本做参数扫描
4. 在相同最佳参数下运行解析参考和 OpenVDB 版本，直接对比
5. 分析 OpenVDB 特有误差来源
