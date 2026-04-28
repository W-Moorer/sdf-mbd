# Multi-body Dynamics Solver

本项目基于 Project Chrono 10.0.0 代码树，面向多体动力学中的复杂几何接触问题，扩展了基于场的接触原语与 OpenVDB 稀疏窄带 SDF 接触后端。当前仓库重点支持网格几何的接触力计算、接触斑块跟踪、切向历史继承、摩擦响应验证，以及论文示例的可复现实验。

## 主要内容

- 保留 Chrono 核心多体动力学框架、碰撞、约束、求解器、FEA、车辆、传感器等原有模块。
- 新增 `chrono::fieldcontact` 场接触原语工具，提供表面采样图、接触斑块提取、法向/切向接触力、黏滑状态和拓扑事件处理。
- 新增场接触运行时状态机，负责持久化接触原语 ID、分裂/合并分类、切向历史传输、接触力与力矩聚合。
- 集成 OpenVDB 稀疏窄带 SDF 示例，用于点接触、面片接触、网格 SDF、非凸网格、齿轮、凸轮、球碰、铰链间隙等案例。
- 提供 Python 脚本生成论文图表、比较 Chrono 原生接触与场接触结果，并执行回归验收。

## 关键源码

- `src/chrono/collision/ChFieldContactPrimitives.h`
  场接触基础数据结构和算法，包括表面图构建、接触原语提取、法向/切向接触设置、黏滑更新与历史继承工具。

- `src/chrono/collision/ChFieldContactRuntime.h`
  场接触运行时跟踪器，包括持久 ID、拓扑事件、切向历史存储、斑块级结果和双向接触体受力汇总。

- `src/demos/core/`
  核心验证与 OpenVDB 示例程序，例如 `demo_CH_field_contact_regression`、`demo_CH_field_contact_feature_switching`、`demo_CH_field_contact_primitives_openvdb`、`demo_CH_simple_gear_openvdb`、`demo_CH_rev_joint_clearance_openvdb` 等。

- `paper_example/`
  论文复现实验入口，包含 `paper_examples_openvdb.cpp`、案例输入、回归清单、绘图脚本和已生成图表。

- `assets/`
  复杂几何案例输入，包括凸轮、偏心滚子、球碰、齿轮、转动副间隙等模型、参考数据和网格文件。

- `scripts/`
  验收、比较和绘图脚本，用于场接触回归、摩擦验证、自由动力学响应、CAD 网格案例和论文结果生成。

## 构建要求

基础构建需要：

- CMake 3.18 或更新版本
- 支持 C++17 的编译器
- Python 3，用于回归检查和图表生成

OpenVDB 相关示例还需要：

- OpenVDB
- TBB
- Imath

在 Windows/MSVC 环境下，当前 CMake 脚本会优先查找 `C:/vcpkg/installed/x64-windows` 中的 OpenVDB、TBB 和 Imath 库。如果依赖安装在其他位置，请通过 CMake 包路径或相应变量让 `find_package(OpenVDB)` 能够发现它们。

## 配置与编译

建议使用独立构建目录：

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

只构建常用场接触目标：

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
cmake --build build --config Release --target demo_CH_field_contact_feature_switching
cmake --build build --config Release --target demo_CH_field_contact_primitives_openvdb
cmake --build build --config Release --target demo_CH_paper_examples_openvdb
```

如果 CMake 未找到 OpenVDB，`*_openvdb` 示例会被跳过，基础 Chrono 目标和非 OpenVDB 场接触验证仍可构建。

## 运行示例

从仓库根目录运行已编译程序。Windows Release 构建的可执行文件通常位于 `build/bin/Release`：

```powershell
build\bin\Release\demo_CH_field_contact_regression.exe
build\bin\Release\demo_CH_field_contact_feature_switching.exe
build\bin\Release\demo_CH_field_contact_primitives_openvdb.exe
build\bin\Release\demo_CH_simple_gear_openvdb.exe
build\bin\Release\demo_CH_rev_joint_clearance_openvdb.exe
```

不同示例会在 `out/` 下写出 CSV、摘要文件或中间结果。绘图和验收脚本默认也从仓库根目录解析输入与输出路径。

## 论文示例复现

`paper_example` 目录锁定了当前稀疏 SDF 接触后端的论文基准输入，包含：

- `eccentric_roller`
- `headon_spheres`
- `headon_spheres_mass_ratio`
- `simple_gear`

运行流程：

```powershell
cmake --build build --config Release --target demo_CH_paper_examples_openvdb
build\bin\Release\demo_CH_paper_examples_openvdb.exe
python paper_example\plot_paper_examples.py --project-root .
python paper_example\check_paper_examples.py --project-root .
```

主要输出：

- `out/paper_example_dynamic_benchmarks/sparse_sdf_frames.csv`
- `out/paper_example_dynamic_benchmarks/comparison_summary.csv`
- `out/paper_example_dynamic_benchmarks/timing_summary.csv`
- `paper_example/figures/*.pdf`

验收阈值记录在 `paper_example/manifest.json` 中。

## 回归测试

仓库默认提供场接触相关 CTest 注册开关：

```powershell
cmake -S . -B build -DCH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release -L field_contact --output-on-failure
```

回归测试覆盖基础场接触、摩擦验证、Chrono 原生接触对比、匹配网格、自由动力学响应、CAD 网格案例、Chrono 原始网格基线和论文示例。部分测试依赖 OpenVDB 和 Python；依赖缺失时，对应测试或目标会被 CMake 跳过。

## 目录结构

```text
.
├── src/                    # Chrono 源码与本项目场接触扩展
├── src/demos/core/          # 场接触、OpenVDB 和核心多体动力学示例
├── assets/                 # 本项目几何模型、参考曲线和案例输入
├── paper_example/           # 论文示例、回归清单、绘图脚本和图表
├── scripts/                 # 验收、比较、绘图和结果生成脚本
├── data/                   # Chrono 原始数据、模型、纹理、YAML 示例
├── paper/                  # 中文论文草稿和编译产物
├── build/                  # 本地构建输出
└── out/                    # 示例和脚本生成的运行结果
```

## 许可证

本仓库继承 Project Chrono 的 BSD 风格许可证。许可证文本见 `LICENSE`。
