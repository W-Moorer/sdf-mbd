# SDF-NCP field_contact 资产 benchmark 结果说明

## 2026-04-28 更新：正式 benchmark 收录顺序

### Benchmark 参考类型约定

当前 benchmark 先按参考来源分类，不能把所有算例都强行要求为 RecurDyn/RMD 对齐：

1. `RecurDyn reference-alignment benchmark`：资产本身包含 RMD、RecurDyn 输出 CSV 或等价 RecurDyn 前端定义。此类 benchmark 需要逐项对齐 part、marker、joint、motion、接触律和输出坐标。`cam_recurdyn_solid_contact` 属于这一类。
2. `Analytic-reference benchmark`：资产由解析参数、程序化几何或解析 envelope 定义，本身没有 RecurDyn 前端映射。此类 benchmark 应对齐解析解或解析几何关系，不应强行补 RMD。`eccentric_roller`、`onset_stress`、`headon_spheres` 和 `headon_spheres_mass_ratio` 属于这一类。
3. `Full-mesh numerical benchmark`：有完整 OBJ/OpenVDB 几何和稳定 rollout，但 reference 可能是数值曲线、field-contact 输出或暂未完整解析的外部模型。`simple_gear` 当前更接近这一类；若使用其 RMD/CSV 作为对照，需要另行完成前端对齐。

### Benchmark 1: `cam_recurdyn_solid_contact`

`cam_recurdyn_solid_contact` 收录为第一项正式 reference-alignment benchmark。该算例使用：

```text
assets/cam/simple_cam.rmd
assets/cam/models/cam_body1.obj
assets/cam/models/cam_body2.obj
assets/cam/data/cam_data.csv
```

建模路径为：

```text
RMD part/marker/joint/motion/surface/contact-law parser
OBJ/OpenVDB full mesh SDF frontend
Chrono MBD constraints: cam rotation motor + follower translational joint
RMD SolidContact1 K/C/KORDER/BPEN force-law baseline
RecurDyn reference CSV comparison
```

该 benchmark 的用途不是声称 frictionless hard SDF-NCP 已经复刻 RecurDyn penalty/damped contact，而是锁定前端建模一致性：几何、marker、joint、驱动方向、输出坐标和 reference 对比链路已经对齐。它是后续评估 `cam_recurdyn_compare` 这类 SDF-NCP hard-contact 方法差异的基准零点。

推荐默认参数：

```text
dt = 0.001
t_end = 3.0
voxel_size = 0.002
```

当前验收统计：

```text
follower_y RMSE       = 0.00059181137613861692
follower_y final err  = 8.1644475163844543e-05
follower_vy RMSE      = 0.030034144970505293
follower_vy final err = 0.0071187760060161323
success_rate          = 1
```

统一 runner 已支持长参考基准：

```text
python scripts\run_sdf_ncp_field_contact_benchmarks.py --include-reference
```

绘图脚本也会读取：

```text
results/sdf_ncp_benchmarks/cam_recurdyn_solid_contact/
```

### Benchmark 2: `eccentric_roller`

`eccentric_roller` 收录为第二项正式 SDF-NCP full-geometry benchmark。该算例使用：

```text
assets/eccentric_roller/eccentric_roller_model.json
assets/eccentric_roller/models/eccentric_disk_cam.obj
assets/eccentric_roller/models/roller_follower.obj
```

建模路径为：

```text
JSON parameter frontend
eccentric_disk_cam OBJ -> OpenVDB SDF
roller_follower OBJ -> surface graph
adaptive active-band connected patch
generic active-patch area-normalized SDF-NCP assembly
follower y-direction free dynamics
analytic envelope comparison
```

该 benchmark 的用途是验证完整 OBJ/OpenVDB 几何下的偏心轮-滚子接触、active patch 形成、gap/normal/lambda 诊断和 follower 响应趋势。它是 analytic-reference full-geometry benchmark：参考来自解析 eccentric envelope 和 JSON 参数，而不是 RecurDyn/RMD。因此它不需要 RecurDyn 前端映射；需要核对的是 OBJ 生成几何、OpenVDB SDF 和解析 envelope 是否一致。

默认参数来自 `eccentric_roller_model.json`：

```text
cam_radius        = 0.03
cam_eccentricity  = 0.006
roller_radius     = 0.01
motor_speed       = -2.0
dt                = 0.002
t_end             = 1.5707963267948966
friction          = 0.0
restitution       = 0.0
```

当前验收统计：

```text
max_penetration              = 1.943490860867314e-05
max_lambda_n                 = 2.4552490347727427
max_ncp_residual             = 4.5511024846491054e-05
max_complementarity_error    = 4.6394615375033133e-05
success_rate                 = 1
analytic envelope final diff = 1.1599401361341022e-05
relative final diff          = 0.00029330347047688534
```

运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe eccentric_roller
python scripts\run_sdf_ncp_field_contact_benchmarks.py --cases eccentric_roller
```

CTest 已注册：

```text
sdf_ncp_benchmark_eccentric_roller_full_geometry
labels: collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

输出路径：

```text
results/sdf_ncp_benchmarks/eccentric_roller/trajectory.csv
results/sdf_ncp_benchmarks/eccentric_roller/summary.csv
results/sdf_ncp_benchmarks/eccentric_roller/comparison.csv
results/sdf_ncp_benchmarks/figures/eccentric_roller/
```

### Benchmark 3: `onset_stress`

`onset_stress` 收录为第三项正式 analytic-reference SDF-NCP benchmark。该算例使用：

```text
assets/onset_stress/onset_stress_model.json
assets/onset_stress/models/onset_cam.obj
assets/onset_stress/models/roller_follower.obj
```

建模路径为：

```text
JSON parameter frontend
onset_cam OBJ -> OpenVDB SDF
roller_follower OBJ -> surface graph
adaptive active-band connected patch
generic active-patch area-normalized SDF-NCP assembly
follower y-direction free dynamics
analytic envelope and onset-time comparison
```

该 benchmark 的用途是验证接触激活时刻、gap 过零、active patch 启动和 SDF-NCP 互补残差。它不是 RecurDyn/RMD reference benchmark；参考来自 JSON 中的解析参数、解析 envelope 和 `target_onset_time`。

默认参数来自 `onset_stress_model.json`：

```text
cam_radius         = 0.03
cam_eccentricity   = 0.006
roller_radius      = 0.01
phase              = 3.141592653589793
motor_speed        = -2.0
dt                 = 0.001
t_end              = 0.45
target_onset_time  = 0.15
gravity_y          = 0.0
friction           = 0.0
restitution        = 0.0
```

当前解析对比统计：

```text
follower_y analytic RMSE       = 0.00066004111206987572
follower_y max abs error       = 0.0018128504933489636
follower_y final abs error     = 0.00069041375550264739
target onset time              = 0.14999999999999999
observed onset time            = 0.15221861543557469
onset time abs error           = 0.0022186154355746945
max_penetration                = 1.469743438065052e-09
max_ncp_residual               = 1.9645705151011818e-08
max_complementarity_error      = 6.5813389256593188e-09
success_rate                   = 1
```

运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe onset_stress
python scripts\run_sdf_ncp_field_contact_benchmarks.py --cases onset_stress
```

CTest 已注册：

```text
sdf_ncp_benchmark_onset_stress_full_geometry
labels: collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

输出路径：

```text
results/sdf_ncp_benchmarks/onset_stress/trajectory.csv
results/sdf_ncp_benchmarks/onset_stress/summary.csv
results/sdf_ncp_benchmarks/onset_stress/comparison.csv
results/sdf_ncp_benchmarks/onset_stress/analytic_comparison.csv
results/sdf_ncp_benchmarks/figures/onset_stress/
```

### Benchmark 4: `headon_spheres`

`headon_spheres` 收录为第四项正式 analytic-reference SDF-NCP benchmark。该算例使用：

```text
assets/headon_spheres/headon_spheres_model.json
assets/headon_spheres/models/ball_sphere.obj
```

建模路径为：

```text
analytic sphere SDF
3D translational two-body dynamics
single normal SDF-NCP contact
frictionless no-restitution impact
analytic common-velocity comparison
```

该 benchmark 的用途是验证解析球 SDF、gap、normal、NCP residual、非穿透和无恢复法向冲量响应。当前 SDF-NCP benchmark 不包含 restitution，因此正式解析参考是动量守恒的 no-restitution common normal velocity，而不是弹性碰撞速度。

当前验收统计：

```text
max_penetration                   = 1.4996337505124302e-13
max_lambda_n                      = 523.5986877637024
max_ncp_residual                  = 9.9426306126246851e-11
max_complementarity_error         = 7.8670589764498038e-11
success_rate                      = 1
analytic common vx sphere_a error = 2.5586043006509129e-07
analytic common vx sphere_b error = 2.5586042995406899e-07
```

运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres
python scripts\run_sdf_ncp_field_contact_benchmarks.py --cases headon_spheres
```

CTest 已注册：

```text
sdf_ncp_benchmark_headon_spheres
labels: collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

输出路径：

```text
results/sdf_ncp_benchmarks/headon_spheres/trajectory.csv
results/sdf_ncp_benchmarks/headon_spheres/summary.csv
results/sdf_ncp_benchmarks/headon_spheres/comparison.csv
results/sdf_ncp_benchmarks/figures/headon_spheres/
```

### Benchmark 5: `headon_spheres_mass_ratio`

`headon_spheres_mass_ratio` 收录为第五项正式 analytic-reference SDF-NCP benchmark。该算例使用：

```text
assets/headon_spheres_mass_ratio/headon_spheres_mass_ratio_model.json
assets/headon_spheres_mass_ratio/models/ball_sphere.obj
```

建模路径与 `headon_spheres` 相同，但改变两球质量比和初始速度，用于验证质量矩阵和法向冲量尺度对 SDF-NCP 响应的影响。

当前验收统计：

```text
max_penetration                   = 6.6557870326278135e-14
max_lambda_n                      = 698.13160995855878
max_ncp_residual                  = 9.9161553285687926e-11
max_complementarity_error         = 4.6532711036623817e-11
success_rate                      = 1
analytic common vx sphere_a error = 2.947859186752666e-07
analytic common vx sphere_b error = 1.4739296017030057e-07
```

运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks.exe headon_spheres_mass_ratio
python scripts\run_sdf_ncp_field_contact_benchmarks.py --cases headon_spheres_mass_ratio
```

CTest 已注册：

```text
sdf_ncp_benchmark_headon_spheres_mass_ratio
labels: collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

输出路径：

```text
results/sdf_ncp_benchmarks/headon_spheres_mass_ratio/trajectory.csv
results/sdf_ncp_benchmarks/headon_spheres_mass_ratio/summary.csv
results/sdf_ncp_benchmarks/headon_spheres_mass_ratio/comparison.csv
results/sdf_ncp_benchmarks/figures/headon_spheres_mass_ratio/
```

### 后续 benchmark 候选优先级

| 优先级 | case | 推荐定位 | 当前状态 | 需要补齐的前端一致性工作 |
| --- | --- | --- | --- | --- |
| 已收录 1 | `cam_recurdyn_solid_contact` | RecurDyn front-end alignment benchmark | 已有 3 秒 reference comparison | 后续可继续解析 RecurDyn 内部 SolidContact 细节 |
| 已收录 2 | `eccentric_roller` | analytic-reference full OBJ/OpenVDB active-patch benchmark | 已有 full OBJ/OpenVDB rollout，summary 稳定 | 无需 RMD 映射；继续核对 OBJ/OpenVDB 与解析 envelope 的一致性 |
| 已收录 3 | `onset_stress` | analytic-reference onset/gap activation benchmark | 已有 full OBJ/OpenVDB rollout、analytic envelope RMSE 和 onset time error | 无需 RMD 映射；继续核对解析 onset 定义、OBJ/OpenVDB 和相位定义 |
| 已收录 4 | `headon_spheres` | analytic-reference sphere SDF and no-restitution impact benchmark | 已有解析共同速度误差、summary 稳定 | 无需 RMD 映射；保留弹性速度作为非目标对照 |
| 已收录 5 | `headon_spheres_mass_ratio` | analytic-reference mass-ratio sphere SDF benchmark | 已有解析共同速度误差、summary 稳定 | 无需 RMD 映射；验证质量比和冲量尺度 |
| 6 | `simple_gear` | 双向齿面 patch、多点接触和旋转自由度 benchmark | 已有 full OBJ/OpenVDB 双向 patch，能稳定运行 | 需要完整对齐 RecurDyn gear RMD 的 joint、驱动、齿轮输出坐标和接触律；当前不应作为高重合 reference benchmark |
| 7 | `rev_joint_clearance` | 关节间隙、多体约束与接触耦合 benchmark | asset 存在，SDF-NCP benchmark 尚未实现 | 需要先建立 RMD/Chrono joint-clearance 前端映射，再接通用 SDF-NCP descriptor 后端 |

建议的正式 benchmark 队列：

1. `cam_recurdyn_solid_contact`：已经收录，作为前端一致性基准。
2. `eccentric_roller`：已经收录，作为 analytic-reference full-geometry active-patch 稳定性基准。
3. `onset_stress`：已经收录，作为 analytic-reference 接触激活时刻 benchmark，重点看 gap/onset 和解析 envelope 误差。
4. `headon_spheres`：已经收录，作为解析球 SDF 和无恢复碰撞基础 benchmark。
5. `headon_spheres_mass_ratio`：已经收录，作为质量比影响和冲量尺度 benchmark。
6. `simple_gear`：适合作为多点 patch 和双向 SDF 查询 benchmark，但应在 RecurDyn 前端完全对齐后再作为 reference-alignment benchmark。
7. `rev_joint_clearance`：工程价值高，但前端约束映射工作量最大，应排在 simple gear 之后。

## 2026-04-28 更新：cam OpenVDB 体素尺寸敏感性

为检查 `cam` 与 RecurDyn reference 之间的剩余微小误差是否主要来自 SDF/OpenVDB 离散化，本轮新增了仅改变 `voxel_size` 的手动 benchmark 模式。所有模式使用相同的 RMD part/marker/joint/motion、相同的 OBJ 几何、相同的 Chrono MBD 约束路径，以及相同的 RMD `SolidContact1` 力律；唯一改变的是 cam OpenVDB level-set 的体素尺寸。

运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact_voxel_004
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact_voxel_002
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact_voxel_001
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact_voxel_0005
```

对比脚本：

```text
python scripts\compare_sdf_ncp_recurdyn_reference.py --case cam --result-case <case_name>
```

最新统计如下：

| result case | voxel size | follower_y RMSE | follower_y max error | follower_y final error | follower_vy RMSE | follower_vy max error | follower_vy final error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `cam_recurdyn_solid_contact_voxel_004` | `0.004` | `0.00061101931906964806` | `0.0053156010063735715` | `0.00017702443583084149` | `0.0309194248833181` | `0.55767604392213199` | `0.0044680395217360641` |
| `cam_recurdyn_solid_contact_voxel_002` | `0.002` | `0.00059181137613861692` | `0.0050523734192073344` | `8.1644475163844543e-05` | `0.030034144970505293` | `0.54944858677160568` | `0.0071187760060161323` |
| `cam_recurdyn_solid_contact_voxel_001` | `0.001` | `0.00060134676363598647` | `0.0050208447963062119` | `2.1740440938444738e-05` | `0.030142885089695558` | `0.54848041699612826` | `0.01025726261932651` |
| `cam_recurdyn_solid_contact_voxel_0005` | `0.0005` | `0.00060961810049083595` | `0.00502176268524962` | `6.1555693597892258e-08` | `0.030317117423723661` | `0.54851239457599466` | `0.010504921546210227` |

结论：

1. 体素尺寸从 `0.004` 细化到 `0.0005` 后，`follower_y final error` 明显下降，说明 SDF/OpenVDB 离散化确实影响局部几何相位和末端位移。
2. `follower_y RMSE` 和 `follower_vy RMSE` 没有随体素尺寸单调下降，并在 `0.002` 到 `0.0005` 之间进入平台区。
3. 当前剩余 RMSE 和速度瞬态误差不主要由 OpenVDB 体素尺寸控制，更可能来自 RecurDyn SolidContact 内部 patch/阻尼/恢复/速度输出细节与当前 SDF force baseline 的实现差异。
4. 对后续正式对比，`voxel_size = 0.002` 仍是较合理的默认值；若关注最终位移误差，可使用 `0.001` 或 `0.0005` 做高分辨率复核，但不应期待它单独消除速度瞬态误差。

## 1. 当前有效 benchmark

当前有效 SDF-NCP benchmark 包括：

```text
headon_spheres
headon_spheres_mass_ratio
cam
eccentric_roller
onset_stress
simple_gear
```

其中：

1. `headon_spheres` 和 `headon_spheres_mass_ratio` 使用解析球 SDF。解析球 SDF 与完整球体几何等价，不属于 reduced 目标轨迹模型。
2. `cam` 使用 `assets/cam` 中的完整 OBJ/OpenVDB SDF，并走 active-band connected patch SDF-NCP 路径；当前窄激活带和几何采样下有效 patch 可能退化为单样本。
3. `eccentric_roller`、`onset_stress` 和 `simple_gear` 使用完整 OBJ/OpenVDB SDF 的多点 patch 接触，不再使用 top-N 压缩候选点。
4. 旧 reduced `eccentric_roller`、`simple_gear` CSV/图表已清理，不再作为正式结果。

## 2. 建模方式摘要

| case | asset | SDF-NCP 建模方式 | 对比来源 | 状态 |
| --- | --- | --- | --- | --- |
| `headon_spheres` | `assets/headon_spheres` | 解析球 SDF，3D 无摩擦平动自由动力学 | analytic no-restitution common velocity | 通过 |
| `headon_spheres_mass_ratio` | `assets/headon_spheres_mass_ratio` | 解析球 SDF，质量比变化 | analytic no-restitution common velocity | 通过 |
| `cam` | `assets/cam` | `cam_body1.obj` OpenVDB SDF + `cam_body2.obj` follower surface graph，active-band connected patch | `cam_data.csv` | 通过 |
| `eccentric_roller` | `assets/eccentric_roller` | `eccentric_disk_cam.obj` OpenVDB SDF + `roller_follower.obj` active-band connected patch | analytic envelope | 通过 |
| `onset_stress` | `assets/onset_stress` | `onset_cam.obj` OpenVDB SDF + `roller_follower.obj` active-band connected patch，记录 onset time | target onset time / analytic envelope | 通过 |
| `simple_gear` | `assets/simple_gear` | RMD 坐标解析，GEAR21 OpenVDB SDF，GEAR22 齿面 active-band connected patch，GEAR22 RX 自由动力学 | `Gear22.csv` | 通过 |

当前 patch 选择规则是局部 active band，而不是 top-N 压缩：

```text
min_phi = min_i phi_i
若 min_phi <= activation_band:
    active = { i | phi_i <= min_phi + activation_band }
    patch = connected_components(active)
```

进入 SDF-NCP residual 的是 patch 内全部 surface samples。若某一时刻只有一个 sample 落入激活带，该 patch 会自然退化为单样本；这不是 top-N 筛选，而是当前网格分辨率、SDF 带宽和激活带参数共同决定的结果。

## 3. 当前结果路径

统一 summary：

```text
results/sdf_ncp_benchmarks/benchmark_summary.csv
results/sdf_ncp_benchmarks/benchmark_comparison_summary.csv
```

各 case：

```text
results/sdf_ncp_benchmarks/headon_spheres/
results/sdf_ncp_benchmarks/headon_spheres_mass_ratio/
results/sdf_ncp_benchmarks/cam/
results/sdf_ncp_benchmarks/eccentric_roller/
results/sdf_ncp_benchmarks/onset_stress/
results/sdf_ncp_benchmarks/simple_gear/
```

图表：

```text
results/sdf_ncp_benchmarks/figures/<case_name>/
```

## 4. 关键数值摘要

最近一次 `benchmark_summary.csv` 中的客观统计：

| case | max_penetration | max_lambda_n | max_ncp_residual | max_complementarity_error | mean_iterations | success_rate |
| --- | --- | --- | --- | --- | --- | --- |
| `headon_spheres` | `1.4996337505124302e-13` | `523.5986877637024` | `9.9426306126246851e-11` | `7.8670589764498038e-11` | `2.7222777222777221` | `1` |
| `headon_spheres_mass_ratio` | `6.6557870326278135e-14` | `698.13160995855878` | `9.9161553285687926e-11` | `4.6532711036623817e-11` | `2.8531468531468533` | `1` |
| `cam` | `1.0000076144933701e-07` | `4047.4349484399163` | `1.0000082598105031e-07` | `7.4211392130154048e-06` | `1.0079365079365079` | `0.99206349206349209` |
| `eccentric_roller` | `9.3801281764172018e-07` | `0.44056665262368833` | `0.0030494781001044768` | `0.0011734927232233626` | `1.5699745547073791` | `0.96055979643765899` |
| `onset_stress` | `1.254193193744868e-09` | `0.23689911162294985` | `1.9740870182034345e-08` | `3.3068568565695699e-09` | `0.042128603104212861` | `1` |
| `simple_gear` | `1.0006110073845775e-07` | `60.882345287614271` | `3.1738614483473959e-07` | `6.1920155857559984e-06` | `3.5098039215686274` | `1` |

这些值只说明当前局部 SDF-NCP residual 在各 case 中保持有限、无 NaN、穿透和互补误差有界。它们不表示与 pressure-field 或 RecurDyn reference 的逐点等价。

## 5. 对比结果解释

`benchmark_comparison_summary.csv` 包含 field-contact/reference 与 SDF-NCP 的同指标对比。

需要谨慎解释以下差异：

1. sphere case 的正式 reference 是无恢复共同法向速度；`comparison.csv` 中保留解析弹性碰撞速度只是 field_contact paper example 的非目标对照，不能作为当前 frictionless SDF-NCP 的验收指标。
2. `cam` 使用完整 mesh SDF 和 active-band connected patch，但窄带内若只有一个 mesh sample 达到最近接触区域，数值上会退化为单样本 patch；它仍无摩擦、无 pressure-field patch pressure distribution。
3. `eccentric_roller` 对比的是解析 envelope。SDF-NCP 使用完整 OBJ sample 查询和 active patch，因此它是几何/响应趋势对比，不是目标轨迹拟合。
4. `onset_stress` 当前 onset time 由 full OBJ/OpenVDB min-gap 过零估计。
5. `simple_gear` 当前与 `Gear22.csv` 差异较大。原因包括：
   - 当前是无摩擦 SDF-NCP；
   - 当前只使用 `GEAR22 surface -> GEAR21 SDF` 单向 active tooth patch；
   - 没有 field-contact pressure patch 分布；
   - 没有 Chrono/RecurDyn 全局约束路径；
   - hard NCP 可能产生比 pressure-field 更强的瞬时法向冲量。

因此 `simple_gear` 当前结果可以作为 full tooth mesh SDF-NCP 后端验证，但不能作为与 field-contact 完全一致的动力学复刻。

## 6. 生成图表

每个 case 至少生成：

```text
gap_vs_time.png
lambda_vs_time.png
penetration_vs_time.png
complementarity_vs_time.png
```

额外图：

```text
headon_spheres/relative_distance_vs_time.png
headon_spheres/velocity_vs_time.png
cam/follower_response_vs_time.png
eccentric_roller/follower_response_vs_time.png
onset_stress/follower_response_vs_time.png
simple_gear/gear_angle_vs_time.png
simple_gear/angular_velocity_vs_time.png
simple_gear/contact_force_vs_time.png
```

## 7. CTest 状态

当前 benchmark CTest：

```text
sdf_ncp_benchmark_headon_spheres
sdf_ncp_benchmark_headon_spheres_mass_ratio
sdf_ncp_benchmark_cam_full_geometry
sdf_ncp_benchmark_eccentric_roller_full_geometry
sdf_ncp_benchmark_onset_stress_full_geometry
sdf_ncp_benchmark_simple_gear_full_geometry
```

统一标签：

```text
collision;field_contact;sdf_ncp;sdf_ncp_benchmark
```

## 8. 复现命令

```text
cmake --build build --config Release --target demo_CH_sdf_ncp_benchmarks
cmake --build build --config Release --target demo_CH_sdf_ncp_benchmarks_openvdb
python scripts\run_sdf_ncp_field_contact_benchmarks.py
python scripts\plot_sdf_ncp_field_contact_benchmarks.py
ctest --test-dir build -C Release -L sdf_ncp_benchmark --output-on-failure
```

## 9. 当前限制和下一步

当前限制：

1. 无摩擦。
2. 未接入 Chrono 主 contact container/descriptor。
3. `simple_gear` 仍是第一阶段单向 tooth-SDF active patch，不是双向一致 contact pair。
4. 当前局部 Newton 使用有限差分 Jacobian。
5. `rev_joint_clearance` 尚未实现。

推荐下一步：

1. 将 `simple_gear` 升级为双向 `GEAR22->GEAR21` 与 `GEAR21->GEAR22` 一致化候选。
2. 为 OpenVDB SDF-NCP 后端补解析/半解析 Jacobian，降低有限差分成本。
3. 将 `cam` 的激活带、采样密度和 patch 宽度做系统参数扫描，避免窄带下长期退化为单样本 patch。
4. 在不修改 Chrono 主 contact path 的前提下，设计 field-contact NCP container 适配层。

## 10. 2026-04-28 修订：通用后端后的 OpenVDB benchmark 结果

本轮将 `cam`、`eccentric_roller`、`onset_stress`、`simple_gear` 的求解路径统一到 `ChSdfNcpGeneralizedProblem` / `SolveSdfNcpGeneralizedProblem`。case 代码只负责 OBJ/OpenVDB/RecurDyn 前端映射、active patch 候选和 `J_i = grad phi dx/dq` 计算。

### 当前 summary

最近一次 `build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe all` 返回 0，并生成以下结果：

| case | max_penetration | max_lambda_n | max_ncp_residual | max_complementarity_error | mean_iterations | success_rate |
| --- | --- | --- | --- | --- | --- | --- |
| `cam` | `5.9059573686681688e-05` | `4047.4250501524621` | `8.2568281365515074e-05` | `0.0054638643699803221` | `0.91269841269841268` | `0.99206349206349209` |
| `eccentric_roller` | `1.943490860867314e-05` | `2.4552490347727427` | `4.5511024846491054e-05` | `4.6394615375033133e-05` | `8.8511450381679388` | `1` |
| `onset_stress` | `1.469743438065052e-09` | `0.4737631710813951` | `1.9645705151011818e-08` | `6.5813389256593188e-09` | `0.16629711751662971` | `1` |
| `simple_gear` | `1.0001883765653474e-07` | `2538.9387857558236` | `2.469618950134651e-07` | `0.00017550745199741425` | `1.392156862745098` | `1` |

这些数值说明当前 generalized SDF-NCP backend 在完整 OBJ/OpenVDB 几何上可以完成 rollout，并保持有限 penetration、NCP residual 和 complementarity diagnostics。

### 与 RecurDyn/reference 的当前对应情况

当前还没有完全对应 RecurDyn 曲线：

| case | reference metric | reference | sdf_ncp | 当前结论 |
| --- | --- | --- | --- | --- |
| `cam` | follower y at nearest final reference | `0.26318247` | `0.40237822700936549` | 偏高约 53%，尚未对应 |
| `cam` | follower vy at nearest final reference | `-0.13720668` | `0.039940748421364698` | 速度方向/幅值尚未对应 |
| `simple_gear` | GEAR22 omega RX at nearest final reference | `-0.92880052999999996` | `0.28271065937761297` | 符号和幅值均未对应 |

这不是因为当前 benchmark 使用 reduced geometry；`cam` 和 `simple_gear` 已经使用完整 OBJ/OpenVDB 几何。主要未对齐项是 RecurDyn 约束、接触律、驱动方向、自由度定义、摩擦/阻尼/恢复参数以及输出坐标定义。

### 后续必须完成的对齐工作

1. 对 `cam` 从 `simple_cam.rmd` 提取完整 joint/motion/contact 设置，并验证 follower 的坐标正方向、初始偏置和输出列定义。
2. 对 `simple_gear` 核对 RMD 中 GEAR21/GEAR22 的 RX 正方向、驱动符号、接触 action/base 定义和 `Gear22.csv` 的角速度符号约定。
3. 确认是否必须纳入 RecurDyn/field-contact 的摩擦、阻尼、恢复系数或 pressure-field 分布；若第一篇仍坚持无摩擦 NCP，应在论文中明确这是模型差异，不能声称逐曲线复刻。
4. 将 `rev_joint_clearance` 按同一 generalized backend 映射，而不是实现新的 case-specific residual。
## 11. 2026-04-28 RecurDyn 前端映射对齐进展

本轮按 RecurDyn 前端建模逐项核对了 `cam` 和 `simple_gear` 的 RMD 信息，并将明确的映射错误修正到 SDF-NCP benchmark 前端。

## 2026-04-29 更新：simple_gear 小时间步稳定性验收

simple_gear 当前不再把 `dt=0.0001` 的发散结果作为有效 benchmark。已新增 dt sweep 验收路径：

```text
demo_CH_sdf_ncp_benchmarks_openvdb.exe simple_gear_dt_001
demo_CH_sdf_ncp_benchmarks_openvdb.exe simple_gear_dt_0005
demo_CH_sdf_ncp_benchmarks_openvdb.exe simple_gear_dt_0001
python scripts/validate_sdf_ncp_simple_gear_dt_sweep.py --no-run
```

新增的后端策略是通用的，不是 simple_gear 曲线拟合：

1. `dt >= 0.001` 保留当前 force-level hard SDF-NCP 基线，保证默认结果不变差。
2. `dt < 0.001` 启用 impulse / velocity-level mixed SDF-NCP。
3. mixed residual/Jacobian 使用前向自动微分计算。
4. 前端向通用后端提供 `normal_velocity_offset`，用于表达 prescribed driver 对 gap-rate 的贡献。
5. force 权重转换为 `dt * force_weight` 的 impulse 权重，避免小时间步下尺度爆炸。

当前验收摘要：

| case | dt | success_rate | max_penetration | MAE vs analytic | RMSE vs analytic | final omega22 |
|---|---:|---:|---:|---:|---:|---:|
| `simple_gear_dt_001` | `0.001` | `1` | `8.2869746620417573e-09` | `0.065929736198861666` | `0.10246898597734395` | `-1.0511661888263963` |
| `simple_gear_dt_0005` | `0.0005` | `1` | `1.0006613138102693e-07` | `0.063336256860321918` | `0.13101143283368741` | `-0.96575938464339262` |
| `simple_gear_dt_0001` | `0.0001` | `1` | `1.0003219585996703e-07` | `0.055549608448871655` | `0.11795684023219731` | `-0.96901542949800457` |

这说明当前实现满足本轮验收标准：`dt=0.001` 不比当前基线差，`dt=0.0005` 与 `dt=0.0001` 不发散，`success_rate = 1`，最大穿透有界，MAE/RMSE 没有因减小时间步而爆炸。

### 2026-04-29 补充：OpenVDB 几何导数与半光滑 active set

simple gear 小时间步路径现在不再用角度有限差分估计 `dJ/dtheta22`。OpenVDB 前端在 `ChSdfContactSampleQuery` 中输出 raw gradient 和 Hessian；gear 前端通过链式法则计算：

```text
dn = (I - n n^T) H dx / ||raw_grad||
dJ/dtheta22 = d(r x n)/dtheta22
```

双向 patch 的两种查询方向都使用同一套几何导数。mixed SDF-NCP 后端继续用轻量前向 AD 组装代数 Jacobian，并把 `contact.jacobian_velocity_derivative` 纳入 normal velocity 和 generalized impulse residual。active samples 在 Newton 内层固定，在外层通过 persistent id、hysteresis 和 warm start 做半光滑 active-set 更新。

重新运行 `python scripts\validate_sdf_ncp_simple_gear_dt_sweep.py` 后得到：

| case | dt | success_rate | max_penetration | MAE vs analytic | RMSE vs analytic | final omega22 |
|---|---:|---:|---:|---:|---:|---:|
| `simple_gear_dt_001` | `0.001` | `1` | `8.2869746620417573e-09` | `0.065929736198861666` | `0.10246898597734395` | `-1.0511661888263963` |
| `simple_gear_dt_0005` | `0.0005` | `1` | `1.0006613138102693e-07` | `0.063319570651797885` | `0.13101056370576455` | `-0.96574456900269146` |
| `simple_gear_dt_0001` | `0.0001` | `1` | `1.0003219585996703e-07` | `0.055550703486665146` | `0.11795585683594351` | `-0.96901368213534234` |

这些数值说明几何导数接入后没有破坏既有验收：`dt=0.001` 保持基线结果；小时间步仍不发散；`success_rate = 1`；误差没有随 `dt` 减小而爆炸。

补充说明：`dt=0.001` 的 force-level hard SDF-NCP 路径不再支持 analytic geometry Jacobian 开关。诊断中强制打开该选项会使该基线的 success rate 和误差显著变差，因此已删除该入口；当前只让小时间步 mixed backend 默认消费 OpenVDB 几何导数，避免破坏既有 benchmark。

### 已对齐项

1. `cam` 的驱动速度不再使用 `cam_model.json` 中的 `3.1415926`，而是解析 `assets/cam/simple_cam.rmd` 中 `RevJoint1.RMotion` 的 `FUNCTION = 3\`，当前 `cam_rmotion_velocity_constant = 3`。
2. `simple_gear` 的 GEAR21 驱动速度解析自 `assets/simple_gear/simple gear.rmd` 中 `RevJoint1.RMotion` 的 `FUNCTION = 1\`，当前 `gear21_rmotion_velocity_constant = 1`。
3. `simple_gear` 修正了 `GEAR22 surface -> GEAR21 SDF` 方向的广义力符号。该方向中，GEAR22 是 action surface，作用在 GEAR22 上的法向力方向应为 `+normal_world`，因此 `J = +(r_22 x n)_x`。反向 `GEAR21 surface -> GEAR22 SDF` 中，GEAR22 是 target body，作用在 GEAR22 上的力为 `-normal_world`，因此仍使用 `J = -(r_22 x n)_x`。
4. `cam` 的 active band 从窄带单点候选扩大到更宽的 active patch 设置，避免长期只输出一个 sample 的退化接触候选。

### 当前结果变化

`simple_gear` 的符号错误修正后，GEAR22 RX 角速度已经从先前的正方向错误变为与参考同号：

```text
reference omega_rx = -0.92880052999999996
sdf_ncp omega_rx   = -0.80187227200161471
relative error     = 0.136658253197148
```

这说明 simple gear 之前最大的问题是 action/base contact 方向下的广义 Jacobian 符号，而不是 NCP 后端公式。

`cam` 在对齐 RMD 驱动速度后仍未对应参考曲线：

```text
reference follower_y = 0.26318247
sdf_ncp follower_y   = 0.40191685834810942
relative error       = 0.52714144809154428

reference follower_vy = -0.13720668
sdf_ncp follower_vy   = 0.090969331811635343
```

因此 cam 当前不对应的主要原因已经从“驱动速度是否对齐”转移到接触模型和约束响应差异：当前 SDF-NCP 是无摩擦、无 RecurDyn penalty/damping 接触律的法向互补模型；RecurDyn `SolidContact1` 中存在 `K = 100000000`、`C = 10000`、`KORDER = 2`、`BPEN = 1e-5` 等接触律参数。若目标是逐曲线复刻 RecurDyn，必须决定是否在 SDF-NCP 外增加等价阻尼/正则化项，或建立专门的 RecurDyn contact-law baseline。纯无摩擦 hard NCP 不应被期待与 penalty/damped SolidContact 曲线完全一致。

### 最新 summary

| case | max_penetration | max_lambda_n | max_ncp_residual | max_complementarity_error | success_rate | 当前曲线对应状态 |
| --- | --- | --- | --- | --- | --- | --- |
| `cam` | `0.0023489876184612513` | `6032.3776152832388` | `0.0032137091529615626` | `8.552245400518693` | `0.96825396825396826` | 未对应 RecurDyn follower 曲线 |
| `simple_gear` | `1.0001349437516183e-07` | `0.0098601987208012006` | `1.0002389870046362e-07` | `1.0025966731214024e-07` | `1` | 已同号，仍有约 13.7% 幅值差异 |

### 下一步必须继续对齐的项

1. `cam`：解析并记录 `SolidContact1` 的 `K/C/KORDER/BPEN/RDF`，决定 SDF-NCP 论文模型是否应保持 hard NCP，还是另设 RecurDyn contact-law 对照。
2. `cam`：核对 `RevJoint1` marker `I=Ground.Marker1`、`J=Body1.Marker1` 和 `TraJoint1` marker `I=Ground.Marker2`、`J=Body2.Marker1` 的 joint frame 方向，确认 follower 输出 `Y:Pos_TY` 是否正是当前 `follower_cm0.y + state.y`。
3. `simple_gear`：继续核对 `GeoSurContact1` 的 `IGGEOMID/JGGEOMID`、action/base marker 和 RecurDyn 输出列符号；当前符号已修正，但幅值差异仍可能来自无摩擦 hard NCP 与 `GGEOMCONTACT` penalty/damping 接触律差异。
4. 将上述 RecurDyn 解析结果写成 machine-readable mapping report，避免后续再靠硬编码或人工记忆维护。
## 2026-04-28 更新：cam descriptor path 结果状态

`cam` 已从 reduced/local generalized residual 切换到 Chrono 多体约束 descriptor path。当前输出仍位于：

```text
results/sdf_ncp_benchmarks/cam/trajectory.csv
results/sdf_ncp_benchmarks/cam/summary.csv
results/sdf_ncp_benchmarks/cam/comparison.csv
```

当前 `summary.csv` 显示该路径可以完成 rollout，并输出有限的 penetration、lambda、NCP residual 和 complementarity diagnostics。最近一次 `cam` descriptor run 的关键统计为：

```text
max_penetration = 9.9550711456686258e-04
max_lambda_n = 1.1009977292603206e+02
max_ncp_residual = 3.1475281310591545e-03
max_complementarity_error = 1.0356163354528472e-01
success_rate = 1
```

与 RecurDyn reference 的当前对比仍未完全对应：

```text
follower_y reference = 0.26318247
follower_y sdf_ncp   = 0.40338595901451391

follower_vy reference = -0.13720668
follower_vy sdf_ncp   = 0.058314895668142318
```

这说明 cam 已进入完整 Chrono 多体约束求解路径，但尚未完成 RecurDyn 前端建模的一一对齐。主要待对齐项是：

1. RMD joint marker frame 和 `REULER` 方向到 Chrono joint frame 的精确映射。
2. follower translational joint 的自由轴是否与 RecurDyn `TraJoint1` 完全一致。
3. RecurDyn `SolidContact1` penalty/damping/contact-law 与当前 frictionless SDF-NCP hard contact 的模型差异。
4. reference CSV 输出坐标是否是 body CM、marker 坐标或其他派生量。

因此当前结果不能写成“已对应 RecurDyn 曲线”。准确表述是：SDF-NCP descriptor 后端已经接入 Chrono 多体约束求解路径，cam 可以用完整 OBJ/OpenVDB 几何运行；曲线差异主要来自 RecurDyn 前端建模和接触律尚未逐项等价。

## 2026-04-28 更新：通用 RecurDyn RMD 映射层用于 cam 对比

为避免 cam benchmark 继续依赖手填坐标，本轮将 cam 的前端建模改为通用 RMD 解析路径。解析内容包括：

```text
PART: id, name, MASS, CM marker id, QG, REULER, IP
MARKER: id, name, PART, QP, REULER
CSURFACE: id, NAME, RM reference marker
JOINT: id, NAME, I marker, J marker, type
MOTION: id, NAME, JOINT, ROTATION/TRANSLATION, VELOCITY, FUNCTION
SOLIDCONTACT: IFLOAT/JFLOAT, ICSURFACEID/JCSURFACEID, BPEN, RDF, K, KORDER, C
GRAVITY: IGRAV/JGRAV/KGRAV
```

cam 现在只通过实体名选择 `Body1`、`Body2`、`RevJoint1`、`TraJoint1`、`RevJoint1.RMotion` 和 `SolidContact1`；实际 CM、joint frame、surface reference marker、gravity、motion speed 和 contact-law audit 均来自 `assets/cam/simple_cam.rmd`。这不是 cam 专用后端，而是 RecurDyn 前端映射层到通用 SDF-NCP descriptor 后端的适配。

新增输出：

```text
results/sdf_ncp_benchmarks/cam/rmd_mapping.csv
```

该文件记录 RMD 中解析出的 part、marker、joint、motion、surface 和 solid contact 参数，用于检查当前 Chrono/SDF-NCP benchmark 与 RecurDyn 输入模型是否逐项一致。

当前 cam benchmark 使用：

```text
ChSystemNSC
ChBodyAuxRef from RMD PART/CM/QG/IP
ChLinkMotorRotationSpeed from RMD RevJoint1 + RevJoint1.RMotion
ChLinkMateGeneric from RMD TraJoint1 marker frame
ChSdfNcpConstraintContactSet descriptor unilateral contact
OBJ/OpenVDB full geometry from RMD CSURFACE reference markers
```

最新短时 benchmark 结果仍未与 RecurDyn 曲线对应：

```text
reference follower_y = 0.26318247
sdf_ncp follower_y   = 0.4057833683628615

reference follower_vy = -0.13720668
sdf_ncp follower_vy   = 0.094847407556277352
```

当前差异不能再归因于 cam 使用 reduced/local residual；cam 已经使用通用 RMD 前端映射和 Chrono descriptor 后端。剩余主要差异是模型层面的：

1. RecurDyn `SolidContact1` 是 penalty/damping contact law，RMD 中已解析到 `K = 100000000`、`C = 10000`、`KORDER = 2`、`BPEN = 1e-5`、`RDF = 0.25`。
2. 当前 SDF-NCP benchmark 是 frictionless hard normal complementarity，不使用 `K/C/KORDER/BPEN/RDF` 作为接触律。
3. 当前 CTest benchmark 仍使用短时间窗；RecurDyn reference CSV 覆盖到 3.0 s。若要做完整曲线对齐，需要增加非 CTest 的 full-duration reference compare runner。
4. 仍需确认 RecurDyn translational joint 的自由轴约定与 Chrono `ChLinkMateGeneric` 轴约定是否完全一致；当前根据 RMD marker frame 使用 local-Z free translational axis。

因此，“通用对比上”的下一步不是给 cam 加专用补丁，而是把 RecurDyn contact-law 和 full-duration reference comparison 做成可选的通用前端/对比模式，并保持 SDF-NCP 后端独立。

## 2026-04-28 更新：完整 RecurDyn 曲线对比输出

新增手动完整对比入口：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_compare
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact
python scripts\compare_sdf_ncp_recurdyn_reference.py --case cam --result-case cam_recurdyn_compare
python scripts\compare_sdf_ncp_recurdyn_reference.py --case cam --result-case cam_recurdyn_solid_contact
```

`cam_recurdyn_compare` 使用完整 RMD/OBJ/OpenVDB/Chrono MBD 前端和 frictionless SDF-NCP descriptor 后端，时间步设置为 RMD `HMAX = 0.001`，输出 3 秒完整曲线。`cam_recurdyn_solid_contact` 使用同一套 RMD/OBJ/OpenVDB/Chrono MBD 前端，但接触力使用从 RMD `SolidContact1` 解析出的 `K/C/KORDER/BPEN`，作为 RecurDyn contact-law 对照基线。两者都不是 reduced cam ODE。

完整对比脚本会按 reference CSV 时间列插值 SDF/Chrono trajectory，并输出：

```text
recurdyn_reference_comparison_timeseries.csv
recurdyn_reference_comparison_summary.csv
figures/reference_comparison/follower_y_overlay.png
figures/reference_comparison/follower_y_error.png
figures/reference_comparison/follower_vy_overlay.png
figures/reference_comparison/follower_vy_error.png
```

图中使用高对比线型：RecurDyn reference 为黑色实线，benchmark 为橙色虚线，误差为绿色曲线。

### 当前完整 3 秒对比统计

frictionless SDF-NCP descriptor path：

```text
follower_y RMSE      = 0.069390710046443149
follower_y final err = 0.057072294926321765
follower_vy RMSE     = 0.5438221097880569
follower_vy final err= 0.46327602078646657
```

RMD SolidContact-law SDF baseline：

```text
follower_y RMSE      = 0.067600173629993446
follower_y final err = 0.058454141571343154
follower_vy RMSE     = 0.52537450106549732
follower_vy final err= 0.4493143881367202
```

这说明当前已经完成了完整时间范围、完整几何、完整 RMD 前端映射、Chrono 多体约束路径和清晰可视化对比；但数值曲线仍未达到“与 RecurDyn 高度重合”。即使把 RMD `K=100000000`、`C=10000`、`KORDER=2`、`BPEN=1e-5` 纳入 SDF force baseline，误差只小幅改善，说明剩余差异很可能来自 RecurDyn SolidContact 内部的接触 patch/压力积分/阻尼限幅/恢复因子实现细节，或者 reference 输出坐标定义仍有未解析项。

因此当前可用于论文/工程报告的准确结论是：对比链路已经完整，图表对比清晰；但 SDF-NCP hard-contact 模型不应被声称已复刻 RecurDyn penalty/damped SolidContact 曲线。下一步若目标是数值曲线高度重合，需要继续解析 RecurDyn SolidContact 的内部力律细节，而不是修改 NCP 后端为 case-specific cam 拟合器。
## 2026-04-28 更新：cam 前端建模一致性修正完成

本轮完成了 `cam` 算例的 RecurDyn 前端建模一致性修正。关键问题不是 SDF-NCP 后端，也不是 OBJ/OpenVDB 几何查询，而是旋转驱动方向的前端映射约定：

```text
RecurDyn rotational motion: J marker relative to I marker
Chrono motor initialization in this benchmark: body_J, body_I
Equivalent Chrono motor speed: - RecurDyn FUNCTION constant
```

代码中新增的通用映射函数：

```text
MapRecurDynRotationVelocityToChronoMotorSpeed(...)
```

该函数只属于 RecurDyn/RMD 前端映射层，不改变 `ChSdfNcpConstraintContactSet`、`ChSdfNcpGeneralizedProblem` 或 Chrono 主 contact container。`comparison.csv` 现在分别记录 RMD 原始速度和 Chrono 等价速度，避免把 `+3 rad/s` 与 `-3 rad/s` 误报为建模差异。

重新运行的默认完整 3 秒对比命令为：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_compare
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe cam_recurdyn_solid_contact
python scripts\compare_sdf_ncp_recurdyn_reference.py --case cam --result-case cam_recurdyn_compare
python scripts\compare_sdf_ncp_recurdyn_reference.py --case cam --result-case cam_recurdyn_solid_contact
```

### 最新 RecurDyn 曲线对比统计

frictionless SDF-NCP descriptor path：

```text
follower_y RMSE       = 0.003405188444209653
follower_y max error  = 0.018874933237232105
follower_y final err  = 0.0010264619191152069
follower_vy RMSE      = 0.082653461054846361
follower_vy max error = 0.74392767919863989
follower_vy final err = 0.02473873209611371
```

RMD SolidContact-law SDF baseline：

```text
follower_y RMSE       = 0.00059181137613861692
follower_y max error  = 0.0050523734192073344
follower_y final err  = 8.1644475163844543e-05
follower_vy RMSE      = 0.030034144970505293
follower_vy max error = 0.54944858677160568
follower_vy final err = 0.0071187760060161323
```

对应输出文件：

```text
results/sdf_ncp_benchmarks/cam_recurdyn_compare/trajectory.csv
results/sdf_ncp_benchmarks/cam_recurdyn_compare/comparison.csv
results/sdf_ncp_benchmarks/cam_recurdyn_compare/recurdyn_reference_comparison_summary.csv
results/sdf_ncp_benchmarks/cam_recurdyn_compare/figures/reference_comparison/

results/sdf_ncp_benchmarks/cam_recurdyn_solid_contact/trajectory.csv
results/sdf_ncp_benchmarks/cam_recurdyn_solid_contact/comparison.csv
results/sdf_ncp_benchmarks/cam_recurdyn_solid_contact/recurdyn_reference_comparison_summary.csv
results/sdf_ncp_benchmarks/cam_recurdyn_solid_contact/figures/reference_comparison/
```

当前结论：

1. `cam` 的 RMD part、marker、joint、motion、gravity、surface reference marker、OBJ/OpenVDB 几何和 follower 输出坐标已经通过同一套前端映射进入 Chrono MBD 路径。
2. 修正旋转 motion 符号后，默认 `cam_recurdyn_compare` 不再依赖诊断用 reverse motor 命令，位移曲线已与 RecurDyn reference 对齐到毫米量级。
3. `cam_recurdyn_solid_contact` 使用同一套前端建模，并启用 RMD `SolidContact1` 的 `K/C/KORDER/BPEN` 力律；该基线与 RecurDyn 曲线重合最好，用于验证前端建模一致性。
4. frictionless SDF-NCP descriptor path 与 RecurDyn SolidContact-law baseline 之间仍存在速度瞬态差异，这是接触律差异，不应再归因于前端建模不一致。
5. SDF-NCP 后端保持通用：后端只接收 gap、normal、weight、body pair 和 unilateral normal constraint 数据；cam 相关内容限制在 RMD/OBJ/OpenVDB 前端映射与 reference 对比层。

## 2026-04-29 更新：AABB BVH 粗检测结果

本轮实现了 backend-neutral AABB BVH 粗检测层。它不改变 SDF-NCP 接触方程，只改变 active sample 候选生成方式：

```text
source surface sample BVH
  -> transformed target SDF local bounds query
  -> candidate sample ids
  -> phi-only active-band scan
  -> full OpenVDB query only for active samples
```

当前接入范围：

1. `cam` descriptor path：follower surface samples 通过 BVH 与 cam OpenVDB bounds 做粗检测。
2. `eccentric_roller` / `onset_stress` local rollout：roller follower samples 通过 BVH 与 cam-like OpenVDB bounds 做粗检测。
3. `simple_gear`：GEAR22->GEAR21 和 GEAR21->GEAR22 两个方向均使用各自 source surface BVH。

为避免精度损失，active band 的 AABB 查询使用 voxel/activation/hysteresis padding；当 BVH 候选不足以满足最小 patch 样本数时，回退到全 surface phi-only 扫描。

### simple_gear dt sweep 结果

重新运行：

```text
python scripts\validate_sdf_ncp_simple_gear_dt_sweep.py
```

得到：

| case | dt | success_rate | max_penetration | MAE vs analytic | RMSE vs analytic | final omega22 | runtime_seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| `simple_gear_dt_001` | `0.001` | `1` | `8.2869746620417573e-09` | `0.065929736198861666` | `0.10246898597734395` | `-1.0511661888263963` | `12.344617400000001` |
| `simple_gear_dt_0005` | `0.0005` | `1` | `1.0006613138102693e-07` | `0.063319570651797885` | `0.13101056370576455` | `-0.96574456900269146` | `11.5995366` |
| `simple_gear_dt_0001` | `0.0001` | `1` | `1.0003219585996703e-07` | `0.055550703486665146` | `0.11795585683594351` | `-0.96901368213534234` | `55.016859500000002` |

`dt=0.001` 的 MAE/RMSE 与 BVH 接入前的基线一致；因此本轮粗检测没有改变该算例的解析误差。墙钟耗时明显下降，主要原因是全 surface scan 不再对每个 sample 都计算 OpenVDB gradient/Hessian。

### 当前限制

1. BVH 仍是静态 local-space sample BVH，每步通过 transformed AABB 查询，不是连续碰撞检测。
2. broad phase 只保证不会改变 active patch 的 intended band；极端高速运动仍需要后续 swept AABB 或时间连续 active set。
3. 真正的 SDF-NCP 收敛精度仍由 active patch、SDF 离散化、接触律和 Newton/AD Jacobian 决定；BVH 只解决候选搜索效率。

## 2026-04-29 更新：rev_joint_clearance SDF-NCP benchmark

`rev_joint_clearance` 已接入完整几何 SDF-NCP benchmark。该算例使用：

```text
assets/rev_joint_clearance/rev_clearance_joint.rmd
assets/rev_joint_clearance/models/body1_subtract1_centered.obj
assets/rev_joint_clearance/models/body3_cylinder1_centered.obj
assets/rev_joint_clearance/data/body2.csv
assets/rev_joint_clearance/data/body3.csv
```

建模路径：

```text
RMD Body1/Body2/Body3 + Fixed1/Fixed2
Body1.Subtract1 OBJ -> OpenVDB SDF
Body3.Cylinder1 OBJ -> surface graph + AABB BVH
Body3 samples -> Body1 SDF query
ChSdfNcpConstraintContactSet descriptor path
RecurDyn body2/body3 CSV comparison
```

当前完整 3 秒运行命令：

```text
build\bin\Release\demo_CH_sdf_ncp_benchmarks_openvdb.exe rev_joint_clearance
python scripts\plot_sdf_ncp_field_contact_benchmarks.py
```

输出：

```text
results/sdf_ncp_benchmarks/rev_joint_clearance/trajectory.csv
results/sdf_ncp_benchmarks/rev_joint_clearance/summary.csv
results/sdf_ncp_benchmarks/rev_joint_clearance/comparison.csv
results/sdf_ncp_benchmarks/rev_joint_clearance/recurdyn_reference_comparison_timeseries.csv
results/sdf_ncp_benchmarks/rev_joint_clearance/recurdyn_reference_comparison_summary.csv
results/sdf_ncp_benchmarks/figures/rev_joint_clearance/
```

当前统计：

| quantity | value |
|---|---:|
| `success_rate` | `1` |
| `max_penetration` | `0` |
| `max_lambda_n` | `46292.474861248957` |
| `max_ncp_residual` | `0.012568139367502403` |
| `max_complementarity_error` | `90.07090329970265` |
| `runtime_seconds` | `15.4644972` |
| `body2_y_rmse_vs_recurdyn` | `0.67076679030986863` |
| `body2_z_rmse_vs_recurdyn` | `0.45840028100800234` |
| `body3_y_rmse_vs_recurdyn` | `0.019253805588352407` |
| `body3_z_rmse_vs_recurdyn` | `0.00456934007894573` |

解释：

1. `Body3` 的位置曲线已经能给出稳定、有限的 SDF-NCP descriptor rollout。
2. `Body2` 的误差仍明显偏大，需要继续核对 RecurDyn 输出坐标是否是 CM、marker 还是 joint frame 坐标。
3. 当前 hard SDF-NCP 与 RecurDyn `GGEOMCONTACT K/C/KORDER/BPEN` penalty/damping 接触律不同，因此不能把曲线差异全部归因于 SDF query。
4. 当前已经进入完整 Chrono MBD 约束路径，不是 reduced local residual。
