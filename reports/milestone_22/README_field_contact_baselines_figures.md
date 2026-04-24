# Milestone 22: Baselines and Paper Figure Generation

This milestone adds the last comparison layer requested for the TeX paper:
baseline variants plus automatic figure generation from benchmark CSV files.

## Baselines

The multi-slip/interlock benchmark now includes five variants:

- `field_primitive`: connected field primitive extraction plus runtime tracking,
  split/merge handling, minimal-rotation transport, and Coulomb return.
- `raw_sample_contacts`: every active surface sample is its own contact.
- `normal_bin_fragments`: connected fragments are split by quantized local
  normals, and the quantized normals are used for force integration.
- `convex_decomposition_proxy`: active contacts are split by fixed body-local
  convex proxy pieces. This approximates the contact bookkeeping produced by a
  convex-decomposition surrogate.
- `chrono_traditional_proxy`: a Chrono-style point-manifold baseline selecting
  a small set of deepest separated contact points with quantized feature
  normals. It is deterministic and shares the analytic benchmark geometry.

The last two are proxy baselines for the paper comparison. They are not a full
replacement for a future native Chrono collision-system experiment, but they
make the intended failure modes measurable in the same deterministic scene set.

Implementation:

`src/demos/core/demo_CH_field_contact_multislip_interlock.cpp`

Run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_multislip_interlock
.\build\bin\Release\demo_CH_field_contact_multislip_interlock.exe
```

Benchmark CSV files:

```text
out/milestone_21/field_contact_multislip_interlock_frames.csv
out/milestone_21/field_contact_multislip_interlock_patches.csv
out/milestone_21/field_contact_multislip_interlock_summary.csv
out/milestone_21/field_contact_multislip_interlock_comparison.csv
```

## Current Baseline Ratios

Ratios are baseline divided by `field_primitive`. Values above 1 mean the
baseline is worse.

```text
guide_rail_sliding,convex_decomposition_proxy:
  mean_patch_count_change_ratio=2.50000000
  rms_normal_jump_ratio=1.96548987
  topology_event_ratio=60.91666667

guide_rail_sliding,chrono_traditional_proxy:
  mean_patch_count_change_ratio=2.00000000
  rms_normal_jump_ratio=2.56431563
  force_oscillation_ratio=3.32208715
  torque_oscillation_ratio=5.56797490

nested_interlock,chrono_traditional_proxy:
  mean_patch_count_change_ratio=4.00000000
  rms_normal_jump_ratio=2.32116370
  force_oscillation_ratio=2.12114201
  torque_oscillation_ratio=3.04096182

multi_patch_rolling_sliding,convex_decomposition_proxy:
  mean_patch_count_change_ratio=1.63636364
  rms_normal_jump_ratio=1.85410253
  topology_event_ratio=23.95000000

multi_patch_rolling_sliding,chrono_traditional_proxy:
  rms_normal_jump_ratio=2.68802800
  force_oscillation_ratio=4.11632485
  torque_oscillation_ratio=3.28031978
```

All reported variants keep:

```text
max_tangential_force_ratio=1.00000000
```

so the baseline comparison does not rely on violating the Coulomb bound.

## Figure Generator

Script:

`scripts/generate_field_contact_figures.py`

Run:

```powershell
python scripts\generate_field_contact_figures.py --input out\milestone_21 --output out\milestone_22\figures
```

The script has no third-party Python dependency. It reads the benchmark CSV
files and writes SVG figures plus an HTML index:

```text
out/milestone_22/figures/index.html
```

Generated per-scene figures:

- `*_force_jump.svg`
- `*_torque_jump.svg`
- `*_patch_count.svg`
- `*_event_timeline.svg`
- `*_coulomb_ratio.svg`
- `*_energy_ratio.svg`

Generated summary figures:

- `summary_force_oscillation.svg`
- `summary_torque_oscillation.svg`
- `summary_patch_churn.svg`
- `summary_normal_jump.svg`
- `summary_energy_ratio.svg`

Current run generated 23 SVG figures.

## Interpretation

Milestones 20-22 now form a reproducible paper pipeline:

- Milestone 20: feature-switching scene with split/merge history
  non-amplification.
- Milestone 21: multi-contact non-convex benchmarks showing compact persistent
  field primitives versus fragmented contacts.
- Milestone 22: convex-decomposition and Chrono-style point-manifold baselines,
  plus automatic paper figure generation for force jump, torque jump, patch
  count, event timeline, Coulomb ratio, and inherited energy ratio.

The strongest paper claim supported by the current numbers is:

```text
field primitives reduce contact-set churn and normal/response oscillation
relative to fragmented feature and point-manifold baselines while preserving
Coulomb feasibility and inherited-history non-amplification.
```
