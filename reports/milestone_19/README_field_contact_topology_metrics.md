# Milestone 19: Contact Topology Consistency Metrics

This milestone implements the TeX contact-topology consistency layer as a
reusable runtime metric accumulator and wires it into the split/merge demo.

## Runtime API

Header:

`src/chrono/collision/ChFieldContactRuntime.h`

New types:

- `FieldContactTopologyMetricsSummary`
- `FieldContactTopologyMetricsAccumulator`

Usage:

```cpp
FieldContactTopologyMetricsAccumulator metrics;

FieldContactStepResult step = tracker.Evaluate(graph, queries, reference, settings);
metrics.Accumulate(step);

FieldContactTopologyMetricsSummary summary = metrics.GetSummary();
```

The accumulator consumes `FieldContactStepResult` and does not alter contact
forces, history, or primitive IDs.

## Metrics

The summary currently reports:

- frames and active frames;
- total and mean primitive count;
- max primitive count;
- total newborn/death/merge/split events;
- mean and max absolute primitive-count change;
- primitive persistence ratio;
- mean inherited source-weight sum;
- mean, RMS, and max matched-normal jump angle;
- mean, RMS, and max force jump;
- mean, RMS, and max torque jump;
- normalized force and torque oscillation indices;
- stick/slip switch count on persistent primitive IDs;
- max Coulomb ratio, energy-gate ratio, and inherited-energy ratio.

## Demo Integration

`demo_CH_field_contact_split_merge_openvdb` now writes these fields to:

`out/milestone_15/field_contact_split_merge_summary.csv`

Current Milestone 15/19 run:

```text
topology_active_frames,1600
topology_mean_patch_count,1.55250000
topology_primitive_persistence_ratio,0.99919485
topology_mean_source_weight_sum,0.99725596
topology_mean_abs_patch_count_change,0.00125078
topology_max_abs_patch_count_change,1
topology_rms_normal_jump_angle,0.00160325
topology_max_normal_jump_angle,0.01398138
topology_mean_force_jump,101.76888121
topology_rms_force_jump,169.74970338
topology_max_force_jump,967.09394872
topology_force_oscillation_index,0.00302145
topology_mean_torque_jump,14.52977227
topology_rms_torque_jump,29.43034318
topology_max_torque_jump,197.11811038
topology_torque_oscillation_index,0.01757526
topology_stick_slip_switches,10
```

The split/merge topology event remains:

```text
2 patches -> 1 patch -> 2 patches
```

with one merge frame and one split frame.

## Regression Coverage

`demo_CH_field_contact_regression` includes a deterministic three-frame fixture
that fixes:

- frame and active-frame counts;
- total, mean, and max primitive count;
- count-variation mean and max;
- primitive persistence ratio;
- source-weight mean;
- mean/RMS/max normal jump;
- mean/RMS/max force jump;
- mean/RMS/max torque jump;
- stick/slip switch count;
- max Coulomb ratio propagation.

Run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
ctest --test-dir build -C Release -R demo_CH_field_contact_regression --output-on-failure
```

## Interpretation

Milestones 15-19 now provide a measurable path from the TeX concept to
implementation:

- Milestone 15: real split/merge event in Chrono dynamics;
- Milestone 16: reusable runtime state machine;
- Milestone 17: objective tangential transport regression;
- Milestone 18: bidirectional symmetric pair conservation;
- Milestone 19: topology consistency metrics for persistence, count variation,
  normal stability, and response continuity.

The next step is to use these metrics as the comparison surface for feature
switching and non-convex benchmark scenes.
