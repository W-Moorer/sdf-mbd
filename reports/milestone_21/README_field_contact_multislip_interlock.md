# Milestone 21: Multi-Point Sliding and Interlock Benchmarks

This milestone adds the first complex non-convex benchmark suite for the TeX
claim that field primitives reduce contact-set churn, normal jumps, and response
oscillation in multi-contact situations.

## Demo

Target:

`src/demos/core/demo_CH_field_contact_multislip_interlock.cpp`

Build and run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_multislip_interlock
.\build\bin\Release\demo_CH_field_contact_multislip_interlock.exe
```

Outputs:

```text
out/milestone_21/field_contact_multislip_interlock_frames.csv
out/milestone_21/field_contact_multislip_interlock_patches.csv
out/milestone_21/field_contact_multislip_interlock_summary.csv
out/milestone_21/field_contact_multislip_interlock_comparison.csv
```

## Scenes

The benchmark remains analytic to keep geometry deterministic while increasing
non-convexity and simultaneous contact complexity.

- `guide_rail_sliding`: a rolling/sliding sphere contacts two raised rails and
  a shallow groove, producing persistent one- and two-patch contact.
- `nested_interlock`: a sphere passes through a slot with collar lips and inner
  shoulders, representing a nested/interlocking transition.
- `multi_patch_rolling_sliding`: a rolling/sliding sphere crosses staggered
  raised pads, producing persistent one- and two-patch contact.

## Compared Variants

- `field_primitive`: connected field-contact primitive extraction plus the
  runtime tracker, split/merge handling, minimal-rotation tangential transport,
  and Coulomb return.
- `raw_sample_contacts`: every active surface sample becomes its own contact
  primitive. This is the contact-set churn baseline.
- `normal_bin_fragments`: active samples are fragmented by quantized local
  normals, and those quantized normals are used in force integration. This is
  the feature/normal-switching baseline.

## Current Run

Field primitive summary:

```text
guide_rail_sliding,field_primitive:
  active_frames=619
  mean_patch_count=1.07845934
  max_patch_count=2
  mean_abs_patch_count_change=0.01714286
  rms_normal_jump_angle=0.02997340
  force_oscillation_index=0.01130530
  torque_oscillation_index=0.01353610

nested_interlock,field_primitive:
  active_frames=152
  mean_patch_count=0.20936639
  max_patch_count=1
  mean_abs_patch_count_change=0.00275862
  rms_normal_jump_angle=0.02669830
  force_oscillation_index=0.00715738
  torque_oscillation_index=0.00831260

multi_patch_rolling_sliding,field_primitive:
  active_frames=740
  mean_patch_count=1.29027963
  max_patch_count=2
  total_merge_patches=18
  mean_abs_patch_count_change=0.02933333
  rms_normal_jump_angle=0.04002010
  force_oscillation_index=0.01510170
  torque_oscillation_index=0.03477840
```

Comparison ratios are reported as baseline divided by field primitive. Values
above 1 mean the baseline is worse.

```text
guide_rail_sliding,raw_sample_contacts:
  mean_patch_count_change_ratio=17.00000000
  rms_normal_jump_ratio=1.41690613
  topology_event_ratio=476.50000000

guide_rail_sliding,normal_bin_fragments:
  mean_patch_count_change_ratio=11.50000000
  rms_normal_jump_ratio=1.80768673
  force_oscillation_ratio=1.04552735
  torque_oscillation_ratio=1.13602314

nested_interlock,raw_sample_contacts:
  mean_patch_count_change_ratio=21.00000000
  rms_normal_jump_ratio=1.09196602
  topology_event_ratio=1044.50000000

nested_interlock,normal_bin_fragments:
  mean_patch_count_change_ratio=14.00000000
  rms_normal_jump_ratio=1.70497080
  force_oscillation_ratio=1.08734630
  torque_oscillation_ratio=1.38616225

multi_patch_rolling_sliding,raw_sample_contacts:
  mean_patch_count_change_ratio=8.27272727
  rms_normal_jump_ratio=1.47989309
  topology_event_ratio=171.70000000

multi_patch_rolling_sliding,normal_bin_fragments:
  mean_patch_count_change_ratio=8.04545455
  rms_normal_jump_ratio=1.97031896
  force_oscillation_ratio=1.34290823
  torque_oscillation_ratio=1.73840578
```

All variants keep:

```text
max_tangential_force_ratio=1.00000000
```

so the comparison does not rely on violating the Coulomb bound.

## Interpretation

The raw-sample baseline demonstrates contact-set instability: the same physical
interaction becomes thousands of newborn/death/merge/split bookkeeping events,
with patch-count-change ratios from 8.27x to 21.0x.

The normal-bin fragmented baseline demonstrates feature-normal instability:
RMS normal jump is 1.70x to 1.97x larger in the interlock and multi-patch cases,
and total force/torque oscillation is also larger:

- interlock: force 1.09x, torque 1.39x;
- multi-patch rolling/sliding: force 1.34x, torque 1.74x.

This is the first benchmark where the TeX claim is visible in the intended
complex setting: field primitives keep the contact collection compact and
persistent while preserving Coulomb feasibility and reducing normal/response
oscillation relative to fragmented feature contacts.
