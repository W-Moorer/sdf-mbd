# Milestone 23: Chrono Native SMC Contact Baseline

This milestone adds a real Chrono native collision/contact baseline for the
field-contact primitive paper experiments.

## What Was Added

Native Chrono benchmark:

`src/demos/core/demo_CH_chrono_native_contact_baseline.cpp`

Comparison and figure script:

`scripts/compare_chrono_native_field.py`

The native benchmark uses Chrono SMC contact with Bullet collision and static
collision bodies built from native Chrono primitives. It covers the same three
benchmark families used in Milestone 21:

- `guide_rail_sliding`: native cylinder over rail/groove boxes.
- `nested_interlock`: native sphere through slot/interlock boxes.
- `multi_patch_rolling_sliding`: native sphere over staggered native box pads.

The moving body is a high-mass kinematic probe driven along the same deterministic
path as the field primitive benchmark. Contact count, force, and torque are
computed by Chrono native SMC. Commanded and post-step actual displacement,
velocity, and acceleration are recorded to check physical consistency of the
trajectory used for comparison.

## Commands

Run native benchmark:

```powershell
cmake --build build --config Release --target demo_CH_chrono_native_contact_baseline
.\build\bin\Release\demo_CH_chrono_native_contact_baseline.exe
```

Align against Milestone 21 field primitive output:

```powershell
python scripts\compare_chrono_native_field.py `
  --native out\milestone_23\chrono_native_baseline_frames.csv `
  --field out\milestone_21\field_contact_multislip_interlock_frames.csv `
  --output out\milestone_23
```

Outputs:

```text
out/milestone_23/chrono_native_baseline_frames.csv
out/milestone_23/chrono_native_field_comparison.csv
out/milestone_23/figures/index.html
```

The script currently generates 18 SVG figures:

- force jump, field vs native;
- torque jump, field vs native;
- patch/contact count timeline;
- displacement norm;
- velocity norm;
- acceleration norm.

## Current Results

Native trajectory tracking remains small, which means the comparison is using a
well-controlled shared kinematic path:

```text
guide_rail_sliding:
  max_position_error=0.00149652
  max_velocity_error=0.00000054
  max_acceleration_error=0.00723141

nested_interlock:
  max_position_error=0.00123345
  max_velocity_error=0.00000100
  max_acceleration_error=0.00823701

multi_patch_rolling_sliding:
  max_position_error=0.00143356
  max_velocity_error=0.00000043
  max_acceleration_error=0.02296619
```

Contact-count churn is consistently larger in the native baseline:

```text
native_to_field_count_change_ratio:
  guide_rail_sliding=3.50000000
  nested_interlock=9.00000000
  multi_patch_rolling_sliding=1.27272727
```

Force/torque response is more nuanced because the native collision geometry is a
piecewise box/cylinder surrogate while the field primitive benchmark uses an
analytic height field. With the current material calibration:

```text
guide_rail_sliding:
  native_to_field_force_oscillation_ratio=0.60401603
  native_to_field_torque_oscillation_ratio=1.06117478
  native_to_field_rms_force_jump_ratio=2.25357035
  native_to_field_rms_torque_jump_ratio=2.03355838

nested_interlock:
  native_to_field_force_oscillation_ratio=0.91927303
  native_to_field_torque_oscillation_ratio=1.03955227
  native_to_field_rms_force_jump_ratio=1.16843296
  native_to_field_rms_torque_jump_ratio=0.90961792

multi_patch_rolling_sliding:
  native_to_field_force_oscillation_ratio=0.70849892
  native_to_field_torque_oscillation_ratio=0.31715374
  native_to_field_rms_force_jump_ratio=0.91023836
  native_to_field_rms_torque_jump_ratio=0.85947233
```

The strongest current conclusion from the real native baseline is therefore
about contact collection stability:

```text
field primitives produce fewer contact-count changes than Chrono native SMC
point/manifold contacts under the same controlled motion.
```

The force/torque comparison is useful but not yet final because native geometric
surrogates and analytic field geometry are not identical. It should be tightened
in a later free-dynamics or matched-mesh SDF/native experiment.

## Meaning For The Paper

Milestone 22 used a deterministic `chrono_traditional_proxy`. Milestone 23 now
adds real Chrono native SMC collision/contact data, which improves baseline
credibility.

For the TeX paper, this should be described carefully:

- the native baseline is real Chrono native SMC, not a proxy;
- the current run is a controlled kinematic-probe comparison, not a free dynamic
  trajectory comparison;
- displacement, velocity, and acceleration curves are generated to show that the
  shared trajectory is physically smooth and accurately tracked;
- contact-count timelines already support the paper claim that field primitives
  reduce contact collection churn;
- force/torque conclusions should be reported with geometry-calibration caveats.

## Next Step

The next validation step should be a matched-geometry experiment:

1. Build the native collision geometry and the field SDF from the same triangle
   mesh.
2. Run both under either the same kinematic probe or the same external drive.
3. Compare displacement, velocity, acceleration, contact count, force jump, and
   torque jump on a truly matched geometry.

That is the point where force/torque superiority can be claimed more strongly.
