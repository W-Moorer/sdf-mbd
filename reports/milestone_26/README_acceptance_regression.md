# Milestone 26: Acceptance-Level Regression

This milestone freezes the key quantitative claims used by the TeX paper into
CI-checkable regression metrics.

## Script

```text
scripts/check_field_contact_acceptance.py
```

Default behavior:

1. run Milestone 20 feature-switching demo;
2. run Milestone 21/22 multi-slip/interlock demo;
3. run Milestone 23 Chrono native SMC baseline;
4. regenerate native-vs-field comparison CSV;
5. assert acceptance thresholds.

Fast CSV-only mode:

```powershell
python scripts\check_field_contact_acceptance.py --check-only
```

Full mode:

```powershell
python scripts\check_field_contact_acceptance.py
```

## CTest Integration

When `CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS=ON`, CMake now registers:

```text
demo_CH_field_contact_acceptance_regression
```

Run:

```powershell
ctest --test-dir build -C Release -R "demo_CH_field_contact_(regression|acceptance_regression)" --output-on-failure
```

## Fixed Acceptance Metrics

Current thresholds:

```text
direct_inherit_energy_threshold = 1.5
normal_bin_normal_jump_threshold = 1.5
native_oscillation_threshold = 1.03
tolerance = 1e-6
```

Checks:

- `narrow_groove_entrance/direct_inherit` must have
  `max_inherited_energy_ratio > 1.5`.
- `narrow_groove_entrance/minimal_rotation_gate_split_merge` must have
  `max_inherited_energy_ratio <= 1 + tol`.
- `normal_bin_fragments` must have
  `rms_normal_jump_ratio > 1.5` for the Milestone 21 comparison rows.
- Chrono native SMC must have at least one force or torque oscillation ratio
  above `1.03`.
- All variants in Milestones 20 and 21 must satisfy
  `max_tangential_force_ratio <= 1 + tol`.
- If Milestone 24 mesh-SDF summary is present, it is also checked for the same
  Coulomb bound.

## Current Result

The current full run passes:

```text
direct_inherit max_inherited_energy_ratio = 1.99992422
minimal_rotation_gate_split_merge max_inherited_energy_ratio = 1.00000000
min normal_bin_fragments rms_normal_jump_ratio = 1.70497080
max native_to_field force/torque oscillation ratio = 1.06117478
max Coulomb ratio = 1.00000000
```

Output summary:

```text
out/milestone_26/acceptance_regression_summary.json
```

## Verification

Validated locally with:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_feature_switching demo_CH_field_contact_multislip_interlock demo_CH_chrono_native_contact_baseline
ctest --test-dir build -C Release -R "demo_CH_field_contact_(regression|acceptance_regression)" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 2
```
