# Milestone 18: Bidirectional Symmetric Pair Contact

This milestone implements the TeX "symmetric treatment" step for field-contact
primitives.

## Runtime API

Header:

`src/chrono/collision/ChFieldContactRuntime.h`

New types:

- `FieldContactBodyWrench`: force plus torque about a body reference point.
- `FieldContactPairSettings`: tolerances used by the pair combiner.
- `FieldContactPairResult`: symmetric pair wrenches, residual diagnostics, raw
  one-sided force magnitudes, deduplication ratio, and max Coulomb ratio.

New function:

```cpp
FieldContactPairResult CombineBidirectionalFieldContactPair(
    const FieldContactStepResult& a_surface_on_b,
    const ChVector3d& reference_a,
    const FieldContactStepResult& b_surface_on_a,
    const ChVector3d& reference_b,
    const FieldContactPairSettings& settings = FieldContactPairSettings());
```

The two `FieldContactStepResult` inputs are produced by separate one-sided
tracker evaluations:

- A surface samples B field;
- B surface samples A field.

The combiner treats the two one-sided forces as duplicate estimates of the same
physical pair interaction. If both directions are active, it uses the averaged
anti-symmetric estimate:

```text
F_A = 0.5 * (F_A_raw - F_B_raw)
F_B = -F_A
```

This prevents bidirectional sampling from doubling the normal response.

## ChBody Accumulation Pattern

The returned torques are expressed about each body reference point. A caller can
apply the pair to two Chrono bodies with the same accumulator pattern used by
the current demos:

```cpp
body_a->AccumulateForce(acc_a, pair.on_a.force, reference_a, false);
body_a->AccumulateTorque(acc_a, pair.on_a.torque, false);
body_b->AccumulateForce(acc_b, pair.on_b.force, reference_b, false);
body_b->AccumulateTorque(acc_b, pair.on_b.torque, false);
```

The pair result also reports force and torque residuals evaluated about the pair
mid-reference.

## Regression Checks

`demo_CH_field_contact_regression` now fixes the Milestone 18 metrics:

- A surface to B field and B surface to A field both produce one-sided runtime
  primitive responses.
- The symmetric pair force residual is zero within `1e-12`.
- The symmetric pair torque residual is zero within `1e-12`.
- Bidirectional sampling does not amplify duplicated normal response:
  `force_amplification_ratio <= 1 + 1e-12`.
- The pair max Coulomb ratio remains `<= 1 + 1e-12`.
- The returned pair wrenches can be accumulated on two `ChBody` instances while
  preserving the same force and torque residual bounds.

Run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
ctest --test-dir build -C Release -R demo_CH_field_contact_regression --output-on-failure
```

## Current Scope

The symmetric combiner is backend-independent and works with OpenVDB, NanoVDB,
analytic SDFs, or cached field-query results. The current regression uses a
small analytic query fixture to keep CI deterministic and independent of
OpenVDB.

The next step is to convert the OpenVDB split/merge demo into a true two-body
pair demo where both bodies own surface graphs and SDF query backends.
