# Milestone 17: Tangential Objectivity Regression

This milestone promotes the TeX tangential-objectivity requirements into
standalone CI-checkable regression metrics.

## Regression Target

Executable:

`demo_CH_field_contact_regression`

CTest entry:

`demo_CH_field_contact_regression`

Run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
ctest --test-dir build -C Release -R demo_CH_field_contact_regression --output-on-failure
```

## Fixed Metrics

### Global rotation covariance

The same tangential update is evaluated in two frames related by a fixed global
rotation `Q`. The regression requires the transported history, gated history,
trial state, final force, and updated elastic history to rotate covariantly:

```text
x_rotated == Q * x_base
F_rotated == Q * F_base
```

Scalar quantities such as force norm, stored energy, and stick/slip state must
remain unchanged. The current tolerance is `1e-12`.

### Rotating normal without slip

A contact normal is slowly rotated while the elastic tangential history and the
relative tangential velocity are both zero. The update must not create
tangential force or nonzero elastic history from frame motion alone.

Current thresholds:

- max tangential force <= `1e-12`;
- max elastic history norm <= `1e-15`.

### Minimal-rotation transport ablation

For a controlled 30 degree normal rotation, the regression compares:

- minimal-rotation transport;
- projection-only transport.

The minimal-rotation result must match the rigidly rotated reference to
relative error <= `1e-12`. The projection-only transport is intentionally kept
as an ablation baseline and must show a visible frame-rotation error:

- projection error ratio >= `0.10`;
- projection energy ratio < `0.90`;
- projection error is at least `1e8` times the guarded minimal-rotation error.

This fixes the TeX claim that direct projection is only a low-order fallback,
not the final objective transport rule.

## Existing Metrics Preserved

The same executable still checks:

- Coulomb disk feasibility;
- inherited history non-amplification;
- synthetic split/merge source classification;
- runtime tracker split/merge classification.

Together, Milestones 15-17 now cover the minimal theoretical contract for the
tangential primitive layer: tangency, objectivity, Coulomb feasibility, history
non-amplification, and explicit topology-event handling.
