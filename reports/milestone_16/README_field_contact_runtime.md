# Milestone 16: Field Contact Runtime API

This milestone moves the field-contact primitive workflow out of the
split/merge demo and into a reusable runtime layer.

## New API

Header:

`src/chrono/collision/ChFieldContactRuntime.h`

Main types:

- `FieldContactRuntimeSettings`: groups extraction, normal contact, tangential
  contact, and history inheritance settings.
- `FieldContactPrimitiveTracker`: owns persistent primitive IDs, previous
  primitive snapshots, tangential history storage, split/merge classification,
  and per-frame state advancement.
- `FieldContactStepResult`: returns per-patch runtime results, frame-level
  event metrics, total contact force, and total contact torque.
- `FieldContactPatchRuntimeResult`: stores the classified primitive, history
  sources, tangential update, Coulomb ratio, inherited-energy ratio, and final
  force contribution.

The runtime API is SDF-backend independent. Callers still provide a
`SurfaceGraph` and a `FieldSampleQuery` array. OpenVDB, NanoVDB, analytic SDFs,
or cached query backends can all feed the same tracker.

## Demo Refactor

`demo_CH_field_contact_split_merge_openvdb` now uses:

```cpp
FieldContactStepResult contact_step =
    contact_tracker.Evaluate(plane_graph, queries, body_pos, contact_settings);
```

The demo keeps only scene generation, OpenVDB SDF sampling, Chrono body
integration, and CSV reporting. The following logic now lives in the runtime
layer:

- active-set extraction;
- primitive extraction;
- persistent ID assignment;
- stable/newborn/split/merge/death classification;
- source-history aggregation;
- normal and tangential primitive force evaluation;
- primitive force/torque summation.

## Regression Coverage

`demo_CH_field_contact_regression` now checks both low-level primitive
invariants and the reusable runtime tracker:

- Coulomb disk feasibility;
- inherited tangential history non-amplification;
- synthetic merge source classification;
- synthetic split source classification;
- runtime tracker merge/split event classification.

Run:

```powershell
cmake --build build --config Release --target demo_CH_field_contact_regression
ctest --test-dir build -C Release -R demo_CH_field_contact_regression --output-on-failure
```

## Current Scope

The runtime API is still single-sided: one body surface graph queries the other
body field. The next engineering step is bidirectional pair evaluation with
deduplication and equal-opposite force application to both bodies.
