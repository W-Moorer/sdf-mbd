# Milestone 14: Field Contact Primitive Core

This milestone moves the field-based contact primitive idea out of a single demo
and into a reusable core header:

`src/chrono/collision/ChFieldContactPrimitives.h`

## Implemented

- Surface graph samples with area weights and explicit adjacency.
- Connected active-set extraction into field contact primitives.
- Area-weighted primitive descriptors: center, normal, tangent basis, mean phi,
  area, and representative velocity.
- Sample-level normal contact integration:
  `dF_i = (k_n max(-phi_i, 0) + c_n max(-v_n, 0)) n_i a_i`.
- Primitive force/torque aggregation from the sample integral.
- Area-Jaccard overlap and geometry fallback for persistent primitive history
  sources.
- Minimal-rotation objective transport for elastic tangential history.
- Non-amplifying history gate.
- Single Coulomb projection for tangential force with elastic state update.

## Smoke Test

The new OpenVDB-backed executable is:

`demo_CH_field_contact_primitives_openvdb`

It runs a kinematic rolling/sliding sphere surface graph against a sparse
OpenVDB sphere SDF and writes:

- `out/milestone_14/field_contact_primitives_core.csv`
- `out/milestone_14/field_contact_primitives_core_summary.csv`

Current summary:

- frames: 420
- surface samples: 4608
- total patch frames: 420
- newborn primitives: 1
- max tangential force ratio: 1.0
- minimum nonzero overlap gate: about 0.974

The force ratio result verifies that the tangential update respects the Coulomb
disk in this smoke test. The overlap gate result verifies persistent history
reuse through the new primitive-source interface.

## Not Yet Complete

- Mesh-to-surface-graph construction from arbitrary triangle meshes.
- Mesh-to-OpenVDB SDF preprocessing and asset loading.
- Bidirectional A-surface-to-B-field and B-surface-to-A-field primitive
  extraction.
- Full split/merge history aggregation in the Chrono dynamics loop.
- Integration into Chrono's collision/contact pipeline beyond explicit external
  force accumulation demos.
- Regression tests for objectivity, non-amplification, Coulomb feasibility, and
  split/merge cases.
