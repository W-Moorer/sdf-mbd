# Milestone 14: Field Contact Primitive Core

This milestone moves the field-based contact primitive idea out of a single demo
and into a reusable core header:

`src/chrono/collision/ChFieldContactPrimitives.h`

## Implemented

- Surface graph samples with area weights and explicit adjacency.
- Triangle-mesh surface graph construction with barycentric vertex area weights.
- Connected active-set extraction into field contact primitives.
- Area-weighted primitive descriptors: center, normal, tangent basis, mean phi,
  area, and representative velocity.
- Sample-level normal contact integration:
  `dF_i = (k_n max(-phi_i, 0) + c_n max(-v_n, 0)) n_i a_i`.
- Primitive force/torque aggregation from the sample integral.
- Area-Jaccard overlap and geometry fallback for persistent primitive history
  sources.
- Multi-source tangential history aggregation for merge-like inheritance, with
  weights normalized so the inherited history is not amplified by topology
  changes.
- Minimal-rotation objective transport for elastic tangential history.
- Non-amplifying history gate.
- Single Coulomb projection for tangential force with elastic state update.

## Smoke Test

The new OpenVDB-backed executable is:

`demo_CH_field_contact_primitives_openvdb`

It runs a kinematic rolling/sliding triangle-mesh sphere surface graph against a
triangle-mesh OpenVDB SDF. The same run also performs reverse extraction
(`target surface -> moving-body SDF`) to exercise the bidirectional primitive
path. It writes:

- `out/milestone_14/field_contact_primitives_core.csv`
- `out/milestone_14/field_contact_primitives_core_summary.csv`

Current summary:

- frames: 420
- target mesh vertices: 8066
- target mesh faces: 16128
- moving mesh vertices: 4514
- moving mesh faces: 9024
- target SDF active voxels: about 999k
- moving SDF active voxels: about 63k
- total patch frames: 420
- total reverse patch frames: 420
- newborn primitives: 1
- max tangential force ratio: 1.0
- minimum nonzero overlap gate: about 0.969
- synthetic merge source count: 2
- synthetic merge weight sum: 1.0

The force ratio result verifies that the tangential update respects the Coulomb
disk in this smoke test. The overlap gate result verifies persistent history
reuse through the new primitive-source interface. The synthetic merge result
verifies that two previous primitives can contribute to one inherited tangential
history without the source weights exceeding one.

## Not Yet Complete

- Asset loading from external mesh files.
- Strong primitive matching and de-duplication between the two bidirectional
  primitive sets.
- Dynamic split/merge scenes that naturally exercise multi-source history in
  the main simulation loop.
- Integration into Chrono's collision/contact pipeline beyond explicit external
  force accumulation demos.
- Regression tests for objectivity, non-amplification, Coulomb feasibility, and
  split/merge cases.
