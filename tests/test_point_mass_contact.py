from __future__ import annotations

import numpy as np

from sdf_mbd.sdf import PlaneSDF
from sdf_mbd.solvers import simulate_ncp_point_mass, simulate_penalty_point_mass


def test_point_mass_plane_penalty_runs() -> None:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    trajectory = simulate_penalty_point_mass(
        ground,
        q0=[0.0, 0.01],
        v0=[0.0, -1.0],
        mass=1.0,
        gravity=9.81,
        dt=1.0e-3,
        t_end=0.02,
        k_n=1.0e4,
    )

    assert len(trajectory) == 21
    assert all(np.all(np.isfinite(state.q)) and np.all(np.isfinite(state.v)) for _, state, _ in trajectory)


def test_point_mass_plane_sdf_ncp_runs() -> None:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    trajectory = simulate_ncp_point_mass(
        ground,
        q0=[0.0, 0.01],
        v0=[0.0, -1.0],
        mass=1.0,
        gravity=9.81,
        dt=1.0e-3,
        t_end=0.02,
        eps=1.0e-6,
    )

    assert len(trajectory) == 21
    assert all(diag.solver_success for _, _, diag in trajectory)
    assert all(np.all(np.isfinite(state.q)) and np.all(np.isfinite(state.v)) for _, state, _ in trajectory)


def test_sdf_ncp_has_small_penetration_in_basic_case() -> None:
    ground = PlaneSDF([0.0, 1.0], 0.0)
    trajectory = simulate_ncp_point_mass(
        ground,
        q0=[0.0, 0.01],
        v0=[0.0, -1.0],
        mass=1.0,
        gravity=9.81,
        dt=1.0e-3,
        t_end=0.05,
        eps=1.0e-6,
    )
    max_penetration = max(diag.penetration for _, _, diag in trajectory)
    max_ncp_residual = max(diag.ncp_residual for _, _, diag in trajectory)

    assert max_penetration < 1.0e-8
    assert max_ncp_residual < 1.0e-8

