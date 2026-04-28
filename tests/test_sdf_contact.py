from __future__ import annotations

import numpy as np

from sdf_mbd.contact import (
    finite_difference_gap_gradient,
    point_contact_jacobian,
    sdf_gap,
    sdf_normal,
)
from sdf_mbd.sdf import CircleSDF, PlaneSDF


def test_sdf_gap_ground_plane() -> None:
    ground = PlaneSDF([0.0, 1.0], 0.0)

    assert np.isclose(sdf_gap(ground, [0.0, 0.2]), 0.2)
    assert np.isclose(sdf_gap(ground, [0.0, -0.1]), -0.1)


def test_plane_sdf_gradient_equals_plane_normal() -> None:
    ground = PlaneSDF([0.0, 2.0], 0.0)

    assert np.allclose(sdf_normal(ground, [0.3, 0.7]), np.array([0.0, 1.0]))


def test_circle_sdf_gradient_has_unit_norm_away_from_center() -> None:
    circle = CircleSDF([0.0, 0.0], 0.5)
    normal = sdf_normal(circle, [0.3, 0.4])

    assert np.isclose(np.linalg.norm(normal), 1.0)
    assert np.allclose(normal, np.array([0.6, 0.8]))


def test_point_contact_jacobian_matches_finite_difference() -> None:
    circle = CircleSDF([0.1, -0.2], 0.75)
    x = np.array([0.7, 0.5])

    jac = point_contact_jacobian(circle, x)
    fd = finite_difference_gap_gradient(circle, x)

    assert jac.shape == (2,)
    assert np.allclose(jac, fd, rtol=1.0e-6, atol=1.0e-8)

