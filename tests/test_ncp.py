from __future__ import annotations

import numpy as np

from sdf_mbd.contact.ncp import (
    complementarity_error,
    smooth_fischer_burmeister,
    smooth_fischer_burmeister_grad,
)


def test_smooth_fischer_burmeister_open_contact() -> None:
    value = smooth_fischer_burmeister(0.2, 0.0, 1.0e-6)

    assert abs(value) < 1.0e-9


def test_smooth_fischer_burmeister_closed_contact() -> None:
    value = smooth_fischer_burmeister(0.0, 5.0, 1.0e-6)

    assert abs(value) < 1.0e-9


def test_smooth_fischer_burmeister_gradient_matches_finite_difference() -> None:
    g = 0.13
    lam = 0.27
    eps = 1.0e-4
    h = 1.0e-6

    d_g, d_lam = smooth_fischer_burmeister_grad(g, lam, eps)
    fd_g = (
        smooth_fischer_burmeister(g + h, lam, eps)
        - smooth_fischer_burmeister(g - h, lam, eps)
    ) / (2.0 * h)
    fd_lam = (
        smooth_fischer_burmeister(g, lam + h, eps)
        - smooth_fischer_burmeister(g, lam - h, eps)
    ) / (2.0 * h)

    assert np.isclose(d_g, fd_g, rtol=1.0e-6, atol=1.0e-8)
    assert np.isclose(d_lam, fd_lam, rtol=1.0e-6, atol=1.0e-8)


def test_complementarity_error_vector_inputs() -> None:
    error = complementarity_error(np.array([0.1, -0.2]), np.array([0.0, 0.3]))

    assert np.isclose(error, 0.26)

