"""SDF contact geometry utilities."""

from __future__ import annotations

import numpy as np

from sdf_mbd.sdf.primitives import grad_phi, phi


def _safe_normalize(v: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    vec = np.asarray(v, dtype=float)
    if fallback is None:
        fallback = np.zeros_like(vec, dtype=float)
        fallback[-1] = 1.0
    norm = np.linalg.norm(vec)
    if norm <= 1.0e-14:
        return np.asarray(fallback, dtype=float).copy()
    return vec / norm


def sdf_gap(sdf: object, x: np.ndarray | list[float] | tuple[float, ...]) -> float:
    """Return the SDF contact gap `g = phi(x)`."""

    return float(phi(sdf, x))


def sdf_normal(sdf: object, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    """Return the normalized SDF normal `n = grad phi(x) / ||grad phi(x)||`."""

    g = grad_phi(sdf, x)
    fallback = np.zeros_like(g, dtype=float)
    fallback[-1] = 1.0
    return _safe_normalize(g, fallback=fallback)


def point_contact_jacobian(sdf: object, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    """Return the point-mass contact Jacobian.

    For a point mass with generalized coordinate `q = x`,

        J_phi(q) = grad phi(x)^T.

    This prototype returns the row as a one-dimensional array with shape
    `(dim,)`. Callers that need matrix form can use `J[None, :]`.
    """

    return np.asarray(grad_phi(sdf, x), dtype=float)


def finite_difference_gap_gradient(
    sdf: object,
    x: np.ndarray | list[float] | tuple[float, ...],
    h: float = 1.0e-6,
) -> np.ndarray:
    """Central finite-difference approximation of `d phi(x) / d x`."""

    point = np.asarray(x, dtype=float)
    out = np.zeros_like(point, dtype=float)
    for i in range(point.size):
        step = np.zeros_like(point, dtype=float)
        step[i] = h
        out[i] = (sdf_gap(sdf, point + step) - sdf_gap(sdf, point - step)) / (2.0 * h)
    return out

