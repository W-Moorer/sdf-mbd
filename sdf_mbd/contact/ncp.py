"""Nonlinear complementarity problem utilities."""

from __future__ import annotations

import numpy as np


def _maybe_scalar(value: np.ndarray) -> float | np.ndarray:
    array = np.asarray(value)
    if array.ndim == 0:
        return float(array)
    return value


def smooth_fischer_burmeister(g: float | np.ndarray, lam: float | np.ndarray, eps: float) -> float | np.ndarray:
    """Evaluate the smoothed Fischer-Burmeister residual.

    For the Signorini condition `0 <= g perpendicular lambda >= 0`, the smoothed
    residual is

        Phi_eps(g, lambda) = sqrt(g^2 + lambda^2 + eps^2) - g - lambda.

    Scalars and NumPy-broadcastable arrays are supported.
    """

    g_arr = np.asarray(g, dtype=float)
    lam_arr = np.asarray(lam, dtype=float)
    eps_value = float(eps)
    root = np.sqrt(g_arr * g_arr + lam_arr * lam_arr + eps_value * eps_value)
    return _maybe_scalar(root - g_arr - lam_arr)


def smooth_fischer_burmeister_grad(
    g: float | np.ndarray,
    lam: float | np.ndarray,
    eps: float,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    """Return derivatives of the smoothed Fischer-Burmeister residual.

    The derivatives are

        dPhi/dg = g / sqrt(g^2 + lambda^2 + eps^2) - 1
        dPhi/dlambda = lambda / sqrt(g^2 + lambda^2 + eps^2) - 1.
    """

    g_arr = np.asarray(g, dtype=float)
    lam_arr = np.asarray(lam, dtype=float)
    eps_value = float(eps)
    denom = np.sqrt(g_arr * g_arr + lam_arr * lam_arr + eps_value * eps_value)
    denom = np.maximum(denom, np.finfo(float).tiny)
    d_g = g_arr / denom - 1.0
    d_lam = lam_arr / denom - 1.0
    return _maybe_scalar(d_g), _maybe_scalar(d_lam)


def complementarity_error(g: float | np.ndarray, lam: float | np.ndarray) -> float:
    """Return a scalar complementarity violation measure.

    For one contact, the diagnostic is

        max(0, -g) + max(0, -lambda) + abs(g * lambda).

    For vector inputs this function sums the per-contact values.
    """

    g_arr = np.asarray(g, dtype=float)
    lam_arr = np.asarray(lam, dtype=float)
    value = np.maximum(0.0, -g_arr) + np.maximum(0.0, -lam_arr) + np.abs(g_arr * lam_arr)
    return float(np.sum(value))

