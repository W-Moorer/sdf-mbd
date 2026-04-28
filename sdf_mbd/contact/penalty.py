"""Penalty-contact baseline assembled from SDF gap and normal."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ncp import complementarity_error
from .sdf_contact import sdf_gap, sdf_normal


@dataclass(frozen=True)
class PenaltyContactResult:
    """Diagnostics for a frictionless normal penalty contact."""

    gap: float
    normal: np.ndarray
    lambda_n: float
    force: np.ndarray
    penetration: float
    complementarity_error: float


def penalty_normal_force_from_sdf(
    sdf: object,
    x: np.ndarray | list[float] | tuple[float, ...],
    k_n: float,
    g_dot: float | None = None,
    c_n: float = 0.0,
) -> PenaltyContactResult:
    """Return normal penalty force from an SDF contact.

    The undamped model is

        lambda_n = k_n * max(0, -g),  force = lambda_n * n.

    If `g_dot` and `c_n` are supplied, a normal damping contribution
    `c_n * max(0, -g_dot)` is added and the multiplier is clamped nonnegative.
    """

    if k_n < 0.0:
        raise ValueError("Normal penalty stiffness k_n must be nonnegative.")
    if c_n < 0.0:
        raise ValueError("Normal damping c_n must be nonnegative.")

    point = np.asarray(x, dtype=float)
    gap = sdf_gap(sdf, point)
    normal = sdf_normal(sdf, point)
    penetration = max(0.0, -gap)
    lambda_n = float(k_n) * penetration
    if g_dot is not None and c_n > 0.0:
        lambda_n += float(c_n) * max(0.0, -float(g_dot))
    lambda_n = max(0.0, lambda_n)
    force = lambda_n * normal
    error = complementarity_error(gap, lambda_n)
    return PenaltyContactResult(
        gap=gap,
        normal=normal,
        lambda_n=lambda_n,
        force=force,
        penetration=penetration,
        complementarity_error=error,
    )

