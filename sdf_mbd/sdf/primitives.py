"""Analytical signed-distance field primitives.

Each primitive exposes `phi(x)` and `grad(x)`. The adapter functions at the end
of this module also accept common alternate method names used by SDF backends.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _as_array(x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _is_single_point(x: np.ndarray) -> bool:
    return x.ndim == 1


def _maybe_scalar(value: np.ndarray) -> float | np.ndarray:
    array = np.asarray(value)
    if array.ndim == 0:
        return float(array)
    return value


def _safe_normalize(v: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    vec = np.asarray(v, dtype=float)
    if fallback is None:
        fallback = np.zeros(vec.shape[-1], dtype=float)
        fallback[-1] = 1.0
    norm = np.linalg.norm(vec, axis=-1, keepdims=True)
    safe = np.where(norm > 1.0e-14, norm, 1.0)
    out = vec / safe
    if np.any(norm <= 1.0e-14):
        out = np.where(norm > 1.0e-14, out, fallback)
    return out


@dataclass(frozen=True)
class PlaneSDF:
    """Signed distance to a plane.

    The plane is represented by `normal dot x - offset = 0`, so
    `phi(x) = normal dot x - offset`. The normal is normalized on construction.
    A 2D ground plane y = 0 is `PlaneSDF([0, 1], 0)`.
    """

    normal: np.ndarray
    offset: float = 0.0

    def __init__(self, normal: np.ndarray | list[float] | tuple[float, ...], offset: float = 0.0):
        n = _as_array(normal)
        norm = np.linalg.norm(n)
        if norm <= 1.0e-14:
            raise ValueError("PlaneSDF normal must be nonzero.")
        object.__setattr__(self, "normal", n / norm)
        object.__setattr__(self, "offset", float(offset))

    def phi(self, x: np.ndarray | list[float] | tuple[float, ...]) -> float | np.ndarray:
        points = _as_array(x)
        value = np.tensordot(points, self.normal, axes=([-1], [0])) - self.offset
        return _maybe_scalar(value)

    def grad(self, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
        points = _as_array(x)
        if _is_single_point(points):
            return self.normal.copy()
        return np.broadcast_to(self.normal, points.shape).copy()

    signed_distance = phi
    gradient = grad


@dataclass(frozen=True)
class CircleSDF:
    """Signed distance to a circle or sphere.

    The dimensionality is determined by `center`. In 2D this is a circle; in 3D
    it is a sphere.
    """

    center: np.ndarray
    radius: float

    def __init__(self, center: np.ndarray | list[float] | tuple[float, ...], radius: float):
        if radius <= 0.0:
            raise ValueError("CircleSDF radius must be positive.")
        object.__setattr__(self, "center", _as_array(center))
        object.__setattr__(self, "radius", float(radius))

    def phi(self, x: np.ndarray | list[float] | tuple[float, ...]) -> float | np.ndarray:
        points = _as_array(x)
        delta = points - self.center
        value = np.linalg.norm(delta, axis=-1) - self.radius
        return _maybe_scalar(value)

    def grad(self, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
        points = _as_array(x)
        delta = points - self.center
        fallback = np.zeros_like(self.center, dtype=float)
        fallback[0] = 1.0
        return _safe_normalize(delta, fallback=fallback)

    signed_distance = phi
    gradient = grad


SphereSDF = CircleSDF


@dataclass(frozen=True)
class AxisAlignedBoxSDF:
    """Signed distance to an axis-aligned box.

    The box is defined by its center and positive half extents. The gradient is
    evaluated by central finite difference, which is robust for smooth points
    away from edges and corners and avoids over-specializing this prototype.
    """

    center: np.ndarray
    half_extents: np.ndarray

    def __init__(
        self,
        center: np.ndarray | list[float] | tuple[float, ...],
        half_extents: np.ndarray | list[float] | tuple[float, ...],
    ):
        c = _as_array(center)
        h = _as_array(half_extents)
        if c.shape != h.shape:
            raise ValueError("Box center and half_extents must have the same shape.")
        if np.any(h <= 0.0):
            raise ValueError("Box half_extents must be positive.")
        object.__setattr__(self, "center", c)
        object.__setattr__(self, "half_extents", h)

    def phi(self, x: np.ndarray | list[float] | tuple[float, ...]) -> float | np.ndarray:
        points = _as_array(x)
        q = np.abs(points - self.center) - self.half_extents
        outside = np.linalg.norm(np.maximum(q, 0.0), axis=-1)
        inside = np.minimum(np.max(q, axis=-1), 0.0)
        return _maybe_scalar(outside + inside)

    def grad(self, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
        points = _as_array(x)
        if not _is_single_point(points):
            return np.stack([self.grad(p) for p in points])

        h = 1.0e-6
        g = np.zeros_like(points, dtype=float)
        for i in range(points.size):
            step = np.zeros_like(points, dtype=float)
            step[i] = h
            g[i] = (self.phi(points + step) - self.phi(points - step)) / (2.0 * h)
        return _safe_normalize(g)

    signed_distance = phi
    gradient = grad


def phi(sdf: object, x: np.ndarray | list[float] | tuple[float, ...]) -> float | np.ndarray:
    """Evaluate a signed distance value using a small compatibility adapter."""

    for name in ("phi", "signed_distance", "distance", "evaluate"):
        method = getattr(sdf, name, None)
        if callable(method):
            return method(x)
    raise AttributeError("SDF object must provide phi(x) or a compatible distance method.")


def grad_phi(sdf: object, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    """Evaluate an SDF gradient using a small compatibility adapter."""

    for name in ("grad", "gradient", "normal"):
        method = getattr(sdf, name, None)
        if callable(method):
            return np.asarray(method(x), dtype=float)
    raise AttributeError("SDF object must provide grad(x) or a compatible gradient method.")

