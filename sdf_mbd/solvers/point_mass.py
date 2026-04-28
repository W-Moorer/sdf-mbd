"""Point-mass SDF contact solvers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from sdf_mbd.contact.ncp import complementarity_error, smooth_fischer_burmeister
from sdf_mbd.contact.penalty import penalty_normal_force_from_sdf
from sdf_mbd.contact.sdf_contact import point_contact_jacobian, sdf_gap, sdf_normal

try:  # SciPy is optional for this small prototype.
    from scipy.optimize import root as scipy_root
except Exception:  # pragma: no cover - exercised only when SciPy is unavailable.
    scipy_root = None


@dataclass(frozen=True)
class PointMassState:
    """2D point-mass state."""

    q: np.ndarray
    v: np.ndarray


@dataclass(frozen=True)
class PointMassDiagnostics:
    """Per-step contact and nonlinear-solve diagnostics."""

    solver_success: bool
    solver_message: str
    solver_iterations: int
    residual_norm: float
    gap: float
    penetration: float
    lambda_n: float
    ncp_residual: float
    complementarity_error: float


@dataclass(frozen=True)
class PointMassStepResult:
    """Point-mass step result."""

    state: PointMassState
    diagnostics: PointMassDiagnostics


def _as_state(q: np.ndarray | list[float], v: np.ndarray | list[float]) -> PointMassState:
    q_arr = np.asarray(q, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    if q_arr.shape != (2,) or v_arr.shape != (2,):
        raise ValueError("This prototype expects 2D q and v arrays with shape (2,).")
    return PointMassState(q=q_arr, v=v_arr)


def _newton_root(
    residual: Callable[[np.ndarray], np.ndarray],
    z0: np.ndarray,
    tol: float,
    max_iter: int,
) -> tuple[np.ndarray, bool, str, int]:
    z = np.asarray(z0, dtype=float).copy()
    h = 1.0e-7
    for iteration in range(1, max_iter + 1):
        r = residual(z)
        if np.linalg.norm(r) <= tol:
            return z, True, "Newton fallback converged.", iteration
        jac = np.zeros((z.size, z.size), dtype=float)
        for j in range(z.size):
            step = np.zeros_like(z)
            step[j] = h
            jac[:, j] = (residual(z + step) - residual(z - step)) / (2.0 * h)
        try:
            dz = np.linalg.solve(jac, -r)
        except np.linalg.LinAlgError:
            return z, False, "Newton fallback Jacobian was singular.", iteration
        z = z + dz
        if np.linalg.norm(dz) <= tol:
            return z, True, "Newton fallback step tolerance reached.", iteration
    return z, False, "Newton fallback reached max iterations.", max_iter


def ncp_point_mass_plane_step(
    sdf: object,
    state: PointMassState,
    mass: float,
    gravity: float,
    dt: float,
    eps: float,
    solver_tol: float = 1.0e-10,
    max_iter: int = 50,
) -> PointMassStepResult:
    """Advance one implicit Euler step for a 2D point mass with SDF-NCP contact.

    Unknowns are `z = [vx_next, vy_next, lambda_n]` with

        q_next = q + dt * v_next
        R_v = M (v_next - v) - dt (Q + J^T lambda_n)
        R_lambda = Phi_eps(phi(q_next), lambda_n).
    """

    if mass <= 0.0:
        raise ValueError("mass must be positive.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if eps < 0.0:
        raise ValueError("eps must be nonnegative.")

    q = np.asarray(state.q, dtype=float)
    v = np.asarray(state.v, dtype=float)
    force_ext = np.array([0.0, -float(mass) * float(gravity)], dtype=float)

    def residual(z: np.ndarray) -> np.ndarray:
        v_next = np.asarray(z[:2], dtype=float)
        lam = float(z[2])
        q_next = q + dt * v_next
        jac = point_contact_jacobian(sdf, q_next)
        gap = sdf_gap(sdf, q_next)
        r_v = mass * (v_next - v) - dt * (force_ext + jac * lam)
        r_lam = smooth_fischer_burmeister(gap, lam, eps)
        return np.array([r_v[0], r_v[1], r_lam], dtype=float)

    v_guess = v + dt * force_ext / mass
    q_guess = q + dt * v_guess
    lam_guess = max(0.0, -sdf_gap(sdf, q_guess)) * mass / max(dt * dt, 1.0e-12)
    z0 = np.array([v_guess[0], v_guess[1], lam_guess], dtype=float)

    if scipy_root is not None:
        solve = scipy_root(residual, z0, method="hybr", options={"xtol": solver_tol, "maxfev": max_iter})
        z = np.asarray(solve.x, dtype=float)
        success = bool(solve.success)
        message = str(solve.message)
        iterations = int(getattr(solve, "nfev", -1))
    else:
        z, success, message, iterations = _newton_root(residual, z0, solver_tol, max_iter)

    lam = float(z[2])
    if lam < -1.0e-8:
        success = False
        message = f"{message} Negative contact multiplier {lam:.6e}."
    elif lam < 0.0:
        z[2] = 0.0
        lam = 0.0

    res = residual(z)
    v_next = z[:2].copy()
    q_next = q + dt * v_next
    gap = sdf_gap(sdf, q_next)
    ncp_residual = float(smooth_fischer_burmeister(gap, lam, eps))
    residual_norm = float(np.linalg.norm(res))
    penetration = max(0.0, -gap)
    diagnostics = PointMassDiagnostics(
        solver_success=success and residual_norm <= max(1.0e-7, 100.0 * solver_tol),
        solver_message=message,
        solver_iterations=iterations,
        residual_norm=residual_norm,
        gap=gap,
        penetration=penetration,
        lambda_n=lam,
        ncp_residual=abs(ncp_residual),
        complementarity_error=complementarity_error(gap, lam),
    )
    return PointMassStepResult(state=PointMassState(q=q_next, v=v_next), diagnostics=diagnostics)


def penalty_point_mass_step(
    sdf: object,
    state: PointMassState,
    mass: float,
    gravity: float,
    dt: float,
    k_n: float,
    c_n: float = 0.0,
) -> PointMassStepResult:
    """Advance one semi-implicit Euler step using an SDF penalty force."""

    if mass <= 0.0:
        raise ValueError("mass must be positive.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")

    q = np.asarray(state.q, dtype=float)
    v = np.asarray(state.v, dtype=float)
    normal = sdf_normal(sdf, q)
    g_dot = float(np.dot(normal, v))
    contact = penalty_normal_force_from_sdf(sdf, q, k_n=k_n, g_dot=g_dot, c_n=c_n)
    force_ext = np.array([0.0, -float(mass) * float(gravity)], dtype=float)
    v_next = v + dt * (force_ext + contact.force) / mass
    q_next = q + dt * v_next
    gap_next = sdf_gap(sdf, q_next)
    ncp_residual = float(smooth_fischer_burmeister(gap_next, contact.lambda_n, 0.0))
    diagnostics = PointMassDiagnostics(
        solver_success=True,
        solver_message="Penalty step.",
        solver_iterations=0,
        residual_norm=0.0,
        gap=gap_next,
        penetration=max(0.0, -gap_next),
        lambda_n=contact.lambda_n,
        ncp_residual=abs(ncp_residual),
        complementarity_error=complementarity_error(gap_next, contact.lambda_n),
    )
    return PointMassStepResult(state=PointMassState(q=q_next, v=v_next), diagnostics=diagnostics)


def simulate_ncp_point_mass(
    sdf: object,
    q0: np.ndarray | list[float],
    v0: np.ndarray | list[float],
    mass: float,
    gravity: float,
    dt: float,
    t_end: float,
    eps: float,
    solver_tol: float = 1.0e-10,
) -> list[tuple[float, PointMassState, PointMassDiagnostics]]:
    """Run a point-mass SDF-NCP trajectory and return time-stamped rows."""

    state = _as_state(q0, v0)
    rows: list[tuple[float, PointMassState, PointMassDiagnostics]] = []
    steps = int(round(t_end / dt))
    for i in range(steps + 1):
        if i == 0:
            gap = sdf_gap(sdf, state.q)
            diag = PointMassDiagnostics(
                solver_success=True,
                solver_message="Initial state.",
                solver_iterations=0,
                residual_norm=0.0,
                gap=gap,
                penetration=max(0.0, -gap),
                lambda_n=0.0,
                ncp_residual=abs(float(smooth_fischer_burmeister(gap, 0.0, eps))),
                complementarity_error=complementarity_error(gap, 0.0),
            )
        else:
            result = ncp_point_mass_plane_step(sdf, state, mass, gravity, dt, eps, solver_tol=solver_tol)
            state = result.state
            diag = result.diagnostics
        rows.append((i * dt, state, diag))
    return rows


def simulate_penalty_point_mass(
    sdf: object,
    q0: np.ndarray | list[float],
    v0: np.ndarray | list[float],
    mass: float,
    gravity: float,
    dt: float,
    t_end: float,
    k_n: float,
    c_n: float = 0.0,
) -> list[tuple[float, PointMassState, PointMassDiagnostics]]:
    """Run a point-mass penalty-contact trajectory and return time-stamped rows."""

    state = _as_state(q0, v0)
    rows: list[tuple[float, PointMassState, PointMassDiagnostics]] = []
    steps = int(round(t_end / dt))
    for i in range(steps + 1):
        if i == 0:
            gap = sdf_gap(sdf, state.q)
            diag = PointMassDiagnostics(
                solver_success=True,
                solver_message="Initial state.",
                solver_iterations=0,
                residual_norm=0.0,
                gap=gap,
                penetration=max(0.0, -gap),
                lambda_n=0.0,
                ncp_residual=abs(float(smooth_fischer_burmeister(gap, 0.0, 0.0))),
                complementarity_error=complementarity_error(gap, 0.0),
            )
        else:
            result = penalty_point_mass_step(sdf, state, mass, gravity, dt, k_n, c_n=c_n)
            state = result.state
            diag = result.diagnostics
        rows.append((i * dt, state, diag))
    return rows

