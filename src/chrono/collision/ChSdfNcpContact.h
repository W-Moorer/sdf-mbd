// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2024 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
//
// Minimal SDF-NCP contact utilities.
//
// This header is intentionally independent from Chrono's contact containers and
// global solver descriptors. It provides the first C++ regression kernel for the
// classical SDF + nonlinear complementarity contact prototype.
//
// =============================================================================

#ifndef CH_SDF_NCP_CONTACT_H
#define CH_SDF_NCP_CONTACT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace chrono {
namespace sdfncp {

inline double RegularizedSmoothingEps(double eps) {
    return std::max(std::abs(eps), 1.0e-15);
}

inline double SmoothFischerBurmeister(double g, double lambda, double eps) {
    const double e = RegularizedSmoothingEps(eps);
    return std::sqrt(g * g + lambda * lambda + e * e) - g - lambda;
}

struct SmoothFischerBurmeisterGradient {
    double dPhi_dg = 0.0;
    double dPhi_dlambda = 0.0;
};

inline SmoothFischerBurmeisterGradient SmoothFischerBurmeisterGrad(double g, double lambda, double eps) {
    const double e = RegularizedSmoothingEps(eps);
    const double denom = std::max(std::sqrt(g * g + lambda * lambda + e * e),
                                  std::numeric_limits<double>::min());

    SmoothFischerBurmeisterGradient grad;
    grad.dPhi_dg = g / denom - 1.0;
    grad.dPhi_dlambda = lambda / denom - 1.0;
    return grad;
}

inline double ComplementarityError(double g, double lambda) {
    return std::max(0.0, -g) + std::max(0.0, -lambda) + std::abs(g * lambda);
}

struct PointMassPlaneState {
    double x = 0.0;
    double y = 0.0;
    double vx = 0.0;
    double vy = 0.0;
};

struct PointMassPlaneSettings {
    double mass = 1.0;
    double gravity = 9.81;
    double dt = 1.0e-3;
    double eps = 1.0e-6;
    double tolerance = 1.0e-10;
    int max_iterations = 40;
};

struct PointMassPlaneResidual {
    std::array<double, 3> value = {{0.0, 0.0, 0.0}};
    double gap = 0.0;
    double lambda_n = 0.0;
    double ncp_residual = 0.0;
};

struct PointMassPlaneDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double gap = 0.0;
    double lambda_n = 0.0;
    double penetration = 0.0;
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
    std::string message;
};

struct PointMassPlaneStepResult {
    PointMassPlaneState state;
    PointMassPlaneDiagnostics diagnostics;
};

inline double Norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline PointMassPlaneResidual ComputePointMassPlaneResidual(const PointMassPlaneState& state,
                                                            const PointMassPlaneSettings& settings,
                                                            const std::array<double, 3>& z) {
    const double vx_next = z[0];
    const double vy_next = z[1];
    const double lambda_n = z[2];
    const double gap_next = state.y + settings.dt * vy_next;

    PointMassPlaneResidual residual;
    residual.value[0] = settings.mass * (vx_next - state.vx);
    residual.value[1] = settings.mass * (vy_next - state.vy) -
                        settings.dt * (-settings.mass * settings.gravity + lambda_n);
    residual.value[2] = SmoothFischerBurmeister(gap_next, lambda_n, settings.eps);
    residual.gap = gap_next;
    residual.lambda_n = lambda_n;
    residual.ncp_residual = residual.value[2];
    return residual;
}

inline std::array<std::array<double, 3>, 3> FiniteDifferencePointMassPlaneJacobian(
    const PointMassPlaneState& state,
    const PointMassPlaneSettings& settings,
    const std::array<double, 3>& z,
    double h = 1.0e-7) {
    std::array<std::array<double, 3>, 3> jac = {{{{0.0, 0.0, 0.0}},
                                                 {{0.0, 0.0, 0.0}},
                                                 {{0.0, 0.0, 0.0}}}};

    for (int col = 0; col < 3; col++) {
        std::array<double, 3> zp = z;
        std::array<double, 3> zm = z;
        zp[col] += h;
        zm[col] -= h;

        const auto rp = ComputePointMassPlaneResidual(state, settings, zp).value;
        const auto rm = ComputePointMassPlaneResidual(state, settings, zm).value;

        for (int row = 0; row < 3; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }

    return jac;
}

inline bool SolveLinear3x3(std::array<std::array<double, 3>, 3> A,
                           std::array<double, 3> b,
                           std::array<double, 3>& x) {
    for (int col = 0; col < 3; col++) {
        int pivot = col;
        double pivot_abs = std::abs(A[col][col]);
        for (int row = col + 1; row < 3; row++) {
            const double candidate = std::abs(A[row][col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }

        if (pivot_abs <= 1.0e-14) {
            return false;
        }

        if (pivot != col) {
            std::swap(A[pivot], A[col]);
            std::swap(b[pivot], b[col]);
        }

        for (int row = col + 1; row < 3; row++) {
            const double factor = A[row][col] / A[col][col];
            A[row][col] = 0.0;
            for (int j = col + 1; j < 3; j++) {
                A[row][j] -= factor * A[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    for (int row = 2; row >= 0; row--) {
        double rhs = b[row];
        for (int col = row + 1; col < 3; col++) {
            rhs -= A[row][col] * x[col];
        }
        if (std::abs(A[row][row]) <= 1.0e-14) {
            return false;
        }
        x[row] = rhs / A[row][row];
    }

    return true;
}

inline PointMassPlaneDiagnostics MakePointMassPlaneDiagnostics(const PointMassPlaneState& state,
                                                              const PointMassPlaneSettings& settings,
                                                              const std::array<double, 3>& z,
                                                              bool success,
                                                              int iterations,
                                                              const std::string& message) {
    const PointMassPlaneResidual residual = ComputePointMassPlaneResidual(state, settings, z);
    PointMassPlaneDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.residual_norm = Norm3(residual.value);
    diagnostics.gap = residual.gap;
    diagnostics.lambda_n = residual.lambda_n;
    diagnostics.penetration = std::max(0.0, -residual.gap);
    diagnostics.ncp_residual = std::abs(residual.ncp_residual);
    diagnostics.complementarity_error = ComplementarityError(residual.gap, residual.lambda_n);
    diagnostics.message = message;
    return diagnostics;
}

inline PointMassPlaneStepResult SolvePointMassPlaneStep(const PointMassPlaneState& state,
                                                        const PointMassPlaneSettings& settings) {
    PointMassPlaneStepResult result;
    result.state = state;

    if (settings.mass <= 0.0 || settings.dt <= 0.0 || settings.max_iterations <= 0) {
        result.diagnostics =
            MakePointMassPlaneDiagnostics(state, settings, {{state.vx, state.vy, 0.0}}, false, 0,
                                          "Invalid point-mass SDF-NCP settings.");
        return result;
    }

    const double vx_guess = state.vx;
    const double vy_guess = state.vy - settings.dt * settings.gravity;
    const double gap_guess = state.y + settings.dt * vy_guess;
    const double lambda_guess =
        std::max(0.0, -gap_guess) * settings.mass / std::max(settings.dt * settings.dt, 1.0e-12);

    std::array<double, 3> z = {{vx_guess, vy_guess, lambda_guess}};
    double residual_norm = Norm3(ComputePointMassPlaneResidual(state, settings, z).value);

    int iterations = 0;
    bool success = residual_norm <= settings.tolerance;
    std::string message = success ? "Initial guess satisfies residual tolerance." : "Newton did not run.";

    for (iterations = 0; !success && iterations < settings.max_iterations; iterations++) {
        const auto residual = ComputePointMassPlaneResidual(state, settings, z).value;
        const auto jac = FiniteDifferencePointMassPlaneJacobian(state, settings, z);
        std::array<double, 3> rhs = {{-residual[0], -residual[1], -residual[2]}};
        std::array<double, 3> dz = {{0.0, 0.0, 0.0}};

        if (!SolveLinear3x3(jac, rhs, dz)) {
            message = "Newton finite-difference Jacobian was singular.";
            break;
        }

        double alpha = 1.0;
        std::array<double, 3> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;

        for (int ls = 0; ls < 12; ls++) {
            candidate = {{z[0] + alpha * dz[0], z[1] + alpha * dz[1], z[2] + alpha * dz[2]}};
            candidate_norm = Norm3(ComputePointMassPlaneResidual(state, settings, candidate).value);
            if (candidate_norm <= residual_norm || alpha <= 1.0e-4) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if (!accepted) {
            message = "Newton line search failed.";
            break;
        }

        z = candidate;
        residual_norm = candidate_norm;

        if (residual_norm <= settings.tolerance) {
            success = true;
            message = "Newton converged.";
            iterations++;
            break;
        }
    }

    if (!success && residual_norm <= std::max(settings.tolerance, 1.0e-10)) {
        success = true;
        message = "Newton reached relaxed residual tolerance.";
    } else if (!success && message == "Newton did not run.") {
        message = "Newton reached max iterations.";
    }

    if (z[2] < -1.0e-8) {
        success = false;
        message = "Newton returned a significantly negative normal multiplier.";
    } else if (z[2] < 0.0) {
        z[2] = 0.0;
    }

    result.state.x = state.x + settings.dt * z[0];
    result.state.y = state.y + settings.dt * z[1];
    result.state.vx = z[0];
    result.state.vy = z[1];
    result.diagnostics = MakePointMassPlaneDiagnostics(state, settings, z, success, iterations, message);
    return result;
}

}  // namespace sdfncp
}  // namespace chrono

#endif

