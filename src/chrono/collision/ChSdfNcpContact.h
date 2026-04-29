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
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "chrono/collision/ChSdfContactGeometry.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChVector3.h"

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

inline ChVector3d NormalizeSdfVector(const ChVector3d& v, const char* message) {
    const double length = v.Length();
    if (length <= 1.0e-14) {
        throw std::invalid_argument(message);
    }
    return v / length;
}

class ChSdfSurface {
  public:
    virtual ~ChSdfSurface() = default;

    virtual double Phi(const ChVector3d& x) const = 0;
    virtual ChVector3d Grad(const ChVector3d& x) const = 0;
};

class ChPlaneSdfSurface : public ChSdfSurface {
  public:
    ChPlaneSdfSurface(const ChVector3d& normal, double offset) : m_offset(offset) {
        m_normal = NormalizeSdfVector(normal, "ChPlaneSdfSurface normal must be nonzero.");
    }

    double Phi(const ChVector3d& x) const override {
        return m_normal.Dot(x) - m_offset;
    }

    ChVector3d Grad(const ChVector3d& x) const override {
        return m_normal;
    }

    const ChVector3d& GetNormal() const {
        return m_normal;
    }

    double GetOffset() const {
        return m_offset;
    }

  private:
    ChVector3d m_normal = ChVector3d(0, 1, 0);
    double m_offset = 0.0;
};

class ChSphereSdfSurface : public ChSdfSurface {
  public:
    ChSphereSdfSurface(const ChVector3d& center, double radius) : m_center(center), m_radius(radius) {
        if (radius <= 0.0) {
            throw std::invalid_argument("ChSphereSdfSurface radius must be positive.");
        }
    }

    double Phi(const ChVector3d& x) const override {
        return (x - m_center).Length() - m_radius;
    }

    ChVector3d Grad(const ChVector3d& x) const override {
        const ChVector3d delta = x - m_center;
        if (delta.Length() <= 1.0e-14) {
            return ChVector3d(1, 0, 0);
        }
        return delta.GetNormalized();
    }

  private:
    ChVector3d m_center = ChVector3d(0, 0, 0);
    double m_radius = 1.0;
};

struct ChSdfPointKinematicsState {
    ChVector3d position = ChVector3d(0, 0, 0);
};

class ChSdfPointKinematics {
  public:
    ChVector3d GetPoint(const ChVector3d& q) const {
        return q;
    }

    ChMatrix33d GetPointJacobianWrtPosition() const {
        return ChMatrix33d(1.0);
    }
};

struct ChSdfNcpContactConstraint {
    double gap = 0.0;
    ChVector3d normal = ChVector3d(0, 1, 0);
    ChVector3d jacobian_position = ChVector3d(0, 1, 0);
    double lambda_n = 0.0;
    double penetration = 0.0;
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
};

inline ChSdfNcpContactConstraint EvaluateSdfNcpContact(const ChSdfSurface& surface,
                                                       const ChVector3d& q,
                                                       double lambda,
                                                       double eps) {
    ChSdfPointKinematics kinematics;
    const ChVector3d x = kinematics.GetPoint(q);
    const ChVector3d normal = NormalizeSdfVector(surface.Grad(x), "SDF gradient must be nonzero.");
    const double gap = surface.Phi(x);

    ChSdfNcpContactConstraint constraint;
    constraint.gap = gap;
    constraint.normal = normal;
    constraint.jacobian_position = kinematics.GetPointJacobianWrtPosition() * normal;
    constraint.lambda_n = lambda;
    constraint.penetration = std::max(0.0, -gap);
    constraint.ncp_residual = SmoothFischerBurmeister(gap, lambda, eps);
    constraint.complementarity_error = ComplementarityError(gap, lambda);
    return constraint;
}

inline ChSdfNcpContactConstraint EvaluateSdfNcpContactQuery(const chrono::sdfcontact::ChSdfContactSampleQuery& query,
                                                            double lambda,
                                                            double eps) {
    const ChVector3d normal = NormalizeSdfVector(query.grad, "SDF query gradient must be nonzero.");

    ChSdfNcpContactConstraint constraint;
    constraint.gap = query.phi;
    constraint.normal = normal;
    constraint.jacobian_position = normal;
    constraint.lambda_n = lambda;
    constraint.penetration = std::max(0.0, -query.phi);
    constraint.ncp_residual = SmoothFischerBurmeister(query.phi, lambda, eps);
    constraint.complementarity_error = ComplementarityError(query.phi, lambda);
    return constraint;
}

struct ChSdfNcpContactPoint {
    std::shared_ptr<const ChSdfSurface> surface;
    double lambda_n = 0.0;
    double gap = 0.0;
    double penetration = 0.0;
    ChVector3d normal = ChVector3d(0, 1, 0);
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
};

class ChSdfNcpContactSet {
  public:
    void AddContact(std::shared_ptr<const ChSdfSurface> surface, double lambda = 0.0) {
        if (!surface) {
            throw std::invalid_argument("ChSdfNcpContactSet contact surface must not be null.");
        }
        ChSdfNcpContactPoint point;
        point.surface = std::move(surface);
        point.lambda_n = lambda;
        m_points.push_back(point);
    }

    size_t GetNumContacts() const {
        return m_points.size();
    }

    bool Empty() const {
        return m_points.empty();
    }

    const ChSdfNcpContactPoint& GetContact(size_t index) const {
        return m_points.at(index);
    }

    ChSdfNcpContactPoint& GetContact(size_t index) {
        return m_points.at(index);
    }

    void SetLambda(size_t index, double lambda) {
        m_points.at(index).lambda_n = lambda;
    }

    double GetLambda(size_t index) const {
        return m_points.at(index).lambda_n;
    }

    void Update(const ChVector3d& q_next, double eps) {
        for (auto& point : m_points) {
            if (!point.surface) {
                throw std::invalid_argument("ChSdfNcpContactSet contains a null surface.");
            }
            const auto constraint = EvaluateSdfNcpContact(*point.surface, q_next, point.lambda_n, eps);
            point.gap = constraint.gap;
            point.penetration = constraint.penetration;
            point.normal = constraint.normal;
            point.ncp_residual = constraint.ncp_residual;
            point.complementarity_error = constraint.complementarity_error;
        }
    }

    ChVector3d ComputeTotalContactForce() const {
        ChVector3d force(0, 0, 0);
        for (const auto& point : m_points) {
            force += point.normal * point.lambda_n;
        }
        return force;
    }

    double MaxPenetration() const {
        double value = 0.0;
        for (const auto& point : m_points) {
            value = std::max(value, point.penetration);
        }
        return value;
    }

    double MaxComplementarityError() const {
        double value = 0.0;
        for (const auto& point : m_points) {
            value = std::max(value, point.complementarity_error);
        }
        return value;
    }

    double MeanComplementarityError() const {
        if (m_points.empty()) {
            return 0.0;
        }
        double sum = 0.0;
        for (const auto& point : m_points) {
            sum += point.complementarity_error;
        }
        return sum / static_cast<double>(m_points.size());
    }

    double NcpResidualNorm() const {
        double sum_sq = 0.0;
        for (const auto& point : m_points) {
            sum_sq += point.ncp_residual * point.ncp_residual;
        }
        return std::sqrt(sum_sq);
    }

  private:
    std::vector<ChSdfNcpContactPoint> m_points;
};

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

struct ChSdfNcpPointMassState {
    ChVector3d position = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
};

struct ChSdfNcpPointMassSettings {
    double mass = 1.0;
    double gravity = 9.81;
    double dt = 1.0e-3;
    double eps = 1.0e-6;
    double tolerance = 1.0e-10;
    int max_iterations = 40;
};

struct ChSdfNcpPointMassResidual {
    std::array<double, 4> value = {{0.0, 0.0, 0.0, 0.0}};
    ChSdfNcpContactConstraint contact;
};

struct ChSdfNcpPointMassDiagnostics {
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

struct ChSdfNcpPointMassStepResult {
    ChSdfNcpPointMassState state;
    ChSdfNcpPointMassDiagnostics diagnostics;
};

inline double Norm4(const std::array<double, 4>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
}

inline ChSdfNcpPointMassResidual ComputeSdfNcpPointMassResidual(const ChSdfSurface& surface,
                                                                const ChSdfNcpPointMassState& state,
                                                                const ChSdfNcpPointMassSettings& settings,
                                                                const std::array<double, 4>& z) {
    const ChVector3d v_next(z[0], z[1], z[2]);
    const double lambda_n = z[3];
    const ChVector3d q_next = state.position + v_next * settings.dt;
    const ChVector3d external_force(0.0, -settings.mass * settings.gravity, 0.0);

    ChSdfNcpPointMassResidual residual;
    residual.contact = EvaluateSdfNcpContact(surface, q_next, lambda_n, settings.eps);

    const ChVector3d r_v =
        settings.mass * (v_next - state.velocity) -
        settings.dt * (external_force + residual.contact.jacobian_position * lambda_n);

    residual.value[0] = r_v.x();
    residual.value[1] = r_v.y();
    residual.value[2] = r_v.z();
    residual.value[3] = residual.contact.ncp_residual;
    return residual;
}

inline std::array<std::array<double, 4>, 4> FiniteDifferenceSdfNcpPointMassJacobian(
    const ChSdfSurface& surface,
    const ChSdfNcpPointMassState& state,
    const ChSdfNcpPointMassSettings& settings,
    const std::array<double, 4>& z,
    double h = 1.0e-7) {
    std::array<std::array<double, 4>, 4> jac = {{{{0.0, 0.0, 0.0, 0.0}},
                                                 {{0.0, 0.0, 0.0, 0.0}},
                                                 {{0.0, 0.0, 0.0, 0.0}},
                                                 {{0.0, 0.0, 0.0, 0.0}}}};

    for (int col = 0; col < 4; col++) {
        std::array<double, 4> zp = z;
        std::array<double, 4> zm = z;
        zp[col] += h;
        zm[col] -= h;

        const auto rp = ComputeSdfNcpPointMassResidual(surface, state, settings, zp).value;
        const auto rm = ComputeSdfNcpPointMassResidual(surface, state, settings, zm).value;

        for (int row = 0; row < 4; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }

    return jac;
}

inline bool SolveLinear4x4(std::array<std::array<double, 4>, 4> A,
                           std::array<double, 4> b,
                           std::array<double, 4>& x) {
    for (int col = 0; col < 4; col++) {
        int pivot = col;
        double pivot_abs = std::abs(A[col][col]);
        for (int row = col + 1; row < 4; row++) {
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

        for (int row = col + 1; row < 4; row++) {
            const double factor = A[row][col] / A[col][col];
            A[row][col] = 0.0;
            for (int j = col + 1; j < 4; j++) {
                A[row][j] -= factor * A[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    for (int row = 3; row >= 0; row--) {
        double rhs = b[row];
        for (int col = row + 1; col < 4; col++) {
            rhs -= A[row][col] * x[col];
        }
        if (std::abs(A[row][row]) <= 1.0e-14) {
            return false;
        }
        x[row] = rhs / A[row][row];
    }

    return true;
}

inline ChSdfNcpPointMassDiagnostics MakeSdfNcpPointMassDiagnostics(const ChSdfSurface& surface,
                                                                   const ChSdfNcpPointMassState& state,
                                                                   const ChSdfNcpPointMassSettings& settings,
                                                                   const std::array<double, 4>& z,
                                                                   bool success,
                                                                   int iterations,
                                                                   const std::string& message) {
    const ChSdfNcpPointMassResidual residual = ComputeSdfNcpPointMassResidual(surface, state, settings, z);

    ChSdfNcpPointMassDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.residual_norm = Norm4(residual.value);
    diagnostics.gap = residual.contact.gap;
    diagnostics.lambda_n = residual.contact.lambda_n;
    diagnostics.penetration = residual.contact.penetration;
    diagnostics.ncp_residual = std::abs(residual.contact.ncp_residual);
    diagnostics.complementarity_error = residual.contact.complementarity_error;
    diagnostics.message = message;
    return diagnostics;
}

inline ChSdfNcpPointMassStepResult SolveSdfNcpPointMassStep(const ChSdfSurface& surface,
                                                            const ChSdfNcpPointMassState& state,
                                                            const ChSdfNcpPointMassSettings& settings) {
    ChSdfNcpPointMassStepResult result;
    result.state = state;

    const std::array<double, 4> invalid_z = {
        {state.velocity.x(), state.velocity.y(), state.velocity.z(), 0.0}};
    if (settings.mass <= 0.0 || settings.dt <= 0.0 || settings.max_iterations <= 0) {
        result.diagnostics = MakeSdfNcpPointMassDiagnostics(surface, state, settings, invalid_z, false, 0,
                                                            "Invalid SDF-NCP point-mass settings.");
        return result;
    }

    const ChVector3d external_accel(0.0, -settings.gravity, 0.0);
    const ChVector3d v_guess = state.velocity + external_accel * settings.dt;
    const ChVector3d q_guess = state.position + v_guess * settings.dt;
    const double gap_guess = surface.Phi(q_guess);
    const double lambda_guess =
        std::max(0.0, -gap_guess) * settings.mass / std::max(settings.dt * settings.dt, 1.0e-12);

    std::array<double, 4> z = {{v_guess.x(), v_guess.y(), v_guess.z(), lambda_guess}};
    double residual_norm = Norm4(ComputeSdfNcpPointMassResidual(surface, state, settings, z).value);

    int iterations = 0;
    bool success = residual_norm <= settings.tolerance;
    std::string message = success ? "Initial guess satisfies residual tolerance." : "Newton did not run.";

    for (iterations = 0; !success && iterations < settings.max_iterations; iterations++) {
        const auto residual = ComputeSdfNcpPointMassResidual(surface, state, settings, z).value;
        const auto jac = FiniteDifferenceSdfNcpPointMassJacobian(surface, state, settings, z);
        std::array<double, 4> rhs = {{-residual[0], -residual[1], -residual[2], -residual[3]}};
        std::array<double, 4> dz = {{0.0, 0.0, 0.0, 0.0}};

        if (!SolveLinear4x4(jac, rhs, dz)) {
            message = "Newton finite-difference Jacobian was singular.";
            break;
        }

        double alpha = 1.0;
        std::array<double, 4> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;

        for (int ls = 0; ls < 12; ls++) {
            candidate = {{z[0] + alpha * dz[0],
                          z[1] + alpha * dz[1],
                          z[2] + alpha * dz[2],
                          z[3] + alpha * dz[3]}};
            candidate_norm = Norm4(ComputeSdfNcpPointMassResidual(surface, state, settings, candidate).value);
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

    if (z[3] < -1.0e-8) {
        success = false;
        message = "Newton returned a significantly negative normal multiplier.";
    } else if (z[3] < 0.0) {
        z[3] = 0.0;
    }

    const ChVector3d v_next(z[0], z[1], z[2]);
    result.state.velocity = v_next;
    result.state.position = state.position + v_next * settings.dt;
    result.diagnostics = MakeSdfNcpPointMassDiagnostics(surface, state, settings, z, success, iterations, message);
    return result;
}

struct ChSdfNcpPointMassContactSetResidual {
    std::vector<double> value;
    ChSdfNcpContactSet contact_set;
};

struct ChSdfNcpPointMassContactSetDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double max_penetration = 0.0;
    double min_lambda_n = 0.0;
    double max_complementarity_error = 0.0;
    double mean_complementarity_error = 0.0;
    double ncp_residual_norm = 0.0;
    std::string message;
};

struct ChSdfNcpPointMassContactSetStepResult {
    ChSdfNcpPointMassState state;
    ChSdfNcpContactSet contact_set;
    ChSdfNcpPointMassContactSetDiagnostics diagnostics;
};

inline double NormDynamic(const std::vector<double>& v) {
    double sum_sq = 0.0;
    for (double value : v) {
        sum_sq += value * value;
    }
    return std::sqrt(sum_sq);
}

inline ChSdfNcpPointMassContactSetResidual ComputeSdfNcpPointMassContactSetResidual(
    const ChSdfNcpContactSet& contact_set,
    const ChSdfNcpPointMassState& state,
    const ChSdfNcpPointMassSettings& settings,
    const std::vector<double>& z) {
    const size_t n_contacts = contact_set.GetNumContacts();
    if (z.size() != 3 + n_contacts) {
        throw std::invalid_argument("SDF-NCP contact-set unknown vector has inconsistent size.");
    }

    const ChVector3d v_next(z[0], z[1], z[2]);
    const ChVector3d q_next = state.position + v_next * settings.dt;
    const ChVector3d external_force(0.0, -settings.mass * settings.gravity, 0.0);

    ChSdfNcpPointMassContactSetResidual residual;
    residual.contact_set = contact_set;
    for (size_t i = 0; i < n_contacts; i++) {
        residual.contact_set.SetLambda(i, z[3 + i]);
    }
    residual.contact_set.Update(q_next, settings.eps);

    const ChVector3d contact_force = residual.contact_set.ComputeTotalContactForce();
    const ChVector3d r_v =
        settings.mass * (v_next - state.velocity) - settings.dt * (external_force + contact_force);

    residual.value.assign(3 + n_contacts, 0.0);
    residual.value[0] = r_v.x();
    residual.value[1] = r_v.y();
    residual.value[2] = r_v.z();
    for (size_t i = 0; i < n_contacts; i++) {
        residual.value[3 + i] = residual.contact_set.GetContact(i).ncp_residual;
    }

    return residual;
}

inline std::vector<std::vector<double>> FiniteDifferenceSdfNcpPointMassContactSetJacobian(
    const ChSdfNcpContactSet& contact_set,
    const ChSdfNcpPointMassState& state,
    const ChSdfNcpPointMassSettings& settings,
    const std::vector<double>& z,
    double h = 1.0e-7) {
    const size_t n = z.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));

    for (size_t col = 0; col < n; col++) {
        std::vector<double> zp = z;
        std::vector<double> zm = z;
        zp[col] += h;
        zm[col] -= h;

        const auto rp = ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, zp).value;
        const auto rm = ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, zm).value;

        for (size_t row = 0; row < n; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }

    return jac;
}

inline bool SolveLinearDynamic(std::vector<std::vector<double>> A,
                               std::vector<double> b,
                               std::vector<double>& x) {
    const size_t n = A.size();
    if (b.size() != n || x.size() != n) {
        return false;
    }

    for (size_t row = 0; row < n; row++) {
        if (A[row].size() != n) {
            return false;
        }
    }

    for (size_t col = 0; col < n; col++) {
        size_t pivot = col;
        double pivot_abs = std::abs(A[col][col]);
        for (size_t row = col + 1; row < n; row++) {
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

        for (size_t row = col + 1; row < n; row++) {
            const double factor = A[row][col] / A[col][col];
            A[row][col] = 0.0;
            for (size_t j = col + 1; j < n; j++) {
                A[row][j] -= factor * A[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    for (int row = static_cast<int>(n) - 1; row >= 0; row--) {
        double rhs = b[static_cast<size_t>(row)];
        for (size_t col = static_cast<size_t>(row) + 1; col < n; col++) {
            rhs -= A[static_cast<size_t>(row)][col] * x[col];
        }
        if (std::abs(A[static_cast<size_t>(row)][static_cast<size_t>(row)]) <= 1.0e-14) {
            return false;
        }
        x[static_cast<size_t>(row)] = rhs / A[static_cast<size_t>(row)][static_cast<size_t>(row)];
    }

    return true;
}

struct ChSdfNcpGeneralizedContact {
    double gap = 0.0;
    std::vector<double> jacobian;
    std::vector<double> jacobian_velocity_derivative;
    double normal_velocity_offset = 0.0;
    double weight = 1.0;
    double lambda_n = 0.0;
    double penetration = 0.0;
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
    int contact_id = -1;
    int patch_id = -1;
};

struct ChSdfNcpGeneralizedProblem {
    double dt = 1.0e-3;
    double eps = 1.0e-6;
    double gap_scale = 1.0;
    double lambda_scale = 1.0;
    double tolerance = 1.0e-10;
    double relaxed_tolerance = 1.0e-8;
    double negative_lambda_tolerance = 1.0e-8;
    int max_iterations = 40;
    bool use_analytic_jacobian = false;
    std::vector<double> current_velocity;
    std::vector<double> mass_diagonal;
    std::vector<double> external_force;
    size_t contact_count = 0;
    std::function<std::vector<ChSdfNcpGeneralizedContact>(const std::vector<double>&)> evaluate_contacts;
};

struct ChSdfNcpGeneralizedResidual {
    std::vector<double> value;
    std::vector<ChSdfNcpGeneralizedContact> contacts;
};

struct ChSdfNcpGeneralizedDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double max_penetration = 0.0;
    double min_lambda_n = 0.0;
    double max_lambda_n = 0.0;
    double max_complementarity_error = 0.0;
    double mean_complementarity_error = 0.0;
    double ncp_residual_norm = 0.0;
    std::string message;
};

struct ChSdfNcpGeneralizedStepResult {
    std::vector<double> next_velocity;
    std::vector<double> lambdas;
    std::vector<ChSdfNcpGeneralizedContact> contacts;
    ChSdfNcpGeneralizedDiagnostics diagnostics;
};

inline void ValidateSdfNcpGeneralizedProblem(const ChSdfNcpGeneralizedProblem& problem) {
    const size_t n_v = problem.current_velocity.size();
    if (n_v == 0) {
        throw std::invalid_argument("Generalized SDF-NCP problem must have at least one velocity DOF.");
    }
    if (problem.mass_diagonal.size() != n_v || problem.external_force.size() != n_v) {
        throw std::invalid_argument("Generalized SDF-NCP mass, velocity, and external-force dimensions differ.");
    }
    for (double mass : problem.mass_diagonal) {
        if (!(mass > 0.0) || !std::isfinite(mass)) {
            throw std::invalid_argument("Generalized SDF-NCP mass diagonal entries must be finite and positive.");
        }
    }
    if (!(problem.gap_scale > 0.0) || !std::isfinite(problem.gap_scale) ||
        !(problem.lambda_scale > 0.0) || !std::isfinite(problem.lambda_scale)) {
        throw std::invalid_argument("Generalized SDF-NCP gap and lambda scales must be finite and positive.");
    }
    if (!problem.evaluate_contacts) {
        throw std::invalid_argument("Generalized SDF-NCP problem needs an evaluate_contacts callback.");
    }
}

inline ChSdfNcpGeneralizedResidual ComputeSdfNcpGeneralizedResidual(const ChSdfNcpGeneralizedProblem& problem,
                                                                    const std::vector<double>& z) {
    ValidateSdfNcpGeneralizedProblem(problem);
    const size_t n_v = problem.current_velocity.size();
    const size_t n_c = problem.contact_count;
    if (z.size() != n_v + n_c) {
        throw std::invalid_argument("Generalized SDF-NCP unknown vector has inconsistent size.");
    }

    std::vector<double> v_next(z.begin(), z.begin() + static_cast<std::ptrdiff_t>(n_v));
    std::vector<ChSdfNcpGeneralizedContact> contacts = problem.evaluate_contacts(v_next);
    if (contacts.size() != n_c) {
        throw std::invalid_argument("Generalized SDF-NCP contact callback returned an inconsistent contact count.");
    }

    std::vector<double> contact_force(n_v, 0.0);
    for (size_t i = 0; i < n_c; i++) {
        auto& contact = contacts[i];
        if (contact.jacobian.size() != n_v) {
            throw std::invalid_argument("Generalized SDF-NCP contact jacobian has inconsistent size.");
        }
        const double lambda = z[n_v + i];
        contact.lambda_n = lambda;
        contact.penetration = std::max(0.0, -contact.gap);
        const double scaled_gap = contact.gap / problem.gap_scale;
        const double scaled_lambda = lambda / problem.lambda_scale;
        contact.ncp_residual = SmoothFischerBurmeister(scaled_gap, scaled_lambda, problem.eps);
        contact.complementarity_error = ComplementarityError(scaled_gap, scaled_lambda);
        for (size_t j = 0; j < n_v; j++) {
            contact_force[j] += contact.jacobian[j] * lambda * contact.weight;
        }
    }

    ChSdfNcpGeneralizedResidual residual;
    residual.contacts = std::move(contacts);
    residual.value.assign(n_v + n_c, 0.0);
    for (size_t j = 0; j < n_v; j++) {
        residual.value[j] = problem.mass_diagonal[j] * (v_next[j] - problem.current_velocity[j]) -
                            problem.dt * (problem.external_force[j] + contact_force[j]);
    }
    for (size_t i = 0; i < n_c; i++) {
        residual.value[n_v + i] = residual.contacts[i].ncp_residual;
    }
    return residual;
}

inline std::vector<std::vector<double>> FiniteDifferenceSdfNcpGeneralizedJacobian(
    const ChSdfNcpGeneralizedProblem& problem,
    const std::vector<double>& z,
    double h = 1.0e-7) {
    const size_t n = z.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));
    for (size_t col = 0; col < n; col++) {
        std::vector<double> zp = z;
        std::vector<double> zm = z;
        zp[col] += h;
        zm[col] -= h;
        const auto rp = ComputeSdfNcpGeneralizedResidual(problem, zp).value;
        const auto rm = ComputeSdfNcpGeneralizedResidual(problem, zm).value;
        for (size_t row = 0; row < n; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }
    return jac;
}

inline std::vector<std::vector<double>> AnalyticSdfNcpGeneralizedJacobian(
    const ChSdfNcpGeneralizedProblem& problem,
    const std::vector<double>& z) {
    const size_t n_v = problem.current_velocity.size();
    const size_t n_c = problem.contact_count;
    const size_t n = n_v + n_c;
    if (z.size() != n) {
        throw std::invalid_argument("Generalized SDF-NCP analytic Jacobian has inconsistent dimension.");
    }

    const auto residual = ComputeSdfNcpGeneralizedResidual(problem, z);
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));

    bool has_contact_jacobian_derivatives = true;
    for (const auto& contact : residual.contacts) {
        if (contact.jacobian_velocity_derivative.size() != n_v * n_v) {
            has_contact_jacobian_derivatives = false;
            break;
        }
    }

    if (has_contact_jacobian_derivatives) {
        for (size_t row = 0; row < n_v; row++) {
            jac[row][row] = problem.mass_diagonal[row];
        }
        for (size_t c = 0; c < n_c; c++) {
            const auto& contact = residual.contacts[c];
            const double lambda = z[n_v + c];
            for (size_t row = 0; row < n_v; row++) {
                for (size_t col = 0; col < n_v; col++) {
                    jac[row][col] -= problem.dt * lambda * contact.weight *
                                     contact.jacobian_velocity_derivative[row * n_v + col];
                }
            }
        }
    } else {
        // Fallback: keep the velocity block as finite differences so SDF query variations are retained.
        constexpr double h = 1.0e-7;
        for (size_t col = 0; col < n_v; col++) {
            std::vector<double> zp = z;
            std::vector<double> zm = z;
            zp[col] += h;
            zm[col] -= h;
            const auto rp = ComputeSdfNcpGeneralizedResidual(problem, zp).value;
            const auto rm = ComputeSdfNcpGeneralizedResidual(problem, zm).value;
            for (size_t row = 0; row < n; row++) {
                jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
            }
        }
    }

    for (size_t c = 0; c < n_c; c++) {
        const auto& contact = residual.contacts[c];
        const size_t lambda_col = n_v + c;
        const size_t n_j = std::min(n_v, contact.jacobian.size());

        for (size_t j = 0; j < n_j; j++) {
            jac[j][lambda_col] += -problem.dt * contact.jacobian[j] * contact.weight;
        }

        const double lambda = z[lambda_col];
        const double scaled_gap = contact.gap / problem.gap_scale;
        const double scaled_lambda = lambda / problem.lambda_scale;
        const auto fb_grad = SmoothFischerBurmeisterGrad(scaled_gap, scaled_lambda, problem.eps);
        const size_t ncp_row = n_v + c;
        for (size_t j = 0; j < n_j; j++) {
            jac[ncp_row][j] = fb_grad.dPhi_dg * problem.dt * contact.jacobian[j] / problem.gap_scale;
        }
        jac[ncp_row][lambda_col] = fb_grad.dPhi_dlambda / problem.lambda_scale;
    }

    return jac;
}

struct ChSdfNcpAdScalar {
    double value = 0.0;
    std::vector<double> derivative;

    ChSdfNcpAdScalar() = default;

    ChSdfNcpAdScalar(double v, size_t n = 0) : value(v), derivative(n, 0.0) {}

    static ChSdfNcpAdScalar Variable(double v, size_t index, size_t n) {
        ChSdfNcpAdScalar out(v, n);
        out.derivative.at(index) = 1.0;
        return out;
    }
};

inline ChSdfNcpAdScalar operator+(const ChSdfNcpAdScalar& a, const ChSdfNcpAdScalar& b) {
    ChSdfNcpAdScalar out(a.value + b.value, a.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] + b.derivative[i];
    }
    return out;
}

inline ChSdfNcpAdScalar operator-(const ChSdfNcpAdScalar& a, const ChSdfNcpAdScalar& b) {
    ChSdfNcpAdScalar out(a.value - b.value, a.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] - b.derivative[i];
    }
    return out;
}

inline ChSdfNcpAdScalar operator*(const ChSdfNcpAdScalar& a, const ChSdfNcpAdScalar& b) {
    ChSdfNcpAdScalar out(a.value * b.value, a.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] * b.value + b.derivative[i] * a.value;
    }
    return out;
}

inline ChSdfNcpAdScalar operator/(const ChSdfNcpAdScalar& a, const ChSdfNcpAdScalar& b) {
    ChSdfNcpAdScalar out(a.value / b.value, a.derivative.size());
    const double denom = b.value * b.value;
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = (a.derivative[i] * b.value - a.value * b.derivative[i]) / denom;
    }
    return out;
}

inline ChSdfNcpAdScalar operator+(const ChSdfNcpAdScalar& a, double b) {
    ChSdfNcpAdScalar out = a;
    out.value += b;
    return out;
}

inline ChSdfNcpAdScalar operator+(double a, const ChSdfNcpAdScalar& b) {
    return b + a;
}

inline ChSdfNcpAdScalar operator-(const ChSdfNcpAdScalar& a, double b) {
    ChSdfNcpAdScalar out = a;
    out.value -= b;
    return out;
}

inline ChSdfNcpAdScalar operator-(double a, const ChSdfNcpAdScalar& b) {
    ChSdfNcpAdScalar out(a - b.value, b.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = -b.derivative[i];
    }
    return out;
}

inline ChSdfNcpAdScalar operator*(const ChSdfNcpAdScalar& a, double b) {
    ChSdfNcpAdScalar out(a.value * b, a.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] * b;
    }
    return out;
}

inline ChSdfNcpAdScalar operator*(double a, const ChSdfNcpAdScalar& b) {
    return b * a;
}

inline ChSdfNcpAdScalar operator/(const ChSdfNcpAdScalar& a, double b) {
    ChSdfNcpAdScalar out(a.value / b, a.derivative.size());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] / b;
    }
    return out;
}

inline ChSdfNcpAdScalar Sqrt(const ChSdfNcpAdScalar& a) {
    const double root = std::sqrt(std::max(a.value, 0.0));
    ChSdfNcpAdScalar out(root, a.derivative.size());
    const double denom = std::max(2.0 * root, std::numeric_limits<double>::min());
    for (size_t i = 0; i < out.derivative.size(); i++) {
        out.derivative[i] = a.derivative[i] / denom;
    }
    return out;
}

inline ChSdfNcpAdScalar SmoothFischerBurmeisterAd(const ChSdfNcpAdScalar& g,
                                                  const ChSdfNcpAdScalar& lambda,
                                                  double eps) {
    const double e = RegularizedSmoothingEps(eps);
    return Sqrt(g * g + lambda * lambda + e * e) - g - lambda;
}

struct ChSdfNcpImpulseMixedSettings {
    double beta = 0.2;
    double cfm = 0.0;
    double velocity_scale = 0.0;
    double impulse_scale = 0.0;
};

struct ChSdfNcpAdResidualJacobian {
    std::vector<double> value;
    std::vector<std::vector<double>> jacobian;
    std::vector<ChSdfNcpGeneralizedContact> contacts;
};

inline ChSdfNcpAdResidualJacobian ComputeSdfNcpGeneralizedImpulseMixedAd(
    const ChSdfNcpGeneralizedProblem& problem,
    const std::vector<double>& z,
    const ChSdfNcpImpulseMixedSettings& settings = ChSdfNcpImpulseMixedSettings()) {
    ValidateSdfNcpGeneralizedProblem(problem);
    const size_t n_v = problem.current_velocity.size();
    const size_t n_c = problem.contact_count;
    const size_t n = n_v + n_c;
    if (z.size() != n) {
        throw std::invalid_argument("Generalized impulse mixed SDF-NCP unknown vector has inconsistent size.");
    }

    std::vector<double> v_next(z.begin(), z.begin() + static_cast<std::ptrdiff_t>(n_v));
    std::vector<ChSdfNcpGeneralizedContact> contacts = problem.evaluate_contacts(v_next);
    if (contacts.size() != n_c) {
        throw std::invalid_argument("Generalized impulse mixed SDF-NCP contact callback returned inconsistent count.");
    }

    const double velocity_scale =
        settings.velocity_scale > 0.0 ? settings.velocity_scale : problem.gap_scale / std::max(problem.dt, 1.0e-12);
    const double impulse_scale = settings.impulse_scale > 0.0 ? settings.impulse_scale : problem.lambda_scale;
    if (!(velocity_scale > 0.0) || !std::isfinite(velocity_scale) || !(impulse_scale > 0.0) ||
        !std::isfinite(impulse_scale)) {
        throw std::invalid_argument("Generalized impulse mixed SDF-NCP scales must be finite and positive.");
    }

    std::vector<ChSdfNcpAdScalar> v_ad;
    std::vector<ChSdfNcpAdScalar> p_ad;
    v_ad.reserve(n_v);
    p_ad.reserve(n_c);
    for (size_t j = 0; j < n_v; j++) {
        v_ad.push_back(ChSdfNcpAdScalar::Variable(z[j], j, n));
    }
    for (size_t i = 0; i < n_c; i++) {
        p_ad.push_back(ChSdfNcpAdScalar::Variable(z[n_v + i], n_v + i, n));
    }

    std::vector<ChSdfNcpAdScalar> residual_ad(n, ChSdfNcpAdScalar(0.0, n));
    for (size_t j = 0; j < n_v; j++) {
        residual_ad[j] = problem.mass_diagonal[j] * (v_ad[j] - problem.current_velocity[j]) -
                         problem.dt * problem.external_force[j];
    }

    for (size_t i = 0; i < n_c; i++) {
        auto& contact = contacts[i];
        if (contact.jacobian.size() != n_v) {
            throw std::invalid_argument("Generalized impulse mixed SDF-NCP contact jacobian has inconsistent size.");
        }

        contact.lambda_n = z[n_v + i];
        contact.penetration = std::max(0.0, -contact.gap);

        ChSdfNcpAdScalar normal_velocity(contact.normal_velocity_offset, n);
        ChSdfNcpAdScalar linearized_gap(contact.gap, n);
        for (size_t j = 0; j < n_v; j++) {
            normal_velocity = normal_velocity + contact.jacobian[j] * v_ad[j];
            linearized_gap = linearized_gap + problem.dt * contact.jacobian[j] * (v_ad[j] - v_next[j]);
            residual_ad[j] = residual_ad[j] - contact.jacobian[j] * p_ad[i] * contact.weight;
        }

        const ChSdfNcpAdScalar mixed_velocity =
            normal_velocity + settings.beta * linearized_gap / std::max(problem.dt, 1.0e-12) + settings.cfm * p_ad[i];
        const ChSdfNcpAdScalar ncp =
            SmoothFischerBurmeisterAd(mixed_velocity / velocity_scale, p_ad[i] / impulse_scale, problem.eps);
        residual_ad[n_v + i] = ncp;
        contact.ncp_residual = ncp.value;
        contact.complementarity_error = ComplementarityError(mixed_velocity.value / velocity_scale,
                                                             contact.lambda_n / impulse_scale);
    }

    ChSdfNcpAdResidualJacobian out;
    out.value.assign(n, 0.0);
    out.jacobian.assign(n, std::vector<double>(n, 0.0));
    for (size_t row = 0; row < n; row++) {
        out.value[row] = residual_ad[row].value;
        for (size_t col = 0; col < n; col++) {
            out.jacobian[row][col] = residual_ad[row].derivative[col];
        }
    }
    out.contacts = std::move(contacts);
    return out;
}

inline ChSdfNcpGeneralizedDiagnostics MakeSdfNcpGeneralizedImpulseMixedDiagnostics(
    const ChSdfNcpGeneralizedProblem& problem,
    const std::vector<double>& z,
    const ChSdfNcpImpulseMixedSettings& settings,
    bool success,
    int iterations,
    const std::string& message) {
    const auto residual = ComputeSdfNcpGeneralizedImpulseMixedAd(problem, z, settings);
    ChSdfNcpGeneralizedDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.message = message;
    diagnostics.residual_norm = NormDynamic(residual.value);
    diagnostics.min_lambda_n = residual.contacts.empty() ? 0.0 : std::numeric_limits<double>::max();
    double sum_comp = 0.0;
    for (const auto& contact : residual.contacts) {
        diagnostics.max_penetration = std::max(diagnostics.max_penetration, contact.penetration);
        diagnostics.min_lambda_n = std::min(diagnostics.min_lambda_n, contact.lambda_n);
        diagnostics.max_lambda_n = std::max(diagnostics.max_lambda_n, std::abs(contact.lambda_n));
        diagnostics.max_complementarity_error =
            std::max(diagnostics.max_complementarity_error, contact.complementarity_error);
        diagnostics.ncp_residual_norm += contact.ncp_residual * contact.ncp_residual;
        sum_comp += contact.complementarity_error;
    }
    diagnostics.ncp_residual_norm = std::sqrt(diagnostics.ncp_residual_norm);
    if (!residual.contacts.empty()) {
        diagnostics.mean_complementarity_error = sum_comp / static_cast<double>(residual.contacts.size());
    }
    return diagnostics;
}

inline ChSdfNcpGeneralizedStepResult SolveSdfNcpGeneralizedImpulseMixedProblem(
    const ChSdfNcpGeneralizedProblem& problem,
    std::vector<double> z,
    const ChSdfNcpImpulseMixedSettings& settings = ChSdfNcpImpulseMixedSettings()) {
    ValidateSdfNcpGeneralizedProblem(problem);
    const size_t n_v = problem.current_velocity.size();
    const size_t n_c = problem.contact_count;
    if (z.empty()) {
        z.assign(n_v + n_c, 0.0);
        for (size_t j = 0; j < n_v; j++) {
            z[j] = problem.current_velocity[j] + problem.dt * problem.external_force[j] / problem.mass_diagonal[j];
        }
    }
    if (z.size() != n_v + n_c) {
        throw std::invalid_argument("Generalized impulse mixed SDF-NCP initial guess has inconsistent size.");
    }
    for (size_t i = 0; i < n_c; i++) {
        if (z[n_v + i] < 0.0) {
            z[n_v + i] = 0.0;
        }
    }

    double residual_norm = NormDynamic(ComputeSdfNcpGeneralizedImpulseMixedAd(problem, z, settings).value);
    bool success = residual_norm <= problem.tolerance;
    int iterations = 0;
    std::string message = success ? "Initial guess satisfies impulse mixed residual tolerance." : "Newton did not run.";

    for (iterations = 0; !success && iterations < problem.max_iterations; iterations++) {
        const auto residual_jacobian = ComputeSdfNcpGeneralizedImpulseMixedAd(problem, z, settings);
        std::vector<double> rhs(residual_jacobian.value.size(), 0.0);
        for (size_t i = 0; i < residual_jacobian.value.size(); i++) {
            rhs[i] = -residual_jacobian.value[i];
        }

        std::vector<double> dz(z.size(), 0.0);
        if (!SolveLinearDynamic(residual_jacobian.jacobian, rhs, dz)) {
            message = "Newton impulse mixed AD Jacobian was singular.";
            break;
        }

        double alpha = 1.0;
        std::vector<double> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;
        for (int ls = 0; ls < 18; ls++) {
            for (size_t i = 0; i < z.size(); i++) {
                candidate[i] = z[i] + alpha * dz[i];
            }
            for (size_t i = 0; i < n_c; i++) {
                double& impulse = candidate[n_v + i];
                if (impulse < 0.0) {
                    impulse = 0.0;
                }
            }

            candidate_norm = NormDynamic(ComputeSdfNcpGeneralizedImpulseMixedAd(problem, candidate, settings).value);
            if (std::isfinite(candidate_norm) && (candidate_norm <= residual_norm || alpha <= 1.0e-5)) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            message = "Newton impulse mixed line search failed.";
            break;
        }

        z = candidate;
        residual_norm = candidate_norm;
        if (residual_norm <= problem.tolerance) {
            success = true;
            iterations++;
            message = "Newton impulse mixed converged.";
            break;
        }
        if (NormDynamic(dz) * alpha <= 1.0e-12 * std::max(1.0, NormDynamic(z))) {
            message = "Newton impulse mixed step reached small update tolerance.";
            break;
        }
    }

    if (!success && residual_norm <= std::max(problem.tolerance, problem.relaxed_tolerance)) {
        success = true;
        message = "Newton impulse mixed reached relaxed residual tolerance.";
    }
    for (size_t i = 0; i < n_c; i++) {
        double& impulse = z[n_v + i];
        if (impulse < -problem.negative_lambda_tolerance) {
            success = false;
            message = "Newton impulse mixed returned a significantly negative normal impulse.";
        } else if (impulse < 0.0) {
            impulse = 0.0;
        }
    }

    const auto final_residual = ComputeSdfNcpGeneralizedImpulseMixedAd(problem, z, settings);
    ChSdfNcpGeneralizedStepResult result;
    result.next_velocity.assign(z.begin(), z.begin() + static_cast<std::ptrdiff_t>(n_v));
    result.lambdas.assign(z.begin() + static_cast<std::ptrdiff_t>(n_v), z.end());
    result.contacts = final_residual.contacts;
    result.diagnostics =
        MakeSdfNcpGeneralizedImpulseMixedDiagnostics(problem, z, settings, success, iterations, message);
    return result;
}

inline ChSdfNcpGeneralizedDiagnostics MakeSdfNcpGeneralizedDiagnostics(
    const ChSdfNcpGeneralizedProblem& problem,
    const std::vector<double>& z,
    bool success,
    int iterations,
    const std::string& message) {
    const auto residual = ComputeSdfNcpGeneralizedResidual(problem, z);
    ChSdfNcpGeneralizedDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.message = message;
    diagnostics.residual_norm = NormDynamic(residual.value);
    diagnostics.min_lambda_n = residual.contacts.empty() ? 0.0 : std::numeric_limits<double>::max();
    double sum_comp = 0.0;
    for (const auto& contact : residual.contacts) {
        diagnostics.max_penetration = std::max(diagnostics.max_penetration, contact.penetration);
        diagnostics.min_lambda_n = std::min(diagnostics.min_lambda_n, contact.lambda_n);
        diagnostics.max_lambda_n = std::max(diagnostics.max_lambda_n, std::abs(contact.lambda_n));
        diagnostics.max_complementarity_error =
            std::max(diagnostics.max_complementarity_error, contact.complementarity_error);
        diagnostics.ncp_residual_norm += contact.ncp_residual * contact.ncp_residual;
        sum_comp += contact.complementarity_error;
    }
    diagnostics.ncp_residual_norm = std::sqrt(diagnostics.ncp_residual_norm);
    if (!residual.contacts.empty()) {
        diagnostics.mean_complementarity_error = sum_comp / static_cast<double>(residual.contacts.size());
    }
    return diagnostics;
}

inline ChSdfNcpGeneralizedStepResult SolveSdfNcpGeneralizedProblem(const ChSdfNcpGeneralizedProblem& problem,
                                                                   std::vector<double> z) {
    ValidateSdfNcpGeneralizedProblem(problem);
    const size_t n_v = problem.current_velocity.size();
    const size_t n_c = problem.contact_count;
    if (z.empty()) {
        z.assign(n_v + n_c, 0.0);
        for (size_t j = 0; j < n_v; j++) {
            z[j] = problem.current_velocity[j] + problem.dt * problem.external_force[j] / problem.mass_diagonal[j];
        }
    }
    if (z.size() != n_v + n_c) {
        throw std::invalid_argument("Generalized SDF-NCP initial guess has inconsistent size.");
    }

    double residual_norm = NormDynamic(ComputeSdfNcpGeneralizedResidual(problem, z).value);
    bool success = residual_norm <= problem.tolerance;
    int iterations = 0;
    std::string message = success ? "Initial guess satisfies residual tolerance." : "Newton did not run.";
    for (iterations = 0; !success && iterations < problem.max_iterations; iterations++) {
        const auto residual = ComputeSdfNcpGeneralizedResidual(problem, z).value;
        const auto jac = problem.use_analytic_jacobian ? AnalyticSdfNcpGeneralizedJacobian(problem, z) :
                                                         FiniteDifferenceSdfNcpGeneralizedJacobian(problem, z);
        std::vector<double> rhs(residual.size(), 0.0);
        for (size_t i = 0; i < residual.size(); i++) {
            rhs[i] = -residual[i];
        }
        std::vector<double> dz(z.size(), 0.0);
        if (!SolveLinearDynamic(jac, rhs, dz)) {
            message = problem.use_analytic_jacobian ? "Newton analytic Jacobian was singular." :
                                                      "Newton finite-difference Jacobian was singular.";
            break;
        }
        double alpha = 1.0;
        std::vector<double> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;
        for (int ls = 0; ls < 14; ls++) {
            for (size_t i = 0; i < z.size(); i++) {
                candidate[i] = z[i] + alpha * dz[i];
            }
            candidate_norm = NormDynamic(ComputeSdfNcpGeneralizedResidual(problem, candidate).value);
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
        if (residual_norm <= problem.tolerance) {
            success = true;
            iterations++;
            message = "Newton converged.";
            break;
        }
    }

    if (!success && residual_norm <= std::max(problem.tolerance, problem.relaxed_tolerance)) {
        success = true;
        message = "Newton reached relaxed residual tolerance.";
    }
    for (size_t i = 0; i < n_c; i++) {
        double& lambda = z[n_v + i];
        if (lambda < -problem.negative_lambda_tolerance) {
            success = false;
            message = "Newton returned a significantly negative normal multiplier.";
        } else if (lambda < 0.0) {
            lambda = 0.0;
        }
    }

    const auto final_residual = ComputeSdfNcpGeneralizedResidual(problem, z);
    ChSdfNcpGeneralizedStepResult result;
    result.next_velocity.assign(z.begin(), z.begin() + static_cast<std::ptrdiff_t>(n_v));
    result.lambdas.assign(z.begin() + static_cast<std::ptrdiff_t>(n_v), z.end());
    result.contacts = final_residual.contacts;
    result.diagnostics = MakeSdfNcpGeneralizedDiagnostics(problem, z, success, iterations, message);
    return result;
}

inline ChSdfNcpPointMassContactSetDiagnostics MakeSdfNcpPointMassContactSetDiagnostics(
    const ChSdfNcpContactSet& contact_set,
    const ChSdfNcpPointMassState& state,
    const ChSdfNcpPointMassSettings& settings,
    const std::vector<double>& z,
    bool success,
    int iterations,
    const std::string& message) {
    const auto residual = ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, z);
    ChSdfNcpPointMassContactSetDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.residual_norm = NormDynamic(residual.value);
    diagnostics.max_penetration = residual.contact_set.MaxPenetration();
    diagnostics.max_complementarity_error = residual.contact_set.MaxComplementarityError();
    diagnostics.mean_complementarity_error = residual.contact_set.MeanComplementarityError();
    diagnostics.ncp_residual_norm = residual.contact_set.NcpResidualNorm();

    diagnostics.min_lambda_n = 0.0;
    for (size_t i = 0; i < residual.contact_set.GetNumContacts(); i++) {
        diagnostics.min_lambda_n = std::min(diagnostics.min_lambda_n, residual.contact_set.GetContact(i).lambda_n);
    }

    diagnostics.message = message;
    return diagnostics;
}

inline ChSdfNcpPointMassContactSetStepResult SolveSdfNcpPointMassContactSetStep(
    const ChSdfNcpContactSet& contact_set,
    const ChSdfNcpPointMassState& state,
    const ChSdfNcpPointMassSettings& settings) {
    ChSdfNcpPointMassContactSetStepResult result;
    result.state = state;
    result.contact_set = contact_set;

    const size_t n_contacts = contact_set.GetNumContacts();
    if (settings.mass <= 0.0 || settings.dt <= 0.0 || settings.max_iterations <= 0) {
        std::vector<double> invalid_z(3 + n_contacts, 0.0);
        invalid_z[0] = state.velocity.x();
        invalid_z[1] = state.velocity.y();
        invalid_z[2] = state.velocity.z();
        result.diagnostics = MakeSdfNcpPointMassContactSetDiagnostics(
            contact_set, state, settings, invalid_z, false, 0, "Invalid SDF-NCP contact-set settings.");
        return result;
    }

    const ChVector3d external_accel(0.0, -settings.gravity, 0.0);
    const ChVector3d v_guess = state.velocity + external_accel * settings.dt;
    const ChVector3d q_guess = state.position + v_guess * settings.dt;
    const double lambda_scale = settings.mass / std::max(settings.dt * settings.dt, 1.0e-12);

    std::vector<double> z(3 + n_contacts, 0.0);
    z[0] = v_guess.x();
    z[1] = v_guess.y();
    z[2] = v_guess.z();
    for (size_t i = 0; i < n_contacts; i++) {
        const auto& point = contact_set.GetContact(i);
        const double gap_guess = point.surface ? point.surface->Phi(q_guess) : 0.0;
        z[3 + i] = std::max(point.lambda_n, std::max(0.0, -gap_guess) * lambda_scale);
    }

    double residual_norm = NormDynamic(ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, z).value);

    int iterations = 0;
    bool success = residual_norm <= settings.tolerance;
    std::string message = success ? "Initial guess satisfies residual tolerance." : "Newton did not run.";

    for (iterations = 0; !success && iterations < settings.max_iterations; iterations++) {
        const auto residual = ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, z).value;
        const auto jac = FiniteDifferenceSdfNcpPointMassContactSetJacobian(contact_set, state, settings, z);

        std::vector<double> rhs(residual.size(), 0.0);
        for (size_t i = 0; i < residual.size(); i++) {
            rhs[i] = -residual[i];
        }

        std::vector<double> dz(z.size(), 0.0);
        if (!SolveLinearDynamic(jac, rhs, dz)) {
            message = "Newton finite-difference Jacobian was singular.";
            break;
        }

        double alpha = 1.0;
        std::vector<double> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;

        for (int ls = 0; ls < 12; ls++) {
            for (size_t i = 0; i < z.size(); i++) {
                candidate[i] = z[i] + alpha * dz[i];
            }
            candidate_norm =
                NormDynamic(ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, candidate).value);
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

    for (size_t i = 0; i < n_contacts; i++) {
        double& lambda = z[3 + i];
        if (lambda < -1.0e-8) {
            success = false;
            message = "Newton returned a significantly negative normal multiplier.";
        } else if (lambda < 0.0) {
            lambda = 0.0;
        }
    }

    const ChVector3d v_next(z[0], z[1], z[2]);
    result.state.velocity = v_next;
    result.state.position = state.position + v_next * settings.dt;

    const auto final_residual = ComputeSdfNcpPointMassContactSetResidual(contact_set, state, settings, z);
    result.contact_set = final_residual.contact_set;
    result.diagnostics = MakeSdfNcpPointMassContactSetDiagnostics(
        contact_set, state, settings, z, success, iterations, message);
    return result;
}

struct ChSdfRigidBody2DState {
    double x = 0.0;
    double y = 0.0;
    double theta = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double omega = 0.0;
};

struct ChSdfRigidBody2DSettings {
    double mass = 1.0;
    double inertia_z = 1.0;
    double gravity = 9.81;
    double dt = 1.0e-3;
    double eps = 1.0e-6;
    double tolerance = 1.0e-10;
    int max_iterations = 40;
};

struct ChSdfRigidBody2DPointJacobian {
    ChVector3d dx_dq = ChVector3d(1, 0, 0);
    ChVector3d dy_dq = ChVector3d(0, 1, 0);
};

inline ChVector3d RigidBody2DCoordinates(const ChSdfRigidBody2DState& state) {
    return ChVector3d(state.x, state.y, state.theta);
}

inline ChVector3d RigidBody2DVelocity(const ChSdfRigidBody2DState& state) {
    return ChVector3d(state.vx, state.vy, state.omega);
}

inline ChVector3d RigidBody2DLocalPointToWorld(const ChVector3d& q, const ChVector3d& r_local) {
    const double c = std::cos(q.z());
    const double s = std::sin(q.z());
    const double x_world = q.x() + c * r_local.x() - s * r_local.y();
    const double y_world = q.y() + s * r_local.x() + c * r_local.y();
    return ChVector3d(x_world, y_world, 0.0);
}

inline ChSdfRigidBody2DPointJacobian RigidBody2DLocalPointJacobian(const ChVector3d& q,
                                                                   const ChVector3d& r_local) {
    const double c = std::cos(q.z());
    const double s = std::sin(q.z());

    ChSdfRigidBody2DPointJacobian jacobian;
    jacobian.dx_dq = ChVector3d(1.0, 0.0, -s * r_local.x() - c * r_local.y());
    jacobian.dy_dq = ChVector3d(0.0, 1.0, c * r_local.x() - s * r_local.y());
    return jacobian;
}

inline ChVector3d RigidBody2DSdfGapJacobian(const ChVector3d& normal,
                                            const ChSdfRigidBody2DPointJacobian& point_jacobian) {
    return point_jacobian.dx_dq * normal.x() + point_jacobian.dy_dq * normal.y();
}

struct ChSdfRigidBody2DContactPoint {
    std::shared_ptr<const ChSdfSurface> surface;
    ChVector3d local_point = ChVector3d(0, 0, 0);
    ChVector3d world_point = ChVector3d(0, 0, 0);
    double lambda_n = 0.0;
    double gap = 0.0;
    double penetration = 0.0;
    ChVector3d normal = ChVector3d(0, 1, 0);
    ChVector3d jacobian = ChVector3d(0, 1, 0);
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
};

inline ChSdfRigidBody2DContactPoint EvaluateRigidBody2DSdfContact(const ChSdfSurface& surface,
                                                                  const ChVector3d& q,
                                                                  const ChVector3d& r_local,
                                                                  double lambda,
                                                                  double eps) {
    const ChVector3d local_point(r_local.x(), r_local.y(), 0.0);
    const ChVector3d world_point = RigidBody2DLocalPointToWorld(q, local_point);
    const ChVector3d normal = NormalizeSdfVector(surface.Grad(world_point), "SDF gradient must be nonzero.");
    const auto point_jacobian = RigidBody2DLocalPointJacobian(q, local_point);
    const ChVector3d sdf_jacobian = RigidBody2DSdfGapJacobian(normal, point_jacobian);
    const double gap = surface.Phi(world_point);

    ChSdfRigidBody2DContactPoint contact;
    contact.local_point = local_point;
    contact.world_point = world_point;
    contact.lambda_n = lambda;
    contact.gap = gap;
    contact.penetration = std::max(0.0, -gap);
    contact.normal = normal;
    contact.jacobian = sdf_jacobian;
    contact.ncp_residual = SmoothFischerBurmeister(gap, lambda, eps);
    contact.complementarity_error = ComplementarityError(gap, lambda);
    return contact;
}

class ChSdfNcpRigidBody2DContactSet {
  public:
    void AddContact(std::shared_ptr<const ChSdfSurface> surface,
                    const ChVector3d& local_point,
                    double lambda = 0.0) {
        if (!surface) {
            throw std::invalid_argument("ChSdfNcpRigidBody2DContactSet contact surface must not be null.");
        }
        ChSdfRigidBody2DContactPoint point;
        point.surface = std::move(surface);
        point.local_point = ChVector3d(local_point.x(), local_point.y(), 0.0);
        point.lambda_n = lambda;
        m_points.push_back(point);
    }

    size_t GetNumContacts() const {
        return m_points.size();
    }

    bool Empty() const {
        return m_points.empty();
    }

    const ChSdfRigidBody2DContactPoint& GetContact(size_t index) const {
        return m_points.at(index);
    }

    ChSdfRigidBody2DContactPoint& GetContact(size_t index) {
        return m_points.at(index);
    }

    void SetLambda(size_t index, double lambda) {
        m_points.at(index).lambda_n = lambda;
    }

    double GetLambda(size_t index) const {
        return m_points.at(index).lambda_n;
    }

    void Update(const ChVector3d& q, double eps) {
        for (auto& point : m_points) {
            if (!point.surface) {
                throw std::invalid_argument("ChSdfNcpRigidBody2DContactSet contains a null surface.");
            }
            auto evaluated =
                EvaluateRigidBody2DSdfContact(*point.surface, q, point.local_point, point.lambda_n, eps);
            evaluated.surface = point.surface;
            point = evaluated;
        }
    }

    ChVector3d ComputeGeneralizedContactForce() const {
        ChVector3d generalized_force(0, 0, 0);
        for (const auto& point : m_points) {
            generalized_force += point.jacobian * point.lambda_n;
        }
        return generalized_force;
    }

    double MaxPenetration() const {
        double value = 0.0;
        for (const auto& point : m_points) {
            value = std::max(value, point.penetration);
        }
        return value;
    }

    double MaxComplementarityError() const {
        double value = 0.0;
        for (const auto& point : m_points) {
            value = std::max(value, point.complementarity_error);
        }
        return value;
    }

    double MeanComplementarityError() const {
        if (m_points.empty()) {
            return 0.0;
        }
        double sum = 0.0;
        for (const auto& point : m_points) {
            sum += point.complementarity_error;
        }
        return sum / static_cast<double>(m_points.size());
    }

    double NcpResidualNorm() const {
        double sum_sq = 0.0;
        for (const auto& point : m_points) {
            sum_sq += point.ncp_residual * point.ncp_residual;
        }
        return std::sqrt(sum_sq);
    }

  private:
    std::vector<ChSdfRigidBody2DContactPoint> m_points;
};

struct ChSdfNcpRigidBody2DResidual {
    std::vector<double> value;
    ChSdfNcpRigidBody2DContactSet contact_set;
};

struct ChSdfNcpRigidBody2DDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double max_penetration = 0.0;
    double min_lambda_n = 0.0;
    double max_complementarity_error = 0.0;
    double mean_complementarity_error = 0.0;
    double ncp_residual_norm = 0.0;
    std::vector<double> lambdas;
    std::vector<double> gaps;
    std::string message;
};

struct ChSdfNcpRigidBody2DStepResult {
    ChSdfRigidBody2DState state;
    ChSdfNcpRigidBody2DContactSet contact_set;
    ChSdfNcpRigidBody2DDiagnostics diagnostics;
};

inline ChSdfNcpRigidBody2DResidual ComputeSdfNcpRigidBody2DResidual(
    const ChSdfNcpRigidBody2DContactSet& contact_set,
    const ChSdfRigidBody2DState& state,
    const ChSdfRigidBody2DSettings& settings,
    const std::vector<double>& z) {
    const size_t n_contacts = contact_set.GetNumContacts();
    if (z.size() != 3 + n_contacts) {
        throw std::invalid_argument("SDF-NCP rigid-body 2D unknown vector has inconsistent size.");
    }

    const ChVector3d v_next(z[0], z[1], z[2]);
    const ChVector3d q_next(state.x + settings.dt * v_next.x(),
                            state.y + settings.dt * v_next.y(),
                            state.theta + settings.dt * v_next.z());
    const ChVector3d external_force(0.0, -settings.mass * settings.gravity, 0.0);

    ChSdfNcpRigidBody2DResidual residual;
    residual.contact_set = contact_set;
    for (size_t i = 0; i < n_contacts; i++) {
        residual.contact_set.SetLambda(i, z[3 + i]);
    }
    residual.contact_set.Update(q_next, settings.eps);

    const ChVector3d contact_force = residual.contact_set.ComputeGeneralizedContactForce();
    const ChVector3d v_current = RigidBody2DVelocity(state);
    const ChVector3d inertia_weighted_delta(settings.mass * (v_next.x() - v_current.x()),
                                            settings.mass * (v_next.y() - v_current.y()),
                                            settings.inertia_z * (v_next.z() - v_current.z()));
    const ChVector3d r_v = inertia_weighted_delta - settings.dt * (external_force + contact_force);

    residual.value.assign(3 + n_contacts, 0.0);
    residual.value[0] = r_v.x();
    residual.value[1] = r_v.y();
    residual.value[2] = r_v.z();
    for (size_t i = 0; i < n_contacts; i++) {
        residual.value[3 + i] = residual.contact_set.GetContact(i).ncp_residual;
    }

    return residual;
}

inline std::vector<std::vector<double>> FiniteDifferenceSdfNcpRigidBody2DJacobian(
    const ChSdfNcpRigidBody2DContactSet& contact_set,
    const ChSdfRigidBody2DState& state,
    const ChSdfRigidBody2DSettings& settings,
    const std::vector<double>& z,
    double h = 1.0e-7) {
    const size_t n = z.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));

    for (size_t col = 0; col < n; col++) {
        std::vector<double> zp = z;
        std::vector<double> zm = z;
        zp[col] += h;
        zm[col] -= h;

        const auto rp = ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, zp).value;
        const auto rm = ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, zm).value;

        for (size_t row = 0; row < n; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }

    return jac;
}

inline ChSdfNcpRigidBody2DDiagnostics MakeSdfNcpRigidBody2DDiagnostics(
    const ChSdfNcpRigidBody2DContactSet& contact_set,
    const ChSdfRigidBody2DState& state,
    const ChSdfRigidBody2DSettings& settings,
    const std::vector<double>& z,
    bool success,
    int iterations,
    const std::string& message) {
    const auto residual = ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, z);

    ChSdfNcpRigidBody2DDiagnostics diagnostics;
    diagnostics.success = success;
    diagnostics.iterations = iterations;
    diagnostics.residual_norm = NormDynamic(residual.value);
    diagnostics.max_penetration = residual.contact_set.MaxPenetration();
    diagnostics.max_complementarity_error = residual.contact_set.MaxComplementarityError();
    diagnostics.mean_complementarity_error = residual.contact_set.MeanComplementarityError();
    diagnostics.ncp_residual_norm = residual.contact_set.NcpResidualNorm();
    diagnostics.min_lambda_n = 0.0;

    if (residual.contact_set.GetNumContacts() > 0) {
        diagnostics.min_lambda_n = std::numeric_limits<double>::infinity();
    }

    for (size_t i = 0; i < residual.contact_set.GetNumContacts(); i++) {
        const auto& contact = residual.contact_set.GetContact(i);
        diagnostics.min_lambda_n = std::min(diagnostics.min_lambda_n, contact.lambda_n);
        diagnostics.lambdas.push_back(contact.lambda_n);
        diagnostics.gaps.push_back(contact.gap);
    }

    diagnostics.message = message;
    return diagnostics;
}

inline ChSdfNcpRigidBody2DStepResult SolveSdfNcpRigidBody2DStep(
    const ChSdfNcpRigidBody2DContactSet& contact_set,
    const ChSdfRigidBody2DState& state,
    const ChSdfRigidBody2DSettings& settings) {
    ChSdfNcpRigidBody2DStepResult result;
    result.state = state;
    result.contact_set = contact_set;

    const size_t n_contacts = contact_set.GetNumContacts();
    if (settings.mass <= 0.0 || settings.inertia_z <= 0.0 || settings.dt <= 0.0 ||
        settings.max_iterations <= 0) {
        std::vector<double> invalid_z(3 + n_contacts, 0.0);
        invalid_z[0] = state.vx;
        invalid_z[1] = state.vy;
        invalid_z[2] = state.omega;
        result.diagnostics = MakeSdfNcpRigidBody2DDiagnostics(
            contact_set, state, settings, invalid_z, false, 0, "Invalid SDF-NCP rigid-body 2D settings.");
        return result;
    }

    const ChVector3d v_guess(state.vx, state.vy - settings.gravity * settings.dt, state.omega);
    const ChVector3d q_guess(state.x + settings.dt * v_guess.x(),
                             state.y + settings.dt * v_guess.y(),
                             state.theta + settings.dt * v_guess.z());
    const double lambda_scale = settings.mass / std::max(settings.dt * settings.dt, 1.0e-12);

    std::vector<double> z(3 + n_contacts, 0.0);
    z[0] = v_guess.x();
    z[1] = v_guess.y();
    z[2] = v_guess.z();
    for (size_t i = 0; i < n_contacts; i++) {
        const auto& contact = contact_set.GetContact(i);
        double gap_guess = 0.0;
        if (contact.surface) {
            gap_guess = contact.surface->Phi(RigidBody2DLocalPointToWorld(q_guess, contact.local_point));
        }
        z[3 + i] = std::max(contact.lambda_n, std::max(0.0, -gap_guess) * lambda_scale);
    }

    double residual_norm = NormDynamic(ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, z).value);

    int iterations = 0;
    bool success = residual_norm <= settings.tolerance;
    std::string message = success ? "Initial guess satisfies residual tolerance." : "Newton did not run.";

    for (iterations = 0; !success && iterations < settings.max_iterations; iterations++) {
        const auto residual = ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, z).value;
        const auto jac = FiniteDifferenceSdfNcpRigidBody2DJacobian(contact_set, state, settings, z);

        std::vector<double> rhs(residual.size(), 0.0);
        for (size_t i = 0; i < residual.size(); i++) {
            rhs[i] = -residual[i];
        }

        std::vector<double> dz(z.size(), 0.0);
        if (!SolveLinearDynamic(jac, rhs, dz)) {
            message = "Newton finite-difference Jacobian was singular.";
            break;
        }

        double alpha = 1.0;
        std::vector<double> candidate = z;
        double candidate_norm = residual_norm;
        bool accepted = false;

        for (int ls = 0; ls < 12; ls++) {
            for (size_t i = 0; i < z.size(); i++) {
                candidate[i] = z[i] + alpha * dz[i];
            }
            candidate_norm = NormDynamic(ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, candidate).value);
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

    for (size_t i = 0; i < n_contacts; i++) {
        double& lambda = z[3 + i];
        if (lambda < -1.0e-8) {
            success = false;
            message = "Newton returned a significantly negative normal multiplier.";
        } else if (lambda < 0.0) {
            lambda = 0.0;
        }
    }

    result.state.x = state.x + settings.dt * z[0];
    result.state.y = state.y + settings.dt * z[1];
    result.state.theta = state.theta + settings.dt * z[2];
    result.state.vx = z[0];
    result.state.vy = z[1];
    result.state.omega = z[2];

    const auto final_residual = ComputeSdfNcpRigidBody2DResidual(contact_set, state, settings, z);
    result.contact_set = final_residual.contact_set;
    result.diagnostics =
        MakeSdfNcpRigidBody2DDiagnostics(contact_set, state, settings, z, success, iterations, message);
    return result;
}

}  // namespace sdfncp
}  // namespace chrono

#endif

