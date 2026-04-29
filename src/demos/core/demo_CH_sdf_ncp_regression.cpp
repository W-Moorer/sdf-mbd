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
// Standalone regression checks for the minimal SDF-NCP C++ prototype.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "chrono/collision/ChSdfNcpContact.h"

using chrono::ChVector3d;
using namespace chrono::sdfncp;

namespace {

int failures = 0;

void Check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << std::endl;
        failures++;
    }
}

void CheckNear(double value, double expected, double tolerance, const std::string& message) {
    if (std::abs(value - expected) > tolerance) {
        std::cerr << "FAILED: " << message << " value=" << value << " expected=" << expected
                  << " tolerance=" << tolerance << std::endl;
        failures++;
    }
}

void CheckVectorNear(const ChVector3d& value,
                     const ChVector3d& expected,
                     double tolerance,
                     const std::string& message) {
    const double error = (value - expected).Length();
    if (error > tolerance) {
        std::cerr << "FAILED: " << message << " error=" << error << " tolerance=" << tolerance << std::endl;
        failures++;
    }
}

std::filesystem::path GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "examples")) {
            return path;
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path();
}

void TestFischerBurmeisterUtilities() {
    const double eps = 1.0e-6;

    Check(std::abs(SmoothFischerBurmeister(0.2, 0.0, eps)) < 1.0e-9,
          "Open contact FB residual should be near zero for small eps");
    Check(std::abs(SmoothFischerBurmeister(0.0, 5.0, eps)) < 1.0e-9,
          "Closed contact FB residual should be near zero for small eps");

    const double g = 0.13;
    const double lambda = 0.27;
    const double h = 1.0e-6;
    const auto grad = SmoothFischerBurmeisterGrad(g, lambda, 1.0e-4);
    const double fd_g = (SmoothFischerBurmeister(g + h, lambda, 1.0e-4) -
                         SmoothFischerBurmeister(g - h, lambda, 1.0e-4)) /
                        (2.0 * h);
    const double fd_lambda = (SmoothFischerBurmeister(g, lambda + h, 1.0e-4) -
                              SmoothFischerBurmeister(g, lambda - h, 1.0e-4)) /
                             (2.0 * h);

    CheckNear(grad.dPhi_dg, fd_g, 1.0e-7, "FB dPhi/dg finite difference");
    CheckNear(grad.dPhi_dlambda, fd_lambda, 1.0e-7, "FB dPhi/dlambda finite difference");
    Check(ComplementarityError(0.2, -0.1) > 0.0, "Negative lambda should produce complementarity error");
    Check(ComplementarityError(-0.1, 0.2) > 0.0, "Negative gap should produce complementarity error");
}

void TestSdfSurfacePlane() {
    ChPlaneSdfSurface plane(ChVector3d(0, 2, 0), 0.0);

    CheckNear(plane.Phi(ChVector3d(0, 1, 0)), 1.0, 1.0e-12, "Plane SDF phi above y=0");
    CheckNear(plane.Phi(ChVector3d(0, -0.5, 0)), -0.5, 1.0e-12, "Plane SDF phi below y=0");
    CheckVectorNear(plane.Grad(ChVector3d(3, 4, 5)), ChVector3d(0, 1, 0), 1.0e-12,
                    "Plane SDF gradient should be normalized normal");
}

void TestSdfNcpContactConstraintEval() {
    ChPlaneSdfSurface plane(ChVector3d(0, 1, 0), 0.0);

    const ChSdfNcpContactConstraint open = EvaluateSdfNcpContact(plane, ChVector3d(0, 0.1, 0), 0.0, 1.0e-6);
    Check(open.gap > 0.0, "Open SDF-NCP contact should have positive gap");
    CheckNear(open.lambda_n, 0.0, 0.0, "Open SDF-NCP contact lambda");
    CheckNear(open.complementarity_error, 0.0, 1.0e-12, "Open SDF-NCP complementarity error");
    CheckVectorNear(open.jacobian_position, ChVector3d(0, 1, 0), 1.0e-12,
                    "Point-contact Jacobian should equal plane normal");

    const ChSdfNcpContactConstraint penetrating =
        EvaluateSdfNcpContact(plane, ChVector3d(0, -0.1, 0), 0.0, 1.0e-6);
    Check(penetrating.penetration > 0.0, "Penetrating SDF-NCP contact should report penetration");
    Check(penetrating.complementarity_error > 0.0, "Penetrating contact should have complementarity error");
}

void TestSdfNcpSharedGeometryQueryEval() {
    chrono::sdfcontact::ChSdfContactSampleQuery query;
    query.phi = -0.02;
    query.grad = ChVector3d(0, 3, 0);
    query.world_pos = ChVector3d(0.1, -0.02, 0.0);
    query.world_vel = ChVector3d(0, -1, 0);

    const ChSdfNcpContactConstraint contact = EvaluateSdfNcpContactQuery(query, 4.0, 1.0e-6);
    CheckNear(contact.gap, query.phi, 1.0e-12, "Shared SDF geometry query should provide NCP gap");
    CheckVectorNear(contact.normal, ChVector3d(0, 1, 0), 1.0e-12,
                    "Shared SDF geometry query gradient should provide NCP normal");
    CheckNear(contact.lambda_n, 4.0, 1.0e-12, "Shared SDF geometry query should preserve lambda");
    Check(contact.complementarity_error > 0.0, "Penetrating shared query should report complementarity error");
}

void TestPointMassPlaneStep() {
    PointMassPlaneSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 40;

    PointMassPlaneState state;
    state.x = 0.0;
    state.y = 1.0e-4;
    state.vx = 0.0;
    state.vy = -1.0;

    const PointMassPlaneStepResult step = SolvePointMassPlaneStep(state, settings);

    Check(step.diagnostics.success, "Point-mass SDF-NCP step should converge");
    Check(step.diagnostics.residual_norm < 1.0e-8, "Solved residual norm should be small");
    Check(step.diagnostics.gap >= -1.0e-10, "Solved step should not significantly penetrate the plane");
    Check(step.diagnostics.lambda_n >= -1.0e-10, "Solved lambda should not be significantly negative");
    Check(step.diagnostics.ncp_residual < 1.0e-8, "Solved NCP residual should be small");
    Check(step.diagnostics.complementarity_error < 1.0e-8, "Solved complementarity error should be small");
}

void TestPointMassPlaneRollout() {
    PointMassPlaneSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 40;

    PointMassPlaneState state;
    state.x = 0.0;
    state.y = 1.0;
    state.vx = 0.0;
    state.vy = 0.0;

    const int steps = 1000;
    double max_penetration = 0.0;
    double min_lambda = 0.0;
    double max_complementarity_error = 0.0;
    double max_residual_norm = 0.0;
    bool all_success = true;
    bool all_finite = true;

    for (int i = 0; i < steps; i++) {
        const PointMassPlaneStepResult step = SolvePointMassPlaneStep(state, settings);
        state = step.state;

        all_success = all_success && step.diagnostics.success;
        all_finite = all_finite && std::isfinite(state.x) && std::isfinite(state.y) &&
                     std::isfinite(state.vx) && std::isfinite(state.vy);
        max_penetration = std::max(max_penetration, step.diagnostics.penetration);
        min_lambda = std::min(min_lambda, step.diagnostics.lambda_n);
        max_complementarity_error =
            std::max(max_complementarity_error, step.diagnostics.complementarity_error);
        max_residual_norm = std::max(max_residual_norm, step.diagnostics.residual_norm);
    }

    Check(all_success, "All SDF-NCP rollout steps should converge");
    Check(all_finite, "SDF-NCP rollout should not produce NaN or Inf");
    Check(max_penetration < 1.0e-8, "SDF-NCP rollout penetration should stay small");
    Check(min_lambda >= -1.0e-8, "SDF-NCP rollout lambda should not be significantly negative");
    Check(max_complementarity_error < 1.0e-8, "SDF-NCP rollout complementarity error should stay small");
    Check(max_residual_norm < 1.0e-8, "SDF-NCP rollout residual norm should stay small");
}

void TestSdfNcpPointMass3DPlaneStep() {
    ChPlaneSdfSurface plane(ChVector3d(0, 1, 0), 0.0);

    ChSdfNcpPointMassSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 40;

    ChSdfNcpPointMassState state;
    state.position = ChVector3d(0.0, 0.1, 0.0);
    state.velocity = ChVector3d(0.0, -1.0, 0.0);

    const ChSdfNcpPointMassStepResult step = SolveSdfNcpPointMassStep(plane, state, settings);

    Check(step.diagnostics.success, "Generic 3D point-mass SDF-NCP step should converge");
    Check(std::isfinite(step.diagnostics.residual_norm), "Generic 3D residual norm should be finite");
    Check(std::isfinite(step.state.position.x()) && std::isfinite(step.state.position.y()) &&
              std::isfinite(step.state.position.z()),
          "Generic 3D point-mass step should not produce NaN position");
    Check(step.diagnostics.lambda_n >= -1.0e-10, "Generic 3D lambda should not be significantly negative");
}

void TestSdfNcpPointMass3DPlaneRollout() {
    ChPlaneSdfSurface plane(ChVector3d(0, 1, 0), 0.0);

    ChSdfNcpPointMassSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 40;

    ChSdfNcpPointMassState state;
    state.position = ChVector3d(0.0, 1.0, 0.0);
    state.velocity = ChVector3d(0.0, 0.0, 0.0);

    const int steps = 1000;
    double max_penetration = 0.0;
    double min_lambda = 0.0;
    double max_complementarity_error = 0.0;
    double max_residual_norm = 0.0;
    bool all_success = true;
    bool all_finite = true;

    for (int i = 0; i < steps; i++) {
        const ChSdfNcpPointMassStepResult step = SolveSdfNcpPointMassStep(plane, state, settings);
        state = step.state;

        all_success = all_success && step.diagnostics.success;
        all_finite = all_finite && std::isfinite(state.position.x()) && std::isfinite(state.position.y()) &&
                     std::isfinite(state.position.z()) && std::isfinite(state.velocity.x()) &&
                     std::isfinite(state.velocity.y()) && std::isfinite(state.velocity.z());
        max_penetration = std::max(max_penetration, step.diagnostics.penetration);
        min_lambda = std::min(min_lambda, step.diagnostics.lambda_n);
        max_complementarity_error =
            std::max(max_complementarity_error, step.diagnostics.complementarity_error);
        max_residual_norm = std::max(max_residual_norm, step.diagnostics.residual_norm);
    }

    Check(all_success, "Generic 3D SDF-NCP rollout steps should converge");
    Check(all_finite, "Generic 3D SDF-NCP rollout should not produce NaN or Inf");
    Check(max_penetration < 1.0e-8, "Generic 3D SDF-NCP rollout penetration should stay small");
    Check(min_lambda >= -1.0e-8, "Generic 3D SDF-NCP rollout lambda should not be significantly negative");
    Check(max_complementarity_error < 1.0e-8,
          "Generic 3D SDF-NCP rollout complementarity error should stay small");
    Check(max_residual_norm < 1.0e-8, "Generic 3D SDF-NCP rollout residual norm should stay small");
}

void TestSdfNcpContactSetEval() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);
    auto wall = std::make_shared<ChPlaneSdfSurface>(ChVector3d(1, 0, 0), 0.0);

    ChSdfNcpContactSet contact_set;
    contact_set.AddContact(ground, 0.0);
    contact_set.AddContact(wall, 0.0);
    contact_set.Update(ChVector3d(0.1, 0.2, 0.0), 1.0e-6);

    Check(contact_set.GetNumContacts() == 2, "Contact set should contain two contacts");
    CheckNear(contact_set.GetContact(0).gap, 0.2, 1.0e-12, "Ground contact-set gap");
    CheckNear(contact_set.GetContact(1).gap, 0.1, 1.0e-12, "Wall contact-set gap");
    CheckVectorNear(contact_set.GetContact(0).normal, ChVector3d(0, 1, 0), 1.0e-12,
                    "Ground contact-set normal");
    CheckVectorNear(contact_set.GetContact(1).normal, ChVector3d(1, 0, 0), 1.0e-12,
                    "Wall contact-set normal");
    CheckVectorNear(contact_set.ComputeTotalContactForce(), ChVector3d(0, 0, 0), 1.0e-12,
                    "Zero multipliers should produce zero total contact force");
}

void TestSdfNcpContactSetForce() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);
    auto wall = std::make_shared<ChPlaneSdfSurface>(ChVector3d(1, 0, 0), 0.0);

    ChSdfNcpContactSet contact_set;
    contact_set.AddContact(ground, 2.0);
    contact_set.AddContact(wall, 3.0);
    contact_set.Update(ChVector3d(0.1, 0.2, 0.0), 1.0e-6);

    CheckVectorNear(contact_set.ComputeTotalContactForce(), ChVector3d(3, 2, 0), 1.0e-12,
                    "Contact set should sum normal forces");
}

void TestSdfNcpPointMassContactSetN1MatchesSingle() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);

    ChSdfNcpPointMassSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 40;

    ChSdfNcpPointMassState state;
    state.position = ChVector3d(0.0, 1.0e-4, 0.0);
    state.velocity = ChVector3d(0.0, -1.0, 0.0);

    ChSdfNcpContactSet contact_set;
    contact_set.AddContact(ground, 0.0);

    const ChSdfNcpPointMassStepResult single = SolveSdfNcpPointMassStep(*ground, state, settings);
    const ChSdfNcpPointMassContactSetStepResult assembled =
        SolveSdfNcpPointMassContactSetStep(contact_set, state, settings);

    Check(single.diagnostics.success, "Single-contact generic residual should converge");
    Check(assembled.diagnostics.success, "N=1 contact-set residual should converge");
    CheckVectorNear(assembled.state.position, single.state.position, 1.0e-9,
                    "N=1 contact-set position should match single residual");
    CheckVectorNear(assembled.state.velocity, single.state.velocity, 1.0e-9,
                    "N=1 contact-set velocity should match single residual");
    CheckNear(assembled.contact_set.GetContact(0).lambda_n, single.diagnostics.lambda_n, 1.0e-8,
              "N=1 contact-set lambda should match single residual");
    CheckNear(assembled.contact_set.GetContact(0).gap, single.diagnostics.gap, 1.0e-9,
              "N=1 contact-set gap should match single residual");
}

void TestSdfNcpPointMassContactSetRolloutTwoPlanes() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);
    auto wall = std::make_shared<ChPlaneSdfSurface>(ChVector3d(1, 0, 0), 0.0);

    ChSdfNcpPointMassSettings settings;
    settings.mass = 1.0;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 50;

    ChSdfNcpPointMassState state;
    state.position = ChVector3d(0.5, 1.0, 0.0);
    state.velocity = ChVector3d(-1.0, 0.0, 0.0);

    ChSdfNcpContactSet contact_set;
    contact_set.AddContact(ground, 0.0);
    contact_set.AddContact(wall, 0.0);

    const int steps = 1000;
    double max_penetration = 0.0;
    double min_lambda = 0.0;
    double max_complementarity_error = 0.0;
    double max_residual_norm = 0.0;
    bool all_success = true;
    bool all_finite = true;

    for (int i = 0; i < steps; i++) {
        const ChSdfNcpPointMassContactSetStepResult step =
            SolveSdfNcpPointMassContactSetStep(contact_set, state, settings);
        state = step.state;
        contact_set = step.contact_set;

        all_success = all_success && step.diagnostics.success;
        all_finite = all_finite && std::isfinite(state.position.x()) && std::isfinite(state.position.y()) &&
                     std::isfinite(state.position.z()) && std::isfinite(state.velocity.x()) &&
                     std::isfinite(state.velocity.y()) && std::isfinite(state.velocity.z());
        max_penetration = std::max(max_penetration, step.diagnostics.max_penetration);
        min_lambda = std::min(min_lambda, step.diagnostics.min_lambda_n);
        max_complementarity_error =
            std::max(max_complementarity_error, step.diagnostics.max_complementarity_error);
        max_residual_norm = std::max(max_residual_norm, step.diagnostics.residual_norm);
    }

    Check(all_success, "Two-plane contact-set rollout steps should converge");
    Check(all_finite, "Two-plane contact-set rollout should not produce NaN or Inf");
    Check(max_penetration < 1.0e-8, "Two-plane contact-set rollout penetration should stay small");
    Check(min_lambda >= -1.0e-8, "Two-plane contact-set rollout lambda should not be significantly negative");
    Check(max_complementarity_error < 1.0e-8,
          "Two-plane contact-set rollout complementarity error should stay small");
    Check(max_residual_norm < 1.0e-8, "Two-plane contact-set rollout residual norm should stay small");
}

void TestSdfNcpRigidBody2DKinematics() {
    const ChVector3d r_local(0.5, -0.25, 0.0);

    const ChVector3d world0 = RigidBody2DLocalPointToWorld(ChVector3d(1.0, 2.0, 0.0), r_local);
    CheckVectorNear(world0, ChVector3d(1.5, 1.75, 0.0), 1.0e-12,
                    "Rigid-body 2D local point should translate without rotation");

    const double pi = std::acos(-1.0);
    const ChVector3d world90 = RigidBody2DLocalPointToWorld(ChVector3d(1.0, 2.0, 0.5 * pi), r_local);
    CheckVectorNear(world90, ChVector3d(1.25, 2.5, 0.0), 1.0e-12,
                    "Rigid-body 2D local point should rotate by pi/2");
}

void TestSdfNcpRigidBody2DJacobianFiniteDifference() {
    ChPlaneSdfSurface ground(ChVector3d(0, 1, 0), 0.0);
    const ChVector3d q(0.2, 0.8, 0.4);
    const ChVector3d r_local(0.5, -0.25, 0.0);

    const auto contact = EvaluateRigidBody2DSdfContact(ground, q, r_local, 0.0, 1.0e-6);

    const auto shift_q = [](const ChVector3d& value, int index, double delta) {
        return ChVector3d(value.x() + (index == 0 ? delta : 0.0),
                          value.y() + (index == 1 ? delta : 0.0),
                          value.z() + (index == 2 ? delta : 0.0));
    };
    const auto gap_at = [&ground, &r_local](const ChVector3d& value) {
        return ground.Phi(RigidBody2DLocalPointToWorld(value, r_local));
    };

    const double h = 1.0e-6;
    ChVector3d fd(0, 0, 0);
    fd.x() = (gap_at(shift_q(q, 0, h)) - gap_at(shift_q(q, 0, -h))) / (2.0 * h);
    fd.y() = (gap_at(shift_q(q, 1, h)) - gap_at(shift_q(q, 1, -h))) / (2.0 * h);
    fd.z() = (gap_at(shift_q(q, 2, h)) - gap_at(shift_q(q, 2, -h))) / (2.0 * h);

    CheckVectorNear(contact.jacobian, fd, 1.0e-8,
                    "Rigid-body 2D SDF contact Jacobian should match finite differences");
}

void TestSdfNcpRigidBody2DContactForce() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);
    const ChVector3d q(0.0, 0.6, 0.2);
    const ChVector3d r_local(0.0, -0.5, 0.0);
    const double lambda = 10.0;

    ChSdfNcpRigidBody2DContactSet contact_set;
    contact_set.AddContact(ground, r_local, lambda);
    contact_set.Update(q, 1.0e-6);

    const ChVector3d generalized_force = contact_set.ComputeGeneralizedContactForce();
    const ChVector3d expected = contact_set.GetContact(0).jacobian * lambda;

    CheckNear(generalized_force.x(), expected.x(), 1.0e-12, "Rigid-body 2D contact generalized Fx");
    Check(generalized_force.y() > 0.0, "Rigid-body 2D contact generalized Fy should be positive");
    CheckNear(generalized_force.y(), expected.y(), 1.0e-12, "Rigid-body 2D contact generalized Fy value");
    CheckNear(generalized_force.z(), expected.z(), 1.0e-12, "Rigid-body 2D contact torque");
}

void TestSdfNcpRigidBody2DStepSinglePoint() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);

    ChSdfRigidBody2DSettings settings;
    settings.mass = 1.0;
    settings.inertia_z = 0.1;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 60;

    ChSdfRigidBody2DState state;
    state.x = 0.0;
    state.y = 0.251;
    state.theta = 0.0;
    state.vx = 0.0;
    state.vy = -1.0;
    state.omega = 0.0;

    ChSdfNcpRigidBody2DContactSet contact_set;
    contact_set.AddContact(ground, ChVector3d(0.0, -0.25, 0.0), 0.0);

    const ChSdfNcpRigidBody2DStepResult step = SolveSdfNcpRigidBody2DStep(contact_set, state, settings);

    Check(step.diagnostics.success, "Rigid-body 2D single-point step should converge");
    Check(std::isfinite(step.diagnostics.residual_norm), "Rigid-body 2D residual norm should be finite");
    Check(std::isfinite(step.state.x) && std::isfinite(step.state.y) && std::isfinite(step.state.theta),
          "Rigid-body 2D single-point step should not produce NaN position");
    Check(std::isfinite(step.state.vx) && std::isfinite(step.state.vy) && std::isfinite(step.state.omega),
          "Rigid-body 2D single-point step should not produce NaN velocity");
    Check(step.diagnostics.min_lambda_n >= -1.0e-8,
          "Rigid-body 2D single-point lambda should not be significantly negative");
}

void TestSdfNcpRigidBody2DRolloutTwoBottomPoints() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);

    ChSdfRigidBody2DSettings settings;
    settings.mass = 1.0;
    settings.inertia_z = 0.10416666666666667;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 80;

    ChSdfRigidBody2DState state;
    state.x = 0.0;
    state.y = 1.0;
    state.theta = 0.1;
    state.vx = 0.0;
    state.vy = 0.0;
    state.omega = 0.0;

    ChSdfNcpRigidBody2DContactSet contact_set;
    contact_set.AddContact(ground, ChVector3d(-0.5, -0.25, 0.0), 0.0);
    contact_set.AddContact(ground, ChVector3d(0.5, -0.25, 0.0), 0.0);

    const int steps = 1000;
    double max_penetration = 0.0;
    double min_lambda = 0.0;
    double max_complementarity_error = 0.0;
    double max_residual_norm = 0.0;
    bool all_success = true;
    bool all_finite = true;

    for (int i = 0; i < steps; i++) {
        const ChSdfNcpRigidBody2DStepResult step = SolveSdfNcpRigidBody2DStep(contact_set, state, settings);
        state = step.state;
        contact_set = step.contact_set;

        all_success = all_success && step.diagnostics.success;
        all_finite = all_finite && std::isfinite(state.x) && std::isfinite(state.y) &&
                     std::isfinite(state.theta) && std::isfinite(state.vx) && std::isfinite(state.vy) &&
                     std::isfinite(state.omega);
        max_penetration = std::max(max_penetration, step.diagnostics.max_penetration);
        min_lambda = std::min(min_lambda, step.diagnostics.min_lambda_n);
        max_complementarity_error =
            std::max(max_complementarity_error, step.diagnostics.max_complementarity_error);
        max_residual_norm = std::max(max_residual_norm, step.diagnostics.residual_norm);
    }

    Check(all_success, "Rigid-body 2D two-point rollout steps should converge");
    Check(all_finite, "Rigid-body 2D rollout should not produce NaN or Inf");
    Check(max_penetration < 1.0e-7, "Rigid-body 2D rollout penetration should stay bounded");
    Check(min_lambda >= -1.0e-8, "Rigid-body 2D rollout lambda should not be significantly negative");
    Check(max_complementarity_error < 1.0e-7,
          "Rigid-body 2D rollout complementarity error should stay bounded");
    Check(max_residual_norm < 1.0e-7, "Rigid-body 2D rollout residual norm should stay bounded");
}

int ExportSdfNcpRigidBody2DRollout() {
    auto ground = std::make_shared<ChPlaneSdfSurface>(ChVector3d(0, 1, 0), 0.0);

    ChSdfRigidBody2DSettings settings;
    settings.mass = 1.0;
    settings.inertia_z = 0.10416666666666667;
    settings.gravity = 9.81;
    settings.dt = 1.0e-3;
    settings.eps = 1.0e-6;
    settings.tolerance = 1.0e-11;
    settings.max_iterations = 80;

    ChSdfRigidBody2DState state;
    state.x = 0.0;
    state.y = 1.0;
    state.theta = 0.1;
    state.vx = 0.0;
    state.vy = 0.0;
    state.omega = 0.0;

    ChSdfNcpRigidBody2DContactSet contact_set;
    contact_set.AddContact(ground, ChVector3d(-0.5, -0.25, 0.0), 0.0);
    contact_set.AddContact(ground, ChVector3d(0.5, -0.25, 0.0), 0.0);

    const auto out_dir = GetProjectRoot() / "results" / "sdf_ncp_cpp";
    std::filesystem::create_directories(out_dir);
    const auto csv_path = out_dir / "rigidbody2d_rollout.csv";

    std::ofstream out(csv_path);
    if (!out) {
        std::cerr << "Failed to open " << csv_path.string() << " for writing." << std::endl;
        return 1;
    }

    out << std::setprecision(17);
    out << "time,x,y,theta,vx,vy,omega,gap0,gap1,lambda0,lambda1,max_penetration,"
           "max_complementarity_error,ncp_residual_norm,residual_norm,iterations,success\n";

    auto write_row = [&out](double time,
                            const ChSdfRigidBody2DState& row_state,
                            const ChSdfNcpRigidBody2DContactSet& row_contacts,
                            const ChSdfNcpRigidBody2DDiagnostics& diag) {
        const auto& c0 = row_contacts.GetContact(0);
        const auto& c1 = row_contacts.GetContact(1);
        out << time << "," << row_state.x << "," << row_state.y << "," << row_state.theta << ","
            << row_state.vx << "," << row_state.vy << "," << row_state.omega << ","
            << c0.gap << "," << c1.gap << "," << c0.lambda_n << "," << c1.lambda_n << ","
            << diag.max_penetration << "," << diag.max_complementarity_error << ","
            << diag.ncp_residual_norm << "," << diag.residual_norm << "," << diag.iterations << ","
            << (diag.success ? 1 : 0) << "\n";
    };

    contact_set.Update(ChVector3d(state.x, state.y, state.theta), settings.eps);
    ChSdfNcpRigidBody2DDiagnostics initial_diag;
    initial_diag.success = true;
    initial_diag.iterations = 0;
    initial_diag.residual_norm = 0.0;
    initial_diag.max_penetration = contact_set.MaxPenetration();
    initial_diag.max_complementarity_error = contact_set.MaxComplementarityError();
    initial_diag.ncp_residual_norm = contact_set.NcpResidualNorm();
    write_row(0.0, state, contact_set, initial_diag);

    const int steps = 1000;
    bool all_success = true;
    for (int i = 1; i <= steps; i++) {
        const ChSdfNcpRigidBody2DStepResult step = SolveSdfNcpRigidBody2DStep(contact_set, state, settings);
        state = step.state;
        contact_set = step.contact_set;
        all_success = all_success && step.diagnostics.success;
        write_row(static_cast<double>(i) * settings.dt, state, contact_set, step.diagnostics);
    }

    std::cout << "Wrote " << csv_path.string() << std::endl;
    if (!all_success) {
        std::cerr << "Warning: at least one rigidbody2d_export step did not satisfy the solver success criterion."
                  << std::endl;
    }
    return all_success ? 0 : 1;
}

}  // namespace

int main(int argc, char* argv[]) {
    const std::string mode = argc > 1 ? argv[1] : "all";

    if (mode == "rigidbody2d_export") {
        return ExportSdfNcpRigidBody2DRollout();
    }

    if (mode == "fb" || mode == "all") {
        TestFischerBurmeisterUtilities();
    }
    if (mode == "surface" || mode == "all") {
        TestSdfSurfacePlane();
    }
    if (mode == "constraint" || mode == "all") {
        TestSdfNcpContactConstraintEval();
    }
    if (mode == "shared_geometry" || mode == "all") {
        TestSdfNcpSharedGeometryQueryEval();
    }
    if (mode == "step" || mode == "all") {
        TestPointMassPlaneStep();
    }
    if (mode == "rollout" || mode == "all") {
        TestPointMassPlaneRollout();
    }
    if (mode == "step3d" || mode == "all") {
        TestSdfNcpPointMass3DPlaneStep();
    }
    if (mode == "rollout3d" || mode == "all") {
        TestSdfNcpPointMass3DPlaneRollout();
    }
    if (mode == "contact_set_eval" || mode == "all") {
        TestSdfNcpContactSetEval();
    }
    if (mode == "contact_set_force" || mode == "all") {
        TestSdfNcpContactSetForce();
    }
    if (mode == "contact_set_n1" || mode == "all") {
        TestSdfNcpPointMassContactSetN1MatchesSingle();
    }
    if (mode == "contact_set_rollout" || mode == "all") {
        TestSdfNcpPointMassContactSetRolloutTwoPlanes();
    }
    if (mode == "rigidbody2d_kinematics" || mode == "all") {
        TestSdfNcpRigidBody2DKinematics();
    }
    if (mode == "rigidbody2d_jacobian" || mode == "all") {
        TestSdfNcpRigidBody2DJacobianFiniteDifference();
    }
    if (mode == "rigidbody2d_force" || mode == "all") {
        TestSdfNcpRigidBody2DContactForce();
    }
    if (mode == "rigidbody2d_step" || mode == "all") {
        TestSdfNcpRigidBody2DStepSinglePoint();
    }
    if (mode == "rigidbody2d_rollout" || mode == "all") {
        TestSdfNcpRigidBody2DRolloutTwoBottomPoints();
    }

    if (mode != "fb" && mode != "surface" && mode != "constraint" && mode != "shared_geometry" &&
        mode != "step" && mode != "rollout" && mode != "step3d" && mode != "rollout3d" &&
        mode != "contact_set_eval" && mode != "contact_set_force" && mode != "contact_set_n1" &&
        mode != "contact_set_rollout" && mode != "rigidbody2d_kinematics" &&
        mode != "rigidbody2d_jacobian" && mode != "rigidbody2d_force" &&
        mode != "rigidbody2d_step" && mode != "rigidbody2d_rollout" &&
        mode != "rigidbody2d_export" && mode != "all") {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 2;
    }

    if (failures > 0) {
        std::cerr << failures << " SDF-NCP regression checks failed." << std::endl;
        return 1;
    }

    std::cout << "SDF-NCP regression checks passed for mode '" << mode << "'." << std::endl;
    return 0;
}

