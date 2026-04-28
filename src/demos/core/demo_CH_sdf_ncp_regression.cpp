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
#include <iostream>
#include <string>

#include "chrono/collision/ChSdfNcpContact.h"

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

}  // namespace

int main(int argc, char* argv[]) {
    const std::string mode = argc > 1 ? argv[1] : "all";

    if (mode == "fb" || mode == "all") {
        TestFischerBurmeisterUtilities();
    }
    if (mode == "step" || mode == "all") {
        TestPointMassPlaneStep();
    }
    if (mode == "rollout" || mode == "all") {
        TestPointMassPlaneRollout();
    }

    if (mode != "fb" && mode != "step" && mode != "rollout" && mode != "all") {
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

