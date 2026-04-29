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
// Standalone SDF-NCP benchmarks reusing existing field_contact asset parameters.
//
// This executable intentionally remains independent from Chrono's production
// contact containers and solver descriptor.  It provides a small reproducible
// benchmark layer for frictionless SDF-NCP contact.
//
// =============================================================================

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "chrono/collision/ChSdfNcpContact.h"
#include "chrono/core/ChVector3.h"

using chrono::ChVector3d;
using namespace chrono::sdfncp;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct BenchmarkStats {
    std::string case_name;
    std::string asset;
    double dt = 0.0;
    double t_end = 0.0;
    int num_steps = 0;
    int num_contacts = 0;
    double max_penetration = 0.0;
    double sum_penetration = 0.0;
    double max_lambda_n = 0.0;
    double max_ncp_residual = 0.0;
    double sum_ncp_residual = 0.0;
    double max_complementarity_error = 0.0;
    double sum_complementarity_error = 0.0;
    double sum_iterations = 0.0;
    double sum_success = 0.0;
    int samples = 0;
    double runtime_seconds = 0.0;
    std::string notes;
};

struct SphereBodyState {
    ChVector3d pos = ChVector3d(0, 0, 0);
    ChVector3d vel = ChVector3d(0, 0, 0);
    ChVector3d omega = ChVector3d(0, 0, 0);
    std::array<double, 4> quat = {{1.0, 0.0, 0.0, 0.0}};
    double mass = 1.0;
    double radius = 0.05;
};

struct SpherePairConfig {
    std::string case_name;
    std::string asset;
    double radius = 0.05;
    double density_a = 1000.0;
    double density_b = 1000.0;
    ChVector3d pos_a0 = ChVector3d(-0.15, 0, 0);
    ChVector3d pos_b0 = ChVector3d(0.15, 0, 0);
    ChVector3d vel_a0 = ChVector3d(1, 0, 0);
    ChVector3d vel_b0 = ChVector3d(0, 0, 0);
    double dt = 5.0e-4;
    double t_end = 0.5;
    double eps = 1.0e-7;
    double tolerance = 1.0e-10;
    int max_iterations = 40;
};

struct SpherePairStepDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double gap = 0.0;
    double lambda_n = 0.0;
    double penetration = 0.0;
    double ncp_residual = 0.0;
    double complementarity_error = 0.0;
};

struct SpherePairStepResult {
    SphereBodyState a;
    SphereBodyState b;
    SpherePairStepDiagnostics diagnostics;
};

std::filesystem::path GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "assets")) {
            return path;
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path();
}

double SphereMass(double radius, double density) {
    return density * 4.0 * kPi * radius * radius * radius / 3.0;
}

double ElasticPostV1(double v1, double v2, double m1, double m2) {
    return ((m1 - m2) / (m1 + m2)) * v1 + (2.0 * m2 / (m1 + m2)) * v2;
}

double ElasticPostV2(double v1, double v2, double m1, double m2) {
    return (2.0 * m1 / (m1 + m2)) * v1 + ((m2 - m1) / (m1 + m2)) * v2;
}

double InelasticCommonVelocity(double v1, double v2, double m1, double m2) {
    return (m1 * v1 + m2 * v2) / (m1 + m2);
}

ChVector3d SafeNormalFromAToB(const ChVector3d& delta) {
    const double length = delta.Length();
    if (length <= 1.0e-14) {
        return ChVector3d(1, 0, 0);
    }
    return delta / length;
}

double SphereGap(const SphereBodyState& a, const SphereBodyState& b) {
    return (b.pos - a.pos).Length() - (a.radius + b.radius);
}

std::vector<double> ComputeSpherePairResidual(const SphereBodyState& a,
                                              const SphereBodyState& b,
                                              double dt,
                                              double eps,
                                              const std::vector<double>& z,
                                              SpherePairStepDiagnostics* diag = nullptr) {
    if (z.size() != 7) {
        throw std::invalid_argument("Sphere-pair SDF-NCP residual expects 7 unknowns.");
    }

    const ChVector3d va_next(z[0], z[1], z[2]);
    const ChVector3d vb_next(z[3], z[4], z[5]);
    const double lambda = z[6];
    const ChVector3d pa_next = a.pos + va_next * dt;
    const ChVector3d pb_next = b.pos + vb_next * dt;
    const ChVector3d normal = SafeNormalFromAToB(pb_next - pa_next);
    const double gap = (pb_next - pa_next).Length() - (a.radius + b.radius);
    const double ncp = SmoothFischerBurmeister(gap, lambda, eps);

    const ChVector3d force_a = normal * (-lambda);
    const ChVector3d force_b = normal * lambda;
    const ChVector3d rv_a = a.mass * (va_next - a.vel) - dt * force_a;
    const ChVector3d rv_b = b.mass * (vb_next - b.vel) - dt * force_b;

    if (diag) {
        diag->gap = gap;
        diag->lambda_n = lambda;
        diag->penetration = std::max(0.0, -gap);
        diag->ncp_residual = std::abs(ncp);
        diag->complementarity_error = ComplementarityError(gap, lambda);
    }

    return {rv_a.x(), rv_a.y(), rv_a.z(), rv_b.x(), rv_b.y(), rv_b.z(), ncp};
}

std::vector<std::vector<double>> FiniteDifferenceJacobian(
    const std::function<std::vector<double>(const std::vector<double>&)>& residual,
    const std::vector<double>& z,
    double h = 1.0e-7) {
    const size_t n = z.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));
    for (size_t col = 0; col < n; col++) {
        std::vector<double> zp = z;
        std::vector<double> zm = z;
        zp[col] += h;
        zm[col] -= h;
        const auto rp = residual(zp);
        const auto rm = residual(zm);
        for (size_t row = 0; row < n; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }
    return jac;
}

SpherePairStepResult SolveSpherePairStep(const SphereBodyState& a,
                                         const SphereBodyState& b,
                                         const SpherePairConfig& config) {
    SpherePairStepResult result;
    result.a = a;
    result.b = b;

    std::vector<double> z = {a.vel.x(), a.vel.y(), a.vel.z(), b.vel.x(), b.vel.y(), b.vel.z(), 0.0};
    const ChVector3d pa_guess = a.pos + a.vel * config.dt;
    const ChVector3d pb_guess = b.pos + b.vel * config.dt;
    const double gap_guess = (pb_guess - pa_guess).Length() - 2.0 * config.radius;
    const double reduced_mass = (a.mass * b.mass) / (a.mass + b.mass);
    z[6] = std::max(0.0, -gap_guess) * reduced_mass / std::max(config.dt * config.dt, 1.0e-12);

    auto residual = [&](const std::vector<double>& zz) {
        return ComputeSpherePairResidual(a, b, config.dt, config.eps, zz, nullptr);
    };

    double residual_norm = NormDynamic(residual(z));
    bool success = residual_norm <= config.tolerance;
    int iterations = 0;

    for (iterations = 0; !success && iterations < config.max_iterations; iterations++) {
        const auto r = residual(z);
        const auto jac = FiniteDifferenceJacobian(residual, z);
        std::vector<double> rhs(r.size(), 0.0);
        for (size_t i = 0; i < r.size(); i++) {
            rhs[i] = -r[i];
        }
        std::vector<double> dz(z.size(), 0.0);
        if (!SolveLinearDynamic(jac, rhs, dz)) {
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
            candidate_norm = NormDynamic(residual(candidate));
            if (candidate_norm <= residual_norm || alpha <= 1.0e-4) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            break;
        }
        z = candidate;
        residual_norm = candidate_norm;
        if (residual_norm <= config.tolerance) {
            success = true;
            iterations++;
            break;
        }
    }

    if (!success && residual_norm <= std::max(config.tolerance, 1.0e-9)) {
        success = true;
    }
    if (z[6] < -1.0e-8) {
        success = false;
    } else if (z[6] < 0.0) {
        z[6] = 0.0;
    }

    result.a.vel = ChVector3d(z[0], z[1], z[2]);
    result.b.vel = ChVector3d(z[3], z[4], z[5]);
    result.a.pos = a.pos + result.a.vel * config.dt;
    result.b.pos = b.pos + result.b.vel * config.dt;

    result.diagnostics.success = success;
    result.diagnostics.iterations = iterations;
    result.diagnostics.residual_norm = residual_norm;
    ComputeSpherePairResidual(a, b, config.dt, config.eps, z, &result.diagnostics);
    return result;
}

void WriteTrajectoryHeader(std::ofstream& out) {
    out << "time,body_id,px,py,pz,q0,q1,q2,q3,vx,vy,vz,wx,wy,wz,contact_id,gap,lambda_n,"
           "penetration,ncp_residual,complementarity_error,residual_norm,iterations,success\n";
}

void WriteSummary(const std::filesystem::path& path, const BenchmarkStats& stats) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "case_name,method,asset,dt,t_end,num_steps,num_contacts,max_penetration,mean_penetration,"
           "max_lambda_n,max_ncp_residual,mean_ncp_residual,max_complementarity_error,"
           "mean_complementarity_error,mean_iterations,success_rate,runtime_seconds,notes\n";
    const double denom = static_cast<double>(std::max(1, stats.samples));
    out << stats.case_name << ",sdf_ncp," << stats.asset << "," << stats.dt << "," << stats.t_end << ","
        << stats.num_steps << "," << stats.num_contacts << "," << stats.max_penetration << ","
        << stats.sum_penetration / denom << "," << stats.max_lambda_n << "," << stats.max_ncp_residual << ","
        << stats.sum_ncp_residual / denom << "," << stats.max_complementarity_error << ","
        << stats.sum_complementarity_error / denom << "," << stats.sum_iterations / denom << ","
        << stats.sum_success / denom << "," << stats.runtime_seconds << "," << stats.notes << "\n";
}

void Accumulate(BenchmarkStats& stats,
                double penetration,
                double lambda_n,
                double ncp_residual,
                double complementarity_error,
                double iterations,
                bool success) {
    stats.max_penetration = std::max(stats.max_penetration, penetration);
    stats.sum_penetration += penetration;
    stats.max_lambda_n = std::max(stats.max_lambda_n, std::abs(lambda_n));
    stats.max_ncp_residual = std::max(stats.max_ncp_residual, ncp_residual);
    stats.sum_ncp_residual += ncp_residual;
    stats.max_complementarity_error = std::max(stats.max_complementarity_error, complementarity_error);
    stats.sum_complementarity_error += complementarity_error;
    stats.sum_iterations += iterations;
    stats.sum_success += success ? 1.0 : 0.0;
    stats.samples++;
}

int RunSpherePairCase(const SpherePairConfig& config) {
    const auto start = std::chrono::steady_clock::now();
    const auto out_dir = GetProjectRoot() / "results" / "sdf_ncp_benchmarks" / config.case_name;
    std::filesystem::create_directories(out_dir);
    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";

    SphereBodyState a;
    SphereBodyState b;
    a.pos = config.pos_a0;
    b.pos = config.pos_b0;
    a.vel = config.vel_a0;
    b.vel = config.vel_b0;
    a.radius = config.radius;
    b.radius = config.radius;
    a.mass = SphereMass(config.radius, config.density_a);
    b.mass = SphereMass(config.radius, config.density_b);

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = config.case_name;
    stats.asset = config.asset;
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    stats.num_contacts = 1;
    stats.notes = "analytic sphere SDF using field_contact asset parameters; no restitution";

    auto write_body = [&](double time, const std::string& id, const SphereBodyState& body, const SpherePairStepDiagnostics& diag) {
        trajectory << time << "," << id << "," << body.pos.x() << "," << body.pos.y() << "," << body.pos.z()
                   << "," << body.quat[0] << "," << body.quat[1] << "," << body.quat[2] << "," << body.quat[3]
                   << "," << body.vel.x() << "," << body.vel.y() << "," << body.vel.z() << ","
                   << body.omega.x() << "," << body.omega.y() << "," << body.omega.z() << ",0,"
                   << diag.gap << "," << diag.lambda_n << "," << diag.penetration << ","
                   << diag.ncp_residual << "," << diag.complementarity_error << "," << diag.residual_norm << ","
                   << diag.iterations << "," << (diag.success ? 1 : 0) << "\n";
    };

    for (int i = 0; i <= steps; i++) {
        SpherePairStepDiagnostics diag;
        diag.success = true;
        diag.gap = SphereGap(a, b);
        diag.penetration = std::max(0.0, -diag.gap);
        diag.lambda_n = 0.0;
        diag.ncp_residual = std::abs(SmoothFischerBurmeister(diag.gap, 0.0, config.eps));
        diag.complementarity_error = ComplementarityError(diag.gap, 0.0);
        if (i > 0) {
            const SpherePairStepResult step = SolveSpherePairStep(a, b, config);
            a = step.a;
            b = step.b;
            diag = step.diagnostics;
        }

        const double time = static_cast<double>(i) * config.dt;
        write_body(time, "sphere_a", a, diag);
        write_body(time, "sphere_b", b, diag);
        Accumulate(stats,
                   diag.penetration,
                   diag.lambda_n,
                   diag.ncp_residual,
                   diag.complementarity_error,
                   static_cast<double>(diag.iterations),
                   diag.success);
    }

    const double elastic_a = ElasticPostV1(config.vel_a0.x(), config.vel_b0.x(), a.mass, b.mass);
    const double elastic_b = ElasticPostV2(config.vel_a0.x(), config.vel_b0.x(), a.mass, b.mass);
    const double inelastic_common =
        InelasticCommonVelocity(config.vel_a0.x(), config.vel_b0.x(), a.mass, b.mass);
    auto write_comparison = [&](const std::string& metric, double reference, double value, const std::string& notes) {
        const double abs_diff = std::abs(value - reference);
        const double rel_diff = std::abs(reference) > 1.0e-14 ? abs_diff / std::abs(reference) : 0.0;
        comparison << config.case_name << "," << metric << "," << reference << "," << value << "," << abs_diff << ","
                   << rel_diff << "," << notes << "\n";
    };
    write_comparison("postimpact_vx_sphere_a", elastic_a, a.vel.x(), "analytic elastic reference from field_contact paper example");
    write_comparison("postimpact_vx_sphere_b", elastic_b, b.vel.x(), "analytic elastic reference from field_contact paper example");
    write_comparison("postimpact_vx_sphere_a_inelastic_common",
                     inelastic_common,
                     a.vel.x(),
                     "analytic no-restitution common normal velocity expected by the frictionless SDF-NCP benchmark");
    write_comparison("postimpact_vx_sphere_b_inelastic_common",
                     inelastic_common,
                     b.vel.x(),
                     "analytic no-restitution common normal velocity expected by the frictionless SDF-NCP benchmark");

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);
    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    return stats.max_penetration < 1.0e-7 && stats.max_complementarity_error < 1.0e-6 &&
                   stats.sum_success / static_cast<double>(std::max(1, stats.samples)) > 0.99
               ? 0
               : 1;
}

SpherePairConfig HeadonSpheresConfig() {
    SpherePairConfig config;
    config.case_name = "headon_spheres";
    config.asset = "assets/headon_spheres";
    config.radius = 0.05;
    config.density_a = 1000.0;
    config.density_b = 1000.0;
    config.pos_a0 = ChVector3d(-0.15, 0, 0);
    config.pos_b0 = ChVector3d(0.15, 0, 0);
    config.vel_a0 = ChVector3d(1.0, 0, 0);
    config.vel_b0 = ChVector3d(0.0, 0, 0);
    config.dt = 5.0e-4;
    config.t_end = 0.5;
    return config;
}

SpherePairConfig HeadonSpheresMassRatioConfig() {
    SpherePairConfig config = HeadonSpheresConfig();
    config.case_name = "headon_spheres_mass_ratio";
    config.asset = "assets/headon_spheres_mass_ratio";
    config.density_a = 1000.0;
    config.density_b = 2000.0;
    config.vel_a0 = ChVector3d(0.0, 0, 0);
    config.vel_b0 = ChVector3d(-1.0, 0, 0);
    return config;
}

int RunMode(const std::string& mode) {
    if (mode == "headon_spheres") {
        return RunSpherePairCase(HeadonSpheresConfig());
    }
    if (mode == "headon_spheres_mass_ratio") {
        return RunSpherePairCase(HeadonSpheresMassRatioConfig());
    }
    if (mode == "eccentric_roller") {
        std::cerr << "SDF-NCP benchmark mode 'eccentric_roller' is disabled until a full asset-geometry SDF "
                     "benchmark is implemented.\n";
        return 2;
    }
    if (mode == "cam") {
        std::cerr << "SDF-NCP benchmark mode 'cam' is disabled until a full asset-geometry SDF benchmark is "
                     "implemented.\n";
        return 2;
    }
    if (mode == "simple_gear") {
        std::cerr << "SDF-NCP benchmark mode 'simple_gear' is disabled until full tooth mesh/SDF geometry is "
                     "implemented.\n";
        return 2;
    }
    if (mode == "all") {
        int rc = 0;
        rc |= RunSpherePairCase(HeadonSpheresConfig());
        rc |= RunSpherePairCase(HeadonSpheresMassRatioConfig());
        return rc;
    }
    std::cerr << "Unknown SDF-NCP benchmark mode: " << mode << std::endl;
    return 2;
}

}  // namespace

int main(int argc, char* argv[]) {
    const std::string mode = argc > 1 ? argv[1] : "all";
    return RunMode(mode);
}
