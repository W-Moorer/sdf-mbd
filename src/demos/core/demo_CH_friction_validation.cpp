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
// Friction validation cases for the field-contact primitive runtime.
//
// The cases intentionally use an analytic plane SDF so that the test isolates
// the patch tangential law, Coulomb clamp, and objective tangential-state
// transport from geometry-discretization errors.
//
// Output:
//   out/friction_validation/inclined_plane.csv
//   out/friction_validation/spring_pull.csv
//   out/friction_validation/objective_transport.csv
//   out/friction_validation/summary.csv
//
// =============================================================================

#define _USE_MATH_DEFINES

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double DegToRad(double deg) {
    return deg * kPi / 180.0;
}

double Clamp(double value, double lo, double hi) {
    return std::max(lo, std::min(hi, value));
}

std::string StickSlipName(StickSlipState state) {
    return state == StickSlipState::Stick ? "stick" : "slip";
}

SurfaceGraph MakeRectPatchGraph(int nu, int nv, double size_u, double size_v) {
    SurfaceGraph graph;
    graph.samples.reserve(static_cast<size_t>(nu * nv));

    double du = size_u / static_cast<double>(nu);
    double dv = size_v / static_cast<double>(nv);
    double area = du * dv;

    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < nu; i++) {
            int id = j * nu + i;
            double u = -0.5 * size_u + (static_cast<double>(i) + 0.5) * du;
            double v = -0.5 * size_v + (static_cast<double>(j) + 0.5) * dv;

            SurfaceSample sample;
            sample.id = id;
            sample.local_pos = ChVector3d(u, 0, v);
            sample.area = area;

            if (i > 0) {
                sample.neighbors.push_back(id - 1);
            }
            if (i + 1 < nu) {
                sample.neighbors.push_back(id + 1);
            }
            if (j > 0) {
                sample.neighbors.push_back(id - nu);
            }
            if (j + 1 < nv) {
                sample.neighbors.push_back(id + nu);
            }

            graph.samples.push_back(sample);
        }
    }

    return graph;
}

std::vector<FieldSampleQuery> BuildPlaneQueries(const SurfaceGraph& graph,
                                                const ChVector3d& face_center_on_plane,
                                                const ChVector3d& axis_u,
                                                const ChVector3d& axis_v,
                                                const ChVector3d& plane_normal,
                                                double penetration,
                                                const ChVector3d& surface_velocity) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(graph.samples.size());

    ChVector3d n = SafeNormalize(plane_normal, ChVector3d(0, 1, 0));
    for (const auto& sample : graph.samples) {
        ChVector3d world_pos = face_center_on_plane + axis_u * sample.local_pos.x() +
                               axis_v * sample.local_pos.z() - n * penetration;

        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = surface_velocity;
        query.grad = n;
        query.phi = world_pos.Dot(n);
        queries.push_back(query);
    }

    return queries;
}

FieldContactRuntimeSettings MakeRuntimeSettings(double dt,
                                                double mu,
                                                double normal_pressure_stiffness,
                                                double tangential_stiffness,
                                                double tangential_damping,
                                                double activation_band) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = activation_band;
    settings.extraction.min_area = 0.0;
    settings.extraction.min_samples = 1;
    settings.extraction.use_penetration_weighted_center = true;

    settings.normal.stiffness = normal_pressure_stiffness;
    settings.normal.damping = 0.0;

    settings.tangential.stiffness = tangential_stiffness;
    settings.tangential.damping = tangential_damping;
    settings.tangential.friction_coefficient = mu;
    settings.tangential.time_step = dt;

    settings.inheritance.min_overlap = 0.01;
    settings.inheritance.min_normal_dot = 0.2;
    settings.inheritance.max_center_distance = 0.05;
    settings.inheritance.geometry_fallback_weight = 0.25;

    return settings;
}

struct InclinedMetrics {
    double stick_max_abs_s = 0.0;
    double stick_max_abs_v = 0.0;
    double stick_max_tangential_ratio = 0.0;
    double slip_mean_acceleration = 0.0;
    double slip_expected_acceleration = 0.0;
    double slip_acceleration_rel_error = 0.0;
    double slip_max_tangential_ratio = 0.0;
    int stick_slip_frames = 0;
    int stick_slip_frames_after_settle = 0;
    int slip_slip_frames = 0;
    double runtime_seconds = 0.0;
};

struct InclinedState {
    double s = 0.0;
    double v = 0.0;
    int slip_frames = 0;
    double max_abs_s = 0.0;
    double max_abs_v = 0.0;
    double max_ratio = 0.0;
    double acceleration_sum = 0.0;
    int acceleration_samples = 0;
};

void RunInclinedPlaneCase(const std::string& out_dir, InclinedMetrics& metrics) {
    const double g = 9.80665;
    const double mass = 1.0;
    const double mu = 0.4;
    const double dt = 1.0e-4;
    const double total_time = 1.0;
    const double normal_pressure_stiffness = 1.0e7;
    const double tangential_stiffness = 5.0e4;
    const double tangential_damping = 6.0;
    const double patch_size = 0.1;
    const double patch_area = patch_size * patch_size;

    SurfaceGraph graph = MakeRectPatchGraph(9, 9, patch_size, patch_size);
    std::ofstream csv(out_dir + "/inclined_plane.csv");
    csv << std::setprecision(12);
    csv << "case,time,s,v,acceleration,expected_acceleration,contact_tangent_force,normal_force,"
           "tangential_ratio,state\n";

    auto simulate = [&](const std::string& name, double angle_deg, bool expect_slip) {
        double theta = DegToRad(angle_deg);
        ChVector3d plane_normal(0, std::cos(theta), std::sin(theta));
        ChVector3d axis_u(1, 0, 0);
        ChVector3d downslope(0, -std::sin(theta), std::cos(theta));
        plane_normal = SafeNormalize(plane_normal, ChVector3d(0, 1, 0));
        downslope = SafeNormalize(downslope, ChVector3d(0, 0, 1));

        double normal_load = mass * g * std::cos(theta);
        double penetration = normal_load / (normal_pressure_stiffness * patch_area);
        double expected_acc = expect_slip ? g * (std::sin(theta) - mu * std::cos(theta)) : 0.0;

        FieldContactPrimitiveTracker tracker;
        tracker.Reset();
        FieldContactRuntimeSettings settings =
            MakeRuntimeSettings(dt, mu, normal_pressure_stiffness, tangential_stiffness,
                                tangential_damping, 5.0 * penetration);

        InclinedState state;
        int steps = static_cast<int>(std::round(total_time / dt));
        for (int step = 0; step <= steps; step++) {
            double time = static_cast<double>(step) * dt;
            ChVector3d face_center = downslope * state.s;
            ChVector3d velocity = downslope * state.v;
            auto queries = BuildPlaneQueries(graph, face_center, axis_u, downslope, plane_normal,
                                             penetration, velocity);
            ChVector3d torque_reference = face_center + plane_normal * 0.05;
            FieldContactStepResult contact = tracker.Evaluate(graph, queries, torque_reference, settings);

            double contact_tangent_force = contact.total_force.Dot(downslope);
            double normal_force = contact.total_force.Dot(plane_normal);
            double gravity_tangent_force = mass * g * std::sin(theta);
            double acceleration = (gravity_tangent_force + contact_tangent_force) / mass;

            StickSlipState patch_state = StickSlipState::Stick;
            double ratio = 0.0;
            if (!contact.patches.empty()) {
                patch_state = contact.patches.front().tangential.state;
                ratio = contact.patches.front().tangential_force_ratio;
            }
            if (patch_state == StickSlipState::Slip) {
                state.slip_frames++;
                if (!expect_slip && time > 0.05) {
                    metrics.stick_slip_frames_after_settle++;
                }
            }
            state.max_ratio = std::max(state.max_ratio, ratio);
            state.max_abs_s = std::max(state.max_abs_s, std::abs(state.s));
            state.max_abs_v = std::max(state.max_abs_v, std::abs(state.v));

            if (expect_slip && time > 0.20) {
                state.acceleration_sum += acceleration;
                state.acceleration_samples++;
            }

            csv << name << "," << time << "," << state.s << "," << state.v << ","
                << acceleration << "," << expected_acc << "," << contact_tangent_force << ","
                << normal_force << "," << ratio << "," << StickSlipName(patch_state) << "\n";

            state.v += acceleration * dt;
            state.s += state.v * dt;
        }

        if (expect_slip) {
            metrics.slip_mean_acceleration =
                state.acceleration_samples > 0 ? state.acceleration_sum / state.acceleration_samples : 0.0;
            metrics.slip_expected_acceleration = expected_acc;
            metrics.slip_acceleration_rel_error =
                std::abs(metrics.slip_mean_acceleration - expected_acc) / std::max(std::abs(expected_acc), 1.0e-12);
            metrics.slip_max_tangential_ratio = state.max_ratio;
            metrics.slip_slip_frames = state.slip_frames;
        } else {
            metrics.stick_max_abs_s = state.max_abs_s;
            metrics.stick_max_abs_v = state.max_abs_v;
            metrics.stick_max_tangential_ratio = state.max_ratio;
            metrics.stick_slip_frames = state.slip_frames;
        }
    };

    auto start = std::chrono::steady_clock::now();
    simulate("stick_15deg", 15.0, false);
    simulate("slip_30deg", 30.0, true);
    auto end = std::chrono::steady_clock::now();
    metrics.runtime_seconds = std::chrono::duration<double>(end - start).count();
}

struct SpringPullMetrics {
    double expected_slip_time = 0.0;
    double first_slip_time = -1.0;
    double first_slip_time_error = 0.0;
    double max_tangential_ratio = 0.0;
    double mean_slip_force_ratio = 0.0;
    double final_x = 0.0;
    int slip_frames = 0;
    double runtime_seconds = 0.0;
};

void RunSpringPullCase(const std::string& out_dir, SpringPullMetrics& metrics) {
    const double mass = 1.0;
    const double mu = 0.4;
    const double normal_load = 10.0;
    const double dt = 1.0e-4;
    const double total_time = 3.0;
    const double driver_velocity = 0.02;
    const double drive_stiffness = 100.0;
    const double normal_pressure_stiffness = 1.0e7;
    const double tangential_stiffness = 2.0e4;
    const double tangential_damping = 2.0;
    const double patch_size = 0.1;
    const double patch_area = patch_size * patch_size;
    const double penetration = normal_load / (normal_pressure_stiffness * patch_area);

    metrics.expected_slip_time = (mu * normal_load) / (drive_stiffness * driver_velocity);

    SurfaceGraph graph = MakeRectPatchGraph(9, 9, patch_size, patch_size);
    ChVector3d normal(0, 1, 0);
    ChVector3d axis_x(1, 0, 0);
    ChVector3d axis_z(0, 0, 1);
    FieldContactPrimitiveTracker tracker;
    tracker.Reset();
    FieldContactRuntimeSettings settings =
        MakeRuntimeSettings(dt, mu, normal_pressure_stiffness, tangential_stiffness,
                            tangential_damping, 5.0 * penetration);

    std::ofstream csv(out_dir + "/spring_pull.csv");
    csv << std::setprecision(12);
    csv << "time,driver_x,x,v,drive_force,contact_x,normal_force,tangential_ratio,state\n";

    double x = 0.0;
    double v = 0.0;
    double slip_force_ratio_sum = 0.0;
    int slip_force_ratio_samples = 0;

    auto start = std::chrono::steady_clock::now();
    int steps = static_cast<int>(std::round(total_time / dt));
    for (int step = 0; step <= steps; step++) {
        double time = static_cast<double>(step) * dt;
        double driver_x = driver_velocity * time;
        double drive_force = drive_stiffness * (driver_x - x);

        ChVector3d face_center(x, 0, 0);
        ChVector3d velocity(v, 0, 0);
        auto queries = BuildPlaneQueries(graph, face_center, axis_x, axis_z, normal, penetration, velocity);
        ChVector3d torque_reference = face_center + normal * 0.05;
        FieldContactStepResult contact = tracker.Evaluate(graph, queries, torque_reference, settings);

        double contact_x = contact.total_force.Dot(axis_x);
        double normal_force = contact.total_force.Dot(normal);
        double acceleration = (drive_force + contact_x) / mass;

        StickSlipState state = StickSlipState::Stick;
        double ratio = 0.0;
        if (!contact.patches.empty()) {
            state = contact.patches.front().tangential.state;
            ratio = contact.patches.front().tangential_force_ratio;
        }

        if (state == StickSlipState::Slip) {
            metrics.slip_frames++;
            if (metrics.first_slip_time < 0.0) {
                metrics.first_slip_time = time;
            }
            if (time > metrics.expected_slip_time + 0.05) {
                slip_force_ratio_sum += std::abs(contact_x) / std::max(mu * normal_force, 1.0e-12);
                slip_force_ratio_samples++;
            }
        }
        metrics.max_tangential_ratio = std::max(metrics.max_tangential_ratio, ratio);

        csv << time << "," << driver_x << "," << x << "," << v << "," << drive_force << ","
            << contact_x << "," << normal_force << "," << ratio << "," << StickSlipName(state) << "\n";

        v += acceleration * dt;
        x += v * dt;
    }
    auto end = std::chrono::steady_clock::now();

    metrics.mean_slip_force_ratio =
        slip_force_ratio_samples > 0 ? slip_force_ratio_sum / static_cast<double>(slip_force_ratio_samples) : 0.0;
    metrics.first_slip_time_error =
        metrics.first_slip_time >= 0.0 ? std::abs(metrics.first_slip_time - metrics.expected_slip_time) : -1.0;
    metrics.final_x = x;
    metrics.runtime_seconds = std::chrono::duration<double>(end - start).count();
}

struct TransportMetrics {
    double max_minimal_norm_error = 0.0;
    double max_minimal_tangent_error = 0.0;
    double final_untransported_tangent_error = 0.0;
    double final_naive_energy_ratio = 0.0;
    double min_naive_energy_ratio = 1.0;
    double runtime_seconds = 0.0;
};

void RunObjectiveTransportCase(const std::string& out_dir, TransportMetrics& metrics) {
    std::ofstream csv(out_dir + "/objective_transport.csv");
    csv << std::setprecision(12);
    csv << "angle_deg,minimal_norm,minimal_norm_error,minimal_tangent_error,"
           "untransported_tangent_error,projected_naive_norm,projected_naive_energy_ratio\n";

    auto start = std::chrono::steady_clock::now();
    ChVector3d n0(0, 1, 0);
    ChVector3d xi0(0.001, 0, 0);
    ChVector3d xi_minimal = xi0;
    ChVector3d xi_naive = xi0;
    double initial_norm = xi0.Length();

    const int steps = 90;
    for (int i = 0; i <= steps; i++) {
        double angle = DegToRad(static_cast<double>(i));
        ChVector3d n_cur(std::sin(angle), std::cos(angle), 0);
        n_cur = SafeNormalize(n_cur, n0);

        if (i > 0) {
            double prev_angle = DegToRad(static_cast<double>(i - 1));
            ChVector3d n_prev(std::sin(prev_angle), std::cos(prev_angle), 0);
            n_prev = SafeNormalize(n_prev, n0);
            xi_minimal = TransportElasticStateMinimalRotation(xi_minimal, n_prev, n_cur);
            xi_naive = ProjectToTangent(xi_naive, n_cur);
        }

        double minimal_norm = xi_minimal.Length();
        double minimal_norm_error = std::abs(minimal_norm - initial_norm) / initial_norm;
        double tangent_error = std::abs(xi_minimal.Dot(n_cur)) / initial_norm;
        double untransported_tangent_error = std::abs(xi0.Dot(n_cur)) / initial_norm;
        double naive_norm = xi_naive.Length();
        double naive_energy_ratio = (naive_norm * naive_norm) / (initial_norm * initial_norm);

        metrics.max_minimal_norm_error = std::max(metrics.max_minimal_norm_error, minimal_norm_error);
        metrics.max_minimal_tangent_error = std::max(metrics.max_minimal_tangent_error, tangent_error);
        metrics.final_untransported_tangent_error = untransported_tangent_error;
        metrics.min_naive_energy_ratio = std::min(metrics.min_naive_energy_ratio, naive_energy_ratio);
        metrics.final_naive_energy_ratio = naive_energy_ratio;

        csv << static_cast<double>(i) << "," << minimal_norm << "," << minimal_norm_error << ","
            << tangent_error << "," << untransported_tangent_error << "," << naive_norm << ","
            << naive_energy_ratio << "\n";
    }
    auto end = std::chrono::steady_clock::now();
    metrics.runtime_seconds = std::chrono::duration<double>(end - start).count();
}

bool WriteSummary(const std::string& out_dir,
                  const InclinedMetrics& inclined,
                  const SpringPullMetrics& spring,
                  const TransportMetrics& transport) {
    std::ofstream summary(out_dir + "/summary.csv");
    if (!summary) {
        return false;
    }

    summary << std::setprecision(12);
    summary << "case,metric,value,target_or_limit,pass\n";

    auto row = [&](const std::string& case_name, const std::string& metric, double value,
                   double target_or_limit, bool pass) {
        summary << case_name << "," << metric << "," << value << "," << target_or_limit << ","
                << (pass ? "true" : "false") << "\n";
    };

    row("inclined_stick", "max_abs_displacement_m", inclined.stick_max_abs_s, 1.5e-3,
        inclined.stick_max_abs_s < 1.5e-3);
    row("inclined_stick", "slip_frames_after_0p05s", static_cast<double>(inclined.stick_slip_frames_after_settle),
        0.0, inclined.stick_slip_frames_after_settle == 0);
    row("inclined_slip", "mean_acceleration_rel_error", inclined.slip_acceleration_rel_error, 0.03,
        inclined.slip_acceleration_rel_error < 0.03);
    row("inclined_slip", "max_tangential_force_ratio", inclined.slip_max_tangential_ratio, 1.0000001,
        inclined.slip_max_tangential_ratio <= 1.0000001);

    row("spring_pull", "first_slip_time_s", spring.first_slip_time, spring.expected_slip_time,
        spring.first_slip_time > 0.0 && spring.first_slip_time_error < 0.18);
    row("spring_pull", "mean_slip_force_ratio", spring.mean_slip_force_ratio, 1.0,
        spring.mean_slip_force_ratio > 0.97 && spring.mean_slip_force_ratio < 1.0000001);
    row("spring_pull", "max_tangential_force_ratio", spring.max_tangential_ratio, 1.0000001,
        spring.max_tangential_ratio <= 1.0000001);

    row("objective_transport", "max_minimal_norm_error", transport.max_minimal_norm_error, 1.0e-12,
        transport.max_minimal_norm_error < 1.0e-12);
    row("objective_transport", "max_minimal_tangent_error", transport.max_minimal_tangent_error, 1.0e-12,
        transport.max_minimal_tangent_error < 1.0e-12);
    row("objective_transport", "final_untransported_tangent_error", transport.final_untransported_tangent_error,
        0.99, transport.final_untransported_tangent_error > 0.99);

    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string out_dir = "out/friction_validation";
    if (argc > 1) {
        out_dir = argv[1];
    }

    std::filesystem::create_directories(out_dir);

    InclinedMetrics inclined;
    SpringPullMetrics spring;
    TransportMetrics transport;

    RunInclinedPlaneCase(out_dir, inclined);
    RunSpringPullCase(out_dir, spring);
    RunObjectiveTransportCase(out_dir, transport);

    if (!WriteSummary(out_dir, inclined, spring, transport)) {
        std::cerr << "Failed to write summary.csv\n";
        return 1;
    }

    std::cout << std::setprecision(8);
    std::cout << "Friction validation output: " << out_dir << "\n";
    std::cout << "Inclined plane runtime: " << inclined.runtime_seconds << " s\n";
    std::cout << "  stick max |s| = " << inclined.stick_max_abs_s
              << ", slip acceleration = " << inclined.slip_mean_acceleration
              << " expected " << inclined.slip_expected_acceleration
              << ", rel err = " << inclined.slip_acceleration_rel_error << "\n";
    std::cout << "Spring pull runtime: " << spring.runtime_seconds << " s\n";
    std::cout << "  first slip time = " << spring.first_slip_time
              << " expected " << spring.expected_slip_time
              << ", mean slip force ratio = " << spring.mean_slip_force_ratio << "\n";
    std::cout << "Objective transport runtime: " << transport.runtime_seconds << " s\n";
    std::cout << "  max minimal norm error = " << transport.max_minimal_norm_error
              << ", final untransported tangent error = " << transport.final_untransported_tangent_error
              << ", final projected-naive energy ratio = " << transport.final_naive_energy_ratio << "\n";

    return 0;
}
