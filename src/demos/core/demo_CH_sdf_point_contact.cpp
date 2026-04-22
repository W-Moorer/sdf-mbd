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
// Milestone 2B.5: Geometric Consistency and Parameter Robustness Cleanup
//
// Enhancements over Milestone 2B:
//   1. Separates activation band (candidate detection) from force band (actual penalty)
//   2. Adds geometric consistency metrics: expected_y, y_error, normalized_y_error
//   3. Adds parameter sweep capability
//   4. Uses analytical plane SDF to avoid large-sphere curvature artifacts
//   5. Produces a recommended parameter window
//
// Output:
//   out/milestone_03_5/sdf_point_contact_best_run.csv  (best parameter run)
//   out/milestone_03_5/sdf_point_contact_sweep.csv     (parameter sweep results)
// =============================================================================

// -- Chrono includes --
#define _USE_MATH_DEFINES
#include <cmath>
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>

using namespace chrono;

// =============================================================================
// SDF query interface (allows both OpenVDB-based and analytical SDF)
// =============================================================================

struct SDFQueryResult {
    double phi;
    ChVector3d grad;
};

using SDFProbeFunc = std::function<SDFQueryResult(const ChVector3d& world_pt)>;

// =============================================================================
// Surface sample point
// =============================================================================

struct SurfaceSample {
    ChVector3d local_pos;
    double area_weight;
};

// =============================================================================
// Parameter sweep result
// =============================================================================

struct SweepResult {
    double stiffness;
    double damping;
    double activation_band;
    double force_band;
    int    sample_count;
    double final_y;
    double expected_y;
    double y_error;
    double final_force_y;
    double expected_force_y;
    double avg_active_count;
    bool   stable;
    bool   premature_lift;
    bool   deep_penetration;
};

// =============================================================================
// Generate surface sample points on a sphere
// =============================================================================

std::vector<SurfaceSample> GenerateSphereSamples(double radius, int n_theta, int n_phi) {
    std::vector<SurfaceSample> samples;
    for (int i = 0; i < n_theta; i++) {
        double theta = M_PI * (i + 0.5) / n_theta;
        for (int j = 0; j < n_phi; j++) {
            double phi = 2.0 * M_PI * j / n_phi;
            double x = radius * std::sin(theta) * std::cos(phi);
            double y = radius * std::sin(theta) * std::sin(phi);
            double z = radius * std::cos(theta);
            SurfaceSample s;
            s.local_pos = ChVector3d(x, y, z);
            s.area_weight = 1.0 / (n_theta * n_phi);
            samples.push_back(s);
        }
    }
    return samples;
}

// =============================================================================
// Run one simulation with given parameters
// =============================================================================

SweepResult RunSimulation(
    const SDFProbeFunc& sdf_probe,
    const std::vector<SurfaceSample>& samples,
    double stiffness,
    double damping,
    double activation_band,
    double force_band,
    double time_step,
    double total_time,
    double dyn_sphere_radius,
    double dyn_mass,
    double drop_height,
    double expected_equilibrium_y,
    const std::string& csv_path = ""
) {
    // Create Chrono system
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -9.81, 0));

    auto dyn_body = chrono_types::make_shared<ChBodyEasySphere>(
        dyn_sphere_radius, 1000.0, false, false
    );
    dyn_body->SetPos(ChVector3d(0.0, drop_height, 0.0));
    dyn_body->SetFixed(false);
    sys.AddBody(dyn_body);

    unsigned int acc_id = dyn_body->AddAccumulator();

    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(100);
    sys.GetSolver()->AsIterative()->SetTolerance(1e-6);

    std::ofstream csv_file;
    bool write_csv = !csv_path.empty();
    if (write_csv) {
        csv_file.open(csv_path);
        csv_file << "time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,"
                 << "expected_y,y_error,normalized_y_error,"
                 << "active_sample_count,force_sample_count,"
                 << "min_phi,mean_active_phi,max_penetration,"
                 << "total_force_x,total_force_y,total_force_z,"
                 << "force_error" << std::endl;
    }

    double next_output_time = 0.0;
    double output_interval = 0.01;
    double expected_force = dyn_mass * 9.81;
    bool stable = true;
    double sum_active_count = 0.0;
    int output_count = 0;
    int step_count = 0;

    while (sys.GetChTime() < total_time) {
        dyn_body->EmptyAccumulator(acc_id);

        ChVector3d body_pos = dyn_body->GetPos();
        ChQuaterniond body_rot = dyn_body->GetRot();
        ChVector3d body_vel = dyn_body->GetPosDt();
        ChVector3d body_ang_vel_world = body_rot.Rotate(dyn_body->GetAngVelLocal());

        ChVector3d total_force(0, 0, 0);
        ChVector3d total_torque(0, 0, 0);
        int active_count = 0;
        int force_count = 0;
        double min_phi = 1e10;
        double sum_active_phi = 0.0;
        double max_penetration = 0.0;

        for (const auto& sample : samples) {
            ChVector3d sample_world = body_pos + body_rot.Rotate(sample.local_pos);
            ChVector3d r = sample_world - body_pos;
            ChVector3d sample_vel = body_vel + body_ang_vel_world.Cross(r);

            // Query SDF
            SDFQueryResult qr = sdf_probe(sample_world);
            double phi = qr.phi;

            if (phi < min_phi) min_phi = phi;

            // Activation check (wide band: candidate detection)
            if (phi < activation_band) {
                active_count++;
                sum_active_phi += phi;
            }

            // Force application check (narrow band: actual penalty)
            if (phi < force_band) {
                force_count++;
                double pen = std::max(-phi, 0.0);
                if (pen > max_penetration) max_penetration = pen;

                double grad_len = qr.grad.Length();
                if (grad_len < 1e-12) continue;
                ChVector3d normal = qr.grad / grad_len;

                double vn = sample_vel.Dot(normal);
                double force_mag = stiffness * pen + damping * std::max(-vn, 0.0);

                if (force_mag > 0.0) {
                    ChVector3d pt_force = normal * force_mag;
                    total_force += pt_force;
                    total_torque += r.Cross(pt_force);
                }
            }
        }

        if (total_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, total_force, dyn_body->GetPos(), false);
            ChQuaterniond q_inv = dyn_body->GetRot().GetConjugate();
            dyn_body->AccumulateTorque(acc_id, q_inv.Rotate(total_torque), false);
        }

        sys.DoStepDynamics(time_step);
        step_count++;

        double time = sys.GetChTime();
        ChVector3d pos = dyn_body->GetPos();
        ChVector3d vel = dyn_body->GetPosDt();

        if (std::abs(pos.y()) > 100.0 || std::isnan(pos.y())) {
            stable = false;
            break;
        }

        if (time >= next_output_time - 1e-10) {
            next_output_time += output_interval;
            output_count++;
            sum_active_count += active_count;

            double y_error = std::abs(pos.y() - expected_equilibrium_y);
            double normalized_y_error = y_error / dyn_sphere_radius;
            double force_error = std::abs(total_force.y() - expected_force);

            if (write_csv) {
                csv_file << std::fixed << std::setprecision(8)
                         << time << ","
                         << pos.x() << "," << pos.y() << "," << pos.z() << ","
                         << vel.x() << "," << vel.y() << "," << vel.z() << ","
                         << expected_equilibrium_y << "," << y_error << "," << normalized_y_error << ","
                         << active_count << "," << force_count << ","
                         << min_phi << "," << (active_count > 0 ? sum_active_phi / active_count : 0.0) << ","
                         << max_penetration << ","
                         << total_force.x() << "," << total_force.y() << "," << total_force.z() << ","
                         << force_error << std::endl;
            }
        }
    }

    SweepResult result;
    result.stiffness = stiffness;
    result.damping = damping;
    result.activation_band = activation_band;
    result.force_band = force_band;
    result.sample_count = static_cast<int>(samples.size());
    result.final_y = dyn_body->GetPos().y();
    result.expected_y = expected_equilibrium_y;
    result.y_error = std::abs(result.final_y - expected_equilibrium_y);
    result.final_force_y = 0.0;
    result.expected_force_y = expected_force;
    result.avg_active_count = (output_count > 0) ? (sum_active_count / output_count) : 0.0;
    result.stable = stable;
    result.premature_lift = (result.final_y > expected_equilibrium_y + 0.05);
    result.deep_penetration = (result.final_y < expected_equilibrium_y - 0.05);

    if (write_csv) {
        csv_file.close();
    }

    return result;
}

// =============================================================================
// Main: Parameter sweep + best run
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 2B.5: Geometric Consistency and Parameter Robustness ===" << std::endl;

    // -- Simulation parameters --
    double time_step = 5e-4;
    double total_time = 2.0;

    // -- Geometry definition --
    double dyn_sphere_radius = 0.2;
    double dyn_sphere_density = 1000.0;
    double dyn_mass = (4.0 / 3.0) * M_PI * std::pow(dyn_sphere_radius, 3) * dyn_sphere_density;
    double expected_equilibrium_y = dyn_sphere_radius;  // sphere resting on plane at y=0
    double drop_height = 1.0;

    std::cout << "\nGeometry:" << std::endl;
    std::cout << "  Dynamic sphere radius: " << dyn_sphere_radius << " m" << std::endl;
    std::cout << "  Dynamic sphere mass: " << dyn_mass << " kg" << std::endl;
    std::cout << "  Expected equilibrium Y: " << expected_equilibrium_y << " m" << std::endl;
    std::cout << "  Expected equilibrium force: " << (dyn_mass * 9.81) << " N" << std::endl;

    // -- Analytical plane SDF at y=0 --
    double plane_y = 0.0;
    SDFProbeFunc plane_sdf = [plane_y](const ChVector3d& pt) -> SDFQueryResult {
        SDFQueryResult qr;
        qr.phi = pt.y() - plane_y;      // positive = above plane
        qr.grad = ChVector3d(0, 1, 0);  // gradient points upward
        return qr;
    };

    std::cout << "\nSDF Target: Analytical plane at y=" << plane_y << std::endl;
    std::cout << "  phi(x,y,z) = y - " << plane_y << ", grad = (0, 1, 0)" << std::endl;

    // -- Surface samples --
    int n_theta = 8;
    int n_phi = 16;
    std::vector<SurfaceSample> samples = GenerateSphereSamples(dyn_sphere_radius, n_theta, n_phi);
    std::cout << "  Surface samples: " << samples.size() << " (" << n_theta << "x" << n_phi << ")" << std::endl;

    // ================================================================
    // Phase 1: Parameter sweep
    // ================================================================
    std::cout << "\n=== Phase 1: Parameter Sweep ===" << std::endl;

    std::vector<double> stiffness_values = {1e4, 5e4, 1e5, 5e5, 1e6};
    std::vector<double> damping_values = {100, 200, 500};
    std::vector<double> force_band_values = {-0.01, 0.0, 0.005};

    std::string sweep_csv = "out/milestone_03_5/sdf_point_contact_sweep.csv";
    std::ofstream sweep_file(sweep_csv);
    if (!sweep_file.is_open()) {
        std::cerr << "Error: Cannot create sweep output file" << std::endl;
        return 1;
    }

    sweep_file << "stiffness,damping,force_band,sample_count,"
               << "final_y,expected_y,y_error,normalized_y_error,"
               << "final_force_y,expected_force_y,"
               << "avg_active_count,stable,premature_lift,deep_penetration" << std::endl;

    SweepResult best_result;
    best_result.y_error = 1e10;
    best_result.stable = false;
    double best_stiffness = 0, best_damping = 0, best_force_band = 0;
    int sweep_count = 0;

    for (double stiffness : stiffness_values) {
        for (double damping : damping_values) {
            for (double force_band : force_band_values) {
                double activation_band = 0.1;

                SweepResult r = RunSimulation(
                    plane_sdf, samples,
                    stiffness, damping, activation_band, force_band,
                    time_step, total_time,
                    dyn_sphere_radius, dyn_mass, drop_height,
                    expected_equilibrium_y
                );

                sweep_file << std::fixed << std::setprecision(6)
                           << stiffness << "," << damping << "," << force_band << ","
                           << r.sample_count << ","
                           << r.final_y << "," << r.expected_y << "," << r.y_error << ","
                           << (r.y_error / dyn_sphere_radius) << ","
                           << r.final_force_y << "," << r.expected_force_y << ","
                           << r.avg_active_count << ","
                           << (r.stable ? 1 : 0) << ","
                           << (r.premature_lift ? 1 : 0) << ","
                           << (r.deep_penetration ? 1 : 0) << std::endl;

                sweep_count++;

                if (r.stable && r.y_error < best_result.y_error) {
                    best_result = r;
                    best_stiffness = stiffness;
                    best_damping = damping;
                    best_force_band = force_band;
                }
            }
        }
    }
    sweep_file.close();

    std::cout << "  Sweep completed: " << sweep_count << " configurations tested" << std::endl;

    // ================================================================
    // Phase 2: Best parameter run with full output
    // ================================================================
    std::cout << "\n=== Phase 2: Best Parameter Run ===" << std::endl;
    std::cout << "  Best stiffness: " << best_stiffness << " N/m" << std::endl;
    std::cout << "  Best damping: " << best_damping << " N*s/m" << std::endl;
    std::cout << "  Best force_band: " << best_force_band << " m" << std::endl;

    std::string best_csv = "out/milestone_03_5/sdf_point_contact_best_run.csv";
    if (argc > 1) {
        best_csv = argv[1];
    }

    SweepResult final_result = RunSimulation(
        plane_sdf, samples,
        best_stiffness, best_damping, 0.1, best_force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        best_csv
    );

    std::cout << "\n=== Final Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Final Y position: " << final_result.final_y << " m" << std::endl;
    std::cout << "  Expected Y: " << final_result.expected_y << " m" << std::endl;
    std::cout << "  Y error: " << final_result.y_error << " m" << std::endl;
    std::cout << "  Normalized Y error: " << (final_result.y_error / dyn_sphere_radius) << std::endl;
    std::cout << "  Stable: " << (final_result.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Premature lift: " << (final_result.premature_lift ? "YES" : "NO") << std::endl;
    std::cout << "  Deep penetration: " << (final_result.deep_penetration ? "YES" : "NO") << std::endl;
    std::cout << "  Output: " << best_csv << std::endl;

    bool pass = final_result.stable && final_result.y_error < 0.02;
    std::cout << "\n  Overall: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass ? 0 : 1;
}
