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
// Milestone 2B.6: Re-integrate 2B.5 cleaned contact logic with OpenVDB SDF
//
// This demo verifies:
//   1. The 2B.5 activation-band / force-band separation works with OpenVDB SDF
//   2. OpenVDB narrow-band truncation and gradient behavior don't break geometry
//   3. Parameter sweep results for OpenVDB are comparable to analytical reference
//
// Scene:
//   - Static target: OpenVDB sphere level set (center at origin, configurable radius)
//   - Dynamic body: small sphere dropped from above, resting on top of static sphere
//
// Output:
//   out/milestone_03_6/sdf_point_contact_openvdb_sweep.csv
//   out/milestone_03_6/sdf_point_contact_openvdb_best_run.csv
//   out/milestone_03_6/sdf_point_contact_reference_vs_openvdb.csv
// =============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <filesystem>

// -- Chrono includes --
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/core/ChFrame.h"

// -- OpenVDB includes --
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/math/Stencils.h>

using namespace chrono;

// =============================================================================
// SDF query interface
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
// Sweep result
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
    double min_phi_final;
    double mean_active_phi_final;
    bool   stable;
    bool   premature_lift;
    bool   deep_penetration;
};

// =============================================================================
// Generate surface samples on a sphere
// =============================================================================

std::vector<SurfaceSample> GenerateSphereSamples(double radius, int n_theta, int n_phi) {
    std::vector<SurfaceSample> samples;
    for (int i = 0; i < n_theta; i++) {
        double theta = M_PI * (i + 0.5) / n_theta;
        for (int j = 0; j < n_phi; j++) {
            double phi_angle = 2.0 * M_PI * j / n_phi;
            double x = radius * std::sin(theta) * std::cos(phi_angle);
            double y = radius * std::cos(theta);
            double z = radius * std::sin(theta) * std::sin(phi_angle);
            SurfaceSample s;
            s.local_pos = ChVector3d(x, y, z);
            s.area_weight = 1.0 / (n_theta * n_phi);
            samples.push_back(s);
        }
    }
    return samples;
}

// =============================================================================
// Analytical sphere SDF (reference)
// =============================================================================

SDFQueryResult AnalyticSphereQuery(const ChVector3d& pt, const ChVector3d& center, double radius) {
    SDFQueryResult qr;
    ChVector3d d = pt - center;
    double dist = d.Length();
    qr.phi = dist - radius;
    if (dist < 1e-12) {
        qr.grad = ChVector3d(1, 0, 0);
    } else {
        qr.grad = d / dist;
    }
    return qr;
}

// =============================================================================
// OpenVDB sphere SDF query
// =============================================================================

struct OpenVDBSphereSDF {
    openvdb::FloatGrid::Ptr grid;
    ChVector3d center_world;
    double sphere_radius;
    double voxel_size;
    double half_width_voxels;
    double band_width_world;

    OpenVDBSphereSDF() : sphere_radius(0.0), voxel_size(0.0), half_width_voxels(0.0), band_width_world(0.0) {}

    SDFQueryResult Query(const ChVector3d& world_pt) const {
        SDFQueryResult qr;

        // Transform world point to SDF local (sphere centered at origin)
        ChVector3d local = world_pt - center_world;

        // Query phi via OpenVDB trilinear interpolation
        openvdb::tools::GridSampler<
            openvdb::FloatGrid::TreeType,
            openvdb::tools::BoxSampler
        > sampler(grid->tree(), grid->transform());
        openvdb::Vec3d local_vec(local.x(), local.y(), local.z());
        qr.phi = static_cast<double>(sampler.wsSample(local_vec));

        // Gradient via central finite differences (in world space)
        double dx = voxel_size * 0.5;
        ChVector3d px = world_pt + ChVector3d(dx, 0, 0);
        ChVector3d mx = world_pt - ChVector3d(dx, 0, 0);
        ChVector3d py = world_pt + ChVector3d(0, dx, 0);
        ChVector3d my = world_pt - ChVector3d(0, dx, 0);
        ChVector3d pz = world_pt + ChVector3d(0, 0, dx);
        ChVector3d mz = world_pt - ChVector3d(0, 0, dx);

        double phi_px = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() + dx, local.y(), local.z())));
        double phi_mx = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() - dx, local.y(), local.z())));
        double phi_py = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() + dx, local.z())));
        double phi_my = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() - dx, local.z())));
        double phi_pz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() + dx)));
        double phi_mz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() - dx)));

        qr.grad = ChVector3d(
            (phi_px - phi_mx) / (2.0 * dx),
            (phi_py - phi_my) / (2.0 * dx),
            (phi_pz - phi_mz) / (2.0 * dx)
        );

        return qr;
    }

    double GetBandWidth() const { return band_width_world; }
};

// =============================================================================
// Helper: ensure output directory exists
// =============================================================================

static std::string GetProjectRoot() {
    auto exe_path = std::filesystem::current_path();
    // We are in build/bin/Release, go up 3 levels
    for (int i = 0; i < 3; i++) {
        exe_path = exe_path.parent_path();
    }
    return exe_path.string();
}

static void EnsureDir(const std::string& dir_path) {
    std::filesystem::create_directories(std::filesystem::path(dir_path));
}

// =============================================================================
// Run one simulation with given parameters (shared core, same as 2B.5)
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
    double final_min_phi = 1e10;
    double final_mean_active_phi = 0.0;
    int final_active_count = 0;

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
            final_min_phi = min_phi;
            final_mean_active_phi = (active_count > 0) ? (sum_active_phi / active_count) : 0.0;
            final_active_count = active_count;

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
    result.min_phi_final = final_min_phi;
    result.mean_active_phi_final = final_mean_active_phi;
    result.stable = stable;
    result.premature_lift = (result.final_y > expected_equilibrium_y + 0.05);
    result.deep_penetration = (result.final_y < expected_equilibrium_y - 0.05);

    if (write_csv) {
        csv_file.close();
    }

    return result;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 2B.6: OpenVDB SDF Re-integration ===" << std::endl;

    // Initialize OpenVDB
    openvdb::initialize();

    // -- Simulation parameters --
    double time_step = 5e-4;
    double total_time = 2.0;

    // -- Geometry definition --
    // Static OpenVDB sphere: center at (0, 0, 0), radius = 1.0m
    // Dynamic sphere: radius = 0.2m, dropped from y=3.0
    // Theoretical contact: dynamic sphere center at y = 1.0 + 0.2 = 1.2
    double static_sphere_radius = 1.0;
    double dyn_sphere_radius = 0.2;
    double dyn_sphere_density = 1000.0;
    double dyn_mass = (4.0 / 3.0) * M_PI * std::pow(dyn_sphere_radius, 3) * dyn_sphere_density;
    double expected_equilibrium_y = static_sphere_radius + dyn_sphere_radius;
    double drop_height = 3.0;

    std::cout << "\nGeometry:" << std::endl;
    std::cout << "  Static OpenVDB sphere radius: " << static_sphere_radius << " m" << std::endl;
    std::cout << "  Dynamic sphere radius: " << dyn_sphere_radius << " m" << std::endl;
    std::cout << "  Dynamic sphere mass: " << dyn_mass << " kg" << std::endl;
    std::cout << "  Expected equilibrium Y (center of dynamic sphere): " << expected_equilibrium_y << " m" << std::endl;
    std::cout << "  Expected equilibrium force: " << (dyn_mass * 9.81) << " N" << std::endl;

    // -- Create OpenVDB sphere level set --
    double voxel_size = 0.05;
    float half_width_vox = 3.0f;
    double band_width_world = half_width_vox * voxel_size;

    std::cout << "\nOpenVDB SDF:" << std::endl;
    std::cout << "  Voxel size: " << voxel_size << " m" << std::endl;
    std::cout << "  Half width: " << half_width_vox << " voxels" << std::endl;
    std::cout << "  Band width: +/- " << band_width_world << " m" << std::endl;

    auto sphere_grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        static_cast<float>(static_sphere_radius),
        openvdb::Vec3f(0.0f, 0.0f, 0.0f),
        static_cast<float>(voxel_size),
        half_width_vox
    );
    sphere_grid->setGridClass(openvdb::GRID_LEVEL_SET);

    OpenVDBSphereSDF openvdb_sdf;
    openvdb_sdf.grid = sphere_grid;
    openvdb_sdf.center_world = ChVector3d(0, 0, 0);
    openvdb_sdf.sphere_radius = static_sphere_radius;
    openvdb_sdf.voxel_size = voxel_size;
    openvdb_sdf.half_width_voxels = half_width_vox;
    openvdb_sdf.band_width_world = band_width_world;

    // -- Analytical sphere SDF (reference baseline) --
    SDFProbeFunc analytic_sdf = [&](const ChVector3d& pt) -> SDFQueryResult {
        return AnalyticSphereQuery(pt, ChVector3d(0, 0, 0), static_sphere_radius);
    };

    // -- OpenVDB SDF probe function --
    SDFProbeFunc openvdb_sdf_probe = [&](const ChVector3d& pt) -> SDFQueryResult {
        return openvdb_sdf.Query(pt);
    };

    // -- Surface samples --
    int n_theta = 8;
    int n_phi = 16;
    std::vector<SurfaceSample> samples = GenerateSphereSamples(dyn_sphere_radius, n_theta, n_phi);
    std::cout << "  Surface samples: " << samples.size() << " (" << n_theta << "x" << n_phi << ")" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_03_6";
    EnsureDir(out_dir);

    // ====================================================================
    // Phase 1: Reference analytical SDF sweep (quick, best params only)
    // ====================================================================
    std::cout << "\n=== Phase 1: Reference Analytical SDF (best params) ===" << std::endl;

    double ref_best_stiffness = 1e5;
    double ref_best_damping = 500;
    double ref_best_force_band = 0.005;
    double activation_band = 0.1;

    SweepResult ref_result = RunSimulation(
        analytic_sdf, samples,
        ref_best_stiffness, ref_best_damping, activation_band, ref_best_force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        ""  // no CSV for reference, just console output
    );

    std::cout << "  Reference final Y: " << ref_result.final_y << " m" << std::endl;
    std::cout << "  Reference y_error: " << ref_result.y_error << " m" << std::endl;
    std::cout << "  Reference stable: " << (ref_result.stable ? "YES" : "NO") << std::endl;

    // ====================================================================
    // Phase 2: OpenVDB parameter sweep
    // ====================================================================
    std::cout << "\n=== Phase 2: OpenVDB Parameter Sweep ===" << std::endl;

    std::vector<double> stiffness_values = {1e4, 5e4, 1e5, 5e5};
    std::vector<double> damping_values = {100, 200, 500};
    std::vector<double> force_band_values = {-0.01, 0.0, 0.005};

    std::string sweep_csv = out_dir + "/sdf_point_contact_openvdb_sweep.csv";
    std::ofstream sweep_file(sweep_csv);
    if (!sweep_file.is_open()) {
        std::cerr << "Error: Cannot create sweep output file: " << sweep_csv << std::endl;
        return 1;
    }

    sweep_file << "stiffness,damping,force_band,sample_count,"
               << "final_y,expected_y,y_error,normalized_y_error,"
               << "final_force_y,expected_force_y,"
               << "avg_active_count,min_phi_final,mean_active_phi_final,"
               << "stable,premature_lift,deep_penetration,"
               << "band_width_world,voxel_size" << std::endl;

    SweepResult openvdb_best_result;
    openvdb_best_result.y_error = 1e10;
    openvdb_best_result.stable = false;
    double openvdb_best_stiffness = 0, openvdb_best_damping = 0, openvdb_best_force_band = 0;
    int sweep_count = 0;

    for (double stiffness : stiffness_values) {
        for (double damping : damping_values) {
            for (double force_band : force_band_values) {
                SweepResult r = RunSimulation(
                    openvdb_sdf_probe, samples,
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
                           << r.min_phi_final << "," << r.mean_active_phi_final << ","
                           << (r.stable ? 1 : 0) << ","
                           << (r.premature_lift ? 1 : 0) << ","
                           << (r.deep_penetration ? 1 : 0) << ","
                           << band_width_world << "," << voxel_size << std::endl;

                sweep_count++;

                if (r.stable && r.y_error < openvdb_best_result.y_error) {
                    openvdb_best_result = r;
                    openvdb_best_stiffness = stiffness;
                    openvdb_best_damping = damping;
                    openvdb_best_force_band = force_band;
                }
            }
        }
    }
    sweep_file.close();

    std::cout << "  Sweep completed: " << sweep_count << " configurations tested" << std::endl;

    // ====================================================================
    // Phase 3: OpenVDB best run with full output
    // ====================================================================
    std::cout << "\n=== Phase 3: OpenVDB Best Parameter Run ===" << std::endl;
    std::cout << "  Best stiffness: " << openvdb_best_stiffness << " N/m" << std::endl;
    std::cout << "  Best damping: " << openvdb_best_damping << " N*s/m" << std::endl;
    std::cout << "  Best force_band: " << openvdb_best_force_band << " m" << std::endl;

    std::string best_csv = out_dir + "/sdf_point_contact_openvdb_best_run.csv";
    SweepResult openvdb_final_result = RunSimulation(
        openvdb_sdf_probe, samples,
        openvdb_best_stiffness, openvdb_best_damping, activation_band, openvdb_best_force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        best_csv
    );

    // ====================================================================
    // Phase 4: Reference vs OpenVDB comparison at best parameters
    // ====================================================================
    std::cout << "\n=== Phase 4: Reference vs OpenVDB Comparison ===" << std::endl;

    // Use the OpenVDB best parameters for both
    double cmp_stiffness = openvdb_best_stiffness;
    double cmp_damping = openvdb_best_damping;
    double cmp_force_band = openvdb_best_force_band;

    SweepResult cmp_ref = RunSimulation(
        analytic_sdf, samples,
        cmp_stiffness, cmp_damping, activation_band, cmp_force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        ""
    );

    SweepResult cmp_ovdb = openvdb_final_result;

    std::string compare_csv = out_dir + "/sdf_point_contact_reference_vs_openvdb.csv";
    std::ofstream compare_file(compare_csv);
    if (compare_file.is_open()) {
        compare_file << "metric,analytical_reference,openvdb_sdf,difference" << std::endl;
        compare_file << std::fixed << std::setprecision(6);
        compare_file << "final_y," << cmp_ref.final_y << "," << cmp_ovdb.final_y << "," << (cmp_ovdb.final_y - cmp_ref.final_y) << std::endl;
        compare_file << "expected_y," << cmp_ref.expected_y << "," << cmp_ovdb.expected_y << ",0" << std::endl;
        compare_file << "y_error," << cmp_ref.y_error << "," << cmp_ovdb.y_error << "," << std::abs(cmp_ovdb.y_error - cmp_ref.y_error) << std::endl;
        compare_file << "normalized_y_error," << (cmp_ref.y_error / dyn_sphere_radius) << "," << (cmp_ovdb.y_error / dyn_sphere_radius) << "," << std::abs(cmp_ovdb.y_error - cmp_ref.y_error) / dyn_sphere_radius << std::endl;
        compare_file << "avg_active_count," << cmp_ref.avg_active_count << "," << cmp_ovdb.avg_active_count << "," << std::abs(cmp_ovdb.avg_active_count - cmp_ref.avg_active_count) << std::endl;
        compare_file << "min_phi_final," << cmp_ref.min_phi_final << "," << cmp_ovdb.min_phi_final << "," << std::abs(cmp_ovdb.min_phi_final - cmp_ref.min_phi_final) << std::endl;
        compare_file << "mean_active_phi," << cmp_ref.mean_active_phi_final << "," << cmp_ovdb.mean_active_phi_final << "," << std::abs(cmp_ovdb.mean_active_phi_final - cmp_ref.mean_active_phi_final) << std::endl;
        compare_file << "stable," << (cmp_ref.stable ? 1 : 0) << "," << (cmp_ovdb.stable ? 1 : 0) << ",0" << std::endl;
        compare_file << "premature_lift," << (cmp_ref.premature_lift ? 1 : 0) << "," << (cmp_ovdb.premature_lift ? 1 : 0) << ",0" << std::endl;
        compare_file << "deep_penetration," << (cmp_ref.deep_penetration ? 1 : 0) << "," << (cmp_ovdb.deep_penetration ? 1 : 0) << ",0" << std::endl;
        compare_file << "stiffness," << cmp_stiffness << "," << cmp_stiffness << ",0" << std::endl;
        compare_file << "damping," << cmp_damping << "," << cmp_damping << ",0" << std::endl;
        compare_file << "force_band," << cmp_force_band << "," << cmp_force_band << ",0" << std::endl;
        compare_file << "activation_band," << activation_band << "," << activation_band << ",0" << std::endl;
        compare_file.close();
    }

    // ====================================================================
    // Summary
    // ====================================================================
    std::cout << "\n=== Final Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n--- Analytical Reference ---" << std::endl;
    std::cout << "  Final Y: " << cmp_ref.final_y << " m" << std::endl;
    std::cout << "  Y error: " << cmp_ref.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (cmp_ref.stable ? "YES" : "NO") << std::endl;

    std::cout << "\n--- OpenVDB SDF ---" << std::endl;
    std::cout << "  Final Y: " << cmp_ovdb.final_y << " m" << std::endl;
    std::cout << "  Y error: " << cmp_ovdb.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (cmp_ovdb.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Min phi (final): " << cmp_ovdb.min_phi_final << " m" << std::endl;
    std::cout << "  Mean active phi: " << cmp_ovdb.mean_active_phi_final << " m" << std::endl;

    std::cout << "\n--- Difference (OpenVDB vs Reference) ---" << std::endl;
    std::cout << "  Y diff: " << (cmp_ovdb.final_y - cmp_ref.final_y) << " m" << std::endl;
    std::cout << "  Y error diff: " << std::abs(cmp_ovdb.y_error - cmp_ref.y_error) << " m" << std::endl;

    bool openvdb_pass = cmp_ovdb.stable && cmp_ovdb.y_error < 0.02;
    bool reference_match = std::abs(cmp_ovdb.y_error - cmp_ref.y_error) < 0.01;

    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "  OpenVDB stable && y_error < 0.02: " << (openvdb_pass ? "PASS" : "FAIL") << std::endl;
    std::cout << "  OpenVDB matches reference within 0.01: " << (reference_match ? "PASS" : "FAIL") << std::endl;

    bool overall_pass = openvdb_pass;
    std::cout << "\n  Overall: " << (overall_pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_point_contact_openvdb_sweep.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_point_contact_openvdb_best_run.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_point_contact_reference_vs_openvdb.csv" << std::endl;

    return overall_pass ? 0 : 1;
}
