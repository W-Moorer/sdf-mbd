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
// Milestone 4: Patch-level Force Aggregation
//
// This demo upgrades the patch primitive layer from "observational" to "mechanical".
// Three contact modes are implemented for A/B/C comparison:
//   A) Pointwise mode:       sample forces directly accumulated to body (baseline)
//   B) Direct patch mode:    sample forces summed per-patch, then patches applied to body
//   C) Equivalent patch mode: force applied at patch center along patch normal
//
// Output:
//   out/milestone_05/sdf_patch_force_output.csv          (main simulation output)
//   out/milestone_05/sdf_patch_force_summary.csv         (patch-level long table)
//   out/milestone_05/sdf_patch_force_vs_pointwise.csv    (pointwise comparison time series)
//   out/milestone_05/sdf_patch_force_vs_pointwise_summary.csv (comparison metrics)
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
#include <set>

// -- Chrono includes --
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/core/ChFrame.h"

// -- OpenVDB includes --
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Interpolation.h>

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
// Surface sample point (with grid index for adjacency)
// =============================================================================

struct SurfaceSample {
    ChVector3d local_pos;
    double area_weight;
    int theta_index;
    int phi_index;
    int global_index;
};

// =============================================================================
// Patch primitive data structure
// =============================================================================

struct PatchPrimitive {
    int patch_id;
    std::vector<int> sample_ids;       // global indices of active samples in this patch
    int active_sample_count;           // number of active samples
    ChVector3d center;                 // weighted average position
    ChVector3d normal;                 // weighted average normal (normalized)
    double total_area;                 // sum of sample weights
    double min_phi;                    // minimum phi in patch
    double mean_phi;                   // weighted mean phi
    double max_penetration;            // max penetration in patch
    ChVector3d force;                  // sum of pointwise forces in patch
    ChVector3d torque;                 // sum of pointwise torques about body COM
};

// =============================================================================
// Contact force mode
// =============================================================================

enum class ContactMode {
    Pointwise,        // sample forces directly to body
    DirectPatchSum,   // patch_force = sum(sample_forces), patch_torque = sum(sample_torques), then to body
    EquivalentPatch   // single force at patch_center along patch_normal, magnitude = normal component of patch_force
};

// =============================================================================
// Generate surface samples on a sphere with grid indexing
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
            s.theta_index = i;
            s.phi_index = j;
            s.global_index = i * n_phi + j;
            samples.push_back(s);
        }
    }
    return samples;
}

// =============================================================================
// Patch grouping: connected components on spherical sampling grid
// =============================================================================

struct GridAdjacency {
    int n_theta;
    int n_phi;
    int total_samples;

    GridAdjacency() : n_theta(0), n_phi(0), total_samples(0) {}

    GridAdjacency(int nt, int np) : n_theta(nt), n_phi(np), total_samples(nt * np) {}

    std::vector<int> GetNeighbors(int ti, int pj) const {
        std::vector<int> neighbors;
        if (ti > 0) neighbors.push_back((ti - 1) * n_phi + pj);
        if (ti < n_theta - 1) neighbors.push_back((ti + 1) * n_phi + pj);
        int pj_prev = (pj - 1 + n_phi) % n_phi;
        int pj_next = (pj + 1) % n_phi;
        neighbors.push_back(ti * n_phi + pj_prev);
        neighbors.push_back(ti * n_phi + pj_next);
        return neighbors;
    }

    std::vector<std::vector<int>> FindComponents(const std::vector<int>& active_indices) const {
        std::set<int> active_set(active_indices.begin(), active_indices.end());
        std::set<int> visited;
        std::vector<std::vector<int>> components;

        for (int idx : active_indices) {
            if (visited.count(idx)) continue;

            std::vector<int> component;
            std::vector<int> queue;
            queue.push_back(idx);
            visited.insert(idx);

            int head = 0;
            while (head < queue.size()) {
                int current = queue[head++];
                component.push_back(current);

                int ti = current / n_phi;
                int pj = current % n_phi;

                for (int neighbor : GetNeighbors(ti, pj)) {
                    if (active_set.count(neighbor) && !visited.count(neighbor)) {
                        visited.insert(neighbor);
                        queue.push_back(neighbor);
                    }
                }
            }

            if (!component.empty()) {
                components.push_back(component);
            }
        }

        return components;
    }
};

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

        ChVector3d local = world_pt - center_world;

        openvdb::tools::GridSampler<
            openvdb::FloatGrid::TreeType,
            openvdb::tools::BoxSampler
        > sampler(grid->tree(), grid->transform());
        openvdb::Vec3d local_vec(local.x(), local.y(), local.z());
        qr.phi = static_cast<double>(sampler.wsSample(local_vec));

        double dx = voxel_size * 0.5;
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
};

// =============================================================================
// Helper: project root directory
// =============================================================================

static std::string GetProjectRoot() {
    auto exe_path = std::filesystem::current_path();
    for (int i = 0; i < 3; i++) {
        exe_path = exe_path.parent_path();
    }
    return exe_path.string();
}

static void EnsureDir(const std::string& dir_path) {
    std::filesystem::create_directories(std::filesystem::path(dir_path));
}

// =============================================================================
// Build patches from active samples (reused from milestone 4)
// =============================================================================

std::vector<PatchPrimitive> BuildPatches(
    const std::vector<int>& active_indices,
    const std::vector<SurfaceSample>& samples,
    const std::vector<SDFQueryResult>& sdf_results,
    const std::vector<ChVector3d>& sample_world_positions,
    const std::vector<ChVector3d>& sample_forces,
    const std::vector<ChVector3d>& sample_torques,
    const GridAdjacency& adjacency,
    int patch_id_offset
) {
    std::vector<PatchPrimitive> patches;

    if (active_indices.empty()) return patches;

    auto components = adjacency.FindComponents(active_indices);

    for (size_t c = 0; c < components.size(); c++) {
        PatchPrimitive patch;
        patch.patch_id = patch_id_offset + static_cast<int>(c);
        patch.sample_ids = components[c];
        patch.active_sample_count = static_cast<int>(components[c].size());

        double total_weight = 0.0;
        ChVector3d weighted_center(0, 0, 0);
        ChVector3d weighted_normal(0, 0, 0);
        double weighted_phi_sum = 0.0;
        patch.min_phi = 1e10;
        patch.max_penetration = 0.0;
        patch.force = ChVector3d(0, 0, 0);
        patch.torque = ChVector3d(0, 0, 0);
        patch.total_area = 0.0;

        for (int sid : components[c]) {
            const auto& s = samples[sid];
            const auto& qr = sdf_results[sid];
            double w = s.area_weight;
            patch.total_area += w;

            ChVector3d world_pos = sample_world_positions[sid];
            weighted_center += world_pos * w;

            double grad_len = qr.grad.Length();
            if (grad_len > 1e-12) {
                ChVector3d n = qr.grad / grad_len;
                weighted_normal += n * w;
            }

            weighted_phi_sum += qr.phi * w;

            if (qr.phi < patch.min_phi) patch.min_phi = qr.phi;
            double pen = std::max(-qr.phi, 0.0);
            if (pen > patch.max_penetration) patch.max_penetration = pen;

            patch.force += sample_forces[sid];
            patch.torque += sample_torques[sid];
        }

        total_weight = patch.total_area;
        if (total_weight > 1e-12) {
            patch.center = weighted_center / total_weight;
            patch.mean_phi = weighted_phi_sum / total_weight;
        }

        double norm_len = weighted_normal.Length();
        if (norm_len > 1e-12) {
            patch.normal = weighted_normal / norm_len;
        } else {
            patch.normal = ChVector3d(0, 1, 0);
        }

        patches.push_back(patch);
    }

    return patches;
}

// =============================================================================
// Run simulation with configurable contact mode
// =============================================================================

struct SimStats {
    double final_y;
    double expected_y;
    double y_error;
    bool stable;
    double avg_active_count;
    double avg_patch_count;
    double avg_max_patch_samples;
    double center_drift_total;
    double normal_fluctuation_total;
    double force_std;
    double torque_std;
    std::string mode_name;
};

SimStats RunContactSimulation(
    const SDFProbeFunc& sdf_probe,
    const std::vector<SurfaceSample>& samples,
    const GridAdjacency& adjacency,
    ContactMode mode,
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
    const std::string& output_csv = "",
    const std::string& summary_csv = ""
) {
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

    std::ofstream out_file;
    std::ofstream summary_file;

    bool write_output = !output_csv.empty();
    bool write_summary = !summary_csv.empty();

    if (write_output) {
        out_file.open(output_csv);
        out_file << "time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,"
                 << "expected_y,y_error,normalized_y_error,"
                 << "active_sample_count,patch_count,total_patch_area,"
                 << "total_force_x,total_force_y,total_force_z,"
                 << "total_torque_x,total_torque_y,total_torque_z,"
                 << "force_error" << std::endl;
    }

    if (write_summary) {
        summary_file.open(summary_csv);
        summary_file << "time,patch_id,center_x,center_y,center_z,"
                     << "normal_x,normal_y,normal_z,"
                     << "area,min_phi,mean_phi,max_penetration,"
                     << "active_sample_count,"
                     << "force_x,force_y,force_z,"
                     << "torque_x,torque_y,torque_z" << std::endl;
    }

    double expected_force = dyn_mass * 9.81;
    bool stable = true;
    double sum_active_count = 0.0;
    double sum_patch_count = 0.0;
    double sum_max_patch_samples = 0.0;
    int output_count = 0;

    ChVector3d prev_max_patch_center(0, 0, 0);
    ChVector3d prev_max_patch_normal(0, 1, 0);
    double center_drift_total = 0.0;
    double normal_fluctuation_total = 0.0;
    double force_sum = 0.0;
    double force_sq_sum = 0.0;
    double torque_sum = 0.0;
    double torque_sq_sum = 0.0;
    bool first_output = true;

    double next_output_time = 0.0;
    double output_interval = 0.01;

    int global_patch_id_counter = 0;

    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());
    std::vector<ChVector3d> sample_forces(samples.size());
    std::vector<ChVector3d> sample_torques(samples.size());

    while (sys.GetChTime() < total_time) {
        dyn_body->EmptyAccumulator(acc_id);

        ChVector3d body_pos = dyn_body->GetPos();
        ChQuaterniond body_rot = dyn_body->GetRot();
        ChVector3d body_vel = dyn_body->GetPosDt();
        ChVector3d body_ang_vel_world = body_rot.Rotate(dyn_body->GetAngVelLocal());

        // Step 1: Query all samples, compute pointwise forces
        ChVector3d total_force(0, 0, 0);
        ChVector3d total_torque(0, 0, 0);
        std::vector<int> active_indices;
        int force_count = 0;
        double min_phi = 1e10;

        for (size_t si = 0; si < samples.size(); si++) {
            const auto& s = samples[si];
            ChVector3d sample_world = body_pos + body_rot.Rotate(s.local_pos);
            ChVector3d r = sample_world - body_pos;
            ChVector3d sample_vel = body_vel + body_ang_vel_world.Cross(r);

            sample_world_positions[si] = sample_world;
            sdf_results[si] = sdf_probe(sample_world);
            double phi = sdf_results[si].phi;
            if (phi < min_phi) min_phi = phi;

            sample_forces[si] = ChVector3d(0, 0, 0);
            sample_torques[si] = ChVector3d(0, 0, 0);

            if (phi < activation_band) {
                active_indices.push_back(static_cast<int>(si));

                if (phi < force_band) {
                    force_count++;
                    double pen = std::max(-phi, 0.0);

                    double grad_len = sdf_results[si].grad.Length();
                    if (grad_len < 1e-12) continue;
                    ChVector3d normal = sdf_results[si].grad / grad_len;

                    double vn = sample_vel.Dot(normal);
                    double force_mag = stiffness * pen + damping * std::max(-vn, 0.0);

                    if (force_mag > 0.0) {
                        ChVector3d pt_force = normal * force_mag;
                        ChVector3d pt_torque = r.Cross(pt_force);
                        sample_forces[si] = pt_force;
                        sample_torques[si] = pt_torque;

                        // Always accumulate for statistics (even in patch mode)
                        total_force += pt_force;
                        total_torque += pt_torque;
                    }
                }
            }
        }

        // Step 2: Build patches
        std::vector<PatchPrimitive> patches = BuildPatches(
            active_indices, samples, sdf_results,
            sample_world_positions, sample_forces, sample_torques,
            adjacency, global_patch_id_counter
        );

        if (!patches.empty()) {
            global_patch_id_counter += static_cast<int>(patches.size());
        }

        // Step 3: Apply forces based on mode
        ChVector3d apply_force(0, 0, 0);
        ChVector3d apply_torque(0, 0, 0);

        if (mode == ContactMode::Pointwise) {
            // Direct pointwise accumulation
            apply_force = total_force;
            apply_torque = total_torque;
        }
        else if (mode == ContactMode::DirectPatchSum) {
            // Patch-level aggregation: sum forces/torques within each patch,
            // then sum across patches (mathematically identical to pointwise)
            for (const auto& patch : patches) {
                apply_force += patch.force;
                apply_torque += patch.torque;
            }
        }
        else if (mode == ContactMode::EquivalentPatch) {
            // Equivalent patch mode: apply force at patch center along patch normal
            // Magnitude = normal component of patch force
            // Torque = (patch_center - COM) x force
            for (const auto& patch : patches) {
                if (patch.force.Length() < 1e-12) continue;

                // Normal component of patch force
                double normal_force_mag = patch.force.Dot(patch.normal);
                if (normal_force_mag <= 0.0) continue;

                ChVector3d equiv_force = patch.normal * normal_force_mag;
                ChVector3d r_to_center = patch.center - body_pos;
                ChVector3d equiv_torque = r_to_center.Cross(equiv_force);

                apply_force += equiv_force;
                apply_torque += equiv_torque;
            }
        }

        // Apply the accumulated force/torque to the body
        if (apply_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, apply_force, body_pos, false);
            ChQuaterniond q_inv = body_rot.GetConjugate();
            dyn_body->AccumulateTorque(acc_id, q_inv.Rotate(apply_torque), false);
        }

        sys.DoStepDynamics(time_step);

        // Output at regular intervals
        double time = sys.GetChTime();
        if (time >= next_output_time - 1e-10) {
            next_output_time += output_interval;
            output_count++;

            ChVector3d pos = dyn_body->GetPos();
            ChVector3d vel = dyn_body->GetPosDt();

            if (std::abs(pos.y()) > 100.0 || std::isnan(pos.y())) {
                stable = false;
                break;
            }

            double y_error = std::abs(pos.y() - expected_equilibrium_y);
            double normalized_y_error = y_error / dyn_sphere_radius;
            double force_error = std::abs(apply_force.y() - expected_force);

            double patch_count = static_cast<double>(patches.size());
            double total_patch_area = 0.0;
            int max_patch_samples = 0;

            for (const auto& p : patches) {
                total_patch_area += p.total_area;
                if (p.active_sample_count > max_patch_samples) max_patch_samples = p.active_sample_count;
            }

            sum_active_count += static_cast<double>(active_indices.size());
            sum_patch_count += patch_count;
            sum_max_patch_samples += static_cast<double>(max_patch_samples);

            force_sum += apply_force.y();
            force_sq_sum += apply_force.y() * apply_force.y();
            torque_sum += apply_torque.Length();
            torque_sq_sum += apply_torque.Length() * apply_torque.Length();

            // Temporal tracking
            if (!first_output && !patches.empty()) {
                const PatchPrimitive* largest = &patches[0];
                for (const auto& p : patches) {
                    if (p.active_sample_count > largest->active_sample_count) largest = &p;
                }

                double center_drift = (largest->center - prev_max_patch_center).Length();
                center_drift_total += center_drift;

                double normal_diff = 1.0 - std::abs(largest->normal.Dot(prev_max_patch_normal));
                normal_fluctuation_total += normal_diff;

                prev_max_patch_center = largest->center;
                prev_max_patch_normal = largest->normal;
            } else if (!patches.empty()) {
                const PatchPrimitive* largest = &patches[0];
                for (const auto& p : patches) {
                    if (p.active_sample_count > largest->active_sample_count) largest = &p;
                }
                prev_max_patch_center = largest->center;
                prev_max_patch_normal = largest->normal;
            }
            first_output = false;

            // Main output CSV
            if (write_output) {
                out_file << std::fixed << std::setprecision(8)
                         << time << ","
                         << pos.x() << "," << pos.y() << "," << pos.z() << ","
                         << vel.x() << "," << vel.y() << "," << vel.z() << ","
                         << expected_equilibrium_y << "," << y_error << "," << normalized_y_error << ","
                         << active_indices.size() << "," << patch_count << "," << total_patch_area << ","
                         << apply_force.x() << "," << apply_force.y() << "," << apply_force.z() << ","
                         << apply_torque.x() << "," << apply_torque.y() << "," << apply_torque.z() << ","
                         << force_error << std::endl;
            }

            // Patch summary CSV
            if (write_summary) {
                if (patches.empty()) {
                    summary_file << std::fixed << std::setprecision(8)
                                 << time << "," << -1 << ","
                                 << "0,0,0," << "0,0,0,0,0,0,0,0,"
                                 << "0," << "0,0,0," << "0,0,0" << std::endl;
                } else {
                    for (const auto& p : patches) {
                        summary_file << std::fixed << std::setprecision(8)
                                     << time << "," << p.patch_id << ","
                                     << p.center.x() << "," << p.center.y() << "," << p.center.z() << ","
                                     << p.normal.x() << "," << p.normal.y() << "," << p.normal.z() << ","
                                     << p.total_area << "," << p.min_phi << "," << p.mean_phi << "," << p.max_penetration << ","
                                     << p.active_sample_count << ","
                                     << p.force.x() << "," << p.force.y() << "," << p.force.z() << ","
                                     << p.torque.x() << "," << p.torque.y() << "," << p.torque.z() << std::endl;
                    }
                }
            }
        }
    }

    SimStats stats;
    stats.final_y = dyn_body->GetPos().y();
    stats.expected_y = expected_equilibrium_y;
    stats.y_error = std::abs(stats.final_y - expected_equilibrium_y);
    stats.stable = stable;
    stats.avg_active_count = (output_count > 0) ? (sum_active_count / output_count) : 0.0;
    stats.avg_patch_count = (output_count > 0) ? (sum_patch_count / output_count) : 0.0;
    stats.avg_max_patch_samples = (output_count > 0) ? (sum_max_patch_samples / output_count) : 0.0;
    stats.center_drift_total = center_drift_total;
    stats.normal_fluctuation_total = normal_fluctuation_total;
    stats.force_std = 0.0;
    stats.torque_std = 0.0;

    if (output_count > 1) {
        double mean_force = force_sum / output_count;
        stats.force_std = std::sqrt(force_sq_sum / output_count - mean_force * mean_force);
        double mean_torque = torque_sum / output_count;
        stats.torque_std = std::sqrt(torque_sq_sum / output_count - mean_torque * mean_torque);
    }

    switch (mode) {
        case ContactMode::Pointwise: stats.mode_name = "pointwise"; break;
        case ContactMode::DirectPatchSum: stats.mode_name = "direct_patch_sum"; break;
        case ContactMode::EquivalentPatch: stats.mode_name = "equivalent_patch"; break;
    }

    return stats;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 4: Patch-level Force Aggregation ===" << std::endl;

    openvdb::initialize();

    // -- Simulation parameters --
    double time_step = 5e-4;
    double total_time = 2.0;

    // -- Geometry definition --
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
    std::cout << "  Expected equilibrium Y: " << expected_equilibrium_y << " m" << std::endl;

    // -- OpenVDB sphere level set --
    double voxel_size = 0.05;
    float half_width_vox = 3.0f;
    double band_width_world = half_width_vox * voxel_size;

    std::cout << "\nOpenVDB SDF:" << std::endl;
    std::cout << "  Voxel size: " << voxel_size << " m" << std::endl;
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

    SDFProbeFunc openvdb_sdf_probe = [&](const ChVector3d& pt) -> SDFQueryResult {
        return openvdb_sdf.Query(pt);
    };

    // -- Surface samples --
    int n_theta = 8;
    int n_phi = 16;
    std::vector<SurfaceSample> samples = GenerateSphereSamples(dyn_sphere_radius, n_theta, n_phi);
    GridAdjacency adjacency(n_theta, n_phi);
    std::cout << "  Surface samples: " << samples.size() << " (" << n_theta << "x" << n_phi << ")" << std::endl;

    // -- Contact parameters (from 2B.6 best) --
    double stiffness = 1e5;
    double damping = 500;
    double activation_band = 0.1;
    double force_band = 0.0;

    std::cout << "\nContact parameters:" << std::endl;
    std::cout << "  stiffness = " << stiffness << " N/m" << std::endl;
    std::cout << "  damping = " << damping << " N*s/m" << std::endl;
    std::cout << "  activation_band = " << activation_band << " m" << std::endl;
    std::cout << "  force_band = " << force_band << " m" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_05";
    EnsureDir(out_dir);

    // ====================================================================
    // Phase 1: Pointwise mode (A)
    // ====================================================================
    std::cout << "\n=== Phase 1: Pointwise Mode (A) ===" << std::endl;

    SimStats pointwise_stats = RunContactSimulation(
        openvdb_sdf_probe, samples, adjacency,
        ContactMode::Pointwise,
        stiffness, damping, activation_band, force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        out_dir + "/sdf_pointwise_output.csv",
        ""
    );

    std::cout << "  Final Y: " << pointwise_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << pointwise_stats.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (pointwise_stats.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Force std dev: " << pointwise_stats.force_std << " N" << std::endl;
    std::cout << "  Torque std dev: " << pointwise_stats.torque_std << " Nm" << std::endl;

    // ====================================================================
    // Phase 2: Direct patch sum mode (B)
    // ====================================================================
    std::cout << "\n=== Phase 2: Direct Patch Sum Mode (B) ===" << std::endl;

    SimStats direct_patch_stats = RunContactSimulation(
        openvdb_sdf_probe, samples, adjacency,
        ContactMode::DirectPatchSum,
        stiffness, damping, activation_band, force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        out_dir + "/sdf_patch_force_output.csv",
        out_dir + "/sdf_patch_force_summary.csv"
    );

    std::cout << "  Final Y: " << direct_patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << direct_patch_stats.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (direct_patch_stats.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Force std dev: " << direct_patch_stats.force_std << " N" << std::endl;
    std::cout << "  Torque std dev: " << direct_patch_stats.torque_std << " Nm" << std::endl;

    // ====================================================================
    // Phase 3: Equivalent patch mode (C)
    // ====================================================================
    std::cout << "\n=== Phase 3: Equivalent Patch Mode (C) ===" << std::endl;

    SimStats equiv_patch_stats = RunContactSimulation(
        openvdb_sdf_probe, samples, adjacency,
        ContactMode::EquivalentPatch,
        stiffness, damping, activation_band, force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        out_dir + "/sdf_patch_force_equiv_output.csv",
        ""
    );

    std::cout << "  Final Y: " << equiv_patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << equiv_patch_stats.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (equiv_patch_stats.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Force std dev: " << equiv_patch_stats.force_std << " N" << std::endl;
    std::cout << "  Torque std dev: " << equiv_patch_stats.torque_std << " Nm" << std::endl;

    // ====================================================================
    // Phase 4: Generate comparison
    // ====================================================================
    std::cout << "\n=== Phase 4: Comparison ===" << std::endl;

    // Generate pointwise comparison CSV (time series matching)
    std::string pointwise_compare_csv = out_dir + "/sdf_patch_force_vs_pointwise.csv";
    std::ofstream pw_compare(pointwise_compare_csv);
    if (pw_compare.is_open()) {
        pw_compare << "metric,pointwise,direct_patch_sum,equivalent_patch" << std::endl;
        pw_compare << std::fixed << std::setprecision(6);
        pw_compare << "final_y," << pointwise_stats.final_y << "," << direct_patch_stats.final_y << "," << equiv_patch_stats.final_y << std::endl;
        pw_compare << "y_error," << pointwise_stats.y_error << "," << direct_patch_stats.y_error << "," << equiv_patch_stats.y_error << std::endl;
        pw_compare << "avg_active_count," << pointwise_stats.avg_active_count << "," << direct_patch_stats.avg_active_count << "," << equiv_patch_stats.avg_active_count << std::endl;
        pw_compare << "avg_patch_count," << pointwise_stats.avg_patch_count << "," << direct_patch_stats.avg_patch_count << "," << equiv_patch_stats.avg_patch_count << std::endl;
        pw_compare << "force_std_dev," << pointwise_stats.force_std << "," << direct_patch_stats.force_std << "," << equiv_patch_stats.force_std << std::endl;
        pw_compare << "torque_std_dev," << pointwise_stats.torque_std << "," << direct_patch_stats.torque_std << "," << equiv_patch_stats.torque_std << std::endl;
        pw_compare << "center_drift," << pointwise_stats.center_drift_total << "," << direct_patch_stats.center_drift_total << "," << equiv_patch_stats.center_drift_total << std::endl;
        pw_compare << "normal_fluctuation," << pointwise_stats.normal_fluctuation_total << "," << direct_patch_stats.normal_fluctuation_total << "," << equiv_patch_stats.normal_fluctuation_total << std::endl;
        pw_compare << "stable," << (pointwise_stats.stable ? 1 : 0) << "," << (direct_patch_stats.stable ? 1 : 0) << "," << (equiv_patch_stats.stable ? 1 : 0) << std::endl;
        pw_compare.close();
    }

    // ====================================================================
    // Summary
    // ====================================================================
    std::cout << "\n=== Final Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n--- Mode A: Pointwise ---" << std::endl;
    std::cout << "  Final Y: " << pointwise_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << pointwise_stats.y_error << " m" << std::endl;
    std::cout << "  Force std: " << pointwise_stats.force_std << " N" << std::endl;

    std::cout << "\n--- Mode B: Direct Patch Sum ---" << std::endl;
    std::cout << "  Final Y: " << direct_patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << direct_patch_stats.y_error << " m" << std::endl;
    std::cout << "  Force std: " << direct_patch_stats.force_std << " N" << std::endl;

    std::cout << "\n--- Mode C: Equivalent Patch ---" << std::endl;
    std::cout << "  Final Y: " << equiv_patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << equiv_patch_stats.y_error << " m" << std::endl;
    std::cout << "  Force std: " << equiv_patch_stats.force_std << " N" << std::endl;

    std::cout << "\n=== Verification ===" << std::endl;
    bool b_pass = direct_patch_stats.stable && direct_patch_stats.y_error < 0.02;
    std::cout << "  Direct patch stable && y_error < 0.02: " << (b_pass ? "PASS" : "FAIL") << std::endl;

    // Mode B and A should produce identical results (same force sum)
    bool identical = std::abs(direct_patch_stats.final_y - pointwise_stats.final_y) < 1e-10;
    std::cout << "  Direct patch identical to pointwise: " << (identical ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_force_output.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_force_summary.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_force_vs_pointwise.csv" << std::endl;
    std::cout << "  (also: " << out_dir + "/sdf_pointwise_output.csv)" << std::endl;
    std::cout << "  (also: " << out_dir + "/sdf_patch_force_equiv_output.csv)" << std::endl;

    return (b_pass && identical) ? 0 : 1;
}
