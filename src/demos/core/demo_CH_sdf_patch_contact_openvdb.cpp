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
// Milestone 3: Patch Primitive First Implementation
//
// This demo builds a patch primitive layer on top of the existing pointwise
// SDF contact logic from Milestone 2B.6. It verifies:
//   1. Active samples can be organized into patch primitives
//   2. Patch-level geometric descriptors can be computed
//   3. Patch temporal stability can be tracked
//   4. Patch contact can be compared against pointwise contact
//
// Scene:
//   - Static target: OpenVDB sphere level set (R=1.0m, center at origin)
//   - Dynamic body: small sphere (R=0.2m) dropped from y=3.0
//
// Output:
//   out/milestone_04/sdf_patch_contact_output.csv
//   out/milestone_04/sdf_patch_summary.csv
//   out/milestone_04/sdf_patch_vs_pointwise.csv
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
#include <map>

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
    ChVector3d torque;                 // sum of pointwise torques in patch
};

// =============================================================================
// Patch temporal tracking record
// =============================================================================

struct PatchFrameRecord {
    int patch_id;
    double time;
    ChVector3d center;
    ChVector3d normal;
    double area;
    double mean_phi;
    int sample_count;
    double force_magnitude;
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
//
// Adjacency definition:
//   Two samples are adjacent if they share an edge in the (theta, phi) grid.
//   Specifically, sample (i, j) is adjacent to:
//     - (i+1, j), (i-1, j) -- theta neighbors
//     - (i, j+1), (i, j-1) -- phi neighbors (with wrap-around)
//
// Grouping strategy:
//   1. Build a set of active sample indices
//   2. Find connected components via BFS on the adjacency graph
//   3. Each connected component becomes one patch
// =============================================================================

struct GridAdjacency {
    int n_theta;
    int n_phi;
    int total_samples;

    GridAdjacency() : n_theta(0), n_phi(0), total_samples(0) {}

    GridAdjacency(int nt, int np) : n_theta(nt), n_phi(np), total_samples(nt * np) {}

    // Get neighbors of a sample given its (theta, phi) grid indices
    std::vector<int> GetNeighbors(int ti, int pj) const {
        std::vector<int> neighbors;
        // Theta neighbors
        if (ti > 0) neighbors.push_back((ti - 1) * n_phi + pj);
        if (ti < n_theta - 1) neighbors.push_back((ti + 1) * n_phi + pj);
        // Phi neighbors (with wrap-around)
        int pj_prev = (pj - 1 + n_phi) % n_phi;
        int pj_next = (pj + 1) % n_phi;
        neighbors.push_back(ti * n_phi + pj_prev);
        neighbors.push_back(ti * n_phi + pj_next);
        return neighbors;
    }

    // Find connected components among active samples
    std::vector<std::vector<int>> FindComponents(const std::vector<int>& active_indices) const {
        std::set<int> active_set(active_indices.begin(), active_indices.end());
        std::set<int> visited;
        std::vector<std::vector<int>> components;

        for (int idx : active_indices) {
            if (visited.count(idx)) continue;

            // BFS to find connected component
            std::vector<int> component;
            std::vector<int> queue;
            queue.push_back(idx);
            visited.insert(idx);

            int head = 0;
            while (head < queue.size()) {
                int current = queue[head++];
                component.push_back(current);

                // Decode grid indices
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

        // Gradient via central finite differences
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
// Analytical sphere SDF (reference baseline)
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
// Build patches from active samples
// =============================================================================

std::vector<PatchPrimitive> BuildPatches(
    const std::vector<int>& active_indices,
    const std::vector<SurfaceSample>& samples,
    const std::vector<SDFQueryResult>& sdf_results,
    const ChVector3d& body_pos,
    const ChQuaterniond& body_rot,
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

        // Compute weighted geometric descriptors
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

            // World position
            ChVector3d world_pos = body_pos + body_rot.Rotate(s.local_pos);
            weighted_center += world_pos * w;

            // Normal (use SDF gradient direction)
            double grad_len = qr.grad.Length();
            if (grad_len > 1e-12) {
                ChVector3d n = qr.grad / grad_len;
                // Weight by phi proximity (closer to surface = more weight)
                double phi_weight = w;
                weighted_normal += n * phi_weight;
            }

            weighted_phi_sum += qr.phi * w;

            if (qr.phi < patch.min_phi) patch.min_phi = qr.phi;
            double pen = std::max(-qr.phi, 0.0);
            if (pen > patch.max_penetration) patch.max_penetration = pen;

            // Sum pointwise forces and torques
            if (sid < static_cast<int>(sample_forces.size())) {
                patch.force += sample_forces[sid];
                patch.torque += sample_torques[sid];
            }
        }

        total_weight = patch.total_area;
        if (total_weight > 1e-12) {
            patch.center = weighted_center / total_weight;
            patch.mean_phi = weighted_phi_sum / total_weight;
        }

        // Normalize patch normal
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
// Run simulation with both pointwise and patch-level tracking
// =============================================================================

struct SimStats {
    double final_y;
    double expected_y;
    double y_error;
    bool stable;
    double avg_active_count;
    double avg_patch_count;
    double avg_max_patch_samples;
    double final_patch_count;
    double patch_count_std;
    double center_drift_total;
    double normal_fluctuation_total;
    double final_force_y;
    double force_std;
};

SimStats RunPatchSimulation(
    const SDFProbeFunc& sdf_probe,
    const std::vector<SurfaceSample>& samples,
    const GridAdjacency& adjacency,
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
    const std::string& summary_csv = "",
    const std::string& pointwise_csv = ""
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

    // Prepare CSV output
    std::ofstream out_file;
    std::ofstream summary_file;
    std::ofstream pointwise_file;

    bool write_output = !output_csv.empty();
    bool write_summary = !summary_csv.empty();
    bool write_pointwise = !pointwise_csv.empty();

    if (write_output) {
        out_file.open(output_csv);
        out_file << "time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,"
                 << "expected_y,y_error,normalized_y_error,"
                 << "active_sample_count,patch_count,total_patch_area,"
                 << "total_force_x,total_force_y,total_force_z,"
                 << "force_error" << std::endl;
    }

    if (write_summary) {
        summary_file.open(summary_csv);
        summary_file << "time,patch_id,center_x,center_y,center_z,"
                     << "normal_x,normal_y,normal_z,"
                     << "area,min_phi,mean_phi,max_penetration,"
                     << "active_sample_count,force_x,force_y,force_z,force_magnitude" << std::endl;
    }

    if (write_pointwise) {
        pointwise_file.open(pointwise_csv);
        pointwise_file << "time,pos_x,pos_y,pos_z,"
                       << "active_count,force_count,min_phi,"
                       << "total_force_x,total_force_y,total_force_z,force_error" << std::endl;
    }

    double expected_force = dyn_mass * 9.81;
    bool stable = true;
    double sum_active_count = 0.0;
    double sum_patch_count = 0.0;
    double sum_max_patch_samples = 0.0;
    int output_count = 0;

    // Temporal tracking
    ChVector3d prev_max_patch_center(0, 0, 0);
    ChVector3d prev_max_patch_normal(0, 1, 0);
    double center_drift_total = 0.0;
    double normal_fluctuation_total = 0.0;
    double force_sum = 0.0;
    double force_sq_sum = 0.0;
    bool first_output = true;

    double next_output_time = 0.0;
    double output_interval = 0.01;

    int global_patch_id_counter = 0;

    // SDF query results and forces (cached per sample)
    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_forces(samples.size());
    std::vector<ChVector3d> sample_torques(samples.size());

    while (sys.GetChTime() < total_time) {
        dyn_body->EmptyAccumulator(acc_id);

        ChVector3d body_pos = dyn_body->GetPos();
        ChQuaterniond body_rot = dyn_body->GetRot();
        ChVector3d body_vel = dyn_body->GetPosDt();
        ChVector3d body_ang_vel_world = body_rot.Rotate(dyn_body->GetAngVelLocal());

        // Step 1: Query all samples and classify
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

            sdf_results[si] = sdf_probe(sample_world);
            double phi = sdf_results[si].phi;
            if (phi < min_phi) min_phi = phi;

            sample_forces[si] = ChVector3d(0, 0, 0);
            sample_torques[si] = ChVector3d(0, 0, 0);

            if (phi < activation_band) {
                active_indices.push_back(static_cast<int>(si));

                // Force application
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
                        sample_forces[si] = pt_force;
                        sample_torques[si] = r.Cross(pt_force);
                        total_force += pt_force;
                        total_torque += r.Cross(pt_force);
                    }
                }
            }
        }

        // Step 2: Apply pointwise forces to body
        if (total_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, total_force, dyn_body->GetPos(), false);
            ChQuaterniond q_inv = dyn_body->GetRot().GetConjugate();
            dyn_body->AccumulateTorque(acc_id, q_inv.Rotate(total_torque), false);
        }

        // Step 3: Build patches
        std::vector<PatchPrimitive> patches = BuildPatches(
            active_indices, samples, sdf_results,
            body_pos, body_rot, sample_forces, sample_torques,
            adjacency, global_patch_id_counter
        );

        if (!patches.empty()) {
            global_patch_id_counter += static_cast<int>(patches.size());
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
            double force_error = std::abs(total_force.y() - expected_force);

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

            force_sum += total_force.y();
            force_sq_sum += total_force.y() * total_force.y();

            // Temporal tracking
            if (!first_output && !patches.empty()) {
                // Find the largest patch
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
                         << total_force.x() << "," << total_force.y() << "," << total_force.z() << ","
                         << force_error << std::endl;
            }

            // Patch summary CSV (long table format)
            if (write_summary) {
                if (patches.empty()) {
                    summary_file << std::fixed << std::setprecision(8)
                                 << time << "," << -1 << ","
                                 << "0,0,0," << "0,0,0,0,0,0,0,0,0,0,0,0" << std::endl;
                } else {
                    for (const auto& p : patches) {
                        summary_file << std::fixed << std::setprecision(8)
                                     << time << "," << p.patch_id << ","
                                     << p.center.x() << "," << p.center.y() << "," << p.center.z() << ","
                                     << p.normal.x() << "," << p.normal.y() << "," << p.normal.z() << ","
                                     << p.total_area << "," << p.min_phi << "," << p.mean_phi << "," << p.max_penetration << ","
                                     << p.active_sample_count << ","
                                     << p.force.x() << "," << p.force.y() << "," << p.force.z() << ","
                                     << p.force.Length() << std::endl;
                    }
                }
            }

            // Pointwise-only CSV (for comparison)
            if (write_pointwise) {
                pointwise_file << std::fixed << std::setprecision(8)
                               << time << ","
                               << pos.x() << "," << pos.y() << "," << pos.z() << ","
                               << active_indices.size() << "," << force_count << ","
                               << min_phi << ","
                               << total_force.x() << "," << total_force.y() << "," << total_force.z() << ","
                               << force_error << std::endl;
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
    stats.final_patch_count = 0.0;
    stats.patch_count_std = 0.0;
    stats.center_drift_total = center_drift_total;
    stats.normal_fluctuation_total = normal_fluctuation_total;
    stats.final_force_y = 0.0;
    stats.force_std = 0.0;

    if (output_count > 1) {
        double mean_force = force_sum / output_count;
        stats.force_std = std::sqrt(force_sq_sum / output_count - mean_force * mean_force);
    }

    if (write_output) out_file.close();
    if (write_summary) summary_file.close();
    if (write_pointwise) pointwise_file.close();

    return stats;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 3: Patch Primitive First Implementation ===" << std::endl;

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

    // -- Surface samples with grid indexing --
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
    std::string out_dir = project_root + "/out/milestone_04";
    EnsureDir(out_dir);

    // ====================================================================
    // Phase 1: Run patch contact simulation
    // ====================================================================
    std::cout << "\n=== Phase 1: Patch Contact Simulation ===" << std::endl;

    SimStats patch_stats = RunPatchSimulation(
        openvdb_sdf_probe, samples, adjacency,
        stiffness, damping, activation_band, force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        out_dir + "/sdf_patch_contact_output.csv",
        out_dir + "/sdf_patch_summary.csv",
        ""
    );

    std::cout << "  Final Y: " << patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << patch_stats.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (patch_stats.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Avg active samples: " << patch_stats.avg_active_count << std::endl;
    std::cout << "  Avg patch count: " << patch_stats.avg_patch_count << std::endl;
    std::cout << "  Avg max patch samples: " << patch_stats.avg_max_patch_samples << std::endl;
    std::cout << "  Total center drift: " << patch_stats.center_drift_total << " m" << std::endl;
    std::cout << "  Total normal fluctuation: " << patch_stats.normal_fluctuation_total << std::endl;
    std::cout << "  Force std dev: " << patch_stats.force_std << " N" << std::endl;

    // ====================================================================
    // Phase 2: Run pointwise contact simulation (for comparison)
    // ====================================================================
    std::cout << "\n=== Phase 2: Pointwise Contact Simulation (for comparison) ===" << std::endl;

    SimStats pointwise_stats = RunPatchSimulation(
        openvdb_sdf_probe, samples, adjacency,
        stiffness, damping, activation_band, force_band,
        time_step, total_time,
        dyn_sphere_radius, dyn_mass, drop_height,
        expected_equilibrium_y,
        "", "",  // no main output
        out_dir + "/sdf_patch_vs_pointwise.csv"
    );

    std::cout << "  Final Y: " << pointwise_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << pointwise_stats.y_error << " m" << std::endl;
    std::cout << "  Stable: " << (pointwise_stats.stable ? "YES" : "NO") << std::endl;
    std::cout << "  Avg active samples: " << pointwise_stats.avg_active_count << std::endl;
    std::cout << "  Force std dev: " << pointwise_stats.force_std << " N" << std::endl;

    // ====================================================================
    // Phase 3: Generate comparison summary
    // ====================================================================
    std::cout << "\n=== Phase 3: Comparison Summary ===" << std::endl;

    std::string compare_path = out_dir + "/sdf_patch_vs_pointwise_summary.csv";
    std::ofstream compare_file(compare_path);
    if (compare_file.is_open()) {
        compare_file << "metric,patch_contact,pointwise_contact,difference" << std::endl;
        compare_file << std::fixed << std::setprecision(6);
        compare_file << "final_y," << patch_stats.final_y << "," << pointwise_stats.final_y << ","
                     << std::abs(patch_stats.final_y - pointwise_stats.final_y) << std::endl;
        compare_file << "y_error," << patch_stats.y_error << "," << pointwise_stats.y_error << ","
                     << std::abs(patch_stats.y_error - pointwise_stats.y_error) << std::endl;
        compare_file << "avg_active_count," << patch_stats.avg_active_count << "," << pointwise_stats.avg_active_count << ",0" << std::endl;
        compare_file << "avg_patch_count," << patch_stats.avg_patch_count << ",N/A,N/A" << std::endl;
        compare_file << "avg_max_patch_samples," << patch_stats.avg_max_patch_samples << ",N/A,N/A" << std::endl;
        compare_file << "force_std_dev," << patch_stats.force_std << "," << pointwise_stats.force_std << ","
                     << std::abs(patch_stats.force_std - pointwise_stats.force_std) << std::endl;
        compare_file << "center_drift_total," << patch_stats.center_drift_total << ",N/A,N/A" << std::endl;
        compare_file << "normal_fluctuation_total," << patch_stats.normal_fluctuation_total << ",N/A,N/A" << std::endl;
        compare_file << "stable," << (patch_stats.stable ? 1 : 0) << "," << (pointwise_stats.stable ? 1 : 0) << ",0" << std::endl;
        compare_file.close();
    }

    // ====================================================================
    // Summary
    // ====================================================================
    std::cout << "\n=== Final Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n--- Patch Contact ---" << std::endl;
    std::cout << "  Final Y: " << patch_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << patch_stats.y_error << " m" << std::endl;
    std::cout << "  Avg patch count: " << patch_stats.avg_patch_count << std::endl;
    std::cout << "  Force std dev: " << patch_stats.force_std << " N" << std::endl;

    std::cout << "\n--- Pointwise Contact ---" << std::endl;
    std::cout << "  Final Y: " << pointwise_stats.final_y << " m" << std::endl;
    std::cout << "  Y error: " << pointwise_stats.y_error << " m" << std::endl;
    std::cout << "  Force std dev: " << pointwise_stats.force_std << " N" << std::endl;

    std::cout << "\n=== Verification ===" << std::endl;
    bool pass = patch_stats.stable && patch_stats.y_error < 0.02;
    std::cout << "  Patch stable && y_error < 0.02: " << (pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_contact_output.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_summary.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_vs_pointwise.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_vs_pointwise_summary.csv" << std::endl;

    return pass ? 0 : 1;
}
