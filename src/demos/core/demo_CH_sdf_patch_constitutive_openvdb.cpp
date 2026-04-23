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
// Milestone 6: Patch Constitutive Law First Version
//
// This demo upgrades patch from "force aggregation container" to "direct force generator".
// It verifies:
//   1. Patch can generate contact force directly from patch geometric descriptors
//   2. Patch constitutive law replaces pointwise sample-by-sample force computation
//   3. Patch constitutive law maintains numerical stability
//   4. Patch constitutive law achieves reasonable equilibrium position and contact force
//   5. Comparison with pointwise and patch-force-sum baselines
//
// Scene: OpenVDB sphere (R=1.0) + dynamic sphere (R=0.2) dropped from y=3.0
//
// Constitutive Law Definition:
//   penetration_eff = area-weighted mean of max(-phi_i, 0)
//   vn_eff = patch center velocity projected to patch normal
//   F_patch = (k_patch * penetration_eff + c_patch * max(-vn_eff, 0)) * patch_normal
//   T_patch = (patch_center - COM) x F_patch
//
// Output:
//   out/milestone_07/sdf_patch_constitutive_output.csv     (frame-level output)
//   out/milestone_07/sdf_patch_constitutive_summary.csv    (summary metrics)
//   out/milestone_07/sdf_patch_constitutive_vs_baselines.csv (baseline comparison)
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
#include <sstream>

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
// Surface sample point
// =============================================================================

struct SurfaceSample {
    ChVector3d local_pos;
    double area_weight;
    int theta_index;
    int phi_index;
    int global_index;
};

// =============================================================================
// Patch primitive (frame-local)
// =============================================================================

struct PatchPrimitive {
    int patch_id;                      // frame-local ID
    std::vector<int> sample_ids;
    int active_sample_count;
    ChVector3d center;
    ChVector3d normal;
    double total_area;
    double min_phi;
    double mean_phi;
    double max_penetration;
    double effective_penetration;     // for constitutive law
    double effective_normal_velocity; // for constitutive law
    ChVector3d force;
    ChVector3d torque;
};

// =============================================================================
// Track status
// =============================================================================

enum class TrackStatus {
    Born,   // first appearance
    Alive,  // continuing from previous frame
    Dead    // no longer present (not written to CSV, used in summary)
};

// =============================================================================
// Persistent patch track (carries state across frames)
// =============================================================================

struct PersistentPatchTrack {
    int persistent_id;
    int current_frame_patch_id;      // frame-local patch ID in current frame
    TrackStatus status;
    ChVector3d center;
    ChVector3d normal;
    double area;
    double mean_phi;
    double min_phi;
    double max_penetration;
    double effective_penetration;
    double effective_normal_velocity;
    ChVector3d force;
    ChVector3d torque;
    int age_in_steps;
    double birth_time;
    double last_seen_time;
    int birth_frame;
};

// =============================================================================
// Tracking summary metrics
// =============================================================================

struct TrackingMetrics {
    int total_born;
    int total_dead;
    int total_frames_with_patches;
    double total_frames;
    double avg_patch_count;
    double max_patch_count;
    double avg_matched_ratio;
    int max_lifetime_steps;
    int avg_lifetime_steps;
    double max_center_drift;
    double avg_center_drift;
    double max_normal_fluctuation;
};

// =============================================================================
// Contact mode enumeration
// =============================================================================

enum class ContactMode {
    Pointwise,        // sample-by-sample force (baseline)
    PatchForceSum,    // patch_force = sum of sample forces
    PatchConstitutive // patch force from constitutive law
};

// =============================================================================
// Mode comparison metrics
// =============================================================================

struct ModeMetrics {
    ContactMode mode;
    std::string mode_name;
    double final_y;
    double y_error;
    double expected_equilibrium_y;
    double total_contact_force_y;
    double total_contact_torque_magnitude;
    double force_std_dev;
    double avg_force_y;
    int active_sample_count;
    int patch_count;
    bool is_stable;
};

// =============================================================================
// Generate surface samples
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
// Grid adjacency and connected components
// =============================================================================

struct GridAdjacency {
    int n_theta;
    int n_phi;

    GridAdjacency() : n_theta(0), n_phi(0) {}
    GridAdjacency(int nt, int np) : n_theta(nt), n_phi(np) {}

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
// OpenVDB sphere SDF
// =============================================================================

struct OpenVDBSphereSDF {
    openvdb::FloatGrid::Ptr grid;
    ChVector3d center_world;
    double sphere_radius;
    double voxel_size;

    OpenVDBSphereSDF() : sphere_radius(0.0), voxel_size(0.0) {}

    SDFQueryResult Query(const ChVector3d& world_pt) const {
        SDFQueryResult qr;
        ChVector3d local = world_pt - center_world;

        openvdb::tools::GridSampler<
            openvdb::FloatGrid::TreeType,
            openvdb::tools::BoxSampler
        > sampler(grid->tree(), grid->transform());
        openvdb::Vec3d lv(local.x(), local.y(), local.z());
        qr.phi = static_cast<double>(sampler.wsSample(lv));

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
// Helpers
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
    const std::vector<ChVector3d>& sample_world_positions,
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
        double weighted_penetration_sum = 0.0;
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
            weighted_penetration_sum += pen * w;
        }

        total_weight = patch.total_area;
        if (total_weight > 1e-12) {
            patch.center = weighted_center / total_weight;
            patch.mean_phi = weighted_phi_sum / total_weight;
            patch.effective_penetration = weighted_penetration_sum / total_weight;
        }

        double norm_len = weighted_normal.Length();
        if (norm_len > 1e-12) {
            patch.normal = weighted_normal / norm_len;
        } else {
            patch.normal = ChVector3d(0, 1, 0);
        }

        patch.effective_normal_velocity = 0.0;

        patches.push_back(patch);
    }

    return patches;
}

// =============================================================================
// Patch matching: current patches to previous persistent tracks
// =============================================================================

struct MatchResult {
    int current_patch_id;
    int persistent_id;
    bool matched;
};

std::vector<MatchResult> MatchPatches(
    const std::vector<PatchPrimitive>& current_patches,
    const std::vector<PersistentPatchTrack>& prev_tracks,
    double distance_threshold,
    double normal_similarity_threshold
) {
    std::vector<MatchResult> matches;
    std::set<int> used_persistent_ids;

    std::vector<size_t> indices(current_patches.size());
    for (size_t i = 0; i < indices.size(); i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(),
              [&current_patches](size_t a, size_t b) {
                  return current_patches[a].force.Length() > current_patches[b].force.Length();
              });

    for (size_t idx : indices) {
        const auto& cp = current_patches[idx];
        MatchResult m;
        m.current_patch_id = cp.patch_id;
        m.persistent_id = -1;
        m.matched = false;

        double best_dist = distance_threshold + 1.0;
        int best_persistent_id = -1;

        for (const auto& track : prev_tracks) {
            double dist = (cp.center - track.center).Length();
            double cos_sim = cp.normal.Dot(track.normal);

            if (cos_sim < normal_similarity_threshold) continue;
            if (dist < best_dist) {
                best_dist = dist;
                best_persistent_id = track.persistent_id;
            }
        }

        if (best_persistent_id >= 0 && !used_persistent_ids.count(best_persistent_id)) {
            m.persistent_id = best_persistent_id;
            m.matched = true;
            used_persistent_ids.insert(best_persistent_id);
        }

        matches.push_back(m);
    }

    std::sort(matches.begin(), matches.end(),
              [](const MatchResult& a, const MatchResult& b) {
                  return a.current_patch_id < b.current_patch_id;
              });

    return matches;
}

// =============================================================================
// Compute patch constitutive force
// =============================================================================

ChVector3d ComputePatchConstitutiveForce(
    const PatchPrimitive& patch,
    const ChVector3d& patch_center_velocity,
    double k_patch,
    double c_patch
) {
    // penetration_eff = area-weighted mean(max(-phi_i, 0))
    double penetration_eff = patch.effective_penetration;

    // vn_eff = patch center velocity projected to patch normal
    double vn_eff = patch_center_velocity.Dot(patch.normal);

    // F_patch = (k_patch * penetration_eff + c_patch * max(-vn_eff, 0)) * patch_normal
    double force_scalar = k_patch * penetration_eff + c_patch * std::max(-vn_eff, 0.0);

    if (force_scalar < 0.0) {
        force_scalar = 0.0;
    }

    return patch.normal * force_scalar;
}

// =============================================================================
// Run single mode simulation
// =============================================================================

struct SimulationResult {
    double final_y;
    double y_error;
    double expected_equilibrium_y;
    double avg_force_y;
    double force_std_dev;
    double total_torque_magnitude;
    int max_patch_count;
    int avg_active_samples;
    bool is_stable;
    std::vector<double> force_y_history;
    std::vector<double> pos_y_history;
    std::vector<int> patch_count_history;
};

SimulationResult RunSimulationMode(
    ContactMode mode,
    const std::string& mode_name,
    double time_step,
    double total_time,
    double static_sphere_radius,
    double dyn_sphere_radius,
    double dyn_mass,
    double drop_height,
    const OpenVDBSphereSDF& openvdb_sdf,
    const SDFProbeFunc& sdf_probe,
    int n_theta,
    int n_phi,
    double stiffness,
    double damping,
    double activation_band,
    double force_band,
    double match_distance_threshold,
    double normal_similarity_threshold,
    const std::string& out_dir
) {
    // -- Chrono setup --
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -9.81, 0));

    auto dyn_body = chrono_types::make_shared<ChBodyEasySphere>(dyn_sphere_radius, 1000.0, false, false);
    dyn_body->SetPos(ChVector3d(0.0, drop_height, 0.0));
    dyn_body->SetFixed(false);
    sys.AddBody(dyn_body);

    unsigned int acc_id = dyn_body->AddAccumulator();
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(100);
    sys.GetSolver()->AsIterative()->SetTolerance(1e-6);

    // -- Persistent state --
    std::vector<PersistentPatchTrack> active_tracks;
    int next_persistent_id = 0;

    // -- Buffers --
    std::vector<SurfaceSample> samples = GenerateSphereSamples(dyn_sphere_radius, n_theta, n_phi);
    GridAdjacency adjacency(n_theta, n_phi);
    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());
    std::vector<ChVector3d> sample_forces(samples.size());
    std::vector<ChVector3d> sample_torques(samples.size());

    // -- Metrics --
    std::vector<double> force_y_values;
    std::vector<double> pos_y_values;
    std::vector<int> patch_counts;
    double total_force_y = 0.0;
    int force_sample_count = 0;

    // -- Expected equilibrium: mg = k * penetration -> penetration = mg/k
    double expected_equilibrium_y = static_sphere_radius + dyn_sphere_radius - (dyn_mass * 9.81) / stiffness;

    while (sys.GetChTime() < total_time) {
        dyn_body->EmptyAccumulator(acc_id);

        ChVector3d body_pos = dyn_body->GetPos();
        ChQuaterniond body_rot = dyn_body->GetRot();
        ChVector3d body_vel = dyn_body->GetPosDt();
        ChVector3d body_ang_vel_world = body_rot.Rotate(dyn_body->GetAngVelLocal());

        // Query samples
        std::vector<int> active_indices;
        ChVector3d total_force(0, 0, 0);
        ChVector3d total_torque(0, 0, 0);

        for (size_t si = 0; si < samples.size(); si++) {
            const auto& s = samples[si];
            ChVector3d sample_world = body_pos + body_rot.Rotate(s.local_pos);
            ChVector3d r = sample_world - body_pos;
            ChVector3d sample_vel = body_vel + body_ang_vel_world.Cross(r);

            sample_world_positions[si] = sample_world;
            sdf_results[si] = sdf_probe(sample_world);
            sample_forces[si] = ChVector3d(0, 0, 0);
            sample_torques[si] = ChVector3d(0, 0, 0);

            double phi = sdf_results[si].phi;
            if (phi < activation_band) {
                active_indices.push_back(static_cast<int>(si));

                // Pointwise mode: compute sample forces
                if (mode == ContactMode::Pointwise) {
                    if (phi < force_band) {
                        double pen = std::max(-phi, 0.0);
                        double grad_len = sdf_results[si].grad.Length();
                        if (grad_len > 1e-12) {
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
            }
        }

        // Patch-based modes
        if (mode == ContactMode::PatchForceSum || mode == ContactMode::PatchConstitutive) {
            int frame_patch_id_offset = next_persistent_id * 1000;
            std::vector<PatchPrimitive> patches = BuildPatches(
                active_indices, samples, sdf_results,
                sample_world_positions, adjacency, frame_patch_id_offset
            );

            // Patch matching
            std::vector<MatchResult> matches = MatchPatches(
                patches, active_tracks,
                match_distance_threshold, normal_similarity_threshold
            );

            std::vector<PersistentPatchTrack> current_tracks;
            std::set<int> matched_persistent_ids;

            for (const auto& m : matches) {
                PersistentPatchTrack track;
                if (m.matched) {
                    for (const auto& prev : active_tracks) {
                        if (prev.persistent_id == m.persistent_id) {
                            track = prev;
                            track.status = TrackStatus::Alive;
                            break;
                        }
                    }
                    track.current_frame_patch_id = m.current_patch_id;
                    matched_persistent_ids.insert(m.persistent_id);
                } else {
                    track.persistent_id = next_persistent_id++;
                    track.current_frame_patch_id = m.current_patch_id;
                    track.status = TrackStatus::Born;
                    track.birth_time = sys.GetChTime();
                    track.birth_frame = 0;
                    track.age_in_steps = 0;
                }

                // Update track data from current patch
                for (const auto& p : patches) {
                    if (p.patch_id == track.current_frame_patch_id) {
                        track.center = p.center;
                        track.normal = p.normal;
                        track.area = p.total_area;
                        track.mean_phi = p.mean_phi;
                        track.min_phi = p.min_phi;
                        track.max_penetration = p.max_penetration;
                        track.effective_penetration = p.effective_penetration;
                        break;
                    }
                }

                track.last_seen_time = sys.GetChTime();
                if (track.status == TrackStatus::Alive) {
                    track.age_in_steps++;
                }

                current_tracks.push_back(track);
            }

            // Compute patch forces
            for (size_t pi = 0; pi < patches.size(); pi++) {
                auto& patch = patches[pi];

                // Compute patch center velocity
                ChVector3d patch_center_vel = body_vel + body_ang_vel_world.Cross(patch.center - body_pos);
                patch.effective_normal_velocity = patch_center_vel.Dot(patch.normal);

                if (mode == ContactMode::PatchForceSum) {
                    // Sum sample forces within patch
                    for (int sid : patch.sample_ids) {
                        double phi = sdf_results[sid].phi;
                        if (phi < force_band) {
                            double pen = std::max(-phi, 0.0);
                            double grad_len = sdf_results[sid].grad.Length();
                            if (grad_len > 1e-12) {
                                ChVector3d normal = sdf_results[sid].grad / grad_len;
                                ChVector3d sample_world = sample_world_positions[sid];
                                ChVector3d r = sample_world - body_pos;
                                ChVector3d sample_vel = body_vel + body_ang_vel_world.Cross(r);
                                double vn = sample_vel.Dot(normal);
                                double force_mag = stiffness * pen + damping * std::max(-vn, 0.0);
                                if (force_mag > 0.0) {
                                    ChVector3d pt_force = normal * force_mag;
                                    patch.force += pt_force;
                                    patch.torque += r.Cross(pt_force);
                                }
                            }
                        }
                    }
                } else if (mode == ContactMode::PatchConstitutive) {
                    // Constitutive law: F_patch = (k * penetration_eff + c * max(-vn_eff, 0)) * normal
                    patch.force = ComputePatchConstitutiveForce(
                        patch, patch_center_vel, stiffness, damping
                    );
                    patch.torque = (patch.center - body_pos).Cross(patch.force);
                }

                total_force += patch.force;
                total_torque += patch.torque;
            }

            active_tracks = current_tracks;
        }

        // Apply forces
        if (total_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, total_force, body_pos, false);
            dyn_body->AccumulateTorque(acc_id, body_rot.GetConjugate().Rotate(total_torque), false);
        }

        sys.DoStepDynamics(time_step);

        // Record metrics
        double current_force_y = total_force.y();
        force_y_values.push_back(current_force_y);
        pos_y_values.push_back(dyn_body->GetPos().y());
        patch_counts.push_back(static_cast<int>(active_indices.size()));

        total_force_y += current_force_y;
        force_sample_count++;
    }

    // Compute statistics
    SimulationResult result;
    result.final_y = pos_y_values.empty() ? 0.0 : pos_y_values.back();
    result.y_error = std::abs(result.final_y - expected_equilibrium_y);
    result.expected_equilibrium_y = expected_equilibrium_y;
    result.avg_force_y = force_sample_count > 0 ? total_force_y / force_sample_count : 0.0;

    // Force std dev
    double sum_sq = 0.0;
    for (double fy : force_y_values) {
        sum_sq += (fy - result.avg_force_y) * (fy - result.avg_force_y);
    }
    result.force_std_dev = force_sample_count > 0 ? std::sqrt(sum_sq / force_sample_count) : 0.0;

    result.total_torque_magnitude = 0.0;
    result.max_patch_count = patch_counts.empty() ? 0 : *std::max_element(patch_counts.begin(), patch_counts.end());
    result.avg_active_samples = patch_counts.empty() ? 0 :
        static_cast<int>(std::accumulate(patch_counts.begin(), patch_counts.end(), 0) / patch_counts.size());

    result.force_y_history = force_y_values;
    result.pos_y_history = pos_y_values;
    result.patch_count_history = patch_counts;

    // Stability: check if position is converging (last 10% of simulation)
    size_t last_10_percent = pos_y_values.size() * 9 / 10;
    double avg_last = 0.0;
    for (size_t i = last_10_percent; i < pos_y_values.size(); i++) {
        avg_last += pos_y_values[i];
    }
    avg_last /= (pos_y_values.size() - last_10_percent);
    result.is_stable = std::abs(avg_last - expected_equilibrium_y) < 0.05;

    return result;
}

// =============================================================================
// Main simulation
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 6: Patch Constitutive Law First Version ===" << std::endl;

    openvdb::initialize();

    // -- Parameters --
    double time_step = 5e-4;
    double total_time = 2.0;

    double static_sphere_radius = 1.0;
    double dyn_sphere_radius = 0.2;
    double dyn_mass = (4.0 / 3.0) * M_PI * std::pow(dyn_sphere_radius, 3) * 1000.0;
    double drop_height = 3.0;

    double voxel_size = 0.05;
    auto sphere_grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        static_cast<float>(static_sphere_radius),
        openvdb::Vec3f(0.0f, 0.0f, 0.0f),
        static_cast<float>(voxel_size),
        3.0f
    );
    sphere_grid->setGridClass(openvdb::GRID_LEVEL_SET);

    OpenVDBSphereSDF openvdb_sdf;
    openvdb_sdf.grid = sphere_grid;
    openvdb_sdf.center_world = ChVector3d(0, 0, 0);
    openvdb_sdf.sphere_radius = static_sphere_radius;
    openvdb_sdf.voxel_size = voxel_size;

    SDFProbeFunc sdf_probe = [&](const ChVector3d& pt) -> SDFQueryResult {
        return openvdb_sdf.Query(pt);
    };

    int n_theta = 8, n_phi = 16;

    double stiffness = 1e5;
    double damping = 500;
    double activation_band = 0.1;
    double force_band = 0.0;

    double match_distance_threshold = 0.5 * dyn_sphere_radius;
    double normal_similarity_threshold = 0.9;

    std::cout << "\nGeometry: static sphere R=" << static_sphere_radius
              << "m, dynamic sphere R=" << dyn_sphere_radius << "m, mass=" << dyn_mass << "kg" << std::endl;
    std::cout << "Contact: stiffness=" << stiffness << " N/m, damping=" << damping << " Ns/m" << std::endl;
    std::cout << "Constitutive Law: F_patch = (k * penetration_eff + c * max(-vn_eff, 0)) * patch_normal" << std::endl;
    std::cout << "  penetration_eff = area-weighted mean(max(-phi_i, 0))" << std::endl;
    std::cout << "  vn_eff = patch center velocity . patch normal" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_07";
    EnsureDir(out_dir);

    // -- Run three modes --
    std::vector<ContactMode> modes = {
        ContactMode::Pointwise,
        ContactMode::PatchForceSum,
        ContactMode::PatchConstitutive
    };
    std::vector<std::string> mode_names = {
        "Pointwise",
        "PatchForceSum",
        "PatchConstitutive"
    };

    std::vector<SimulationResult> results;

    for (size_t mi = 0; mi < modes.size(); mi++) {
        std::cout << "\n=== Running Mode: " << mode_names[mi] << " ===" << std::endl;

        SimulationResult result = RunSimulationMode(
            modes[mi],
            mode_names[mi],
            time_step,
            total_time,
            static_sphere_radius,
            dyn_sphere_radius,
            dyn_mass,
            drop_height,
            openvdb_sdf,
            sdf_probe,
            n_theta,
            n_phi,
            stiffness,
            damping,
            activation_band,
            force_band,
            match_distance_threshold,
            normal_similarity_threshold,
            out_dir
        );

        results.push_back(result);

        std::cout << "  Final Y: " << std::fixed << std::setprecision(6) << result.final_y << " m" << std::endl;
        std::cout << "  Expected Y: " << result.expected_equilibrium_y << " m" << std::endl;
        std::cout << "  Y Error: " << result.y_error << " m" << std::endl;
        std::cout << "  Avg Force Y: " << result.avg_force_y << " N" << std::endl;
        std::cout << "  Force Std Dev: " << result.force_std_dev << " N" << std::endl;
        std::cout << "  Stable: " << (result.is_stable ? "YES" : "NO") << std::endl;
    }

    // -- Write output CSV --
    std::ofstream output_file(out_dir + "/sdf_patch_constitutive_output.csv");
    output_file << "mode,final_y,expected_y,y_error,avg_force_y,force_std_dev,torque_magnitude,stable" << std::endl;
    for (size_t mi = 0; mi < modes.size(); mi++) {
        output_file << mode_names[mi] << ","
                    << std::fixed << std::setprecision(8)
                    << results[mi].final_y << ","
                    << results[mi].expected_equilibrium_y << ","
                    << results[mi].y_error << ","
                    << results[mi].avg_force_y << ","
                    << results[mi].force_std_dev << ","
                    << results[mi].total_torque_magnitude << ","
                    << (results[mi].is_stable ? 1 : 0) << std::endl;
    }
    output_file.close();

    // -- Write summary CSV --
    std::ofstream summary_file(out_dir + "/sdf_patch_constitutive_summary.csv");
    summary_file << "metric,pointwise,patch_force_sum,patch_constitutive" << std::endl;
    summary_file << "final_y,"
                 << std::fixed << std::setprecision(6)
                 << results[0].final_y << ","
                 << results[1].final_y << ","
                 << results[2].final_y << std::endl;
    summary_file << "y_error,"
                 << results[0].y_error << ","
                 << results[1].y_error << ","
                 << results[2].y_error << std::endl;
    summary_file << "avg_force_y,"
                 << results[0].avg_force_y << ","
                 << results[1].avg_force_y << ","
                 << results[2].avg_force_y << std::endl;
    summary_file << "force_std_dev,"
                 << results[0].force_std_dev << ","
                 << results[1].force_std_dev << ","
                 << results[2].force_std_dev << std::endl;
    summary_file << "stable,"
                 << (results[0].is_stable ? 1 : 0) << ","
                 << (results[1].is_stable ? 1 : 0) << ","
                 << (results[2].is_stable ? 1 : 0) << std::endl;
    summary_file.close();

    // -- Write comparison CSV --
    std::ofstream compare_file(out_dir + "/sdf_patch_constitutive_vs_baselines.csv");
    compare_file << "comparison_metric,pointwise,patch_force_sum,patch_constitutive,notes" << std::endl;
    compare_file << "final_y,"
                 << std::fixed << std::setprecision(6)
                 << results[0].final_y << ","
                 << results[1].final_y << ","
                 << results[2].final_y
                 << ",Equilibrium position" << std::endl;
    compare_file << "y_error,"
                 << results[0].y_error << ","
                 << results[1].y_error << ","
                 << results[2].y_error
                 << ",Deviation from expected equilibrium" << std::endl;
    compare_file << "avg_force_y,"
                 << results[0].avg_force_y << ","
                 << results[1].avg_force_y << ","
                 << results[2].avg_force_y
                 << ",Expected mg=" << (dyn_mass * 9.81) << " N" << std::endl;
    compare_file << "force_std_dev,"
                 << results[0].force_std_dev << ","
                 << results[1].force_std_dev << ","
                 << results[2].force_std_dev
                 << ",Force oscillation magnitude" << std::endl;
    compare_file << "stable,"
                 << (results[0].is_stable ? 1 : 0) << ","
                 << (results[1].is_stable ? 1 : 0) << ","
                 << (results[2].is_stable ? 1 : 0)
                 << ",Position convergence check" << std::endl;
    compare_file.close();

    // -- Console summary --
    std::cout << "\n=== Comparison Summary ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Metric"
              << std::setw(20) << "Pointwise"
              << std::setw(20) << "PatchForceSum"
              << std::setw(20) << "PatchConstitutive" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::left << std::setw(20) << "Final Y (m)"
              << std::fixed << std::setprecision(6)
              << std::setw(20) << results[0].final_y
              << std::setw(20) << results[1].final_y
              << std::setw(20) << results[2].final_y << std::endl;
    std::cout << std::left << std::setw(20) << "Y Error (m)"
              << std::setw(20) << results[0].y_error
              << std::setw(20) << results[1].y_error
              << std::setw(20) << results[2].y_error << std::endl;
    std::cout << std::left << std::setw(20) << "Avg Force (N)"
              << std::setw(20) << results[0].avg_force_y
              << std::setw(20) << results[1].avg_force_y
              << std::setw(20) << results[2].avg_force_y << std::endl;
    std::cout << std::left << std::setw(20) << "Force Std (N)"
              << std::setw(20) << results[0].force_std_dev
              << std::setw(20) << results[1].force_std_dev
              << std::setw(20) << results[2].force_std_dev << std::endl;
    std::cout << std::left << std::setw(20) << "Stable"
              << std::setw(20) << (results[0].is_stable ? "YES" : "NO")
              << std::setw(20) << (results[1].is_stable ? "YES" : "NO")
              << std::setw(20) << (results[2].is_stable ? "YES" : "NO") << std::endl;

    // -- Verification --
    bool all_run_success = true;
    bool constitutive_stable = results[2].is_stable;
    bool constitutive_reasonable = results[2].y_error < 0.05;

    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "  All modes completed: " << (all_run_success ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Constitutive mode stable: " << (constitutive_stable ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Constitutive mode reasonable y_error: " << (constitutive_reasonable ? "PASS" : "FAIL") << std::endl;

    bool pass = all_run_success && constitutive_stable && constitutive_reasonable;
    std::cout << "  Overall: " << (pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_constitutive_output.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_constitutive_summary.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_constitutive_vs_baselines.csv" << std::endl;

    return pass ? 0 : 1;
}
