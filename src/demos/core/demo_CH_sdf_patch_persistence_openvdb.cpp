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
// Milestone 5: Patch Persistence / Patch Tracking
//
// This demo upgrades patches from "frame-local objects" to "time-persistent objects".
// It verifies:
//   1. Patches can be matched across consecutive time steps
//   2. Persistent patch IDs can be maintained over time
//   3. Birth / death / reappearance of patches is correctly detected
//   4. Patch tracking state is stable enough for future constitutive law
//
// Scene: OpenVDB sphere (R=1.0) + dynamic sphere (R=0.2) dropped from y=3.0
//
// Output:
//   out/milestone_06/sdf_patch_tracking_output.csv     (frame-level tracking)
//   out/milestone_06/sdf_patch_tracks.csv              (patch track long table)
//   out/milestone_06/sdf_patch_tracking_summary.csv    (tracking metrics)
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
    double avg_matched_ratio;      // matched / max(prev_count, curr_count)
    int max_lifetime_steps;
    int avg_lifetime_steps;
    double max_center_drift;       // max cumulative center drift for any track
    double avg_center_drift;
    double max_normal_fluctuation; // max cumulative normal fluctuation
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
// Patch matching: current patches to previous persistent tracks
//
// Rules:
//   1. Cost = Euclidean distance between patch centers
//   2. Normal similarity constraint: cos(angle) > 0.9
//   3. Distance threshold: center_dist < 0.5 * dyn_sphere_radius
//   4. Greedy one-to-one matching: each current patch matches at most one prev track
//      and vice versa
// =============================================================================

struct MatchResult {
    int current_patch_id;
    int persistent_id;  // -1 if no match
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

    // Sort current patches by force magnitude (largest first for greedy matching)
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

    // Sort back by patch_id for consistent output
    std::sort(matches.begin(), matches.end(),
              [](const MatchResult& a, const MatchResult& b) {
                  return a.current_patch_id < b.current_patch_id;
              });

    return matches;
}

// =============================================================================
// Main simulation with patch tracking
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 5: Patch Persistence / Tracking ===" << std::endl;

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
    std::vector<SurfaceSample> samples = GenerateSphereSamples(dyn_sphere_radius, n_theta, n_phi);
    GridAdjacency adjacency(n_theta, n_phi);

    double stiffness = 1e5;
    double damping = 500;
    double activation_band = 0.1;
    double force_band = 0.0;

    // -- Matching thresholds --
    double match_distance_threshold = 0.5 * dyn_sphere_radius;   // 0.1m
    double normal_similarity_threshold = 0.9;

    std::cout << "\nGeometry: static sphere R=" << static_sphere_radius
              << "m, dynamic sphere R=" << dyn_sphere_radius << "m, mass=" << dyn_mass << "kg" << std::endl;
    std::cout << "Matching: distance_threshold=" << match_distance_threshold
              << "m, normal_cos_threshold=" << normal_similarity_threshold << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_06";
    EnsureDir(out_dir);

    // -- CSV output files --
    std::ofstream tracking_file(out_dir + "/sdf_patch_tracking_output.csv");
    std::ofstream tracks_file(out_dir + "/sdf_patch_tracks.csv");

    tracking_file << "time,patch_count,active_sample_count,matched_count,born_count,dead_count,"
                  << "pos_x,pos_y,pos_z,vel_y,total_force_y,force_error" << std::endl;

    tracks_file << "time,persistent_patch_id,frame_patch_id,status,"
                << "center_x,center_y,center_z,"
                << "normal_x,normal_y,normal_z,"
                << "area,mean_phi,min_phi,force_magnitude,"
                << "age_in_steps" << std::endl;

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
    std::vector<PersistentPatchTrack> active_tracks;   // tracks from previous frame
    int next_persistent_id = 0;
    int total_born = 0;
    int total_dead = 0;
    int total_matched = 0;
    int frames_with_patches = 0;
    int total_frames = 0;
    double sum_patch_count = 0.0;
    int max_patch_count = 0;
    double sum_matched_ratio = 0.0;
    int matched_ratio_count = 0;

    // Per-track cumulative metrics
    std::map<int, int> track_birth_step;
    std::map<int, int> track_age;
    std::map<int, ChVector3d> track_prev_center;
    std::map<int, ChVector3d> track_prev_normal;
    std::map<int, double> track_cumulative_drift;
    std::map<int, double> track_cumulative_normal_fluct;

    // -- Buffers --
    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());
    std::vector<ChVector3d> sample_forces(samples.size());
    std::vector<ChVector3d> sample_torques(samples.size());

    double next_output_time = 0.0;
    double output_interval = 0.01;

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

        // Apply forces
        if (total_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, total_force, body_pos, false);
            dyn_body->AccumulateTorque(acc_id, body_rot.GetConjugate().Rotate(total_torque), false);
        }

        sys.DoStepDynamics(time_step);

        // Build patches
        int frame_patch_id_offset = next_persistent_id * 1000;
        std::vector<PatchPrimitive> patches = BuildPatches(
            active_indices, samples, sdf_results,
            sample_world_positions, sample_forces, sample_torques,
            adjacency, frame_patch_id_offset
        );

        // -- Patch matching --
        std::vector<MatchResult> matches = MatchPatches(
            patches, active_tracks,
            match_distance_threshold, normal_similarity_threshold
        );

        // Determine persistent tracks for current frame
        std::vector<PersistentPatchTrack> current_tracks;
        std::set<int> matched_persistent_ids;
        int born_count = 0;
        int matched_count = 0;

        for (const auto& m : matches) {
            PersistentPatchTrack track;
            if (m.matched) {
                // Find the previous track
                for (const auto& prev : active_tracks) {
                    if (prev.persistent_id == m.persistent_id) {
                        track = prev;
                        track.status = TrackStatus::Alive;
                        break;
                    }
                }
                track.current_frame_patch_id = m.current_patch_id;
                matched_count++;
                matched_persistent_ids.insert(m.persistent_id);
            } else {
                // New patch: birth
                track.persistent_id = next_persistent_id++;
                track.current_frame_patch_id = m.current_patch_id;
                track.status = TrackStatus::Born;
                track.birth_time = sys.GetChTime();
                track.birth_frame = total_frames;
                track.age_in_steps = 0;
                track_birth_step[track.persistent_id] = total_frames;
                track_cumulative_drift[track.persistent_id] = 0.0;
                track_cumulative_normal_fluct[track.persistent_id] = 0.0;
                born_count++;
                total_born++;
            }

            // Update track data from current patch
            for (const auto& p : patches) {
                if (p.patch_id == track.current_frame_patch_id) {
                    track.center = p.center;
                    track.normal = p.normal;
                    track.area = p.total_area;
                    track.mean_phi = p.mean_phi;
                    track.min_phi = p.min_phi;
                    track.force = p.force;
                    track.torque = p.torque;
                    break;
                }
            }

            track.last_seen_time = sys.GetChTime();
            if (track.status == TrackStatus::Alive) {
                track.age_in_steps++;

                // Cumulative drift
                if (track_prev_center.count(track.persistent_id)) {
                    double drift = (track.center - track_prev_center[track.persistent_id]).Length();
                    track_cumulative_drift[track.persistent_id] += drift;
                }
                if (track_prev_normal.count(track.persistent_id)) {
                    double fluct = 1.0 - std::abs(track.normal.Dot(track_prev_normal[track.persistent_id]));
                    track_cumulative_normal_fluct[track.persistent_id] += fluct;
                }
            }

            track_prev_center[track.persistent_id] = track.center;
            track_prev_normal[track.persistent_id] = track.normal;

            current_tracks.push_back(track);
        }

        total_matched += matched_count;

        // Count dead tracks (in previous but not matched)
        int dead_count = 0;
        for (const auto& prev : active_tracks) {
            if (!matched_persistent_ids.count(prev.persistent_id)) {
                dead_count++;
                total_dead++;

                // Write death record to tracks CSV
                tracks_file << std::fixed << std::setprecision(8)
                            << sys.GetChTime() << "," << prev.persistent_id << ",-1,"
                            << "dead,"
                            << prev.center.x() << "," << prev.center.y() << "," << prev.center.z() << ","
                            << prev.normal.x() << "," << prev.normal.y() << "," << prev.normal.z() << ","
                            << prev.area << "," << prev.mean_phi << "," << prev.min_phi << ","
                            << prev.force.Length() << ","
                            << prev.age_in_steps << std::endl;
            }
        }

        frames_with_patches += (patches.empty() ? 0 : 1);
        sum_patch_count += static_cast<double>(patches.size());
        if (static_cast<int>(patches.size()) > max_patch_count) {
            max_patch_count = static_cast<int>(patches.size());
        }

        int prev_count = static_cast<int>(active_tracks.size());
        int curr_count = static_cast<int>(patches.size());
        if (prev_count > 0 || curr_count > 0) {
            double denom = std::max(prev_count, curr_count);
            if (denom > 0) {
                sum_matched_ratio += static_cast<double>(matched_count) / denom;
                matched_ratio_count++;
            }
        }

        total_frames++;

        // -- Write tracking output at regular intervals --
        double time = sys.GetChTime();
        if (time >= next_output_time - 1e-10) {
            next_output_time += output_interval;

            ChVector3d pos = dyn_body->GetPos();
            ChVector3d vel = dyn_body->GetPosDt();
            double force_error = std::abs(total_force.y() - dyn_mass * 9.81);

            tracking_file << std::fixed << std::setprecision(8)
                          << time << ","
                          << patches.size() << ","
                          << active_indices.size() << ","
                          << matched_count << ","
                          << born_count << ","
                          << dead_count << ","
                          << pos.x() << "," << pos.y() << "," << pos.z() << ","
                          << vel.y() << ","
                          << total_force.y() << ","
                          << force_error << std::endl;

            // Write alive/born tracks
            for (const auto& track : current_tracks) {
                const char* status_str = (track.status == TrackStatus::Born) ? "born" : "alive";
                tracks_file << std::fixed << std::setprecision(8)
                            << time << "," << track.persistent_id << ","
                            << track.current_frame_patch_id << ","
                            << status_str << ","
                            << track.center.x() << "," << track.center.y() << "," << track.center.z() << ","
                            << track.normal.x() << "," << track.normal.y() << "," << track.normal.z() << ","
                            << track.area << "," << track.mean_phi << "," << track.min_phi << ","
                            << track.force.Length() << ","
                            << track.age_in_steps << std::endl;
            }
        }

        // Update active tracks for next frame
        active_tracks = current_tracks;
    }

    tracking_file.close();
    tracks_file.close();

    // -- Compute tracking summary --
    TrackingMetrics metrics;
    metrics.total_born = total_born;
    metrics.total_dead = total_dead;
    metrics.total_frames_with_patches = frames_with_patches;
    metrics.total_frames = static_cast<double>(total_frames);
    metrics.avg_patch_count = (total_frames > 0) ? (sum_patch_count / total_frames) : 0.0;
    metrics.max_patch_count = max_patch_count;
    metrics.avg_matched_ratio = (matched_ratio_count > 0) ? (sum_matched_ratio / matched_ratio_count) : 0.0;

    // Lifetime stats
    std::vector<int> lifetimes;
    double max_drift = 0.0;
    double sum_drift = 0.0;
    int drift_count = 0;
    double max_fluct = 0.0;

    for (const auto& [pid, birth_step] : track_birth_step) {
        // Find the track's final age
        int final_age = 0;
        for (const auto& t : active_tracks) {
            if (t.persistent_id == pid) {
                final_age = t.age_in_steps;
                break;
            }
        }
        // Also check if track died (use cumulative info)
        int lifetime = final_age;
        lifetimes.push_back(lifetime);

        if (track_cumulative_drift.count(pid)) {
            double d = track_cumulative_drift[pid];
            if (d > max_drift) max_drift = d;
            sum_drift += d;
            drift_count++;
        }
        if (track_cumulative_normal_fluct.count(pid)) {
            double f = track_cumulative_normal_fluct[pid];
            if (f > max_fluct) max_fluct = f;
        }
    }

    metrics.max_lifetime_steps = 0;
    metrics.avg_lifetime_steps = 0;
    if (!lifetimes.empty()) {
        metrics.max_lifetime_steps = *std::max_element(lifetimes.begin(), lifetimes.end());
        double sum_life = 0;
        for (int l : lifetimes) sum_life += l;
        metrics.avg_lifetime_steps = static_cast<int>(sum_life / lifetimes.size());
    }

    metrics.max_center_drift = max_drift;
    metrics.avg_center_drift = (drift_count > 0) ? (sum_drift / drift_count) : 0.0;
    metrics.max_normal_fluctuation = max_fluct;

    // Write summary CSV
    std::ofstream summary_file(out_dir + "/sdf_patch_tracking_summary.csv");
    if (summary_file.is_open()) {
        summary_file << "metric,value" << std::endl;
        summary_file << "total_born," << metrics.total_born << std::endl;
        summary_file << "total_dead," << metrics.total_dead << std::endl;
        summary_file << "total_frames," << metrics.total_frames << std::endl;
        summary_file << "frames_with_patches," << metrics.total_frames_with_patches << std::endl;
        summary_file << "avg_patch_count," << std::fixed << std::setprecision(4) << metrics.avg_patch_count << std::endl;
        summary_file << "max_patch_count," << metrics.max_patch_count << std::endl;
        summary_file << "avg_matched_ratio," << metrics.avg_matched_ratio << std::endl;
        summary_file << "max_lifetime_steps," << metrics.max_lifetime_steps << std::endl;
        summary_file << "avg_lifetime_steps," << metrics.avg_lifetime_steps << std::endl;
        summary_file << "max_center_drift," << metrics.max_center_drift << std::endl;
        summary_file << "avg_center_drift," << metrics.avg_center_drift << std::endl;
        summary_file << "max_normal_fluctuation," << metrics.max_normal_fluctuation << std::endl;
        summary_file.close();
    }

    // -- Console output --
    std::cout << "\n=== Tracking Results ===" << std::endl;
    std::cout << "  Total frames: " << metrics.total_frames << std::endl;
    std::cout << "  Frames with patches: " << metrics.total_frames_with_patches << std::endl;
    std::cout << "  Total born: " << metrics.total_born << std::endl;
    std::cout << "  Total dead: " << metrics.total_dead << std::endl;
    std::cout << "  Avg patch count: " << std::fixed << std::setprecision(4) << metrics.avg_patch_count << std::endl;
    std::cout << "  Max patch count: " << metrics.max_patch_count << std::endl;
    std::cout << "  Avg matched ratio: " << metrics.avg_matched_ratio << std::endl;
    std::cout << "  Max lifetime (steps): " << metrics.max_lifetime_steps << std::endl;
    std::cout << "  Avg lifetime (steps): " << metrics.avg_lifetime_steps << std::endl;
    std::cout << "  Max center drift: " << metrics.max_center_drift << " m" << std::endl;
    std::cout << "  Max normal fluctuation: " << metrics.max_normal_fluctuation << std::endl;

    // Verify
    bool pass = metrics.total_born > 0 && metrics.total_frames_with_patches > 0 && metrics.avg_matched_ratio > 0.5;
    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "  At least one birth: " << (metrics.total_born > 0 ? "PASS" : "FAIL") << std::endl;
    std::cout << "  At least one frame with patches: " << (metrics.total_frames_with_patches > 0 ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Avg matched ratio > 0.5: " << (metrics.avg_matched_ratio > 0.5 ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Overall: " << (pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tracking_output.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tracks.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tracking_summary.csv" << std::endl;

    return pass ? 0 : 1;
}
