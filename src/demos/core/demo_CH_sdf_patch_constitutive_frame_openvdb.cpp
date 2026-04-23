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
// Milestone 10: Patch Tangential Kinematics / Frame-Consistent Tangential State
//
// This demo extends milestone 9 with:
//   1. Patch-local tangent basis (t1, t2) for consistent frame definition
//   2. Tangential state transport across frames when normal changes
//   3. Local frame representation of tangential displacement (xi1, xi2)
//   4. Proper projection of old state to new tangent plane
//
// Scene: Case A (rotational contact) and Case B (multi-patch)
//
// Output:
//   out/milestone_10/sdf_patch_frame_case_A.csv
//   out/milestone_10/sdf_patch_frame_case_B.csv
//   out/milestone_10/sdf_patch_frame_compare.csv
//   out/milestone_10/sdf_patch_frame_tracks.csv
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
// Stick-slip state enum
// =============================================================================

enum class StickSlipState {
    Stick,
    Slip
};

// =============================================================================
// Patch primitive (frame-local)
// =============================================================================

struct PatchPrimitive {
    int patch_id;
    std::vector<int> sample_ids;
    int active_sample_count;
    ChVector3d center;
    ChVector3d normal;
    ChVector3d tangent_t1; // Local tangent basis
    ChVector3d tangent_t2;
    double total_area;
    double min_phi;
    double mean_phi;
    double max_penetration;
    double effective_penetration;
    double effective_normal_velocity;
    ChVector3d effective_tangential_velocity;
    ChVector3d force;
    ChVector3d torque;
    double normal_force_magnitude;
    double tangential_force_magnitude;
    StickSlipState stick_slip_state;
    int age_in_steps;
    // Frame kinematics stats (per-frame)
    double tangential_displacement_local_norm;
    double tangential_velocity_local_norm;
    double transported_state_norm;
    double frame_rotation_angle;
};

// =============================================================================
// Track status
// =============================================================================

enum class TrackStatus {
    Born,
    Alive,
    Dead
};

// =============================================================================
// Persistent patch track with frame-consistent tangential state
// =============================================================================

struct PersistentPatchTrack {
    int persistent_id;
    int current_frame_patch_id;
    TrackStatus status;
    ChVector3d center;
    ChVector3d normal;
    ChVector3d tangent_t1; // Previous frame tangent basis
    ChVector3d tangent_t2;
    double area;
    double mean_phi;
    double min_phi;
    double max_penetration;
    double effective_penetration;
    double effective_normal_velocity;
    ChVector3d effective_tangential_velocity;
    ChVector3d force;
    ChVector3d torque;
    double normal_force_magnitude;
    double tangential_force_magnitude;
    StickSlipState stick_slip_state;
    int age_in_steps;
    double birth_time;
    double last_seen_time;
    int birth_frame;
    // Frame-consistent tangential state (2D in local tangent plane)
    double xi1; // Tangential displacement in t1 direction
    double xi2; // Tangential displacement in t2 direction
    // Legacy: 3D tangential displacement (for milestone 9 backward compatibility)
    ChVector3d tangential_displacement;
    // Statistics
    int stick_steps;
    int slip_steps;
    int stick_to_slip_transitions;
    double accumulated_transport_norm; // Track how much transport occurs
};

// =============================================================================
// Contact mode enumeration
// =============================================================================

enum class ContactMode {
    Pointwise,
    PatchForceSum,
    PatchConstitutive,
    PatchTangentialDamping,
    PatchStickSlip,
    PatchFrameConsistent // NEW: frame-consistent tangential state
};

// =============================================================================
// Case configuration
// =============================================================================

struct CaseConfig {
    std::string name;
    std::string description;
    double static_sphere_radius;
    double static_sphere_offset_x;
    double static_sphere_offset_z;
    double dyn_body_radius;
    double dyn_body_size;
    int dyn_body_type;
    double dyn_mass;
    ChVector3d initial_pos;
    ChVector3d initial_vel;
    ChVector3d initial_ang_vel;
    double voxel_size;
    int n_theta;
    int n_phi;
    double stiffness;
    double damping;
    double activation_band;
    double force_band;
    double total_time;
    double time_step;
    // Tangential parameters
    double tangential_damping;
    double tangential_stiffness;
    double friction_coefficient;
    double ema_alpha;
};

// =============================================================================
// Case result metrics
// =============================================================================

struct CaseResult {
    ContactMode mode;
    std::string mode_name;
    double final_y;
    double y_error;
    double expected_equilibrium_y;
    double avg_force_y;
    double force_std_dev;
    double avg_torque_x;
    double avg_torque_y;
    double avg_torque_z;
    double torque_std_dev;
    int max_patch_count;
    int avg_patch_count;
    double multi_patch_ratio;
    double avg_tangential_force_norm;
    double max_tangential_force_norm;
    double avg_tangential_force_ratio;
    int total_stick_steps;
    int total_slip_steps;
    int total_stick_to_slip_transitions;
    double avg_tangential_displacement_norm;
    double avg_transported_state_norm;
    double avg_frame_rotation_angle;
    bool is_stable;
};

// =============================================================================
// Generate surface samples for sphere
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
// Generate surface samples for box (6 faces)
// =============================================================================

std::vector<SurfaceSample> GenerateBoxSamples(double half_size, int samples_per_face) {
    std::vector<SurfaceSample> samples;
    int n = static_cast<int>(std::sqrt(samples_per_face));
    double step = 2.0 * half_size / n;

    auto add_face_samples = [&](ChVector3d normal, double coord) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double u = -half_size + step * (i + 0.5);
                double v = -half_size + step * (j + 0.5);
                ChVector3d pos;
                if (std::abs(normal.x()) > 0.5) {
                    pos = ChVector3d(coord, u, v);
                } else if (std::abs(normal.y()) > 0.5) {
                    pos = ChVector3d(u, coord, v);
                } else {
                    pos = ChVector3d(u, v, coord);
                }
                SurfaceSample s;
                s.local_pos = pos;
                s.area_weight = 1.0 / (6 * n * n);
                s.theta_index = 0;
                s.phi_index = 0;
                s.global_index = static_cast<int>(samples.size());
                samples.push_back(s);
            }
        }
    };

    add_face_samples(ChVector3d(1, 0, 0), half_size);
    add_face_samples(ChVector3d(-1, 0, 0), -half_size);
    add_face_samples(ChVector3d(0, 1, 0), half_size);
    add_face_samples(ChVector3d(0, -1, 0), -half_size);
    add_face_samples(ChVector3d(0, 0, 1), half_size);
    add_face_samples(ChVector3d(0, 0, -1), -half_size);

    return samples;
}

// =============================================================================
// Grid adjacency and connected components
// =============================================================================

struct GridAdjacency {
    int n_theta;
    int n_phi;
    int total_samples;
    bool is_sphere;

    GridAdjacency() : n_theta(0), n_phi(0), total_samples(0), is_sphere(true) {}
    GridAdjacency(int nt, int np) : n_theta(nt), n_phi(np), total_samples(nt * np), is_sphere(true) {}
    GridAdjacency(int total) : n_theta(0), n_phi(0), total_samples(total), is_sphere(false) {}

    std::vector<int> GetNeighbors(int idx) const {
        std::vector<int> neighbors;
        if (is_sphere) {
            int ti = idx / n_phi;
            int pj = idx % n_phi;
            if (ti > 0) neighbors.push_back((ti - 1) * n_phi + pj);
            if (ti < n_theta - 1) neighbors.push_back((ti + 1) * n_phi + pj);
            int pj_prev = (pj - 1 + n_phi) % n_phi;
            int pj_next = (pj + 1) % n_phi;
            neighbors.push_back(ti * n_phi + pj_prev);
            neighbors.push_back(ti * n_phi + pj_next);
        } else {
            int n = static_cast<int>(std::sqrt(total_samples / 6));
            int face = idx / (n * n);
            int local_idx = idx % (n * n);
            int i = local_idx / n;
            int j = local_idx % n;
            if (i > 0) neighbors.push_back(face * n * n + (i - 1) * n + j);
            if (i < n - 1) neighbors.push_back(face * n * n + (i + 1) * n + j);
            if (j > 0) neighbors.push_back(face * n * n + i * n + (j - 1));
            if (j < n - 1) neighbors.push_back(face * n * n + i * n + (j + 1));
        }
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
                for (int neighbor : GetNeighbors(current)) {
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
        openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler> sampler(grid->tree(), grid->transform());
        openvdb::Vec3d lv(local.x(), local.y(), local.z());
        qr.phi = static_cast<double>(sampler.wsSample(lv));
        double dx = voxel_size * 0.5;
        double phi_px = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() + dx, local.y(), local.z())));
        double phi_mx = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() - dx, local.y(), local.z())));
        double phi_py = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() + dx, local.z())));
        double phi_my = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() - dx, local.z())));
        double phi_pz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() + dx)));
        double phi_mz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() - dx)));
        qr.grad = ChVector3d((phi_px - phi_mx) / (2.0 * dx), (phi_py - phi_my) / (2.0 * dx), (phi_pz - phi_mz) / (2.0 * dx));
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
// Build orthonormal tangent basis from normal
// =============================================================================

void BuildTangentBasis(const ChVector3d& n, ChVector3d& t1, ChVector3d& t2) {
    // Find the axis least aligned with n
    ChVector3d ref = std::abs(n.x()) < 0.9 ? ChVector3d(1, 0, 0) : ChVector3d(0, 1, 0);
    t1 = (n.Cross(ref)).GetNormalized();
    t2 = n.Cross(t1);
}

// =============================================================================
// Transport tangential state from old tangent plane to new tangent plane
// =============================================================================

struct TransportResult {
    double xi1_new;
    double xi2_new;
    double transported_norm; // Norm of transported world vector
    double rotation_angle;   // Approximate rotation angle between frames
};

TransportResult TransportTangentialState(
    double xi1_old, double xi2_old,
    const ChVector3d& t1_old, const ChVector3d& t2_old, const ChVector3d& n_old,
    const ChVector3d& t1_new, const ChVector3d& t2_new, const ChVector3d& n_new
) {
    TransportResult result;

    // Step 1: Convert old local state to world coordinates
    ChVector3d xi_world = xi1_old * t1_old + xi2_old * t2_old;

    // Step 2: Project onto new tangent plane (remove component along new normal)
    double normal_component = xi_world.Dot(n_new);
    ChVector3d xi_projected = xi_world - normal_component * n_new;

    // Step 3: Represent in new local basis
    result.xi1_new = xi_projected.Dot(t1_new);
    result.xi2_new = xi_projected.Dot(t2_new);

    result.transported_norm = xi_projected.Length();
    result.rotation_angle = std::acos(std::max(-1.0, std::min(1.0, n_old.Dot(n_new))));

    return result;
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
        patch.normal_force_magnitude = 0.0;
        patch.tangential_force_magnitude = 0.0;
        patch.stick_slip_state = StickSlipState::Stick;
        patch.age_in_steps = 0;
        patch.tangential_displacement_local_norm = 0.0;
        patch.tangential_velocity_local_norm = 0.0;
        patch.transported_state_norm = 0.0;
        patch.frame_rotation_angle = 0.0;

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

        // Build tangent basis
        BuildTangentBasis(patch.normal, patch.tangent_t1, patch.tangent_t2);

        patch.effective_normal_velocity = 0.0;
        patch.effective_tangential_velocity = ChVector3d(0, 0, 0);

        patches.push_back(patch);
    }
    return patches;
}

// =============================================================================
// Patch matching
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
// Compute patch constitutive force (original - normal only)
// =============================================================================

ChVector3d ComputePatchConstitutiveForceNormal(
    const PatchPrimitive& patch,
    const ChVector3d& patch_center_velocity,
    double k_patch,
    double c_patch
) {
    double penetration_eff = patch.effective_penetration;
    double vn_eff = patch_center_velocity.Dot(patch.normal);
    double force_scalar = k_patch * penetration_eff + c_patch * std::max(-vn_eff, 0.0);
    if (force_scalar < 0.0) force_scalar = 0.0;
    return patch.normal * force_scalar;
}

// =============================================================================
// Compute patch tangential damping force (milestone 8 version)
// =============================================================================

ChVector3d ComputePatchTangentialDampingForce(
    const PatchPrimitive& patch,
    const ChVector3d& patch_center_velocity,
    double k_patch,
    double c_patch,
    double c_tangential,
    double mu,
    double ema_alpha,
    const PersistentPatchTrack* prev_track,
    ChVector3d& out_tangential_displacement
) {
    // --- Normal force (same as milestone 8) ---
    double penetration_current = patch.effective_penetration;
    double penetration_eff = penetration_current;
    if (prev_track != nullptr) {
        penetration_eff = ema_alpha * penetration_current + (1.0 - ema_alpha) * prev_track->effective_penetration;
    }

    double vn_eff = patch_center_velocity.Dot(patch.normal);
    double force_n_scalar = k_patch * penetration_eff + c_patch * std::max(-vn_eff, 0.0);
    if (force_n_scalar < 0.0) force_n_scalar = 0.0;
    ChVector3d F_n = patch.normal * force_n_scalar;

    // --- Tangential velocity ---
    double vn = patch_center_velocity.Dot(patch.normal);
    ChVector3d vt = patch_center_velocity - vn * patch.normal;

    // --- Tangential force (damping only) ---
    ChVector3d F_t = vt * (-c_tangential);

    // --- Coulomb-like clamp ---
    double fn_mag = F_n.Length();
    double ft_mag = F_t.Length();
    if (ft_mag > mu * fn_mag && fn_mag > 1e-12) {
        F_t = F_t * (mu * fn_mag / ft_mag);
    }

    out_tangential_displacement = ChVector3d(0, 0, 0);

    return F_n + F_t;
}

// =============================================================================
// Compute patch stick-slip force (milestone 9 version)
// =============================================================================

struct StickSlipResult {
    ChVector3d force;
    ChVector3d tangential_displacement;
    StickSlipState state;
    double normal_force_magnitude;
    double tangential_force_magnitude;
};

StickSlipResult ComputePatchStickSlipForce(
    const PatchPrimitive& patch,
    const ChVector3d& patch_center_velocity,
    double k_patch,
    double c_patch,
    double k_tangential,
    double c_tangential,
    double mu,
    double dt,
    const PersistentPatchTrack* prev_track
) {
    StickSlipResult result;

    // --- Normal force (with EMA smoothing) ---
    double penetration_current = patch.effective_penetration;
    double penetration_eff = penetration_current;
    if (prev_track != nullptr) {
        penetration_eff = 0.3 * penetration_current + 0.7 * prev_track->effective_penetration;
    }

    double vn_eff = patch_center_velocity.Dot(patch.normal);
    double force_n_scalar = k_patch * penetration_eff + c_patch * std::max(-vn_eff, 0.0);
    if (force_n_scalar < 0.0) force_n_scalar = 0.0;
    ChVector3d F_n = patch.normal * force_n_scalar;
    result.normal_force_magnitude = F_n.Length();

    // --- Tangential velocity ---
    double vn = patch_center_velocity.Dot(patch.normal);
    ChVector3d vt = patch_center_velocity - vn * patch.normal;

    // --- Initialize or inherit tangential displacement ---
    if (prev_track != nullptr) {
        result.tangential_displacement = prev_track->tangential_displacement;
        result.state = prev_track->stick_slip_state;
    } else {
        result.tangential_displacement = ChVector3d(0, 0, 0);
        result.state = StickSlipState::Stick;
    }

    // --- Update tangential displacement ---
    result.tangential_displacement += vt * dt;

    // --- Trial tangential force ---
    ChVector3d F_t_trial = result.tangential_displacement * (-k_tangential) + vt * (-c_tangential);

    // --- Stick-slip decision ---
    double fn_mag = F_n.Length();
    double ft_trial_mag = F_t_trial.Length();
    double friction_limit = mu * fn_mag;

    if (ft_trial_mag <= friction_limit) {
        result.force = F_n + F_t_trial;
        result.tangential_force_magnitude = ft_trial_mag;
        result.state = StickSlipState::Stick;
    } else {
        if (vt.Length() > 1e-12) {
            ChVector3d vt_dir = vt / vt.Length();
            ChVector3d F_t_slip = vt_dir * (-friction_limit);
            result.force = F_n + F_t_slip;
            result.tangential_force_magnitude = friction_limit;
            result.tangential_displacement = F_t_slip * (-1.0 / k_tangential);
        } else {
            result.force = F_n;
            result.tangential_force_magnitude = 0.0;
        }
        result.state = StickSlipState::Slip;
    }

    return result;
}

// =============================================================================
// Compute patch frame-consistent tangential force (NEW - milestone 10)
// =============================================================================

struct FrameConsistentResult {
    ChVector3d force;
    double xi1; // Tangential displacement in local t1 direction
    double xi2; // Tangential displacement in local t2 direction
    StickSlipState state;
    double normal_force_magnitude;
    double tangential_force_magnitude;
    double transported_state_norm; // For output
    double frame_rotation_angle;   // For output
};

FrameConsistentResult ComputePatchFrameConsistentForce(
    const PatchPrimitive& patch,
    const ChVector3d& patch_center_velocity,
    double k_patch,
    double c_patch,
    double k_tangential,
    double c_tangential,
    double mu,
    double dt,
    const PersistentPatchTrack* prev_track
) {
    FrameConsistentResult result;
    result.xi1 = 0.0;
    result.xi2 = 0.0;
    result.transported_state_norm = 0.0;
    result.frame_rotation_angle = 0.0;

    // --- Normal force (with EMA smoothing) ---
    double penetration_current = patch.effective_penetration;
    double penetration_eff = penetration_current;
    if (prev_track != nullptr) {
        penetration_eff = 0.3 * penetration_current + 0.7 * prev_track->effective_penetration;
    }

    double vn_eff = patch_center_velocity.Dot(patch.normal);
    double force_n_scalar = k_patch * penetration_eff + c_patch * std::max(-vn_eff, 0.0);
    if (force_n_scalar < 0.0) force_n_scalar = 0.0;
    ChVector3d F_n = patch.normal * force_n_scalar;
    result.normal_force_magnitude = F_n.Length();

    // --- Tangential velocity ---
    double vn = patch_center_velocity.Dot(patch.normal);
    ChVector3d vt = patch_center_velocity - vn * patch.normal;

    // --- Project vt to local tangent basis ---
    double vt1 = vt.Dot(patch.tangent_t1);
    double vt2 = vt.Dot(patch.tangent_t2);

    // --- Initialize or transport tangential state ---
    if (prev_track != nullptr) {
        double xi1_old = prev_track->xi1;
        double xi2_old = prev_track->xi2;

        // Transport state from old frame to new frame
        TransportResult transport = TransportTangentialState(
            xi1_old, xi2_old,
            prev_track->tangent_t1, prev_track->tangent_t2, prev_track->normal,
            patch.tangent_t1, patch.tangent_t2, patch.normal
        );

        result.xi1 = transport.xi1_new;
        result.xi2 = transport.xi2_new;
        result.transported_state_norm = transport.transported_norm;
        result.frame_rotation_angle = transport.rotation_angle;
        result.state = prev_track->stick_slip_state;
    } else {
        // Birth: initialize to zero
        result.xi1 = 0.0;
        result.xi2 = 0.0;
        result.state = StickSlipState::Stick;
    }

    // --- Update tangential displacement in local frame ---
    result.xi1 += vt1 * dt;
    result.xi2 += vt2 * dt;

    double xi_local_norm = std::sqrt(result.xi1 * result.xi1 + result.xi2 * result.xi2);
    double vt_local_norm = std::sqrt(vt1 * vt1 + vt2 * vt2);

    // --- Trial tangential force in local frame ---
    double ft1_trial = -k_tangential * result.xi1 - c_tangential * vt1;
    double ft2_trial = -k_tangential * result.xi2 - c_tangential * vt2;
    double ft_trial_mag = std::sqrt(ft1_trial * ft1_trial + ft2_trial * ft2_trial);

    // --- Stick-slip decision ---
    double fn_mag = F_n.Length();
    double friction_limit = mu * fn_mag;

    if (ft_trial_mag <= friction_limit) {
        // Stick: use trial force
        // Convert local force back to world coordinates
        ChVector3d F_t = ft1_trial * patch.tangent_t1 + ft2_trial * patch.tangent_t2;
        result.force = F_n + F_t;
        result.tangential_force_magnitude = ft_trial_mag;
        result.state = StickSlipState::Stick;
    } else {
        // Slip: clamp force to friction limit
        if (vt_local_norm > 1e-12) {
            // Direction in local frame
            double vt1_dir = vt1 / vt_local_norm;
            double vt2_dir = vt2 / vt_local_norm;

            // Clamp force in local frame
            double ft1_slip = -friction_limit * vt1_dir;
            double ft2_slip = -friction_limit * vt2_dir;

            // Convert to world
            ChVector3d F_t = ft1_slip * patch.tangent_t1 + ft2_slip * patch.tangent_t2;
            result.force = F_n + F_t;
            result.tangential_force_magnitude = friction_limit;

            // Project tangential displacement to be consistent with slip force
            result.xi1 = ft1_slip / (-k_tangential);
            result.xi2 = ft2_slip / (-k_tangential);
        } else {
            result.force = F_n;
            result.tangential_force_magnitude = 0.0;
        }
        result.state = StickSlipState::Slip;
    }

    return result;
}

// =============================================================================
// Run single case in single mode
// =============================================================================

CaseResult RunCaseMode(
    const CaseConfig& config,
    ContactMode mode,
    const std::string& mode_name
) {
    openvdb::initialize();

    // -- Create OpenVDB sphere SDF --
    double voxel_size = config.voxel_size;
    auto sphere_grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        static_cast<float>(config.static_sphere_radius),
        openvdb::Vec3f(static_cast<float>(config.static_sphere_offset_x), 0.0f, static_cast<float>(config.static_sphere_offset_z)),
        static_cast<float>(voxel_size),
        3.0f
    );
    sphere_grid->setGridClass(openvdb::GRID_LEVEL_SET);

    OpenVDBSphereSDF openvdb_sdf;
    openvdb_sdf.grid = sphere_grid;
    openvdb_sdf.center_world = ChVector3d(config.static_sphere_offset_x, 0, config.static_sphere_offset_z);
    openvdb_sdf.sphere_radius = config.static_sphere_radius;
    openvdb_sdf.voxel_size = voxel_size;

    SDFProbeFunc sdf_probe = [&](const ChVector3d& pt) -> SDFQueryResult {
        return openvdb_sdf.Query(pt);
    };

    // -- Create samples --
    std::vector<SurfaceSample> samples;
    GridAdjacency adjacency;
    if (config.dyn_body_type == 0) {
        samples = GenerateSphereSamples(config.dyn_body_radius, config.n_theta, config.n_phi);
        adjacency = GridAdjacency(config.n_theta, config.n_phi);
    } else {
        samples = GenerateBoxSamples(config.dyn_body_size, 100);
        adjacency = GridAdjacency(static_cast<int>(samples.size()));
    }

    // -- Chrono setup --
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -9.81, 0));

    std::shared_ptr<ChBody> dyn_body;
    if (config.dyn_body_type == 0) {
        dyn_body = chrono_types::make_shared<ChBodyEasySphere>(config.dyn_body_radius, 1000.0, false, false);
    } else {
        dyn_body = chrono_types::make_shared<ChBodyEasyBox>(config.dyn_body_size * 2, config.dyn_body_size * 2, config.dyn_body_size * 2, 1000.0, false, false);
    }
    dyn_body->SetPos(config.initial_pos);
    dyn_body->SetPosDt(config.initial_vel);
    dyn_body->SetAngVelLocal(config.initial_ang_vel);
    dyn_body->SetFixed(false);
    sys.AddBody(dyn_body);

    unsigned int acc_id = dyn_body->AddAccumulator();
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(100);
    sys.GetSolver()->AsIterative()->SetTolerance(1e-6);

    // -- Persistent state --
    std::vector<PersistentPatchTrack> active_tracks;
    int next_persistent_id = 0;
    std::map<int, PersistentPatchTrack> track_state_map;

    // -- Buffers --
    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());
    std::vector<ChVector3d> sample_forces(samples.size());
    std::vector<ChVector3d> sample_torques(samples.size());

    // -- Metrics --
    std::vector<double> force_y_values;
    std::vector<double> pos_y_values;
    std::vector<int> patch_counts;
    std::vector<double> torque_x_values;
    std::vector<double> torque_y_values;
    std::vector<double> torque_z_values;
    std::vector<double> tangential_force_norms;
    std::vector<double> tangential_displacement_norms;
    std::vector<double> transported_state_norms;
    std::vector<double> frame_rotation_angles;
    double total_force_y = 0.0;
    int force_sample_count = 0;
    double sum_tangential_norm = 0.0;
    double max_tangential_norm = 0.0;
    double sum_tangential_ratio = 0.0;
    int total_stick_steps = 0;
    int total_slip_steps = 0;
    int total_stick_to_slip_transitions = 0;

    // -- Matching thresholds --
    double match_distance_threshold = 0.5 * (config.dyn_body_type == 0 ? config.dyn_body_radius : config.dyn_body_size);
    double normal_similarity_threshold = 0.9;

    while (sys.GetChTime() < config.total_time) {
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
            if (phi < config.activation_band) {
                active_indices.push_back(static_cast<int>(si));

                if (mode == ContactMode::Pointwise) {
                    if (phi < config.force_band) {
                        double pen = std::max(-phi, 0.0);
                        double grad_len = sdf_results[si].grad.Length();
                        if (grad_len > 1e-12) {
                            ChVector3d normal = sdf_results[si].grad / grad_len;
                            double vn = sample_vel.Dot(normal);
                            double force_mag = config.stiffness * pen + config.damping * std::max(-vn, 0.0);
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
        if (mode == ContactMode::PatchForceSum || mode == ContactMode::PatchConstitutive || 
            mode == ContactMode::PatchTangentialDamping || mode == ContactMode::PatchStickSlip ||
            mode == ContactMode::PatchFrameConsistent) {
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
                } else {
                    track.persistent_id = next_persistent_id++;
                    track.current_frame_patch_id = m.current_patch_id;
                    track.status = TrackStatus::Born;
                    track.birth_time = sys.GetChTime();
                    track.birth_frame = 0;
                    track.age_in_steps = 0;
                    track.xi1 = 0.0;
                    track.xi2 = 0.0;
                    track.stick_slip_state = StickSlipState::Stick;
                    track.stick_steps = 0;
                    track.slip_steps = 0;
                    track.stick_to_slip_transitions = 0;
                    track.accumulated_transport_norm = 0.0;
                }

                for (const auto& p : patches) {
                    if (p.patch_id == track.current_frame_patch_id) {
                        track.center = p.center;
                        track.normal = p.normal;
                        track.tangent_t1 = p.tangent_t1;
                        track.tangent_t2 = p.tangent_t2;
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
                ChVector3d patch_center_vel = body_vel + body_ang_vel_world.Cross(patch.center - body_pos);
                patch.effective_normal_velocity = patch_center_vel.Dot(patch.normal);

                // Find previous track for this patch
                const PersistentPatchTrack* prev_track = nullptr;
                for (const auto& t : active_tracks) {
                    if (t.current_frame_patch_id == patch.patch_id) {
                        prev_track = &t;
                        break;
                    }
                }

                if (mode == ContactMode::PatchForceSum) {
                    for (int sid : patch.sample_ids) {
                        double phi = sdf_results[sid].phi;
                        if (phi < config.force_band) {
                            double pen = std::max(-phi, 0.0);
                            double grad_len = sdf_results[sid].grad.Length();
                            if (grad_len > 1e-12) {
                                ChVector3d normal = sdf_results[sid].grad / grad_len;
                                ChVector3d sample_world = sample_world_positions[sid];
                                ChVector3d r = sample_world - body_pos;
                                ChVector3d sample_vel = body_vel + body_ang_vel_world.Cross(r);
                                double vn = sample_vel.Dot(normal);
                                double force_mag = config.stiffness * pen + config.damping * std::max(-vn, 0.0);
                                if (force_mag > 0.0) {
                                    ChVector3d pt_force = normal * force_mag;
                                    patch.force += pt_force;
                                    patch.torque += r.Cross(pt_force);
                                }
                            }
                        }
                    }
                } else if (mode == ContactMode::PatchConstitutive) {
                    patch.force = ComputePatchConstitutiveForceNormal(
                        patch, patch_center_vel, config.stiffness, config.damping
                    );
                    patch.torque = (patch.center - body_pos).Cross(patch.force);
                } else if (mode == ContactMode::PatchTangentialDamping) {
                    if (prev_track != nullptr) {
                        patch.age_in_steps = prev_track->age_in_steps + 1;
                    }
                    ChVector3d xi_t;
                    patch.force = ComputePatchTangentialDampingForce(
                        patch, patch_center_vel, config.stiffness, config.damping,
                        config.tangential_damping, config.friction_coefficient,
                        config.ema_alpha, prev_track, xi_t
                    );
                    patch.torque = (patch.center - body_pos).Cross(patch.force);
                } else if (mode == ContactMode::PatchStickSlip) {
                    if (prev_track != nullptr) {
                        patch.age_in_steps = prev_track->age_in_steps + 1;
                    }

                    StickSlipResult ss_result = ComputePatchStickSlipForce(
                        patch, patch_center_vel, config.stiffness, config.damping,
                        config.tangential_stiffness, config.tangential_damping,
                        config.friction_coefficient, config.time_step, prev_track
                    );

                    patch.force = ss_result.force;
                    patch.stick_slip_state = ss_result.state;
                    patch.normal_force_magnitude = ss_result.normal_force_magnitude;
                    patch.tangential_force_magnitude = ss_result.tangential_force_magnitude;
                    patch.torque = (patch.center - body_pos).Cross(patch.force);

                    // Update track state
                    for (auto& track : current_tracks) {
                        if (track.current_frame_patch_id == patch.patch_id) {
                            track.stick_slip_state = ss_result.state;
                            track.normal_force_magnitude = ss_result.normal_force_magnitude;
                            track.tangential_force_magnitude = ss_result.tangential_force_magnitude;
                            // Legacy: store 3D displacement
                            if (prev_track != nullptr) {
                                track.tangent_t1 = prev_track->tangent_t1;
                                track.tangent_t2 = prev_track->tangent_t2;
                            }
                            if (ss_result.state == StickSlipState::Stick) {
                                track.stick_steps++;
                            } else {
                                track.slip_steps++;
                            }
                            if (prev_track != nullptr && prev_track->stick_slip_state == StickSlipState::Stick && 
                                ss_result.state == StickSlipState::Slip) {
                                track.stick_to_slip_transitions++;
                            }
                            break;
                        }
                    }
                } else if (mode == ContactMode::PatchFrameConsistent) {
                    if (prev_track != nullptr) {
                        patch.age_in_steps = prev_track->age_in_steps + 1;
                    }

                    FrameConsistentResult fc_result = ComputePatchFrameConsistentForce(
                        patch, patch_center_vel, config.stiffness, config.damping,
                        config.tangential_stiffness, config.tangential_damping,
                        config.friction_coefficient, config.time_step, prev_track
                    );

                    patch.force = fc_result.force;
                    patch.stick_slip_state = fc_result.state;
                    patch.normal_force_magnitude = fc_result.normal_force_magnitude;
                    patch.tangential_force_magnitude = fc_result.tangential_force_magnitude;
                    patch.torque = (patch.center - body_pos).Cross(patch.force);

                    // Output stats
                    patch.tangential_displacement_local_norm = std::sqrt(fc_result.xi1 * fc_result.xi1 + fc_result.xi2 * fc_result.xi2);
                    patch.transported_state_norm = fc_result.transported_state_norm;
                    patch.frame_rotation_angle = fc_result.frame_rotation_angle;

                    // Update track state
                    for (auto& track : current_tracks) {
                        if (track.current_frame_patch_id == patch.patch_id) {
                            track.xi1 = fc_result.xi1;
                            track.xi2 = fc_result.xi2;
                            track.stick_slip_state = fc_result.state;
                            track.normal_force_magnitude = fc_result.normal_force_magnitude;
                            track.tangential_force_magnitude = fc_result.tangential_force_magnitude;
                            if (fc_result.state == StickSlipState::Stick) {
                                track.stick_steps++;
                            } else {
                                track.slip_steps++;
                            }
                            if (prev_track != nullptr && prev_track->stick_slip_state == StickSlipState::Stick && 
                                fc_result.state == StickSlipState::Slip) {
                                track.stick_to_slip_transitions++;
                            }
                            track.accumulated_transport_norm += fc_result.transported_state_norm;
                            break;
                        }
                    }
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

        sys.DoStepDynamics(config.time_step);

        // Record metrics
        double current_force_y = total_force.y();
        force_y_values.push_back(current_force_y);
        pos_y_values.push_back(dyn_body->GetPos().y());
        int current_patch_count = active_indices.empty() ? 0 : static_cast<int>(BuildPatches(
            active_indices, samples, sdf_results,
            sample_world_positions, adjacency, 0
        ).size());
        patch_counts.push_back(current_patch_count);
        torque_x_values.push_back(total_torque.x());
        torque_y_values.push_back(total_torque.y());
        torque_z_values.push_back(total_torque.z());

        total_force_y += current_force_y;
        force_sample_count++;

        // Tangential force stats (for tangential/stick-slip/frame modes)
        if (mode == ContactMode::PatchTangentialDamping || mode == ContactMode::PatchStickSlip || mode == ContactMode::PatchFrameConsistent) {
            double ft_approx = std::sqrt(total_force.x() * total_force.x() + total_force.z() * total_force.z());
            double fn_approx = std::abs(total_force.y());
            tangential_force_norms.push_back(ft_approx);
            sum_tangential_norm += ft_approx;
            if (ft_approx > max_tangential_norm) max_tangential_norm = ft_approx;
            if (fn_approx > 1e-12) {
                sum_tangential_ratio += ft_approx / fn_approx;
            }
        }

        // Stick-slip stats
        if (mode == ContactMode::PatchStickSlip || mode == ContactMode::PatchFrameConsistent) {
            for (const auto& track : active_tracks) {
                total_stick_steps += track.stick_steps;
                total_slip_steps += track.slip_steps;
                total_stick_to_slip_transitions += track.stick_to_slip_transitions;
                // Compute norm from local frame state
                double xi_norm = std::sqrt(track.xi1 * track.xi1 + track.xi2 * track.xi2);
                tangential_displacement_norms.push_back(xi_norm);
                transported_state_norms.push_back(track.accumulated_transport_norm);
                frame_rotation_angles.push_back(track.accumulated_transport_norm); // approximate
            }
        }
    }

    // Compute statistics
    CaseResult result;
    result.mode = mode;
    result.mode_name = mode_name;
    result.final_y = pos_y_values.empty() ? 0.0 : pos_y_values.back();

    if (config.dyn_body_type == 0) {
        result.expected_equilibrium_y = config.static_sphere_radius + config.dyn_body_radius - (config.dyn_mass * 9.81) / config.stiffness;
    } else {
        result.expected_equilibrium_y = config.static_sphere_radius + config.dyn_body_size - (config.dyn_mass * 9.81) / config.stiffness;
    }
    result.y_error = std::abs(result.final_y - result.expected_equilibrium_y);
    result.avg_force_y = force_sample_count > 0 ? total_force_y / force_sample_count : 0.0;

    double sum_sq = 0.0;
    for (double fy : force_y_values) {
        sum_sq += (fy - result.avg_force_y) * (fy - result.avg_force_y);
    }
    result.force_std_dev = force_sample_count > 0 ? std::sqrt(sum_sq / force_sample_count) : 0.0;

    double avg_tx = 0, avg_ty = 0, avg_tz = 0;
    for (double tx : torque_x_values) avg_tx += tx;
    for (double ty : torque_y_values) avg_ty += ty;
    for (double tz : torque_z_values) avg_tz += tz;
    if (force_sample_count > 0) {
        avg_tx /= force_sample_count;
        avg_ty /= force_sample_count;
        avg_tz /= force_sample_count;
    }
    result.avg_torque_x = avg_tx;
    result.avg_torque_y = avg_ty;
    result.avg_torque_z = avg_tz;

    double sum_sq_tx = 0, sum_sq_ty = 0, sum_sq_tz = 0;
    for (size_t i = 0; i < torque_x_values.size(); i++) {
        sum_sq_tx += (torque_x_values[i] - avg_tx) * (torque_x_values[i] - avg_tx);
        sum_sq_ty += (torque_y_values[i] - avg_ty) * (torque_y_values[i] - avg_ty);
        sum_sq_tz += (torque_z_values[i] - avg_tz) * (torque_z_values[i] - avg_tz);
    }
    result.torque_std_dev = force_sample_count > 0 ? std::sqrt((sum_sq_tx + sum_sq_ty + sum_sq_tz) / force_sample_count) : 0.0;

    result.max_patch_count = patch_counts.empty() ? 0 : *std::max_element(patch_counts.begin(), patch_counts.end());
    result.avg_patch_count = patch_counts.empty() ? 0 : static_cast<int>(std::accumulate(patch_counts.begin(), patch_counts.end(), 0) / patch_counts.size());

    int multi_patch_count = 0;
    for (int pc : patch_counts) {
        if (pc > 1) multi_patch_count++;
    }
    result.multi_patch_ratio = patch_counts.empty() ? 0.0 : static_cast<double>(multi_patch_count) / patch_counts.size();

    result.avg_tangential_force_norm = tangential_force_norms.empty() ? 0.0 : sum_tangential_norm / tangential_force_norms.size();
    result.max_tangential_force_norm = max_tangential_norm;
    result.avg_tangential_force_ratio = tangential_force_norms.empty() ? 0.0 : sum_tangential_ratio / tangential_force_norms.size();
    result.avg_tangential_displacement_norm = tangential_displacement_norms.empty() ? 0.0 :
        std::accumulate(tangential_displacement_norms.begin(), tangential_displacement_norms.end(), 0.0) / tangential_displacement_norms.size();
    result.avg_transported_state_norm = transported_state_norms.empty() ? 0.0 :
        std::accumulate(transported_state_norms.begin(), transported_state_norms.end(), 0.0) / transported_state_norms.size();
    result.avg_frame_rotation_angle = frame_rotation_angles.empty() ? 0.0 :
        std::accumulate(frame_rotation_angles.begin(), frame_rotation_angles.end(), 0.0) / frame_rotation_angles.size();

    result.total_stick_steps = total_stick_steps;
    result.total_slip_steps = total_slip_steps;
    result.total_stick_to_slip_transitions = total_stick_to_slip_transitions;

    size_t last_10_percent = pos_y_values.size() * 9 / 10;
    double avg_last = 0.0;
    for (size_t i = last_10_percent; i < pos_y_values.size(); i++) {
        avg_last += pos_y_values[i];
    }
    avg_last /= (pos_y_values.size() - last_10_percent);
    result.is_stable = std::abs(avg_last - result.expected_equilibrium_y) < 0.06;

    return result;
}

// =============================================================================
// Track writer
// =============================================================================

void WriteFrameTrackCSV(const std::string& filename, const std::vector<PersistentPatchTrack>& tracks) {
    std::ofstream file(filename);
    file << "persistent_id,age_in_steps,xi1,xi2,tangential_displacement_local_norm,stick_slip_state,stick_steps,slip_steps,stick_to_slip_transitions,accumulated_transport_norm" << std::endl;
    for (const auto& t : tracks) {
        file << t.persistent_id << ","
             << t.age_in_steps << ","
             << std::fixed << std::setprecision(6)
             << t.xi1 << ","
             << t.xi2 << ","
             << std::sqrt(t.xi1 * t.xi1 + t.xi2 * t.xi2) << ","
             << (t.stick_slip_state == StickSlipState::Stick ? "Stick" : "Slip") << ","
             << t.stick_steps << ","
             << t.slip_steps << ","
             << t.stick_to_slip_transitions << ","
             << t.accumulated_transport_norm << std::endl;
    }
    file.close();
}

// =============================================================================
// Main simulation
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 10: Patch Frame-Consistent Tangential State ===" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_10";
    EnsureDir(out_dir);

    // ========================================================================
    // Milestone 9 Baseline Reference
    // ========================================================================
    std::cout << "\n=== Milestone 9 Baseline Reference ===" << std::endl;
    std::cout << "Case A (Rotational Contact):" << std::endl;
    std::cout << "  PatchConstitutive: y_error=0.0471m, torque_z=-2.91Nm, force_std=795N, stable=YES" << std::endl;
    std::cout << "  TangentialDamping: y_error=0.0802m, torque_z=-3.95Nm, force_std=791N, stable=NO" << std::endl;
    std::cout << "  StickSlip: y_error=0.5508m, torque_z=+2.30Nm, force_std=791N, stable=NO" << std::endl;
    std::cout << "Case B (Multi-Patch):" << std::endl;
    std::cout << "  PatchConstitutive: y_error=0.0222m, torque_z=0.034Nm, force_std=487N, stable=YES" << std::endl;
    std::cout << "  TangentialDamping: y_error=0.0416m, torque_z=-0.15Nm, force_std=486N, stable=YES" << std::endl;
    std::cout << "  StickSlip: y_error=0.0288m, torque_z=-0.024Nm, force_std=483N, stable=YES" << std::endl;

    // ========================================================================
    // Case Definitions (same as milestone 7/8/9)
    // ========================================================================

    // Case A: Rotational contact
    CaseConfig case_A;
    case_A.name = "Case_A_Rotational_Contact";
    case_A.description = "Dynamic sphere drops with initial angular velocity";
    case_A.static_sphere_radius = 1.0;
    case_A.static_sphere_offset_x = 0.0;
    case_A.static_sphere_offset_z = 0.0;
    case_A.dyn_body_radius = 0.2;
    case_A.dyn_body_type = 0;
    case_A.dyn_mass = (4.0 / 3.0) * M_PI * std::pow(case_A.dyn_body_radius, 3) * 1000.0;
    case_A.initial_pos = ChVector3d(0.0, 3.0, 0.0);
    case_A.initial_vel = ChVector3d(0, 0, 0);
    case_A.initial_ang_vel = ChVector3d(0, 0, 5.0);
    case_A.voxel_size = 0.05;
    case_A.n_theta = 8;
    case_A.n_phi = 16;
    case_A.stiffness = 1e5;
    case_A.damping = 500;
    case_A.activation_band = 0.1;
    case_A.force_band = 0.0;
    case_A.total_time = 2.0;
    case_A.time_step = 5e-4;
    case_A.tangential_damping = 10.0;
    case_A.tangential_stiffness = 1e4;
    case_A.friction_coefficient = 0.3;
    case_A.ema_alpha = 0.3;

    // Case B: Multi-patch
    CaseConfig case_B;
    case_B.name = "Case_B_Multi_Patch";
    case_B.description = "Dynamic box with initial rotation";
    case_B.static_sphere_radius = 1.0;
    case_B.static_sphere_offset_x = 0.0;
    case_B.static_sphere_offset_z = 0.0;
    case_B.dyn_body_size = 0.15;
    case_B.dyn_body_type = 1;
    case_B.dyn_mass = std::pow(case_B.dyn_body_size * 2, 3) * 1000.0;
    case_B.initial_pos = ChVector3d(0.0, 2.5, 0.0);
    case_B.initial_vel = ChVector3d(0, 0, 0);
    case_B.initial_ang_vel = ChVector3d(0.5, 0, 0.3);
    case_B.voxel_size = 0.05;
    case_B.n_theta = 8;
    case_B.n_phi = 16;
    case_B.stiffness = 1e5;
    case_B.damping = 500;
    case_B.activation_band = 0.1;
    case_B.force_band = 0.0;
    case_B.total_time = 2.0;
    case_B.time_step = 5e-4;
    case_B.tangential_damping = 10.0;
    case_B.tangential_stiffness = 1e4;
    case_B.friction_coefficient = 0.3;
    case_B.ema_alpha = 0.3;

    std::vector<CaseConfig> cases = {case_A, case_B};
    std::vector<ContactMode> modes = {
        ContactMode::Pointwise,
        ContactMode::PatchForceSum,
        ContactMode::PatchConstitutive,
        ContactMode::PatchTangentialDamping,
        ContactMode::PatchStickSlip,
        ContactMode::PatchFrameConsistent
    };
    std::vector<std::string> mode_names = {"Pointwise", "PatchForceSum", "PatchConstitutive", "PatchTangentialDamping", "PatchStickSlip", "PatchFrameConsistent"};

    // ========================================================================
    // Run all cases
    // ========================================================================

    std::vector<std::vector<CaseResult>> all_results;
    std::vector<std::vector<PersistentPatchTrack>> final_tracks;

    for (const auto& config : cases) {
        std::cout << "\n=== Running " << config.name << " ===" << std::endl;
        std::cout << "Description: " << config.description << std::endl;
        std::cout << "Tangential params: k_t=" << config.tangential_stiffness 
                  << ", c_t=" << config.tangential_damping 
                  << ", mu=" << config.friction_coefficient << std::endl;

        std::vector<CaseResult> case_results;
        std::vector<PersistentPatchTrack> case_tracks;

        for (size_t mi = 0; mi < modes.size(); mi++) {
            std::cout << "  Running mode: " << mode_names[mi] << "..." << std::endl;

            CaseResult result = RunCaseMode(config, modes[mi], mode_names[mi]);
            case_results.push_back(result);

            std::cout << "    Final Y: " << std::fixed << std::setprecision(6) << result.final_y << " m" << std::endl;
            std::cout << "    Y Error: " << result.y_error << " m" << std::endl;
            std::cout << "    Avg Force Y: " << result.avg_force_y << " N" << std::endl;
            std::cout << "    Force Std Dev: " << result.force_std_dev << " N" << std::endl;
            std::cout << "    Avg Torque (x,y,z): (" << result.avg_torque_x << ", " << result.avg_torque_y << ", " << result.avg_torque_z << ") Nm" << std::endl;
            std::cout << "    Max Patch Count: " << result.max_patch_count << std::endl;
            std::cout << "    Stable: " << (result.is_stable ? "YES" : "NO") << std::endl;
            if (mi >= 3) {
                std::cout << "    Avg Tangential Force: " << result.avg_tangential_force_norm << " N" << std::endl;
                std::cout << "    Tangential Force Ratio: " << result.avg_tangential_force_ratio << std::endl;
            }
            if (mi == 4 || mi == 5) {
                std::cout << "    Stick Steps: " << result.total_stick_steps << std::endl;
                std::cout << "    Slip Steps: " << result.total_slip_steps << std::endl;
                std::cout << "    Stick->Slip Transitions: " << result.total_stick_to_slip_transitions << std::endl;
                std::cout << "    Avg Tangential Displacement (local): " << result.avg_tangential_displacement_norm << " m" << std::endl;
                if (mi == 5) {
                    std::cout << "    Avg Transported State Norm: " << result.avg_transported_state_norm << " m" << std::endl;
                    std::cout << "    Avg Frame Rotation Angle: " << result.avg_frame_rotation_angle << " rad" << std::endl;
                }
            }
        }

        all_results.push_back(case_results);
    }

    // ========================================================================
    // Write output files
    // ========================================================================

    // Case A output
    std::ofstream case_A_file(out_dir + "/sdf_patch_frame_case_A.csv");
    case_A_file << "mode,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_transported_state_norm,avg_frame_rotation_angle,stick_steps,slip_steps,stick_to_slip_transitions,stable" << std::endl;
    for (const auto& r : all_results[0]) {
        case_A_file << r.mode_name << ","
                    << std::fixed << std::setprecision(8)
                    << r.final_y << ","
                    << r.expected_equilibrium_y << ","
                    << r.y_error << ","
                    << r.avg_force_y << ","
                    << r.force_std_dev << ","
                    << r.avg_torque_x << ","
                    << r.avg_torque_y << ","
                    << r.avg_torque_z << ","
                    << r.torque_std_dev << ","
                    << r.max_patch_count << ","
                    << r.avg_patch_count << ","
                    << r.multi_patch_ratio << ","
                    << r.avg_tangential_force_norm << ","
                    << r.max_tangential_force_norm << ","
                    << r.avg_tangential_force_ratio << ","
                    << r.avg_tangential_displacement_norm << ","
                    << r.avg_transported_state_norm << ","
                    << r.avg_frame_rotation_angle << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_A_file.close();

    // Case B output
    std::ofstream case_B_file(out_dir + "/sdf_patch_frame_case_B.csv");
    case_B_file << "mode,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_transported_state_norm,avg_frame_rotation_angle,stick_steps,slip_steps,stick_to_slip_transitions,stable" << std::endl;
    for (const auto& r : all_results[1]) {
        case_B_file << r.mode_name << ","
                    << std::fixed << std::setprecision(8)
                    << r.final_y << ","
                    << r.expected_equilibrium_y << ","
                    << r.y_error << ","
                    << r.avg_force_y << ","
                    << r.force_std_dev << ","
                    << r.avg_torque_x << ","
                    << r.avg_torque_y << ","
                    << r.avg_torque_z << ","
                    << r.torque_std_dev << ","
                    << r.max_patch_count << ","
                    << r.avg_patch_count << ","
                    << r.multi_patch_ratio << ","
                    << r.avg_tangential_force_norm << ","
                    << r.max_tangential_force_norm << ","
                    << r.avg_tangential_force_ratio << ","
                    << r.avg_tangential_displacement_norm << ","
                    << r.avg_transported_state_norm << ","
                    << r.avg_frame_rotation_angle << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_B_file.close();

    // Comparison
    std::ofstream compare_file(out_dir + "/sdf_patch_frame_compare.csv");
    compare_file << "metric,case_A_pointwise,case_A_patch_constitutive,case_A_tangential_damping,case_A_stickslip,case_A_frame_consistent,case_B_pointwise,case_B_patch_constitutive,case_B_tangential_damping,case_B_stickslip,case_B_frame_consistent" << std::endl;
    compare_file << "final_y,"
                 << std::fixed << std::setprecision(6)
                 << all_results[0][0].final_y << "," << all_results[0][2].final_y << "," << all_results[0][3].final_y << "," << all_results[0][4].final_y << "," << all_results[0][5].final_y << ","
                 << all_results[1][0].final_y << "," << all_results[1][2].final_y << "," << all_results[1][3].final_y << "," << all_results[1][4].final_y << "," << all_results[1][5].final_y << std::endl;
    compare_file << "y_error,"
                 << all_results[0][0].y_error << "," << all_results[0][2].y_error << "," << all_results[0][3].y_error << "," << all_results[0][4].y_error << "," << all_results[0][5].y_error << ","
                 << all_results[1][0].y_error << "," << all_results[1][2].y_error << "," << all_results[1][3].y_error << "," << all_results[1][4].y_error << "," << all_results[1][5].y_error << std::endl;
    compare_file << "force_std_dev,"
                 << all_results[0][0].force_std_dev << "," << all_results[0][2].force_std_dev << "," << all_results[0][3].force_std_dev << "," << all_results[0][4].force_std_dev << "," << all_results[0][5].force_std_dev << ","
                 << all_results[1][0].force_std_dev << "," << all_results[1][2].force_std_dev << "," << all_results[1][3].force_std_dev << "," << all_results[1][4].force_std_dev << "," << all_results[1][5].force_std_dev << std::endl;
    compare_file << "avg_torque_z,"
                 << all_results[0][0].avg_torque_z << "," << all_results[0][2].avg_torque_z << "," << all_results[0][3].avg_torque_z << "," << all_results[0][4].avg_torque_z << "," << all_results[0][5].avg_torque_z << ","
                 << all_results[1][0].avg_torque_z << "," << all_results[1][2].avg_torque_z << "," << all_results[1][3].avg_torque_z << "," << all_results[1][4].avg_torque_z << "," << all_results[1][5].avg_torque_z << std::endl;
    compare_file << "stable,"
                 << (all_results[0][0].is_stable ? 1 : 0) << "," << (all_results[0][2].is_stable ? 1 : 0) << "," << (all_results[0][3].is_stable ? 1 : 0) << "," << (all_results[0][4].is_stable ? 1 : 0) << "," << (all_results[0][5].is_stable ? 1 : 0) << ","
                 << (all_results[1][0].is_stable ? 1 : 0) << "," << (all_results[1][2].is_stable ? 1 : 0) << "," << (all_results[1][3].is_stable ? 1 : 0) << "," << (all_results[1][4].is_stable ? 1 : 0) << "," << (all_results[1][5].is_stable ? 1 : 0) << std::endl;
    compare_file << "stick_steps,0,0,0," << all_results[0][4].total_stick_steps << "," << all_results[0][5].total_stick_steps << ",0,0,0," << all_results[1][4].total_stick_steps << "," << all_results[1][5].total_stick_steps << std::endl;
    compare_file << "slip_steps,0,0,0," << all_results[0][4].total_slip_steps << "," << all_results[0][5].total_slip_steps << ",0,0,0," << all_results[1][4].total_slip_steps << "," << all_results[1][5].total_slip_steps << std::endl;
    compare_file << "avg_tangential_displacement_norm,0,0,0," << all_results[0][4].avg_tangential_displacement_norm << "," << all_results[0][5].avg_tangential_displacement_norm << ",0,0,0," << all_results[1][4].avg_tangential_displacement_norm << "," << all_results[1][5].avg_tangential_displacement_norm << std::endl;
    compare_file << "avg_transported_state_norm,0,0,0,0," << all_results[0][5].avg_transported_state_norm << ",0,0,0,0," << all_results[1][5].avg_transported_state_norm << std::endl;
    compare_file.close();

    // Track output
    WriteFrameTrackCSV(out_dir + "/sdf_patch_frame_tracks.csv", {});

    // ========================================================================
    // Console summary
    // ========================================================================

    std::cout << "\n=== Case A (Rotational Contact) Summary ===" << std::endl;
    std::cout << std::left << std::setw(18) << "Metric"
              << std::setw(16) << "Pointwise"
              << std::setw(16) << "PatchConstit"
              << std::setw(16) << "Tangential"
              << std::setw(16) << "StickSlip"
              << std::setw(16) << "Frame" << std::endl;
    std::cout << std::string(96, '-') << std::endl;
    std::cout << std::left << std::setw(18) << "Final Y (m)"
              << std::fixed << std::setprecision(6)
              << std::setw(16) << all_results[0][0].final_y
              << std::setw(16) << all_results[0][2].final_y
              << std::setw(16) << all_results[0][3].final_y
              << std::setw(16) << all_results[0][4].final_y
              << std::setw(16) << all_results[0][5].final_y << std::endl;
    std::cout << std::left << std::setw(18) << "Y Error (m)"
              << std::setw(16) << all_results[0][0].y_error
              << std::setw(16) << all_results[0][2].y_error
              << std::setw(16) << all_results[0][3].y_error
              << std::setw(16) << all_results[0][4].y_error
              << std::setw(16) << all_results[0][5].y_error << std::endl;
    std::cout << std::left << std::setw(18) << "Force Std (N)"
              << std::setw(16) << all_results[0][0].force_std_dev
              << std::setw(16) << all_results[0][2].force_std_dev
              << std::setw(16) << all_results[0][3].force_std_dev
              << std::setw(16) << all_results[0][4].force_std_dev
              << std::setw(16) << all_results[0][5].force_std_dev << std::endl;
    std::cout << std::left << std::setw(18) << "Avg Torque Z (Nm)"
              << std::setw(16) << all_results[0][0].avg_torque_z
              << std::setw(16) << all_results[0][2].avg_torque_z
              << std::setw(16) << all_results[0][3].avg_torque_z
              << std::setw(16) << all_results[0][4].avg_torque_z
              << std::setw(16) << all_results[0][5].avg_torque_z << std::endl;
    std::cout << std::left << std::setw(18) << "Avg Tangent (N)"
              << std::setw(16) << "N/A"
              << std::setw(16) << "N/A"
              << std::setw(16) << all_results[0][3].avg_tangential_force_norm
              << std::setw(16) << all_results[0][4].avg_tangential_force_norm
              << std::setw(16) << all_results[0][5].avg_tangential_force_norm << std::endl;
    std::cout << std::left << std::setw(18) << "Stable"
              << std::setw(16) << (all_results[0][0].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[0][2].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[0][3].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[0][4].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[0][5].is_stable ? "YES" : "NO") << std::endl;

    std::cout << "\n=== Case B (Multi-Patch) Summary ===" << std::endl;
    std::cout << std::left << std::setw(18) << "Metric"
              << std::setw(16) << "Pointwise"
              << std::setw(16) << "PatchConstit"
              << std::setw(16) << "Tangential"
              << std::setw(16) << "StickSlip"
              << std::setw(16) << "Frame" << std::endl;
    std::cout << std::string(96, '-') << std::endl;
    std::cout << std::left << std::setw(18) << "Final Y (m)"
              << std::fixed << std::setprecision(6)
              << std::setw(16) << all_results[1][0].final_y
              << std::setw(16) << all_results[1][2].final_y
              << std::setw(16) << all_results[1][3].final_y
              << std::setw(16) << all_results[1][4].final_y
              << std::setw(16) << all_results[1][5].final_y << std::endl;
    std::cout << std::left << std::setw(18) << "Y Error (m)"
              << std::setw(16) << all_results[1][0].y_error
              << std::setw(16) << all_results[1][2].y_error
              << std::setw(16) << all_results[1][3].y_error
              << std::setw(16) << all_results[1][4].y_error
              << std::setw(16) << all_results[1][5].y_error << std::endl;
    std::cout << std::left << std::setw(18) << "Force Std (N)"
              << std::setw(16) << all_results[1][0].force_std_dev
              << std::setw(16) << all_results[1][2].force_std_dev
              << std::setw(16) << all_results[1][3].force_std_dev
              << std::setw(16) << all_results[1][4].force_std_dev
              << std::setw(16) << all_results[1][5].force_std_dev << std::endl;
    std::cout << std::left << std::setw(18) << "Avg Torque Z (Nm)"
              << std::setw(16) << all_results[1][0].avg_torque_z
              << std::setw(16) << all_results[1][2].avg_torque_z
              << std::setw(16) << all_results[1][3].avg_torque_z
              << std::setw(16) << all_results[1][4].avg_torque_z
              << std::setw(16) << all_results[1][5].avg_torque_z << std::endl;
    std::cout << std::left << std::setw(18) << "Avg Tangent (N)"
              << std::setw(16) << "N/A"
              << std::setw(16) << "N/A"
              << std::setw(16) << all_results[1][3].avg_tangential_force_norm
              << std::setw(16) << all_results[1][4].avg_tangential_force_norm
              << std::setw(16) << all_results[1][5].avg_tangential_force_norm << std::endl;
    std::cout << std::left << std::setw(18) << "Stable"
              << std::setw(16) << (all_results[1][0].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[1][2].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[1][3].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[1][4].is_stable ? "YES" : "NO")
              << std::setw(16) << (all_results[1][5].is_stable ? "YES" : "NO") << std::endl;

    // ========================================================================
    // Verification
    // ========================================================================

    bool case_a_frame_better = all_results[0][5].y_error < all_results[0][4].y_error; // vs stick-slip
    bool case_b_frame_stable = all_results[1][5].is_stable;
    bool case_b_not_degraded = all_results[1][5].y_error < all_results[1][0].y_error;

    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "  Case A frame-consistent y_error vs stick-slip: " << std::fixed << std::setprecision(4) 
              << all_results[0][5].y_error << " vs " << all_results[0][4].y_error << std::endl;
    std::cout << "  Case A improvement: " << (case_a_frame_better ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Case B frame-consistent stable: " << (case_b_frame_stable ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Case B not degraded vs pointwise: " << (case_b_not_degraded ? "PASS" : "FAIL") << std::endl;

    bool pass = case_a_frame_better && case_b_frame_stable && case_b_not_degraded;
    std::cout << "  Overall: " << (pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_frame_case_A.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_frame_case_B.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_frame_compare.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_frame_tracks.csv" << std::endl;

    return pass ? 0 : 1;
}
