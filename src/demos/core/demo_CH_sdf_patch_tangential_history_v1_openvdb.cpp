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
// Milestone 14: Tangential History V1 - Geometry-Consistent Transfer
//
// This demo replaces fixed carry-over with v1 geometry-consistent pipeline:
//   1. History state bound to persistent patch (world-coordinate xi_elastic)
//   2. Project history onto current tangent plane
//   3. Compute adaptive alpha from geometry/match quality/topology
//   4. xi_transfer = alpha * projected_history
//   5. Add current-step tangential increment -> xi_trial
//   6. Coulomb return mapping: stick if within limit, slip otherwise
//
// Scene: Case A (rotational contact) and Case B (multi-patch)
//
// Output:
//   out/milestone_14/sdf_patch_tangential_history_v1_case_A.csv
//   out/milestone_14/sdf_patch_tangential_history_v1_case_B.csv
//   out/milestone_14/sdf_patch_tangential_history_v1_summary.csv
//   out/milestone_14/sdf_patch_tangential_history_v1_diagnostics.csv
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
    int matched_persistent_id;
    std::vector<int> sample_ids;
    int active_sample_count;
    ChVector3d center;
    ChVector3d normal;
    ChVector3d tangent_t1;
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
    double patch_radius_estimate;
    double tangential_displacement_local_norm;
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
// Persistent patch track
// =============================================================================

struct PersistentPatchTrack {
    int persistent_id;
    int current_frame_patch_id;
    TrackStatus status;
    ChVector3d center;
    ChVector3d normal;
    ChVector3d tangent_t1;
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
    double patch_radius_estimate;
    double xi1;
    double xi2;
    ChVector3d tangential_displacement;
    int stick_steps;
    int slip_steps;
    int stick_to_slip_transitions;
    double accumulated_transport_norm;
    int total_transport_clamp_hits;
    double sum_attenuation;
    double sum_carry_over;
    int topology_birth_event;
};

// =============================================================================
// Tangential history state bound to persistent patch (V1)
// =============================================================================

struct PatchTangentialState {
    int persistent_id;
    ChVector3d xi_elastic_world;
    ChVector3d normal_prev;
    ChVector3d t1_prev;
    ChVector3d t2_prev;
    double match_quality_prev;
    double overlap_prev;
    double area_prev;
    int age;
    bool valid;
};

// =============================================================================
// Carry-over strategy
// =============================================================================

enum class CarryOverStrategy {
    Fixed,
    AdaptiveV1
};

// =============================================================================
// Transport configuration
// =============================================================================

struct TransportConfig {
    std::string name;
    bool enable_clamp;
    bool enable_attenuation;
    CarryOverStrategy carry_strategy;
    double xi_max;
    double attenuation_beta;
    double fixed_carry_alpha;
    double adaptive_alpha_base;
    double adaptive_theta_scale;
    double adaptive_match_scale;
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
    double tangential_damping;
    double tangential_stiffness;
    double friction_coefficient;
    double ema_alpha;
};

// =============================================================================
// Result metrics
// =============================================================================

struct RunResult {
    std::string config_name;
    bool enable_clamp;
    bool enable_attenuation;
    CarryOverStrategy carry_strategy;
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
    double avg_raw_transport_norm;
    double avg_limited_transport_norm;
    double avg_transport_attenuation;
    double avg_carry_over_factor;
    double min_carry_over_factor;
    double max_carry_over_factor;
    double avg_adaptive_alpha;
    double min_adaptive_alpha;
    double max_adaptive_alpha;
    int total_transport_clamp_hits;
    double avg_theta;
    double max_theta;
    double avg_match_distance;
    int total_birth_events;
    int total_reset_events;
    double avg_overlap;
    double min_overlap;
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
    ChVector3d ref = std::abs(n.x()) < 0.9 ? ChVector3d(1, 0, 0) : ChVector3d(0, 1, 0);
    t1 = (n.Cross(ref)).GetNormalized();
    t2 = n.Cross(t1);
}

// =============================================================================
// Project vector onto tangent plane (remove normal component)
// =============================================================================

ChVector3d ProjectToTangentPlane(const ChVector3d& vec_world, const ChVector3d& normal_current) {
    double normal_component = vec_world.Dot(normal_current);
    return vec_world - normal_component * normal_current;
}

// =============================================================================
// Compute adaptive alpha for V1 history transfer
// =============================================================================

double ComputeAdaptiveAlphaV1(
    double theta,
    double match_distance,
    bool is_newborn,
    const TransportConfig& config
) {
    double alpha = config.adaptive_alpha_base;

    double theta_factor = std::exp(-theta / config.adaptive_theta_scale);
    alpha *= theta_factor;

    double match_factor = std::exp(-match_distance / config.adaptive_match_scale);
    alpha *= match_factor;

    if (is_newborn) {
        alpha = 0.0;
    }

    alpha = std::max(0.0, std::min(1.0, alpha));

    return alpha;
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
        patch.matched_persistent_id = -1;
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
        patch.patch_radius_estimate = 0.0;
        patch.tangential_displacement_local_norm = 0.0;

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

        patch.patch_radius_estimate = std::sqrt(patch.total_area / M_PI);
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
    double match_distance;
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
        m.match_distance = distance_threshold + 1.0;
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
            m.match_distance = best_dist;
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
// Run single case with single transport config
// =============================================================================

RunResult RunCaseWithConfig(
    const CaseConfig& config,
    const TransportConfig& transport_cfg,
    const std::string& mode_name
) {
    openvdb::initialize();

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

    std::vector<SurfaceSample> samples;
    GridAdjacency adjacency;
    if (config.dyn_body_type == 0) {
        samples = GenerateSphereSamples(config.dyn_body_radius, config.n_theta, config.n_phi);
        adjacency = GridAdjacency(config.n_theta, config.n_phi);
    } else {
        samples = GenerateBoxSamples(config.dyn_body_size, 100);
        adjacency = GridAdjacency(static_cast<int>(samples.size()));
    }

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

    std::vector<PersistentPatchTrack> active_tracks;
    int next_persistent_id = 0;

    // V1 tangential history state map (keyed by persistent_id)
    std::map<int, PatchTangentialState> tangential_history;

    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());

    std::vector<double> force_y_values;
    std::vector<double> pos_y_values;
    std::vector<int> patch_counts;
    std::vector<double> torque_x_values;
    std::vector<double> torque_y_values;
    std::vector<double> torque_z_values;
    std::vector<double> tangential_force_norms;
    std::vector<double> tangential_displacement_norms;
    std::vector<double> raw_transport_norms;
    std::vector<double> limited_transport_norms;
    std::vector<double> attenuation_values;
    std::vector<double> carry_values;
    std::vector<double> adaptive_alpha_values;
    std::vector<double> theta_values;
    std::vector<double> match_distance_values;
    std::vector<double> overlap_values;
    int total_birth_events = 0;
    int total_reset_events = 0;
    int total_transport_clamp_hits = 0;
    double total_force_y = 0.0;
    int force_sample_count = 0;
    double sum_tangential_norm = 0.0;
    double max_tangential_norm = 0.0;
    double sum_tangential_ratio = 0.0;
    int total_stick_steps = 0;
    int total_slip_steps = 0;
    int total_stick_to_slip_transitions = 0;

    double match_distance_threshold = 0.5 * (config.dyn_body_type == 0 ? config.dyn_body_radius : config.dyn_body_size);
    double normal_similarity_threshold = 0.9;

    int frame_offset = 0;

    while (sys.GetChTime() < config.total_time) {
        dyn_body->EmptyAccumulator(acc_id);

        ChVector3d body_pos = dyn_body->GetPos();
        ChQuaterniond body_rot = dyn_body->GetRot();
        ChVector3d body_vel = dyn_body->GetPosDt();
        ChVector3d body_ang_vel_world = body_rot.Rotate(dyn_body->GetAngVelLocal());

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

            double phi = sdf_results[si].phi;
            if (phi < config.activation_band) {
                active_indices.push_back(static_cast<int>(si));
            }
        }

        int frame_patch_id_offset = frame_offset * 1000;
        frame_offset++;
        std::vector<PatchPrimitive> patches = BuildPatches(
            active_indices, samples, sdf_results,
            sample_world_positions, adjacency, frame_patch_id_offset
        );

        std::vector<MatchResult> matches = MatchPatches(
            patches, active_tracks,
            match_distance_threshold, normal_similarity_threshold
        );

        std::vector<PersistentPatchTrack> current_tracks;

        for (const auto& m : matches) {
            PersistentPatchTrack track;
            track.topology_birth_event = 0;
            if (m.matched) {
                for (const auto& prev : active_tracks) {
                    if (prev.persistent_id == m.persistent_id) {
                        track = prev;
                        track.status = TrackStatus::Alive;
                        track.topology_birth_event = 0;
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
                track.total_transport_clamp_hits = 0;
                track.sum_attenuation = 0.0;
                track.sum_carry_over = 0.0;
                track.topology_birth_event = 1;
                total_birth_events++;
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
                    track.patch_radius_estimate = p.patch_radius_estimate;
                    for (auto& qp : patches) {
                        if (qp.patch_id == track.current_frame_patch_id) {
                            qp.matched_persistent_id = track.persistent_id;
                            break;
                        }
                    }
                    break;
                }
            }

            track.last_seen_time = sys.GetChTime();
            if (track.status == TrackStatus::Alive) {
                track.age_in_steps++;
            }

            current_tracks.push_back(track);
        }

        // Compute patch forces with V1 tangential history transfer
        for (size_t pi = 0; pi < patches.size(); pi++) {
            auto& patch = patches[pi];
            ChVector3d patch_center_vel = body_vel + body_ang_vel_world.Cross(patch.center - body_pos);

            const PersistentPatchTrack* prev_track = nullptr;
            if (patch.matched_persistent_id >= 0) {
                for (const auto& t : active_tracks) {
                    if (t.persistent_id == patch.matched_persistent_id) {
                        prev_track = &t;
                        break;
                    }
                }
            }

            // --- Normal force (with EMA smoothing) ---
            double penetration_current = patch.effective_penetration;
            double penetration_eff = penetration_current;
            if (prev_track != nullptr) {
                penetration_eff = 0.3 * penetration_current + 0.7 * prev_track->effective_penetration;
            }

            double vn_eff = patch_center_vel.Dot(patch.normal);
            double force_n_scalar = config.stiffness * penetration_eff + config.damping * std::max(-vn_eff, 0.0);
            if (force_n_scalar < 0.0) force_n_scalar = 0.0;
            ChVector3d F_n = patch.normal * force_n_scalar;

            // --- Tangential velocity ---
            double vn = patch_center_vel.Dot(patch.normal);
            ChVector3d vt = patch_center_vel - vn * patch.normal;

            double vt1 = vt.Dot(patch.tangent_t1);
            double vt2 = vt.Dot(patch.tangent_t2);

            // --- V1 Tangential History Transfer ---
            double xi1 = 0.0, xi2 = 0.0;
            double raw_trans_norm = 0.0, limited_trans_norm = 0.0;
            double attenuation = 1.0, carry_factor = 1.0;
            double adaptive_alpha = 1.0;
            double frame_rot_angle = 0.0;
            double match_dist = 0.0;
            bool is_newborn = (prev_track == nullptr || prev_track->topology_birth_event == 1);

            if (prev_track != nullptr && tangential_history.count(prev_track->persistent_id)) {
                auto& hist = tangential_history[prev_track->persistent_id];
                if (hist.valid) {
                    frame_rot_angle = std::acos(std::max(-1.0, std::min(1.0, hist.normal_prev.Dot(patch.normal))));
                    match_dist = (patch.center - prev_track->center).Length();

                    // Step 3: Project history onto current tangent plane
                    ChVector3d xi_proj = ProjectToTangentPlane(hist.xi_elastic_world, patch.normal);

                    raw_trans_norm = xi_proj.Length();

                    // Step 4: Compute adaptive alpha
                    adaptive_alpha = ComputeAdaptiveAlphaV1(
                        frame_rot_angle,
                        match_dist,
                        is_newborn,
                        transport_cfg
                    );

                    // Step 5: xi_transfer = alpha * projected_history
                    ChVector3d xi_transfer = xi_proj * adaptive_alpha;

                    // Step 6: Add current step increment
                    double delta_xi1 = vt1 * config.time_step;
                    double delta_xi2 = vt2 * config.time_step;
                    ChVector3d delta_xi = delta_xi1 * patch.tangent_t1 + delta_xi2 * patch.tangent_t2;
                    ChVector3d xi_trial = xi_transfer + delta_xi;

                    // Step 7: Coulomb return mapping
                    double fn_mag = F_n.Length();
                    double friction_limit = config.friction_coefficient * fn_mag;

                    // Trial tangential force
                    ChVector3d tau_trial = config.tangential_stiffness * xi_trial + config.tangential_damping * vt;
                    double tau_trial_mag = tau_trial.Length();

                    if (tau_trial_mag <= friction_limit) {
                        // Stick
                        xi1 = xi_trial.Dot(patch.tangent_t1);
                        xi2 = xi_trial.Dot(patch.tangent_t2);
                        patch.stick_slip_state = StickSlipState::Stick;
                    } else {
                        // Slip: project back to friction limit
                        ChVector3d tau_new = tau_trial * (friction_limit / tau_trial_mag);
                        ChVector3d xi_elastic_new = tau_new / config.tangential_stiffness;
                        xi1 = xi_elastic_new.Dot(patch.tangent_t1);
                        xi2 = xi_elastic_new.Dot(patch.tangent_t2);
                        patch.stick_slip_state = StickSlipState::Slip;
                    }

                    limited_trans_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
                    carry_factor = adaptive_alpha;
                } else {
                    // Invalid history: reset
                    total_reset_events++;
                    xi1 = vt1 * config.time_step;
                    xi2 = vt2 * config.time_step;
                    patch.stick_slip_state = StickSlipState::Stick;
                    limited_trans_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
                }
            } else {
                // No history: newborn or first contact
                total_reset_events++;
                xi1 = vt1 * config.time_step;
                xi2 = vt2 * config.time_step;
                patch.stick_slip_state = StickSlipState::Stick;
                limited_trans_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
            }

            // Compute tangential force from xi
            ChVector3d F_t = -config.tangential_stiffness * (xi1 * patch.tangent_t1 + xi2 * patch.tangent_t2)
                           - config.tangential_damping * vt;

            // Apply Coulomb clamp
            double ft_mag = F_t.Length();
            double fn_mag = F_n.Length();
            double friction_limit = config.friction_coefficient * fn_mag;
            if (ft_mag > friction_limit && fn_mag > 1e-12) {
                F_t = F_t * (friction_limit / ft_mag);
                ft_mag = friction_limit;
            }

            patch.force = F_n + F_t;
            patch.normal_force_magnitude = fn_mag;
            patch.tangential_force_magnitude = ft_mag;
            patch.torque = (patch.center - body_pos).Cross(patch.force);
            patch.tangential_displacement_local_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);

            // Update tangential history state
            ChVector3d xi_elastic_world_new = xi1 * patch.tangent_t1 + xi2 * patch.tangent_t2;
            if (patch.matched_persistent_id >= 0) {
                PatchTangentialState new_hist;
                new_hist.persistent_id = patch.matched_persistent_id;
                new_hist.xi_elastic_world = xi_elastic_world_new;
                new_hist.normal_prev = patch.normal;
                new_hist.t1_prev = patch.tangent_t1;
                new_hist.t2_prev = patch.tangent_t2;
                new_hist.match_quality_prev = std::exp(-match_dist / config.dyn_body_radius);
                new_hist.overlap_prev = patch.effective_penetration;
                new_hist.area_prev = patch.total_area;
                new_hist.age = (prev_track != nullptr) ? prev_track->age_in_steps : 0;
                new_hist.valid = true;
                tangential_history[patch.matched_persistent_id] = new_hist;
            }

            // Update track state
            for (auto& track : current_tracks) {
                if (track.persistent_id == patch.matched_persistent_id) {
                    track.xi1 = xi1;
                    track.xi2 = xi2;
                    track.stick_slip_state = patch.stick_slip_state;
                    track.normal_force_magnitude = fn_mag;
                    track.tangential_force_magnitude = ft_mag;
                    if (patch.stick_slip_state == StickSlipState::Stick) track.stick_steps++;
                    else track.slip_steps++;
                    if (prev_track != nullptr && prev_track->stick_slip_state == StickSlipState::Stick && patch.stick_slip_state == StickSlipState::Slip) {
                        track.stick_to_slip_transitions++;
                    }
                    break;
                }
            }

            total_force += patch.force;
            total_torque += patch.torque;

            // Record diagnostics
            theta_values.push_back(frame_rot_angle);
            match_distance_values.push_back(match_dist);
            overlap_values.push_back(patch.effective_penetration);
            if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV1) {
                adaptive_alpha_values.push_back(adaptive_alpha);
            }
        }

        active_tracks = current_tracks;

        if (total_force.Length() > 1e-12) {
            dyn_body->AccumulateForce(acc_id, total_force, body_pos, false);
            dyn_body->AccumulateTorque(acc_id, body_rot.GetConjugate().Rotate(total_torque), false);
        }

        sys.DoStepDynamics(config.time_step);

        double current_force_y = total_force.y();
        force_y_values.push_back(current_force_y);
        pos_y_values.push_back(dyn_body->GetPos().y());
        int current_patch_count = patches.size();
        patch_counts.push_back(current_patch_count);
        torque_x_values.push_back(total_torque.x());
        torque_y_values.push_back(total_torque.y());
        torque_z_values.push_back(total_torque.z());

        total_force_y += current_force_y;
        force_sample_count++;

        double ft_approx = std::sqrt(total_force.x() * total_force.x() + total_force.z() * total_force.z());
        double fn_approx = std::abs(total_force.y());
        tangential_force_norms.push_back(ft_approx);
        sum_tangential_norm += ft_approx;
        if (ft_approx > max_tangential_norm) max_tangential_norm = ft_approx;
        if (fn_approx > 1e-12) {
            sum_tangential_ratio += ft_approx / fn_approx;
        }

        for (const auto& track : active_tracks) {
            total_stick_steps += track.stick_steps;
            total_slip_steps += track.slip_steps;
            total_stick_to_slip_transitions += track.stick_to_slip_transitions;
            double xi_norm = std::sqrt(track.xi1 * track.xi1 + track.xi2 * track.xi2);
            tangential_displacement_norms.push_back(xi_norm);
            raw_transport_norms.push_back(track.accumulated_transport_norm);
            limited_transport_norms.push_back(track.accumulated_transport_norm);
            attenuation_values.push_back(track.sum_attenuation / std::max(1, track.age_in_steps));
            carry_values.push_back(track.sum_carry_over / std::max(1, track.age_in_steps));
        }
        total_transport_clamp_hits += active_tracks.empty() ? 0 : active_tracks.back().total_transport_clamp_hits;
    }

    // Compute statistics
    RunResult result;
    result.config_name = mode_name;
    result.enable_clamp = transport_cfg.enable_clamp;
    result.enable_attenuation = transport_cfg.enable_attenuation;
    result.carry_strategy = transport_cfg.carry_strategy;
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
    result.avg_raw_transport_norm = raw_transport_norms.empty() ? 0.0 :
        std::accumulate(raw_transport_norms.begin(), raw_transport_norms.end(), 0.0) / raw_transport_norms.size();
    result.avg_limited_transport_norm = limited_transport_norms.empty() ? 0.0 :
        std::accumulate(limited_transport_norms.begin(), limited_transport_norms.end(), 0.0) / limited_transport_norms.size();
    result.avg_transport_attenuation = attenuation_values.empty() ? 0.0 :
        std::accumulate(attenuation_values.begin(), attenuation_values.end(), 0.0) / attenuation_values.size();
    result.avg_carry_over_factor = carry_values.empty() ? 0.0 :
        std::accumulate(carry_values.begin(), carry_values.end(), 0.0) / carry_values.size();
    result.min_carry_over_factor = carry_values.empty() ? 0.0 : *std::min_element(carry_values.begin(), carry_values.end());
    result.max_carry_over_factor = carry_values.empty() ? 0.0 : *std::max_element(carry_values.begin(), carry_values.end());

    result.avg_adaptive_alpha = adaptive_alpha_values.empty() ? 0.0 :
        std::accumulate(adaptive_alpha_values.begin(), adaptive_alpha_values.end(), 0.0) / adaptive_alpha_values.size();
    result.min_adaptive_alpha = adaptive_alpha_values.empty() ? 0.0 : *std::min_element(adaptive_alpha_values.begin(), adaptive_alpha_values.end());
    result.max_adaptive_alpha = adaptive_alpha_values.empty() ? 0.0 : *std::max_element(adaptive_alpha_values.begin(), adaptive_alpha_values.end());

    result.avg_theta = theta_values.empty() ? 0.0 :
        std::accumulate(theta_values.begin(), theta_values.end(), 0.0) / theta_values.size();
    result.max_theta = theta_values.empty() ? 0.0 : *std::max_element(theta_values.begin(), theta_values.end());
    result.avg_match_distance = match_distance_values.empty() ? 0.0 :
        std::accumulate(match_distance_values.begin(), match_distance_values.end(), 0.0) / match_distance_values.size();
    result.total_birth_events = total_birth_events;
    result.total_reset_events = total_reset_events;

    result.avg_overlap = overlap_values.empty() ? 0.0 :
        std::accumulate(overlap_values.begin(), overlap_values.end(), 0.0) / overlap_values.size();
    result.min_overlap = overlap_values.empty() ? 0.0 : *std::min_element(overlap_values.begin(), overlap_values.end());

    result.total_stick_steps = total_stick_steps;
    result.total_slip_steps = total_slip_steps;
    result.total_stick_to_slip_transitions = total_stick_to_slip_transitions;
    result.total_transport_clamp_hits = total_transport_clamp_hits;

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
// Main simulation
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 14: Tangential History V1 ===" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_14";
    EnsureDir(out_dir);

    // ========================================================================
    // Baseline Reference
    // ========================================================================
    std::cout << "\n=== Baseline Reference (from Milestone 13) ===" << std::endl;
    std::cout << "Case A:" << std::endl;
    std::cout << "  A-only: y_error=0.380m, torque_z=+1.45Nm, trans_norm=3.56m, stable=NO" << std::endl;
    std::cout << "  C-fixed: y_error=0.075m, torque_z=+0.53Nm, trans_norm=0.12m, stable=NO" << std::endl;
    std::cout << "Case B:" << std::endl;
    std::cout << "  A-only: y_error=0.031m, torque_z=-0.025Nm, trans_norm=2.01m, stable=YES" << std::endl;
    std::cout << "  C-fixed: y_error=0.074m, torque_z=-0.121Nm, trans_norm=0.13m, stable=YES" << std::endl;

    // ========================================================================
    // Case Definitions
    // ========================================================================

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

    // ========================================================================
    // Experimental Matrix
    // ========================================================================

    double adaptive_alpha_base = 0.8;
    double adaptive_theta_scale = 0.3;
    double adaptive_match_scale = 0.05;

    TransportConfig cfg_A = {"A-only", true, false, CarryOverStrategy::Fixed, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale};
    TransportConfig cfg_C_fixed = {"C-fixed", false, false, CarryOverStrategy::Fixed, 0.01, 5.0, 0.5, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale};
    TransportConfig cfg_C_adaptive_v1 = {"C-adaptive-v1", false, false, CarryOverStrategy::AdaptiveV1, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale};
    TransportConfig cfg_A_C_adaptive_v1 = {"A+C-adaptive-v1", true, false, CarryOverStrategy::AdaptiveV1, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale};
    TransportConfig cfg_A_C_fixed = {"A+C-fixed", true, false, CarryOverStrategy::Fixed, 0.01, 5.0, 0.5, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale};

    std::vector<TransportConfig> configs = {cfg_A, cfg_C_fixed, cfg_C_adaptive_v1, cfg_A_C_adaptive_v1, cfg_A_C_fixed};

    // ========================================================================
    // Run all cases
    // ========================================================================

    std::vector<std::vector<RunResult>> all_results;

    for (const auto& config : cases) {
        std::cout << "\n=== Running " << config.name << " ===" << std::endl;
        std::cout << "Description: " << config.description << std::endl;

        std::vector<RunResult> case_results;

        for (size_t ci = 0; ci < configs.size(); ci++) {
            std::cout << "  Running config: " << configs[ci].name << "..." << std::endl;

            RunResult result = RunCaseWithConfig(config, configs[ci], configs[ci].name);
            case_results.push_back(result);

            std::cout << "    Final Y: " << std::fixed << std::setprecision(6) << result.final_y << " m" << std::endl;
            std::cout << "    Y Error: " << result.y_error << " m" << std::endl;
            std::cout << "    Force Std Dev: " << result.force_std_dev << " N" << std::endl;
            std::cout << "    Avg Torque Z: " << result.avg_torque_z << " Nm" << std::endl;
            std::cout << "    Avg Adaptive Alpha: " << result.avg_adaptive_alpha << std::endl;
            std::cout << "    Min/Max Adaptive Alpha: " << result.min_adaptive_alpha << " / " << result.max_adaptive_alpha << std::endl;
            std::cout << "    Avg Theta: " << result.avg_theta << " rad" << std::endl;
            std::cout << "    Max Theta: " << result.max_theta << " rad" << std::endl;
            std::cout << "    Total Reset Events: " << result.total_reset_events << std::endl;
            std::cout << "    Stable: " << (result.is_stable ? "YES" : "NO") << std::endl;
        }

        all_results.push_back(case_results);
    }

    // ========================================================================
    // Write output files
    // ========================================================================

    std::ofstream case_A_file(out_dir + "/sdf_patch_tangential_history_v1_case_A.csv");
    case_A_file << "config,clamp,attenuation,carry_strategy,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_raw_transport_norm,avg_limited_transport_norm,avg_transport_attenuation,avg_carry_over_factor,min_carry_over_factor,max_carry_over_factor,avg_adaptive_alpha,min_adaptive_alpha,max_adaptive_alpha,avg_theta,max_theta,avg_match_distance,total_birth_events,total_reset_events,total_transport_clamp_hits,stick_steps,slip_steps,stick_to_slip_transitions,avg_overlap,min_overlap,stable" << std::endl;
    for (const auto& r : all_results[0]) {
        std::string carry_str = (r.carry_strategy == CarryOverStrategy::AdaptiveV1) ? "adaptive-v1" : "fixed";
        case_A_file << r.config_name << ","
                    << (r.enable_clamp ? 1 : 0) << ","
                    << (r.enable_attenuation ? 1 : 0) << ","
                    << carry_str << ","
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
                    << r.avg_raw_transport_norm << ","
                    << r.avg_limited_transport_norm << ","
                    << r.avg_transport_attenuation << ","
                    << r.avg_carry_over_factor << ","
                    << r.min_carry_over_factor << ","
                    << r.max_carry_over_factor << ","
                    << r.avg_adaptive_alpha << ","
                    << r.min_adaptive_alpha << ","
                    << r.max_adaptive_alpha << ","
                    << r.avg_theta << ","
                    << r.max_theta << ","
                    << r.avg_match_distance << ","
                    << r.total_birth_events << ","
                    << r.total_reset_events << ","
                    << r.total_transport_clamp_hits << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << r.avg_overlap << ","
                    << r.min_overlap << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_A_file.close();

    std::ofstream case_B_file(out_dir + "/sdf_patch_tangential_history_v1_case_B.csv");
    case_B_file << "config,clamp,attenuation,carry_strategy,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_raw_transport_norm,avg_limited_transport_norm,avg_transport_attenuation,avg_carry_over_factor,min_carry_over_factor,max_carry_over_factor,avg_adaptive_alpha,min_adaptive_alpha,max_adaptive_alpha,avg_theta,max_theta,avg_match_distance,total_birth_events,total_reset_events,total_transport_clamp_hits,stick_steps,slip_steps,stick_to_slip_transitions,avg_overlap,min_overlap,stable" << std::endl;
    for (const auto& r : all_results[1]) {
        std::string carry_str = (r.carry_strategy == CarryOverStrategy::AdaptiveV1) ? "adaptive-v1" : "fixed";
        case_B_file << r.config_name << ","
                    << (r.enable_clamp ? 1 : 0) << ","
                    << (r.enable_attenuation ? 1 : 0) << ","
                    << carry_str << ","
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
                    << r.avg_raw_transport_norm << ","
                    << r.avg_limited_transport_norm << ","
                    << r.avg_transport_attenuation << ","
                    << r.avg_carry_over_factor << ","
                    << r.min_carry_over_factor << ","
                    << r.max_carry_over_factor << ","
                    << r.avg_adaptive_alpha << ","
                    << r.min_adaptive_alpha << ","
                    << r.max_adaptive_alpha << ","
                    << r.avg_theta << ","
                    << r.max_theta << ","
                    << r.avg_match_distance << ","
                    << r.total_birth_events << ","
                    << r.total_reset_events << ","
                    << r.total_transport_clamp_hits << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << r.avg_overlap << ","
                    << r.min_overlap << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_B_file.close();

    // Summary
    std::ofstream summary_file(out_dir + "/sdf_patch_tangential_history_v1_summary.csv");
    summary_file << "config,clamp,attenuation,carry_strategy,case_A_y_error,case_A_torque_z,case_A_trans_norm,case_A_alpha_mean,case_A_alpha_min,case_A_stable,case_B_y_error,case_B_torque_z,case_B_trans_norm,case_B_alpha_mean,case_B_alpha_min,case_B_stable" << std::endl;
    for (size_t ci = 0; ci < configs.size(); ci++) {
        std::string carry_str = (configs[ci].carry_strategy == CarryOverStrategy::AdaptiveV1) ? "adaptive-v1" : "fixed";
        summary_file << configs[ci].name << ","
                     << (configs[ci].enable_clamp ? 1 : 0) << ","
                     << (configs[ci].enable_attenuation ? 1 : 0) << ","
                     << carry_str << ","
                     << std::fixed << std::setprecision(6)
                     << all_results[0][ci].y_error << ","
                     << all_results[0][ci].avg_torque_z << ","
                     << all_results[0][ci].avg_limited_transport_norm << ","
                     << all_results[0][ci].avg_adaptive_alpha << ","
                     << all_results[0][ci].min_adaptive_alpha << ","
                     << (all_results[0][ci].is_stable ? 1 : 0) << ","
                     << all_results[1][ci].y_error << ","
                     << all_results[1][ci].avg_torque_z << ","
                     << all_results[1][ci].avg_limited_transport_norm << ","
                     << all_results[1][ci].avg_adaptive_alpha << ","
                     << all_results[1][ci].min_adaptive_alpha << ","
                     << (all_results[1][ci].is_stable ? 1 : 0) << std::endl;
    }
    summary_file.close();

    // Diagnostics
    std::ofstream diag_file(out_dir + "/sdf_patch_tangential_history_v1_diagnostics.csv");
    diag_file << "config,case_A_avg_theta,case_A_max_theta,case_A_avg_match_dist,case_A_birth_events,case_A_reset_events,case_A_avg_alpha,case_A_min_alpha,case_A_max_alpha,case_B_avg_theta,case_B_max_theta,case_B_avg_match_dist,case_B_birth_events,case_B_reset_events,case_B_avg_alpha,case_B_min_alpha,case_B_max_alpha" << std::endl;
    for (size_t ci = 0; ci < configs.size(); ci++) {
        diag_file << configs[ci].name << ","
                  << std::fixed << std::setprecision(6)
                  << all_results[0][ci].avg_theta << "," << all_results[0][ci].max_theta << "," << all_results[0][ci].avg_match_distance << "," << all_results[0][ci].total_birth_events << "," << all_results[0][ci].total_reset_events << ","
                  << all_results[0][ci].avg_adaptive_alpha << "," << all_results[0][ci].min_adaptive_alpha << "," << all_results[0][ci].max_adaptive_alpha << ","
                  << all_results[1][ci].avg_theta << "," << all_results[1][ci].max_theta << "," << all_results[1][ci].avg_match_distance << "," << all_results[1][ci].total_birth_events << "," << all_results[1][ci].total_reset_events << ","
                  << all_results[1][ci].avg_adaptive_alpha << "," << all_results[1][ci].min_adaptive_alpha << "," << all_results[1][ci].max_adaptive_alpha << std::endl;
    }
    diag_file.close();

    // ========================================================================
    // Console summary
    // ========================================================================

    std::cout << "\n=== Case A Tangential History V1 Summary ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Config"
              << std::setw(12) << "y_error"
              << std::setw(12) << "torque_z"
              << std::setw(12) << "trans_norm"
              << std::setw(12) << "alpha_avg"
              << std::setw(12) << "alpha_min"
              << std::setw(8) << "stable" << std::endl;
    std::cout << std::string(76, '-') << std::endl;
    for (const auto& r : all_results[0]) {
        std::cout << std::left << std::setw(20) << r.config_name
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << r.y_error
                  << std::setw(12) << r.avg_torque_z
                  << std::setw(12) << r.avg_limited_transport_norm
                  << std::setw(12) << r.avg_adaptive_alpha
                  << std::setw(12) << r.min_adaptive_alpha
                  << std::setw(8) << (r.is_stable ? "YES" : "NO") << std::endl;
    }

    std::cout << "\n=== Case B Tangential History V1 Summary ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Config"
              << std::setw(12) << "y_error"
              << std::setw(12) << "torque_z"
              << std::setw(12) << "trans_norm"
              << std::setw(12) << "alpha_avg"
              << std::setw(12) << "alpha_min"
              << std::setw(8) << "stable" << std::endl;
    std::cout << std::string(76, '-') << std::endl;
    for (const auto& r : all_results[1]) {
        std::cout << std::left << std::setw(20) << r.config_name
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << r.y_error
                  << std::setw(12) << r.avg_torque_z
                  << std::setw(12) << r.avg_limited_transport_norm
                  << std::setw(12) << r.avg_adaptive_alpha
                  << std::setw(12) << r.min_adaptive_alpha
                  << std::setw(8) << (r.is_stable ? "YES" : "NO") << std::endl;
    }

    // ========================================================================
    // Verification
    // ========================================================================

    bool case_a_v1_vs_fixed = all_results[0][2].y_error <= all_results[0][1].y_error + 0.01;
    bool case_b_v1_stable = all_results[1][2].is_stable;
    bool case_b_v1_not_degraded = all_results[1][2].y_error < all_results[1][0].y_error + 0.01;

    std::cout << "\n=== Verification ===" << std::endl;
    std::cout << "  Case A C-adaptive-v1 vs C-fixed: " << std::fixed << std::setprecision(4) 
              << all_results[0][2].y_error << " vs " << all_results[0][1].y_error << std::endl;
    std::cout << "  Case A V1 improvement (or not worse): " << (case_a_v1_vs_fixed ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Case B V1 stable: " << (case_b_v1_stable ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Case B V1 not degraded: " << (case_b_v1_not_degraded ? "PASS" : "FAIL") << std::endl;

    bool pass = case_a_v1_vs_fixed && case_b_v1_stable && case_b_v1_not_degraded;
    std::cout << "  Overall: " << (pass ? "PASS" : "FAIL") << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_history_v1_case_A.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_history_v1_case_B.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_history_v1_summary.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_history_v1_diagnostics.csv" << std::endl;

    return pass ? 0 : 1;
}
