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
// Milestone 13: Adaptive Carry-Over v1
//
// This demo replaces fixed carry-over (alpha = constant) with adaptive alpha:
//   alpha = f(theta, match_quality, topology_event)
//
// Rules:
//   1. theta (normal change angle): larger theta -> smaller alpha
//   2. match quality (centroid distance): worse match -> smaller alpha  
//   3. topology events: newborn/unmatched -> alpha = 0
//
// Scene: Case A (rotational contact) and Case B (multi-patch)
//
// Output:
//   out/milestone_13/sdf_patch_adaptive_carryover_case_A.csv
//   out/milestone_13/sdf_patch_adaptive_carryover_case_B.csv
//   out/milestone_13/sdf_patch_adaptive_carryover_summary.csv
//   out/milestone_13/sdf_patch_adaptive_carryover_diagnostics.csv
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
#include <numeric>
#include <limits>

// -- Chrono includes --
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/core/ChFrame.h"

// -- OpenVDB includes --
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Interpolation.h>

using namespace chrono;

static double ClampSigned(double value) {
    return std::max(-1.0, std::min(1.0, value));
}

static ChVector3d SafeNormalize(const ChVector3d& v, const ChVector3d& fallback) {
    double len = v.Length();
    return len > 1.0e-14 ? v / len : fallback;
}

static ChVector3d ProjectToTangent(const ChVector3d& v, const ChVector3d& n) {
    return v - v.Dot(n) * n;
}

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
    int topology_birth_event; // 1 if born this frame
};

// =============================================================================
// Carry-over strategy
// =============================================================================

enum class CarryOverStrategy {
    Fixed,      // alpha = constant
    Adaptive,   // alpha = f(theta, match_quality, topology)
    AdaptiveV1,  // V1: geometry-consistent transfer with return mapping
    AdaptiveV21,  // V2.1: transported elastic state + validity gate
    AdaptiveV21b,  // V2.1b: V2.1 with normalized tail-aware validity gate
    AdaptiveV22,  // V2.2: V2.1b with direction-consistency gate
    AdaptiveV23,  // V2.3: V2.1b with direction-pollution gate
    AdaptiveV24  // V2.4: V2.1b with basis-distortion + direction-pollution gate
};

static std::string CarryOverStrategyName(CarryOverStrategy strategy) {
    if (strategy == CarryOverStrategy::AdaptiveV24) {
        return "adaptive-v2.4";
    }
    if (strategy == CarryOverStrategy::AdaptiveV23) {
        return "adaptive-v2.3";
    }
    if (strategy == CarryOverStrategy::AdaptiveV22) {
        return "adaptive-v2.2";
    }
    if (strategy == CarryOverStrategy::AdaptiveV21b) {
        return "adaptive-v2.1b";
    }
    if (strategy == CarryOverStrategy::AdaptiveV21) {
        return "adaptive-v2.1";
    }
    if (strategy == CarryOverStrategy::AdaptiveV1) {
        return "adaptive-v1";
    }
    if (strategy == CarryOverStrategy::Adaptive) {
        return "adaptive";
    }
    return "fixed";
}

// =============================================================================
// Transport limit mechanism flags
// =============================================================================

struct TransportConfig {
    std::string name;
    bool enable_clamp;
    bool enable_attenuation;
    CarryOverStrategy carry_strategy;
    double xi_max;
    double attenuation_beta;
    double fixed_carry_alpha;
    // Adaptive carry-over parameters
    double adaptive_alpha_base;    // Base alpha when theta=0
    double adaptive_theta_scale;   // Theta scale for decay (radians)
    double adaptive_match_scale;   // Match distance scale for decay (meters)
    double adaptive_slip_risk_scale; // V2.1: slip risk scale for weak correction (0=disabled)
};

// =============================================================================
// Adaptive alpha result
// =============================================================================

struct AdaptiveAlphaResult {
    double alpha;
    double theta;
    double match_distance;
    int topology_event; // 0=stable, 1=born, 2=drifted
};

// =============================================================================
// Transport result with diagnostics
// =============================================================================

struct TransportResult {
    double xi1_new;
    double xi2_new;
    double raw_transport_norm;
    double limited_transport_norm;
    double transport_attenuation;
    double carry_over_factor;
    int transport_clamp_hit;
    double frame_rotation_angle;
    double adaptive_alpha;
    double theta_ref = 0.0;
    double match_ref = 0.0;
    double theta_rel = 0.0;
    double match_rel = 0.0;
    double w_theta = 1.0;
    double w_match = 1.0;
    double w_topo = 1.0;
    double w_slip = 1.0;
    double risk_combo = 0.0;
    double w_tail = 1.0;
    double direction_consistency = 1.0;
    double w_dir = 1.0;
    double gate_before_direction = 1.0;
    double direction_pollution = 0.0;
    double w_pollution = 1.0;
    double gate_before_pollution = 1.0;
    double gate_after_pollution = 1.0;
    double tau_hist_z = 0.0;
    double tau_ref_z = 0.0;
    double basis_distortion = 0.0;
    double w_unified = 1.0;
    double gate_before_unified = 1.0;
    double gate_after_unified = 1.0;
    double xi_prev_basis_t1 = 0.0;
    double xi_prev_basis_t2 = 0.0;
    double xi_cur_basis_t1 = 0.0;
    double xi_cur_basis_t2 = 0.0;
    double xi_transport_dir_t1 = 0.0;
    double xi_transport_dir_t2 = 0.0;
    double delta_xi_dir_t1 = 0.0;
    double delta_xi_dir_t2 = 0.0;
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
    double avg_patch_count;
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
    double avg_frame_rotation_angle;
    double avg_theta;
    double max_theta;
    double avg_match_distance;
    double theta_p50;
    double theta_p95;
    double match_distance_p50;
    double match_distance_p95;
    double avg_theta_rel;
    double theta_rel_p50;
    double theta_rel_p95;
    double theta_rel_p99;
    double max_theta_rel;
    double avg_match_rel;
    double match_rel_p50;
    double match_rel_p95;
    double match_rel_p99;
    double max_match_rel;
    double avg_risk_combo;
    double risk_combo_p50;
    double risk_combo_p95;
    double risk_combo_p99;
    double max_risk_combo;
    double avg_gate_w_theta;
    double avg_gate_w_match;
    double avg_gate_w_topo;
    double avg_gate_w_slip;
    double avg_gate_w_tail;
    double avg_direction_consistency;
    double direction_consistency_p50;
    double direction_consistency_p95;
    double direction_consistency_p99;
    double avg_gate_w_dir;
    double avg_basis_distortion;
    double basis_distortion_p50;
    double basis_distortion_p95;
    double basis_distortion_p99;
    double avg_gate_w_unified;
    double basis_gate_before_avg;
    double basis_gate_after_avg;
    double avg_direction_pollution;
    double direction_pollution_p50;
    double direction_pollution_p95;
    double direction_pollution_p99;
    double avg_gate_w_pollution;
    double pollution_gate_before_avg;
    double pollution_gate_after_avg;
    double direction_bad_gate_before_avg;
    double direction_bad_gate_after_avg;
    double high_risk_gate_avg;
    double high_risk_v21_gate_avg;
    double median_risk_gate_avg;
    double median_risk_v21_gate_avg;
    int total_birth_events;

    // -- Clamp diagnostics --
    int clamp_trigger_count;
    double clamp_pre_norm_avg;
    double clamp_post_norm_avg;
    double clamp_max_pre_norm;
    double clamp_max_post_norm;

    // -- Patch traces --
    struct RunResultTrace {
        int frame_id;
        int persistent_id;
        double xi_max;
        bool clamp_enabled;
        double clamp_condition_value;
        bool clamp_triggered;
        double xi_history_norm;
        double xi_proj_norm;
        double xi_transfer_norm;
        double xi_before_increment_norm;
        double xi_after_increment_norm;
        double xi_before_clamp_norm;
        double xi_after_clamp_norm;
        double ft_trial_mag;
        double friction_limit;
        std::string stick_or_slip;
        double ft_final_mag;
    };
    std::vector<RunResultTrace> patch_traces;

    bool is_stable;
};

// =============================================================================
// Generate surface samples for sphere
// =============================================================================

std::vector<SurfaceSample> GenerateSphereSamples(double radius, int n_theta, int n_phi) {
    std::vector<SurfaceSample> samples;
    for (int i = 0; i < n_theta; i++) {
        double theta = M_PI * (i + 0.5) / n_theta;
        double theta_min = M_PI * i / n_theta;
        double theta_max = M_PI * (i + 1) / n_theta;
        double ring_area = 2.0 * M_PI * radius * radius * (std::cos(theta_min) - std::cos(theta_max));
        double sample_area = ring_area / n_phi;
        for (int j = 0; j < n_phi; j++) {
            double phi_angle = 2.0 * M_PI * j / n_phi;
            double x = radius * std::sin(theta) * std::cos(phi_angle);
            double y = radius * std::cos(theta);
            double z = radius * std::sin(theta) * std::sin(phi_angle);
            SurfaceSample s;
            s.local_pos = ChVector3d(x, y, z);
            s.area_weight = sample_area;
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
                s.area_weight = step * step;
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

static double Percentile(std::vector<double> values, double q) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    q = std::max(0.0, std::min(1.0, q));
    size_t idx = static_cast<size_t>(std::floor(q * static_cast<double>(values.size() - 1)));
    idx = std::min(idx, values.size() - 1);
    return values[idx];
}

static double MeanValue(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
}

static double MeanGateForRiskRange(
    const std::vector<double>& risks,
    const std::vector<double>& gates,
    double min_risk,
    double max_risk
) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < risks.size() && i < gates.size(); i++) {
        if (risks[i] >= min_risk && risks[i] <= max_risk) {
            sum += gates[i];
            count++;
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

static double MeanGateForDirectionMax(
    const std::vector<double>& direction_consistencies,
    const std::vector<double>& gates,
    double max_direction_consistency
) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < direction_consistencies.size() && i < gates.size(); i++) {
        if (direction_consistencies[i] <= max_direction_consistency) {
            sum += gates[i];
            count++;
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

static double MeanGateForPollutionMin(
    const std::vector<double>& pollutions,
    const std::vector<double>& gates,
    double min_pollution
) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < pollutions.size() && i < gates.size(); i++) {
        if (pollutions[i] >= min_pollution) {
            sum += gates[i];
            count++;
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

static double MeanGateForBasisMin(
    const std::vector<double>& distortions,
    const std::vector<double>& gates,
    double min_distortion
) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < distortions.size() && i < gates.size(); i++) {
        if (distortions[i] >= min_distortion) {
            sum += gates[i];
            count++;
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

// =============================================================================
// V2.1a gate normalization state
// =============================================================================

struct GateNormalizationStats {
    std::vector<double> theta_samples;
    std::vector<double> match_samples;

    double ThetaRef() const {
        return std::max(1.0e-4, MeanValue(theta_samples));
    }

    double MatchRef() const {
        return std::max(1.0e-4, MeanValue(match_samples));
    }

    void Observe(double theta, double match_distance) {
        if (std::isfinite(theta)) {
            theta_samples.push_back(std::max(0.0, theta));
        }
        if (std::isfinite(match_distance)) {
            match_samples.push_back(std::max(0.0, match_distance));
        }
    }
};

// =============================================================================
// Build orthonormal tangent basis from normal
// =============================================================================

void BuildTangentBasis(const ChVector3d& n, ChVector3d& t1, ChVector3d& t2) {
    ChVector3d nn = SafeNormalize(n, ChVector3d(0, 1, 0));
    ChVector3d ref = std::abs(nn.x()) < 0.9 ? ChVector3d(1, 0, 0) : ChVector3d(0, 1, 0);
    t1 = SafeNormalize(nn.Cross(ref), ChVector3d(0, 0, 1));
    t2 = SafeNormalize(nn.Cross(t1), ChVector3d(1, 0, 0));
}

// =============================================================================
// Compute adaptive carry-over alpha
// =============================================================================

AdaptiveAlphaResult ComputeAdaptiveAlpha(
    double theta,
    double match_distance,
    int topology_event,
    const TransportConfig& config
) {
    AdaptiveAlphaResult result;
    result.theta = theta;
    result.match_distance = match_distance;
    result.topology_event = topology_event;

    // Base alpha
    double alpha = config.adaptive_alpha_base;

    // Factor 1: Theta decay (rotation)
    // alpha *= exp(-theta / theta_scale)
    double theta_factor = std::exp(-theta / config.adaptive_theta_scale);
    alpha *= theta_factor;

    // Factor 2: Match quality decay
    // alpha *= exp(-dist / match_scale)
    double match_factor = std::exp(-match_distance / config.adaptive_match_scale);
    alpha *= match_factor;

    // Factor 3: Topology event override
    if (topology_event == 1) {
        // Newborn: zero carry-over
        alpha = 0.0;
    } else if (topology_event == 2) {
        // Drifted: significant reduction
        alpha *= 0.1;
    }

    // Clamp to [0, 1]
    alpha = std::max(0.0, std::min(1.0, alpha));

    result.alpha = alpha;
    return result;
}

// =============================================================================
// V1 Tangential History State (world-coordinate elastic state)
// =============================================================================

struct TangentialHistoryState {
    int persistent_id;
    ChVector3d xi_elastic_world;
    ChVector3d normal_prev;
    ChVector3d t1_prev;
    ChVector3d t2_prev;
    double match_quality_prev;
    double area_prev;
    int topology_state_prev;    // V2.1: 0=stable, 1=split/merge, 2=newborn/death
    int age;
    bool valid;
};

// =============================================================================
// V2.1 Transport operator: project elastic state to current tangent plane
// =============================================================================

ChVector3d TransportElasticStateToCurrentTangent(
    const ChVector3d& xi_elastic_world_prev,
    const ChVector3d& n_prev,
    const ChVector3d& n_cur
) {
    ChVector3d a = SafeNormalize(n_prev, ChVector3d(0, 1, 0));
    ChVector3d b = SafeNormalize(n_cur, a);
    ChVector3d xi_prev_tangent = ProjectToTangent(xi_elastic_world_prev, a);

    ChVector3d v = a.Cross(b);
    double c = ClampSigned(a.Dot(b));
    double s2 = v.Dot(v);

    ChVector3d xi_rotated = xi_prev_tangent;
    if (s2 > 1.0e-24) {
        ChVector3d vx = v.Cross(xi_prev_tangent);
        ChVector3d vvx = v.Cross(vx);
        xi_rotated = xi_prev_tangent + vx + vvx * ((1.0 - c) / s2);
    } else if (c < 0.0) {
        // Anti-parallel normals make the shortest rotation axis underdetermined.
        // Keep the non-amplifying part only; topology/gate should usually reset it.
        xi_rotated = ProjectToTangent(xi_prev_tangent, b);
    }

    return ProjectToTangent(xi_rotated, b);
}

// =============================================================================
// V2.1 Validity gate: compute history retention weight w in [0, 1]
// =============================================================================

double ValidityGateV21(
    double theta,
    double match_distance,
    int topology_event,
    double slip_risk,
    const TransportConfig& config
) {
    double w_theta = std::exp(-theta / config.adaptive_theta_scale);
    double w_match = std::exp(-match_distance / config.adaptive_match_scale);
    
    double w_topo;
    if (topology_event == 1) {
        w_topo = 0.0;
    } else if (topology_event == 2) {
        w_topo = 0.3;
    } else {
        w_topo = 1.0;
    }
    
    double w_slip = 1.0;
    if (config.adaptive_slip_risk_scale > 1e-12) {
        double lambda = 0.15;
        double excess = std::max(0.0, slip_risk - 1.0);
        w_slip = 1.0 / (1.0 + lambda * excess);
    }
    
    double w = w_theta * w_match * w_topo * w_slip;
    return std::max(0.0, std::min(1.0, w));
}

// =============================================================================
// V2.1b Validity gate: normalized history validity gate with tail suppression
// =============================================================================

double ValidityGateV21b(
    double theta,
    double match_distance,
    int topology_event,
    double slip_risk,
    const TransportConfig& config,
    const GateNormalizationStats& gate_stats,
    TransportResult& result
) {
    result.theta_ref = gate_stats.ThetaRef();
    result.match_ref = gate_stats.MatchRef();
    result.theta_rel = theta / (result.theta_ref + 1.0e-12);
    result.match_rel = match_distance / (result.match_ref + 1.0e-12);
    result.risk_combo = std::max(result.theta_rel, result.match_rel);

    const double theta_rel_scale = std::max(1.0e-12, config.adaptive_theta_scale);
    const double match_rel_scale = std::max(1.0e-12, config.adaptive_match_scale);
    result.w_theta = std::exp(-result.theta_rel / theta_rel_scale);
    result.w_match = std::exp(-result.match_rel / match_rel_scale);

    if (topology_event == 1) {
        result.w_topo = 0.0;
    } else if (topology_event == 2) {
        result.w_topo = 0.3;
    } else {
        result.w_topo = 1.0;
    }

    result.w_slip = 1.0;
    if (config.adaptive_slip_risk_scale > 1.0e-12) {
        double excess = std::max(0.0, slip_risk - 1.0);
        result.w_slip = std::exp(-config.adaptive_slip_risk_scale * excess);
        result.w_slip = std::max(0.85, std::min(1.0, result.w_slip));
    }

    const double risk_knee = 1.5;
    const double tail_scale = 0.5;
    result.w_tail = 1.0;
    if (result.risk_combo > risk_knee) {
        result.w_tail = std::exp(-(result.risk_combo - risk_knee) / tail_scale);
    }

    double w = result.w_theta * result.w_match * result.w_topo * result.w_slip * result.w_tail;
    return std::max(0.0, std::min(1.0, w));
}

double DirectionConsistencyGate(
    const ChVector3d& xi_transport,
    const ChVector3d& delta_xi,
    TransportResult& result
) {
    double xi_norm = xi_transport.Length();
    double delta_norm = delta_xi.Length();
    result.direction_consistency = 1.0;
    result.w_dir = 1.0;

    if (xi_norm <= 1.0e-12 || delta_norm <= 1.0e-12) {
        return 1.0;
    }

    double cos_dir = xi_transport.Dot(delta_xi) / (xi_norm * delta_norm);
    cos_dir = std::max(-1.0, std::min(1.0, cos_dir));
    double directional_alignment = std::max(0.0, cos_dir);
    double direction_evidence = std::min(1.0, delta_norm / (0.5 * xi_norm + 1.0e-12));
    result.direction_consistency = directional_alignment * direction_evidence;

    const double dir_knee = 0.75;
    const double dir_floor = 0.90;
    const double theta_rel_activation = 0.65;
    if (result.theta_rel < theta_rel_activation) {
        return 1.0;
    }
    if (result.direction_consistency < dir_knee) {
        result.w_dir = dir_floor + (1.0 - dir_floor) * (result.direction_consistency / dir_knee);
    }

    return result.w_dir;
}

double DirectionPollutionGate(
    const ChVector3d& xi_transport,
    const ChVector3d& delta_xi,
    const ChVector3d& tangential_velocity,
    const ChVector3d& normal_force,
    const ChVector3d& patch_moment_arm,
    double tangential_stiffness,
    double tangential_damping,
    TransportResult& result
) {
    result.direction_pollution = 0.0;
    result.w_pollution = 1.0;
    result.gate_before_pollution = result.gate_before_direction;
    result.gate_after_pollution = result.gate_before_direction;
    result.tau_hist_z = 0.0;
    result.tau_ref_z = 0.0;

    ChVector3d ft_hist = -tangential_stiffness * xi_transport;
    ChVector3d ft_current = -tangential_stiffness * delta_xi - tangential_damping * tangential_velocity;
    ChVector3d ref_force = normal_force + ft_current;

    result.tau_hist_z = patch_moment_arm.Cross(ft_hist).z();
    result.tau_ref_z = patch_moment_arm.Cross(ref_force).z();

    double hist_abs = std::abs(result.tau_hist_z);
    double ref_abs = std::abs(result.tau_ref_z);
    if (hist_abs <= 1.0e-12 || ref_abs <= 1.0e-12) {
        return 1.0;
    }

    double opposition = std::max(
        0.0,
        -(result.tau_hist_z * result.tau_ref_z) / (hist_abs * ref_abs + 1.0e-18)
    );
    double hist_share = hist_abs / (hist_abs + ref_abs + 1.0e-18);
    result.direction_pollution = std::max(0.0, std::min(1.0, opposition * hist_share));

    const double theta_rel_activation = 0.65;
    const double gate_activation = 0.05;
    const double pollution_knee = 0.10;
    const double pollution_strength = 0.30;
    const double pollution_floor = 0.70;

    if (result.theta_rel < theta_rel_activation ||
        result.gate_before_direction < gate_activation ||
        result.direction_pollution <= pollution_knee) {
        return 1.0;
    }

    double excess = (result.direction_pollution - pollution_knee) / (1.0 - pollution_knee);
    excess = std::max(0.0, std::min(1.0, excess));
    result.w_pollution = std::max(pollution_floor, 1.0 - pollution_strength * excess);
    result.gate_after_pollution = result.gate_before_direction * result.w_pollution;
    return result.w_pollution;
}

double BasisDistortionPollutionGate(
    double xi1_old,
    double xi2_old,
    const ChVector3d& t1_old,
    const ChVector3d& t2_old,
    const ChVector3d& xi_transport,
    const ChVector3d& t1_new,
    const ChVector3d& t2_new,
    const ChVector3d& delta_xi,
    const ChVector3d& tangential_velocity,
    const ChVector3d& normal_force,
    const ChVector3d& patch_moment_arm,
    double tangential_stiffness,
    double tangential_damping,
    TransportResult& result
) {
    result.basis_distortion = 0.0;
    result.direction_pollution = 0.0;
    result.w_unified = 1.0;
    result.w_pollution = 1.0;
    result.gate_before_unified = result.gate_before_direction;
    result.gate_after_unified = result.gate_before_direction;
    result.gate_before_pollution = result.gate_before_direction;
    result.gate_after_pollution = result.gate_before_direction;
    result.tau_hist_z = 0.0;
    result.tau_ref_z = 0.0;
    result.xi_prev_basis_t1 = xi1_old;
    result.xi_prev_basis_t2 = xi2_old;
    result.xi_cur_basis_t1 = xi_transport.Dot(t1_new);
    result.xi_cur_basis_t2 = xi_transport.Dot(t2_new);
    ChVector3d xi_prev_world_from_basis = xi1_old * t1_old + xi2_old * t2_old;

    double prev_norm = std::sqrt(xi1_old * xi1_old + xi2_old * xi2_old);
    double cur_norm = std::sqrt(result.xi_cur_basis_t1 * result.xi_cur_basis_t1 +
                                result.xi_cur_basis_t2 * result.xi_cur_basis_t2);
    if (prev_norm > 1.0e-12 && cur_norm > 1.0e-12) {
        double projection_change = (xi_prev_world_from_basis - xi_transport).Length() / prev_norm;
        double projection_distortion = std::max(0.0, std::min(1.0, projection_change / 0.01));
        double u_prev1 = xi1_old / prev_norm;
        double u_prev2 = xi2_old / prev_norm;
        double cur1 = result.xi_cur_basis_t1 / cur_norm;
        double cur2 = result.xi_cur_basis_t2 / cur_norm;
        double pred1 = t1_old.Dot(t1_new) * u_prev1 + t2_old.Dot(t1_new) * u_prev2;
        double pred2 = t1_old.Dot(t2_new) * u_prev1 + t2_old.Dot(t2_new) * u_prev2;
        double pred_norm = std::sqrt(pred1 * pred1 + pred2 * pred2);
        double component_dot = std::max(-1.0, std::min(1.0, u_prev1 * cur1 + u_prev2 * cur2));
        double component_distortion = 0.5 * (1.0 - component_dot);
        double axis_distortion = std::abs(u_prev1) * 0.5 * (1.0 - std::max(-1.0, std::min(1.0, t1_old.Dot(t1_new)))) +
                                 std::abs(u_prev2) * 0.5 * (1.0 - std::max(-1.0, std::min(1.0, t2_old.Dot(t2_new))));
        axis_distortion = std::max(0.0, std::min(1.0, axis_distortion));
        double prediction_distortion = 0.0;
        if (pred_norm > 1.0e-12) {
            pred1 /= pred_norm;
            pred2 /= pred_norm;
            double dir_dot = std::max(-1.0, std::min(1.0, pred1 * cur1 + pred2 * cur2));
            prediction_distortion = 0.5 * (1.0 - dir_dot);
        }
        bool prev_t1_dominant = std::abs(xi1_old) >= std::abs(xi2_old);
        bool cur_t1_dominant = std::abs(result.xi_cur_basis_t1) >= std::abs(result.xi_cur_basis_t2);
        double dominant_swap = (prev_t1_dominant != cur_t1_dominant) ? 1.0 : 0.0;
        double prev_dom = prev_t1_dominant ? u_prev1 : u_prev2;
        double cur_dom = prev_t1_dominant ? cur1 : cur2;
        double dominant_sign_flip = (prev_dom * cur_dom < 0.0) ? 1.0 : 0.0;
        result.basis_distortion = std::max(
            0.0,
            std::min(1.0,
                     0.50 * component_distortion +
                     0.25 * axis_distortion +
                     0.15 * dominant_swap +
                     0.05 * dominant_sign_flip +
                     0.05 * prediction_distortion)
        );
        result.basis_distortion = std::max(result.basis_distortion, projection_distortion);
    }

    ChVector3d ft_hist = -tangential_stiffness * xi_transport;
    ChVector3d ft_current = -tangential_stiffness * delta_xi - tangential_damping * tangential_velocity;
    ChVector3d ref_force = normal_force + ft_current;
    result.tau_hist_z = patch_moment_arm.Cross(ft_hist).z();
    result.tau_ref_z = patch_moment_arm.Cross(ref_force).z();

    double hist_abs = std::abs(result.tau_hist_z);
    double ref_abs = std::abs(result.tau_ref_z);
    if (hist_abs > 1.0e-12 && ref_abs > 1.0e-12) {
        double opposition = std::max(
            0.0,
            -(result.tau_hist_z * result.tau_ref_z) / (hist_abs * ref_abs + 1.0e-18)
        );
        double hist_share = hist_abs / (hist_abs + ref_abs + 1.0e-18);
        result.direction_pollution = std::max(0.0, std::min(1.0, opposition * hist_share));
    }

    const double theta_rel_activation = 0.65;
    const double gate_activation = 0.05;
    const double pollution_knee = 0.10;
    const double pollution_strength = 0.30;
    const double pollution_floor = 0.70;
    if (result.theta_rel < theta_rel_activation ||
        result.gate_before_direction < gate_activation ||
        result.direction_pollution <= pollution_knee) {
        return 1.0;
    }

    double pollution_excess = (result.direction_pollution - pollution_knee) / (1.0 - pollution_knee);
    pollution_excess = std::max(0.0, std::min(1.0, pollution_excess));
    result.w_pollution = std::max(pollution_floor, 1.0 - pollution_strength * pollution_excess);
    result.gate_after_pollution = result.gate_before_direction * result.w_pollution;

    const double basis_knee = 0.08;
    const double basis_strength = 0.50;
    const double basis_floor = 0.70;
    double basis_excess = 0.0;
    if (result.basis_distortion > basis_knee) {
        basis_excess = (result.basis_distortion - basis_knee) / (1.0 - basis_knee);
        basis_excess = std::max(0.0, std::min(1.0, basis_excess));
    }
    double w_basis_extra = std::max(basis_floor, 1.0 - basis_strength * basis_excess * pollution_excess);
    result.w_unified = result.w_pollution * w_basis_extra;
    result.gate_after_unified = result.gate_before_direction * result.w_unified;
    result.gate_after_pollution = result.gate_after_unified;
    return result.w_unified;
}

// =============================================================================
// V2.1: Transported elastic state + validity gate + unchanged return mapping
// =============================================================================

TransportResult TransportWithCarryOverV21(
    double xi1_old, double xi2_old,
    const ChVector3d& t1_old, const ChVector3d& t2_old, const ChVector3d& n_old,
    const ChVector3d& t1_new, const ChVector3d& t2_new, const ChVector3d& n_new,
    double match_distance,
    int topology_event,
    double vt1, double vt2,
    double time_step,
    double tangential_stiffness,
    double tangential_damping,
    double friction_coefficient,
    double fn_magnitude,
    bool is_newborn,
    const TangentialHistoryState* history_state,
    const TransportConfig& config
) {
    TransportResult result;
    result.raw_transport_norm = 0.0;
    result.limited_transport_norm = 0.0;
    result.transport_attenuation = 1.0;
    result.carry_over_factor = 1.0;
    result.transport_clamp_hit = 0;
    result.frame_rotation_angle = 0.0;
    result.adaptive_alpha = 0.0;
    result.xi1_new = 0.0;
    result.xi2_new = 0.0;

    if (history_state == nullptr || !history_state->valid || is_newborn) {
        result.xi1_new = vt1 * time_step;
        result.xi2_new = vt2 * time_step;
        result.limited_transport_norm = std::sqrt(result.xi1_new * result.xi1_new + result.xi2_new * result.xi2_new);
        return result;
    }

    // Step 1: Read old xi_elastic from persistent patch
    ChVector3d xi_elastic_world_prev = history_state->xi_elastic_world;

    // Step 2: Transport xi_elastic to current tangent plane
    ChVector3d xi_transport = TransportElasticStateToCurrentTangent(
        xi_elastic_world_prev, history_state->normal_prev, n_new
    );
    result.raw_transport_norm = xi_transport.Length();
    result.frame_rotation_angle = std::acos(std::max(-1.0, std::min(1.0, history_state->normal_prev.Dot(n_new))));

    // Step 3: Validity gate
    double theta = result.frame_rotation_angle;
    double delta_xi1 = vt1 * time_step;
    double delta_xi2 = vt2 * time_step;
    ChVector3d xi_trial_estimate = xi_transport + delta_xi1 * t1_new + delta_xi2 * t2_new;
    ChVector3d tau_trial_estimate = tangential_stiffness * xi_trial_estimate + tangential_damping * (vt1 * t1_new + vt2 * t2_new);
    double tau_trial_mag_estimate = tau_trial_estimate.Length();
    double friction_limit_estimate = friction_coefficient * fn_magnitude;
    double slip_risk = (friction_limit_estimate > 1e-12) ? (tau_trial_mag_estimate / friction_limit_estimate) : 0.0;

    double w = ValidityGateV21(theta, match_distance, topology_event, slip_risk, config);
    result.adaptive_alpha = w;
    result.carry_over_factor = w;

    ChVector3d xi_valid = xi_transport * w;

    // Step 4: Add current step increment
    ChVector3d delta_xi = delta_xi1 * t1_new + delta_xi2 * t2_new;
    ChVector3d xi_trial = xi_valid + delta_xi;

    // The common force path below performs the single Coulomb return mapping.
    result.xi1_new = xi_trial.Dot(t1_new);
    result.xi2_new = xi_trial.Dot(t2_new);
    result.limited_transport_norm = xi_valid.Length();
    return result;
}

// =============================================================================
// V2.1b: transported elastic state + normalized tail-aware validity gate + unchanged return mapping
// =============================================================================

TransportResult TransportWithCarryOverV21b(
    double xi1_old, double xi2_old,
    const ChVector3d& t1_old, const ChVector3d& t2_old, const ChVector3d& n_old,
    const ChVector3d& t1_new, const ChVector3d& t2_new, const ChVector3d& n_new,
    double match_distance,
    int topology_event,
    double vt1, double vt2,
    double time_step,
    double tangential_stiffness,
    double tangential_damping,
    double friction_coefficient,
    double fn_magnitude,
    bool is_newborn,
    const TangentialHistoryState* history_state,
    const TransportConfig& config,
    GateNormalizationStats& gate_stats,
    bool enable_direction_gate,
    bool enable_pollution_gate,
    bool enable_unified_gate,
    const ChVector3d& patch_moment_arm,
    const ChVector3d& normal_force
) {
    TransportResult result;
    result.raw_transport_norm = 0.0;
    result.limited_transport_norm = 0.0;
    result.transport_attenuation = 1.0;
    result.carry_over_factor = 1.0;
    result.transport_clamp_hit = 0;
    result.frame_rotation_angle = 0.0;
    result.adaptive_alpha = 0.0;
    result.xi1_new = 0.0;
    result.xi2_new = 0.0;

    if (history_state == nullptr || !history_state->valid || is_newborn) {
        result.xi1_new = vt1 * time_step;
        result.xi2_new = vt2 * time_step;
        result.limited_transport_norm = std::sqrt(result.xi1_new * result.xi1_new + result.xi2_new * result.xi2_new);
        return result;
    }

    // Step 1: Read old xi_elastic from persistent patch
    ChVector3d xi_elastic_world_prev = history_state->xi_elastic_world;

    // Step 2: Transport xi_elastic to current tangent plane
    ChVector3d xi_transport = TransportElasticStateToCurrentTangent(
        xi_elastic_world_prev, history_state->normal_prev, n_new
    );
    result.raw_transport_norm = xi_transport.Length();
    result.xi_prev_basis_t1 = xi1_old;
    result.xi_prev_basis_t2 = xi2_old;
    result.xi_cur_basis_t1 = xi_transport.Dot(t1_new);
    result.xi_cur_basis_t2 = xi_transport.Dot(t2_new);
    if (result.raw_transport_norm > 1.0e-12) {
        result.xi_transport_dir_t1 = xi_transport.Dot(t1_new) / result.raw_transport_norm;
        result.xi_transport_dir_t2 = xi_transport.Dot(t2_new) / result.raw_transport_norm;
    }
    result.frame_rotation_angle = std::acos(std::max(-1.0, std::min(1.0, history_state->normal_prev.Dot(n_new))));

    // Step 3: Normalized tail-aware validity gate
    double theta = result.frame_rotation_angle;
    double delta_xi1 = vt1 * time_step;
    double delta_xi2 = vt2 * time_step;
    double delta_xi_norm = std::sqrt(delta_xi1 * delta_xi1 + delta_xi2 * delta_xi2);
    if (delta_xi_norm > 1.0e-12) {
        result.delta_xi_dir_t1 = delta_xi1 / delta_xi_norm;
        result.delta_xi_dir_t2 = delta_xi2 / delta_xi_norm;
    }
    ChVector3d xi_trial_estimate = xi_transport + delta_xi1 * t1_new + delta_xi2 * t2_new;
    ChVector3d tau_trial_estimate = tangential_stiffness * xi_trial_estimate + tangential_damping * (vt1 * t1_new + vt2 * t2_new);
    double tau_trial_mag_estimate = tau_trial_estimate.Length();
    double friction_limit_estimate = friction_coefficient * fn_magnitude;
    double slip_risk = (friction_limit_estimate > 1e-12) ? (tau_trial_mag_estimate / friction_limit_estimate) : 0.0;

    double w = ValidityGateV21b(theta, match_distance, topology_event, slip_risk, config, gate_stats, result);
    gate_stats.Observe(theta, match_distance);
    result.gate_before_direction = w;
    ChVector3d delta_xi = delta_xi1 * t1_new + delta_xi2 * t2_new;
    ChVector3d tangential_velocity = vt1 * t1_new + vt2 * t2_new;
    if (enable_unified_gate) {
        w *= BasisDistortionPollutionGate(
            xi1_old,
            xi2_old,
            t1_old,
            t2_old,
            xi_transport,
            t1_new,
            t2_new,
            delta_xi,
            tangential_velocity,
            normal_force,
            patch_moment_arm,
            tangential_stiffness,
            tangential_damping,
            result
        );
    } else if (enable_pollution_gate) {
        w *= DirectionPollutionGate(
            xi_transport,
            delta_xi,
            tangential_velocity,
            normal_force,
            patch_moment_arm,
            tangential_stiffness,
            tangential_damping,
            result
        );
    } else if (enable_direction_gate) {
        w *= DirectionConsistencyGate(xi_transport, delta_xi1 * t1_new + delta_xi2 * t2_new, result);
    }
    result.adaptive_alpha = w;
    result.carry_over_factor = w;

    ChVector3d xi_valid = xi_transport * w;

    // Step 4: Add current step increment
    ChVector3d xi_trial = xi_valid + delta_xi;

    // The common force path below performs the single Coulomb return mapping.
    result.xi1_new = xi_trial.Dot(t1_new);
    result.xi2_new = xi_trial.Dot(t2_new);
    result.limited_transport_norm = xi_valid.Length();
    return result;
}

// =============================================================================
// V1: Geometry-consistent tangential history transfer with return mapping
// =============================================================================

TransportResult TransportWithCarryOverV1(
    double xi1_old, double xi2_old,
    const ChVector3d& t1_old, const ChVector3d& t2_old, const ChVector3d& n_old,
    const ChVector3d& t1_new, const ChVector3d& t2_new, const ChVector3d& n_new,
    double match_distance,
    int topology_event,
    double vt1, double vt2,
    double time_step,
    double tangential_stiffness,
    double tangential_damping,
    double friction_coefficient,
    double fn_magnitude,
    bool is_newborn,
    const TangentialHistoryState* history_state,
    const TransportConfig& config
) {
    TransportResult result;
    result.raw_transport_norm = 0.0;
    result.limited_transport_norm = 0.0;
    result.transport_attenuation = 1.0;
    result.carry_over_factor = 1.0;
    result.transport_clamp_hit = 0;
    result.frame_rotation_angle = 0.0;
    result.adaptive_alpha = 0.0;
    result.xi1_new = 0.0;
    result.xi2_new = 0.0;

    if (history_state == nullptr || !history_state->valid || is_newborn) {
        result.xi1_new = vt1 * time_step;
        result.xi2_new = vt2 * time_step;
        result.limited_transport_norm = std::sqrt(result.xi1_new * result.xi1_new + result.xi2_new * result.xi2_new);
        return result;
    }

    result.frame_rotation_angle = std::acos(std::max(-1.0, std::min(1.0, history_state->normal_prev.Dot(n_new))));

    ChVector3d xi_proj = history_state->xi_elastic_world;
    double normal_component = xi_proj.Dot(n_new);
    xi_proj = xi_proj - normal_component * n_new;

    result.raw_transport_norm = xi_proj.Length();

    double theta = result.frame_rotation_angle;
    double alpha = config.adaptive_alpha_base;
    alpha *= std::exp(-theta / config.adaptive_theta_scale);
    alpha *= std::exp(-match_distance / config.adaptive_match_scale);
    alpha = std::max(0.0, std::min(1.0, alpha));

    result.adaptive_alpha = alpha;

    ChVector3d xi_transfer = xi_proj * alpha;

    double delta_xi1 = vt1 * time_step;
    double delta_xi2 = vt2 * time_step;
    ChVector3d delta_xi = delta_xi1 * t1_new + delta_xi2 * t2_new;
    ChVector3d xi_trial = xi_transfer + delta_xi;

    result.xi1_new = xi_trial.Dot(t1_new);
    result.xi2_new = xi_trial.Dot(t2_new);
    result.limited_transport_norm = xi_transfer.Length();
    result.carry_over_factor = alpha;

    return result;
}

// =============================================================================
// Transport with carry-over
// =============================================================================

TransportResult TransportWithCarryOver(
    double xi1_old, double xi2_old,
    const ChVector3d& t1_old, const ChVector3d& t2_old, const ChVector3d& n_old,
    const ChVector3d& t1_new, const ChVector3d& t2_new, const ChVector3d& n_new,
    double match_distance,
    int topology_event,
    const TransportConfig& config
) {
    TransportResult result;
    result.raw_transport_norm = 0.0;
    result.limited_transport_norm = 0.0;
    result.transport_attenuation = 1.0;
    result.carry_over_factor = 1.0;
    result.transport_clamp_hit = 0;
    result.frame_rotation_angle = 0.0;
    result.adaptive_alpha = 1.0;

    // Step 1: Convert old local state to world coordinates
    ChVector3d xi_world = xi1_old * t1_old + xi2_old * t2_old;

    // Step 2: Project onto new tangent plane
    double normal_component = xi_world.Dot(n_new);
    ChVector3d xi_projected = xi_world - normal_component * n_new;

    result.raw_transport_norm = xi_projected.Length();
    result.frame_rotation_angle = std::acos(std::max(-1.0, std::min(1.0, n_old.Dot(n_new))));

    // Step 3: Apply rotation-aware attenuation if enabled
    if (config.enable_attenuation) {
        double theta = result.frame_rotation_angle;
        double attenuation = std::exp(-config.attenuation_beta * theta);
        result.transport_attenuation = attenuation;
        xi_projected *= attenuation;
    }

    // Step 4: Apply magnitude clamp if enabled
    if (config.enable_clamp) {
        double xi_proj_norm = xi_projected.Length();
        if (xi_proj_norm > config.xi_max) {
            xi_projected = xi_projected * (config.xi_max / xi_proj_norm);
            result.transport_clamp_hit = 1;
        }
    }

    // Step 5: Apply carry-over (fixed or adaptive)
    // NOTE: When fixed_carry_alpha <= 0, skip carry-over entirely to match
    // MS12 semantics where enable_carry_over=false means no carry-over multiplication
    double carry_factor = 1.0;
    if (config.carry_strategy == CarryOverStrategy::Fixed) {
        if (config.fixed_carry_alpha > 0.0) {
            carry_factor = config.fixed_carry_alpha;
            xi_projected *= carry_factor;
        }
        // else: carry_factor = 1.0, skip multiplication (equivalent to enable_carry_over=false)
    } else if (config.carry_strategy == CarryOverStrategy::Adaptive) {
        AdaptiveAlphaResult alpha_result = ComputeAdaptiveAlpha(
            result.frame_rotation_angle,
            match_distance,
            topology_event,
            config
        );
        carry_factor = alpha_result.alpha;
        result.adaptive_alpha = carry_factor;
        if (carry_factor > 0.0) {
            xi_projected *= carry_factor;
        }
    }
    result.carry_over_factor = carry_factor;

    result.limited_transport_norm = xi_projected.Length();

    // Step 6: Represent in new local basis
    result.xi1_new = xi_projected.Dot(t1_new);
    result.xi2_new = xi_projected.Dot(t2_new);

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
    const std::string& mode_name,
    const std::string& out_dir
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

    // -- V1 tangential history state map --
    std::map<int, TangentialHistoryState> tangential_history;
    GateNormalizationStats gate_stats;

    // -- Buffers --
    std::vector<SDFQueryResult> sdf_results(samples.size());
    std::vector<ChVector3d> sample_world_positions(samples.size());

    // -- Metrics --
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
    std::vector<double> theta_rel_values;
    std::vector<double> match_rel_values;
    std::vector<double> risk_combo_values;
    std::vector<double> v21_abs_gate_values;
    std::vector<double> final_gate_values;
    std::vector<double> gate_w_theta_values;
    std::vector<double> gate_w_match_values;
    std::vector<double> gate_w_topo_values;
    std::vector<double> gate_w_slip_values;
    std::vector<double> gate_w_tail_values;
    std::vector<double> direction_consistency_values;
    std::vector<double> gate_w_dir_values;
    std::vector<double> gate_before_direction_values;
    std::vector<double> direction_pollution_values;
    std::vector<double> gate_w_pollution_values;
    std::vector<double> gate_before_pollution_values;
    std::vector<double> basis_distortion_values;
    std::vector<double> gate_w_unified_values;
    std::vector<double> gate_before_unified_values;
    int total_birth_events = 0;
    int total_transport_clamp_hits = 0;
    double total_force_y = 0.0;
    int force_sample_count = 0;
    double sum_tangential_norm = 0.0;
    double max_tangential_norm = 0.0;
    double sum_tangential_ratio = 0.0;
    int total_stick_steps = 0;
    int total_slip_steps = 0;
    int total_stick_to_slip_transitions = 0;

    // -- Clamp diagnostics --
    int clamp_trigger_count = 0;
    double clamp_pre_norm_sum = 0.0;
    double clamp_post_norm_sum = 0.0;
    double clamp_max_pre_norm = 0.0;
    double clamp_max_post_norm = 0.0;
    int clamp_sample_count = 0;

    // -- Global xi_before_clamp distribution --
    std::vector<double> all_xi_before_clamp_values;
    
    // -- Max xi_before_clamp tracker for automatic trace selection --
    struct MaxXiTracker {
        double max_xi = 0.0;
        int max_frame = -1;
        int max_patch_idx = -1;
        int max_persistent_id = -1;
        double xi_history_norm = 0.0;
        double xi_proj_norm = 0.0;
        double xi_transfer_norm = 0.0;
        double xi_after_increment_norm = 0.0;
        double xi_before_clamp_norm = 0.0;
        double xi_after_clamp_norm = 0.0;
        double xi_before_return_mapping_norm = 0.0;
        double xi_after_return_mapping_norm = 0.0;
        double ft_trial_mag = 0.0;
        double friction_limit = 0.0;
        std::string stick_or_slip = "not_applicable";
        double ft_final_mag = 0.0;
        double fn_mag = 0.0;
        bool is_valid = false;  // True only if full force computation completed
    };
    MaxXiTracker max_xi_tracker;

    // -- Patch trace for diagnostic --
    struct PatchTrace {
        int frame_id;
        int persistent_id;
        std::string config_name;
        double xi_max;
        bool clamp_enabled;
        double clamp_condition_value;
        bool clamp_triggered;
        double xi_history_norm;
        double xi_proj_norm;
        double xi_transfer_norm;
        double xi_before_increment_norm;
        double xi_after_increment_norm;
        double xi_before_clamp_norm;
        double xi_after_clamp_norm;
        double xi_before_return_mapping_norm;
        double xi_after_return_mapping_norm;
        double ft_trial_mag;
        double friction_limit;
        std::string stick_or_slip;
        double ft_final_mag;
    };
    std::vector<PatchTrace> patch_traces;

    struct GateSampleTrace {
        bool valid = false;
        int frame_id = -1;
        int patch_index = -1;
        int persistent_id = -1;
        double theta = 0.0;
        double theta_ref = 0.0;
        double theta_rel = 0.0;
        double match_distance = 0.0;
        double match_ref = 0.0;
        double match_rel = 0.0;
        double risk_combo = 0.0;
        double xi_transport_norm = 0.0;
        double delta_xi_norm = 0.0;
        double xi_prev_basis_t1 = 0.0;
        double xi_prev_basis_t2 = 0.0;
        double xi_cur_basis_t1 = 0.0;
        double xi_cur_basis_t2 = 0.0;
        double xi_transport_dir_t1 = 0.0;
        double xi_transport_dir_t2 = 0.0;
        double delta_xi_dir_t1 = 0.0;
        double delta_xi_dir_t2 = 0.0;
        double direction_consistency = 1.0;
        double gate_before_direction = 1.0;
        double gate_after_direction = 1.0;
        double direction_pollution = 0.0;
        double gate_before_pollution = 1.0;
        double gate_after_pollution = 1.0;
        double basis_distortion = 0.0;
        double gate_before_unified = 1.0;
        double gate_after_unified = 1.0;
        double w_theta = 1.0;
        double w_match = 1.0;
        double w_topo = 1.0;
        double w_slip = 1.0;
        double w_tail = 1.0;
        double w_dir = 1.0;
        double w_pollution = 1.0;
        double w_unified = 1.0;
        double tau_hist_z = 0.0;
        double tau_ref_z = 0.0;
        double validity_gate = 1.0;
        double v21_abs_w_theta = 1.0;
        double v21_abs_w_match = 1.0;
        double v21_abs_gate = 1.0;
        double raw_transport_norm = 0.0;
        double limited_transport_norm = 0.0;
        double ft_trial_mag = 0.0;
        double friction_limit = 0.0;
        std::string stick_or_slip = "unknown";
        double ft_final_mag = 0.0;
        double score = -1.0;
    };
    GateSampleTrace gate_sample_trace;
    GateSampleTrace gate_normal_sample_trace;

    // Trace target: select specific frame/patch for detailed output
    // Use frame numbers that exist in simulation (~4000 frames for 2s at dt=5e-4)
    const int trace_target_frame = 100;  // Trace at frame 100 (early steady state)

    // -- Matching thresholds --
    double match_distance_threshold = 0.5 * (config.dyn_body_type == 0 ? config.dyn_body_radius : config.dyn_body_size);
    double normal_similarity_threshold = 0.9;

    // -- Frame patch ID offset for uniqueness --
    int frame_offset = 0;

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

            double phi = sdf_results[si].phi;
            if (phi < config.activation_band) {
                active_indices.push_back(static_cast<int>(si));
            }
        }

        // Patch-based force computation
        int frame_patch_id_offset = frame_offset * 1000;
        frame_offset++;
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
            track.topology_birth_event = 0;
            if (m.matched) {
                for (const auto& prev : active_tracks) {
                    if (prev.persistent_id == m.persistent_id) {
                        track = prev;
                        track.status = TrackStatus::Alive;
                        track.topology_birth_event = 0; // Reset birth flag for existing track
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

        // Compute patch forces
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

            // --- Transport with carry-over ---
            double xi1 = 0.0, xi2 = 0.0;
            double raw_trans_norm = 0.0, limited_trans_norm = 0.0;
            double attenuation = 1.0, carry_factor = 1.0;
            double adaptive_alpha = 1.0;
            int clamp_hit = 0;
            double frame_rot_angle = 0.0;
            double match_dist = 0.0;
            int topology_event = 0;
            bool is_newborn = (prev_track == nullptr || prev_track->topology_birth_event == 1);
            bool is_v22 = (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22);
            bool is_v23 = (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23);
            bool is_v24 = (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24);

            // Find V1 history state
            const TangentialHistoryState* v1_history = nullptr;
            if ((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV1 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) &&
                patch.matched_persistent_id >= 0) {
                auto it = tangential_history.find(patch.matched_persistent_id);
                if (it != tangential_history.end() && it->second.valid) {
                    v1_history = &it->second;
                }
            }

            if (prev_track != nullptr) {
                // Determine topology event
                if (prev_track->topology_birth_event) {
                    topology_event = 1; // born
                }
                double center_drift = (patch.center - prev_track->center).Length();
                if (center_drift > patch.patch_radius_estimate * 2.0) {
                    topology_event = 2; // drifted
                }

                match_dist = center_drift;

                if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV1 ||
                    transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21 ||
                    transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                    transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                    transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                    transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) {
                    TransportResult trans;
                    if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) {
                        trans = TransportWithCarryOverV21b(
                            prev_track->xi1, prev_track->xi2,
                            prev_track->tangent_t1, prev_track->tangent_t2, prev_track->normal,
                            patch.tangent_t1, patch.tangent_t2, patch.normal,
                            match_dist, topology_event,
                            vt1, vt2,
                            config.time_step,
                            config.tangential_stiffness,
                            config.tangential_damping,
                            config.friction_coefficient,
                            F_n.Length(),
                            is_newborn,
                            v1_history,
                            transport_cfg,
                            gate_stats,
                            is_v22,
                            is_v23,
                            is_v24,
                            patch.center - body_pos,
                            F_n
                        );
                    } else if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21) {
                        trans = TransportWithCarryOverV21(
                            prev_track->xi1, prev_track->xi2,
                            prev_track->tangent_t1, prev_track->tangent_t2, prev_track->normal,
                            patch.tangent_t1, patch.tangent_t2, patch.normal,
                            match_dist, topology_event,
                            vt1, vt2,
                            config.time_step,
                            config.tangential_stiffness,
                            config.tangential_damping,
                            config.friction_coefficient,
                            F_n.Length(),
                            is_newborn,
                            v1_history,
                            transport_cfg
                        );
                    } else {
                        trans = TransportWithCarryOverV1(
                            prev_track->xi1, prev_track->xi2,
                            prev_track->tangent_t1, prev_track->tangent_t2, prev_track->normal,
                            patch.tangent_t1, patch.tangent_t2, patch.normal,
                            match_dist, topology_event,
                            vt1, vt2,
                            config.time_step,
                            config.tangential_stiffness,
                            config.tangential_damping,
                            config.friction_coefficient,
                            F_n.Length(),
                            is_newborn,
                            v1_history,
                            transport_cfg
                        );
                    }
                    xi1 = trans.xi1_new;
                    xi2 = trans.xi2_new;
                    raw_trans_norm = trans.raw_transport_norm;
                    limited_trans_norm = trans.limited_transport_norm;
                    attenuation = trans.transport_attenuation;
                    carry_factor = trans.carry_over_factor;
                    adaptive_alpha = trans.adaptive_alpha;
                    clamp_hit = trans.transport_clamp_hit;
                    frame_rot_angle = trans.frame_rotation_angle;
                    if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) {
                        theta_rel_values.push_back(trans.theta_rel);
                        match_rel_values.push_back(trans.match_rel);
                        risk_combo_values.push_back(trans.risk_combo);
                        double v21_abs_gate = std::exp(-frame_rot_angle / 0.3) *
                                              std::exp(-match_dist / 0.05) *
                                              trans.w_topo *
                                              trans.w_slip;
                        v21_abs_gate_values.push_back(v21_abs_gate);
                        final_gate_values.push_back(trans.adaptive_alpha);
                        gate_w_theta_values.push_back(trans.w_theta);
                        gate_w_match_values.push_back(trans.w_match);
                        gate_w_topo_values.push_back(trans.w_topo);
                        gate_w_slip_values.push_back(trans.w_slip);
                        gate_w_tail_values.push_back(trans.w_tail);
                        direction_consistency_values.push_back(trans.direction_consistency);
                        gate_w_dir_values.push_back(trans.w_dir);
                        gate_before_direction_values.push_back(trans.gate_before_direction);
                        direction_pollution_values.push_back(trans.direction_pollution);
                        gate_w_pollution_values.push_back(trans.w_pollution);
                        gate_before_pollution_values.push_back(trans.gate_before_pollution);
                        basis_distortion_values.push_back(trans.basis_distortion);
                        gate_w_unified_values.push_back(trans.w_unified);
                        gate_before_unified_values.push_back(trans.gate_before_unified);

                        auto fill_gate_trace = [&](GateSampleTrace& trace, double score) {
                            trace.valid = true;
                            trace.frame_id = frame_offset;
                            trace.patch_index = static_cast<int>(pi);
                            trace.persistent_id = patch.matched_persistent_id;
                            trace.theta = frame_rot_angle;
                            trace.theta_ref = trans.theta_ref;
                            trace.theta_rel = trans.theta_rel;
                            trace.match_distance = match_dist;
                            trace.match_ref = trans.match_ref;
                            trace.match_rel = trans.match_rel;
                            trace.risk_combo = trans.risk_combo;
                            trace.xi_transport_norm = raw_trans_norm;
                            trace.delta_xi_norm = std::sqrt(std::pow(vt1 * config.time_step, 2) + std::pow(vt2 * config.time_step, 2));
                            trace.xi_prev_basis_t1 = trans.xi_prev_basis_t1;
                            trace.xi_prev_basis_t2 = trans.xi_prev_basis_t2;
                            trace.xi_cur_basis_t1 = trans.xi_cur_basis_t1;
                            trace.xi_cur_basis_t2 = trans.xi_cur_basis_t2;
                            trace.xi_transport_dir_t1 = trans.xi_transport_dir_t1;
                            trace.xi_transport_dir_t2 = trans.xi_transport_dir_t2;
                            trace.delta_xi_dir_t1 = trans.delta_xi_dir_t1;
                            trace.delta_xi_dir_t2 = trans.delta_xi_dir_t2;
                            trace.direction_consistency = trans.direction_consistency;
                            trace.gate_before_direction = trans.gate_before_direction;
                            trace.gate_after_direction = trans.adaptive_alpha;
                            trace.direction_pollution = trans.direction_pollution;
                            trace.gate_before_pollution = trans.gate_before_pollution;
                            trace.gate_after_pollution = trans.gate_after_pollution;
                            trace.basis_distortion = trans.basis_distortion;
                            trace.gate_before_unified = trans.gate_before_unified;
                            trace.gate_after_unified = trans.gate_after_unified;
                            trace.w_theta = trans.w_theta;
                            trace.w_match = trans.w_match;
                            trace.w_topo = trans.w_topo;
                            trace.w_slip = trans.w_slip;
                            trace.w_tail = trans.w_tail;
                            trace.w_dir = trans.w_dir;
                            trace.w_pollution = trans.w_pollution;
                            trace.w_unified = trans.w_unified;
                            trace.tau_hist_z = trans.tau_hist_z;
                            trace.tau_ref_z = trans.tau_ref_z;
                            trace.validity_gate = trans.adaptive_alpha;
                            trace.v21_abs_w_theta = std::exp(-frame_rot_angle / 0.3);
                            trace.v21_abs_w_match = std::exp(-match_dist / 0.05);
                            trace.v21_abs_gate = v21_abs_gate;
                            trace.raw_transport_norm = raw_trans_norm;
                            trace.limited_transport_norm = limited_trans_norm;
                            trace.score = score;
                        };

                        double high_score = is_v24 ? ((1.0 - trans.w_unified) *
                                                       trans.basis_distortion *
                                                       trans.direction_pollution *
                                                       std::max(1.0e-12, std::abs(trans.tau_hist_z)))
                                                   : (is_v23 ? ((1.0 - trans.w_pollution) *
                                                                trans.direction_pollution *
                                                                std::max(1.0e-12, std::abs(trans.tau_hist_z)))
                                                             : (is_v22 ? ((1.0 - trans.direction_consistency) * std::max(1.0e-12, raw_trans_norm))
                                                                       : trans.risk_combo));
                        if (high_score > gate_sample_trace.score) {
                            fill_gate_trace(gate_sample_trace, high_score);
                        }

                        double normal_score = raw_trans_norm;
                        bool normal_sample = is_v24 ? (trans.basis_distortion <= 0.05 &&
                                                        trans.direction_pollution <= 0.05 &&
                                                        trans.risk_combo <= 1.0)
                                                    : (is_v23 ? (trans.direction_pollution <= 0.05 && trans.risk_combo <= 1.0)
                                                              : (is_v22 ? (trans.direction_consistency >= 0.8 && trans.risk_combo <= 1.0)
                                                                        : (trans.risk_combo <= 1.0)));
                        if (normal_sample && normal_score > gate_normal_sample_trace.score) {
                            fill_gate_trace(gate_normal_sample_trace, normal_score);
                        }
                    }

                    // V1: increment already included in transport, no extra addition needed
                    // But if clamp is enabled, apply clamp to final value
                    double xi_norm_v1 = std::sqrt(xi1 * xi1 + xi2 * xi2);
                    
                    // Always update max tracker for trace sampling (regardless of clamp)
                    if (xi_norm_v1 > max_xi_tracker.max_xi) {
                        max_xi_tracker.max_xi = xi_norm_v1;
                        max_xi_tracker.max_frame = frame_offset;
                        max_xi_tracker.max_patch_idx = static_cast<int>(pi);
                        max_xi_tracker.max_persistent_id = patch.matched_persistent_id;
                        max_xi_tracker.xi_history_norm = (prev_track != nullptr) ? std::sqrt(prev_track->xi1 * prev_track->xi1 + prev_track->xi2 * prev_track->xi2) : 0.0;
                        max_xi_tracker.xi_proj_norm = raw_trans_norm;
                        max_xi_tracker.xi_transfer_norm = xi_norm_v1;
                        max_xi_tracker.xi_after_increment_norm = xi_norm_v1;
                        max_xi_tracker.xi_before_clamp_norm = xi_norm_v1;
                    }
                    
                    if (transport_cfg.enable_clamp) {
                        // Record global distribution (ALL samples, not just triggered)
                        all_xi_before_clamp_values.push_back(xi_norm_v1);
                        
                        if (xi_norm_v1 > transport_cfg.xi_max) {
                            double scale = transport_cfg.xi_max / xi_norm_v1;
                            xi1 *= scale;
                            xi2 *= scale;
                            clamp_hit = 1;
                            clamp_trigger_count++;
                            clamp_pre_norm_sum += xi_norm_v1;
                            clamp_post_norm_sum += std::sqrt(xi1 * xi1 + xi2 * xi2);
                            clamp_max_pre_norm = std::max(clamp_max_pre_norm, xi_norm_v1);
                            clamp_max_post_norm = std::max(clamp_max_post_norm, std::sqrt(xi1 * xi1 + xi2 * xi2));
                            clamp_sample_count++;
                            
                            // Update max tracker with post-clamp info
                            if (max_xi_tracker.max_frame == frame_offset && max_xi_tracker.max_persistent_id == patch.matched_persistent_id) {
                                max_xi_tracker.xi_after_clamp_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
                            }
                        } else {
                            // Update max tracker with post-clamp info (no change)
                            if (max_xi_tracker.max_frame == frame_offset && max_xi_tracker.max_persistent_id == patch.matched_persistent_id) {
                                max_xi_tracker.xi_after_clamp_norm = xi_norm_v1;
                            }
                        }
                    }
                } else {
                    // Original: transport + incremental update
                    TransportResult trans = TransportWithCarryOver(
                        prev_track->xi1, prev_track->xi2,
                        prev_track->tangent_t1, prev_track->tangent_t2, prev_track->normal,
                        patch.tangent_t1, patch.tangent_t2, patch.normal,
                        match_dist, topology_event,
                        transport_cfg
                    );
                    xi1 = trans.xi1_new;
                    xi2 = trans.xi2_new;
                    raw_trans_norm = trans.raw_transport_norm;
                    limited_trans_norm = trans.limited_transport_norm;
                    attenuation = trans.transport_attenuation;
                    carry_factor = trans.carry_over_factor;
                    adaptive_alpha = trans.adaptive_alpha;
                    clamp_hit = trans.transport_clamp_hit;
                    frame_rot_angle = trans.frame_rotation_angle;

                    // Record clamp diagnostics BEFORE increment (this is already clamped by TransportWithCarryOver)
                    // For accurate diagnostics, we need to record the value AFTER clamp but BEFORE increment
                    if (transport_cfg.enable_clamp && clamp_hit) {
                        double xi_norm_after_clamp = std::sqrt(xi1 * xi1 + xi2 * xi2);
                        clamp_trigger_count++;
                        // Since TransportWithCarryOver already clamped, we record the clamped value
                        // pre_norm = value after internal clamp (from TransportWithCarryOver)
                        // post_norm = same as pre_norm here (no additional clamp applied yet)
                        clamp_pre_norm_sum += xi_norm_after_clamp;
                        clamp_post_norm_sum += xi_norm_after_clamp;
                        clamp_max_pre_norm = std::max(clamp_max_pre_norm, xi_norm_after_clamp);
                        clamp_max_post_norm = std::max(clamp_max_post_norm, xi_norm_after_clamp);
                        clamp_sample_count++;
                    }

                    // Update tangential displacement with current step increment
                    xi1 += vt1 * config.time_step;
                    xi2 += vt2 * config.time_step;

                    // After increment: record xi norm for trace sampling (regardless of clamp)
                    double xi_norm_after_increment = std::sqrt(xi1 * xi1 + xi2 * xi2);
                    
                    // Always update max tracker for trace sampling (regardless of clamp)
                    if (xi_norm_after_increment > max_xi_tracker.max_xi) {
                        max_xi_tracker.max_xi = xi_norm_after_increment;
                        max_xi_tracker.max_frame = frame_offset;
                        max_xi_tracker.max_patch_idx = static_cast<int>(pi);
                        max_xi_tracker.max_persistent_id = patch.matched_persistent_id;
                        max_xi_tracker.xi_history_norm = (prev_track != nullptr) ? std::sqrt(prev_track->xi1 * prev_track->xi1 + prev_track->xi2 * prev_track->xi2) : 0.0;
                        max_xi_tracker.xi_proj_norm = raw_trans_norm;
                        max_xi_tracker.xi_transfer_norm = std::sqrt((xi1 - vt1 * config.time_step) * (xi1 - vt1 * config.time_step) + (xi2 - vt2 * config.time_step) * (xi2 - vt2 * config.time_step));
                        max_xi_tracker.xi_after_increment_norm = xi_norm_after_increment;
                        max_xi_tracker.xi_before_clamp_norm = xi_norm_after_increment;
                    }
                    
                    if (transport_cfg.enable_clamp) {
                        // Record global distribution (ALL samples, not just triggered)
                        all_xi_before_clamp_values.push_back(xi_norm_after_increment);
                        
                        // Only count as trigger if clamp actually modifies the value
                        if (xi_norm_after_increment > transport_cfg.xi_max) {
                            double xi_before = xi_norm_after_increment;
                            double scale = transport_cfg.xi_max / xi_norm_after_increment;
                            xi1 *= scale;
                            xi2 *= scale;
                            double xi_after = std::sqrt(xi1 * xi1 + xi2 * xi2);
                            
                            clamp_trigger_count++;
                            clamp_pre_norm_sum += xi_before;
                            clamp_post_norm_sum += xi_after;
                            clamp_max_pre_norm = std::max(clamp_max_pre_norm, xi_before);
                            clamp_max_post_norm = std::max(clamp_max_post_norm, xi_after);
                            clamp_sample_count++;
                            
                            // Update max tracker with post-clamp info
                            if (max_xi_tracker.max_frame == frame_offset && max_xi_tracker.max_persistent_id == patch.matched_persistent_id) {
                                max_xi_tracker.xi_after_clamp_norm = xi_after;
                            }
                        } else {
                            // Update max tracker with post-clamp info (no change)
                            if (max_xi_tracker.max_frame == frame_offset && max_xi_tracker.max_persistent_id == patch.matched_persistent_id) {
                                max_xi_tracker.xi_after_clamp_norm = xi_norm_after_increment;
                            }
                        }
                    }
                }

                for (auto& track : current_tracks) {
                    if (track.persistent_id == patch.matched_persistent_id) {
                        track.stick_slip_state = prev_track->stick_slip_state;
                        track.sum_attenuation += attenuation;
                        track.sum_carry_over += carry_factor;
                        track.total_transport_clamp_hits += clamp_hit;
                        track.accumulated_transport_norm += limited_trans_norm;
                        break;
                    }
                }
            } else {
                // Newborn patch: initialize
                xi1 = vt1 * config.time_step;
                xi2 = vt2 * config.time_step;
            }

            // Trial tangential force in local frame. Elastic history and damping are
            // separated so the stored state remains the recoverable elastic part.
            double ft1_el_trial = -config.tangential_stiffness * xi1;
            double ft2_el_trial = -config.tangential_stiffness * xi2;
            double ft_el_trial_mag = std::sqrt(ft1_el_trial * ft1_el_trial + ft2_el_trial * ft2_el_trial);
            double ft1_damp = -config.tangential_damping * vt1;
            double ft2_damp = -config.tangential_damping * vt2;
            double ft1_trial = ft1_el_trial + ft1_damp;
            double ft2_trial = ft2_el_trial + ft2_damp;
            double ft_trial_mag = std::sqrt(ft1_trial * ft1_trial + ft2_trial * ft2_trial);

            // Coulomb return mapping: project the elastic trial force first, then
            // project the elastic+damping force to the final friction disk.
            double fn_mag = F_n.Length();
            double friction_limit = config.friction_coefficient * fn_mag;
            StickSlipState state;
            ChVector3d F_t;
            double ft_mag;

            if (friction_limit <= 1.0e-12) {
                F_t = ChVector3d(0, 0, 0);
                ft_mag = 0.0;
                xi1 = 0.0;
                xi2 = 0.0;
                state = StickSlipState::Slip;
            } else {
                bool elastic_slip = ft_el_trial_mag > friction_limit;
                double ft1_el = ft1_el_trial;
                double ft2_el = ft2_el_trial;
                if (elastic_slip) {
                    double elastic_scale = friction_limit / ft_el_trial_mag;
                    ft1_el *= elastic_scale;
                    ft2_el *= elastic_scale;
                    xi1 = -ft1_el / config.tangential_stiffness;
                    xi2 = -ft2_el / config.tangential_stiffness;
                }

                double ft1_combined = ft1_el + ft1_damp;
                double ft2_combined = ft2_el + ft2_damp;
                double ft_combined_mag = std::sqrt(ft1_combined * ft1_combined + ft2_combined * ft2_combined);
                bool damping_projection = ft_combined_mag > friction_limit;

                if (!elastic_slip && !damping_projection) {
                    F_t = ft1_combined * patch.tangent_t1 + ft2_combined * patch.tangent_t2;
                    ft_mag = ft_combined_mag;
                    state = StickSlipState::Stick;
                } else {
                    if (damping_projection && ft_combined_mag > 1.0e-12) {
                        double combined_scale = friction_limit / ft_combined_mag;
                        ft1_combined *= combined_scale;
                        ft2_combined *= combined_scale;
                    } else {
                        ft1_combined = ft1_el;
                        ft2_combined = ft2_el;
                    }
                    F_t = ft1_combined * patch.tangent_t1 + ft2_combined * patch.tangent_t2;
                    ft_mag = F_t.Length();
                    state = StickSlipState::Slip;
                }
            }

            patch.force = F_n + F_t;
            patch.stick_slip_state = state;
            patch.normal_force_magnitude = fn_mag;
            patch.tangential_force_magnitude = ft_mag;
            patch.torque = (patch.center - body_pos).Cross(patch.force);
            patch.tangential_displacement_local_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);

            if ((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) &&
                patch.matched_persistent_id >= 0) {
                auto update_gate_trace_force_fields = [&](GateSampleTrace& trace) {
                    if (trace.frame_id == frame_offset &&
                        trace.patch_index == static_cast<int>(pi) &&
                        trace.persistent_id == patch.matched_persistent_id) {
                        trace.ft_trial_mag = ft_trial_mag;
                        trace.friction_limit = friction_limit;
                        trace.stick_or_slip = (state == StickSlipState::Stick) ? "stick" : "slip";
                        trace.ft_final_mag = ft_mag;
                    }
                };
                update_gate_trace_force_fields(gate_sample_trace);
                update_gate_trace_force_fields(gate_normal_sample_trace);
            }

            // -- Update max_xi_tracker with force fields (AFTER force computation) --
            if (max_xi_tracker.max_frame == frame_offset && max_xi_tracker.max_patch_idx == static_cast<int>(pi)) {
                max_xi_tracker.ft_trial_mag = ft_trial_mag;
                max_xi_tracker.friction_limit = friction_limit;
                max_xi_tracker.stick_or_slip = (state == StickSlipState::Stick) ? "stick" : "slip";
                max_xi_tracker.ft_final_mag = ft_mag;
                max_xi_tracker.fn_mag = fn_mag;
                max_xi_tracker.xi_before_return_mapping_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
                max_xi_tracker.xi_after_return_mapping_norm = std::sqrt(xi1 * xi1 + xi2 * xi2);
                max_xi_tracker.is_valid = true;
            }

            // Update track state
            for (auto& track : current_tracks) {
                if (track.persistent_id == patch.matched_persistent_id) {
                    track.xi1 = xi1;
                    track.xi2 = xi2;
                    track.stick_slip_state = state;
                    track.normal_force_magnitude = fn_mag;
                    track.tangential_force_magnitude = ft_mag;
                    if (state == StickSlipState::Stick) track.stick_steps++;
                    else track.slip_steps++;
                    if (prev_track != nullptr && prev_track->stick_slip_state == StickSlipState::Stick && state == StickSlipState::Slip) {
                        track.stick_to_slip_transitions++;
                    }
                    break;
                }
            }

            // Update V1/V2.1 tangential history state
            if ((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV1 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                 transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) &&
                patch.matched_persistent_id >= 0) {
                ChVector3d xi_elastic_world_new = xi1 * patch.tangent_t1 + xi2 * patch.tangent_t2;
                TangentialHistoryState new_hist;
                new_hist.persistent_id = patch.matched_persistent_id;
                new_hist.xi_elastic_world = xi_elastic_world_new;
                new_hist.normal_prev = patch.normal;
                new_hist.t1_prev = patch.tangent_t1;
                new_hist.t2_prev = patch.tangent_t2;
                new_hist.match_quality_prev = std::exp(-match_dist / (config.dyn_body_type == 0 ? config.dyn_body_radius : config.dyn_body_size));
                new_hist.area_prev = patch.total_area;
                new_hist.topology_state_prev = 0;
                new_hist.age = (prev_track != nullptr) ? prev_track->age_in_steps : 0;
                new_hist.valid = true;
                tangential_history[patch.matched_persistent_id] = new_hist;
            }

            total_force += patch.force;
            total_torque += patch.torque;

            // Record diagnostics
            theta_values.push_back(frame_rot_angle);
            match_distance_values.push_back(match_dist);
            if (transport_cfg.carry_strategy == CarryOverStrategy::Adaptive || 
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV1 ||
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21 ||
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
                transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) {
                adaptive_alpha_values.push_back(adaptive_alpha);
            }

            // -- Record patch trace for target frame --
            if (frame_offset == trace_target_frame && patch.matched_persistent_id >= 0 && patch_traces.size() < 5) {
                // Select patches with highest tangential force for tracing
                PatchTrace trace;
                trace.frame_id = frame_offset;
                trace.persistent_id = patch.matched_persistent_id;
                trace.config_name = mode_name;
                trace.xi_max = transport_cfg.xi_max;
                trace.clamp_enabled = transport_cfg.enable_clamp;
                trace.clamp_condition_value = std::sqrt(xi1 * xi1 + xi2 * xi2);
                trace.clamp_triggered = (clamp_hit > 0);
                
                // Estimate xi_history_norm from prev_track if available
                trace.xi_history_norm = (prev_track != nullptr) ? std::sqrt(prev_track->xi1 * prev_track->xi1 + prev_track->xi2 * prev_track->xi2) : 0.0;
                trace.xi_proj_norm = raw_trans_norm;
                trace.xi_transfer_norm = (prev_track != nullptr) ? std::sqrt(xi1 * xi1 + xi2 * xi2) : 0.0;
                trace.xi_before_increment_norm = (prev_track != nullptr) ? std::sqrt(xi1 * xi1 + xi2 * xi2) : 0.0;
                trace.xi_after_increment_norm = patch.tangential_displacement_local_norm;
                trace.xi_before_clamp_norm = trace.xi_after_increment_norm;
                trace.xi_after_clamp_norm = patch.tangential_displacement_local_norm;
                trace.xi_before_return_mapping_norm = trace.xi_after_clamp_norm;
                trace.xi_after_return_mapping_norm = patch.tangential_displacement_local_norm;
                trace.ft_trial_mag = ft_trial_mag;
                trace.friction_limit = friction_limit;
                trace.stick_or_slip = (state == StickSlipState::Stick) ? "stick" : "slip";
                trace.ft_final_mag = ft_mag;
                
                patch_traces.push_back(trace);
            }
        }

        active_tracks = current_tracks;

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
        int current_patch_count = patches.size();
        patch_counts.push_back(current_patch_count);
        torque_x_values.push_back(total_torque.x());
        torque_y_values.push_back(total_torque.y());
        torque_z_values.push_back(total_torque.z());

        total_force_y += current_force_y;
        force_sample_count++;

        // Tangential force stats
        double ft_approx = std::sqrt(total_force.x() * total_force.x() + total_force.z() * total_force.z());
        double fn_approx = std::abs(total_force.y());
        tangential_force_norms.push_back(ft_approx);
        sum_tangential_norm += ft_approx;
        if (ft_approx > max_tangential_norm) max_tangential_norm = ft_approx;
        if (fn_approx > 1e-12) {
            sum_tangential_ratio += ft_approx / fn_approx;
        }

        // Stick-slip and transport stats
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
    result.avg_patch_count = patch_counts.empty() ? 0.0 :
        static_cast<double>(std::accumulate(patch_counts.begin(), patch_counts.end(), 0)) / patch_counts.size();

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
    result.theta_p50 = Percentile(theta_values, 0.50);
    result.theta_p95 = Percentile(theta_values, 0.95);
    result.match_distance_p50 = Percentile(match_distance_values, 0.50);
    result.match_distance_p95 = Percentile(match_distance_values, 0.95);
    result.avg_theta_rel = MeanValue(theta_rel_values);
    result.theta_rel_p50 = Percentile(theta_rel_values, 0.50);
    result.theta_rel_p95 = Percentile(theta_rel_values, 0.95);
    result.theta_rel_p99 = Percentile(theta_rel_values, 0.99);
    result.max_theta_rel = theta_rel_values.empty() ? 0.0 : *std::max_element(theta_rel_values.begin(), theta_rel_values.end());
    result.avg_match_rel = MeanValue(match_rel_values);
    result.match_rel_p50 = Percentile(match_rel_values, 0.50);
    result.match_rel_p95 = Percentile(match_rel_values, 0.95);
    result.match_rel_p99 = Percentile(match_rel_values, 0.99);
    result.max_match_rel = match_rel_values.empty() ? 0.0 : *std::max_element(match_rel_values.begin(), match_rel_values.end());
    result.avg_risk_combo = MeanValue(risk_combo_values);
    result.risk_combo_p50 = Percentile(risk_combo_values, 0.50);
    result.risk_combo_p95 = Percentile(risk_combo_values, 0.95);
    result.risk_combo_p99 = Percentile(risk_combo_values, 0.99);
    result.max_risk_combo = risk_combo_values.empty() ? 0.0 : *std::max_element(risk_combo_values.begin(), risk_combo_values.end());
    result.avg_gate_w_theta = MeanValue(gate_w_theta_values);
    result.avg_gate_w_match = MeanValue(gate_w_match_values);
    result.avg_gate_w_topo = MeanValue(gate_w_topo_values);
    result.avg_gate_w_slip = MeanValue(gate_w_slip_values);
    result.avg_gate_w_tail = MeanValue(gate_w_tail_values);
    result.avg_gate_w_dir = MeanValue(gate_w_dir_values);
    result.avg_direction_consistency = MeanValue(direction_consistency_values);
    result.direction_consistency_p50 = Percentile(direction_consistency_values, 0.50);
    result.direction_consistency_p95 = Percentile(direction_consistency_values, 0.95);
    result.direction_consistency_p99 = Percentile(direction_consistency_values, 0.99);
    result.direction_bad_gate_before_avg = MeanGateForDirectionMax(direction_consistency_values, gate_before_direction_values, 0.50);
    result.direction_bad_gate_after_avg = MeanGateForDirectionMax(direction_consistency_values, final_gate_values, 0.50);
    result.avg_basis_distortion = MeanValue(basis_distortion_values);
    result.basis_distortion_p50 = Percentile(basis_distortion_values, 0.50);
    result.basis_distortion_p95 = Percentile(basis_distortion_values, 0.95);
    result.basis_distortion_p99 = Percentile(basis_distortion_values, 0.99);
    result.avg_gate_w_unified = MeanValue(gate_w_unified_values);
    result.basis_gate_before_avg = MeanGateForBasisMin(basis_distortion_values, gate_before_unified_values, result.basis_distortion_p95);
    result.basis_gate_after_avg = MeanGateForBasisMin(basis_distortion_values, final_gate_values, result.basis_distortion_p95);
    result.avg_direction_pollution = MeanValue(direction_pollution_values);
    result.direction_pollution_p50 = Percentile(direction_pollution_values, 0.50);
    result.direction_pollution_p95 = Percentile(direction_pollution_values, 0.95);
    result.direction_pollution_p99 = Percentile(direction_pollution_values, 0.99);
    result.avg_gate_w_pollution = MeanValue(gate_w_pollution_values);
    result.pollution_gate_before_avg = MeanGateForPollutionMin(direction_pollution_values, gate_before_pollution_values, result.direction_pollution_p95);
    result.pollution_gate_after_avg = MeanGateForPollutionMin(direction_pollution_values, final_gate_values, result.direction_pollution_p95);
    result.high_risk_gate_avg = MeanGateForRiskRange(risk_combo_values, final_gate_values, result.risk_combo_p95, std::numeric_limits<double>::infinity());
    result.high_risk_v21_gate_avg = MeanGateForRiskRange(risk_combo_values, v21_abs_gate_values, result.risk_combo_p95, std::numeric_limits<double>::infinity());
    result.median_risk_gate_avg = MeanGateForRiskRange(risk_combo_values, final_gate_values, result.risk_combo_p50 * 0.95, result.risk_combo_p50 * 1.05);
    result.median_risk_v21_gate_avg = MeanGateForRiskRange(risk_combo_values, v21_abs_gate_values, result.risk_combo_p50 * 0.95, result.risk_combo_p50 * 1.05);
    result.total_birth_events = total_birth_events;

    result.total_stick_steps = total_stick_steps;
    result.total_slip_steps = total_slip_steps;
    result.total_stick_to_slip_transitions = total_stick_to_slip_transitions;
    result.total_transport_clamp_hits = total_transport_clamp_hits;
    result.avg_frame_rotation_angle = result.avg_theta;

    // -- Clamp diagnostics --
    result.clamp_trigger_count = clamp_trigger_count;
    result.clamp_pre_norm_avg = clamp_sample_count > 0 ? clamp_pre_norm_sum / clamp_sample_count : 0.0;
    result.clamp_post_norm_avg = clamp_sample_count > 0 ? clamp_post_norm_sum / clamp_sample_count : 0.0;
    result.clamp_max_pre_norm = clamp_max_pre_norm;
    result.clamp_max_post_norm = clamp_max_post_norm;

    // -- Compute global xi_before_clamp distribution --
    double xi_before_clamp_mean = 0.0;
    double xi_before_clamp_p50 = 0.0;
    double xi_before_clamp_p95 = 0.0;
    double xi_before_clamp_p99 = 0.0;
    double xi_before_clamp_max = 0.0;
    
    if (!all_xi_before_clamp_values.empty()) {
        // Sort for percentiles
        std::vector<double> sorted_xi = all_xi_before_clamp_values;
        std::sort(sorted_xi.begin(), sorted_xi.end());
        
        size_t n = sorted_xi.size();
        xi_before_clamp_max = sorted_xi.back();
        xi_before_clamp_mean = std::accumulate(sorted_xi.begin(), sorted_xi.end(), 0.0) / n;
        xi_before_clamp_p50 = sorted_xi[n / 2];
        xi_before_clamp_p95 = sorted_xi[static_cast<size_t>(n * 0.95)];
        xi_before_clamp_p99 = sorted_xi[static_cast<size_t>(n * 0.99)];
    }

    // -- Write global distribution file --
    {
        std::string case_name = (config.dyn_body_type == 0) ? "Case_A" : "Case_B";
        std::string dist_file = out_dir + "/sdf_patch_adaptive_carryover_xi_dist_" + case_name + "_" + mode_name + ".csv";
        std::ofstream df(dist_file);
        df << "config,case_name,sample_count,xi_before_clamp_mean,xi_before_clamp_p50,xi_before_clamp_p95,xi_before_clamp_p99,xi_before_clamp_max,xi_max,clamp_trigger_count" << std::endl;
        df << mode_name << "," << case_name << ","
           << all_xi_before_clamp_values.size() << ","
           << std::fixed << std::setprecision(8)
           << xi_before_clamp_mean << ","
           << xi_before_clamp_p50 << ","
           << xi_before_clamp_p95 << ","
           << xi_before_clamp_p99 << ","
           << xi_before_clamp_max << ","
           << transport_cfg.xi_max << ","
           << clamp_trigger_count << std::endl;
        df.close();
    }

    // -- Write max xi_before_clamp trace file --
    {
        std::string case_name = (config.dyn_body_type == 0) ? "Case_A" : "Case_B";
        std::string trace_file = out_dir + "/sdf_patch_adaptive_carryover_max_trace_" + case_name + "_" + mode_name + ".csv";
        std::ofstream tf(trace_file);
        tf << "case_name,config_name,frame_id,patch_index,persistent_id,xi_max,clamp_enabled,"
           << "clamp_condition_variable_name,clamp_condition_value,clamp_triggered,"
           << "xi_history_norm,xi_proj_norm,xi_transfer_norm,xi_after_increment_norm,xi_before_clamp_norm,xi_after_clamp_norm,"
           << "xi_before_return_mapping_norm,xi_after_return_mapping_norm,"
           << "ft_trial_mag,friction_limit,stick_or_slip,ft_final_mag,"
           << "sample_rank_by_xi_before_clamp,is_global_max_sample" << std::endl;
        
        std::string cond_var_name = "xi_norm_after_increment";
        
        // Determine sample rank (1 = global max)
        int sample_rank = (max_xi_tracker.max_frame >= 0) ? 1 : 0;
        std::string is_global_max = (sample_rank == 1) ? "yes" : "no";
        
        tf << case_name << "," << mode_name << ","
           << max_xi_tracker.max_frame << "," << max_xi_tracker.max_patch_idx << ","
           << max_xi_tracker.max_persistent_id << ","
           << std::fixed << std::setprecision(8)
           << transport_cfg.xi_max << "," << (transport_cfg.enable_clamp ? 1 : 0) << ","
           << cond_var_name << "," << max_xi_tracker.max_xi << ","
           << (max_xi_tracker.max_xi > transport_cfg.xi_max ? 1 : 0) << ","
           << max_xi_tracker.xi_history_norm << "," << max_xi_tracker.xi_proj_norm << ","
           << max_xi_tracker.xi_transfer_norm << "," << max_xi_tracker.xi_after_increment_norm << ","
           << max_xi_tracker.xi_before_clamp_norm << "," << max_xi_tracker.xi_after_clamp_norm << ","
           << max_xi_tracker.xi_before_return_mapping_norm << "," << max_xi_tracker.xi_after_return_mapping_norm << ","
           << max_xi_tracker.ft_trial_mag << "," << max_xi_tracker.friction_limit << ","
           << max_xi_tracker.stick_or_slip << "," << max_xi_tracker.ft_final_mag << ","
           << sample_rank << "," << is_global_max << std::endl;
        tf.close();
    }

    // -- Write V2.1b/V2.2 gate representative samples --
    if (transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV21b ||
        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22 ||
        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23 ||
        transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) {
        std::string case_name = (config.dyn_body_type == 0) ? "Case_A" : "Case_B";
        std::string gate_file = out_dir + "/sdf_patch_adaptive_carryover_gate_sample_" + case_name + "_" + mode_name + ".csv";
        std::ofstream gf(gate_file);
        gf << "trace_type,case_name,config_name,frame_id,patch_index,persistent_id,theta,theta_ref,theta_rel,match_distance,match_ref,match_rel,"
           << "risk_combo,xi_transport_norm,delta_xi_norm,xi_prev_basis_t1,xi_prev_basis_t2,xi_cur_basis_t1,xi_cur_basis_t2,"
           << "xi_transport_dir_t1,xi_transport_dir_t2,delta_xi_dir_t1,delta_xi_dir_t2,"
           << "basis_distortion,direction_consistency,direction_pollution,gate_before_direction,gate_after_direction,"
           << "gate_before_pollution,gate_after_pollution,gate_before_unified,gate_after_unified,"
           << "tau_hist_z,tau_ref_z,w_theta,w_match,w_topo,w_slip,w_tail,w_dir,w_pollution,w_unified,validity_gate,v21_abs_w_theta,v21_abs_w_match,v21_abs_gate,"
           << "raw_transport_norm,limited_transport_norm,ft_trial_mag,friction_limit,stick_or_slip,ft_final_mag,score" << std::endl;
        auto write_gate_trace = [&](const std::string& trace_type, const GateSampleTrace& trace) {
            gf << trace_type << "," << case_name << "," << mode_name << ","
               << trace.frame_id << "," << trace.patch_index << "," << trace.persistent_id << ","
               << std::fixed << std::setprecision(8)
               << trace.theta << "," << trace.theta_ref << "," << trace.theta_rel << ","
               << trace.match_distance << "," << trace.match_ref << "," << trace.match_rel << ","
               << trace.risk_combo << ","
               << trace.xi_transport_norm << "," << trace.delta_xi_norm << ","
               << trace.xi_prev_basis_t1 << "," << trace.xi_prev_basis_t2 << ","
               << trace.xi_cur_basis_t1 << "," << trace.xi_cur_basis_t2 << ","
               << trace.xi_transport_dir_t1 << "," << trace.xi_transport_dir_t2 << ","
               << trace.delta_xi_dir_t1 << "," << trace.delta_xi_dir_t2 << ","
               << trace.basis_distortion << "," << trace.direction_consistency << "," << trace.direction_pollution << ","
               << trace.gate_before_direction << "," << trace.gate_after_direction << ","
               << trace.gate_before_pollution << "," << trace.gate_after_pollution << ","
               << trace.gate_before_unified << "," << trace.gate_after_unified << ","
               << trace.tau_hist_z << "," << trace.tau_ref_z << ","
               << trace.w_theta << "," << trace.w_match << "," << trace.w_topo << ","
               << trace.w_slip << "," << trace.w_tail << "," << trace.w_dir << "," << trace.w_pollution << "," << trace.w_unified << "," << trace.validity_gate << ","
               << trace.v21_abs_w_theta << "," << trace.v21_abs_w_match << "," << trace.v21_abs_gate << ","
               << trace.raw_transport_norm << "," << trace.limited_transport_norm << ","
               << trace.ft_trial_mag << "," << trace.friction_limit << "," << trace.stick_or_slip << "," << trace.ft_final_mag << ","
               << trace.score << std::endl;
        };
        write_gate_trace((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV24) ? "high_basis_distortion_pollution" :
                         ((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV23) ? "high_direction_pollution" :
                          ((transport_cfg.carry_strategy == CarryOverStrategy::AdaptiveV22) ? "high_direction_inconsistency" : "high_risk")),
                         gate_sample_trace);
        write_gate_trace("normal_retention", gate_normal_sample_trace);
        gf.close();
    }

    size_t last_10_percent = pos_y_values.size() * 9 / 10;
    double avg_last = 0.0;
    for (size_t i = last_10_percent; i < pos_y_values.size(); i++) {
        avg_last += pos_y_values[i];
    }
    avg_last /= (pos_y_values.size() - last_10_percent);
    result.is_stable = std::abs(avg_last - result.expected_equilibrium_y) < 0.06;

    return result;
}

static ChVector3d RotateAroundAxis(const ChVector3d& value, const ChVector3d& axis, double angle) {
    ChVector3d a = SafeNormalize(axis, ChVector3d(0, 0, 1));
    double c = std::cos(angle);
    double s = std::sin(angle);
    return value * c + a.Cross(value) * s + a * (a.Dot(value) * (1.0 - c));
}

static void WriteTangentialTheoryDiagnostics(const std::string& out_dir) {
    const ChVector3d n0(0, 1, 0);
    const ChVector3d xi0(0.01, 0, 0);
    const ChVector3d axis(0, 0, 1);
    const ChVector3d q_axis(1, 2, 3);
    const double q_angle = 0.73;

    std::string csv_path = out_dir + "/sdf_patch_tangential_transport_theory.csv";
    std::ofstream csv(csv_path);
    csv << "angle_rad,projection_norm,minrot_norm,expected_norm,"
        << "projection_parallel_error,minrot_parallel_error,minrot_objectivity_error" << std::endl;

    double max_projection_error = 0.0;
    double max_minrot_error = 0.0;
    double max_objectivity_error = 0.0;

    for (int i = 0; i <= 24; i++) {
        double angle = (0.5 * M_PI) * static_cast<double>(i) / 24.0;
        ChVector3d n1 = RotateAroundAxis(n0, axis, angle);
        ChVector3d expected = RotateAroundAxis(xi0, axis, angle);
        ChVector3d projection_transport = ProjectToTangent(xi0, n1);
        ChVector3d minrot_transport = TransportElasticStateToCurrentTangent(xi0, n0, n1);

        double projection_error = (projection_transport - expected).Length();
        double minrot_error = (minrot_transport - expected).Length();

        ChVector3d q_n0 = RotateAroundAxis(n0, q_axis, q_angle);
        ChVector3d q_n1 = RotateAroundAxis(n1, q_axis, q_angle);
        ChVector3d q_xi0 = RotateAroundAxis(xi0, q_axis, q_angle);
        ChVector3d q_transport_expected = RotateAroundAxis(minrot_transport, q_axis, q_angle);
        ChVector3d q_transport_actual = TransportElasticStateToCurrentTangent(q_xi0, q_n0, q_n1);
        double objectivity_error = (q_transport_actual - q_transport_expected).Length();

        max_projection_error = std::max(max_projection_error, projection_error);
        max_minrot_error = std::max(max_minrot_error, minrot_error);
        max_objectivity_error = std::max(max_objectivity_error, objectivity_error);

        csv << std::fixed << std::setprecision(10)
            << angle << ","
            << projection_transport.Length() << ","
            << minrot_transport.Length() << ","
            << expected.Length() << ","
            << projection_error << ","
            << minrot_error << ","
            << objectivity_error << std::endl;
    }
    csv.close();

    std::ofstream summary(out_dir + "/sdf_patch_tangential_theory_summary.csv");
    summary << "metric,value" << std::endl;
    summary << std::fixed << std::setprecision(10)
            << "max_projection_parallel_error," << max_projection_error << std::endl
            << "max_minrot_parallel_error," << max_minrot_error << std::endl
            << "max_minrot_objectivity_error," << max_objectivity_error << std::endl;
    summary.close();
}

// =============================================================================
// Main simulation
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 13: Adaptive Carry-Over v1 ===" << std::endl;

    // -- Output directory --
    std::string project_root = GetProjectRoot();
    std::string out_dir = project_root + "/out/milestone_13";
    EnsureDir(out_dir);
    WriteTangentialTheoryDiagnostics(out_dir);

    // ========================================================================
    // Baseline Reference
    // ========================================================================
    std::cout << "\n=== Baseline Reference (from Milestone 12) ===" << std::endl;
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

    // Adaptive carry-over parameters
    double adaptive_alpha_base = 0.8;     // Base alpha when theta=0 (higher than fixed 0.5)
    double adaptive_theta_scale = 0.3;     // V1/V2.1 absolute theta decay scale (radians)
    double adaptive_match_scale = 0.05;    // V1/V2.1 absolute match distance decay scale (meters)
    double v21b_theta_rel_scale = 0.7;     // V2.1b relative theta decay scale (dimensionless)
    double v21b_match_rel_scale = 8.0;     // V2.1b relative match distance decay scale (dimensionless)

    TransportConfig cfg_A = {"A-only", true, false, CarryOverStrategy::Fixed, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale, 0.0};
    TransportConfig cfg_C_fixed = {"C-fixed", false, false, CarryOverStrategy::Fixed, 0.01, 5.0, 0.5, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale, 0.0};
    TransportConfig cfg_C_adaptive_v1 = {"C-adaptive-v1", false, false, CarryOverStrategy::AdaptiveV1, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale, 0.0};
    TransportConfig cfg_C_adaptive_v21 = {"C-adaptive-v2.1", false, false, CarryOverStrategy::AdaptiveV21, 0.01, 5.0, 0.0, adaptive_alpha_base, adaptive_theta_scale, adaptive_match_scale, 0.5};
    TransportConfig cfg_C_adaptive_v21b = {"C-adaptive-v2.1b", false, false, CarryOverStrategy::AdaptiveV21b, 0.01, 5.0, 0.0, adaptive_alpha_base, v21b_theta_rel_scale, v21b_match_rel_scale, 0.05};
    TransportConfig cfg_C_adaptive_v22 = {"C-adaptive-v2.2", false, false, CarryOverStrategy::AdaptiveV22, 0.01, 5.0, 0.0, adaptive_alpha_base, v21b_theta_rel_scale, v21b_match_rel_scale, 0.05};
    TransportConfig cfg_C_adaptive_v23 = {"C-adaptive-v2.3", false, false, CarryOverStrategy::AdaptiveV23, 0.01, 5.0, 0.0, adaptive_alpha_base, v21b_theta_rel_scale, v21b_match_rel_scale, 0.05};
    TransportConfig cfg_C_adaptive_v24 = {"C-adaptive-v2.4", false, false, CarryOverStrategy::AdaptiveV24, 0.01, 5.0, 0.0, adaptive_alpha_base, v21b_theta_rel_scale, v21b_match_rel_scale, 0.05};

    std::vector<TransportConfig> configs = {
        cfg_A,
        cfg_C_fixed,
        cfg_C_adaptive_v1,
        cfg_C_adaptive_v21b,
        cfg_C_adaptive_v23,
        cfg_C_adaptive_v24
    };

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

            RunResult result = RunCaseWithConfig(config, configs[ci], configs[ci].name, out_dir);
            case_results.push_back(result);

            std::cout << "    Final Y: " << std::fixed << std::setprecision(6) << result.final_y << " m" << std::endl;
            std::cout << "    Y Error: " << result.y_error << " m" << std::endl;
            std::cout << "    Force Std Dev: " << result.force_std_dev << " N" << std::endl;
            std::cout << "    Avg Torque Z: " << result.avg_torque_z << " Nm" << std::endl;
            std::cout << "    Avg Carry-Over: " << result.avg_carry_over_factor << std::endl;
            std::cout << "    Avg Adaptive Alpha: " << result.avg_adaptive_alpha << std::endl;
            std::cout << "    Min/Max Adaptive Alpha: " << result.min_adaptive_alpha << " / " << result.max_adaptive_alpha << std::endl;
            std::cout << "    Avg Theta: " << result.avg_theta << " rad" << std::endl;
            std::cout << "    Max Theta: " << result.max_theta << " rad" << std::endl;
            std::cout << "    Avg Match Distance: " << result.avg_match_distance << " m" << std::endl;
            std::cout << "    Avg Theta Rel / Match Rel: " << result.avg_theta_rel << " / " << result.avg_match_rel << std::endl;
            std::cout << "    Avg Gate Components (theta, match, topo, slip, tail, dir, pollution, unified): "
                      << result.avg_gate_w_theta << ", "
                      << result.avg_gate_w_match << ", "
                      << result.avg_gate_w_topo << ", "
                      << result.avg_gate_w_slip << ", "
                      << result.avg_gate_w_tail << ", "
                      << result.avg_gate_w_dir << ", "
                      << result.avg_gate_w_pollution << ", "
                      << result.avg_gate_w_unified << std::endl;
            std::cout << "    Avg Direction Consistency: " << result.avg_direction_consistency << std::endl;
            std::cout << "    Avg Basis Distortion: " << result.avg_basis_distortion << std::endl;
            std::cout << "    Avg Direction Pollution: " << result.avg_direction_pollution << std::endl;
            std::cout << "    Birth Events: " << result.total_birth_events << std::endl;
            std::cout << "    Stable: " << (result.is_stable ? "YES" : "NO") << std::endl;
        }

        all_results.push_back(case_results);
    }

    // ========================================================================
    // Write output files
    // ========================================================================

    // Case A output
    std::ofstream case_A_file(out_dir + "/sdf_patch_adaptive_carryover_case_A.csv");
    case_A_file << "config,clamp,attenuation,carry_strategy,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_raw_transport_norm,avg_limited_transport_norm,avg_transport_attenuation,avg_carry_over_factor,min_carry_over_factor,max_carry_over_factor,avg_adaptive_alpha,min_adaptive_alpha,max_adaptive_alpha,avg_theta,max_theta,avg_match_distance,theta_p50,theta_p95,match_distance_p50,match_distance_p95,avg_theta_rel,theta_rel_p50,theta_rel_p95,theta_rel_p99,max_theta_rel,avg_match_rel,match_rel_p50,match_rel_p95,match_rel_p99,max_match_rel,avg_risk_combo,risk_combo_p50,risk_combo_p95,risk_combo_p99,max_risk_combo,avg_gate_w_theta,avg_gate_w_match,avg_gate_w_topo,avg_gate_w_slip,avg_gate_w_tail,avg_gate_w_dir,avg_direction_consistency,direction_consistency_p50,direction_consistency_p95,direction_consistency_p99,avg_basis_distortion,basis_distortion_p50,basis_distortion_p95,basis_distortion_p99,avg_gate_w_unified,basis_gate_before_avg,basis_gate_after_avg,avg_direction_pollution,direction_pollution_p50,direction_pollution_p95,direction_pollution_p99,avg_gate_w_pollution,pollution_gate_before_avg,pollution_gate_after_avg,direction_bad_gate_before_avg,direction_bad_gate_after_avg,high_risk_gate_avg,high_risk_v21_gate_avg,median_risk_gate_avg,median_risk_v21_gate_avg,total_birth_events,total_transport_clamp_hits,clamp_trigger_count,clamp_pre_norm_avg,clamp_post_norm_avg,clamp_max_pre_norm,clamp_max_post_norm,stick_steps,slip_steps,stick_to_slip_transitions,stable" << std::endl;
    for (const auto& r : all_results[0]) {
        std::string carry_str = CarryOverStrategyName(r.carry_strategy);
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
                    << r.theta_p50 << ","
                    << r.theta_p95 << ","
                    << r.match_distance_p50 << ","
                    << r.match_distance_p95 << ","
                    << r.avg_theta_rel << ","
                    << r.theta_rel_p50 << ","
                    << r.theta_rel_p95 << ","
                    << r.theta_rel_p99 << ","
                    << r.max_theta_rel << ","
                    << r.avg_match_rel << ","
                    << r.match_rel_p50 << ","
                    << r.match_rel_p95 << ","
                    << r.match_rel_p99 << ","
                    << r.max_match_rel << ","
                    << r.avg_risk_combo << ","
                    << r.risk_combo_p50 << ","
                    << r.risk_combo_p95 << ","
                    << r.risk_combo_p99 << ","
                    << r.max_risk_combo << ","
                    << r.avg_gate_w_theta << ","
                    << r.avg_gate_w_match << ","
                    << r.avg_gate_w_topo << ","
                    << r.avg_gate_w_slip << ","
                    << r.avg_gate_w_tail << ","
                    << r.avg_gate_w_dir << ","
                    << r.avg_direction_consistency << ","
                    << r.direction_consistency_p50 << ","
                    << r.direction_consistency_p95 << ","
                    << r.direction_consistency_p99 << ","
                    << r.avg_basis_distortion << ","
                    << r.basis_distortion_p50 << ","
                    << r.basis_distortion_p95 << ","
                    << r.basis_distortion_p99 << ","
                    << r.avg_gate_w_unified << ","
                    << r.basis_gate_before_avg << ","
                    << r.basis_gate_after_avg << ","
                    << r.avg_direction_pollution << ","
                    << r.direction_pollution_p50 << ","
                    << r.direction_pollution_p95 << ","
                    << r.direction_pollution_p99 << ","
                    << r.avg_gate_w_pollution << ","
                    << r.pollution_gate_before_avg << ","
                    << r.pollution_gate_after_avg << ","
                    << r.direction_bad_gate_before_avg << ","
                    << r.direction_bad_gate_after_avg << ","
                    << r.high_risk_gate_avg << ","
                    << r.high_risk_v21_gate_avg << ","
                    << r.median_risk_gate_avg << ","
                    << r.median_risk_v21_gate_avg << ","
                    << r.total_birth_events << ","
                    << r.total_transport_clamp_hits << ","
                    << r.clamp_trigger_count << ","
                    << r.clamp_pre_norm_avg << ","
                    << r.clamp_post_norm_avg << ","
                    << r.clamp_max_pre_norm << ","
                    << r.clamp_max_post_norm << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_A_file.close();

    // Case B output
    std::ofstream case_B_file(out_dir + "/sdf_patch_adaptive_carryover_case_B.csv");
    case_B_file << "config,clamp,attenuation,carry_strategy,final_y,expected_y,y_error,avg_force_y,force_std_dev,avg_torque_x,avg_torque_y,avg_torque_z,torque_std_dev,max_patch_count,avg_patch_count,multi_patch_ratio,avg_tangential_force_norm,max_tangential_force_norm,avg_tangential_force_ratio,avg_tangential_displacement_norm,avg_raw_transport_norm,avg_limited_transport_norm,avg_transport_attenuation,avg_carry_over_factor,min_carry_over_factor,max_carry_over_factor,avg_adaptive_alpha,min_adaptive_alpha,max_adaptive_alpha,avg_theta,max_theta,avg_match_distance,theta_p50,theta_p95,match_distance_p50,match_distance_p95,avg_theta_rel,theta_rel_p50,theta_rel_p95,theta_rel_p99,max_theta_rel,avg_match_rel,match_rel_p50,match_rel_p95,match_rel_p99,max_match_rel,avg_risk_combo,risk_combo_p50,risk_combo_p95,risk_combo_p99,max_risk_combo,avg_gate_w_theta,avg_gate_w_match,avg_gate_w_topo,avg_gate_w_slip,avg_gate_w_tail,avg_gate_w_dir,avg_direction_consistency,direction_consistency_p50,direction_consistency_p95,direction_consistency_p99,avg_basis_distortion,basis_distortion_p50,basis_distortion_p95,basis_distortion_p99,avg_gate_w_unified,basis_gate_before_avg,basis_gate_after_avg,avg_direction_pollution,direction_pollution_p50,direction_pollution_p95,direction_pollution_p99,avg_gate_w_pollution,pollution_gate_before_avg,pollution_gate_after_avg,direction_bad_gate_before_avg,direction_bad_gate_after_avg,high_risk_gate_avg,high_risk_v21_gate_avg,median_risk_gate_avg,median_risk_v21_gate_avg,total_birth_events,total_transport_clamp_hits,clamp_trigger_count,clamp_pre_norm_avg,clamp_post_norm_avg,clamp_max_pre_norm,clamp_max_post_norm,stick_steps,slip_steps,stick_to_slip_transitions,stable" << std::endl;
    for (const auto& r : all_results[1]) {
        std::string carry_str = CarryOverStrategyName(r.carry_strategy);
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
                    << r.theta_p50 << ","
                    << r.theta_p95 << ","
                    << r.match_distance_p50 << ","
                    << r.match_distance_p95 << ","
                    << r.avg_theta_rel << ","
                    << r.theta_rel_p50 << ","
                    << r.theta_rel_p95 << ","
                    << r.theta_rel_p99 << ","
                    << r.max_theta_rel << ","
                    << r.avg_match_rel << ","
                    << r.match_rel_p50 << ","
                    << r.match_rel_p95 << ","
                    << r.match_rel_p99 << ","
                    << r.max_match_rel << ","
                    << r.avg_risk_combo << ","
                    << r.risk_combo_p50 << ","
                    << r.risk_combo_p95 << ","
                    << r.risk_combo_p99 << ","
                    << r.max_risk_combo << ","
                    << r.avg_gate_w_theta << ","
                    << r.avg_gate_w_match << ","
                    << r.avg_gate_w_topo << ","
                    << r.avg_gate_w_slip << ","
                    << r.avg_gate_w_tail << ","
                    << r.avg_gate_w_dir << ","
                    << r.avg_direction_consistency << ","
                    << r.direction_consistency_p50 << ","
                    << r.direction_consistency_p95 << ","
                    << r.direction_consistency_p99 << ","
                    << r.avg_basis_distortion << ","
                    << r.basis_distortion_p50 << ","
                    << r.basis_distortion_p95 << ","
                    << r.basis_distortion_p99 << ","
                    << r.avg_gate_w_unified << ","
                    << r.basis_gate_before_avg << ","
                    << r.basis_gate_after_avg << ","
                    << r.avg_direction_pollution << ","
                    << r.direction_pollution_p50 << ","
                    << r.direction_pollution_p95 << ","
                    << r.direction_pollution_p99 << ","
                    << r.avg_gate_w_pollution << ","
                    << r.pollution_gate_before_avg << ","
                    << r.pollution_gate_after_avg << ","
                    << r.direction_bad_gate_before_avg << ","
                    << r.direction_bad_gate_after_avg << ","
                    << r.high_risk_gate_avg << ","
                    << r.high_risk_v21_gate_avg << ","
                    << r.median_risk_gate_avg << ","
                    << r.median_risk_v21_gate_avg << ","
                    << r.total_birth_events << ","
                    << r.total_transport_clamp_hits << ","
                    << r.clamp_trigger_count << ","
                    << r.clamp_pre_norm_avg << ","
                    << r.clamp_post_norm_avg << ","
                    << r.clamp_max_pre_norm << ","
                    << r.clamp_max_post_norm << ","
                    << r.total_stick_steps << ","
                    << r.total_slip_steps << ","
                    << r.total_stick_to_slip_transitions << ","
                    << (r.is_stable ? 1 : 0) << std::endl;
    }
    case_B_file.close();

    // Summary
    std::ofstream summary_file(out_dir + "/sdf_patch_adaptive_carryover_summary.csv");
    summary_file << "config,clamp,attenuation,carry_strategy,case_A_y_error,case_A_torque_z,case_A_trans_norm,case_A_alpha_mean,case_A_alpha_min,case_A_stable,case_B_y_error,case_B_torque_z,case_B_trans_norm,case_B_alpha_mean,case_B_alpha_min,case_B_stable" << std::endl;
    for (size_t ci = 0; ci < configs.size(); ci++) {
        std::string carry_str = CarryOverStrategyName(configs[ci].carry_strategy);
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
    std::ofstream diag_file(out_dir + "/sdf_patch_adaptive_carryover_diagnostics.csv");
    diag_file << "config,case_A_avg_theta,case_A_theta_p50,case_A_theta_p95,case_A_avg_theta_rel,case_A_theta_rel_p50,case_A_theta_rel_p95,case_A_max_theta,case_A_avg_match_dist,case_A_match_p50,case_A_match_p95,case_A_avg_match_rel,case_A_match_rel_p50,case_A_match_rel_p95,case_A_w_theta,case_A_w_match,case_A_w_topo,case_A_w_slip,case_A_birth_events,case_A_avg_alpha,case_A_min_alpha,case_A_max_alpha,case_B_avg_theta,case_B_theta_p50,case_B_theta_p95,case_B_avg_theta_rel,case_B_theta_rel_p50,case_B_theta_rel_p95,case_B_max_theta,case_B_avg_match_dist,case_B_match_p50,case_B_match_p95,case_B_avg_match_rel,case_B_match_rel_p50,case_B_match_rel_p95,case_B_w_theta,case_B_w_match,case_B_w_topo,case_B_w_slip,case_B_birth_events,case_B_avg_alpha,case_B_min_alpha,case_B_max_alpha" << std::endl;
    for (size_t ci = 0; ci < configs.size(); ci++) {
        diag_file << configs[ci].name << ","
                  << std::fixed << std::setprecision(6)
                  << all_results[0][ci].avg_theta << "," << all_results[0][ci].theta_p50 << "," << all_results[0][ci].theta_p95 << ","
                  << all_results[0][ci].avg_theta_rel << "," << all_results[0][ci].theta_rel_p50 << "," << all_results[0][ci].theta_rel_p95 << "," << all_results[0][ci].max_theta << ","
                  << all_results[0][ci].avg_match_distance << "," << all_results[0][ci].match_distance_p50 << "," << all_results[0][ci].match_distance_p95 << ","
                  << all_results[0][ci].avg_match_rel << "," << all_results[0][ci].match_rel_p50 << "," << all_results[0][ci].match_rel_p95 << ","
                  << all_results[0][ci].avg_gate_w_theta << "," << all_results[0][ci].avg_gate_w_match << "," << all_results[0][ci].avg_gate_w_topo << "," << all_results[0][ci].avg_gate_w_slip << ","
                  << all_results[0][ci].total_birth_events << "," << all_results[0][ci].avg_adaptive_alpha << "," << all_results[0][ci].min_adaptive_alpha << "," << all_results[0][ci].max_adaptive_alpha << ","
                  << all_results[1][ci].avg_theta << "," << all_results[1][ci].theta_p50 << "," << all_results[1][ci].theta_p95 << ","
                  << all_results[1][ci].avg_theta_rel << "," << all_results[1][ci].theta_rel_p50 << "," << all_results[1][ci].theta_rel_p95 << "," << all_results[1][ci].max_theta << ","
                  << all_results[1][ci].avg_match_distance << "," << all_results[1][ci].match_distance_p50 << "," << all_results[1][ci].match_distance_p95 << ","
                  << all_results[1][ci].avg_match_rel << "," << all_results[1][ci].match_rel_p50 << "," << all_results[1][ci].match_rel_p95 << ","
                  << all_results[1][ci].avg_gate_w_theta << "," << all_results[1][ci].avg_gate_w_match << "," << all_results[1][ci].avg_gate_w_topo << "," << all_results[1][ci].avg_gate_w_slip << ","
                  << all_results[1][ci].total_birth_events << "," << all_results[1][ci].avg_adaptive_alpha << "," << all_results[1][ci].min_adaptive_alpha << "," << all_results[1][ci].max_adaptive_alpha << std::endl;
    }
    diag_file.close();

    // ========================================================================
    // Console summary
    // ========================================================================

    std::cout << "\n=== Case A Adaptive Carry-Over Summary ===" << std::endl;
    std::cout << std::left << std::setw(18) << "Config"
              << std::setw(12) << "y_error"
              << std::setw(12) << "torque_z"
              << std::setw(12) << "trans_norm"
              << std::setw(12) << "alpha_avg"
              << std::setw(12) << "alpha_min"
              << std::setw(8) << "stable" << std::endl;
    std::cout << std::string(74, '-') << std::endl;
    for (const auto& r : all_results[0]) {
        std::cout << std::left << std::setw(18) << r.config_name
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << r.y_error
                  << std::setw(12) << r.avg_torque_z
                  << std::setw(12) << r.avg_limited_transport_norm
                  << std::setw(12) << r.avg_adaptive_alpha
                  << std::setw(12) << r.min_adaptive_alpha
                  << std::setw(8) << (r.is_stable ? "YES" : "NO") << std::endl;
    }

    std::cout << "\n=== Case B Adaptive Carry-Over Summary ===" << std::endl;
    std::cout << std::left << std::setw(18) << "Config"
              << std::setw(12) << "y_error"
              << std::setw(12) << "torque_z"
              << std::setw(12) << "trans_norm"
              << std::setw(12) << "alpha_avg"
              << std::setw(12) << "alpha_min"
              << std::setw(8) << "stable" << std::endl;
    std::cout << std::string(74, '-') << std::endl;
    for (const auto& r : all_results[1]) {
        std::cout << std::left << std::setw(18) << r.config_name
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

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_adaptive_carryover_case_A.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_adaptive_carryover_case_B.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_adaptive_carryover_summary.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_adaptive_carryover_diagnostics.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_transport_theory.csv" << std::endl;
    std::cout << "  " << out_dir + "/sdf_patch_tangential_theory_summary.csv" << std::endl;

    return 0;
}
