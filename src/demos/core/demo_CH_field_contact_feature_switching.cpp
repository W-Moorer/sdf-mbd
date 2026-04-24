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
// Milestone 20: feature-switching-sensitive field-contact scenes.
//
// This demo keeps the SDF backend analytic so the comparison isolates contact
// primitive history, transport, and split/merge logic. A kinematic sphere slides
// across non-smooth/non-convex height fields and each frame is evaluated with:
//   - current runtime field primitive;
//   - no tangential history;
//   - direct nearest-source history inheritance;
//   - projection-only transport;
//   - minimal-rotation transport with gate and split/merge aggregation.
//
// Outputs:
//   out/milestone_20/field_contact_feature_switching_frames.csv
//   out/milestone_20/field_contact_feature_switching_patches.csv
//   out/milestone_20/field_contact_feature_switching_summary.csv
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

enum class FeatureScenario {
    FoldedSeam,
    ConcavePolylineCorner,
    NarrowGrooveEntrance
};

enum class TransportVariant {
    RuntimeFieldPrimitive,
    NoHistory,
    DirectInherit,
    ProjectionTransport,
    MinimalRotationGateSplitMerge
};

struct ScenarioConfig {
    FeatureScenario scenario;
    std::string name;
    double sphere_radius = 0.16;
    double bottom_y = 0.0;
    double x0 = -0.32;
    double x1 = 0.32;
    double z0 = 0.0;
    double z1 = 0.0;
    double total_time = 1.20;
    int frames = 601;
};

struct VariantState {
    TransportVariant variant;
    FieldContactPrimitiveTracker runtime_tracker;
    std::vector<PrimitiveSnapshot> previous_snapshots;
    std::map<int, TangentialHistory> history_store;
    FieldContactTopologyMetricsAccumulator metrics;
    int next_persistent_id = 0;
};

struct ScenarioVariantSummary {
    std::string scenario;
    std::string variant;
    FieldContactTopologyMetricsSummary metrics;
};

static const char* ScenarioName(FeatureScenario scenario) {
    switch (scenario) {
        case FeatureScenario::FoldedSeam:
            return "folded_seam";
        case FeatureScenario::ConcavePolylineCorner:
            return "concave_polyline_corner";
        case FeatureScenario::NarrowGrooveEntrance:
            return "narrow_groove_entrance";
    }
    return "unknown";
}

static const char* VariantName(TransportVariant variant) {
    switch (variant) {
        case TransportVariant::RuntimeFieldPrimitive:
            return "current_field_primitive";
        case TransportVariant::NoHistory:
            return "no_history";
        case TransportVariant::DirectInherit:
            return "direct_inherit";
        case TransportVariant::ProjectionTransport:
            return "projection_transport";
        case TransportVariant::MinimalRotationGateSplitMerge:
            return "minimal_rotation_gate_split_merge";
    }
    return "unknown";
}

static std::string GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "paper")) {
            return path.string();
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path().string();
}

static double SmoothStep(double edge0, double edge1, double x) {
    if (edge0 == edge1) {
        return x >= edge1 ? 1.0 : 0.0;
    }
    double t = Clamp01((x - edge0) / (edge1 - edge0));
    return t * t * (3.0 - 2.0 * t);
}

static double SmoothAbs(double x, double eps) {
    return std::sqrt(x * x + eps * eps);
}

class AnalyticFeatureHeightField {
  public:
    explicit AnalyticFeatureHeightField(FeatureScenario scenario) : m_scenario(scenario) {}

    FieldSampleQuery Query(const ChVector3d& world_pos, const ChVector3d& world_vel) const {
        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = world_vel;

        query.phi = world_pos.y() - Height(world_pos.x(), world_pos.z());

        const double h = 2.0e-4;
        const double dhdx = (Height(world_pos.x() + h, world_pos.z()) -
                             Height(world_pos.x() - h, world_pos.z())) /
                            (2.0 * h);
        const double dhdz = (Height(world_pos.x(), world_pos.z() + h) -
                             Height(world_pos.x(), world_pos.z() - h)) /
                            (2.0 * h);
        query.grad = SafeNormalize(ChVector3d(-dhdx, 1.0, -dhdz), ChVector3d(0, 1, 0));
        return query;
    }

  private:
    double Height(double x, double z) const {
        switch (m_scenario) {
            case FeatureScenario::FoldedSeam:
                return FoldedSeamHeight(x, z);
            case FeatureScenario::ConcavePolylineCorner:
                return ConcaveCornerHeight(x, z);
            case FeatureScenario::NarrowGrooveEntrance:
                return NarrowGrooveHeight(x, z);
        }
        return 0.0;
    }

    static double FoldedSeamHeight(double x, double z) {
        const double slope = 0.30;
        const double seam_softening = 0.0010;
        const double crown = 0.010;
        const double lateral_camber = -0.006 * std::exp(-(z * z) / (2.0 * 0.24 * 0.24));
        return crown + slope * SmoothAbs(x, seam_softening) + lateral_camber;
    }

    static double ConcaveCornerHeight(double x, double z) {
        const double half_width = 0.145;
        const double depth = 0.040;
        const double floor = 0.010;
        const double side_slope = 0.22;
        double triangular_notch = std::max(0.0, 1.0 - SmoothAbs(x, 0.0008) / half_width);
        double local_relief = -depth * triangular_notch;
        double outside_ramp = side_slope * std::max(0.0, SmoothAbs(x, 0.0008) - half_width);
        double z_taper = 1.0 - 0.12 * SmoothStep(0.08, 0.24, std::abs(z));
        return (floor + local_relief + outside_ramp) * z_taper;
    }

    static double NarrowGrooveHeight(double x, double z) {
        double constriction = SmoothStep(-0.30, -0.08, x) * (1.0 - SmoothStep(0.08, 0.30, x));
        double half_width = 0.105 + (0.032 - 0.105) * constriction;
        double depth = 0.038;
        double shoulder = 0.006 * constriction;
        double normalized_z = z / std::max(half_width, 1.0e-6);
        double groove_profile = std::exp(-std::pow(normalized_z, 4.0));
        return shoulder - depth * groove_profile;
    }

    FeatureScenario m_scenario;
};

static FieldContactRuntimeSettings MakeRuntimeSettings(double time_step) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = 0.004;
    settings.extraction.min_area = 1.0e-8;
    settings.extraction.min_samples = 3;
    settings.extraction.use_penetration_weighted_center = true;

    settings.normal.stiffness = 1.6e6;
    settings.normal.damping = 0.0;

    settings.tangential.stiffness = 2.5e3;
    settings.tangential.damping = 0.0;
    settings.tangential.friction_coefficient = 0.45;
    settings.tangential.time_step = time_step;

    settings.inheritance.min_overlap = 0.006;
    settings.inheritance.min_normal_dot = 0.20;
    settings.inheritance.max_center_distance = 0.075;
    settings.inheritance.geometry_fallback_weight = 0.20;
    return settings;
}

static ChVector3d Interpolate(const ChVector3d& a, const ChVector3d& b, double u) {
    return a + (b - a) * u;
}

static std::vector<FieldSampleQuery> BuildSphereQueries(const SurfaceGraph& graph,
                                                        const AnalyticFeatureHeightField& field,
                                                        const ChVector3d& center,
                                                        const ChVector3d& velocity) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(graph.samples.size());
    for (const auto& sample : graph.samples) {
        ChVector3d world_pos = center + sample.local_pos;
        queries.push_back(field.Query(world_pos, velocity));
    }
    return queries;
}

static void ClassifyPatches(std::map<int, std::vector<int>>& primary_users,
                            FieldContactStepResult& result,
                            int& next_persistent_id) {
    for (size_t pi = 0; pi < result.patches.size(); pi++) {
        auto& patch_result = result.patches[pi];
        if (patch_result.sources.empty()) {
            patch_result.persistent_id = next_persistent_id++;
            patch_result.event = FieldContactPrimitiveEvent::Newborn;
            result.stats.newborn_count++;
            continue;
        }

        int primary = patch_result.sources.front().persistent_id;
        auto users_it = primary_users.find(primary);
        if (users_it == primary_users.end()) {
            patch_result.persistent_id = next_persistent_id++;
            patch_result.event = FieldContactPrimitiveEvent::Newborn;
            result.stats.newborn_count++;
            continue;
        }

        const auto& users = users_it->second;
        bool split_child = users.size() > 1 && users.front() != static_cast<int>(pi);
        bool merge_patch = patch_result.sources.size() > 1;

        if (split_child) {
            patch_result.persistent_id = next_persistent_id++;
            patch_result.event = merge_patch ? FieldContactPrimitiveEvent::SplitMerge :
                                               FieldContactPrimitiveEvent::Split;
            result.stats.split_count++;
        } else {
            patch_result.persistent_id = primary;
            if (merge_patch) {
                patch_result.event = FieldContactPrimitiveEvent::Merge;
                result.stats.merge_count++;
            } else if (users.size() > 1) {
                patch_result.event = FieldContactPrimitiveEvent::SplitPrimary;
                result.stats.split_count++;
            } else {
                patch_result.event = FieldContactPrimitiveEvent::Stable;
            }
        }
    }
}

static TangentialUpdateResult UpdateTangentialContactFromPreparedXi(const ChVector3d& prepared_xi,
                                                                    bool has_prepared_history,
                                                                    const ChVector3d& current_normal,
                                                                    const ChVector3d& current_tangential_velocity,
                                                                    double normal_force_magnitude,
                                                                    double history_gate,
                                                                    const TangentialContactSettings& settings) {
    TangentialUpdateResult result;
    ChVector3d normal = SafeNormalize(current_normal, ChVector3d(0, 1, 0));
    ChVector3d vt = ProjectToTangent(current_tangential_velocity, normal);
    result.gate = Clamp01(history_gate);
    result.friction_limit = std::max(0.0, settings.friction_coefficient * normal_force_magnitude);

    if (has_prepared_history) {
        result.transported_xi = prepared_xi;
    }

    result.stored_energy_before_gate = 0.5 * settings.stiffness * result.transported_xi.Dot(result.transported_xi);
    result.gated_xi = result.transported_xi * result.gate;
    result.stored_energy_after_gate = 0.5 * settings.stiffness * result.gated_xi.Dot(result.gated_xi);

    result.trial_xi = ProjectToTangent(result.gated_xi + vt * settings.time_step, normal);
    ChVector3d elastic_force_trial = -settings.stiffness * result.trial_xi;
    result.elastic_trial_force_norm = elastic_force_trial.Length();

    ChVector3d xi_elastic = result.trial_xi;
    ChVector3d elastic_force = elastic_force_trial;

    if (settings.stiffness <= 0.0 || result.friction_limit <= 1.0e-14) {
        xi_elastic = ChVector3d(0, 0, 0);
        elastic_force = ChVector3d(0, 0, 0);
        result.state = StickSlipState::Slip;
    } else if (result.elastic_trial_force_norm > result.friction_limit) {
        double scale = result.friction_limit / result.elastic_trial_force_norm;
        elastic_force = elastic_force_trial * scale;
        xi_elastic = -elastic_force / settings.stiffness;
        result.state = StickSlipState::Slip;
    } else {
        result.state = StickSlipState::Stick;
    }

    ChVector3d damping_force = -settings.damping * vt;
    ChVector3d total_force = elastic_force + damping_force;
    double total_norm = total_force.Length();
    if (total_norm > result.friction_limit && result.friction_limit > 1.0e-14) {
        total_force *= result.friction_limit / total_norm;
        result.state = StickSlipState::Slip;
    }

    result.force = ProjectToTangent(total_force, normal);
    result.final_force_norm = result.force.Length();
    result.history.valid = true;
    result.history.xi_elastic_world = ProjectToTangent(xi_elastic, normal);
    result.history.normal = normal;
    return result;
}

static void ApplyTangentialVariant(FieldContactPatchRuntimeResult& patch_result,
                                   TransportVariant variant,
                                   const ChVector3d& torque_reference,
                                   const FieldContactRuntimeSettings& settings,
                                   const std::map<int, TangentialHistory>& old_history_store,
                                   std::map<int, TangentialHistory>& new_history_store) {
    ChVector3d prepared_xi(0, 0, 0);
    bool has_prepared_history = false;
    double history_gate = 0.0;

    if (variant == TransportVariant::MinimalRotationGateSplitMerge) {
        std::vector<WeightedTangentialHistorySource> weighted_history_sources;
        for (const auto& source : patch_result.sources) {
            auto hist_it = old_history_store.find(source.persistent_id);
            if (hist_it == old_history_store.end() || !hist_it->second.valid) {
                continue;
            }

            WeightedTangentialHistorySource weighted_source;
            weighted_source.history = hist_it->second;
            weighted_source.weight = source.weight;
            weighted_history_sources.push_back(weighted_source);

            ChVector3d source_xi = TransportElasticStateMinimalRotation(hist_it->second.xi_elastic_world,
                                                                        hist_it->second.normal,
                                                                        patch_result.patch.normal);
            patch_result.source_energy_bound += source.weight * 0.5 * settings.tangential.stiffness *
                                                source_xi.Dot(source_xi);
        }

        TangentialHistory aggregated_history =
            AggregateTangentialHistorySources(weighted_history_sources,
                                              patch_result.patch.normal,
                                              patch_result.persistent_id);
        if (aggregated_history.valid) {
            prepared_xi = aggregated_history.xi_elastic_world;
            has_prepared_history = true;
            history_gate = 1.0;
            patch_result.inherited_energy = 0.5 * settings.tangential.stiffness * prepared_xi.Dot(prepared_xi);
        }
    } else if (variant != TransportVariant::NoHistory && !patch_result.sources.empty()) {
        for (const auto& source : patch_result.sources) {
            auto hist_it = old_history_store.find(source.persistent_id);
            if (hist_it == old_history_store.end() || !hist_it->second.valid) {
                continue;
            }

            has_prepared_history = true;
            const TangentialHistory& old_history = hist_it->second;

            if (variant == TransportVariant::DirectInherit) {
                prepared_xi += old_history.xi_elastic_world;
                patch_result.source_energy_bound += 0.5 * settings.tangential.stiffness *
                                                    old_history.xi_elastic_world.Dot(old_history.xi_elastic_world);
            } else if (variant == TransportVariant::ProjectionTransport) {
                ChVector3d projected_xi = ProjectToTangent(old_history.xi_elastic_world, patch_result.patch.normal);
                prepared_xi += projected_xi * source.weight;
                patch_result.source_energy_bound += source.weight * 0.5 * settings.tangential.stiffness *
                                                    projected_xi.Dot(projected_xi);
            }
        }

        if (has_prepared_history) {
            history_gate = 1.0;
            patch_result.inherited_energy = 0.5 * settings.tangential.stiffness * prepared_xi.Dot(prepared_xi);
        }
    }

    if (patch_result.source_energy_bound > 1.0e-16) {
        patch_result.inherited_energy_ratio = patch_result.inherited_energy / patch_result.source_energy_bound;
    }

    ChVector3d vt = ProjectToTangent(patch_result.patch.representative_velocity, patch_result.patch.normal);
    patch_result.tangential =
        UpdateTangentialContactFromPreparedXi(prepared_xi,
                                              has_prepared_history,
                                              patch_result.patch.normal,
                                              vt,
                                              patch_result.patch.normal_force.Length(),
                                              history_gate,
                                              settings.tangential);
    patch_result.tangential.history.persistent_id = patch_result.persistent_id;
    new_history_store[patch_result.persistent_id] = patch_result.tangential.history;

    patch_result.patch.tangential_force = patch_result.tangential.force;
    patch_result.patch.force = patch_result.patch.normal_force + patch_result.patch.tangential_force;
    patch_result.patch.torque += (patch_result.patch.center - torque_reference).Cross(patch_result.patch.tangential_force);

    patch_result.tangential_force_ratio = patch_result.tangential.friction_limit > 1.0e-12 ?
                                              patch_result.tangential.final_force_norm /
                                                  patch_result.tangential.friction_limit :
                                              0.0;
    patch_result.energy_gate_ratio = patch_result.tangential.stored_energy_before_gate > 1.0e-16 ?
                                         patch_result.tangential.stored_energy_after_gate /
                                             patch_result.tangential.stored_energy_before_gate :
                                         0.0;
}

static FieldContactStepResult EvaluateVariant(VariantState& state,
                                              const SurfaceGraph& graph,
                                              const std::vector<FieldSampleQuery>& queries,
                                              const ChVector3d& torque_reference,
                                              const FieldContactRuntimeSettings& settings) {
    if (state.variant == TransportVariant::RuntimeFieldPrimitive) {
        return state.runtime_tracker.Evaluate(graph, queries, torque_reference, settings);
    }

    std::vector<int> active = BuildActiveSet(queries, settings.extraction.activation_band);
    std::vector<PrimitivePatch> primitive_patches =
        ExtractPrimitives(graph, queries, active, settings.extraction);

    FieldContactStepResult result;
    result.stats.patch_count = static_cast<int>(primitive_patches.size());
    result.patches.resize(primitive_patches.size());

    std::map<int, std::vector<int>> primary_users;
    std::set<int> referenced_previous;

    for (size_t pi = 0; pi < primitive_patches.size(); pi++) {
        auto& patch_result = result.patches[pi];
        patch_result.patch = primitive_patches[pi];
        patch_result.sources =
            ComputeHistorySources(patch_result.patch, state.previous_snapshots, graph, settings.inheritance);

        result.stats.max_source_count =
            std::max(result.stats.max_source_count, static_cast<int>(patch_result.sources.size()));

        for (const auto& source : patch_result.sources) {
            patch_result.source_weight_sum += source.weight;
            referenced_previous.insert(source.persistent_id);
        }

        if (!patch_result.sources.empty()) {
            primary_users[patch_result.sources.front().persistent_id].push_back(static_cast<int>(pi));
        }
    }

    for (auto& entry : primary_users) {
        auto& users = entry.second;
        std::sort(users.begin(), users.end(), [&](int a, int b) {
            double wa = result.patches[a].sources.empty() ? 0.0 : result.patches[a].sources.front().weight;
            double wb = result.patches[b].sources.empty() ? 0.0 : result.patches[b].sources.front().weight;
            return wa > wb;
        });
        result.stats.max_previous_reuse =
            std::max(result.stats.max_previous_reuse, static_cast<int>(users.size()));
    }

    ClassifyPatches(primary_users, result, state.next_persistent_id);

    for (const auto& prev : state.previous_snapshots) {
        if (!referenced_previous.count(prev.persistent_id)) {
            result.stats.death_count++;
        }
    }

    std::vector<PrimitiveSnapshot> current_snapshots;
    std::map<int, TangentialHistory> new_history_store;

    for (auto& patch_result : result.patches) {
        ApplyNormalContactIntegral(patch_result.patch, graph, queries, torque_reference, settings.normal);
        ApplyTangentialVariant(patch_result,
                               state.variant,
                               torque_reference,
                               settings,
                               state.history_store,
                               new_history_store);

        result.total_force += patch_result.patch.force;
        result.total_torque += patch_result.patch.torque;
        current_snapshots.push_back(MakeSnapshot(patch_result.patch, patch_result.persistent_id));

        result.stats.max_tangential_force_ratio =
            std::max(result.stats.max_tangential_force_ratio, patch_result.tangential_force_ratio);
        result.stats.max_energy_gate_ratio =
            std::max(result.stats.max_energy_gate_ratio, patch_result.energy_gate_ratio);
        result.stats.max_inherited_energy_ratio =
            std::max(result.stats.max_inherited_energy_ratio, patch_result.inherited_energy_ratio);
    }

    state.previous_snapshots = current_snapshots;
    state.history_store = new_history_store;
    return result;
}

static std::string JoinSources(const std::vector<HistorySource>& sources) {
    std::ostringstream ss;
    for (size_t i = 0; i < sources.size(); i++) {
        if (i > 0) {
            ss << ";";
        }
        ss << sources[i].persistent_id << ":" << std::fixed << std::setprecision(4) << sources[i].weight;
    }
    return ss.str();
}

static void WriteFrameRow(std::ofstream& out,
                          const std::string& scenario,
                          const std::string& variant,
                          int frame,
                          double time,
                          const ChVector3d& center,
                          const ChVector3d& velocity,
                          const FieldContactStepResult& step) {
    out << scenario << "," << variant << "," << frame << "," << time << "," << center.x() << "," << center.y()
        << "," << center.z() << "," << velocity.x() << "," << velocity.y() << "," << velocity.z() << ","
        << step.stats.patch_count << "," << step.stats.newborn_count << "," << step.stats.merge_count << ","
        << step.stats.split_count << "," << step.stats.death_count << "," << step.stats.max_source_count << ","
        << step.stats.max_previous_reuse << "," << step.total_force.x() << "," << step.total_force.y() << ","
        << step.total_force.z() << "," << step.total_force.Length() << "," << step.total_torque.x() << ","
        << step.total_torque.y() << "," << step.total_torque.z() << "," << step.total_torque.Length() << ","
        << step.stats.max_tangential_force_ratio << "," << step.stats.max_inherited_energy_ratio << "\n";
}

static void WritePatchRows(std::ofstream& out,
                           const std::string& scenario,
                           const std::string& variant,
                           int frame,
                           double time,
                           const FieldContactStepResult& step) {
    for (size_t i = 0; i < step.patches.size(); i++) {
        const auto& patch = step.patches[i];
        out << scenario << "," << variant << "," << frame << "," << time << "," << i << ","
            << patch.persistent_id << "," << FieldContactPrimitiveEventName(patch.event) << ","
            << patch.sources.size() << "," << patch.source_weight_sum << "," << patch.patch.sample_ids.size()
            << "," << patch.patch.area << "," << patch.patch.mean_phi << "," << patch.patch.max_penetration
            << "," << patch.patch.center.x() << "," << patch.patch.center.y() << "," << patch.patch.center.z()
            << "," << patch.patch.normal.x() << "," << patch.patch.normal.y() << "," << patch.patch.normal.z()
            << "," << patch.patch.normal_force.Length() << "," << patch.patch.tangential_force.Length() << ","
            << patch.tangential_force_ratio << "," << patch.energy_gate_ratio << ","
            << patch.inherited_energy_ratio << "," << (patch.tangential.state == StickSlipState::Stick ? "stick" : "slip")
            << ",\"" << JoinSources(patch.sources) << "\"\n";
    }
}

static void WriteSummaryRow(std::ofstream& out, const ScenarioVariantSummary& summary) {
    const auto& m = summary.metrics;
    out << summary.scenario << "," << summary.variant << "," << m.frames << "," << m.active_frames << ","
        << m.total_patch_count << "," << m.mean_patch_count << "," << m.max_patch_count << ","
        << m.total_newborn << "," << m.total_death << "," << m.total_merge_patches << ","
        << m.total_split_patches << "," << m.max_source_count << "," << m.max_previous_reuse << ","
        << m.primitive_persistence_ratio << "," << m.mean_source_weight_sum << ","
        << m.mean_abs_patch_count_change << "," << m.max_abs_patch_count_change << ","
        << m.mean_normal_jump_angle << "," << m.rms_normal_jump_angle << "," << m.max_normal_jump_angle << ","
        << m.mean_force_jump << "," << m.rms_force_jump << "," << m.max_force_jump << ","
        << m.max_force_norm << "," << m.force_oscillation_index << "," << m.mean_torque_jump << ","
        << m.rms_torque_jump << "," << m.max_torque_jump << "," << m.max_torque_norm << ","
        << m.torque_oscillation_index << "," << m.total_stick_slip_switches << ","
        << m.max_tangential_force_ratio << "," << m.max_energy_gate_ratio << ","
        << m.max_inherited_energy_ratio << "\n";
}

static void RunScenario(const ScenarioConfig& scenario,
                        const SurfaceGraph& sphere_graph,
                        const std::vector<TransportVariant>& variants,
                        std::ofstream& frames_out,
                        std::ofstream& patches_out,
                        std::vector<ScenarioVariantSummary>& summaries) {
    const double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    FieldContactRuntimeSettings settings = MakeRuntimeSettings(dt);
    AnalyticFeatureHeightField field(scenario.scenario);

    std::vector<VariantState> states;
    states.reserve(variants.size());
    for (auto variant : variants) {
        VariantState state;
        state.variant = variant;
        state.runtime_tracker.Reset();
        state.metrics.Reset();
        states.push_back(state);
    }

    ChVector3d start(scenario.x0, scenario.bottom_y + scenario.sphere_radius, scenario.z0);
    ChVector3d end(scenario.x1, scenario.bottom_y + scenario.sphere_radius, scenario.z1);
    ChVector3d velocity = (end - start) / scenario.total_time;

    std::cout << "Scenario " << scenario.name << " (" << scenario.frames << " frames)" << std::endl;

    for (int frame = 0; frame < scenario.frames; frame++) {
        double u = static_cast<double>(frame) / static_cast<double>(scenario.frames - 1);
        double time = u * scenario.total_time;
        ChVector3d center = Interpolate(start, end, u);
        std::vector<FieldSampleQuery> queries = BuildSphereQueries(sphere_graph, field, center, velocity);

        for (auto& state : states) {
            FieldContactStepResult step = EvaluateVariant(state, sphere_graph, queries, center, settings);
            state.metrics.Accumulate(step);

            std::string scenario_name = scenario.name;
            std::string variant_name = VariantName(state.variant);
            WriteFrameRow(frames_out, scenario_name, variant_name, frame, time, center, velocity, step);
            WritePatchRows(patches_out, scenario_name, variant_name, frame, time, step);
        }
    }

    for (auto& state : states) {
        ScenarioVariantSummary summary;
        summary.scenario = scenario.name;
        summary.variant = VariantName(state.variant);
        summary.metrics = state.metrics.GetSummary();
        summaries.push_back(summary);

        std::cout << "  " << summary.variant << ": active=" << summary.metrics.active_frames
                  << " persistence=" << summary.metrics.primitive_persistence_ratio
                  << " force_osc=" << summary.metrics.force_oscillation_index
                  << " max_ratio=" << summary.metrics.max_tangential_force_ratio
                  << " stick_slip=" << summary.metrics.total_stick_slip_switches << std::endl;
    }
}

}  // namespace

int main() {
    std::cout << "Milestone 20: field contact feature-switching scenes" << std::endl;

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_20";
    std::filesystem::create_directories(out_dir);

    std::ofstream frames_out(out_dir / "field_contact_feature_switching_frames.csv");
    std::ofstream patches_out(out_dir / "field_contact_feature_switching_patches.csv");
    std::ofstream summary_out(out_dir / "field_contact_feature_switching_summary.csv");

    frames_out << std::fixed << std::setprecision(8);
    patches_out << std::fixed << std::setprecision(8);
    summary_out << std::fixed << std::setprecision(8);

    frames_out << "scenario,variant,frame,time,center_x,center_y,center_z,velocity_x,velocity_y,velocity_z,"
                  "patch_count,newborn_count,merge_count,split_count,death_count,max_source_count,"
                  "max_previous_reuse,total_force_x,total_force_y,total_force_z,total_force_norm,"
                  "total_torque_x,total_torque_y,total_torque_z,total_torque_norm,"
                  "max_tangential_force_ratio,max_inherited_energy_ratio\n";

    patches_out << "scenario,variant,frame,time,patch_index,persistent_id,event,source_count,source_weight_sum,"
                   "sample_count,area,mean_phi,max_penetration,center_x,center_y,center_z,normal_x,normal_y,"
                   "normal_z,normal_force_norm,tangential_force_norm,tangential_force_ratio,energy_gate_ratio,"
                   "inherited_energy_ratio,stick_slip,sources\n";

    summary_out << "scenario,variant,frames,active_frames,total_patch_count,mean_patch_count,max_patch_count,"
                   "total_newborn,total_death,total_merge_patches,total_split_patches,max_source_count,"
                   "max_previous_reuse,primitive_persistence_ratio,mean_source_weight_sum,"
                   "mean_abs_patch_count_change,max_abs_patch_count_change,mean_normal_jump_angle,"
                   "rms_normal_jump_angle,max_normal_jump_angle,mean_force_jump,rms_force_jump,max_force_jump,"
                   "max_force_norm,force_oscillation_index,mean_torque_jump,rms_torque_jump,max_torque_jump,"
                   "max_torque_norm,torque_oscillation_index,total_stick_slip_switches,"
                   "max_tangential_force_ratio,max_energy_gate_ratio,max_inherited_energy_ratio\n";

    SurfaceGraph sphere_graph = MakeSphereSurfaceGraph(0.16, 24, 48);

    std::vector<ScenarioConfig> scenarios = {
        {FeatureScenario::FoldedSeam, "folded_seam", 0.16, 0.004, -0.34, 0.34, 0.0, 0.0, 1.20, 601},
        {FeatureScenario::ConcavePolylineCorner, "concave_polyline_corner", 0.16, -0.016, -0.31, 0.31, 0.0, 0.0, 1.20, 601},
        {FeatureScenario::NarrowGrooveEntrance, "narrow_groove_entrance", 0.16, -0.031, -0.38, 0.38, 0.0, 0.0, 1.35, 676},
    };

    std::vector<TransportVariant> variants = {
        TransportVariant::RuntimeFieldPrimitive,
        TransportVariant::NoHistory,
        TransportVariant::DirectInherit,
        TransportVariant::ProjectionTransport,
        TransportVariant::MinimalRotationGateSplitMerge,
    };

    std::vector<ScenarioVariantSummary> summaries;
    for (const auto& scenario : scenarios) {
        RunScenario(scenario, sphere_graph, variants, frames_out, patches_out, summaries);
    }

    for (const auto& summary : summaries) {
        WriteSummaryRow(summary_out, summary);
    }

    std::cout << "Wrote:" << std::endl;
    std::cout << "  " << (out_dir / "field_contact_feature_switching_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "field_contact_feature_switching_patches.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "field_contact_feature_switching_summary.csv").string() << std::endl;
    return 0;
}
