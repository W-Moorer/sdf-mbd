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
// Milestone 21: multi-point sliding and interlock benchmarks.
//
// This executable extends the Milestone 20 analytic feature-switching demo to
// harder non-convex scenes:
//   - groove/guide-rail sliding;
//   - nested/interlocking lips;
//   - multi-patch rolling/sliding over staggered pads.
//
// It compares complete field primitives against two deliberately fragmented
// baselines. The comparison surface is the Milestone 19 topology/response
// metric summary: contact-set churn, normal jumps, force/torque oscillation,
// Coulomb feasibility, and inherited-history energy.
//
// Outputs:
//   out/milestone_21/field_contact_multislip_interlock_frames.csv
//   out/milestone_21/field_contact_multislip_interlock_patches.csv
//   out/milestone_21/field_contact_multislip_interlock_summary.csv
//   out/milestone_21/field_contact_multislip_interlock_comparison.csv
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

enum class NonconvexScenario {
    GuideRailSliding,
    NestedInterlock,
    MultiPatchRollingSliding
};

enum class BenchmarkVariant {
    FieldPrimitive,
    RawSampleContacts,
    NormalBinFragments
};

struct ScenarioConfig {
    NonconvexScenario scenario;
    std::string name;
    double sphere_radius = 0.16;
    double bottom_y = 0.0;
    ChVector3d start = ChVector3d(0, 0, 0);
    ChVector3d end = ChVector3d(0, 0, 0);
    double z_oscillation = 0.0;
    double z_oscillation_cycles = 1.0;
    double total_time = 1.4;
    int frames = 701;
    double rolling_scale = 1.0;
};

struct BodyKinematics {
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
    double roll_angle_x = 0.0;
    double roll_angle_z = 0.0;
    ChVector3d angular_velocity = ChVector3d(0, 0, 0);
};

struct VariantState {
    BenchmarkVariant variant;
    FieldContactPrimitiveTracker tracker;
    FieldContactTopologyMetricsAccumulator metrics;
};

struct ScenarioVariantSummary {
    std::string scenario;
    std::string variant;
    FieldContactTopologyMetricsSummary metrics;
};

static const char* ScenarioName(NonconvexScenario scenario) {
    switch (scenario) {
        case NonconvexScenario::GuideRailSliding:
            return "guide_rail_sliding";
        case NonconvexScenario::NestedInterlock:
            return "nested_interlock";
        case NonconvexScenario::MultiPatchRollingSliding:
            return "multi_patch_rolling_sliding";
    }
    return "unknown";
}

static const char* VariantName(BenchmarkVariant variant) {
    switch (variant) {
        case BenchmarkVariant::FieldPrimitive:
            return "field_primitive";
        case BenchmarkVariant::RawSampleContacts:
            return "raw_sample_contacts";
        case BenchmarkVariant::NormalBinFragments:
            return "normal_bin_fragments";
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

static double Gaussian(double x, double z, double sx, double sz) {
    double qx = x / std::max(sx, 1.0e-9);
    double qz = z / std::max(sz, 1.0e-9);
    return std::exp(-0.5 * (qx * qx + qz * qz));
}

class AnalyticNonconvexHeightField {
  public:
    explicit AnalyticNonconvexHeightField(NonconvexScenario scenario) : m_scenario(scenario) {}

    FieldSampleQuery Query(const ChVector3d& world_pos, const ChVector3d& world_vel) const {
        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = world_vel;
        query.phi = world_pos.y() - Height(world_pos.x(), world_pos.z());

        const double h = 2.0e-4;
        double dhdx = (Height(world_pos.x() + h, world_pos.z()) -
                       Height(world_pos.x() - h, world_pos.z())) /
                      (2.0 * h);
        double dhdz = (Height(world_pos.x(), world_pos.z() + h) -
                       Height(world_pos.x(), world_pos.z() - h)) /
                      (2.0 * h);
        query.grad = SafeNormalize(ChVector3d(-dhdx, 1.0, -dhdz), ChVector3d(0, 1, 0));
        return query;
    }

  private:
    double Height(double x, double z) const {
        switch (m_scenario) {
            case NonconvexScenario::GuideRailSliding:
                return GuideRailHeight(x, z);
            case NonconvexScenario::NestedInterlock:
                return NestedInterlockHeight(x, z);
            case NonconvexScenario::MultiPatchRollingSliding:
                return MultiPatchPadHeight(x, z);
        }
        return 0.0;
    }

    static double GuideRailHeight(double x, double z) {
        double left_rail = Gaussian(0.0, z - 0.086, 1.0, 0.020);
        double right_rail = Gaussian(0.0, z + 0.086, 1.0, 0.020);
        double rail = left_rail + right_rail;
        double center_trough = std::exp(-std::pow(z / 0.060, 4.0));
        double rail_waviness = 0.0035 * std::sin(2.0 * kPi * x / 0.22) * rail;
        double local_gate = 0.006 * Gaussian(x, z, 0.095, 0.14);
        return 0.004 + 0.029 * rail - 0.016 * center_trough + rail_waviness + local_gate;
    }

    static double NestedInterlockHeight(double x, double z) {
        double slot = std::exp(-std::pow(z / 0.070, 4.0));
        double collar = Gaussian(x, z - 0.095, 0.065, 0.023) + Gaussian(x, z + 0.095, 0.065, 0.023);
        double inner_lip = (Gaussian(x - 0.105, z, 0.050, 0.032) + Gaussian(x + 0.105, z, 0.050, 0.032));
        double stagger = 0.004 * std::sin(2.0 * kPi * x / 0.18) * std::exp(-std::pow(z / 0.12, 2.0));
        return 0.002 - 0.024 * slot + 0.030 * collar + 0.014 * inner_lip + stagger;
    }

    static double MultiPatchPadHeight(double x, double z) {
        double height = 0.002;
        for (int i = -4; i <= 4; i++) {
            double xi = 0.105 * static_cast<double>(i);
            double zi = (i % 2 == 0) ? 0.072 : -0.072;
            height += 0.033 * Gaussian(x - xi, z - zi, 0.023, 0.024);
            height += 0.027 * Gaussian(x - xi, z + 0.62 * zi, 0.026, 0.022);
        }
        height += 0.002 * std::sin(2.0 * kPi * x / 0.17) * std::cos(2.0 * kPi * z / 0.24);
        return height;
    }

    NonconvexScenario m_scenario;
};

static FieldContactRuntimeSettings MakeRuntimeSettings(double time_step) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = 0.004;
    settings.extraction.min_area = 1.0e-8;
    settings.extraction.min_samples = 3;
    settings.extraction.use_penetration_weighted_center = true;

    settings.normal.stiffness = 1.8e6;
    settings.normal.damping = 0.0;

    settings.tangential.stiffness = 3.0e3;
    settings.tangential.damping = 0.0;
    settings.tangential.friction_coefficient = 0.45;
    settings.tangential.time_step = time_step;

    settings.inheritance.min_overlap = 0.006;
    settings.inheritance.min_normal_dot = 0.18;
    settings.inheritance.max_center_distance = 0.085;
    settings.inheritance.geometry_fallback_weight = 0.20;
    return settings;
}

static ChVector3d RotateX(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(v.x(), c * v.y() - s * v.z(), s * v.y() + c * v.z());
}

static ChVector3d RotateZ(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z());
}

static ChVector3d RotateLocal(const ChVector3d& v, double angle_x, double angle_z) {
    return RotateX(RotateZ(v, angle_z), angle_x);
}

static BodyKinematics KinematicsAtFrame(const ScenarioConfig& scenario, int frame) {
    double u = static_cast<double>(frame) / static_cast<double>(scenario.frames - 1);
    double time_step = scenario.total_time / static_cast<double>(scenario.frames - 1);
    ChVector3d start = scenario.start + ChVector3d(0, scenario.sphere_radius + scenario.bottom_y, 0);
    ChVector3d end = scenario.end + ChVector3d(0, scenario.sphere_radius + scenario.bottom_y, 0);

    auto center_at = [&](double uu) {
        double osc = scenario.z_oscillation * std::sin(2.0 * kPi * scenario.z_oscillation_cycles * uu);
        return start + (end - start) * uu + ChVector3d(0, 0, osc);
    };

    double u0 = std::max(0.0, u - 1.0 / static_cast<double>(scenario.frames - 1));
    double u1 = std::min(1.0, u + 1.0 / static_cast<double>(scenario.frames - 1));

    BodyKinematics kin;
    kin.center = center_at(u);
    kin.velocity = (center_at(u1) - center_at(u0)) / ((u1 - u0) * scenario.total_time);
    kin.roll_angle_z = -scenario.rolling_scale * (kin.center.x() - start.x()) / scenario.sphere_radius;
    kin.roll_angle_x = scenario.rolling_scale * (kin.center.z() - start.z()) / scenario.sphere_radius;
    kin.angular_velocity = ChVector3d(kin.velocity.z() / scenario.sphere_radius,
                                      0.0,
                                      -kin.velocity.x() / scenario.sphere_radius) *
                           scenario.rolling_scale;
    (void)time_step;
    return kin;
}

static std::vector<FieldSampleQuery> BuildQueries(const SurfaceGraph& graph,
                                                  const AnalyticNonconvexHeightField& field,
                                                  const BodyKinematics& kin) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(graph.samples.size());
    for (const auto& sample : graph.samples) {
        ChVector3d rotated = RotateLocal(sample.local_pos, kin.roll_angle_x, kin.roll_angle_z);
        ChVector3d world_pos = kin.center + rotated;
        ChVector3d world_vel = kin.velocity + kin.angular_velocity.Cross(rotated);
        queries.push_back(field.Query(world_pos, world_vel));
    }
    return queries;
}

static PrimitivePatch BuildPatchFromSampleIds(const SurfaceGraph& graph,
                                              const std::vector<FieldSampleQuery>& queries,
                                              const std::vector<int>& sample_ids,
                                              int primitive_id,
                                              bool penetration_weighted_center) {
    PrimitivePatch primitive;
    primitive.primitive_id = primitive_id;
    primitive.sample_ids = sample_ids;

    double area = 0.0;
    double area_weighted_phi = 0.0;
    double area_weighted_penetration = 0.0;
    double center_weight = 0.0;
    ChVector3d area_center(0, 0, 0);
    ChVector3d penetration_center(0, 0, 0);
    ChVector3d velocity_sum(0, 0, 0);
    ChVector3d normal_sum(0, 0, 0);
    double min_phi = 1.0e100;
    double max_penetration = 0.0;

    for (int sid : sample_ids) {
        if (sid < 0 || sid >= static_cast<int>(graph.samples.size()) || sid >= static_cast<int>(queries.size())) {
            continue;
        }

        const auto& sample = graph.samples[sid];
        const auto& query = queries[sid];
        double sample_area = std::max(0.0, sample.area);
        double penetration = std::max(-query.phi, 0.0);
        ChVector3d sample_normal = SafeNormalize(query.grad, ChVector3d(0, 1, 0));

        area += sample_area;
        area_weighted_phi += query.phi * sample_area;
        area_weighted_penetration += penetration * sample_area;
        area_center += query.world_pos * sample_area;
        velocity_sum += query.world_vel * sample_area;
        normal_sum += sample_normal * sample_area;
        min_phi = std::min(min_phi, query.phi);
        max_penetration = std::max(max_penetration, penetration);

        if (penetration_weighted_center) {
            double w = sample_area * penetration;
            penetration_center += query.world_pos * w;
            center_weight += w;
        }
    }

    if (area > 1.0e-16) {
        primitive.area = area;
        primitive.center = center_weight > 1.0e-16 ? penetration_center / center_weight : area_center / area;
        primitive.representative_velocity = velocity_sum / area;
        primitive.normal = SafeNormalize(normal_sum, ChVector3d(0, 1, 0));
        BuildTangentBasis(primitive.normal, primitive.tangent_t1, primitive.tangent_t2);
        primitive.mean_phi = area_weighted_phi / area;
        primitive.min_phi = min_phi;
        primitive.max_penetration = max_penetration;
        primitive.mean_penetration = area_weighted_penetration / area;
    }

    return primitive;
}

static std::vector<PrimitivePatch> BuildRawSamplePatches(const SurfaceGraph& graph,
                                                         const std::vector<FieldSampleQuery>& queries,
                                                         const PatchExtractionSettings& settings) {
    std::vector<PrimitivePatch> patches;
    std::vector<int> active = BuildActiveSet(queries, settings.activation_band);
    patches.reserve(active.size());
    int primitive_id = 0;
    for (int sid : active) {
        patches.push_back(BuildPatchFromSampleIds(graph, queries, std::vector<int>{sid}, primitive_id++, false));
    }
    return patches;
}

static int QuantizeNormalComponent(double value) {
    const int bins = 5;
    const double tangent_range = 0.65;
    double t = Clamp01(0.5 * (value / tangent_range + 1.0));
    int id = static_cast<int>(std::floor(t * bins));
    return std::max(0, std::min(bins - 1, id));
}

static double NormalBinCenter(int id) {
    const int bins = 5;
    const double tangent_range = 0.65;
    return (((static_cast<double>(id) + 0.5) / static_cast<double>(bins)) * 2.0 - 1.0) * tangent_range;
}

static ChVector3d SnapNormalToBin(const ChVector3d& normal) {
    ChVector3d n = SafeNormalize(normal, ChVector3d(0, 1, 0));
    int bx = QuantizeNormalComponent(n.x());
    int bz = QuantizeNormalComponent(n.z());
    double nx = NormalBinCenter(bx);
    double nz = NormalBinCenter(bz);
    double tangent_sq = nx * nx + nz * nz;
    if (tangent_sq > 0.72) {
        double scale = std::sqrt(0.72 / tangent_sq);
        nx *= scale;
        nz *= scale;
        tangent_sq = nx * nx + nz * nz;
    }
    double ny = std::sqrt(std::max(0.0, 1.0 - tangent_sq));
    return SafeNormalize(ChVector3d(nx, ny, nz), ChVector3d(0, 1, 0));
}

static int NormalBinKey(const ChVector3d& normal) {
    const int bins = 5;
    ChVector3d n = SafeNormalize(normal, ChVector3d(0, 1, 0));
    int bx = QuantizeNormalComponent(n.x());
    int bz = QuantizeNormalComponent(n.z());
    return bx + bins * bz;
}

static std::vector<PrimitivePatch> BuildNormalBinPatches(const SurfaceGraph& graph,
                                                         const std::vector<FieldSampleQuery>& queries,
                                                         const PatchExtractionSettings& settings) {
    std::vector<PrimitivePatch> patches;
    std::vector<int> active = BuildActiveSet(queries, settings.activation_band);
    std::vector<char> active_mask(graph.samples.size(), 0);
    std::vector<int> bin(graph.samples.size(), -1);
    for (int sid : active) {
        if (sid >= 0 && sid < static_cast<int>(graph.samples.size())) {
            active_mask[sid] = 1;
            bin[sid] = NormalBinKey(SafeNormalize(queries[sid].grad, ChVector3d(0, 1, 0)));
        }
    }

    std::vector<char> visited(graph.samples.size(), 0);
    int primitive_id = 0;
    for (int seed : active) {
        if (seed < 0 || seed >= static_cast<int>(graph.samples.size()) || visited[seed] || !active_mask[seed]) {
            continue;
        }

        std::vector<int> component;
        std::vector<int> queue;
        queue.push_back(seed);
        visited[seed] = 1;
        size_t head = 0;
        while (head < queue.size()) {
            int current = queue[head++];
            component.push_back(current);
            for (int neighbor : graph.samples[current].neighbors) {
                if (neighbor < 0 || neighbor >= static_cast<int>(graph.samples.size())) {
                    continue;
                }
                if (!visited[neighbor] && active_mask[neighbor] && bin[neighbor] == bin[seed]) {
                    visited[neighbor] = 1;
                    queue.push_back(neighbor);
                }
            }
        }

        if (!component.empty()) {
            patches.push_back(BuildPatchFromSampleIds(graph, queries, component, primitive_id++, true));
        }
    }

    return patches;
}

static FieldContactStepResult EvaluateVariant(VariantState& state,
                                              const SurfaceGraph& graph,
                                              const std::vector<FieldSampleQuery>& queries,
                                              const ChVector3d& torque_reference,
                                              const FieldContactRuntimeSettings& settings) {
    if (state.variant == BenchmarkVariant::FieldPrimitive) {
        return state.tracker.Evaluate(graph, queries, torque_reference, settings);
    }

    std::vector<PrimitivePatch> patches;
    if (state.variant == BenchmarkVariant::RawSampleContacts) {
        patches = BuildRawSamplePatches(graph, queries, settings.extraction);
    } else if (state.variant == BenchmarkVariant::NormalBinFragments) {
        std::vector<FieldSampleQuery> binned_queries = queries;
        for (auto& query : binned_queries) {
            if (query.phi < settings.extraction.activation_band) {
                query.grad = SnapNormalToBin(query.grad);
            }
        }
        patches = BuildNormalBinPatches(graph, binned_queries, settings.extraction);
        return state.tracker.EvaluatePatches(graph, binned_queries, patches, torque_reference, settings);
    } else {
        patches = BuildNormalBinPatches(graph, queries, settings.extraction);
    }

    return state.tracker.EvaluatePatches(graph, queries, patches, torque_reference, settings);
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
                          const BodyKinematics& kin,
                          const FieldContactStepResult& step) {
    out << scenario << "," << variant << "," << frame << "," << time << "," << kin.center.x() << ","
        << kin.center.y() << "," << kin.center.z() << "," << kin.velocity.x() << "," << kin.velocity.y()
        << "," << kin.velocity.z() << "," << step.stats.patch_count << "," << step.stats.newborn_count
        << "," << step.stats.merge_count << "," << step.stats.split_count << "," << step.stats.death_count
        << "," << step.stats.max_source_count << "," << step.stats.max_previous_reuse << ","
        << step.total_force.x() << "," << step.total_force.y() << "," << step.total_force.z() << ","
        << step.total_force.Length() << "," << step.total_torque.x() << "," << step.total_torque.y()
        << "," << step.total_torque.z() << "," << step.total_torque.Length() << ","
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

static double SafeRatio(double numerator, double denominator) {
    return denominator > 1.0e-14 ? numerator / denominator : 0.0;
}

static void WriteComparisonRows(std::ofstream& out,
                                const std::string& scenario,
                                const FieldContactTopologyMetricsSummary& field,
                                const ScenarioVariantSummary& baseline) {
    const auto& b = baseline.metrics;
    out << scenario << "," << baseline.variant << ","
        << SafeRatio(b.mean_abs_patch_count_change, field.mean_abs_patch_count_change) << ","
        << SafeRatio(static_cast<double>(b.max_abs_patch_count_change),
                     static_cast<double>(field.max_abs_patch_count_change))
        << "," << SafeRatio(b.rms_normal_jump_angle, field.rms_normal_jump_angle) << ","
        << SafeRatio(b.force_oscillation_index, field.force_oscillation_index) << ","
        << SafeRatio(b.torque_oscillation_index, field.torque_oscillation_index) << ","
        << SafeRatio(static_cast<double>(b.total_newborn + b.total_death + b.total_merge_patches + b.total_split_patches),
                     static_cast<double>(field.total_newborn + field.total_death + field.total_merge_patches + field.total_split_patches))
        << "," << b.max_tangential_force_ratio << "," << b.max_inherited_energy_ratio << "\n";
}

static void RunScenario(const ScenarioConfig& scenario,
                        const SurfaceGraph& graph,
                        const std::vector<BenchmarkVariant>& variants,
                        std::ofstream& frames_out,
                        std::ofstream& patches_out,
                        std::vector<ScenarioVariantSummary>& summaries) {
    double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    FieldContactRuntimeSettings settings = MakeRuntimeSettings(dt);
    AnalyticNonconvexHeightField field(scenario.scenario);

    std::vector<VariantState> states;
    states.reserve(variants.size());
    for (auto variant : variants) {
        VariantState state;
        state.variant = variant;
        state.tracker.Reset();
        state.metrics.Reset();
        states.push_back(state);
    }

    std::cout << "Scenario " << scenario.name << " (" << scenario.frames << " frames)" << std::endl;

    for (int frame = 0; frame < scenario.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(scenario, frame);
        std::vector<FieldSampleQuery> queries = BuildQueries(graph, field, kin);

        for (auto& state : states) {
            FieldContactStepResult step = EvaluateVariant(state, graph, queries, kin.center, settings);
            state.metrics.Accumulate(step);

            std::string variant_name = VariantName(state.variant);
            WriteFrameRow(frames_out, scenario.name, variant_name, frame, time, kin, step);
            WritePatchRows(patches_out, scenario.name, variant_name, frame, time, step);
        }
    }

    for (auto& state : states) {
        ScenarioVariantSummary summary;
        summary.scenario = scenario.name;
        summary.variant = VariantName(state.variant);
        summary.metrics = state.metrics.GetSummary();
        summaries.push_back(summary);

        std::cout << "  " << summary.variant << ": mean_patch=" << summary.metrics.mean_patch_count
                  << " mean_dcount=" << summary.metrics.mean_abs_patch_count_change
                  << " rms_normal=" << summary.metrics.rms_normal_jump_angle
                  << " force_osc=" << summary.metrics.force_oscillation_index
                  << " torque_osc=" << summary.metrics.torque_oscillation_index
                  << " max_ratio=" << summary.metrics.max_tangential_force_ratio << std::endl;
    }
}

}  // namespace

int main() {
    std::cout << "Milestone 21: field contact multi-slip and interlock benchmarks" << std::endl;

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_21";
    std::filesystem::create_directories(out_dir);

    std::ofstream frames_out(out_dir / "field_contact_multislip_interlock_frames.csv");
    std::ofstream patches_out(out_dir / "field_contact_multislip_interlock_patches.csv");
    std::ofstream summary_out(out_dir / "field_contact_multislip_interlock_summary.csv");
    std::ofstream comparison_out(out_dir / "field_contact_multislip_interlock_comparison.csv");

    frames_out << std::fixed << std::setprecision(8);
    patches_out << std::fixed << std::setprecision(8);
    summary_out << std::fixed << std::setprecision(8);
    comparison_out << std::fixed << std::setprecision(8);

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

    comparison_out << "scenario,baseline_variant,mean_patch_count_change_ratio,max_patch_count_change_ratio,"
                      "rms_normal_jump_ratio,force_oscillation_ratio,torque_oscillation_ratio,"
                      "topology_event_ratio,max_tangential_force_ratio,max_inherited_energy_ratio\n";

    SurfaceGraph sphere_graph = MakeSphereSurfaceGraph(0.16, 26, 52);

    std::vector<ScenarioConfig> scenarios = {
        {NonconvexScenario::GuideRailSliding,
         "guide_rail_sliding",
         0.16,
         0.012,
         ChVector3d(-0.52, 0, 0.0),
         ChVector3d(0.52, 0, 0.0),
         0.010,
         2.0,
         1.40,
         701,
         1.0},
        {NonconvexScenario::NestedInterlock,
         "nested_interlock",
         0.16,
         0.000,
         ChVector3d(-0.42, 0, -0.012),
         ChVector3d(0.42, 0, 0.012),
         0.030,
         1.5,
         1.45,
         726,
         0.8},
        {NonconvexScenario::MultiPatchRollingSliding,
         "multi_patch_rolling_sliding",
         0.16,
         0.016,
         ChVector3d(-0.48, 0, -0.085),
         ChVector3d(0.48, 0, 0.085),
         0.020,
         2.5,
         1.50,
         751,
         1.0},
    };

    std::vector<BenchmarkVariant> variants = {
        BenchmarkVariant::FieldPrimitive,
        BenchmarkVariant::RawSampleContacts,
        BenchmarkVariant::NormalBinFragments,
    };

    std::vector<ScenarioVariantSummary> summaries;
    for (const auto& scenario : scenarios) {
        RunScenario(scenario, sphere_graph, variants, frames_out, patches_out, summaries);
    }

    for (const auto& summary : summaries) {
        WriteSummaryRow(summary_out, summary);
    }

    for (const auto& scenario : scenarios) {
        const FieldContactTopologyMetricsSummary* field_summary = nullptr;
        for (const auto& summary : summaries) {
            if (summary.scenario == scenario.name && summary.variant == "field_primitive") {
                field_summary = &summary.metrics;
                break;
            }
        }
        if (!field_summary) {
            continue;
        }

        for (const auto& summary : summaries) {
            if (summary.scenario == scenario.name && summary.variant != "field_primitive") {
                WriteComparisonRows(comparison_out, scenario.name, *field_summary, summary);
            }
        }
    }

    std::cout << "Wrote:" << std::endl;
    std::cout << "  " << (out_dir / "field_contact_multislip_interlock_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "field_contact_multislip_interlock_patches.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "field_contact_multislip_interlock_summary.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "field_contact_multislip_interlock_comparison.csv").string() << std::endl;
    return 0;
}
