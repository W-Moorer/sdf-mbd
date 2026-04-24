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
// Regression tests for field-based contact primitive invariants.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"
#include "chrono/core/ChRotation.h"
#include "chrono/physics/ChBody.h"

#include "gtest/gtest.h"

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

SurfaceGraph MakeIndexedGraph(int count) {
    SurfaceGraph graph;
    graph.samples.resize(count);
    for (int i = 0; i < count; i++) {
        graph.samples[i].id = i;
        graph.samples[i].area = 1.0;
        graph.samples[i].local_pos = ChVector3d(static_cast<double>(i), 0, 0);
        if (i > 0) {
            graph.samples[i].neighbors.push_back(i - 1);
        }
        if (i + 1 < count) {
            graph.samples[i].neighbors.push_back(i + 1);
        }
    }
    return graph;
}

PrimitivePatch MakePatch(std::vector<int> sample_ids, const ChVector3d& center = ChVector3d(0, 0, 0)) {
    PrimitivePatch patch;
    patch.sample_ids = sample_ids;
    patch.center = center;
    patch.normal = ChVector3d(0, 1, 0);
    patch.area = static_cast<double>(sample_ids.size());
    return patch;
}

PrimitiveSnapshot MakeSnapshot(int persistent_id, std::vector<int> sample_ids) {
    PrimitiveSnapshot snapshot;
    snapshot.persistent_id = persistent_id;
    snapshot.sample_ids = sample_ids;
    snapshot.center = ChVector3d(0, 0, 0);
    snapshot.normal = ChVector3d(0, 1, 0);
    snapshot.area = static_cast<double>(sample_ids.size());
    return snapshot;
}

std::vector<FieldSampleQuery> MakeQueries(const SurfaceGraph& graph) {
    std::vector<FieldSampleQuery> queries(graph.samples.size());
    for (size_t i = 0; i < graph.samples.size(); i++) {
        queries[i].phi = 0.01;
        queries[i].grad = ChVector3d(0, 1, 0);
        queries[i].world_pos = graph.samples[i].local_pos;
        queries[i].world_vel = ChVector3d(0, 0, 0);
    }
    return queries;
}

FieldContactStepResult MakeOneSidedPairStep(const ChVector3d& force,
                                            const ChVector3d& contact_point,
                                            double tangent_ratio) {
    FieldContactStepResult step;
    FieldContactPatchRuntimeResult patch_result;
    patch_result.patch.center = contact_point;
    patch_result.patch.force = force;
    patch_result.patch.normal_force = force;
    patch_result.patch.area = 1.0;
    patch_result.tangential_force_ratio = tangent_ratio;

    step.patches.push_back(patch_result);
    step.stats.patch_count = 1;
    step.stats.max_tangential_force_ratio = tangent_ratio;
    step.total_force = force;
    return step;
}

ChVector3d ProjectionOnlyTransport(const ChVector3d& xi_elastic_world_prev, const ChVector3d& normal_cur) {
    return ProjectToTangent(xi_elastic_world_prev, SafeNormalize(normal_cur, ChVector3d(0, 1, 0)));
}

double SourceWeightSum(const std::vector<HistorySource>& sources) {
    double sum = 0.0;
    for (const auto& source : sources) {
        sum += source.weight;
    }
    return sum;
}

struct ClassifiedPatch {
    std::string event;
    int persistent_id;
    std::vector<HistorySource> sources;
};

std::vector<ClassifiedPatch> ClassifyEvents(const std::vector<PrimitivePatch>& current,
                                            const std::vector<PrimitiveSnapshot>& previous,
                                            const SurfaceGraph& graph) {
    HistoryInheritanceSettings settings;
    settings.min_overlap = 0.01;
    settings.min_normal_dot = 0.1;
    settings.max_center_distance = 10.0;
    settings.geometry_fallback_weight = 0.0;

    std::vector<ClassifiedPatch> classified(current.size());
    std::map<int, std::vector<int>> primary_users;
    for (size_t i = 0; i < current.size(); i++) {
        classified[i].sources = ComputeHistorySources(current[i], previous, graph, settings);
        if (!classified[i].sources.empty()) {
            primary_users[classified[i].sources.front().persistent_id].push_back(static_cast<int>(i));
        }
    }

    int next_persistent_id = 1000;
    for (size_t i = 0; i < current.size(); i++) {
        auto& out = classified[i];
        if (out.sources.empty()) {
            out.event = "newborn";
            out.persistent_id = next_persistent_id++;
            continue;
        }

        int primary = out.sources.front().persistent_id;
        const auto& users = primary_users[primary];
        bool split_child = users.size() > 1 && users.front() != static_cast<int>(i);
        bool merge_patch = out.sources.size() > 1;

        if (split_child) {
            out.event = merge_patch ? "split_merge" : "split";
            out.persistent_id = next_persistent_id++;
        } else {
            out.persistent_id = primary;
            if (merge_patch) {
                out.event = "merge";
            } else if (users.size() > 1) {
                out.event = "split_primary";
            } else {
                out.event = "stable";
            }
        }
    }

    return classified;
}

}  // namespace

TEST(FieldContactPrimitives, TangentialUpdateRespectsCoulombDisk) {
    TangentialContactSettings settings;
    settings.stiffness = 100.0;
    settings.damping = 5.0;
    settings.friction_coefficient = 0.4;
    settings.time_step = 1.0;

    TangentialUpdateResult result = UpdateTangentialContact(
        nullptr, ChVector3d(0, 1, 0), ChVector3d(10, 0, 0), 100.0, 0.0, settings);

    ASSERT_EQ(result.state, StickSlipState::Slip);
    ASSERT_LE(result.final_force_norm, result.friction_limit + 1.0e-12);
    ASSERT_NEAR(result.friction_limit, 40.0, 1.0e-12);
    ASSERT_NEAR(result.force.Dot(ChVector3d(0, 1, 0)), 0.0, 1.0e-12);
}

TEST(FieldContactPrimitives, HistoryAggregationDoesNotAmplifyWeightedSourceEnergy) {
    const double stiffness = 2.0e5;
    ChVector3d normal(0, 1, 0);

    TangentialHistory h0;
    h0.valid = true;
    h0.normal = normal;
    h0.xi_elastic_world = ChVector3d(0.002, 0, 0);

    TangentialHistory h1;
    h1.valid = true;
    h1.normal = normal;
    h1.xi_elastic_world = ChVector3d(0, 0, 0.002);

    std::vector<WeightedTangentialHistorySource> sources = {
        {h0, 0.5},
        {h1, 0.5}
    };

    TangentialHistory aggregated = AggregateTangentialHistorySources(sources, normal, 10);

    double source_bound = 0.0;
    for (const auto& source : sources) {
        source_bound += source.weight * 0.5 * stiffness *
                        source.history.xi_elastic_world.Dot(source.history.xi_elastic_world);
    }
    double inherited_energy = 0.5 * stiffness *
                              aggregated.xi_elastic_world.Dot(aggregated.xi_elastic_world);

    ASSERT_TRUE(aggregated.valid);
    ASSERT_LE(inherited_energy, source_bound + 1.0e-12);
    ASSERT_NEAR(inherited_energy / source_bound, 0.5, 1.0e-12);
    ASSERT_NEAR(aggregated.xi_elastic_world.Dot(normal), 0.0, 1.0e-12);
}

TEST(FieldContactPrimitives, MergeEventHasTwoBoundedHistorySources) {
    SurfaceGraph graph = MakeIndexedGraph(4);

    std::vector<PrimitiveSnapshot> previous = {
        MakeSnapshot(0, {0, 1}),
        MakeSnapshot(1, {2, 3})
    };
    std::vector<PrimitivePatch> current = {
        MakePatch({0, 1, 2, 3})
    };

    auto classified = ClassifyEvents(current, previous, graph);

    ASSERT_EQ(classified.size(), 1);
    ASSERT_EQ(classified[0].event, "merge");
    ASSERT_EQ(classified[0].persistent_id, 0);
    ASSERT_EQ(classified[0].sources.size(), 2);
    ASSERT_NEAR(SourceWeightSum(classified[0].sources), 1.0, 1.0e-12);
    ASSERT_LE(SourceWeightSum(classified[0].sources), 1.0 + 1.0e-12);
}

TEST(FieldContactPrimitives, SplitEventReusesOnePreviousPrimitiveWithoutAmplifyingWeights) {
    SurfaceGraph graph = MakeIndexedGraph(4);

    std::vector<PrimitiveSnapshot> previous = {
        MakeSnapshot(7, {0, 1, 2, 3})
    };
    std::vector<PrimitivePatch> current = {
        MakePatch({0, 1}, ChVector3d(-1, 0, 0)),
        MakePatch({2, 3}, ChVector3d(1, 0, 0))
    };

    auto classified = ClassifyEvents(current, previous, graph);

    ASSERT_EQ(classified.size(), 2);
    ASSERT_EQ(classified[0].event, "split_primary");
    ASSERT_EQ(classified[1].event, "split");
    ASSERT_EQ(classified[0].persistent_id, 7);
    ASSERT_NE(classified[1].persistent_id, 7);
    ASSERT_EQ(classified[0].sources.size(), 1);
    ASSERT_EQ(classified[1].sources.size(), 1);
    ASSERT_EQ(classified[0].sources[0].persistent_id, 7);
    ASSERT_EQ(classified[1].sources[0].persistent_id, 7);
    ASSERT_NEAR(classified[0].sources[0].weight, 0.5, 1.0e-12);
    ASSERT_NEAR(classified[1].sources[0].weight, 0.5, 1.0e-12);
    ASSERT_LE(classified[0].sources[0].weight + classified[1].sources[0].weight, 1.0 + 1.0e-12);
}

TEST(FieldContactPrimitives, RuntimeTrackerClassifiesMergeAndSplitEvents) {
    SurfaceGraph graph = MakeIndexedGraph(4);
    std::vector<FieldSampleQuery> queries = MakeQueries(graph);

    FieldContactRuntimeSettings settings;
    settings.inheritance.min_overlap = 0.01;
    settings.inheritance.min_normal_dot = 0.1;
    settings.inheritance.max_center_distance = 10.0;
    settings.inheritance.geometry_fallback_weight = 0.0;

    FieldContactPrimitiveTracker tracker;
    tracker.EvaluatePatches(graph,
                            queries,
                            {MakePatch({0, 1, 2}), MakePatch({3})},
                            ChVector3d(0, 0, 0),
                            settings);
    FieldContactStepResult merge =
        tracker.EvaluatePatches(graph, queries, {MakePatch({0, 1, 2, 3})}, ChVector3d(0, 0, 0), settings);

    ASSERT_EQ(merge.stats.patch_count, 1);
    ASSERT_EQ(merge.stats.merge_count, 1);
    ASSERT_EQ(merge.patches[0].event, FieldContactPrimitiveEvent::Merge);
    ASSERT_EQ(merge.patches[0].sources.size(), 2);
    ASSERT_LE(merge.patches[0].source_weight_sum, 1.0 + 1.0e-12);

    tracker.Reset();
    tracker.EvaluatePatches(graph, queries, {MakePatch({0, 1, 2, 3})}, ChVector3d(0, 0, 0), settings);
    FieldContactStepResult split =
        tracker.EvaluatePatches(graph,
                                queries,
                                {MakePatch({0, 1, 2}), MakePatch({3})},
                                ChVector3d(0, 0, 0),
                                settings);

    ASSERT_EQ(split.stats.patch_count, 2);
    ASSERT_EQ(split.stats.split_count, 2);
    ASSERT_EQ(split.patches[0].event, FieldContactPrimitiveEvent::SplitPrimary);
    ASSERT_EQ(split.patches[1].event, FieldContactPrimitiveEvent::Split);
    ASSERT_EQ(split.patches[0].sources[0].persistent_id, split.patches[1].sources[0].persistent_id);
    ASSERT_LE(split.patches[0].source_weight_sum + split.patches[1].source_weight_sum, 1.0 + 1.0e-12);
}

TEST(FieldContactPrimitives, TangentialUpdateIsGloballyRotationCovariant) {
    TangentialContactSettings settings;
    settings.stiffness = 250.0;
    settings.damping = 3.0;
    settings.friction_coefficient = 0.8;
    settings.time_step = 0.01;

    ChVector3d previous_normal = SafeNormalize(ChVector3d(0.15, 0.96, -0.23), ChVector3d(0, 1, 0));
    ChVector3d current_normal = SafeNormalize(ChVector3d(-0.32, 0.89, 0.18), ChVector3d(0, 1, 0));

    TangentialHistory previous;
    previous.valid = true;
    previous.normal = previous_normal;
    previous.xi_elastic_world = ProjectToTangent(ChVector3d(0.012, -0.004, 0.006), previous_normal);

    ChVector3d vt = ProjectToTangent(ChVector3d(0.08, -0.03, 0.05), current_normal);
    TangentialUpdateResult base =
        UpdateTangentialContact(&previous, current_normal, vt, 100.0, 0.85, settings);

    ChQuaterniond q = QuatFromAngleAxis(0.73, SafeNormalize(ChVector3d(0.31, -0.72, 0.54), ChVector3d(1, 0, 0)));

    TangentialHistory rotated_previous;
    rotated_previous.valid = true;
    rotated_previous.normal = q.Rotate(previous.normal);
    rotated_previous.xi_elastic_world = q.Rotate(previous.xi_elastic_world);

    TangentialUpdateResult rotated =
        UpdateTangentialContact(&rotated_previous, q.Rotate(current_normal), q.Rotate(vt), 100.0, 0.85, settings);

    ASSERT_NEAR((rotated.transported_xi - q.Rotate(base.transported_xi)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR((rotated.gated_xi - q.Rotate(base.gated_xi)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR((rotated.trial_xi - q.Rotate(base.trial_xi)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR((rotated.force - q.Rotate(base.force)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR((rotated.history.xi_elastic_world - q.Rotate(base.history.xi_elastic_world)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(rotated.final_force_norm, base.final_force_norm, 1.0e-12);
    ASSERT_NEAR(rotated.stored_energy_after_gate, base.stored_energy_after_gate, 1.0e-12);
    ASSERT_EQ(rotated.state, base.state);
}

TEST(FieldContactPrimitives, RotatingNormalWithoutSlipDoesNotCreateTangentialForce) {
    TangentialContactSettings settings;
    settings.stiffness = 1000.0;
    settings.damping = 0.0;
    settings.friction_coefficient = 0.5;
    settings.time_step = 0.002;

    TangentialHistory history;
    history.valid = true;
    history.normal = ChVector3d(0, 1, 0);
    history.xi_elastic_world = ChVector3d(0, 0, 0);

    double max_force = 0.0;
    double max_history_norm = 0.0;
    for (int i = 1; i <= 120; i++) {
        double angle = 0.35 * static_cast<double>(i) / 120.0;
        ChQuaterniond q = QuatFromAngleAxis(angle, ChVector3d(0, 0, 1));
        ChVector3d normal = q.Rotate(ChVector3d(0, 1, 0));
        TangentialUpdateResult result =
            UpdateTangentialContact(&history, normal, ChVector3d(0, 0, 0), 20.0, 1.0, settings);

        max_force = std::max(max_force, result.final_force_norm);
        max_history_norm = std::max(max_history_norm, result.history.xi_elastic_world.Length());
        history = result.history;
    }

    ASSERT_LE(max_force, 1.0e-12);
    ASSERT_LE(max_history_norm, 1.0e-15);
}

TEST(FieldContactPrimitives, MinimalRotationTransportBeatsProjectionOnlyTransport) {
    ChVector3d normal_prev(0, 1, 0);
    ChVector3d xi_prev(0.001, 0, 0);
    ChQuaterniond q = QuatFromAngleAxis(30.0 * CH_DEG_TO_RAD, ChVector3d(0, 0, 1));
    ChVector3d normal_cur = q.Rotate(normal_prev);
    ChVector3d expected = q.Rotate(xi_prev);

    ChVector3d minimal = TransportElasticStateMinimalRotation(xi_prev, normal_prev, normal_cur);
    ChVector3d projected = ProjectionOnlyTransport(xi_prev, normal_cur);

    double xi_norm = expected.Length();
    double minimal_error_ratio = (minimal - expected).Length() / xi_norm;
    double projection_error_ratio = (projected - expected).Length() / xi_norm;
    double projection_energy_ratio = projected.Dot(projected) / expected.Dot(expected);

    ASSERT_LE(minimal_error_ratio, 1.0e-12);
    ASSERT_GE(projection_error_ratio, 0.10);
    ASSERT_LT(projection_energy_ratio, 0.90);
    ASSERT_GT(projection_error_ratio, 1.0e8 * std::max(minimal_error_ratio, 1.0e-16));
}

TEST(FieldContactPrimitives, BidirectionalPairConservesAndDeduplicatesForces) {
    ChVector3d reference_a(-0.4, 0.1, -0.2);
    ChVector3d reference_b(0.6, -0.15, 0.25);
    ChVector3d contact_point(0.1, 0.05, 0.2);

    FieldContactStepResult a_to_b = MakeOneSidedPairStep(ChVector3d(12, 30, -4), contact_point, 0.72);
    FieldContactStepResult b_to_a = MakeOneSidedPairStep(ChVector3d(-10, -34, 3), contact_point, 0.91);

    FieldContactPairResult pair =
        CombineBidirectionalFieldContactPair(a_to_b, reference_a, b_to_a, reference_b);

    ASSERT_NEAR(pair.force_residual.Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(pair.torque_residual.Length(), 0.0, 1.0e-12);
    ASSERT_LE(pair.force_amplification_ratio, 1.0 + 1.0e-12);
    ASSERT_LT(pair.symmetric_force_norm, a_to_b.total_force.Length() + b_to_a.total_force.Length());
    ASSERT_NEAR(pair.max_tangential_force_ratio, 0.91, 1.0e-12);
    ASSERT_LE(pair.max_tangential_force_ratio, 1.0 + 1.0e-12);
}

TEST(FieldContactPrimitives, BidirectionalRuntimeEvaluatesBothSurfaceFieldDirections) {
    SurfaceGraph graph_a = MakeIndexedGraph(1);
    SurfaceGraph graph_b = MakeIndexedGraph(1);
    ChVector3d reference_a(-0.5, 0, 0);
    ChVector3d reference_b(0.5, 0, 0);
    ChVector3d contact_point(0, 0.1, 0);

    std::vector<FieldSampleQuery> queries_a(1);
    queries_a[0].phi = -0.01;
    queries_a[0].grad = ChVector3d(1, 0, 0);
    queries_a[0].world_pos = contact_point;
    queries_a[0].world_vel = ChVector3d(0, 0, 0);

    std::vector<FieldSampleQuery> queries_b(1);
    queries_b[0].phi = -0.01;
    queries_b[0].grad = ChVector3d(-1, 0, 0);
    queries_b[0].world_pos = contact_point;
    queries_b[0].world_vel = ChVector3d(0, 0, 0);

    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = 0.0;
    settings.extraction.min_samples = 1;
    settings.extraction.min_area = 0.0;
    settings.normal.stiffness = 1000.0;
    settings.normal.damping = 0.0;
    settings.tangential.stiffness = 100.0;
    settings.tangential.damping = 0.0;
    settings.tangential.friction_coefficient = 0.5;
    settings.tangential.time_step = 0.001;

    FieldContactPrimitiveTracker tracker_a;
    FieldContactPrimitiveTracker tracker_b;
    FieldContactStepResult a_to_b = tracker_a.Evaluate(graph_a, queries_a, reference_a, settings);
    FieldContactStepResult b_to_a = tracker_b.Evaluate(graph_b, queries_b, reference_b, settings);
    FieldContactPairResult pair =
        CombineBidirectionalFieldContactPair(a_to_b, reference_a, b_to_a, reference_b);

    ASSERT_EQ(a_to_b.stats.patch_count, 1);
    ASSERT_EQ(b_to_a.stats.patch_count, 1);
    ASSERT_NEAR((a_to_b.total_force - ChVector3d(10, 0, 0)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR((b_to_a.total_force - ChVector3d(-10, 0, 0)).Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(pair.force_residual.Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(pair.torque_residual.Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(pair.force_amplification_ratio, 1.0, 1.0e-12);
    ASSERT_LE(pair.max_tangential_force_ratio, 1.0 + 1.0e-12);
}

TEST(FieldContactPrimitives, SymmetricPairWrenchCanBeAppliedToTwoBodies) {
    ChVector3d reference_a(-0.25, 0, 0);
    ChVector3d reference_b(0.25, 0, 0);
    ChVector3d contact_point(0.0, 0.2, -0.1);

    FieldContactStepResult a_to_b = MakeOneSidedPairStep(ChVector3d(0, 18, 4), contact_point, 0.4);
    FieldContactStepResult b_to_a = MakeOneSidedPairStep(ChVector3d(0, -18, -4), contact_point, 0.5);
    FieldContactPairResult pair =
        CombineBidirectionalFieldContactPair(a_to_b, reference_a, b_to_a, reference_b);

    ChBody body_a;
    ChBody body_b;
    body_a.SetRot(ChQuaterniond(1, 0, 0, 0));
    body_b.SetRot(ChQuaterniond(1, 0, 0, 0));
    body_a.SetPos(reference_a);
    body_b.SetPos(reference_b);

    unsigned int acc_a = body_a.AddAccumulator();
    unsigned int acc_b = body_b.AddAccumulator();
    body_a.EmptyAccumulator(acc_a);
    body_b.EmptyAccumulator(acc_b);

    body_a.AccumulateForce(acc_a, pair.on_a.force, reference_a, false);
    body_a.AccumulateTorque(acc_a, pair.on_a.torque, false);
    body_b.AccumulateForce(acc_b, pair.on_b.force, reference_b, false);
    body_b.AccumulateTorque(acc_b, pair.on_b.torque, false);

    ChVector3d accumulated_force_residual = body_a.GetAccumulatedForce(acc_a) + body_b.GetAccumulatedForce(acc_b);
    ChVector3d pair_reference = 0.5 * (reference_a + reference_b);
    ChVector3d accumulated_torque_residual =
        body_a.GetAccumulatedTorque(acc_a) + (reference_a - pair_reference).Cross(body_a.GetAccumulatedForce(acc_a)) +
        body_b.GetAccumulatedTorque(acc_b) + (reference_b - pair_reference).Cross(body_b.GetAccumulatedForce(acc_b));

    ASSERT_NEAR(accumulated_force_residual.Length(), 0.0, 1.0e-12);
    ASSERT_NEAR(accumulated_torque_residual.Length(), 0.0, 1.0e-12);
}
