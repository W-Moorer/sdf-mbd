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

#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactPrimitives.h"

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
