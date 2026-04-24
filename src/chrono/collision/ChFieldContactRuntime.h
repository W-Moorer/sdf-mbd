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
// Runtime state machine for field-based contact primitives.
//
// This layer turns per-sample field queries into persistent primitive contact
// responses. It owns topology classification, persistent IDs, tangential
// history transfer, and primitive-level force aggregation. Concrete SDF
// backends remain outside this header and provide FieldSampleQuery arrays.
//
// =============================================================================

#ifndef CH_FIELD_CONTACT_RUNTIME_H
#define CH_FIELD_CONTACT_RUNTIME_H

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "chrono/collision/ChFieldContactPrimitives.h"

namespace chrono {
namespace fieldcontact {

enum class FieldContactPrimitiveEvent {
    Stable,
    Newborn,
    SplitPrimary,
    Split,
    Merge,
    SplitMerge
};

inline const char* FieldContactPrimitiveEventName(FieldContactPrimitiveEvent event) {
    switch (event) {
        case FieldContactPrimitiveEvent::Stable:
            return "stable";
        case FieldContactPrimitiveEvent::Newborn:
            return "newborn";
        case FieldContactPrimitiveEvent::SplitPrimary:
            return "split_primary";
        case FieldContactPrimitiveEvent::Split:
            return "split";
        case FieldContactPrimitiveEvent::Merge:
            return "merge";
        case FieldContactPrimitiveEvent::SplitMerge:
            return "split_merge";
    }
    return "unknown";
}

struct FieldContactRuntimeSettings {
    PatchExtractionSettings extraction;
    NormalContactSettings normal;
    TangentialContactSettings tangential;
    HistoryInheritanceSettings inheritance;
};

struct FieldContactPatchRuntimeResult {
    PrimitivePatch patch;
    std::vector<HistorySource> sources;
    int persistent_id = -1;
    FieldContactPrimitiveEvent event = FieldContactPrimitiveEvent::Stable;
    TangentialUpdateResult tangential;
    double source_weight_sum = 0.0;
    double source_energy_bound = 0.0;
    double inherited_energy = 0.0;
    double inherited_energy_ratio = 0.0;
    double tangential_force_ratio = 0.0;
    double energy_gate_ratio = 0.0;
};

struct FieldContactFrameStats {
    int patch_count = 0;
    int newborn_count = 0;
    int merge_count = 0;
    int split_count = 0;
    int death_count = 0;
    int max_source_count = 0;
    int max_previous_reuse = 0;
    double max_tangential_force_ratio = 0.0;
    double max_energy_gate_ratio = 0.0;
    double max_inherited_energy_ratio = 0.0;
};

struct FieldContactStepResult {
    std::vector<FieldContactPatchRuntimeResult> patches;
    FieldContactFrameStats stats;
    ChVector3d total_force = ChVector3d(0, 0, 0);
    ChVector3d total_torque = ChVector3d(0, 0, 0);
};

class FieldContactPrimitiveTracker {
  public:
    void Reset(int next_persistent_id = 0) {
        m_next_persistent_id = next_persistent_id;
        m_previous_snapshots.clear();
        m_history_store.clear();
    }

    int GetNextPersistentId() const {
        return m_next_persistent_id;
    }

    const std::vector<PrimitiveSnapshot>& GetPreviousSnapshots() const {
        return m_previous_snapshots;
    }

    const std::map<int, TangentialHistory>& GetHistoryStore() const {
        return m_history_store;
    }

    FieldContactStepResult Evaluate(const SurfaceGraph& graph,
                                    const std::vector<FieldSampleQuery>& queries,
                                    const ChVector3d& torque_reference,
                                    const FieldContactRuntimeSettings& settings) {
        std::vector<int> active = BuildActiveSet(queries, settings.extraction.activation_band);
        std::vector<PrimitivePatch> primitive_patches =
            ExtractPrimitives(graph, queries, active, settings.extraction);
        return EvaluatePatches(graph, queries, primitive_patches, torque_reference, settings);
    }

    FieldContactStepResult EvaluatePatches(const SurfaceGraph& graph,
                                           const std::vector<FieldSampleQuery>& queries,
                                           const std::vector<PrimitivePatch>& primitive_patches,
                                           const ChVector3d& torque_reference,
                                           const FieldContactRuntimeSettings& settings) {
        FieldContactStepResult result;
        result.stats.patch_count = static_cast<int>(primitive_patches.size());
        result.patches.resize(primitive_patches.size());

        std::map<int, std::vector<int>> primary_users;
        std::set<int> referenced_previous;

        for (size_t pi = 0; pi < primitive_patches.size(); pi++) {
            auto& patch_result = result.patches[pi];
            patch_result.patch = primitive_patches[pi];
            patch_result.sources =
                ComputeHistorySources(patch_result.patch, m_previous_snapshots, graph, settings.inheritance);
            result.stats.max_source_count =
                std::max(result.stats.max_source_count, static_cast<int>(patch_result.sources.size()));

            for (const auto& source : patch_result.sources) {
                patch_result.source_weight_sum += source.weight;
                referenced_previous.insert(source.persistent_id);
            }

            if (!patch_result.sources.empty()) {
                int primary = patch_result.sources.front().persistent_id;
                primary_users[primary].push_back(static_cast<int>(pi));
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

        ClassifyCurrentPatches(primary_users, result);

        for (const auto& prev : m_previous_snapshots) {
            if (!referenced_previous.count(prev.persistent_id)) {
                result.stats.death_count++;
            }
        }

        std::vector<PrimitiveSnapshot> current_snapshots;
        std::map<int, TangentialHistory> new_history_store;

        for (auto& patch_result : result.patches) {
            ApplyNormalContactIntegral(patch_result.patch,
                                       graph,
                                       queries,
                                       torque_reference,
                                       settings.normal);
            ApplyTangentialContact(patch_result, torque_reference, settings, new_history_store);

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

        m_previous_snapshots = current_snapshots;
        m_history_store = new_history_store;
        return result;
    }

  private:
    void ClassifyCurrentPatches(const std::map<int, std::vector<int>>& primary_users,
                                FieldContactStepResult& result) {
        for (size_t pi = 0; pi < result.patches.size(); pi++) {
            auto& patch_result = result.patches[pi];
            if (patch_result.sources.empty()) {
                patch_result.persistent_id = m_next_persistent_id++;
                patch_result.event = FieldContactPrimitiveEvent::Newborn;
                result.stats.newborn_count++;
                continue;
            }

            int primary = patch_result.sources.front().persistent_id;
            auto users_it = primary_users.find(primary);
            if (users_it == primary_users.end()) {
                patch_result.persistent_id = m_next_persistent_id++;
                patch_result.event = FieldContactPrimitiveEvent::Newborn;
                result.stats.newborn_count++;
                continue;
            }
            const auto& users = users_it->second;
            bool split_child = users.size() > 1 && users.front() != static_cast<int>(pi);
            bool merge_patch = patch_result.sources.size() > 1;

            if (split_child) {
                patch_result.persistent_id = m_next_persistent_id++;
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

    void ApplyTangentialContact(FieldContactPatchRuntimeResult& patch_result,
                                const ChVector3d& torque_reference,
                                const FieldContactRuntimeSettings& settings,
                                std::map<int, TangentialHistory>& new_history_store) const {
        std::vector<WeightedTangentialHistorySource> weighted_history_sources;
        for (const auto& source : patch_result.sources) {
            auto hist_it = m_history_store.find(source.persistent_id);
            if (hist_it == m_history_store.end() || !hist_it->second.valid) {
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
        const TangentialHistory* previous_history = aggregated_history.valid ? &aggregated_history : nullptr;

        patch_result.inherited_energy = 0.5 * settings.tangential.stiffness *
                                        aggregated_history.xi_elastic_world.Dot(aggregated_history.xi_elastic_world);
        patch_result.inherited_energy_ratio = patch_result.source_energy_bound > 1.0e-16 ?
                                                  patch_result.inherited_energy / patch_result.source_energy_bound :
                                                  0.0;

        ChVector3d vt = ProjectToTangent(patch_result.patch.representative_velocity, patch_result.patch.normal);
        patch_result.tangential =
            UpdateTangentialContact(previous_history,
                                    patch_result.patch.normal,
                                    vt,
                                    patch_result.patch.normal_force.Length(),
                                    previous_history ? 1.0 : 0.0,
                                    settings.tangential);
        patch_result.tangential.history.persistent_id = patch_result.persistent_id;
        new_history_store[patch_result.persistent_id] = patch_result.tangential.history;

        patch_result.patch.tangential_force = patch_result.tangential.force;
        patch_result.patch.force = patch_result.patch.normal_force + patch_result.patch.tangential_force;
        patch_result.patch.torque +=
            (patch_result.patch.center - torque_reference).Cross(patch_result.patch.tangential_force);

        patch_result.tangential_force_ratio = patch_result.tangential.friction_limit > 1.0e-12 ?
                                                  patch_result.tangential.final_force_norm /
                                                      patch_result.tangential.friction_limit :
                                                  0.0;
        patch_result.energy_gate_ratio = patch_result.tangential.stored_energy_before_gate > 1.0e-16 ?
                                             patch_result.tangential.stored_energy_after_gate /
                                                 patch_result.tangential.stored_energy_before_gate :
                                             0.0;
    }

    int m_next_persistent_id = 0;
    std::vector<PrimitiveSnapshot> m_previous_snapshots;
    std::map<int, TangentialHistory> m_history_store;
};

}  // namespace fieldcontact
}  // namespace chrono

#endif
