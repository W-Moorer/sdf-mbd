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
// Field-based contact primitive utilities.
//
// This header contains the geometry- and constitutive-law parts of the
// field-based contact primitive prototype. It is intentionally independent of a
// concrete SDF backend such as OpenVDB or NanoVDB: callers provide surface
// samples, adjacency, world-space field query results, and sample velocities.
//
// =============================================================================

#ifndef CH_FIELD_CONTACT_PRIMITIVES_H
#define CH_FIELD_CONTACT_PRIMITIVES_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "chrono/core/ChVector3.h"

namespace chrono {
namespace fieldcontact {

inline double Clamp01(double value) {
    return std::max(0.0, std::min(1.0, value));
}

inline double ClampSigned(double value) {
    return std::max(-1.0, std::min(1.0, value));
}

inline ChVector3d SafeNormalize(const ChVector3d& v, const ChVector3d& fallback) {
    double len = v.Length();
    return len > 1.0e-14 ? v / len : fallback;
}

inline ChVector3d ProjectToTangent(const ChVector3d& v, const ChVector3d& n) {
    return v - v.Dot(n) * n;
}

inline void BuildTangentBasis(const ChVector3d& normal, ChVector3d& t1, ChVector3d& t2) {
    ChVector3d n = SafeNormalize(normal, ChVector3d(0, 1, 0));
    ChVector3d ref = std::abs(n.x()) < 0.9 ? ChVector3d(1, 0, 0) : ChVector3d(0, 1, 0);
    t1 = SafeNormalize(n.Cross(ref), ChVector3d(0, 0, 1));
    t2 = SafeNormalize(n.Cross(t1), ChVector3d(1, 0, 0));
}

inline ChVector3d TransportElasticStateMinimalRotation(const ChVector3d& xi_elastic_world_prev,
                                                       const ChVector3d& normal_prev,
                                                       const ChVector3d& normal_cur) {
    ChVector3d a = SafeNormalize(normal_prev, ChVector3d(0, 1, 0));
    ChVector3d b = SafeNormalize(normal_cur, a);
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
        xi_rotated = ProjectToTangent(xi_prev_tangent, b);
    }

    return ProjectToTangent(xi_rotated, b);
}

enum class StickSlipState {
    Stick,
    Slip
};

enum class PrimitiveTopologyEvent {
    Stable,
    Newborn,
    Split,
    Merge,
    Death
};

struct SurfaceSample {
    int id = -1;
    ChVector3d local_pos = ChVector3d(0, 0, 0);
    double area = 0.0;
    std::vector<int> neighbors;
};

struct SurfaceGraph {
    std::vector<SurfaceSample> samples;

    std::vector<std::vector<int>> FindConnectedComponents(const std::vector<int>& active_indices) const {
        std::vector<std::vector<int>> components;
        if (samples.empty() || active_indices.empty()) {
            return components;
        }

        std::vector<char> active(samples.size(), 0);
        for (int idx : active_indices) {
            if (idx >= 0 && idx < static_cast<int>(samples.size())) {
                active[idx] = 1;
            }
        }

        std::vector<char> visited(samples.size(), 0);
        for (int seed : active_indices) {
            if (seed < 0 || seed >= static_cast<int>(samples.size()) || visited[seed] || !active[seed]) {
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

                for (int neighbor : samples[current].neighbors) {
                    if (neighbor < 0 || neighbor >= static_cast<int>(samples.size())) {
                        continue;
                    }
                    if (active[neighbor] && !visited[neighbor]) {
                        visited[neighbor] = 1;
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

struct TriangleFace {
    int v0 = -1;
    int v1 = -1;
    int v2 = -1;
};

inline void AddUniqueNeighbor(std::vector<int>& neighbors, int value) {
    if (std::find(neighbors.begin(), neighbors.end(), value) == neighbors.end()) {
        neighbors.push_back(value);
    }
}

inline SurfaceGraph MakeTriangleMeshSurfaceGraph(const std::vector<ChVector3d>& vertices,
                                                 const std::vector<TriangleFace>& faces,
                                                 double min_triangle_area = 1.0e-16) {
    SurfaceGraph graph;
    graph.samples.resize(vertices.size());

    for (size_t i = 0; i < vertices.size(); i++) {
        graph.samples[i].id = static_cast<int>(i);
        graph.samples[i].local_pos = vertices[i];
        graph.samples[i].area = 0.0;
    }

    auto valid_index = [&](int idx) {
        return idx >= 0 && idx < static_cast<int>(vertices.size());
    };

    auto add_edge = [&](int a, int b) {
        if (a == b || !valid_index(a) || !valid_index(b)) {
            return;
        }
        AddUniqueNeighbor(graph.samples[a].neighbors, b);
        AddUniqueNeighbor(graph.samples[b].neighbors, a);
    };

    for (const auto& face : faces) {
        if (!valid_index(face.v0) || !valid_index(face.v1) || !valid_index(face.v2)) {
            continue;
        }

        ChVector3d a = vertices[face.v0];
        ChVector3d b = vertices[face.v1];
        ChVector3d c = vertices[face.v2];
        double area = 0.5 * (b - a).Cross(c - a).Length();
        if (area <= min_triangle_area) {
            continue;
        }

        double vertex_area = area / 3.0;
        graph.samples[face.v0].area += vertex_area;
        graph.samples[face.v1].area += vertex_area;
        graph.samples[face.v2].area += vertex_area;

        add_edge(face.v0, face.v1);
        add_edge(face.v1, face.v2);
        add_edge(face.v2, face.v0);
    }

    for (auto& sample : graph.samples) {
        std::sort(sample.neighbors.begin(), sample.neighbors.end());
        sample.neighbors.erase(std::unique(sample.neighbors.begin(), sample.neighbors.end()),
                               sample.neighbors.end());
    }

    return graph;
}

inline SurfaceGraph MakeSphereSurfaceGraph(double radius, int n_theta, int n_phi) {
    SurfaceGraph graph;
    if (radius <= 0.0 || n_theta <= 0 || n_phi <= 2) {
        return graph;
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    graph.samples.reserve(static_cast<size_t>(n_theta * n_phi));

    for (int i = 0; i < n_theta; i++) {
        double theta = pi * (i + 0.5) / n_theta;
        double theta_min = pi * i / n_theta;
        double theta_max = pi * (i + 1) / n_theta;
        double ring_area = 2.0 * pi * radius * radius * (std::cos(theta_min) - std::cos(theta_max));
        double sample_area = ring_area / n_phi;

        for (int j = 0; j < n_phi; j++) {
            double phi = 2.0 * pi * j / n_phi;
            int idx = i * n_phi + j;

            SurfaceSample sample;
            sample.id = idx;
            sample.local_pos = ChVector3d(radius * std::sin(theta) * std::cos(phi),
                                          radius * std::cos(theta),
                                          radius * std::sin(theta) * std::sin(phi));
            sample.area = sample_area;

            if (i > 0) {
                sample.neighbors.push_back((i - 1) * n_phi + j);
            }
            if (i < n_theta - 1) {
                sample.neighbors.push_back((i + 1) * n_phi + j);
            }
            sample.neighbors.push_back(i * n_phi + ((j - 1 + n_phi) % n_phi));
            sample.neighbors.push_back(i * n_phi + ((j + 1) % n_phi));

            graph.samples.push_back(sample);
        }
    }

    return graph;
}

struct FieldSampleQuery {
    double phi = 0.0;
    ChVector3d grad = ChVector3d(0, 1, 0);
    ChVector3d world_pos = ChVector3d(0, 0, 0);
    ChVector3d world_vel = ChVector3d(0, 0, 0);
};

struct PrimitivePatch {
    int primitive_id = -1;
    std::vector<int> sample_ids;
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d representative_velocity = ChVector3d(0, 0, 0);
    ChVector3d normal = ChVector3d(0, 1, 0);
    ChVector3d tangent_t1 = ChVector3d(1, 0, 0);
    ChVector3d tangent_t2 = ChVector3d(0, 0, 1);
    double area = 0.0;
    double mean_phi = 0.0;
    double min_phi = 0.0;
    double max_penetration = 0.0;
    double mean_penetration = 0.0;
    ChVector3d normal_force = ChVector3d(0, 0, 0);
    ChVector3d tangential_force = ChVector3d(0, 0, 0);
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
};

struct PatchExtractionSettings {
    double activation_band = 0.0;
    double min_area = 0.0;
    int min_samples = 1;
    bool use_penetration_weighted_center = true;
};

inline std::vector<int> BuildActiveSet(const std::vector<FieldSampleQuery>& queries,
                                       double activation_band) {
    std::vector<int> active;
    for (size_t i = 0; i < queries.size(); i++) {
        if (queries[i].phi < activation_band) {
            active.push_back(static_cast<int>(i));
        }
    }
    return active;
}

inline std::vector<PrimitivePatch> ExtractPrimitives(const SurfaceGraph& graph,
                                                     const std::vector<FieldSampleQuery>& queries,
                                                     const std::vector<int>& active_indices,
                                                     const PatchExtractionSettings& settings) {
    std::vector<PrimitivePatch> primitives;
    if (graph.samples.empty() || queries.size() != graph.samples.size()) {
        return primitives;
    }

    auto components = graph.FindConnectedComponents(active_indices);
    int primitive_id = 0;

    for (const auto& component : components) {
        if (static_cast<int>(component.size()) < settings.min_samples) {
            continue;
        }

        double area = 0.0;
        double area_weighted_phi = 0.0;
        double area_weighted_penetration = 0.0;
        double center_weight = 0.0;
        ChVector3d area_center(0, 0, 0);
        ChVector3d penetration_center(0, 0, 0);
        ChVector3d velocity_sum(0, 0, 0);
        ChVector3d normal_sum(0, 0, 0);
        double min_phi = std::numeric_limits<double>::max();
        double max_penetration = 0.0;

        for (int sid : component) {
            if (sid < 0 || sid >= static_cast<int>(graph.samples.size())) {
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

            if (settings.use_penetration_weighted_center) {
                double w = sample_area * penetration;
                penetration_center += query.world_pos * w;
                center_weight += w;
            }
        }

        if (area < settings.min_area || area <= 1.0e-16) {
            continue;
        }

        PrimitivePatch primitive;
        primitive.primitive_id = primitive_id++;
        primitive.sample_ids = component;
        primitive.area = area;
        primitive.center = center_weight > 1.0e-16 ? penetration_center / center_weight : area_center / area;
        primitive.representative_velocity = velocity_sum / area;
        primitive.normal = SafeNormalize(normal_sum, ChVector3d(0, 1, 0));
        BuildTangentBasis(primitive.normal, primitive.tangent_t1, primitive.tangent_t2);
        primitive.mean_phi = area_weighted_phi / area;
        primitive.min_phi = min_phi;
        primitive.max_penetration = max_penetration;
        primitive.mean_penetration = area_weighted_penetration / area;

        primitives.push_back(primitive);
    }

    return primitives;
}

struct NormalContactSettings {
    // This stiffness is a pressure stiffness because sample force is integrated
    // as (stiffness * penetration) * sample_area.
    double stiffness = 1.0e7;
    double damping = 1.0e4;
};

struct NormalContactResult {
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
    double active_area = 0.0;
    double mean_penetration = 0.0;
    double max_penetration = 0.0;
};

inline NormalContactResult ComputeNormalContactIntegral(const PrimitivePatch& primitive,
                                                        const SurfaceGraph& graph,
                                                        const std::vector<FieldSampleQuery>& queries,
                                                        const ChVector3d& torque_reference,
                                                        const NormalContactSettings& settings) {
    NormalContactResult result;

    double penetration_area_sum = 0.0;
    for (int sid : primitive.sample_ids) {
        if (sid < 0 || sid >= static_cast<int>(graph.samples.size()) ||
            sid >= static_cast<int>(queries.size())) {
            continue;
        }

        const auto& sample = graph.samples[sid];
        const auto& query = queries[sid];
        double penetration = std::max(-query.phi, 0.0);
        if (penetration <= 0.0 || sample.area <= 0.0) {
            continue;
        }

        ChVector3d normal = SafeNormalize(query.grad, primitive.normal);
        double vn = query.world_vel.Dot(normal);
        double pressure = settings.stiffness * penetration + settings.damping * std::max(-vn, 0.0);
        pressure = std::max(0.0, pressure);

        ChVector3d dF = normal * (pressure * sample.area);
        result.force += dF;
        result.torque += (query.world_pos - torque_reference).Cross(dF);
        result.active_area += sample.area;
        penetration_area_sum += penetration * sample.area;
        result.max_penetration = std::max(result.max_penetration, penetration);
    }

    if (result.active_area > 1.0e-16) {
        result.mean_penetration = penetration_area_sum / result.active_area;
    }

    return result;
}

inline void ApplyNormalContactIntegral(PrimitivePatch& primitive,
                                       const SurfaceGraph& graph,
                                       const std::vector<FieldSampleQuery>& queries,
                                       const ChVector3d& torque_reference,
                                       const NormalContactSettings& settings) {
    NormalContactResult response =
        ComputeNormalContactIntegral(primitive, graph, queries, torque_reference, settings);
    primitive.normal_force = response.force;
    primitive.force = primitive.normal_force + primitive.tangential_force;
    primitive.torque = response.torque + (primitive.center - torque_reference).Cross(primitive.tangential_force);
}

struct PrimitiveSnapshot {
    int persistent_id = -1;
    std::vector<int> sample_ids;
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d normal = ChVector3d(0, 1, 0);
    double area = 0.0;
};

inline PrimitiveSnapshot MakeSnapshot(const PrimitivePatch& primitive, int persistent_id) {
    PrimitiveSnapshot snapshot;
    snapshot.persistent_id = persistent_id;
    snapshot.sample_ids = primitive.sample_ids;
    snapshot.center = primitive.center;
    snapshot.normal = primitive.normal;
    snapshot.area = primitive.area;
    return snapshot;
}

inline double ComputeSampleArea(const SurfaceGraph& graph, const std::vector<int>& sample_ids) {
    double area = 0.0;
    for (int sid : sample_ids) {
        if (sid >= 0 && sid < static_cast<int>(graph.samples.size())) {
            area += std::max(0.0, graph.samples[sid].area);
        }
    }
    return area;
}

inline double ComputeAreaIntersection(const SurfaceGraph& graph,
                                      const std::vector<int>& a,
                                      const std::vector<int>& b) {
    std::unordered_set<int> bset;
    bset.reserve(b.size());
    for (int sid : b) {
        bset.insert(sid);
    }

    double intersection = 0.0;
    for (int sid : a) {
        if (bset.count(sid) && sid >= 0 && sid < static_cast<int>(graph.samples.size())) {
            intersection += std::max(0.0, graph.samples[sid].area);
        }
    }
    return intersection;
}

inline double ComputeAreaJaccard(const SurfaceGraph& graph,
                                 const std::vector<int>& a,
                                 const std::vector<int>& b) {
    double area_a = ComputeSampleArea(graph, a);
    double area_b = ComputeSampleArea(graph, b);
    double intersection = ComputeAreaIntersection(graph, a, b);
    double uni = area_a + area_b - intersection;
    return uni > 1.0e-16 ? intersection / uni : 0.0;
}

struct HistorySource {
    int persistent_id = -1;
    int previous_index = -1;
    double overlap = 0.0;
    double normal_dot = 0.0;
    double center_distance = 0.0;
    double weight = 0.0;
};

struct HistoryInheritanceSettings {
    double min_overlap = 0.02;
    double min_normal_dot = 0.5;
    double max_center_distance = 0.25;
    double geometry_fallback_weight = 0.25;
};

inline std::vector<HistorySource> ComputeHistorySources(const PrimitivePatch& current,
                                                        const std::vector<PrimitiveSnapshot>& previous,
                                                        const SurfaceGraph& graph,
                                                        const HistoryInheritanceSettings& settings) {
    std::vector<HistorySource> sources;
    double raw_sum = 0.0;

    for (size_t i = 0; i < previous.size(); i++) {
        const auto& prev = previous[i];
        double overlap = ComputeAreaJaccard(graph, current.sample_ids, prev.sample_ids);
        double normal_dot = current.normal.Dot(prev.normal);
        double distance = (current.center - prev.center).Length();

        if (normal_dot < settings.min_normal_dot) {
            continue;
        }

        double raw_weight = 0.0;
        if (overlap >= settings.min_overlap) {
            raw_weight = overlap;
        } else if (distance <= settings.max_center_distance) {
            double proximity = 1.0 - distance / std::max(settings.max_center_distance, 1.0e-12);
            raw_weight = settings.geometry_fallback_weight * std::max(0.0, proximity) *
                         std::max(0.0, normal_dot);
        }

        if (raw_weight <= 0.0) {
            continue;
        }

        HistorySource source;
        source.persistent_id = prev.persistent_id;
        source.previous_index = static_cast<int>(i);
        source.overlap = overlap;
        source.normal_dot = normal_dot;
        source.center_distance = distance;
        source.weight = raw_weight;
        sources.push_back(source);
        raw_sum += raw_weight;
    }

    if (raw_sum > 1.0e-16) {
        double scale = std::min(1.0, raw_sum) / raw_sum;
        for (auto& source : sources) {
            source.weight *= scale;
        }
    }

    std::sort(sources.begin(), sources.end(), [](const HistorySource& a, const HistorySource& b) {
        return a.weight > b.weight;
    });

    return sources;
}

struct TangentialHistory {
    int persistent_id = -1;
    bool valid = false;
    ChVector3d xi_elastic_world = ChVector3d(0, 0, 0);
    ChVector3d normal = ChVector3d(0, 1, 0);
};

struct WeightedTangentialHistorySource {
    TangentialHistory history;
    double weight = 0.0;
};

inline TangentialHistory AggregateTangentialHistorySources(
    const std::vector<WeightedTangentialHistorySource>& sources,
    const ChVector3d& current_normal,
    int persistent_id) {
    TangentialHistory result;
    result.persistent_id = persistent_id;
    result.normal = SafeNormalize(current_normal, ChVector3d(0, 1, 0));
    result.xi_elastic_world = ChVector3d(0, 0, 0);
    result.valid = false;

    double weight_sum = 0.0;
    for (const auto& source : sources) {
        double w = Clamp01(source.weight);
        if (!source.history.valid || w <= 0.0) {
            continue;
        }

        ChVector3d transported = TransportElasticStateMinimalRotation(source.history.xi_elastic_world,
                                                                      source.history.normal,
                                                                      result.normal);
        result.xi_elastic_world += transported * w;
        weight_sum += w;
    }

    result.xi_elastic_world = ProjectToTangent(result.xi_elastic_world, result.normal);
    result.valid = weight_sum > 1.0e-16;
    return result;
}

struct TangentialContactSettings {
    double stiffness = 1.0e5;
    double damping = 0.0;
    double friction_coefficient = 0.5;
    double time_step = 1.0e-3;
};

struct TangentialUpdateResult {
    TangentialHistory history;
    ChVector3d transported_xi = ChVector3d(0, 0, 0);
    ChVector3d gated_xi = ChVector3d(0, 0, 0);
    ChVector3d trial_xi = ChVector3d(0, 0, 0);
    ChVector3d force = ChVector3d(0, 0, 0);
    double friction_limit = 0.0;
    double elastic_trial_force_norm = 0.0;
    double final_force_norm = 0.0;
    double gate = 0.0;
    double stored_energy_before_gate = 0.0;
    double stored_energy_after_gate = 0.0;
    StickSlipState state = StickSlipState::Stick;
};

inline TangentialUpdateResult UpdateTangentialContact(const TangentialHistory* previous_history,
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

    if (previous_history && previous_history->valid) {
        result.transported_xi = TransportElasticStateMinimalRotation(previous_history->xi_elastic_world,
                                                                     previous_history->normal,
                                                                     normal);
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

}  // namespace fieldcontact
}  // namespace chrono

#endif
