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
// Backend-neutral SDF contact geometry frontend.
//
// This header owns only geometric data exchanged between an SDF source and a
// contact backend: surface samples, sample adjacency, triangle faces, and SDF
// query values. It intentionally does not define pressure/field-contact
// constitutive laws or SDF-NCP complementarity residuals.
//
// Mesh-to-SDF frontends such as OpenVDB/NanoVDB demos can produce these types.
// Pressure-field and SDF-NCP backends can then consume the same geometry/query
// stream without depending on each other's solver or force logic.
//
// =============================================================================

#ifndef CH_SDF_CONTACT_GEOMETRY_H
#define CH_SDF_CONTACT_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "chrono/core/ChVector3.h"

namespace chrono {
namespace sdfcontact {

inline ChVector3d NormalizeOrFallback(const ChVector3d& v, const ChVector3d& fallback) {
    const double len = v.Length();
    return len > 1.0e-14 ? v / len : fallback;
}

struct ChSdfTriangleFace {
    int v0 = -1;
    int v1 = -1;
    int v2 = -1;
};

struct ChSdfContactSurfaceSample {
    int id = -1;
    ChVector3d local_pos = ChVector3d(0, 0, 0);
    double area = 0.0;
    std::vector<int> neighbors;
};

struct ChSdfContactSurfaceGraph {
    std::vector<ChSdfContactSurfaceSample> samples;

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

inline double ChSdfAxisValue(const ChVector3d& v, int axis) {
    return axis == 0 ? v.x() : (axis == 1 ? v.y() : v.z());
}

struct ChSdfAabb {
    ChVector3d min = ChVector3d(std::numeric_limits<double>::max(),
                                std::numeric_limits<double>::max(),
                                std::numeric_limits<double>::max());
    ChVector3d max = ChVector3d(-std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max());

    bool IsValid() const {
        return min.x() <= max.x() && min.y() <= max.y() && min.z() <= max.z();
    }

    void Reset() {
        min = ChVector3d(std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
        max = ChVector3d(-std::numeric_limits<double>::max(),
                         -std::numeric_limits<double>::max(),
                         -std::numeric_limits<double>::max());
    }

    void Include(const ChVector3d& p) {
        min = ChVector3d(std::min(min.x(), p.x()), std::min(min.y(), p.y()), std::min(min.z(), p.z()));
        max = ChVector3d(std::max(max.x(), p.x()), std::max(max.y(), p.y()), std::max(max.z(), p.z()));
    }

    void Include(const ChSdfAabb& other) {
        if (!other.IsValid()) {
            return;
        }
        Include(other.min);
        Include(other.max);
    }

    void Expand(double padding) {
        if (!IsValid()) {
            return;
        }
        const double p = std::max(0.0, padding);
        min -= ChVector3d(p, p, p);
        max += ChVector3d(p, p, p);
    }

    bool Overlaps(const ChSdfAabb& other) const {
        if (!IsValid() || !other.IsValid()) {
            return false;
        }
        return min.x() <= other.max.x() && max.x() >= other.min.x() && min.y() <= other.max.y() &&
               max.y() >= other.min.y() && min.z() <= other.max.z() && max.z() >= other.min.z();
    }

    ChVector3d Center() const {
        return (min + max) * 0.5;
    }

    ChVector3d Extent() const {
        return max - min;
    }

    int MaxExtentAxis() const {
        const ChVector3d e = Extent();
        if (e.x() >= e.y() && e.x() >= e.z()) {
            return 0;
        }
        return e.y() >= e.z() ? 1 : 2;
    }
};

inline ChSdfAabb MakeAabbFromPoints(const std::vector<ChVector3d>& points) {
    ChSdfAabb bounds;
    for (const auto& p : points) {
        bounds.Include(p);
    }
    return bounds;
}

inline ChSdfAabb MakeAabbFromSurfaceGraph(const ChSdfContactSurfaceGraph& graph) {
    ChSdfAabb bounds;
    for (const auto& sample : graph.samples) {
        if (sample.area > 0.0) {
            bounds.Include(sample.local_pos);
        }
    }
    return bounds;
}

class ChSdfContactSampleBvh {
  public:
    explicit ChSdfContactSampleBvh(const ChSdfContactSurfaceGraph& graph, int leaf_size = 16) {
        Build(graph, leaf_size);
    }

    ChSdfContactSampleBvh() = default;

    void Build(const ChSdfContactSurfaceGraph& graph, int leaf_size = 16) {
        m_entries.clear();
        m_nodes.clear();
        m_leaf_size = std::max(1, leaf_size);
        m_entries.reserve(graph.samples.size());
        for (const auto& sample : graph.samples) {
            if (sample.area <= 0.0) {
                continue;
            }
            m_entries.push_back(Entry{sample.id, sample.local_pos});
        }
        if (!m_entries.empty()) {
            BuildRecursive(0, static_cast<int>(m_entries.size()));
        }
    }

    bool Empty() const {
        return m_nodes.empty();
    }

    size_t Size() const {
        return m_entries.size();
    }

    void Query(const ChSdfAabb& bounds, std::vector<int>& sample_ids) const {
        if (m_nodes.empty() || !bounds.IsValid()) {
            return;
        }
        QueryNode(0, bounds, sample_ids);
    }

    std::vector<int> Query(const ChSdfAabb& bounds) const {
        std::vector<int> sample_ids;
        Query(bounds, sample_ids);
        return sample_ids;
    }

  private:
    struct Entry {
        int sample_id = -1;
        ChVector3d pos = ChVector3d(0, 0, 0);
    };

    struct Node {
        ChSdfAabb bounds;
        int begin = 0;
        int end = 0;
        int left = -1;
        int right = -1;

        bool IsLeaf() const {
            return left < 0 && right < 0;
        }
    };

    int BuildRecursive(int begin, int end) {
        Node node;
        node.begin = begin;
        node.end = end;
        for (int i = begin; i < end; i++) {
            node.bounds.Include(m_entries[static_cast<size_t>(i)].pos);
        }

        const int node_index = static_cast<int>(m_nodes.size());
        m_nodes.push_back(node);

        if (end - begin <= m_leaf_size) {
            return node_index;
        }

        const int axis = node.bounds.MaxExtentAxis();
        const int mid = begin + (end - begin) / 2;
        std::nth_element(m_entries.begin() + begin,
                         m_entries.begin() + mid,
                         m_entries.begin() + end,
                         [axis](const Entry& a, const Entry& b) {
                             return ChSdfAxisValue(a.pos, axis) < ChSdfAxisValue(b.pos, axis);
                         });

        m_nodes[static_cast<size_t>(node_index)].left = BuildRecursive(begin, mid);
        m_nodes[static_cast<size_t>(node_index)].right = BuildRecursive(mid, end);
        return node_index;
    }

    void QueryNode(int node_index, const ChSdfAabb& bounds, std::vector<int>& sample_ids) const {
        const Node& node = m_nodes[static_cast<size_t>(node_index)];
        if (!node.bounds.Overlaps(bounds)) {
            return;
        }
        if (node.IsLeaf()) {
            for (int i = node.begin; i < node.end; i++) {
                sample_ids.push_back(m_entries[static_cast<size_t>(i)].sample_id);
            }
            return;
        }
        QueryNode(node.left, bounds, sample_ids);
        QueryNode(node.right, bounds, sample_ids);
    }

    std::vector<Entry> m_entries;
    std::vector<Node> m_nodes;
    int m_leaf_size = 16;
};

struct ChSdfContactSampleQuery {
    double phi = 0.0;
    ChVector3d grad = ChVector3d(0, 1, 0);
    ChVector3d raw_grad = ChVector3d(0, 1, 0);
    double raw_grad_norm = 1.0;
    bool has_hessian = false;
    // Columns of the local SDF Hessian d(raw_grad)/d(local_pos).
    ChVector3d hessian_x = ChVector3d(0, 0, 0);
    ChVector3d hessian_y = ChVector3d(0, 0, 0);
    ChVector3d hessian_z = ChVector3d(0, 0, 0);
    ChVector3d world_pos = ChVector3d(0, 0, 0);
    ChVector3d world_vel = ChVector3d(0, 0, 0);

    ChVector3d HessianTimes(const ChVector3d& local_direction) const {
        return hessian_x * local_direction.x() + hessian_y * local_direction.y() +
               hessian_z * local_direction.z();
    }

    ChVector3d UnitGradDirectionalDerivative(const ChVector3d& local_direction) const {
        if (!has_hessian) {
            return ChVector3d(0, 0, 0);
        }

        const double grad_len = raw_grad_norm > 1.0e-14 && std::isfinite(raw_grad_norm) ? raw_grad_norm : 1.0;
        const ChVector3d unit_grad = NormalizeOrFallback(raw_grad, grad);
        const ChVector3d d_raw_grad = HessianTimes(local_direction);
        return (d_raw_grad - unit_grad * unit_grad.Dot(d_raw_grad)) / grad_len;
    }
};

inline void AddUniqueNeighbor(std::vector<int>& neighbors, int value) {
    if (std::find(neighbors.begin(), neighbors.end(), value) == neighbors.end()) {
        neighbors.push_back(value);
    }
}

inline ChSdfContactSurfaceGraph MakeTriangleMeshSurfaceGraph(const std::vector<ChVector3d>& vertices,
                                                             const std::vector<ChSdfTriangleFace>& faces,
                                                             double min_triangle_area = 1.0e-16) {
    ChSdfContactSurfaceGraph graph;
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

        const ChVector3d a = vertices[face.v0];
        const ChVector3d b = vertices[face.v1];
        const ChVector3d c = vertices[face.v2];
        const double area = 0.5 * (b - a).Cross(c - a).Length();
        if (area <= min_triangle_area) {
            continue;
        }

        const double vertex_area = area / 3.0;
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

inline ChSdfContactSurfaceGraph MakeSphereSurfaceGraph(double radius, int n_theta, int n_phi) {
    ChSdfContactSurfaceGraph graph;
    if (radius <= 0.0 || n_theta <= 0 || n_phi <= 2) {
        return graph;
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    graph.samples.reserve(static_cast<size_t>(n_theta * n_phi));

    for (int i = 0; i < n_theta; i++) {
        const double theta = pi * (i + 0.5) / n_theta;
        const double theta_min = pi * i / n_theta;
        const double theta_max = pi * (i + 1) / n_theta;
        const double ring_area = 2.0 * pi * radius * radius * (std::cos(theta_min) - std::cos(theta_max));
        const double sample_area = ring_area / n_phi;

        for (int j = 0; j < n_phi; j++) {
            const double phi = 2.0 * pi * j / n_phi;
            const int idx = i * n_phi + j;

            ChSdfContactSurfaceSample sample;
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

}  // namespace sdfcontact
}  // namespace chrono

#endif
