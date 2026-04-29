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

struct ChSdfContactSampleQuery {
    double phi = 0.0;
    ChVector3d grad = ChVector3d(0, 1, 0);
    ChVector3d world_pos = ChVector3d(0, 0, 0);
    ChVector3d world_vel = ChVector3d(0, 0, 0);
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
