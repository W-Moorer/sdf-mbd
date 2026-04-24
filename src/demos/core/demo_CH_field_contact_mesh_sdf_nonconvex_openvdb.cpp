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
// Milestone 24: non-convex mesh -> OpenVDB SDF field-contact benchmarks.
//
// This executable moves the Milestone 21/22 deterministic primitive benchmark
// away from analytic height fields.  Each static target is an explicit triangle
// mesh, converted to an OpenVDB level-set SDF, and then queried by the existing
// field-contact primitive runtime.
//
// Outputs:
//   out/milestone_24/mesh_sdf_nonconvex_frames.csv
//   out/milestone_24/mesh_sdf_nonconvex_patches.csv
//   out/milestone_24/mesh_sdf_nonconvex_summary.csv
//   out/milestone_24/mesh_sdf_nonconvex_comparison.csv
//   out/milestone_24/mesh_sdf_nonconvex_geometry.csv
//
// =============================================================================

#define _USE_MATH_DEFINES

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
#include <tuple>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct DemoTriangleMesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct OpenVDBSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 0.006;
};

enum class MeshSdfScenario {
    UChannel,
    ToothedRail,
    KeySlotInsertion,
    StaggeredPadTrack
};

enum class BenchmarkVariant {
    FieldPrimitive,
    RawSampleContacts,
    NormalBinFragments,
    ConvexDecompositionProxy,
    ChronoTraditionalProxy
};

struct ScenarioConfig {
    MeshSdfScenario scenario;
    std::string name;
    ChVector3d start = ChVector3d(0, 0, 0);
    ChVector3d end = ChVector3d(0, 0, 0);
    double probe_radius = 0.12;
    bool use_key_probe = false;
    double z_oscillation = 0.0;
    double z_oscillation_cycles = 1.0;
    double total_time = 1.2;
    int frames = 501;
    double rolling_scale = 1.0;
    double voxel_size = 0.006;
};

struct SceneData {
    DemoTriangleMesh target_mesh;
    SurfaceGraph moving_graph;
    OpenVDBSDF sdf;
};

struct BodyKinematics {
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
    double roll_angle_x = 0.0;
    double roll_angle_z = 0.0;
    ChVector3d angular_velocity = ChVector3d(0, 0, 0);
};

struct VariantState {
    BenchmarkVariant variant = BenchmarkVariant::FieldPrimitive;
    FieldContactPrimitiveTracker tracker;
    FieldContactTopologyMetricsAccumulator metrics;
};

struct ScenarioVariantSummary {
    std::string scenario;
    std::string variant;
    FieldContactTopologyMetricsSummary metrics;
};

static const char* VariantName(BenchmarkVariant variant) {
    switch (variant) {
        case BenchmarkVariant::FieldPrimitive:
            return "field_primitive";
        case BenchmarkVariant::RawSampleContacts:
            return "raw_sample_contacts";
        case BenchmarkVariant::NormalBinFragments:
            return "normal_bin_fragments";
        case BenchmarkVariant::ConvexDecompositionProxy:
            return "convex_decomposition_proxy";
        case BenchmarkVariant::ChronoTraditionalProxy:
            return "chrono_traditional_proxy";
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

static int AddDedupVertex(DemoTriangleMesh& mesh,
                          std::map<std::tuple<long long, long long, long long>, int>& index,
                          const ChVector3d& p) {
    const double scale = 1.0e9;
    auto key = std::make_tuple(static_cast<long long>(std::llround(p.x() * scale)),
                               static_cast<long long>(std::llround(p.y() * scale)),
                               static_cast<long long>(std::llround(p.z() * scale)));
    auto it = index.find(key);
    if (it != index.end()) {
        return it->second;
    }
    int id = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(p);
    index[key] = id;
    return id;
}

static void AddGridQuad(DemoTriangleMesh& mesh,
                        std::map<std::tuple<long long, long long, long long>, int>& index,
                        const ChVector3d& p00,
                        const ChVector3d& p10,
                        const ChVector3d& p01,
                        int nu,
                        int nv) {
    nu = std::max(1, nu);
    nv = std::max(1, nv);
    std::vector<int> ids(static_cast<size_t>((nu + 1) * (nv + 1)), -1);

    auto id = [nu](int iu, int iv) {
        return iv * (nu + 1) + iu;
    };

    for (int iv = 0; iv <= nv; iv++) {
        double v = static_cast<double>(iv) / static_cast<double>(nv);
        for (int iu = 0; iu <= nu; iu++) {
            double u = static_cast<double>(iu) / static_cast<double>(nu);
            ChVector3d p = p00 + (p10 - p00) * u + (p01 - p00) * v;
            ids[id(iu, iv)] = AddDedupVertex(mesh, index, p);
        }
    }

    for (int iv = 0; iv < nv; iv++) {
        for (int iu = 0; iu < nu; iu++) {
            int v00 = ids[id(iu, iv)];
            int v10 = ids[id(iu + 1, iv)];
            int v01 = ids[id(iu, iv + 1)];
            int v11 = ids[id(iu + 1, iv + 1)];
            mesh.faces.push_back({v00, v10, v11});
            mesh.faces.push_back({v00, v11, v01});
        }
    }
}

static void AddBoxSurface(DemoTriangleMesh& mesh,
                          const ChVector3d& center,
                          const ChVector3d& half,
                          int nx,
                          int ny,
                          int nz) {
    std::map<std::tuple<long long, long long, long long>, int> index;
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
        const ChVector3d& p = mesh.vertices[i];
        const double scale = 1.0e9;
        index[std::make_tuple(static_cast<long long>(std::llround(p.x() * scale)),
                              static_cast<long long>(std::llround(p.y() * scale)),
                              static_cast<long long>(std::llround(p.z() * scale)))] = static_cast<int>(i);
    }

    const double x0 = center.x() - half.x();
    const double x1 = center.x() + half.x();
    const double y0 = center.y() - half.y();
    const double y1 = center.y() + half.y();
    const double z0 = center.z() - half.z();
    const double z1 = center.z() + half.z();

    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x0, y1, z0), ChVector3d(x1, y0, z0), ny, nx);  // -Z
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z1), ChVector3d(x1, y0, z1), ChVector3d(x0, y1, z1), nx, ny);  // +Z
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x0, y0, z1), ChVector3d(x0, y1, z0), nz, ny);  // -X
    AddGridQuad(mesh, index, ChVector3d(x1, y0, z0), ChVector3d(x1, y1, z0), ChVector3d(x1, y0, z1), ny, nz);  // +X
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x1, y0, z0), ChVector3d(x0, y0, z1), nx, nz);  // -Y
    AddGridQuad(mesh, index, ChVector3d(x0, y1, z0), ChVector3d(x0, y1, z1), ChVector3d(x1, y1, z0), nz, nx);  // +Y
}

static void AppendMesh(DemoTriangleMesh& dst, const DemoTriangleMesh& src) {
    int vertex_offset = static_cast<int>(dst.vertices.size());
    dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
    for (const auto& face : src.faces) {
        dst.faces.push_back({face.v0 + vertex_offset, face.v1 + vertex_offset, face.v2 + vertex_offset});
    }
}

static void AddTriangularPrismToothedRidge(DemoTriangleMesh& mesh,
                                           const ChVector3d& center,
                                           double base_width_x,
                                           double height_y,
                                           double half_z) {
    int base = static_cast<int>(mesh.vertices.size());
    double x0 = center.x() - 0.5 * base_width_x;
    double x1 = center.x() + 0.5 * base_width_x;
    double xp = center.x();
    double y0 = center.y();
    double yp = center.y() + height_y;
    double z0 = center.z() - half_z;
    double z1 = center.z() + half_z;

    mesh.vertices.push_back(ChVector3d(x0, y0, z0));  // 0
    mesh.vertices.push_back(ChVector3d(x1, y0, z0));  // 1
    mesh.vertices.push_back(ChVector3d(xp, yp, z0));  // 2
    mesh.vertices.push_back(ChVector3d(x0, y0, z1));  // 3
    mesh.vertices.push_back(ChVector3d(x1, y0, z1));  // 4
    mesh.vertices.push_back(ChVector3d(xp, yp, z1));  // 5

    auto f = [&](int a, int b, int c) {
        mesh.faces.push_back({base + a, base + b, base + c});
    };

    f(0, 2, 1);  // -Z
    f(3, 4, 5);  // +Z
    f(0, 1, 4);
    f(0, 4, 3);  // base
    f(0, 3, 5);
    f(0, 5, 2);  // left slope
    f(1, 2, 5);
    f(1, 5, 4);  // right slope
}

static DemoTriangleMesh BuildUChannelTarget() {
    DemoTriangleMesh mesh;
    AddBoxSurface(mesh, ChVector3d(0.0, -0.025, 0.0), ChVector3d(0.66, 0.025, 0.26), 12, 1, 5);
    AddBoxSurface(mesh, ChVector3d(0.0, 0.055, -0.170), ChVector3d(0.66, 0.055, 0.030), 12, 3, 1);
    AddBoxSurface(mesh, ChVector3d(0.0, 0.055, 0.170), ChVector3d(0.66, 0.055, 0.030), 12, 3, 1);
    AddBoxSurface(mesh, ChVector3d(-0.36, 0.025, 0.0), ChVector3d(0.045, 0.025, 0.125), 1, 1, 3);
    AddBoxSurface(mesh, ChVector3d(0.36, 0.025, 0.0), ChVector3d(0.045, 0.025, 0.125), 1, 1, 3);
    return mesh;
}

static DemoTriangleMesh BuildToothedRailTarget() {
    DemoTriangleMesh mesh;
    AddBoxSurface(mesh, ChVector3d(0.0, -0.018, 0.0), ChVector3d(0.64, 0.018, 0.18), 14, 1, 4);
    for (int i = -5; i <= 5; i++) {
        double x = 0.105 * static_cast<double>(i);
        AddTriangularPrismToothedRidge(mesh, ChVector3d(x, 0.0, 0.0), 0.070, 0.048, 0.160);
    }
    AddBoxSurface(mesh, ChVector3d(0.0, 0.018, -0.145), ChVector3d(0.62, 0.018, 0.012), 12, 1, 1);
    AddBoxSurface(mesh, ChVector3d(0.0, 0.018, 0.145), ChVector3d(0.62, 0.018, 0.012), 12, 1, 1);
    return mesh;
}

static DemoTriangleMesh BuildKeySlotTarget() {
    DemoTriangleMesh mesh;
    AddBoxSurface(mesh, ChVector3d(0.0, 0.060, -0.125), ChVector3d(0.45, 0.060, 0.028), 12, 3, 1);
    AddBoxSurface(mesh, ChVector3d(0.0, 0.060, 0.125), ChVector3d(0.45, 0.060, 0.028), 12, 3, 1);
    AddBoxSurface(mesh, ChVector3d(0.0, -0.018, 0.0), ChVector3d(0.47, 0.018, 0.165), 12, 1, 4);
    AddBoxSurface(mesh, ChVector3d(-0.24, 0.125, -0.073), ChVector3d(0.065, 0.022, 0.018), 2, 1, 1);
    AddBoxSurface(mesh, ChVector3d(-0.24, 0.125, 0.073), ChVector3d(0.065, 0.022, 0.018), 2, 1, 1);
    AddBoxSurface(mesh, ChVector3d(0.24, 0.125, -0.073), ChVector3d(0.065, 0.022, 0.018), 2, 1, 1);
    AddBoxSurface(mesh, ChVector3d(0.24, 0.125, 0.073), ChVector3d(0.065, 0.022, 0.018), 2, 1, 1);
    return mesh;
}

static DemoTriangleMesh BuildStaggeredPadTarget() {
    DemoTriangleMesh mesh;
    for (int i = -5; i <= 5; i++) {
        double x = 0.105 * static_cast<double>(i);
        double z = (i % 2 == 0) ? 0.060 : -0.060;
        AddBoxSurface(mesh, ChVector3d(x, 0.014, z), ChVector3d(0.034, 0.014, 0.040), 2, 1, 2);
        AddBoxSurface(mesh, ChVector3d(x + 0.050, 0.012, -0.62 * z), ChVector3d(0.030, 0.012, 0.036), 2, 1, 2);
    }
    AddBoxSurface(mesh, ChVector3d(0.0, -0.014, 0.0), ChVector3d(0.64, 0.014, 0.17), 14, 1, 4);
    return mesh;
}

static DemoTriangleMesh BuildKeyProbeMesh() {
    DemoTriangleMesh mesh;
    DemoTriangleMesh stem;
    DemoTriangleMesh head;
    AddBoxSurface(stem, ChVector3d(0.0, -0.010, 0.0), ChVector3d(0.100, 0.028, 0.054), 6, 2, 4);
    AddBoxSurface(head, ChVector3d(0.0, 0.034, 0.0), ChVector3d(0.100, 0.020, 0.103), 6, 2, 8);
    AppendMesh(mesh, stem);
    AppendMesh(mesh, head);
    return mesh;
}

static openvdb::FloatGrid::Ptr BuildLevelSetFromTriangleMesh(const DemoTriangleMesh& mesh,
                                                             double voxel_size,
                                                             float half_width_voxels) {
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec3I> triangles;
    points.reserve(mesh.vertices.size());
    triangles.reserve(mesh.faces.size());

    for (const auto& v : mesh.vertices) {
        points.emplace_back(static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z()));
    }
    for (const auto& face : mesh.faces) {
        triangles.emplace_back(face.v0, face.v1, face.v2);
    }

    auto transform = openvdb::math::Transform::createLinearTransform(voxel_size);
    auto grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*transform, points, triangles, half_width_voxels);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
    return grid;
}

static SceneData BuildScene(const ScenarioConfig& config) {
    SceneData scene;
    scene.sdf.voxel_size = config.voxel_size;

    switch (config.scenario) {
        case MeshSdfScenario::UChannel:
            scene.target_mesh = BuildUChannelTarget();
            break;
        case MeshSdfScenario::ToothedRail:
            scene.target_mesh = BuildToothedRailTarget();
            break;
        case MeshSdfScenario::KeySlotInsertion:
            scene.target_mesh = BuildKeySlotTarget();
            break;
        case MeshSdfScenario::StaggeredPadTrack:
            scene.target_mesh = BuildStaggeredPadTarget();
            break;
    }

    if (config.use_key_probe) {
        DemoTriangleMesh key = BuildKeyProbeMesh();
        scene.moving_graph = MakeTriangleMeshSurfaceGraph(key.vertices, key.faces);
    } else {
        scene.moving_graph = MakeSphereSurfaceGraph(config.probe_radius, 20, 40);
    }

    scene.sdf.grid = BuildLevelSetFromTriangleMesh(scene.target_mesh, config.voxel_size, 7.0f);
    return scene;
}

static FieldContactRuntimeSettings MakeRuntimeSettings(double time_step, double activation_band) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = activation_band;
    settings.extraction.min_area = 1.0e-8;
    settings.extraction.min_samples = 3;
    settings.extraction.use_penetration_weighted_center = true;

    settings.normal.stiffness = 1.8e6;
    settings.normal.damping = 0.0;

    settings.tangential.stiffness = 3.0e3;
    settings.tangential.damping = 0.0;
    settings.tangential.friction_coefficient = 0.45;
    settings.tangential.time_step = time_step;

    settings.inheritance.min_overlap = 0.005;
    settings.inheritance.min_normal_dot = 0.16;
    settings.inheritance.max_center_distance = 0.075;
    settings.inheritance.geometry_fallback_weight = 0.18;
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

    auto center_at = [&](double uu) {
        double osc = scenario.z_oscillation * std::sin(2.0 * kPi * scenario.z_oscillation_cycles * uu);
        return scenario.start + (scenario.end - scenario.start) * uu + ChVector3d(0, 0, osc);
    };

    double du = 1.0 / static_cast<double>(scenario.frames - 1);
    double u0 = std::max(0.0, u - du);
    double u1 = std::min(1.0, u + du);

    BodyKinematics kin;
    kin.center = center_at(u);
    kin.velocity = (center_at(u1) - center_at(u0)) / ((u1 - u0) * scenario.total_time);
    if (scenario.rolling_scale > 0.0 && scenario.probe_radius > 1.0e-9) {
        kin.roll_angle_z = -scenario.rolling_scale * (kin.center.x() - scenario.start.x()) / scenario.probe_radius;
        kin.roll_angle_x = scenario.rolling_scale * (kin.center.z() - scenario.start.z()) / scenario.probe_radius;
        kin.angular_velocity = ChVector3d(kin.velocity.z() / scenario.probe_radius,
                                          0.0,
                                          -kin.velocity.x() / scenario.probe_radius) *
                               scenario.rolling_scale;
    }
    return kin;
}

using SdfSampler = openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler>;

static FieldSampleQuery QueryOpenVDB(const OpenVDBSDF& sdf,
                                     SdfSampler& sampler,
                                     const ChVector3d& world_pos,
                                     const ChVector3d& world_vel) {
    FieldSampleQuery query;
    query.world_pos = world_pos;
    query.world_vel = world_vel;

    query.phi = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y(), world_pos.z())));

    const double h = 0.5 * sdf.voxel_size;
    double phi_px = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x() + h, world_pos.y(), world_pos.z())));
    double phi_mx = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x() - h, world_pos.y(), world_pos.z())));
    double phi_py = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y() + h, world_pos.z())));
    double phi_my = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y() - h, world_pos.z())));
    double phi_pz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y(), world_pos.z() + h)));
    double phi_mz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y(), world_pos.z() - h)));

    query.grad = SafeNormalize(ChVector3d((phi_px - phi_mx) / (2.0 * h),
                                          (phi_py - phi_my) / (2.0 * h),
                                          (phi_pz - phi_mz) / (2.0 * h)),
                               ChVector3d(0, 1, 0));
    return query;
}

static std::vector<FieldSampleQuery> BuildQueries(const SurfaceGraph& graph,
                                                  const OpenVDBSDF& sdf,
                                                  const BodyKinematics& kin) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(graph.samples.size());
    SdfSampler sampler(sdf.grid->tree(), sdf.grid->transform());

    for (const auto& sample : graph.samples) {
        ChVector3d rotated = RotateLocal(sample.local_pos, kin.roll_angle_x, kin.roll_angle_z);
        ChVector3d world_pos = kin.center + rotated;
        ChVector3d world_vel = kin.velocity + kin.angular_velocity.Cross(rotated);
        queries.push_back(QueryOpenVDB(sdf, sampler, world_pos, world_vel));
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

static int ConvexProxyPieceKey(const ChVector3d& local_pos) {
    double angle = std::atan2(local_pos.z(), local_pos.x());
    int sector = static_cast<int>(std::floor((angle + kPi) / (kPi / 3.0)));
    sector = std::max(0, std::min(5, sector));
    int lower_band = local_pos.y() < 0.0 ? 0 : 1;
    return sector + 6 * lower_band;
}

static std::vector<PrimitivePatch> BuildConvexDecompositionProxyPatches(
    const SurfaceGraph& graph,
    const std::vector<FieldSampleQuery>& queries,
    const PatchExtractionSettings& settings) {
    std::vector<PrimitivePatch> patches;
    std::vector<int> active = BuildActiveSet(queries, settings.activation_band);
    std::vector<char> active_mask(graph.samples.size(), 0);
    std::vector<int> piece(graph.samples.size(), -1);
    for (int sid : active) {
        if (sid >= 0 && sid < static_cast<int>(graph.samples.size())) {
            active_mask[sid] = 1;
            piece[sid] = ConvexProxyPieceKey(graph.samples[sid].local_pos);
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
                if (!visited[neighbor] && active_mask[neighbor] && piece[neighbor] == piece[seed]) {
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

static std::vector<PrimitivePatch> BuildChronoTraditionalProxyPatches(
    const SurfaceGraph& graph,
    const std::vector<FieldSampleQuery>& queries,
    const PatchExtractionSettings& settings) {
    std::vector<int> active = BuildActiveSet(queries, settings.activation_band);
    std::sort(active.begin(), active.end(), [&](int a, int b) {
        return queries[a].phi < queries[b].phi;
    });

    std::vector<int> selected;
    selected.reserve(10);
    const int max_contacts = 10;
    const double min_spacing = 0.035;

    for (int sid : active) {
        if (sid < 0 || sid >= static_cast<int>(graph.samples.size())) {
            continue;
        }
        if (queries[sid].phi >= 0.0 && !selected.empty()) {
            continue;
        }

        bool separated = true;
        for (int previous : selected) {
            if ((graph.samples[sid].local_pos - graph.samples[previous].local_pos).Length() < min_spacing) {
                separated = false;
                break;
            }
        }
        if (!separated) {
            continue;
        }

        selected.push_back(sid);
        if (static_cast<int>(selected.size()) >= max_contacts) {
            break;
        }
    }

    std::vector<PrimitivePatch> patches;
    patches.reserve(selected.size());
    int primitive_id = 0;
    for (int sid : selected) {
        patches.push_back(BuildPatchFromSampleIds(graph, queries, std::vector<int>{sid}, primitive_id++, false));
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
    } else if (state.variant == BenchmarkVariant::ConvexDecompositionProxy) {
        patches = BuildConvexDecompositionProxyPatches(graph, queries, settings.extraction);
    } else if (state.variant == BenchmarkVariant::ChronoTraditionalProxy) {
        std::vector<FieldSampleQuery> chrono_queries = queries;
        for (auto& query : chrono_queries) {
            if (query.phi < settings.extraction.activation_band) {
                query.grad = SnapNormalToBin(query.grad);
            }
        }
        patches = BuildChronoTraditionalProxyPatches(graph, chrono_queries, settings.extraction);
        return state.tracker.EvaluatePatches(graph, chrono_queries, patches, torque_reference, settings);
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
                        const std::vector<BenchmarkVariant>& variants,
                        std::ofstream& frames_out,
                        std::ofstream& patches_out,
                        std::ofstream& geometry_out,
                        std::vector<ScenarioVariantSummary>& summaries) {
    SceneData scene = BuildScene(scenario);
    const double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    FieldContactRuntimeSettings settings = MakeRuntimeSettings(dt, 1.35 * scenario.voxel_size);

    geometry_out << scenario.name << "," << scene.target_mesh.vertices.size() << "," << scene.target_mesh.faces.size()
                 << "," << scene.moving_graph.samples.size() << "," << scene.sdf.voxel_size << ","
                 << scene.sdf.grid->activeVoxelCount() << "\n";

    std::vector<VariantState> states;
    states.reserve(variants.size());
    for (auto variant : variants) {
        VariantState state;
        state.variant = variant;
        state.tracker.Reset();
        state.metrics.Reset();
        states.push_back(state);
    }

    std::cout << "Scenario " << scenario.name << ": target vertices=" << scene.target_mesh.vertices.size()
              << " faces=" << scene.target_mesh.faces.size()
              << " active_voxels=" << scene.sdf.grid->activeVoxelCount()
              << " moving_samples=" << scene.moving_graph.samples.size() << std::endl;

    for (int frame = 0; frame < scenario.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(scenario, frame);
        std::vector<FieldSampleQuery> queries = BuildQueries(scene.moving_graph, scene.sdf, kin);

        for (auto& state : states) {
            FieldContactStepResult step = EvaluateVariant(state, scene.moving_graph, queries, kin.center, settings);
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

        std::cout << "  " << summary.variant << ": active_frames=" << summary.metrics.active_frames
                  << " mean_patch=" << summary.metrics.mean_patch_count
                  << " mean_dcount=" << summary.metrics.mean_abs_patch_count_change
                  << " rms_normal=" << summary.metrics.rms_normal_jump_angle
                  << " force_osc=" << summary.metrics.force_oscillation_index
                  << " torque_osc=" << summary.metrics.torque_oscillation_index
                  << " max_ratio=" << summary.metrics.max_tangential_force_ratio << std::endl;
    }
}

}  // namespace

int main() {
    std::cout << "Milestone 24: mesh-to-OpenVDB SDF non-convex field-contact benchmarks" << std::endl;

    openvdb::initialize();

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_24";
    std::filesystem::create_directories(out_dir);

    std::ofstream frames_out(out_dir / "mesh_sdf_nonconvex_frames.csv");
    std::ofstream patches_out(out_dir / "mesh_sdf_nonconvex_patches.csv");
    std::ofstream summary_out(out_dir / "mesh_sdf_nonconvex_summary.csv");
    std::ofstream comparison_out(out_dir / "mesh_sdf_nonconvex_comparison.csv");
    std::ofstream geometry_out(out_dir / "mesh_sdf_nonconvex_geometry.csv");

    frames_out << std::fixed << std::setprecision(8);
    patches_out << std::fixed << std::setprecision(8);
    summary_out << std::fixed << std::setprecision(8);
    comparison_out << std::fixed << std::setprecision(8);
    geometry_out << std::fixed << std::setprecision(8);

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

    geometry_out << "scenario,target_vertices,target_faces,moving_samples,voxel_size,active_voxels\n";

    std::vector<ScenarioConfig> scenarios = {
        {MeshSdfScenario::UChannel,
         "u_channel_sphere_slide",
         ChVector3d(-0.54, 0.116, 0.000),
         ChVector3d(0.54, 0.116, 0.000),
         0.120,
         false,
         0.045,
         2.0,
         1.25,
         501,
         1.0,
         0.006},
        {MeshSdfScenario::ToothedRail,
         "toothed_rail_sphere_slide",
         ChVector3d(-0.54, 0.145, -0.010),
         ChVector3d(0.54, 0.145, 0.010),
         0.110,
         false,
         0.022,
         2.0,
         1.25,
         501,
         1.0,
         0.006},
        {MeshSdfScenario::KeySlotInsertion,
         "key_slot_insertion",
         ChVector3d(-0.54, 0.070, 0.000),
         ChVector3d(0.54, 0.070, 0.012),
         0.100,
         true,
         0.018,
         1.5,
         1.30,
         521,
         0.0,
         0.0055},
        {MeshSdfScenario::StaggeredPadTrack,
         "staggered_pad_track",
         ChVector3d(-0.54, 0.145, -0.074),
         ChVector3d(0.54, 0.145, 0.074),
         0.120,
         false,
         0.036,
         2.5,
         1.35,
         541,
         1.0,
         0.006},
    };

    std::vector<BenchmarkVariant> variants = {
        BenchmarkVariant::FieldPrimitive,
        BenchmarkVariant::RawSampleContacts,
        BenchmarkVariant::NormalBinFragments,
        BenchmarkVariant::ConvexDecompositionProxy,
        BenchmarkVariant::ChronoTraditionalProxy,
    };

    std::vector<ScenarioVariantSummary> summaries;
    for (const auto& scenario : scenarios) {
        RunScenario(scenario, variants, frames_out, patches_out, geometry_out, summaries);
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
    std::cout << "  " << (out_dir / "mesh_sdf_nonconvex_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "mesh_sdf_nonconvex_patches.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "mesh_sdf_nonconvex_summary.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "mesh_sdf_nonconvex_comparison.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "mesh_sdf_nonconvex_geometry.csv").string() << std::endl;
    return 0;
}
