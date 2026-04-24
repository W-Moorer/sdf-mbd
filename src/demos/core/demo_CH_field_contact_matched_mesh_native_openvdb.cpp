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
// Milestone 27: matched native mesh baseline.
//
// Each benchmark creates one explicit triangle mesh and reuses that exact mesh
// for two backends:
//   1. Chrono native SMC collision through ChCollisionShapeTriangleMesh.
//   2. Field-contact primitive queries through an OpenVDB level-set SDF.
//
// This removes the geometry mismatch between the previous Chrono-native box
// proxy baseline and the OpenVDB field primitive experiments. The output is
// aligned per frame and includes contact/patch counts, force/torque jumps,
// displacement/velocity/acceleration tracking, Coulomb ratio, and stability
// ratios suitable for regression checks.
//
// Outputs:
//   out/milestone_27/matched_mesh_field_frames.csv
//   out/milestone_27/matched_mesh_native_frames.csv
//   out/milestone_27/matched_mesh_comparison.csv
//   out/milestone_27/matched_mesh_geometry.csv
//
// =============================================================================

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/collision/ChFieldContactRuntime.h"
#include "chrono/core/ChRotation.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChSystemSMC.h"

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

enum class MatchedScenario {
    UChannel,
    StaggeredPadTrack
};

struct ScenarioConfig {
    MatchedScenario scenario;
    std::string name;
    ChVector3d start = ChVector3d(0, 0, 0);
    ChVector3d end = ChVector3d(0, 0, 0);
    double probe_radius = 0.12;
    double z_oscillation = 0.0;
    double z_oscillation_cycles = 1.0;
    double total_time = 1.25;
    int frames = 501;
    double rolling_scale = 1.0;
    double voxel_size = 0.006;
    double native_mesh_radius = 0.0;
};

struct SceneData {
    DemoTriangleMesh target_mesh;
    SurfaceGraph moving_graph;
    OpenVDBSDF sdf;
};

struct BodyKinematics {
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
    ChVector3d acceleration = ChVector3d(0, 0, 0);
    double roll_angle_x = 0.0;
    double roll_angle_z = 0.0;
    ChVector3d angular_velocity = ChVector3d(0, 0, 0);
};

struct FieldFrameRecord {
    std::string scenario;
    int frame = 0;
    double time = 0.0;
    BodyKinematics kin;
    int patch_count = 0;
    int newborn_count = 0;
    int merge_count = 0;
    int split_count = 0;
    int death_count = 0;
    int max_source_count = 0;
    int max_previous_reuse = 0;
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
    double max_tangential_force_ratio = 0.0;
    double max_inherited_energy_ratio = 0.0;
};

struct NativeFrameRecord {
    std::string scenario;
    int frame = 0;
    double time = 0.0;
    BodyKinematics cmd;
    ChVector3d actual_pos = ChVector3d(0, 0, 0);
    ChVector3d actual_vel = ChVector3d(0, 0, 0);
    ChVector3d actual_acc = ChVector3d(0, 0, 0);
    ChVector3d pos_error = ChVector3d(0, 0, 0);
    ChVector3d vel_error = ChVector3d(0, 0, 0);
    ChVector3d acc_error = ChVector3d(0, 0, 0);
    int contact_count = 0;
    ChVector3d contact_normal = ChVector3d(0, 0, 0);
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
};

struct NativeMetricSummary {
    int frames = 0;
    int active_frames = 0;
    double mean_contact_count = 0.0;
    int max_contact_count = 0;
    double mean_abs_contact_count_change = 0.0;
    int max_abs_contact_count_change = 0;
    double mean_force_jump = 0.0;
    double rms_force_jump = 0.0;
    double max_force_jump = 0.0;
    double max_force_norm = 0.0;
    double force_oscillation_index = 0.0;
    double rms_contact_normal_jump_angle = 0.0;
    double max_contact_normal_jump_angle = 0.0;
    double mean_torque_jump = 0.0;
    double rms_torque_jump = 0.0;
    double max_torque_jump = 0.0;
    double max_torque_norm = 0.0;
    double torque_oscillation_index = 0.0;
    double rms_force_direction_jump_angle = 0.0;
    double max_force_direction_jump_angle = 0.0;
    double max_position_error = 0.0;
    double max_velocity_error = 0.0;
    double max_acceleration_error = 0.0;
    double rms_position_error = 0.0;
    double rms_velocity_error = 0.0;
    double rms_acceleration_error = 0.0;
};

struct ScenarioResult {
    ScenarioConfig config;
    FieldContactTopologyMetricsSummary field_summary;
    NativeMetricSummary native_summary;
};

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

    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x0, y1, z0), ChVector3d(x1, y0, z0), ny, nx);
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z1), ChVector3d(x1, y0, z1), ChVector3d(x0, y1, z1), nx, ny);
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x0, y0, z1), ChVector3d(x0, y1, z0), nz, ny);
    AddGridQuad(mesh, index, ChVector3d(x1, y0, z0), ChVector3d(x1, y1, z0), ChVector3d(x1, y0, z1), ny, nz);
    AddGridQuad(mesh, index, ChVector3d(x0, y0, z0), ChVector3d(x1, y0, z0), ChVector3d(x0, y0, z1), nx, nz);
    AddGridQuad(mesh, index, ChVector3d(x0, y1, z0), ChVector3d(x0, y1, z1), ChVector3d(x1, y1, z0), nz, nx);
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
    scene.target_mesh =
        config.scenario == MatchedScenario::UChannel ? BuildUChannelTarget() : BuildStaggeredPadTarget();
    scene.moving_graph = MakeSphereSurfaceGraph(config.probe_radius, 20, 40);
    scene.sdf.grid = BuildLevelSetFromTriangleMesh(scene.target_mesh, config.voxel_size, 7.0f);
    return scene;
}

static std::shared_ptr<ChTriangleMeshConnected> ToChronoTriangleMesh(const DemoTriangleMesh& mesh) {
    auto out = chrono_types::make_shared<ChTriangleMeshConnected>();
    out->GetCoordsVertices() = mesh.vertices;
    auto& indices = out->GetIndicesVertices();
    indices.reserve(mesh.faces.size());
    for (const auto& face : mesh.faces) {
        indices.push_back(ChVector3i(face.v0, face.v1, face.v2));
    }
    return out;
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
    double du = 1.0 / static_cast<double>(scenario.frames - 1);
    double dt = scenario.total_time * du;

    auto center_at = [&](double uu) {
        double clamped = std::max(0.0, std::min(1.0, uu));
        double osc = scenario.z_oscillation * std::sin(2.0 * kPi * scenario.z_oscillation_cycles * clamped);
        return scenario.start + (scenario.end - scenario.start) * clamped + ChVector3d(0, 0, osc);
    };

    double u0 = std::max(0.0, u - du);
    double u1 = std::min(1.0, u + du);

    BodyKinematics kin;
    kin.center = center_at(u);
    kin.velocity = (center_at(u1) - center_at(u0)) / ((u1 - u0) * scenario.total_time);
    if (frame > 0 && frame + 1 < scenario.frames) {
        kin.acceleration = (center_at(u + du) - 2.0 * center_at(u) + center_at(u - du)) / (dt * dt);
    }

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

static std::shared_ptr<ChContactMaterialSMC> MakeMaterial() {
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetFriction(0.45f);
    mat->SetRestitution(0.02f);
    mat->SetYoungModulus(1.2e4f);
    mat->SetPoissonRatio(0.30f);
    return mat;
}

class NativeContactNormalReporter : public ChContactContainer::ReportContactCallback {
  public:
    bool OnReportContact(const ChVector3d& pA,
                         const ChVector3d& pB,
                         const ChMatrix33<>& plane_coord,
                         double distance,
                         double eff_radius,
                         const ChVector3d& react_forces,
                         const ChVector3d& react_torques,
                         ChContactable* contactobjA,
                         ChContactable* contactobjB,
                         int constraint_offset) override {
        (void)pA;
        (void)pB;
        (void)distance;
        (void)eff_radius;
        (void)react_torques;
        (void)contactobjA;
        (void)contactobjB;
        (void)constraint_offset;

        ChVector3d normal = SafeNormalize(plane_coord.GetAxisX(), ChVector3d(0, 1, 0));
        if (react_forces.x() < 0.0) {
            normal *= -1.0;
        }
        double weight = std::max(1.0e-12, std::abs(react_forces.x()));
        weighted_normal += normal * weight;
        weight_sum += weight;
        count++;
        return true;
    }

    ChVector3d GetAverageNormal(const ChVector3d& total_force) const {
        if (weight_sum <= 1.0e-14 || count == 0) {
            return ChVector3d(0, 0, 0);
        }
        ChVector3d normal = SafeNormalize(weighted_normal / weight_sum, ChVector3d(0, 1, 0));
        if (total_force.Length() > 1.0e-12 && normal.Dot(total_force) < 0.0) {
            normal *= -1.0;
        }
        return normal;
    }

    ChVector3d weighted_normal = ChVector3d(0, 0, 0);
    double weight_sum = 0.0;
    int count = 0;
};

static std::shared_ptr<ChBody> AddMatchedMeshTarget(ChSystemSMC& sys,
                                                    const std::shared_ptr<ChContactMaterialSMC>& mat,
                                                    const DemoTriangleMesh& mesh,
                                                    double native_mesh_radius) {
    auto target = chrono_types::make_shared<ChBody>();
    target->SetFixed(true);
    target->EnableCollision(true);
    auto chrono_mesh = ToChronoTriangleMesh(mesh);
    auto shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat, chrono_mesh, true, false, native_mesh_radius);
    target->AddCollisionShape(shape);
    sys.AddBody(target);
    return target;
}

static std::shared_ptr<ChBody> AddMovingSphere(ChSystemSMC& sys,
                                               const std::shared_ptr<ChContactMaterialSMC>& mat,
                                               double radius) {
    const double kinematic_probe_density = 1.0e7;
    auto body = chrono_types::make_shared<ChBodyEasySphere>(radius, kinematic_probe_density, false, true, mat);
    body->SetFixed(false);
    body->EnableCollision(true);
    sys.AddBody(body);
    return body;
}

static void ApplyCommandedKinematics(const std::shared_ptr<ChBody>& body, const BodyKinematics& kin) {
    body->SetPos(kin.center);
    body->SetPosDt(kin.velocity);
    body->SetPosDt2(kin.acceleration);
    body->SetRot(QuatFromAngleX(kin.roll_angle_x) * QuatFromAngleZ(kin.roll_angle_z));
    body->SetAngVelLocal(kin.angular_velocity);
}

static double SafeRatio(double numerator, double denominator) {
    return denominator > 1.0e-14 ? numerator / denominator : 0.0;
}

static double DirectionJumpAngle(const ChVector3d& a, const ChVector3d& b) {
    double la = a.Length();
    double lb = b.Length();
    if (la <= 1.0e-12 || lb <= 1.0e-12) {
        return 0.0;
    }
    double dot = (a / la).Dot(b / lb);
    dot = std::max(-1.0, std::min(1.0, dot));
    return std::acos(dot);
}

static NativeMetricSummary SummarizeNative(const std::vector<NativeFrameRecord>& frames) {
    NativeMetricSummary out;
    out.frames = static_cast<int>(frames.size());
    if (frames.empty()) {
        return out;
    }

    double count_sum = 0.0;
    double count_change_sum = 0.0;
    int count_change_samples = 0;
    double force_jump_sum = 0.0;
    double force_jump_sq_sum = 0.0;
    int force_jump_samples = 0;
    double torque_jump_sum = 0.0;
    double torque_jump_sq_sum = 0.0;
    int torque_jump_samples = 0;
    double contact_normal_jump_sq_sum = 0.0;
    int contact_normal_jump_samples = 0;
    double direction_jump_sq_sum = 0.0;
    int direction_jump_samples = 0;
    double pos_err_sq_sum = 0.0;
    double vel_err_sq_sum = 0.0;
    double acc_err_sq_sum = 0.0;

    for (size_t i = 0; i < frames.size(); i++) {
        const auto& row = frames[i];
        count_sum += row.contact_count;
        out.max_contact_count = std::max(out.max_contact_count, row.contact_count);
        if (row.contact_count > 0) {
            out.active_frames++;
        }
        out.max_force_norm = std::max(out.max_force_norm, row.force.Length());
        out.max_torque_norm = std::max(out.max_torque_norm, row.torque.Length());
        out.max_position_error = std::max(out.max_position_error, row.pos_error.Length());
        out.max_velocity_error = std::max(out.max_velocity_error, row.vel_error.Length());
        out.max_acceleration_error = std::max(out.max_acceleration_error, row.acc_error.Length());
        pos_err_sq_sum += row.pos_error.Length2();
        vel_err_sq_sum += row.vel_error.Length2();
        acc_err_sq_sum += row.acc_error.Length2();

        if (i > 0) {
            const auto& prev = frames[i - 1];
            int count_change = std::abs(row.contact_count - prev.contact_count);
            count_change_sum += count_change;
            out.max_abs_contact_count_change = std::max(out.max_abs_contact_count_change, count_change);
            count_change_samples++;

            double force_jump = (row.force - prev.force).Length();
            force_jump_sum += force_jump;
            force_jump_sq_sum += force_jump * force_jump;
            out.max_force_jump = std::max(out.max_force_jump, force_jump);
            force_jump_samples++;

            double torque_jump = (row.torque - prev.torque).Length();
            torque_jump_sum += torque_jump;
            torque_jump_sq_sum += torque_jump * torque_jump;
            out.max_torque_jump = std::max(out.max_torque_jump, torque_jump);
            torque_jump_samples++;

            double direction_jump = DirectionJumpAngle(row.force, prev.force);
            if (row.force.Length() > 1.0e-12 && prev.force.Length() > 1.0e-12) {
                direction_jump_sq_sum += direction_jump * direction_jump;
                out.max_force_direction_jump_angle = std::max(out.max_force_direction_jump_angle, direction_jump);
                direction_jump_samples++;
            }

            double contact_normal_jump = DirectionJumpAngle(row.contact_normal, prev.contact_normal);
            if (row.contact_normal.Length() > 1.0e-12 && prev.contact_normal.Length() > 1.0e-12) {
                contact_normal_jump_sq_sum += contact_normal_jump * contact_normal_jump;
                out.max_contact_normal_jump_angle = std::max(out.max_contact_normal_jump_angle, contact_normal_jump);
                contact_normal_jump_samples++;
            }
        }
    }

    out.mean_contact_count = count_sum / static_cast<double>(frames.size());
    if (count_change_samples > 0) {
        out.mean_abs_contact_count_change = count_change_sum / static_cast<double>(count_change_samples);
    }
    if (force_jump_samples > 0) {
        out.mean_force_jump = force_jump_sum / static_cast<double>(force_jump_samples);
        out.rms_force_jump = std::sqrt(force_jump_sq_sum / static_cast<double>(force_jump_samples));
        out.force_oscillation_index = SafeRatio(out.rms_force_jump, out.max_force_norm);
    }
    if (torque_jump_samples > 0) {
        out.mean_torque_jump = torque_jump_sum / static_cast<double>(torque_jump_samples);
        out.rms_torque_jump = std::sqrt(torque_jump_sq_sum / static_cast<double>(torque_jump_samples));
        out.torque_oscillation_index = SafeRatio(out.rms_torque_jump, out.max_torque_norm);
    }
    if (direction_jump_samples > 0) {
        out.rms_force_direction_jump_angle =
            std::sqrt(direction_jump_sq_sum / static_cast<double>(direction_jump_samples));
    }
    if (contact_normal_jump_samples > 0) {
        out.rms_contact_normal_jump_angle =
            std::sqrt(contact_normal_jump_sq_sum / static_cast<double>(contact_normal_jump_samples));
    }
    out.rms_position_error = std::sqrt(pos_err_sq_sum / static_cast<double>(frames.size()));
    out.rms_velocity_error = std::sqrt(vel_err_sq_sum / static_cast<double>(frames.size()));
    out.rms_acceleration_error = std::sqrt(acc_err_sq_sum / static_cast<double>(frames.size()));
    return out;
}

static FieldContactTopologyMetricsSummary RunFieldScenario(const ScenarioConfig& scenario,
                                                           const SceneData& scene,
                                                           std::ofstream& frames_out) {
    const double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    FieldContactRuntimeSettings settings = MakeRuntimeSettings(dt, 1.35 * scenario.voxel_size);
    FieldContactPrimitiveTracker tracker;
    FieldContactTopologyMetricsAccumulator metrics;
    tracker.Reset();
    metrics.Reset();

    for (int frame = 0; frame < scenario.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(scenario, frame);
        auto queries = BuildQueries(scene.moving_graph, scene.sdf, kin);
        FieldContactStepResult step = tracker.Evaluate(scene.moving_graph, queries, kin.center, settings);
        metrics.Accumulate(step);

        FieldFrameRecord row;
        row.scenario = scenario.name;
        row.frame = frame;
        row.time = time;
        row.kin = kin;
        row.patch_count = step.stats.patch_count;
        row.newborn_count = step.stats.newborn_count;
        row.merge_count = step.stats.merge_count;
        row.split_count = step.stats.split_count;
        row.death_count = step.stats.death_count;
        row.max_source_count = step.stats.max_source_count;
        row.max_previous_reuse = step.stats.max_previous_reuse;
        row.force = step.total_force;
        row.torque = step.total_torque;
        row.max_tangential_force_ratio = step.stats.max_tangential_force_ratio;
        row.max_inherited_energy_ratio = step.stats.max_inherited_energy_ratio;

        frames_out << row.scenario << "," << row.frame << "," << row.time << ","
                   << kin.center.x() << "," << kin.center.y() << "," << kin.center.z() << ","
                   << kin.velocity.x() << "," << kin.velocity.y() << "," << kin.velocity.z() << ","
                   << kin.acceleration.x() << "," << kin.acceleration.y() << "," << kin.acceleration.z() << ","
                   << row.patch_count << "," << row.newborn_count << "," << row.merge_count << ","
                   << row.split_count << "," << row.death_count << "," << row.max_source_count << ","
                   << row.max_previous_reuse << "," << row.force.x() << "," << row.force.y() << ","
                   << row.force.z() << "," << row.force.Length() << "," << row.torque.x() << ","
                   << row.torque.y() << "," << row.torque.z() << "," << row.torque.Length() << ","
                   << row.max_tangential_force_ratio << "," << row.max_inherited_energy_ratio << "\n";
    }

    return metrics.GetSummary();
}

static NativeMetricSummary RunNativeScenario(const ScenarioConfig& scenario,
                                             const SceneData& scene,
                                             std::ofstream& frames_out) {
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(120);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-7);

    auto mat = MakeMaterial();
    ChCollisionInfo::SetDefaultEffectiveCurvatureRadius(scenario.probe_radius);
    AddMatchedMeshTarget(sys, mat, scene.target_mesh, scenario.native_mesh_radius);
    auto moving = AddMovingSphere(sys, mat, scenario.probe_radius);

    const double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    ChVector3d previous_actual_vel(0, 0, 0);
    bool has_previous_actual_vel = false;
    std::vector<NativeFrameRecord> rows;
    rows.reserve(static_cast<size_t>(scenario.frames));

    for (int frame = 0; frame < scenario.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(scenario, frame);
        ApplyCommandedKinematics(moving, kin);

        sys.DoStepDynamics(dt);

        NativeFrameRecord row;
        row.scenario = scenario.name;
        row.frame = frame;
        row.time = time;
        row.cmd = kin;
        row.actual_pos = moving->GetPos();
        row.actual_vel = moving->GetPosDt();
        row.actual_acc = has_previous_actual_vel ? (row.actual_vel - previous_actual_vel) / dt : kin.acceleration;
        previous_actual_vel = row.actual_vel;
        has_previous_actual_vel = true;
        row.pos_error = row.actual_pos - kin.center;
        row.vel_error = row.actual_vel - kin.velocity;
        row.acc_error = row.actual_acc - kin.acceleration;
        row.contact_count = static_cast<int>(sys.GetNumContacts());
        row.force = moving->GetContactForce();
        row.torque = moving->GetContactTorque();
        auto normal_reporter = chrono_types::make_shared<NativeContactNormalReporter>();
        sys.GetContactContainer()->ReportAllContacts(normal_reporter);
        row.contact_normal = normal_reporter->GetAverageNormal(row.force);
        rows.push_back(row);

        frames_out << row.scenario << "," << row.frame << "," << row.time << ","
                   << kin.center.x() << "," << kin.center.y() << "," << kin.center.z() << ","
                   << kin.velocity.x() << "," << kin.velocity.y() << "," << kin.velocity.z() << ","
                   << kin.acceleration.x() << "," << kin.acceleration.y() << "," << kin.acceleration.z() << ","
                   << row.actual_pos.x() << "," << row.actual_pos.y() << "," << row.actual_pos.z() << ","
                   << row.actual_vel.x() << "," << row.actual_vel.y() << "," << row.actual_vel.z() << ","
                   << row.actual_acc.x() << "," << row.actual_acc.y() << "," << row.actual_acc.z() << ","
                   << row.pos_error.Length() << "," << row.vel_error.Length() << "," << row.acc_error.Length() << ","
                   << row.contact_count << "," << row.contact_normal.x() << "," << row.contact_normal.y() << ","
                   << row.contact_normal.z() << "," << row.force.x() << "," << row.force.y() << "," << row.force.z()
                   << "," << row.force.Length() << "," << row.torque.x() << "," << row.torque.y() << ","
                   << row.torque.z() << "," << row.torque.Length() << "\n";
    }

    return SummarizeNative(rows);
}

static void WriteComparisonHeader(std::ofstream& out) {
    out << "scenario,"
           "field_frames,field_active_frames,field_mean_patch_count,field_max_patch_count,"
           "field_mean_abs_patch_count_change,field_max_abs_patch_count_change,"
           "field_rms_normal_jump_angle,field_max_normal_jump_angle,"
           "field_force_oscillation_index,field_torque_oscillation_index,"
           "field_rms_force_jump,field_rms_torque_jump,field_max_tangential_force_ratio,"
           "native_frames,native_active_frames,native_mean_contact_count,native_max_contact_count,"
           "native_mean_abs_contact_count_change,native_max_abs_contact_count_change,"
           "native_rms_contact_normal_jump_angle,native_max_contact_normal_jump_angle,"
           "native_rms_force_direction_jump_angle,native_max_force_direction_jump_angle,"
           "native_force_oscillation_index,native_torque_oscillation_index,"
           "native_rms_force_jump,native_rms_torque_jump,"
           "native_max_position_error,native_max_velocity_error,native_max_acceleration_error,"
           "native_rms_position_error,native_rms_velocity_error,native_rms_acceleration_error,"
           "native_to_field_count_churn_ratio,native_to_field_max_count_jump_ratio,"
           "native_contact_normal_to_field_normal_jump_ratio,"
           "native_force_direction_to_field_normal_jump_ratio,"
           "native_to_field_force_oscillation_ratio,native_to_field_torque_oscillation_ratio\n";
}

static void WriteComparisonRow(std::ofstream& out, const ScenarioResult& result) {
    const auto& f = result.field_summary;
    const auto& n = result.native_summary;
    out << result.config.name << ","
        << f.frames << "," << f.active_frames << "," << f.mean_patch_count << "," << f.max_patch_count << ","
        << f.mean_abs_patch_count_change << "," << f.max_abs_patch_count_change << ","
        << f.rms_normal_jump_angle << "," << f.max_normal_jump_angle << ","
        << f.force_oscillation_index << "," << f.torque_oscillation_index << ","
        << f.rms_force_jump << "," << f.rms_torque_jump << "," << f.max_tangential_force_ratio << ","
        << n.frames << "," << n.active_frames << "," << n.mean_contact_count << "," << n.max_contact_count << ","
        << n.mean_abs_contact_count_change << "," << n.max_abs_contact_count_change << ","
        << n.rms_contact_normal_jump_angle << "," << n.max_contact_normal_jump_angle << ","
        << n.rms_force_direction_jump_angle << "," << n.max_force_direction_jump_angle << ","
        << n.force_oscillation_index << "," << n.torque_oscillation_index << ","
        << n.rms_force_jump << "," << n.rms_torque_jump << ","
        << n.max_position_error << "," << n.max_velocity_error << "," << n.max_acceleration_error << ","
        << n.rms_position_error << "," << n.rms_velocity_error << "," << n.rms_acceleration_error << ","
        << SafeRatio(n.mean_abs_contact_count_change, f.mean_abs_patch_count_change) << ","
        << SafeRatio(static_cast<double>(n.max_abs_contact_count_change),
                     static_cast<double>(f.max_abs_patch_count_change))
        << "," << SafeRatio(n.rms_contact_normal_jump_angle, f.rms_normal_jump_angle) << ","
        << SafeRatio(n.rms_force_direction_jump_angle, f.rms_normal_jump_angle) << ","
        << SafeRatio(n.force_oscillation_index, f.force_oscillation_index) << ","
        << SafeRatio(n.torque_oscillation_index, f.torque_oscillation_index) << "\n";
}

}  // namespace

int main() {
    std::cout << "Milestone 27: matched native mesh/OpenVDB field baseline" << std::endl;
    openvdb::initialize();

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_27";
    std::filesystem::create_directories(out_dir);

    std::ofstream field_frames(out_dir / "matched_mesh_field_frames.csv");
    std::ofstream native_frames(out_dir / "matched_mesh_native_frames.csv");
    std::ofstream comparison(out_dir / "matched_mesh_comparison.csv");
    std::ofstream geometry(out_dir / "matched_mesh_geometry.csv");

    field_frames << std::fixed << std::setprecision(8);
    native_frames << std::fixed << std::setprecision(8);
    comparison << std::fixed << std::setprecision(8);
    geometry << std::fixed << std::setprecision(8);

    field_frames << "scenario,frame,time,center_x,center_y,center_z,velocity_x,velocity_y,velocity_z,"
                    "acceleration_x,acceleration_y,acceleration_z,patch_count,newborn_count,merge_count,"
                    "split_count,death_count,max_source_count,max_previous_reuse,total_force_x,total_force_y,"
                    "total_force_z,total_force_norm,total_torque_x,total_torque_y,total_torque_z,total_torque_norm,"
                    "max_tangential_force_ratio,max_inherited_energy_ratio\n";
    native_frames << "scenario,frame,time,cmd_pos_x,cmd_pos_y,cmd_pos_z,cmd_vel_x,cmd_vel_y,cmd_vel_z,"
                     "cmd_acc_x,cmd_acc_y,cmd_acc_z,actual_pos_x,actual_pos_y,actual_pos_z,"
                     "actual_vel_x,actual_vel_y,actual_vel_z,actual_acc_x,actual_acc_y,actual_acc_z,"
                     "pos_error_norm,vel_error_norm,acc_error_norm,contact_count,"
                     "contact_normal_x,contact_normal_y,contact_normal_z,"
                     "contact_force_x,contact_force_y,contact_force_z,contact_force_norm,"
                     "contact_torque_x,contact_torque_y,contact_torque_z,contact_torque_norm\n";
    WriteComparisonHeader(comparison);
    geometry << "scenario,target_vertices,target_faces,moving_samples,voxel_size,active_voxels,native_mesh_radius\n";

    std::vector<ScenarioConfig> scenarios = {
        {MatchedScenario::UChannel,
         "matched_u_channel_sphere_slide",
         ChVector3d(-0.54, 0.116, 0.000),
         ChVector3d(0.54, 0.116, 0.000),
         0.120,
         0.045,
         1.0,
         1.25,
         501,
         1.0,
         0.006,
         0.0},
        {MatchedScenario::StaggeredPadTrack,
         "matched_staggered_pad_track",
         ChVector3d(-0.54, 0.145, -0.074),
         ChVector3d(0.54, 0.145, 0.074),
         0.120,
         0.036,
         2.0,
         1.35,
         541,
         1.0,
         0.006,
         0.0},
    };

    std::vector<ScenarioResult> results;
    for (const auto& scenario : scenarios) {
        std::cout << "Scenario " << scenario.name << " (" << scenario.frames << " frames)" << std::endl;
        SceneData scene = BuildScene(scenario);
        geometry << scenario.name << "," << scene.target_mesh.vertices.size() << "," << scene.target_mesh.faces.size()
                 << "," << scene.moving_graph.samples.size() << "," << scene.sdf.voxel_size << ","
                 << scene.sdf.grid->activeVoxelCount() << "," << scenario.native_mesh_radius << "\n";

        ScenarioResult result;
        result.config = scenario;
        result.field_summary = RunFieldScenario(scenario, scene, field_frames);
        result.native_summary = RunNativeScenario(scenario, scene, native_frames);
        results.push_back(result);

        std::cout << "  field mean_patch=" << result.field_summary.mean_patch_count
                  << " mean_dcount=" << result.field_summary.mean_abs_patch_count_change
                  << " normal_rms=" << result.field_summary.rms_normal_jump_angle
                  << " coulomb=" << result.field_summary.max_tangential_force_ratio << std::endl;
        std::cout << "  native mean_contacts=" << result.native_summary.mean_contact_count
                  << " mean_dcount=" << result.native_summary.mean_abs_contact_count_change
                  << " normal_rms=" << result.native_summary.rms_contact_normal_jump_angle
                  << " force_dir_rms=" << result.native_summary.rms_force_direction_jump_angle
                  << " pos_err_max=" << result.native_summary.max_position_error << std::endl;
    }

    for (const auto& result : results) {
        WriteComparisonRow(comparison, result);
    }

    std::cout << "Wrote:" << std::endl;
    std::cout << "  " << (out_dir / "matched_mesh_field_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "matched_mesh_native_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "matched_mesh_comparison.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "matched_mesh_geometry.csv").string() << std::endl;
    return 0;
}
