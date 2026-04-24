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
// Milestone 29: real CAD/mesh case study.
//
// This benchmark uses existing Chrono CAD/mesh-derived OBJ assets:
//   data/models/bulldozer/shoe_collision.obj
//   data/models/semicapsule.obj
//
// The track-shoe mesh is converted to an OpenVDB SDF and the semicapsule mesh is
// used as the moving field-contact surface graph. The same OBJ meshes are also
// used for a Chrono native mesh-collision baseline.
//
// Outputs:
//   out/milestone_29/cad_mesh_case_study_frames.csv
//   out/milestone_29/cad_mesh_case_study_comparison.csv
//   out/milestone_29/cad_mesh_case_study_geometry.csv
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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "chrono/collision/ChCollisionShapeConvexHull.h"
#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/collision/ChFieldContactRuntime.h"
#include "chrono/core/ChRotation.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChSystemSMC.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kFrictionCoefficient = 0.45;

struct DemoTriangleMesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct MeshBounds {
    ChVector3d min = ChVector3d(0, 0, 0);
    ChVector3d max = ChVector3d(0, 0, 0);
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d size = ChVector3d(0, 0, 0);
};

struct OpenVDBSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 0.004;
};

struct CadCaseStudyConfig {
    std::string scenario = "cad_track_shoe_semicapsule_slide";
    std::filesystem::path target_obj;
    std::filesystem::path moving_obj;
    int frames = 701;
    double total_time = 1.40;
    double voxel_size = 0.0015;
    double moving_scale = 0.22;
    double native_target_radius = 0.0;
    double native_moving_radius = 0.0;
};

struct SceneData {
    DemoTriangleMesh target_mesh;
    DemoTriangleMesh moving_mesh_local;
    SurfaceGraph moving_graph;
    MeshBounds target_bounds;
    MeshBounds moving_bounds_original;
    MeshBounds moving_bounds_local;
    OpenVDBSDF sdf;
};

struct BodyKinematics {
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
    double yaw = 0.0;
    double roll = 0.0;
    ChVector3d angular_velocity = ChVector3d(0, 0, 0);
};

struct NativeFrame {
    int contact_count = 0;
    ChVector3d contact_normal = ChVector3d(0, 0, 0);
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
    double max_penetration = 0.0;
    double max_coulomb_ratio = 0.0;
};

struct NativeSummary {
    int frames = 0;
    int active_frames = 0;
    double mean_contact_count = 0.0;
    int max_contact_count = 0;
    double mean_abs_contact_count_change = 0.0;
    int max_abs_contact_count_change = 0;
    double rms_contact_normal_jump_angle = 0.0;
    double max_contact_normal_jump_angle = 0.0;
    double rms_force_jump = 0.0;
    double max_force_jump = 0.0;
    double max_force_norm = 0.0;
    double force_oscillation_index = 0.0;
    double rms_torque_jump = 0.0;
    double max_torque_jump = 0.0;
    double max_torque_norm = 0.0;
    double torque_oscillation_index = 0.0;
    double max_penetration = 0.0;
    double max_coulomb_ratio = 0.0;
};

struct NativeContactReporter : public ChContactContainer::ReportContactCallback {
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
        (void)eff_radius;
        (void)react_torques;
        (void)contactobjA;
        (void)contactobjB;
        (void)constraint_offset;

        count++;
        max_penetration = std::max(max_penetration, std::max(0.0, -distance));

        ChVector3d normal = SafeNormalize(plane_coord.GetAxisX(), ChVector3d(0, 1, 0));
        if (react_forces.x() < 0.0) {
            normal *= -1.0;
        }
        double normal_weight = std::max(1.0e-12, std::abs(react_forces.x()));
        weighted_normal += normal * normal_weight;
        normal_weight_sum += normal_weight;

        double fn = std::abs(react_forces.x());
        double ft = std::sqrt(react_forces.y() * react_forces.y() + react_forces.z() * react_forces.z());
        if (fn > 1.0e-12) {
            max_coulomb_ratio = std::max(max_coulomb_ratio, ft / (kFrictionCoefficient * fn));
        }
        return true;
    }

    ChVector3d AverageNormal(const ChVector3d& total_force) const {
        if (normal_weight_sum <= 1.0e-14 || count == 0) {
            return ChVector3d(0, 0, 0);
        }
        ChVector3d normal = SafeNormalize(weighted_normal / normal_weight_sum, ChVector3d(0, 1, 0));
        if (total_force.Length() > 1.0e-12 && normal.Dot(total_force) < 0.0) {
            normal *= -1.0;
        }
        return normal;
    }

    int count = 0;
    ChVector3d weighted_normal = ChVector3d(0, 0, 0);
    double normal_weight_sum = 0.0;
    double max_penetration = 0.0;
    double max_coulomb_ratio = 0.0;
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

static MeshBounds ComputeBounds(const DemoTriangleMesh& mesh) {
    MeshBounds b;
    if (mesh.vertices.empty()) {
        return b;
    }
    b.min = mesh.vertices.front();
    b.max = mesh.vertices.front();
    for (const auto& v : mesh.vertices) {
        b.min.x() = std::min(b.min.x(), v.x());
        b.min.y() = std::min(b.min.y(), v.y());
        b.min.z() = std::min(b.min.z(), v.z());
        b.max.x() = std::max(b.max.x(), v.x());
        b.max.y() = std::max(b.max.y(), v.y());
        b.max.z() = std::max(b.max.z(), v.z());
    }
    b.center = 0.5 * (b.min + b.max);
    b.size = b.max - b.min;
    return b;
}

static int ParseObjIndex(const std::string& token, int vertex_count) {
    std::string head = token;
    size_t slash = head.find('/');
    if (slash != std::string::npos) {
        head = head.substr(0, slash);
    }
    if (head.empty()) {
        return -1;
    }
    int raw = std::stoi(head);
    if (raw > 0) {
        return raw - 1;
    }
    if (raw < 0) {
        return vertex_count + raw;
    }
    return -1;
}

static DemoTriangleMesh LoadObjTriangleMesh(const std::filesystem::path& path) {
    DemoTriangleMesh mesh;
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open OBJ mesh: " + path.string());
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        std::string tag;
        iss >> tag;
        if (tag == "v") {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            iss >> x >> y >> z;
            mesh.vertices.push_back(ChVector3d(x, y, z));
        } else if (tag == "f") {
            std::vector<int> ids;
            std::string token;
            while (iss >> token) {
                int id = ParseObjIndex(token, static_cast<int>(mesh.vertices.size()));
                if (id >= 0 && id < static_cast<int>(mesh.vertices.size())) {
                    ids.push_back(id);
                }
            }
            for (size_t i = 1; i + 1 < ids.size(); i++) {
                mesh.faces.push_back({ids[0], ids[i], ids[i + 1]});
            }
        }
    }
    return mesh;
}

static DemoTriangleMesh RecenterMesh(const DemoTriangleMesh& mesh, const ChVector3d& center, double scale) {
    DemoTriangleMesh out = mesh;
    for (auto& v : out.vertices) {
        v = (v - center) * scale;
    }
    return out;
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

static SceneData BuildScene(const CadCaseStudyConfig& config) {
    SceneData scene;
    scene.sdf.voxel_size = config.voxel_size;
    scene.target_mesh = LoadObjTriangleMesh(config.target_obj);
    DemoTriangleMesh moving_original = LoadObjTriangleMesh(config.moving_obj);
    scene.target_bounds = ComputeBounds(scene.target_mesh);
    scene.moving_bounds_original = ComputeBounds(moving_original);
    scene.moving_mesh_local = RecenterMesh(moving_original, scene.moving_bounds_original.center, config.moving_scale);
    scene.moving_bounds_local = ComputeBounds(scene.moving_mesh_local);
    scene.moving_graph = MakeTriangleMeshSurfaceGraph(scene.moving_mesh_local.vertices, scene.moving_mesh_local.faces);
    scene.sdf.grid = BuildLevelSetFromTriangleMesh(scene.target_mesh, config.voxel_size, 8.0f);
    return scene;
}

static FieldContactRuntimeSettings MakeRuntimeSettings(double dt, double activation_band) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = activation_band;
    settings.extraction.min_area = 1.0e-9;
    settings.extraction.min_samples = 2;
    settings.extraction.use_penetration_weighted_center = true;

    settings.normal.stiffness = 1.2e6;
    settings.normal.damping = 0.0;

    settings.tangential.stiffness = 2.6e3;
    settings.tangential.damping = 0.0;
    settings.tangential.friction_coefficient = kFrictionCoefficient;
    settings.tangential.time_step = dt;

    settings.inheritance.min_overlap = 0.004;
    settings.inheritance.min_normal_dot = 0.14;
    settings.inheritance.max_center_distance = 0.050;
    settings.inheritance.geometry_fallback_weight = 0.18;
    return settings;
}

static ChVector3d RotateY(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(c * v.x() + s * v.z(), v.y(), -s * v.x() + c * v.z());
}

static ChVector3d RotateZ(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z());
}

static ChVector3d RotateLocal(const ChVector3d& v, double yaw, double roll) {
    return RotateY(RotateZ(v, roll), yaw);
}

static BodyKinematics KinematicsAtFrame(const CadCaseStudyConfig& config, int frame) {
    double u = static_cast<double>(frame) / static_cast<double>(config.frames - 1);
    double du = 1.0 / static_cast<double>(config.frames - 1);

    auto center_at = [](double uu) {
        double clamped = std::max(0.0, std::min(1.0, uu));
        double x = 0.055 + 0.145 * clamped;
        double z = 0.074 * std::sin(2.0 * kPi * 1.25 * clamped) + 0.020 * std::sin(2.0 * kPi * 3.0 * clamped);
        double y = 0.1088 + 0.0015 * std::sin(2.0 * kPi * 0.75 * clamped);
        return ChVector3d(x, y, z);
    };

    auto yaw_at = [](double uu) {
        double clamped = std::max(0.0, std::min(1.0, uu));
        return 0.55 * std::sin(2.0 * kPi * 1.7 * clamped) + 0.35 * clamped;
    };

    auto roll_at = [](double uu) {
        double clamped = std::max(0.0, std::min(1.0, uu));
        return 0.22 * std::sin(2.0 * kPi * 2.3 * clamped);
    };

    double u0 = std::max(0.0, u - du);
    double u1 = std::min(1.0, u + du);
    double denom = (u1 - u0) * config.total_time;

    BodyKinematics kin;
    kin.center = center_at(u);
    kin.velocity = denom > 0.0 ? (center_at(u1) - center_at(u0)) / denom : ChVector3d(0, 0, 0);
    kin.yaw = yaw_at(u);
    kin.roll = roll_at(u);
    double yaw_rate = denom > 0.0 ? (yaw_at(u1) - yaw_at(u0)) / denom : 0.0;
    double roll_rate = denom > 0.0 ? (roll_at(u1) - roll_at(u0)) / denom : 0.0;
    kin.angular_velocity = ChVector3d(0.0, yaw_rate, roll_rate);
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
        ChVector3d rotated = RotateLocal(sample.local_pos, kin.yaw, kin.roll);
        ChVector3d world_pos = kin.center + rotated;
        ChVector3d world_vel = kin.velocity + kin.angular_velocity.Cross(rotated);
        queries.push_back(QueryOpenVDB(sdf, sampler, world_pos, world_vel));
    }
    return queries;
}

static double MaxPatchPenetration(const FieldContactStepResult& step) {
    double out = 0.0;
    for (const auto& patch : step.patches) {
        out = std::max(out, patch.patch.max_penetration);
    }
    return out;
}

static int CountSlipPatches(const FieldContactStepResult& step) {
    int count = 0;
    for (const auto& patch : step.patches) {
        if (patch.tangential.state == StickSlipState::Slip) {
            count++;
        }
    }
    return count;
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

static double SafeRatio(double numerator, double denominator) {
    return denominator > 1.0e-14 ? numerator / denominator : 0.0;
}

static std::shared_ptr<ChContactMaterialSMC> MakeMaterial() {
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetFriction(static_cast<float>(kFrictionCoefficient));
    mat->SetRestitution(0.02f);
    mat->SetYoungModulus(1.2e4f);
    mat->SetPoissonRatio(0.30f);
    return mat;
}

static NativeSummary SummarizeNative(const std::vector<NativeFrame>& frames) {
    NativeSummary out;
    out.frames = static_cast<int>(frames.size());
    if (frames.empty()) {
        return out;
    }

    double count_sum = 0.0;
    double count_change_sum = 0.0;
    int count_change_samples = 0;
    double normal_jump_sq_sum = 0.0;
    int normal_jump_samples = 0;
    double force_jump_sq_sum = 0.0;
    int force_jump_samples = 0;
    double torque_jump_sq_sum = 0.0;
    int torque_jump_samples = 0;

    for (size_t i = 0; i < frames.size(); i++) {
        const auto& row = frames[i];
        count_sum += row.contact_count;
        out.max_contact_count = std::max(out.max_contact_count, row.contact_count);
        out.max_force_norm = std::max(out.max_force_norm, row.force.Length());
        out.max_torque_norm = std::max(out.max_torque_norm, row.torque.Length());
        out.max_penetration = std::max(out.max_penetration, row.max_penetration);
        out.max_coulomb_ratio = std::max(out.max_coulomb_ratio, row.max_coulomb_ratio);
        if (row.contact_count > 0) {
            out.active_frames++;
        }

        if (i > 0) {
            const auto& prev = frames[i - 1];
            int count_change = std::abs(row.contact_count - prev.contact_count);
            count_change_sum += count_change;
            out.max_abs_contact_count_change = std::max(out.max_abs_contact_count_change, count_change);
            count_change_samples++;

            double normal_jump = DirectionJumpAngle(row.contact_normal, prev.contact_normal);
            if (row.contact_normal.Length() > 1.0e-12 && prev.contact_normal.Length() > 1.0e-12) {
                normal_jump_sq_sum += normal_jump * normal_jump;
                out.max_contact_normal_jump_angle = std::max(out.max_contact_normal_jump_angle, normal_jump);
                normal_jump_samples++;
            }

            double force_jump = (row.force - prev.force).Length();
            force_jump_sq_sum += force_jump * force_jump;
            out.max_force_jump = std::max(out.max_force_jump, force_jump);
            force_jump_samples++;

            double torque_jump = (row.torque - prev.torque).Length();
            torque_jump_sq_sum += torque_jump * torque_jump;
            out.max_torque_jump = std::max(out.max_torque_jump, torque_jump);
            torque_jump_samples++;
        }
    }

    out.mean_contact_count = count_sum / static_cast<double>(frames.size());
    if (count_change_samples > 0) {
        out.mean_abs_contact_count_change = count_change_sum / static_cast<double>(count_change_samples);
    }
    if (normal_jump_samples > 0) {
        out.rms_contact_normal_jump_angle = std::sqrt(normal_jump_sq_sum / static_cast<double>(normal_jump_samples));
    }
    if (force_jump_samples > 0) {
        out.rms_force_jump = std::sqrt(force_jump_sq_sum / static_cast<double>(force_jump_samples));
        out.force_oscillation_index = SafeRatio(out.rms_force_jump, out.max_force_norm);
    }
    if (torque_jump_samples > 0) {
        out.rms_torque_jump = std::sqrt(torque_jump_sq_sum / static_cast<double>(torque_jump_samples));
        out.torque_oscillation_index = SafeRatio(out.rms_torque_jump, out.max_torque_norm);
    }
    return out;
}

static std::shared_ptr<ChBody> AddMeshBody(ChSystemSMC& sys,
                                           const std::shared_ptr<ChContactMaterialSMC>& mat,
                                           const DemoTriangleMesh& mesh,
                                           bool fixed,
                                           bool is_static,
                                           bool is_convex,
                                           double radius) {
    auto body = chrono_types::make_shared<ChBody>();
    body->SetFixed(fixed);
    body->SetMass(1.0);
    body->SetInertiaXX(ChVector3d(1.0e-3, 1.0e-3, 1.0e-3));
    body->EnableCollision(true);
    if (is_convex && !is_static) {
        body->AddCollisionShape(chrono_types::make_shared<ChCollisionShapeConvexHull>(mat, mesh.vertices));
    } else {
        auto shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat,
                                                                             ToChronoTriangleMesh(mesh),
                                                                             is_static,
                                                                             is_convex,
                                                                             radius);
        body->AddCollisionShape(shape);
    }
    sys.AddBody(body);
    return body;
}

static NativeSummary RunNativeBaseline(const CadCaseStudyConfig& config,
                                       const SceneData& scene,
                                       std::ofstream& frames_out) {
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(120);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-7);

    auto mat = MakeMaterial();
    AddMeshBody(sys, mat, scene.target_mesh, true, true, false, config.native_target_radius);
    auto moving = AddMeshBody(sys, mat, scene.moving_mesh_local, false, false, true, config.native_moving_radius);

    double dt = config.total_time / static_cast<double>(config.frames - 1);
    std::vector<NativeFrame> native_frames;
    native_frames.reserve(static_cast<size_t>(config.frames));

    for (int frame = 0; frame < config.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(config, frame);
        moving->SetPos(kin.center);
        moving->SetRot(QuatFromAngleY(kin.yaw) * QuatFromAngleZ(kin.roll));
        moving->SetPosDt(kin.velocity);
        moving->SetAngVelParent(kin.angular_velocity);

        sys.DoStepDynamics(dt);

        NativeFrame native;
        native.contact_count = static_cast<int>(sys.GetNumContacts());
        native.force = moving->GetContactForce();
        native.torque = moving->GetContactTorque();
        auto reporter = chrono_types::make_shared<NativeContactReporter>();
        sys.GetContactContainer()->ReportAllContacts(reporter);
        native.max_penetration = reporter->max_penetration;
        native.max_coulomb_ratio = reporter->max_coulomb_ratio;
        native.contact_normal = reporter->AverageNormal(native.force);
        native_frames.push_back(native);

        frames_out << config.scenario << ",native_mesh_baseline," << frame << "," << time << ","
                   << kin.center.x() << "," << kin.center.y() << "," << kin.center.z() << ","
                   << kin.velocity.x() << "," << kin.velocity.y() << "," << kin.velocity.z() << ","
                   << kin.yaw << "," << kin.roll << ","
                   << native.contact_count << ",0,0,0,0,0,0,0,"
                   << native.force.x() << "," << native.force.y() << "," << native.force.z() << ","
                   << native.force.Length() << "," << native.torque.x() << "," << native.torque.y() << ","
                   << native.torque.z() << "," << native.torque.Length() << ","
                   << native.max_penetration << "," << native.max_coulomb_ratio << "\n";
    }

    return SummarizeNative(native_frames);
}

static FieldContactTopologyMetricsSummary RunFieldCase(const CadCaseStudyConfig& config,
                                                       const SceneData& scene,
                                                       std::ofstream& frames_out,
                                                       double& max_penetration) {
    double dt = config.total_time / static_cast<double>(config.frames - 1);
    FieldContactRuntimeSettings settings = MakeRuntimeSettings(dt, 1.35 * config.voxel_size);
    FieldContactPrimitiveTracker tracker;
    FieldContactTopologyMetricsAccumulator metrics;
    tracker.Reset();
    metrics.Reset();
    max_penetration = 0.0;

    for (int frame = 0; frame < config.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(config, frame);
        auto queries = BuildQueries(scene.moving_graph, scene.sdf, kin);
        FieldContactStepResult step = tracker.Evaluate(scene.moving_graph, queries, kin.center, settings);
        metrics.Accumulate(step);
        double frame_penetration = MaxPatchPenetration(step);
        max_penetration = std::max(max_penetration, frame_penetration);

        frames_out << config.scenario << ",field_primitive_openvdb," << frame << "," << time << ","
                   << kin.center.x() << "," << kin.center.y() << "," << kin.center.z() << ","
                   << kin.velocity.x() << "," << kin.velocity.y() << "," << kin.velocity.z() << ","
                   << kin.yaw << "," << kin.roll << ","
                   << "0," << step.stats.patch_count << "," << step.stats.newborn_count << ","
                   << step.stats.merge_count << "," << step.stats.split_count << "," << step.stats.death_count << ","
                   << step.stats.max_source_count << "," << CountSlipPatches(step) << ","
                   << step.total_force.x() << "," << step.total_force.y() << "," << step.total_force.z() << ","
                   << step.total_force.Length() << "," << step.total_torque.x() << "," << step.total_torque.y() << ","
                   << step.total_torque.z() << "," << step.total_torque.Length() << ","
                   << frame_penetration << "," << step.stats.max_tangential_force_ratio << "\n";
    }

    return metrics.GetSummary();
}

}  // namespace

int main() {
    std::cout << "Milestone 29: CAD mesh case study" << std::endl;
    openvdb::initialize();

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_29";
    std::filesystem::create_directories(out_dir);

    CadCaseStudyConfig config;
    config.target_obj = std::filesystem::path(project_root) / "data" / "models" / "bulldozer" / "shoe_collision.obj";
    config.moving_obj = std::filesystem::path(project_root) / "data" / "models" / "semicapsule.obj";

    SceneData scene = BuildScene(config);

    std::ofstream frames_out(out_dir / "cad_mesh_case_study_frames.csv");
    std::ofstream comparison_out(out_dir / "cad_mesh_case_study_comparison.csv");
    std::ofstream geometry_out(out_dir / "cad_mesh_case_study_geometry.csv");
    frames_out << std::fixed << std::setprecision(8);
    comparison_out << std::fixed << std::setprecision(8);
    geometry_out << std::fixed << std::setprecision(8);

    frames_out << "scenario,variant,frame,time,center_x,center_y,center_z,velocity_x,velocity_y,velocity_z,"
                  "yaw,roll,contact_count,patch_count,newborn_count,merge_count,split_count,death_count,"
                  "max_source_count,slip_patch_count,total_force_x,total_force_y,total_force_z,total_force_norm,"
                  "total_torque_x,total_torque_y,total_torque_z,total_torque_norm,max_penetration,coulomb_ratio\n";
    comparison_out << "scenario,target_obj,moving_obj,target_vertices,target_faces,moving_vertices,moving_faces,"
                      "field_active_frames,native_active_frames,field_mean_patch_count,native_mean_contact_count,"
                      "field_mean_abs_count_change,native_mean_abs_contact_count_change,native_to_field_count_churn_ratio,"
                      "field_rms_normal_jump_angle,native_rms_contact_normal_jump_angle,"
                      "native_to_field_normal_jump_ratio,field_force_oscillation_index,native_force_oscillation_index,"
                      "native_to_field_force_oscillation_ratio,field_torque_oscillation_index,native_torque_oscillation_index,"
                      "native_to_field_torque_oscillation_ratio,field_max_penetration,native_max_penetration,"
                      "field_max_coulomb_ratio,native_max_coulomb_ratio,openvdb_active_voxels\n";
    geometry_out << "scenario,target_obj,moving_obj,target_vertices,target_faces,moving_vertices,moving_faces,"
                    "moving_scale,"
                    "target_min_x,target_min_y,target_min_z,target_max_x,target_max_y,target_max_z,"
                    "moving_original_min_x,moving_original_min_y,moving_original_min_z,"
                    "moving_original_max_x,moving_original_max_y,moving_original_max_z,"
                    "moving_local_min_x,moving_local_min_y,moving_local_min_z,"
                    "moving_local_max_x,moving_local_max_y,moving_local_max_z,"
                    "moving_graph_samples,voxel_size,openvdb_active_voxels\n";

    geometry_out << config.scenario << "," << config.target_obj.generic_string() << ","
                 << config.moving_obj.generic_string() << ","
                 << scene.target_mesh.vertices.size() << "," << scene.target_mesh.faces.size() << ","
                 << scene.moving_mesh_local.vertices.size() << "," << scene.moving_mesh_local.faces.size() << ","
                 << config.moving_scale << ","
                 << scene.target_bounds.min.x() << "," << scene.target_bounds.min.y() << "," << scene.target_bounds.min.z()
                 << "," << scene.target_bounds.max.x() << "," << scene.target_bounds.max.y() << ","
                 << scene.target_bounds.max.z() << ","
                 << scene.moving_bounds_original.min.x() << "," << scene.moving_bounds_original.min.y() << ","
                 << scene.moving_bounds_original.min.z() << "," << scene.moving_bounds_original.max.x() << ","
                 << scene.moving_bounds_original.max.y() << "," << scene.moving_bounds_original.max.z() << ","
                 << scene.moving_bounds_local.min.x() << "," << scene.moving_bounds_local.min.y() << ","
                 << scene.moving_bounds_local.min.z() << "," << scene.moving_bounds_local.max.x() << ","
                 << scene.moving_bounds_local.max.y() << "," << scene.moving_bounds_local.max.z() << ","
                 << scene.moving_graph.samples.size() << "," << config.voxel_size << ","
                 << scene.sdf.grid->activeVoxelCount() << "\n";

    double field_max_penetration = 0.0;
    FieldContactTopologyMetricsSummary field = RunFieldCase(config, scene, frames_out, field_max_penetration);
    NativeSummary native = RunNativeBaseline(config, scene, frames_out);

    comparison_out << config.scenario << "," << config.target_obj.generic_string() << ","
                   << config.moving_obj.generic_string() << ","
                   << scene.target_mesh.vertices.size() << "," << scene.target_mesh.faces.size() << ","
                   << scene.moving_mesh_local.vertices.size() << "," << scene.moving_mesh_local.faces.size() << ","
                   << field.active_frames << "," << native.active_frames << ","
                   << field.mean_patch_count << "," << native.mean_contact_count << ","
                   << field.mean_abs_patch_count_change << "," << native.mean_abs_contact_count_change << ","
                   << SafeRatio(native.mean_abs_contact_count_change, field.mean_abs_patch_count_change) << ","
                   << field.rms_normal_jump_angle << "," << native.rms_contact_normal_jump_angle << ","
                   << SafeRatio(native.rms_contact_normal_jump_angle, field.rms_normal_jump_angle) << ","
                   << field.force_oscillation_index << "," << native.force_oscillation_index << ","
                   << SafeRatio(native.force_oscillation_index, field.force_oscillation_index) << ","
                   << field.torque_oscillation_index << "," << native.torque_oscillation_index << ","
                   << SafeRatio(native.torque_oscillation_index, field.torque_oscillation_index) << ","
                   << field_max_penetration << "," << native.max_penetration << ","
                   << field.max_tangential_force_ratio << "," << native.max_coulomb_ratio << ","
                   << scene.sdf.grid->activeVoxelCount() << "\n";

    std::cout << "  target vertices/faces=" << scene.target_mesh.vertices.size() << "/" << scene.target_mesh.faces.size()
              << " moving vertices/faces=" << scene.moving_mesh_local.vertices.size() << "/"
              << scene.moving_mesh_local.faces.size() << std::endl;
    std::cout << "  field mean_patch=" << field.mean_patch_count
              << " normal_rms=" << field.rms_normal_jump_angle
              << " force_osc=" << field.force_oscillation_index
              << " coulomb=" << field.max_tangential_force_ratio
              << " max_pen=" << field_max_penetration << std::endl;
    std::cout << "  native mean_contacts=" << native.mean_contact_count
              << " normal_rms=" << native.rms_contact_normal_jump_angle
              << " force_osc=" << native.force_oscillation_index
              << " coulomb=" << native.max_coulomb_ratio
              << " max_pen=" << native.max_penetration << std::endl;
    std::cout << "Wrote:" << std::endl;
    std::cout << "  " << (out_dir / "cad_mesh_case_study_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "cad_mesh_case_study_comparison.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "cad_mesh_case_study_geometry.csv").string() << std::endl;
    return 0;
}
