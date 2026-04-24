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
// Milestone 28: free dynamics response for field-contact primitives.
//
// Unlike the kinematic path benchmarks, the moving sphere here has real mass and
// inertia.  A velocity/preload actuator applies external forces, while contact is
// supplied either by Chrono native SMC mesh collision or by field-contact
// primitive forces accumulated directly on a ChBody.
//
// Outputs:
//   out/milestone_28/free_dynamics_frames.csv
//   out/milestone_28/free_dynamics_summary.csv
//   out/milestone_28/free_dynamics_comparison.csv
//   out/milestone_28/free_dynamics_geometry.csv
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
constexpr double kFrictionCoefficient = 0.45;

struct DemoTriangleMesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct OpenVDBSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 0.006;
};

enum class FreeDynamicsScenarioType {
    UChannel,
    StaggeredPadTrack
};

enum class DynamicsVariant {
    FieldPrimitiveAccumulator,
    NativeMeshSMC
};

struct FreeDynamicsScenario {
    FreeDynamicsScenarioType scenario;
    std::string name;
    ChVector3d initial_pos = ChVector3d(0, 0, 0);
    ChVector3d initial_vel = ChVector3d(0, 0, 0);
    double radius = 0.12;
    double density = 950.0;
    double dt = 5.0e-4;
    int steps = 2000;
    double drive_speed = 0.45;
    double drive_gain = 70.0;
    double max_drive_force = 18.0;
    double normal_preload = 28.0;
    double vertical_damping = 8.0;
    double z_reference_start = 0.0;
    double z_reference_end = 0.0;
    double z_oscillation = 0.0;
    double z_oscillation_cycles = 1.0;
    double lateral_gain = 65.0;
    double lateral_damping = 7.0;
    double max_lateral_force = 12.0;
    double voxel_size = 0.006;
    double native_mesh_radius = 0.0;
    std::filesystem::path target_obj;
    std::string source_demo;
};

struct SceneData {
    DemoTriangleMesh target_mesh;
    SurfaceGraph moving_graph;
    OpenVDBSDF sdf;
};

struct FrameSample {
    std::string scenario;
    std::string variant;
    int frame = 0;
    double time = 0.0;
    ChVector3d pos = ChVector3d(0, 0, 0);
    ChVector3d vel = ChVector3d(0, 0, 0);
    ChVector3d acc = ChVector3d(0, 0, 0);
    ChVector3d ang_vel = ChVector3d(0, 0, 0);
    ChVector3d external_force = ChVector3d(0, 0, 0);
    ChVector3d contact_force = ChVector3d(0, 0, 0);
    ChVector3d contact_torque = ChVector3d(0, 0, 0);
    int contact_count = 0;
    int patch_count = 0;
    int stick_count = 0;
    int slip_count = 0;
    double kinetic_energy = 0.0;
    double external_power = 0.0;
    double contact_power = 0.0;
    double external_work = 0.0;
    double contact_work = 0.0;
    double max_penetration = 0.0;
    double coulomb_ratio = 0.0;
};

struct DynamicsSummary {
    std::string scenario;
    std::string variant;
    int frames = 0;
    int active_frames = 0;
    double max_kinetic_energy = 0.0;
    double final_kinetic_energy = 0.0;
    double external_work = 0.0;
    double contact_work = 0.0;
    double positive_contact_work = 0.0;
    double energy_budget_ratio = 0.0;
    double max_penetration = 0.0;
    double max_coulomb_ratio = 0.0;
    double max_speed = 0.0;
    double rms_acceleration = 0.0;
    double rms_acceleration_jump = 0.0;
    double rms_lateral_velocity = 0.0;
    double final_x = 0.0;
    double final_y = 0.0;
    double final_z = 0.0;
    double mean_count = 0.0;
    double mean_abs_count_change = 0.0;
    int max_count = 0;
    int total_stick_frames = 0;
    int total_slip_frames = 0;
    int stick_slip_switches = 0;
};

struct NativeContactReport : public ChContactContainer::ReportContactCallback {
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
        (void)plane_coord;
        (void)eff_radius;
        (void)react_torques;
        (void)contactobjA;
        (void)contactobjB;
        (void)constraint_offset;

        count++;
        max_penetration = std::max(max_penetration, std::max(0.0, -distance));
        double fn = std::abs(react_forces.x());
        double ft = std::sqrt(react_forces.y() * react_forces.y() + react_forces.z() * react_forces.z());
        if (fn > 1.0e-12) {
            max_coulomb_ratio = std::max(max_coulomb_ratio, ft / (kFrictionCoefficient * fn));
            if (ft > 0.98 * kFrictionCoefficient * fn) {
                slip_like_count++;
            } else {
                stick_like_count++;
            }
        }
        return true;
    }

    int count = 0;
    int stick_like_count = 0;
    int slip_like_count = 0;
    double max_penetration = 0.0;
    double max_coulomb_ratio = 0.0;
};

static const char* VariantName(DynamicsVariant variant) {
    switch (variant) {
        case DynamicsVariant::FieldPrimitiveAccumulator:
            return "field_primitive_accumulator";
        case DynamicsVariant::NativeMeshSMC:
            return "native_mesh_smc";
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

static double Clamp(double value, double lo, double hi) {
    return std::max(lo, std::min(hi, value));
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
    if (mesh.vertices.empty() || mesh.faces.empty()) {
        throw std::runtime_error("OBJ mesh has no usable triangles: " + path.string());
    }
    return mesh;
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

static SceneData BuildScene(const FreeDynamicsScenario& config) {
    SceneData scene;
    scene.sdf.voxel_size = config.voxel_size;
    if (!config.target_obj.empty()) {
        scene.target_mesh = LoadObjTriangleMesh(config.target_obj);
    } else {
        scene.target_mesh = config.scenario == FreeDynamicsScenarioType::UChannel ? BuildUChannelTarget() :
                                                                                 BuildStaggeredPadTarget();
    }
    scene.moving_graph = MakeSphereSurfaceGraph(config.radius, 20, 40);
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

    settings.normal.stiffness = 8.0e5;
    settings.normal.damping = 4.0e4;

    settings.tangential.stiffness = 2.2e3;
    settings.tangential.damping = 12.0;
    settings.tangential.friction_coefficient = kFrictionCoefficient;
    settings.tangential.time_step = time_step;

    settings.inheritance.min_overlap = 0.005;
    settings.inheritance.min_normal_dot = 0.16;
    settings.inheritance.max_center_distance = 0.075;
    settings.inheritance.geometry_fallback_weight = 0.18;
    return settings;
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
                                                  const std::shared_ptr<ChBody>& body) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(graph.samples.size());
    SdfSampler sampler(sdf.grid->tree(), sdf.grid->transform());

    for (const auto& sample : graph.samples) {
        ChVector3d world_pos = body->TransformPointLocalToParent(sample.local_pos);
        ChVector3d world_vel = body->PointSpeedLocalToParent(sample.local_pos);
        queries.push_back(QueryOpenVDB(sdf, sampler, world_pos, world_vel));
    }
    return queries;
}

static std::shared_ptr<ChContactMaterialSMC> MakeMaterial() {
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetFriction(static_cast<float>(kFrictionCoefficient));
    mat->SetRestitution(0.02f);
    mat->SetYoungModulus(1.2e4f);
    mat->SetPoissonRatio(0.30f);
    return mat;
}

static void ConfigureSystem(ChSystemSMC& sys) {
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(120);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-7);
}

static void AddMatchedMeshTarget(ChSystemSMC& sys,
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
}

static std::shared_ptr<ChBody> AddMovingSphere(ChSystemSMC& sys,
                                               const std::shared_ptr<ChContactMaterialSMC>& mat,
                                               const FreeDynamicsScenario& scenario,
                                               bool native_collision) {
    auto body = chrono_types::make_shared<ChBodyEasySphere>(scenario.radius,
                                                            scenario.density,
                                                            false,
                                                            native_collision,
                                                            mat);
    body->SetPos(scenario.initial_pos);
    body->SetPosDt(scenario.initial_vel);
    body->SetFixed(false);
    body->SetSleepingAllowed(false);
    body->EnableCollision(native_collision);
    sys.AddBody(body);
    return body;
}

static double ZReference(const FreeDynamicsScenario& scenario, double time) {
    double total_time = scenario.dt * static_cast<double>(scenario.steps);
    double u = total_time > 0.0 ? Clamp(time / total_time, 0.0, 1.0) : 0.0;
    double z_line = scenario.z_reference_start + (scenario.z_reference_end - scenario.z_reference_start) * u;
    return z_line + scenario.z_oscillation * std::sin(2.0 * kPi * scenario.z_oscillation_cycles * u);
}

static ChVector3d ComputeExternalForce(const FreeDynamicsScenario& scenario,
                                       const std::shared_ptr<ChBody>& body,
                                       double time) {
    ChVector3d pos = body->GetPos();
    ChVector3d vel = body->GetPosDt();
    double fx = Clamp(scenario.drive_gain * (scenario.drive_speed - vel.x()),
                      -scenario.max_drive_force,
                      scenario.max_drive_force);
    double fy = -scenario.normal_preload - scenario.vertical_damping * vel.y();
    double z_ref = ZReference(scenario, time);
    double fz = Clamp(scenario.lateral_gain * (z_ref - pos.z()) - scenario.lateral_damping * vel.z(),
                      -scenario.max_lateral_force,
                      scenario.max_lateral_force);
    return ChVector3d(fx, fy, fz);
}

static double KineticEnergy(const std::shared_ptr<ChBody>& body) {
    ChVector3d v = body->GetPosDt();
    ChVector3d w = body->GetAngVelLocal();
    ChVector3d Iw = body->GetInertia() * w;
    return 0.5 * body->GetMass() * v.Dot(v) + 0.5 * w.Dot(Iw);
}

static double MaxPatchPenetration(const FieldContactStepResult& step) {
    double out = 0.0;
    for (const auto& patch : step.patches) {
        out = std::max(out, patch.patch.max_penetration);
    }
    return out;
}

static void CountStickSlip(const FieldContactStepResult& step, int& stick_count, int& slip_count) {
    stick_count = 0;
    slip_count = 0;
    for (const auto& patch : step.patches) {
        if (patch.tangential.state == StickSlipState::Slip) {
            slip_count++;
        } else {
            stick_count++;
        }
    }
}

static FrameSample MakeFrameSample(const FreeDynamicsScenario& scenario,
                                   DynamicsVariant variant,
                                   int frame,
                                   double time,
                                   const std::shared_ptr<ChBody>& body,
                                   const ChVector3d& previous_velocity,
                                   const ChVector3d& external_force,
                                   const ChVector3d& contact_force,
                                   const ChVector3d& contact_torque,
                                   int count,
                                   int patch_count,
                                   int stick_count,
                                   int slip_count,
                                   double max_penetration,
                                   double coulomb_ratio,
                                   double external_work,
                                   double contact_work,
                                   double dt) {
    FrameSample row;
    row.scenario = scenario.name;
    row.variant = VariantName(variant);
    row.frame = frame;
    row.time = time;
    row.pos = body->GetPos();
    row.vel = body->GetPosDt();
    row.acc = frame > 0 ? (row.vel - previous_velocity) / dt : ChVector3d(0, 0, 0);
    row.ang_vel = body->GetAngVelParent();
    row.external_force = external_force;
    row.contact_force = contact_force;
    row.contact_torque = contact_torque;
    row.contact_count = count;
    row.patch_count = patch_count;
    row.stick_count = stick_count;
    row.slip_count = slip_count;
    row.kinetic_energy = KineticEnergy(body);
    row.external_power = external_force.Dot(row.vel);
    row.contact_power = contact_force.Dot(row.vel) + contact_torque.Dot(row.ang_vel);
    row.external_work = external_work;
    row.contact_work = contact_work;
    row.max_penetration = max_penetration;
    row.coulomb_ratio = coulomb_ratio;
    return row;
}

static DynamicsSummary SummarizeFrames(const std::vector<FrameSample>& rows) {
    DynamicsSummary out;
    if (rows.empty()) {
        return out;
    }

    out.scenario = rows.front().scenario;
    out.variant = rows.front().variant;
    out.frames = static_cast<int>(rows.size());
    out.external_work = rows.back().external_work;
    out.contact_work = rows.back().contact_work;
    out.final_kinetic_energy = rows.back().kinetic_energy;
    out.final_x = rows.back().pos.x();
    out.final_y = rows.back().pos.y();
    out.final_z = rows.back().pos.z();

    double count_sum = 0.0;
    double count_change_sum = 0.0;
    int count_change_samples = 0;
    double acc_sq_sum = 0.0;
    double acc_jump_sq_sum = 0.0;
    int acc_jump_samples = 0;
    double lateral_vel_sq_sum = 0.0;

    for (size_t i = 0; i < rows.size(); i++) {
        const auto& row = rows[i];
        int count = row.variant == "field_primitive_accumulator" ? row.patch_count : row.contact_count;
        if (count > 0) {
            out.active_frames++;
        }
        count_sum += count;
        out.max_count = std::max(out.max_count, count);
        out.max_kinetic_energy = std::max(out.max_kinetic_energy, row.kinetic_energy);
        out.max_penetration = std::max(out.max_penetration, row.max_penetration);
        out.max_coulomb_ratio = std::max(out.max_coulomb_ratio, row.coulomb_ratio);
        out.max_speed = std::max(out.max_speed, row.vel.Length());
        acc_sq_sum += row.acc.Length2();
        lateral_vel_sq_sum += row.vel.z() * row.vel.z();
        if (row.contact_power > 0.0) {
            out.positive_contact_work += row.contact_power * (i > 0 ? (row.time - rows[i - 1].time) : 0.0);
        }
        if (row.stick_count > 0) {
            out.total_stick_frames++;
        }
        if (row.slip_count > 0) {
            out.total_slip_frames++;
        }
        if (i > 0) {
            const auto& prev = rows[i - 1];
            int prev_count = prev.variant == "field_primitive_accumulator" ? prev.patch_count : prev.contact_count;
            count_change_sum += std::abs(count - prev_count);
            count_change_samples++;
            ChVector3d da = row.acc - prev.acc;
            acc_jump_sq_sum += da.Length2();
            acc_jump_samples++;
            if ((prev.slip_count > 0) != (row.slip_count > 0)) {
                out.stick_slip_switches++;
            }
        }
    }

    out.mean_count = count_sum / static_cast<double>(rows.size());
    if (count_change_samples > 0) {
        out.mean_abs_count_change = count_change_sum / static_cast<double>(count_change_samples);
    }
    out.rms_acceleration = std::sqrt(acc_sq_sum / static_cast<double>(rows.size()));
    if (acc_jump_samples > 0) {
        out.rms_acceleration_jump = std::sqrt(acc_jump_sq_sum / static_cast<double>(acc_jump_samples));
    }
    out.rms_lateral_velocity = std::sqrt(lateral_vel_sq_sum / static_cast<double>(rows.size()));

    double energy_budget = std::max(1.0e-8, std::abs(out.external_work) + 1.0);
    out.energy_budget_ratio = out.max_kinetic_energy / energy_budget;
    return out;
}

static void WriteFrame(std::ofstream& out, const FrameSample& row) {
    out << row.scenario << "," << row.variant << "," << row.frame << "," << row.time << ","
        << row.pos.x() << "," << row.pos.y() << "," << row.pos.z() << ","
        << row.vel.x() << "," << row.vel.y() << "," << row.vel.z() << ","
        << row.acc.x() << "," << row.acc.y() << "," << row.acc.z() << ","
        << row.ang_vel.x() << "," << row.ang_vel.y() << "," << row.ang_vel.z() << ","
        << row.external_force.x() << "," << row.external_force.y() << "," << row.external_force.z() << ","
        << row.contact_force.x() << "," << row.contact_force.y() << "," << row.contact_force.z() << ","
        << row.contact_force.Length() << "," << row.contact_torque.x() << "," << row.contact_torque.y() << ","
        << row.contact_torque.z() << "," << row.contact_torque.Length() << ","
        << row.contact_count << "," << row.patch_count << "," << row.stick_count << "," << row.slip_count << ","
        << row.kinetic_energy << "," << row.external_power << "," << row.contact_power << ","
        << row.external_work << "," << row.contact_work << "," << row.max_penetration << ","
        << row.coulomb_ratio << "\n";
}

static void WriteSummary(std::ofstream& out, const DynamicsSummary& s) {
    out << s.scenario << "," << s.variant << "," << s.frames << "," << s.active_frames << ","
        << s.mean_count << "," << s.mean_abs_count_change << "," << s.max_count << ","
        << s.max_kinetic_energy << "," << s.final_kinetic_energy << "," << s.external_work << ","
        << s.contact_work << "," << s.positive_contact_work << "," << s.energy_budget_ratio << ","
        << s.max_penetration << "," << s.max_coulomb_ratio << "," << s.max_speed << ","
        << s.rms_acceleration << "," << s.rms_acceleration_jump << "," << s.rms_lateral_velocity << ","
        << s.final_x << "," << s.final_y << "," << s.final_z << ","
        << s.total_stick_frames << "," << s.total_slip_frames << "," << s.stick_slip_switches << "\n";
}

static DynamicsSummary RunFieldScenario(const FreeDynamicsScenario& scenario,
                                        const SceneData& scene,
                                        std::ofstream& frames_out) {
    ChSystemSMC sys;
    ConfigureSystem(sys);
    auto mat = MakeMaterial();
    auto body = AddMovingSphere(sys, mat, scenario, false);
    unsigned int accumulator = body->AddAccumulator();

    FieldContactRuntimeSettings settings = MakeRuntimeSettings(scenario.dt, 1.35 * scenario.voxel_size);
    FieldContactPrimitiveTracker tracker;
    tracker.Reset();

    std::vector<FrameSample> rows;
    rows.reserve(static_cast<size_t>(scenario.steps + 1));
    ChVector3d previous_velocity = body->GetPosDt();
    double external_work = 0.0;
    double contact_work = 0.0;

    for (int frame = 0; frame <= scenario.steps; frame++) {
        double time = static_cast<double>(frame) * scenario.dt;
        ChVector3d velocity_before = body->GetPosDt();
        ChVector3d angular_velocity_before = body->GetAngVelParent();
        ChVector3d external_force = ComputeExternalForce(scenario, body, time);

        auto queries = BuildQueries(scene.moving_graph, scene.sdf, body);
        FieldContactStepResult step = tracker.Evaluate(scene.moving_graph, queries, body->GetPos(), settings);
        int stick_count = 0;
        int slip_count = 0;
        CountStickSlip(step, stick_count, slip_count);

        double external_power = external_force.Dot(velocity_before);
        double contact_power = step.total_force.Dot(velocity_before) + step.total_torque.Dot(angular_velocity_before);
        if (frame > 0) {
            external_work += external_power * scenario.dt;
            contact_work += contact_power * scenario.dt;
        }

        body->EmptyAccumulator(accumulator);
        body->AccumulateForce(accumulator, external_force, body->GetPos(), false);
        body->AccumulateForce(accumulator, step.total_force, body->GetPos(), false);
        body->AccumulateTorque(accumulator, step.total_torque, false);
        sys.ForceUpdate();
        sys.DoStepDynamics(scenario.dt);

        FrameSample row = MakeFrameSample(scenario,
                                          DynamicsVariant::FieldPrimitiveAccumulator,
                                          frame,
                                          time,
                                          body,
                                          previous_velocity,
                                          external_force,
                                          step.total_force,
                                          step.total_torque,
                                          0,
                                          step.stats.patch_count,
                                          stick_count,
                                          slip_count,
                                          MaxPatchPenetration(step),
                                          step.stats.max_tangential_force_ratio,
                                          external_work,
                                          contact_work,
                                          scenario.dt);
        previous_velocity = row.vel;
        WriteFrame(frames_out, row);
        rows.push_back(row);
    }

    return SummarizeFrames(rows);
}

static DynamicsSummary RunNativeScenario(const FreeDynamicsScenario& scenario,
                                         const SceneData& scene,
                                         std::ofstream& frames_out) {
    ChSystemSMC sys;
    ConfigureSystem(sys);
    auto mat = MakeMaterial();
    ChCollisionInfo::SetDefaultEffectiveCurvatureRadius(scenario.radius);
    AddMatchedMeshTarget(sys, mat, scene.target_mesh, scenario.native_mesh_radius);
    auto body = AddMovingSphere(sys, mat, scenario, true);
    unsigned int accumulator = body->AddAccumulator();

    std::vector<FrameSample> rows;
    rows.reserve(static_cast<size_t>(scenario.steps + 1));
    ChVector3d previous_velocity = body->GetPosDt();
    double external_work = 0.0;
    double contact_work = 0.0;

    for (int frame = 0; frame <= scenario.steps; frame++) {
        double time = static_cast<double>(frame) * scenario.dt;
        ChVector3d velocity_before = body->GetPosDt();
        ChVector3d angular_velocity_before = body->GetAngVelParent();
        ChVector3d external_force = ComputeExternalForce(scenario, body, time);

        body->EmptyAccumulator(accumulator);
        body->AccumulateForce(accumulator, external_force, body->GetPos(), false);
        sys.ForceUpdate();
        sys.DoStepDynamics(scenario.dt);

        ChVector3d contact_force = body->GetContactForce();
        ChVector3d contact_torque = body->GetContactTorque();
        double external_power = external_force.Dot(velocity_before);
        double native_contact_power = contact_force.Dot(velocity_before) + contact_torque.Dot(angular_velocity_before);
        if (frame > 0) {
            external_work += external_power * scenario.dt;
            contact_work += native_contact_power * scenario.dt;
        }

        auto report = chrono_types::make_shared<NativeContactReport>();
        sys.GetContactContainer()->ReportAllContacts(report);

        FrameSample row = MakeFrameSample(scenario,
                                          DynamicsVariant::NativeMeshSMC,
                                          frame,
                                          time,
                                          body,
                                          previous_velocity,
                                          external_force,
                                          contact_force,
                                          contact_torque,
                                          report->count,
                                          0,
                                          report->stick_like_count,
                                          report->slip_like_count,
                                          report->max_penetration,
                                          report->max_coulomb_ratio,
                                          external_work,
                                          contact_work,
                                          scenario.dt);
        previous_velocity = row.vel;
        WriteFrame(frames_out, row);
        rows.push_back(row);
    }

    return SummarizeFrames(rows);
}

static double SafeRatio(double numerator, double denominator) {
    return denominator > 1.0e-14 ? numerator / denominator : 0.0;
}

static void WriteComparisonHeader(std::ofstream& out) {
    out << "scenario,"
           "field_max_kinetic_energy,native_max_kinetic_energy,"
           "field_energy_budget_ratio,native_energy_budget_ratio,"
           "field_contact_work,native_contact_work,field_positive_contact_work,native_positive_contact_work,"
           "field_max_penetration,native_max_penetration,"
           "field_max_coulomb_ratio,native_max_coulomb_ratio,"
           "field_rms_acceleration_jump,native_rms_acceleration_jump,"
           "field_to_native_acceleration_jump_ratio,"
           "field_mean_count,native_mean_count,field_mean_abs_count_change,native_mean_abs_count_change,"
           "field_final_x,native_final_x,field_final_y,native_final_y,field_final_z,native_final_z,"
           "trajectory_final_position_delta,field_stick_slip_switches,native_stick_slip_switches\n";
}

static void WriteComparisonRow(std::ofstream& out,
                               const std::string& scenario,
                               const DynamicsSummary& field,
                               const DynamicsSummary& native) {
    ChVector3d field_final(field.final_x, field.final_y, field.final_z);
    ChVector3d native_final(native.final_x, native.final_y, native.final_z);
    out << scenario << ","
        << field.max_kinetic_energy << "," << native.max_kinetic_energy << ","
        << field.energy_budget_ratio << "," << native.energy_budget_ratio << ","
        << field.contact_work << "," << native.contact_work << ","
        << field.positive_contact_work << "," << native.positive_contact_work << ","
        << field.max_penetration << "," << native.max_penetration << ","
        << field.max_coulomb_ratio << "," << native.max_coulomb_ratio << ","
        << field.rms_acceleration_jump << "," << native.rms_acceleration_jump << ","
        << SafeRatio(field.rms_acceleration_jump, native.rms_acceleration_jump) << ","
        << field.mean_count << "," << native.mean_count << ","
        << field.mean_abs_count_change << "," << native.mean_abs_count_change << ","
        << field.final_x << "," << native.final_x << ","
        << field.final_y << "," << native.final_y << ","
        << field.final_z << "," << native.final_z << ","
        << (field_final - native_final).Length() << ","
        << field.stick_slip_switches << "," << native.stick_slip_switches << "\n";
}

}  // namespace

int main() {
    std::cout << "Milestone 28: free dynamics field-contact response" << std::endl;
    openvdb::initialize();

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_28";
    std::filesystem::create_directories(out_dir);

    std::ofstream frames_out(out_dir / "free_dynamics_frames.csv");
    std::ofstream summary_out(out_dir / "free_dynamics_summary.csv");
    std::ofstream comparison_out(out_dir / "free_dynamics_comparison.csv");
    std::ofstream geometry_out(out_dir / "free_dynamics_geometry.csv");

    frames_out << std::fixed << std::setprecision(8);
    summary_out << std::fixed << std::setprecision(8);
    comparison_out << std::fixed << std::setprecision(8);
    geometry_out << std::fixed << std::setprecision(8);

    frames_out << "scenario,variant,frame,time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z,"
                  "ang_vel_x,ang_vel_y,ang_vel_z,external_force_x,external_force_y,external_force_z,"
                  "contact_force_x,contact_force_y,contact_force_z,contact_force_norm,"
                  "contact_torque_x,contact_torque_y,contact_torque_z,contact_torque_norm,"
                  "contact_count,patch_count,stick_count,slip_count,kinetic_energy,external_power,contact_power,"
                  "external_work,contact_work,max_penetration,coulomb_ratio\n";
    summary_out << "scenario,variant,frames,active_frames,mean_count,mean_abs_count_change,max_count,"
                   "max_kinetic_energy,final_kinetic_energy,external_work,contact_work,positive_contact_work,"
                   "energy_budget_ratio,max_penetration,max_coulomb_ratio,max_speed,rms_acceleration,"
                   "rms_acceleration_jump,rms_lateral_velocity,final_x,final_y,final_z,total_stick_frames,"
                   "total_slip_frames,stick_slip_switches\n";
    WriteComparisonHeader(comparison_out);
    geometry_out << "scenario,source_demo,target_obj,target_vertices,target_faces,moving_samples,voxel_size,"
                    "active_voxels,native_mesh_radius\n";

    const auto data_dir = std::filesystem::path(project_root) / "data" / "models";

    std::vector<FreeDynamicsScenario> scenarios = {
        {FreeDynamicsScenarioType::UChannel,
         "free_u_channel_sphere_slide",
         ChVector3d(-0.54, 0.116, 0.000),
         ChVector3d(0.18, 0.0, 0.0),
         0.120,
         950.0,
         5.0e-4,
         2400,
         0.42,
         72.0,
         16.0,
         30.0,
         8.0,
         0.000,
         0.000,
         0.018,
         1.0,
         60.0,
         6.0,
         10.0,
         0.006,
         0.0},
        {FreeDynamicsScenarioType::StaggeredPadTrack,
         "free_staggered_pad_track",
         ChVector3d(-0.54, 0.145, -0.074),
         ChVector3d(0.16, 0.0, 0.0),
         0.120,
         950.0,
         5.0e-4,
         2500,
         0.38,
         65.0,
         14.0,
         16.0,
         7.0,
         -0.074,
         0.074,
         0.014,
         2.0,
         55.0,
         6.0,
         10.0,
         0.006,
         0.0},
        {FreeDynamicsScenarioType::StaggeredPadTrack,
         "free_chrono_original_shoe_collision_slide",
         ChVector3d(0.050, 0.128, -0.055),
         ChVector3d(0.060, 0.0, 0.0),
         0.055,
         950.0,
         5.0e-4,
         1800,
         0.120,
         45.0,
         8.0,
         16.0,
         5.0,
         -0.055,
         0.055,
         0.012,
         1.0,
         35.0,
         4.0,
         7.0,
         0.002,
         0.0,
         data_dir / "bulldozer" / "shoe_collision.obj",
         "src/demos/mbs/demo_MBS_tracks.cpp"},
    };

    for (const auto& scenario : scenarios) {
        std::cout << "Scenario " << scenario.name << " (" << scenario.steps << " steps)" << std::endl;
        SceneData scene = BuildScene(scenario);
        geometry_out << scenario.name << "," << scenario.source_demo << "," << scenario.target_obj.generic_string() << ","
                     << scene.target_mesh.vertices.size() << "," << scene.target_mesh.faces.size() << ","
                     << scene.moving_graph.samples.size() << "," << scene.sdf.voxel_size << ","
                     << scene.sdf.grid->activeVoxelCount() << "," << scenario.native_mesh_radius << "\n";

        DynamicsSummary field = RunFieldScenario(scenario, scene, frames_out);
        DynamicsSummary native = RunNativeScenario(scenario, scene, frames_out);
        WriteSummary(summary_out, field);
        WriteSummary(summary_out, native);
        WriteComparisonRow(comparison_out, scenario.name, field, native);

        std::cout << "  field: maxK=" << field.max_kinetic_energy
                  << " max_pen=" << field.max_penetration
                  << " coulomb=" << field.max_coulomb_ratio
                  << " acc_jump=" << field.rms_acceleration_jump << std::endl;
        std::cout << "  native: maxK=" << native.max_kinetic_energy
                  << " max_pen=" << native.max_penetration
                  << " coulomb=" << native.max_coulomb_ratio
                  << " acc_jump=" << native.rms_acceleration_jump << std::endl;
    }

    std::cout << "Wrote:" << std::endl;
    std::cout << "  " << (out_dir / "free_dynamics_frames.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "free_dynamics_summary.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "free_dynamics_comparison.csv").string() << std::endl;
    std::cout << "  " << (out_dir / "free_dynamics_geometry.csv").string() << std::endl;
    return 0;
}
