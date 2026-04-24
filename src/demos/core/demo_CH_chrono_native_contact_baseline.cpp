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
// Milestone 23: Chrono native SMC contact baseline.
//
// This demo runs real Chrono native collision/contact for the same benchmark
// families used by the field-contact primitive experiments:
//   - cylinder over rail/groove;
//   - sphere through a nested/interlock slot;
//   - sphere over staggered multi-pad terrain.
//
// The moving body is driven along the same deterministic kinematic path used by
// the field primitive benchmark. Chrono's native SMC solver supplies contact
// count, force, and torque. The CSV also records commanded and post-step actual
// displacement/velocity/acceleration so the trajectory physicality can be
// checked when comparing against field primitive outputs.
//
// Output:
//   out/milestone_23/chrono_native_baseline_frames.csv
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "chrono/core/ChRotation.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"

using namespace chrono;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

enum class NativeScenario {
    GuideRailSliding,
    NestedInterlock,
    MultiPatchRollingSliding
};

enum class MovingBodyType {
    Sphere,
    CylinderZ
};

struct ScenarioConfig {
    NativeScenario scenario;
    std::string name;
    MovingBodyType moving_type = MovingBodyType::Sphere;
    double radius = 0.16;
    double cylinder_length = 0.30;
    double bottom_y = 0.0;
    ChVector3d start = ChVector3d(0, 0, 0);
    ChVector3d end = ChVector3d(0, 0, 0);
    double z_oscillation = 0.0;
    double z_oscillation_cycles = 1.0;
    double total_time = 1.4;
    int frames = 701;
    double rolling_scale = 1.0;
};

struct BodyKinematics {
    ChVector3d center = ChVector3d(0, 0, 0);
    ChVector3d velocity = ChVector3d(0, 0, 0);
    ChVector3d acceleration = ChVector3d(0, 0, 0);
    double roll_angle_x = 0.0;
    double roll_angle_z = 0.0;
    ChVector3d angular_velocity = ChVector3d(0, 0, 0);
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

static BodyKinematics KinematicsAtFrame(const ScenarioConfig& scenario, int frame) {
    double u = static_cast<double>(frame) / static_cast<double>(scenario.frames - 1);
    double du = 1.0 / static_cast<double>(scenario.frames - 1);
    double dt = scenario.total_time * du;
    ChVector3d start = scenario.start + ChVector3d(0, scenario.radius + scenario.bottom_y, 0);
    ChVector3d end = scenario.end + ChVector3d(0, scenario.radius + scenario.bottom_y, 0);

    auto center_at = [&](double uu) {
        double clamped = std::max(0.0, std::min(1.0, uu));
        double osc = scenario.z_oscillation * std::sin(2.0 * kPi * scenario.z_oscillation_cycles * clamped);
        return start + (end - start) * clamped + ChVector3d(0, 0, osc);
    };

    double u0 = std::max(0.0, u - du);
    double u1 = std::min(1.0, u + du);

    BodyKinematics kin;
    kin.center = center_at(u);
    kin.velocity = (center_at(u1) - center_at(u0)) / ((u1 - u0) * scenario.total_time);

    if (frame > 0 && frame + 1 < scenario.frames) {
        kin.acceleration = (center_at(u + du) - 2.0 * center_at(u) + center_at(u - du)) / (dt * dt);
    } else {
        kin.acceleration = ChVector3d(0, 0, 0);
    }

    kin.roll_angle_z = -scenario.rolling_scale * (kin.center.x() - start.x()) / scenario.radius;
    kin.roll_angle_x = scenario.rolling_scale * (kin.center.z() - start.z()) / scenario.radius;
    kin.angular_velocity = ChVector3d(kin.velocity.z() / scenario.radius,
                                      0.0,
                                      -kin.velocity.x() / scenario.radius) *
                           scenario.rolling_scale;
    return kin;
}

static std::shared_ptr<ChContactMaterialSMC> MakeMaterial() {
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetFriction(0.45f);
    mat->SetRestitution(0.02f);
    mat->SetYoungModulus(1.2e4f);
    mat->SetPoissonRatio(0.30f);
    return mat;
}

static void AddBox(ChSystemSMC& sys,
                   const std::shared_ptr<ChContactMaterialSMC>& mat,
                   double sx,
                   double sy,
                   double sz,
                   const ChVector3d& pos) {
    auto box = chrono_types::make_shared<ChBodyEasyBox>(sx, sy, sz, 1000.0, false, true, mat);
    box->SetPos(pos);
    box->SetFixed(true);
    box->EnableCollision(true);
    sys.AddBody(box);
}

static void AddGuideRailGeometry(ChSystemSMC& sys, const std::shared_ptr<ChContactMaterialSMC>& mat) {
    AddBox(sys, mat, 1.35, 0.018, 0.055, ChVector3d(0.0, -0.009, 0.0));
    AddBox(sys, mat, 1.35, 0.052, 0.026, ChVector3d(0.0, 0.026, 0.086));
    AddBox(sys, mat, 1.35, 0.052, 0.026, ChVector3d(0.0, 0.026, -0.086));
    for (int i = -3; i <= 3; i++) {
        double x = 0.16 * static_cast<double>(i);
        AddBox(sys, mat, 0.035, 0.020, 0.19, ChVector3d(x, 0.010, 0.0));
    }
}

static void AddNestedInterlockGeometry(ChSystemSMC& sys, const std::shared_ptr<ChContactMaterialSMC>& mat) {
    AddBox(sys, mat, 1.05, 0.018, 0.050, ChVector3d(0.0, -0.009, 0.0));
    AddBox(sys, mat, 1.00, 0.050, 0.026, ChVector3d(0.0, 0.025, 0.106));
    AddBox(sys, mat, 1.00, 0.050, 0.026, ChVector3d(0.0, 0.025, -0.106));
    AddBox(sys, mat, 0.080, 0.052, 0.170, ChVector3d(-0.115, 0.026, 0.0));
    AddBox(sys, mat, 0.080, 0.052, 0.170, ChVector3d(0.115, 0.026, 0.0));
    AddBox(sys, mat, 0.110, 0.036, 0.052, ChVector3d(0.0, 0.018, 0.060));
    AddBox(sys, mat, 0.110, 0.036, 0.052, ChVector3d(0.0, 0.018, -0.060));
}

static void AddMultiPadGeometry(ChSystemSMC& sys, const std::shared_ptr<ChContactMaterialSMC>& mat) {
    AddBox(sys, mat, 1.15, 0.012, 0.34, ChVector3d(0.0, -0.006, 0.0));
    for (int i = -4; i <= 4; i++) {
        double x = 0.105 * static_cast<double>(i);
        double z = (i % 2 == 0) ? 0.072 : -0.072;
        AddBox(sys, mat, 0.052, 0.050, 0.050, ChVector3d(x, 0.025, z));
        AddBox(sys, mat, 0.058, 0.042, 0.046, ChVector3d(x, 0.021, -0.62 * z));
    }
}

static std::shared_ptr<ChBody> AddMovingBody(ChSystemSMC& sys,
                                             const std::shared_ptr<ChContactMaterialSMC>& mat,
                                             const ScenarioConfig& scenario) {
    std::shared_ptr<ChBody> body;
    const double kinematic_probe_density = 1.0e7;
    if (scenario.moving_type == MovingBodyType::CylinderZ) {
        body = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Z,
                                                             scenario.radius,
                                                             scenario.cylinder_length,
                                                             kinematic_probe_density,
                                                             false,
                                                             true,
                                                             mat);
    } else {
        body = chrono_types::make_shared<ChBodyEasySphere>(scenario.radius, kinematic_probe_density, false, true, mat);
    }

    body->SetFixed(false);
    body->EnableCollision(true);
    sys.AddBody(body);
    return body;
}

static void ApplyCommandedKinematics(const std::shared_ptr<ChBody>& body,
                                     const BodyKinematics& kin,
                                     MovingBodyType body_type) {
    body->SetPos(kin.center);
    body->SetPosDt(kin.velocity);
    body->SetPosDt2(kin.acceleration);

    if (body_type == MovingBodyType::CylinderZ) {
        body->SetRot(QuatFromAngleZ(kin.roll_angle_z));
        body->SetAngVelLocal(ChVector3d(0, 0, kin.angular_velocity.z()));
    } else {
        body->SetRot(QuatFromAngleX(kin.roll_angle_x) * QuatFromAngleZ(kin.roll_angle_z));
        body->SetAngVelLocal(kin.angular_velocity);
    }
}

static void AddScenarioGeometry(ChSystemSMC& sys,
                                const std::shared_ptr<ChContactMaterialSMC>& mat,
                                NativeScenario scenario) {
    switch (scenario) {
        case NativeScenario::GuideRailSliding:
            AddGuideRailGeometry(sys, mat);
            break;
        case NativeScenario::NestedInterlock:
            AddNestedInterlockGeometry(sys, mat);
            break;
        case NativeScenario::MultiPatchRollingSliding:
            AddMultiPadGeometry(sys, mat);
            break;
    }
}

static void WriteFrameHeader(std::ofstream& out) {
    out << "scenario,frame,time,"
           "cmd_pos_x,cmd_pos_y,cmd_pos_z,cmd_vel_x,cmd_vel_y,cmd_vel_z,cmd_acc_x,cmd_acc_y,cmd_acc_z,"
           "actual_pos_x,actual_pos_y,actual_pos_z,actual_vel_x,actual_vel_y,actual_vel_z,"
           "actual_acc_x,actual_acc_y,actual_acc_z,pos_error_norm,vel_error_norm,acc_error_norm,"
           "contact_count,contact_force_x,contact_force_y,contact_force_z,contact_force_norm,"
           "contact_torque_x,contact_torque_y,contact_torque_z,contact_torque_norm\n";
}

static void RunScenario(const ScenarioConfig& scenario, std::ofstream& out) {
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(120);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-7);

    auto mat = MakeMaterial();
    ChCollisionInfo::SetDefaultEffectiveCurvatureRadius(scenario.radius);
    AddScenarioGeometry(sys, mat, scenario.scenario);
    auto moving = AddMovingBody(sys, mat, scenario);

    double dt = scenario.total_time / static_cast<double>(scenario.frames - 1);
    ChVector3d previous_actual_vel(0, 0, 0);
    bool has_previous_actual_vel = false;

    std::cout << "Scenario " << scenario.name << " (" << scenario.frames << " frames)" << std::endl;

    for (int frame = 0; frame < scenario.frames; frame++) {
        double time = static_cast<double>(frame) * dt;
        BodyKinematics kin = KinematicsAtFrame(scenario, frame);
        ApplyCommandedKinematics(moving, kin, scenario.moving_type);

        sys.DoStepDynamics(dt);

        ChVector3d actual_pos = moving->GetPos();
        ChVector3d actual_vel = moving->GetPosDt();
        ChVector3d actual_acc = has_previous_actual_vel ? (actual_vel - previous_actual_vel) / dt : kin.acceleration;
        previous_actual_vel = actual_vel;
        has_previous_actual_vel = true;

        ChVector3d contact_force = moving->GetContactForce();
        ChVector3d contact_torque = moving->GetContactTorque();
        int contact_count = static_cast<int>(sys.GetNumContacts());

        ChVector3d pos_error = actual_pos - kin.center;
        ChVector3d vel_error = actual_vel - kin.velocity;
        ChVector3d acc_error = actual_acc - kin.acceleration;

        out << scenario.name << "," << frame << "," << time << ","
            << kin.center.x() << "," << kin.center.y() << "," << kin.center.z() << ","
            << kin.velocity.x() << "," << kin.velocity.y() << "," << kin.velocity.z() << ","
            << kin.acceleration.x() << "," << kin.acceleration.y() << "," << kin.acceleration.z() << ","
            << actual_pos.x() << "," << actual_pos.y() << "," << actual_pos.z() << ","
            << actual_vel.x() << "," << actual_vel.y() << "," << actual_vel.z() << ","
            << actual_acc.x() << "," << actual_acc.y() << "," << actual_acc.z() << ","
            << pos_error.Length() << "," << vel_error.Length() << "," << acc_error.Length() << ","
            << contact_count << ","
            << contact_force.x() << "," << contact_force.y() << "," << contact_force.z() << ","
            << contact_force.Length() << ","
            << contact_torque.x() << "," << contact_torque.y() << "," << contact_torque.z() << ","
            << contact_torque.Length() << "\n";
    }
}

}  // namespace

int main() {
    std::cout << "Milestone 23: Chrono native SMC contact baseline" << std::endl;

    const std::string project_root = GetProjectRoot();
    const std::filesystem::path out_dir = std::filesystem::path(project_root) / "out" / "milestone_23";
    std::filesystem::create_directories(out_dir);

    std::ofstream out(out_dir / "chrono_native_baseline_frames.csv");
    out << std::fixed << std::setprecision(8);
    WriteFrameHeader(out);

    std::vector<ScenarioConfig> scenarios = {
        {NativeScenario::GuideRailSliding,
         "guide_rail_sliding",
         MovingBodyType::CylinderZ,
         0.16,
         0.30,
         0.012,
         ChVector3d(-0.52, 0, 0.0),
         ChVector3d(0.52, 0, 0.0),
         0.010,
         2.0,
         1.40,
         701,
         1.0},
        {NativeScenario::NestedInterlock,
         "nested_interlock",
         MovingBodyType::Sphere,
         0.16,
         0.30,
         0.000,
         ChVector3d(-0.42, 0, -0.012),
         ChVector3d(0.42, 0, 0.012),
         0.030,
         1.5,
         1.45,
         726,
         0.8},
        {NativeScenario::MultiPatchRollingSliding,
         "multi_patch_rolling_sliding",
         MovingBodyType::Sphere,
         0.16,
         0.30,
         0.016,
         ChVector3d(-0.48, 0, -0.085),
         ChVector3d(0.48, 0, 0.085),
         0.020,
         2.5,
         1.50,
         751,
         1.0},
    };

    for (const auto& scenario : scenarios) {
        RunScenario(scenario, out);
    }

    std::cout << "Wrote " << (out_dir / "chrono_native_baseline_frames.csv").string() << std::endl;
    return 0;
}
