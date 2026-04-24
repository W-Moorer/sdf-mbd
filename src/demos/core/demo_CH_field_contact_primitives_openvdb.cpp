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
// Field-based contact primitive implementation smoke test.
//
// This demo is deliberately smaller than the milestone research demos. Its
// purpose is to exercise the reusable implementation in
// chrono/collision/ChFieldContactPrimitives.h:
//   - surface graph active-set extraction
//   - connected field contact primitive construction
//   - sample-level area integral for normal force and torque
//   - overlap-based persistent primitive history sources
//   - objective tangential history transport and Coulomb projection
//
// =============================================================================

#define _USE_MATH_DEFINES

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactPrimitives.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetSphere.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

struct OpenVDBSphereSDF {
    openvdb::FloatGrid::Ptr grid;
    ChVector3d center = ChVector3d(0, 0, 0);
    double voxel_size = 0.02;

    FieldSampleQuery Query(const ChVector3d& world_pos, const ChVector3d& world_vel) const {
        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = world_vel;

        openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler> sampler(
            grid->tree(), grid->transform());

        ChVector3d local = world_pos - center;
        openvdb::Vec3d p(local.x(), local.y(), local.z());
        query.phi = static_cast<double>(sampler.wsSample(p));

        double h = 0.5 * voxel_size;
        double phi_px = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() + h, local.y(), local.z())));
        double phi_mx = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x() - h, local.y(), local.z())));
        double phi_py = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() + h, local.z())));
        double phi_my = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y() - h, local.z())));
        double phi_pz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() + h)));
        double phi_mz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(local.x(), local.y(), local.z() - h)));

        query.grad = ChVector3d((phi_px - phi_mx) / (2.0 * h),
                                (phi_py - phi_my) / (2.0 * h),
                                (phi_pz - phi_mz) / (2.0 * h));
        return query;
    }
};

static ChVector3d RotateZ(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z());
}

static std::string GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 3; i++) {
        path = path.parent_path();
    }
    return path.string();
}

static double SumHistoryWeights(const std::vector<HistorySource>& sources) {
    double sum = 0.0;
    for (const auto& source : sources) {
        sum += source.weight;
    }
    return sum;
}

}  // namespace

int main(int argc, char* argv[]) {
    std::cout << "=== Field-Based Contact Primitive Core Smoke Test ===" << std::endl;

    openvdb::initialize();

    const double target_radius = 1.0;
    const double moving_radius = 0.25;
    const double voxel_size = 0.01;
    const double dt = 0.002;
    const int frame_count = 420;

    OpenVDBSphereSDF sdf;
    sdf.center = ChVector3d(0, 0, 0);
    sdf.voxel_size = voxel_size;
    sdf.grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        static_cast<float>(target_radius),
        openvdb::Vec3f(0.0f, 0.0f, 0.0f),
        static_cast<float>(voxel_size),
        4.0f);
    sdf.grid->setGridClass(openvdb::GRID_LEVEL_SET);

    SurfaceGraph graph = MakeSphereSurfaceGraph(moving_radius, 48, 96);

    PatchExtractionSettings extract_settings;
    extract_settings.activation_band = 0.02;
    extract_settings.min_area = 1.0e-6;
    extract_settings.min_samples = 2;

    NormalContactSettings normal_settings;
    normal_settings.stiffness = 8.0e6;
    normal_settings.damping = 1.0e4;

    TangentialContactSettings tangential_settings;
    tangential_settings.stiffness = 2.0e5;
    tangential_settings.damping = 30.0;
    tangential_settings.friction_coefficient = 0.35;
    tangential_settings.time_step = dt;

    HistoryInheritanceSettings inheritance_settings;
    inheritance_settings.min_overlap = 0.02;
    inheritance_settings.min_normal_dot = 0.4;
    inheritance_settings.max_center_distance = 0.08;
    inheritance_settings.geometry_fallback_weight = 0.1;

    std::string out_dir = GetProjectRoot() + "/out/milestone_14";
    std::filesystem::create_directories(out_dir);
    std::ofstream csv(out_dir + "/field_contact_primitives_core.csv");
    csv << "frame,time,patch_count,active_area,total_force_x,total_force_y,total_force_z,"
        << "total_torque_x,total_torque_y,total_torque_z,main_persistent_id,main_area,"
        << "main_overlap_gate,main_normal_force,main_tangential_force,main_friction_limit,"
        << "main_stick_slip,history_energy_before_gate,history_energy_after_gate" << std::endl;

    std::vector<PrimitiveSnapshot> previous_snapshots;
    std::map<int, TangentialHistory> history_store;
    int next_persistent_id = 0;

    double max_tangent_ratio = 0.0;
    double min_overlap_gate = 1.0;
    int total_newborn = 0;
    int total_patch_frames = 0;

    for (int frame = 0; frame < frame_count; frame++) {
        double time = frame * dt;
        double roll_rate = 1.4;
        double angle = roll_rate * time;

        ChVector3d body_pos(0.08 * std::sin(0.9 * time),
                            target_radius + moving_radius - 0.014 + 0.004 * std::cos(0.6 * time),
                            0.04 * std::cos(0.7 * time));
        ChVector3d body_vel(0.08 * 0.9 * std::cos(0.9 * time),
                            -0.004 * 0.6 * std::sin(0.6 * time),
                            -0.04 * 0.7 * std::sin(0.7 * time));
        ChVector3d body_omega(0, 0, roll_rate);

        std::vector<FieldSampleQuery> queries(graph.samples.size());
        for (size_t si = 0; si < graph.samples.size(); si++) {
            ChVector3d local_rotated = RotateZ(graph.samples[si].local_pos, angle);
            ChVector3d world_pos = body_pos + local_rotated;
            ChVector3d world_vel = body_vel + body_omega.Cross(local_rotated);
            queries[si] = sdf.Query(world_pos, world_vel);
        }

        std::vector<int> active = BuildActiveSet(queries, extract_settings.activation_band);
        std::vector<PrimitivePatch> patches = ExtractPrimitives(graph, queries, active, extract_settings);

        ChVector3d total_force(0, 0, 0);
        ChVector3d total_torque(0, 0, 0);
        double active_area = 0.0;
        int main_persistent_id = -1;
        double main_area = 0.0;
        double main_overlap_gate = 0.0;
        double main_normal_force = 0.0;
        double main_tangential_force = 0.0;
        double main_friction_limit = 0.0;
        double main_energy_before = 0.0;
        double main_energy_after = 0.0;
        StickSlipState main_state = StickSlipState::Stick;

        std::vector<PrimitiveSnapshot> current_snapshots;
        current_snapshots.reserve(patches.size());

        for (auto& patch : patches) {
            ApplyNormalContactIntegral(patch, graph, queries, body_pos, normal_settings);

            std::vector<HistorySource> sources =
                ComputeHistorySources(patch, previous_snapshots, graph, inheritance_settings);
            double gate = SumHistoryWeights(sources);
            int persistent_id = -1;
            if (!sources.empty()) {
                persistent_id = sources.front().persistent_id;
            } else {
                persistent_id = next_persistent_id++;
                total_newborn++;
                gate = 0.0;
            }

            const TangentialHistory* previous_history = nullptr;
            auto hist_it = history_store.find(persistent_id);
            if (hist_it != history_store.end()) {
                previous_history = &hist_it->second;
            }

            ChVector3d vt = ProjectToTangent(patch.representative_velocity, patch.normal);
            TangentialUpdateResult tangential =
                UpdateTangentialContact(previous_history,
                                        patch.normal,
                                        vt,
                                        patch.normal_force.Length(),
                                        gate,
                                        tangential_settings);
            tangential.history.persistent_id = persistent_id;
            history_store[persistent_id] = tangential.history;

            patch.tangential_force = tangential.force;
            patch.force = patch.normal_force + patch.tangential_force;
            patch.torque += (patch.center - body_pos).Cross(patch.tangential_force);

            current_snapshots.push_back(MakeSnapshot(patch, persistent_id));

            total_force += patch.force;
            total_torque += patch.torque;
            active_area += patch.area;
            total_patch_frames++;

            double ratio = tangential.friction_limit > 1.0e-12 ?
                               tangential.final_force_norm / tangential.friction_limit :
                               0.0;
            max_tangent_ratio = std::max(max_tangent_ratio, ratio);
            if (!sources.empty()) {
                min_overlap_gate = std::min(min_overlap_gate, gate);
            }

            if (main_persistent_id < 0 || patch.normal_force.Length() > main_normal_force) {
                main_persistent_id = persistent_id;
                main_area = patch.area;
                main_overlap_gate = gate;
                main_normal_force = patch.normal_force.Length();
                main_tangential_force = tangential.final_force_norm;
                main_friction_limit = tangential.friction_limit;
                main_energy_before = tangential.stored_energy_before_gate;
                main_energy_after = tangential.stored_energy_after_gate;
                main_state = tangential.state;
            }
        }

        csv << frame << ","
            << std::fixed << std::setprecision(8)
            << time << ","
            << patches.size() << ","
            << active_area << ","
            << total_force.x() << ","
            << total_force.y() << ","
            << total_force.z() << ","
            << total_torque.x() << ","
            << total_torque.y() << ","
            << total_torque.z() << ","
            << main_persistent_id << ","
            << main_area << ","
            << main_overlap_gate << ","
            << main_normal_force << ","
            << main_tangential_force << ","
            << main_friction_limit << ","
            << (main_state == StickSlipState::Stick ? "stick" : "slip") << ","
            << main_energy_before << ","
            << main_energy_after << std::endl;

        previous_snapshots = current_snapshots;
    }

    csv.close();

    std::ofstream summary(out_dir + "/field_contact_primitives_core_summary.csv");
    summary << "metric,value" << std::endl;
    summary << std::fixed << std::setprecision(8)
            << "frames," << frame_count << std::endl
            << "surface_samples," << graph.samples.size() << std::endl
            << "total_patch_frames," << total_patch_frames << std::endl
            << "newborn_primitives," << total_newborn << std::endl
            << "max_tangential_force_ratio," << max_tangent_ratio << std::endl
            << "min_nonzero_overlap_gate," << (min_overlap_gate == 1.0 ? 0.0 : min_overlap_gate) << std::endl;
    summary.close();

    std::cout << "Output:" << std::endl;
    std::cout << "  " << out_dir + "/field_contact_primitives_core.csv" << std::endl;
    std::cout << "  " << out_dir + "/field_contact_primitives_core_summary.csv" << std::endl;
    std::cout << "Max tangential force / Coulomb limit: " << max_tangent_ratio << std::endl;
    std::cout << "Newborn primitives: " << total_newborn << std::endl;

    return 0;
}
