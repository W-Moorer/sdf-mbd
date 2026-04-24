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
// Milestone 15: real split/merge field-contact primitive demo.
//
// A triangulated surface plane moves against an OpenVDB SDF generated from two
// overlapping sphere meshes. As the plane height changes, the active surface
// region naturally evolves from two disconnected patches to one merged patch
// and then splits again. The demo records overlap-based history sources,
// split/merge classifications, tangential-history energy, and force jumps.
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
#include <vector>

#include "chrono/collision/ChFieldContactPrimitives.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

struct DemoTriangleMesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct OpenVDBSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 0.01;

    FieldSampleQuery Query(const ChVector3d& world_pos, const ChVector3d& world_vel) const {
        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = world_vel;

        openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler> sampler(
            grid->tree(), grid->transform());

        openvdb::Vec3d p(world_pos.x(), world_pos.y(), world_pos.z());
        query.phi = static_cast<double>(sampler.wsSample(p));

        double h = 0.5 * voxel_size;
        double phi_px = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x() + h, world_pos.y(), world_pos.z())));
        double phi_mx = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x() - h, world_pos.y(), world_pos.z())));
        double phi_py = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y() + h, world_pos.z())));
        double phi_my = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y() - h, world_pos.z())));
        double phi_pz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y(), world_pos.z() + h)));
        double phi_mz = static_cast<double>(sampler.wsSample(openvdb::Vec3d(world_pos.x(), world_pos.y(), world_pos.z() - h)));

        query.grad = ChVector3d((phi_px - phi_mx) / (2.0 * h),
                                (phi_py - phi_my) / (2.0 * h),
                                (phi_pz - phi_mz) / (2.0 * h));
        return query;
    }
};

static DemoTriangleMesh BuildTriangulatedSphereMesh(double radius,
                                                    int n_latitude,
                                                    int n_longitude,
                                                    const ChVector3d& offset) {
    DemoTriangleMesh mesh;
    if (radius <= 0.0 || n_latitude < 3 || n_longitude < 3) {
        return mesh;
    }

    int north = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(offset + ChVector3d(0, radius, 0));

    auto ring_index = [n_longitude](int ring, int j) {
        return 1 + (ring - 1) * n_longitude + ((j + n_longitude) % n_longitude);
    };

    for (int i = 1; i < n_latitude; i++) {
        double theta = M_PI * static_cast<double>(i) / static_cast<double>(n_latitude);
        for (int j = 0; j < n_longitude; j++) {
            double phi = 2.0 * M_PI * static_cast<double>(j) / static_cast<double>(n_longitude);
            mesh.vertices.push_back(offset + ChVector3d(radius * std::sin(theta) * std::cos(phi),
                                                        radius * std::cos(theta),
                                                        radius * std::sin(theta) * std::sin(phi)));
        }
    }

    int south = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(offset + ChVector3d(0, -radius, 0));

    for (int j = 0; j < n_longitude; j++) {
        mesh.faces.push_back({north, ring_index(1, j + 1), ring_index(1, j)});
    }

    for (int i = 1; i < n_latitude - 1; i++) {
        for (int j = 0; j < n_longitude; j++) {
            int a = ring_index(i, j);
            int b = ring_index(i, j + 1);
            int c = ring_index(i + 1, j);
            int d = ring_index(i + 1, j + 1);
            mesh.faces.push_back({a, b, d});
            mesh.faces.push_back({a, d, c});
        }
    }

    for (int j = 0; j < n_longitude; j++) {
        mesh.faces.push_back({ring_index(n_latitude - 1, j), ring_index(n_latitude - 1, j + 1), south});
    }

    return mesh;
}

static void AppendMesh(DemoTriangleMesh& dst, const DemoTriangleMesh& src) {
    int vertex_offset = static_cast<int>(dst.vertices.size());
    dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
    for (const auto& face : src.faces) {
        dst.faces.push_back({face.v0 + vertex_offset, face.v1 + vertex_offset, face.v2 + vertex_offset});
    }
}

static DemoTriangleMesh BuildDoubleSphereTarget(double radius, double separation) {
    DemoTriangleMesh target;
    AppendMesh(target, BuildTriangulatedSphereMesh(radius, 56, 112, ChVector3d(-0.5 * separation, 0, 0)));
    AppendMesh(target, BuildTriangulatedSphereMesh(radius, 56, 112, ChVector3d(0.5 * separation, 0, 0)));
    return target;
}

static DemoTriangleMesh BuildPlaneMesh(double half_x, double half_z, int nx, int nz) {
    DemoTriangleMesh mesh;
    if (half_x <= 0.0 || half_z <= 0.0 || nx < 1 || nz < 1) {
        return mesh;
    }

    mesh.vertices.reserve(static_cast<size_t>((nx + 1) * (nz + 1)));
    for (int iz = 0; iz <= nz; iz++) {
        double z = -half_z + 2.0 * half_z * static_cast<double>(iz) / static_cast<double>(nz);
        for (int ix = 0; ix <= nx; ix++) {
            double x = -half_x + 2.0 * half_x * static_cast<double>(ix) / static_cast<double>(nx);
            mesh.vertices.push_back(ChVector3d(x, 0, z));
        }
    }

    auto vid = [nx](int ix, int iz) {
        return iz * (nx + 1) + ix;
    };

    mesh.faces.reserve(static_cast<size_t>(2 * nx * nz));
    for (int iz = 0; iz < nz; iz++) {
        for (int ix = 0; ix < nx; ix++) {
            int v00 = vid(ix, iz);
            int v10 = vid(ix + 1, iz);
            int v01 = vid(ix, iz + 1);
            int v11 = vid(ix + 1, iz + 1);
            mesh.faces.push_back({v00, v10, v11});
            mesh.faces.push_back({v00, v11, v01});
        }
    }

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

static std::string GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 3; i++) {
        path = path.parent_path();
    }
    return path.string();
}

static std::string JoinSourceIds(const std::vector<HistorySource>& sources) {
    std::ostringstream ss;
    for (size_t i = 0; i < sources.size(); i++) {
        if (i > 0) {
            ss << ";";
        }
        ss << sources[i].persistent_id << ":" << std::fixed << std::setprecision(4) << sources[i].weight;
    }
    return ss.str();
}

static double PlaneHeight(double phase, double y_high, double y_low) {
    if (phase < 0.20) {
        return y_high;
    }
    if (phase < 0.45) {
        double u = (phase - 0.20) / 0.25;
        double s = 0.5 - 0.5 * std::cos(M_PI * u);
        return y_high + (y_low - y_high) * s;
    }
    if (phase < 0.65) {
        return y_low;
    }
    if (phase < 0.90) {
        double u = (phase - 0.65) / 0.25;
        double s = 0.5 - 0.5 * std::cos(M_PI * u);
        return y_low + (y_high - y_low) * s;
    }
    return y_high;
}

static double PlaneHeightDt(double phase, double total_time, double y_high, double y_low) {
    const double eps = 1.0e-4;
    double p0 = std::max(0.0, phase - eps);
    double p1 = std::min(1.0, phase + eps);
    return (PlaneHeight(p1, y_high, y_low) - PlaneHeight(p0, y_high, y_low)) /
           ((p1 - p0) * total_time);
}

}  // namespace

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 15: Field Contact Primitive Split/Merge ===" << std::endl;

    openvdb::initialize();

    const double target_radius = 0.45;
    const double target_separation = 0.56;
    const double voxel_size = 0.01;
    const double y_high = 0.431;
    const double y_low = 0.305;
    const double total_time = 1.6;
    const double dt = 0.002;
    const int frame_count = static_cast<int>(total_time / dt);

    DemoTriangleMesh target_mesh = BuildDoubleSphereTarget(target_radius, target_separation);
    DemoTriangleMesh plane_mesh = BuildPlaneMesh(0.85, 0.45, 95, 55);
    SurfaceGraph plane_graph = MakeTriangleMeshSurfaceGraph(plane_mesh.vertices, plane_mesh.faces);

    OpenVDBSDF target_sdf;
    target_sdf.voxel_size = voxel_size;
    target_sdf.grid = BuildLevelSetFromTriangleMesh(target_mesh, voxel_size, 4.0f);

    PatchExtractionSettings extract_settings;
    extract_settings.activation_band = 0.012;
    extract_settings.min_area = 1.0e-5;
    extract_settings.min_samples = 3;

    NormalContactSettings normal_settings;
    normal_settings.stiffness = 5.0e6;
    normal_settings.damping = 1.0e4;

    TangentialContactSettings tangential_settings;
    tangential_settings.stiffness = 1.5e5;
    tangential_settings.damping = 20.0;
    tangential_settings.friction_coefficient = 0.35;
    tangential_settings.time_step = dt;

    HistoryInheritanceSettings inheritance_settings;
    inheritance_settings.min_overlap = 0.01;
    inheritance_settings.min_normal_dot = 0.35;
    inheritance_settings.max_center_distance = 0.18;
    inheritance_settings.geometry_fallback_weight = 0.05;

    std::string out_dir = GetProjectRoot() + "/out/milestone_15";
    std::filesystem::create_directories(out_dir);

    std::ofstream frame_csv(out_dir + "/field_contact_split_merge_frames.csv");
    frame_csv << "frame,time,plane_y,patch_count,newborn_count,merge_count,split_count,death_count,"
              << "max_source_count,max_previous_reuse,total_force_x,total_force_y,total_force_z,"
              << "total_torque_x,total_torque_y,total_torque_z,force_jump,torque_jump,max_tangent_ratio,"
              << "max_energy_gate_ratio,max_inherited_energy_ratio" << std::endl;

    std::ofstream patch_csv(out_dir + "/field_contact_split_merge_patches.csv");
    patch_csv << "frame,time,patch_index,persistent_id,event,area,center_x,center_y,center_z,"
              << "source_count,source_weight_sum,sources,normal_force,tangential_force,friction_limit,"
              << "tangent_ratio,energy_before_gate,energy_after_gate,energy_gate_ratio,"
              << "source_energy_bound,inherited_energy,inherited_energy_ratio,stick_slip" << std::endl;

    std::vector<PrimitiveSnapshot> previous_snapshots;
    std::map<int, TangentialHistory> history_store;
    int next_persistent_id = 0;

    int frames_with_two = 0;
    int frames_with_one = 0;
    int frames_with_merge = 0;
    int frames_with_split = 0;
    int total_newborn = 0;
    int total_death = 0;
    int total_merge_patches = 0;
    int total_split_patches = 0;
    int max_patch_count = 0;
    int max_source_count_all = 0;
    int max_previous_reuse_all = 0;
    double max_tangent_ratio_all = 0.0;
    double max_energy_gate_ratio_all = 0.0;
    double max_inherited_energy_ratio_all = 0.0;
    double max_force_jump = 0.0;
    double max_torque_jump = 0.0;
    ChVector3d previous_total_force(0, 0, 0);
    ChVector3d previous_total_torque(0, 0, 0);

    for (int frame = 0; frame < frame_count; frame++) {
        double time = frame * dt;
        double phase = time / total_time;
        double y = PlaneHeight(phase, y_high, y_low);
        double ydot = PlaneHeightDt(phase, total_time, y_high, y_low);
        double x = 0.035 * std::sin(2.0 * M_PI * phase);
        double xdot = 0.035 * (2.0 * M_PI / total_time) * std::cos(2.0 * M_PI * phase);
        ChVector3d body_pos(x, y, 0);
        ChVector3d body_vel(xdot, ydot, 0);

        std::vector<FieldSampleQuery> queries(plane_graph.samples.size());
        for (size_t si = 0; si < plane_graph.samples.size(); si++) {
            ChVector3d world_pos = body_pos + plane_graph.samples[si].local_pos;
            queries[si] = target_sdf.Query(world_pos, body_vel);
        }

        std::vector<int> active = BuildActiveSet(queries, extract_settings.activation_band);
        std::vector<PrimitivePatch> patches = ExtractPrimitives(plane_graph, queries, active, extract_settings);

        struct PatchCandidate {
            std::vector<HistorySource> sources;
            int persistent_id = -1;
            std::string event = "stable";
        };

        std::vector<PatchCandidate> candidates(patches.size());
        std::map<int, std::vector<int>> primary_users;
        std::set<int> referenced_previous;

        for (size_t pi = 0; pi < patches.size(); pi++) {
            candidates[pi].sources = ComputeHistorySources(patches[pi], previous_snapshots, plane_graph, inheritance_settings);
            if (!candidates[pi].sources.empty()) {
                int primary = candidates[pi].sources.front().persistent_id;
                primary_users[primary].push_back(static_cast<int>(pi));
                for (const auto& source : candidates[pi].sources) {
                    referenced_previous.insert(source.persistent_id);
                }
            }
        }

        for (auto& entry : primary_users) {
            auto& users = entry.second;
            std::sort(users.begin(), users.end(), [&](int a, int b) {
                double wa = candidates[a].sources.empty() ? 0.0 : candidates[a].sources.front().weight;
                double wb = candidates[b].sources.empty() ? 0.0 : candidates[b].sources.front().weight;
                return wa > wb;
            });
        }

        int newborn_count = 0;
        int merge_count = 0;
        int split_count = 0;
        int max_source_count = 0;
        int max_previous_reuse = 0;
        for (const auto& entry : primary_users) {
            max_previous_reuse = std::max(max_previous_reuse, static_cast<int>(entry.second.size()));
        }

        for (size_t pi = 0; pi < patches.size(); pi++) {
            auto& candidate = candidates[pi];
            max_source_count = std::max(max_source_count, static_cast<int>(candidate.sources.size()));

            if (candidate.sources.empty()) {
                candidate.persistent_id = next_persistent_id++;
                candidate.event = "newborn";
                newborn_count++;
                total_newborn++;
                continue;
            }

            int primary = candidate.sources.front().persistent_id;
            const auto& users = primary_users[primary];
            bool split_child = users.size() > 1 && users.front() != static_cast<int>(pi);
            bool merge_patch = candidate.sources.size() > 1;

            if (split_child) {
                candidate.persistent_id = next_persistent_id++;
                candidate.event = merge_patch ? "split_merge" : "split";
                split_count++;
                total_split_patches++;
            } else {
                candidate.persistent_id = primary;
                if (merge_patch) {
                    candidate.event = "merge";
                    merge_count++;
                    total_merge_patches++;
                } else if (users.size() > 1) {
                    candidate.event = "split_primary";
                    split_count++;
                    total_split_patches++;
                } else {
                    candidate.event = "stable";
                }
            }
        }

        int death_count = 0;
        for (const auto& prev : previous_snapshots) {
            if (!referenced_previous.count(prev.persistent_id)) {
                death_count++;
            }
        }
        total_death += death_count;

        std::vector<PrimitiveSnapshot> current_snapshots;
        std::map<int, TangentialHistory> new_history_store;
        ChVector3d total_force(0, 0, 0);
        ChVector3d total_torque(0, 0, 0);
        double max_tangent_ratio = 0.0;
        double max_energy_gate_ratio = 0.0;
        double max_inherited_energy_ratio = 0.0;

        for (size_t pi = 0; pi < patches.size(); pi++) {
            auto& patch = patches[pi];
            auto& candidate = candidates[pi];

            ApplyNormalContactIntegral(patch, plane_graph, queries, body_pos, normal_settings);

            std::vector<WeightedTangentialHistorySource> weighted_history_sources;
            double source_weight_sum = 0.0;
            double source_energy_bound = 0.0;
            for (const auto& source : candidate.sources) {
                source_weight_sum += source.weight;
                auto hist_it = history_store.find(source.persistent_id);
                if (hist_it != history_store.end() && hist_it->second.valid) {
                    WeightedTangentialHistorySource weighted_source;
                    weighted_source.history = hist_it->second;
                    weighted_source.weight = source.weight;
                    weighted_history_sources.push_back(weighted_source);

                    ChVector3d source_xi = TransportElasticStateMinimalRotation(hist_it->second.xi_elastic_world,
                                                                                hist_it->second.normal,
                                                                                patch.normal);
                    source_energy_bound += source.weight * 0.5 * tangential_settings.stiffness * source_xi.Dot(source_xi);
                }
            }

            TangentialHistory aggregated_history =
                AggregateTangentialHistorySources(weighted_history_sources, patch.normal, candidate.persistent_id);
            const TangentialHistory* previous_history = aggregated_history.valid ? &aggregated_history : nullptr;
            double inherited_energy = 0.5 * tangential_settings.stiffness *
                                      aggregated_history.xi_elastic_world.Dot(aggregated_history.xi_elastic_world);
            double inherited_energy_ratio = source_energy_bound > 1.0e-16 ?
                                                inherited_energy / source_energy_bound :
                                                0.0;

            ChVector3d vt = ProjectToTangent(patch.representative_velocity, patch.normal);
            TangentialUpdateResult tangential =
                UpdateTangentialContact(previous_history,
                                        patch.normal,
                                        vt,
                                        patch.normal_force.Length(),
                                        previous_history ? 1.0 : 0.0,
                                        tangential_settings);
            tangential.history.persistent_id = candidate.persistent_id;
            new_history_store[candidate.persistent_id] = tangential.history;

            patch.tangential_force = tangential.force;
            patch.force = patch.normal_force + patch.tangential_force;
            patch.torque += (patch.center - body_pos).Cross(patch.tangential_force);

            total_force += patch.force;
            total_torque += patch.torque;
            current_snapshots.push_back(MakeSnapshot(patch, candidate.persistent_id));

            double tangent_ratio = tangential.friction_limit > 1.0e-12 ?
                                       tangential.final_force_norm / tangential.friction_limit :
                                       0.0;
            double energy_ratio = tangential.stored_energy_before_gate > 1.0e-16 ?
                                      tangential.stored_energy_after_gate / tangential.stored_energy_before_gate :
                                      0.0;
            max_tangent_ratio = std::max(max_tangent_ratio, tangent_ratio);
            max_energy_gate_ratio = std::max(max_energy_gate_ratio, energy_ratio);
            max_inherited_energy_ratio = std::max(max_inherited_energy_ratio, inherited_energy_ratio);

            patch_csv << frame << ","
                      << std::fixed << std::setprecision(8)
                      << time << ","
                      << pi << ","
                      << candidate.persistent_id << ","
                      << candidate.event << ","
                      << patch.area << ","
                      << patch.center.x() << ","
                      << patch.center.y() << ","
                      << patch.center.z() << ","
                      << candidate.sources.size() << ","
                      << source_weight_sum << ","
                      << JoinSourceIds(candidate.sources) << ","
                      << patch.normal_force.Length() << ","
                      << tangential.final_force_norm << ","
                      << tangential.friction_limit << ","
                      << tangent_ratio << ","
                      << tangential.stored_energy_before_gate << ","
                      << tangential.stored_energy_after_gate << ","
                      << energy_ratio << ","
                      << source_energy_bound << ","
                      << inherited_energy << ","
                      << inherited_energy_ratio << ","
                      << (tangential.state == StickSlipState::Stick ? "stick" : "slip") << std::endl;
        }

        double force_jump = frame > 0 ? (total_force - previous_total_force).Length() : 0.0;
        double torque_jump = frame > 0 ? (total_torque - previous_total_torque).Length() : 0.0;
        previous_total_force = total_force;
        previous_total_torque = total_torque;

        max_patch_count = std::max(max_patch_count, static_cast<int>(patches.size()));
        max_source_count_all = std::max(max_source_count_all, max_source_count);
        max_previous_reuse_all = std::max(max_previous_reuse_all, max_previous_reuse);
        max_tangent_ratio_all = std::max(max_tangent_ratio_all, max_tangent_ratio);
        max_energy_gate_ratio_all = std::max(max_energy_gate_ratio_all, max_energy_gate_ratio);
        max_inherited_energy_ratio_all = std::max(max_inherited_energy_ratio_all, max_inherited_energy_ratio);
        max_force_jump = std::max(max_force_jump, force_jump);
        max_torque_jump = std::max(max_torque_jump, torque_jump);
        if (patches.size() == 1) {
            frames_with_one++;
        }
        if (patches.size() == 2) {
            frames_with_two++;
        }
        if (merge_count > 0) {
            frames_with_merge++;
        }
        if (split_count > 0) {
            frames_with_split++;
        }

        frame_csv << frame << ","
                  << std::fixed << std::setprecision(8)
                  << time << ","
                  << y << ","
                  << patches.size() << ","
                  << newborn_count << ","
                  << merge_count << ","
                  << split_count << ","
                  << death_count << ","
                  << max_source_count << ","
                  << max_previous_reuse << ","
                  << total_force.x() << ","
                  << total_force.y() << ","
                  << total_force.z() << ","
                  << total_torque.x() << ","
                  << total_torque.y() << ","
                  << total_torque.z() << ","
                  << force_jump << ","
                  << torque_jump << ","
                  << max_tangent_ratio << ","
                  << max_energy_gate_ratio << ","
                  << max_inherited_energy_ratio << std::endl;

        previous_snapshots = current_snapshots;
        history_store = new_history_store;
    }

    frame_csv.close();
    patch_csv.close();

    std::ofstream summary(out_dir + "/field_contact_split_merge_summary.csv");
    summary << "metric,value" << std::endl;
    summary << std::fixed << std::setprecision(8)
            << "frames," << frame_count << std::endl
            << "plane_vertices," << plane_mesh.vertices.size() << std::endl
            << "plane_faces," << plane_mesh.faces.size() << std::endl
            << "target_vertices," << target_mesh.vertices.size() << std::endl
            << "target_faces," << target_mesh.faces.size() << std::endl
            << "target_sdf_active_voxels," << target_sdf.grid->activeVoxelCount() << std::endl
            << "frames_with_one_patch," << frames_with_one << std::endl
            << "frames_with_two_patches," << frames_with_two << std::endl
            << "frames_with_merge," << frames_with_merge << std::endl
            << "frames_with_split," << frames_with_split << std::endl
            << "total_newborn," << total_newborn << std::endl
            << "total_death," << total_death << std::endl
            << "total_merge_patches," << total_merge_patches << std::endl
            << "total_split_patches," << total_split_patches << std::endl
            << "max_patch_count," << max_patch_count << std::endl
            << "max_source_count," << max_source_count_all << std::endl
            << "max_previous_reuse," << max_previous_reuse_all << std::endl
            << "max_tangential_force_ratio," << max_tangent_ratio_all << std::endl
            << "max_energy_gate_ratio," << max_energy_gate_ratio_all << std::endl
            << "max_inherited_energy_ratio," << max_inherited_energy_ratio_all << std::endl
            << "max_force_jump," << max_force_jump << std::endl
            << "max_torque_jump," << max_torque_jump << std::endl;
    summary.close();

    std::cout << "Output:" << std::endl;
    std::cout << "  " << out_dir + "/field_contact_split_merge_frames.csv" << std::endl;
    std::cout << "  " << out_dir + "/field_contact_split_merge_patches.csv" << std::endl;
    std::cout << "  " << out_dir + "/field_contact_split_merge_summary.csv" << std::endl;
    std::cout << "Frames with merge: " << frames_with_merge << std::endl;
    std::cout << "Frames with split: " << frames_with_split << std::endl;
    std::cout << "Max source count: " << max_source_count_all << std::endl;
    std::cout << "Max previous reuse: " << max_previous_reuse_all << std::endl;
    std::cout << "Max tangential force / Coulomb limit: " << max_tangent_ratio_all << std::endl;

    return 0;
}
