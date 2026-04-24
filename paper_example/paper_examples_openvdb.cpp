// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Paper example sparse-SDF contact-force dynamics benchmark runner.
//
// This executable exercises the current field-contact runtime on the four
// paper_example benchmark cases:
//   - eccentric_roller
//   - headon_spheres
//   - headon_spheres_mass_ratio
//   - onset_stress
//
// The reported trajectories are advanced from field-contact forces.  The
// reference curves are used only for error measurement.
//
// Outputs:
//   out/paper_example_dynamic_benchmarks/sparse_sdf_frames.csv
//   out/paper_example_dynamic_benchmarks/comparison_summary.csv
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
#include <sstream>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct Mesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct SparseSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 1.0e-4;
};

struct SummaryRow {
    std::string case_name;
    std::string quantity;
    double max_abs_error = 0.0;
    double rms_error = 0.0;
    double tolerance = 0.0;
    bool passed = false;
};

using SdfSampler = openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler>;

static std::string GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "paper_example")) {
            return path.string();
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path().string();
}

static Mesh LoadObj(const std::filesystem::path& path) {
    Mesh mesh;
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open OBJ: " + path.string());
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string tag;
        iss >> tag;
        if (tag == "v") {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            iss >> x >> y >> z;
            mesh.vertices.emplace_back(x, y, z);
        } else if (tag == "f") {
            std::vector<int> ids;
            std::string token;
            while (iss >> token) {
                size_t slash = token.find('/');
                std::string v = slash == std::string::npos ? token : token.substr(0, slash);
                int idx = std::stoi(v);
                if (idx < 0) {
                    idx = static_cast<int>(mesh.vertices.size()) + idx + 1;
                }
                ids.push_back(idx - 1);
            }
            for (size_t i = 1; i + 1 < ids.size(); i++) {
                mesh.faces.push_back(TriangleFace{ids[0], ids[i], ids[i + 1]});
            }
        }
    }
    return mesh;
}

static Mesh BuildClosedSphereMesh(double radius, int rings, int sectors) {
    Mesh mesh;
    mesh.vertices.emplace_back(0.0, radius, 0.0);
    for (int i = 1; i < rings; i++) {
        double theta = kPi * static_cast<double>(i) / static_cast<double>(rings);
        for (int j = 0; j < sectors; j++) {
            double phi = 2.0 * kPi * static_cast<double>(j) / static_cast<double>(sectors);
            mesh.vertices.emplace_back(radius * std::sin(theta) * std::cos(phi),
                                       radius * std::cos(theta),
                                       radius * std::sin(theta) * std::sin(phi));
        }
    }
    int south = static_cast<int>(mesh.vertices.size());
    mesh.vertices.emplace_back(0.0, -radius, 0.0);

    auto vid = [&](int ring, int sector) {
        return 1 + (ring - 1) * sectors + ((sector + sectors) % sectors);
    };

    for (int j = 0; j < sectors; j++) {
        mesh.faces.push_back(TriangleFace{0, vid(1, j + 1), vid(1, j)});
    }
    for (int i = 1; i < rings - 1; i++) {
        for (int j = 0; j < sectors; j++) {
            int a = vid(i, j);
            int b = vid(i, j + 1);
            int c = vid(i + 1, j);
            int d = vid(i + 1, j + 1);
            mesh.faces.push_back(TriangleFace{a, b, c});
            mesh.faces.push_back(TriangleFace{b, d, c});
        }
    }
    for (int j = 0; j < sectors; j++) {
        mesh.faces.push_back(TriangleFace{vid(rings - 1, j), vid(rings - 1, j + 1), south});
    }
    return mesh;
}

static SparseSDF BuildSDF(const Mesh& mesh, double voxel_size, float half_width_voxels) {
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
    return SparseSDF{grid, voxel_size};
}

static ChVector3d RotateZ(const ChVector3d& v, double angle) {
    double c = std::cos(angle);
    double s = std::sin(angle);
    return ChVector3d(c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z());
}

static double SampleSDF(const SparseSDF& sdf, SdfSampler& sampler, const ChVector3d& p) {
    return static_cast<double>(sampler.wsSample(openvdb::Vec3d(p.x(), p.y(), p.z())));
}

static ChVector3d SDFGradient(const SparseSDF& sdf, SdfSampler& sampler, const ChVector3d& p) {
    const double h = 0.5 * sdf.voxel_size;
    ChVector3d grad((SampleSDF(sdf, sampler, p + ChVector3d(h, 0, 0)) -
                         SampleSDF(sdf, sampler, p - ChVector3d(h, 0, 0))) /
                        (2.0 * h),
                    (SampleSDF(sdf, sampler, p + ChVector3d(0, h, 0)) -
                         SampleSDF(sdf, sampler, p - ChVector3d(0, h, 0))) /
                        (2.0 * h),
                    (SampleSDF(sdf, sampler, p + ChVector3d(0, 0, h)) -
                         SampleSDF(sdf, sampler, p - ChVector3d(0, 0, h))) /
                        (2.0 * h));
    return SafeNormalize(grad, ChVector3d(0, 1, 0));
}

static FieldContactRuntimeSettings MakeSettings(double dt,
                                                double activation_band,
                                                double normal_stiffness,
                                                double normal_damping) {
    FieldContactRuntimeSettings settings;
    settings.extraction.activation_band = activation_band;
    settings.extraction.min_area = 0.0;
    settings.extraction.min_samples = 1;
    settings.extraction.use_penetration_weighted_center = true;
    settings.normal.stiffness = normal_stiffness;
    settings.normal.damping = normal_damping;
    settings.tangential.time_step = dt;
    settings.tangential.friction_coefficient = 0.0;
    settings.tangential.stiffness = 0.0;
    settings.tangential.damping = 0.0;
    return settings;
}

static double CylinderMass(double radius, double thickness, double density) {
    return kPi * radius * radius * thickness * density;
}

static double CamEnvelopeY(double eccentricity, double cam_radius, double roller_radius, double theta) {
    double cx = eccentricity * std::cos(theta);
    double cy = eccentricity * std::sin(theta);
    double reach = cam_radius + roller_radius;
    return cy + std::sqrt(std::max(0.0, reach * reach - cx * cx));
}

static std::vector<FieldSampleQuery> BuildCamQueries(const SurfaceGraph& roller_graph,
                                                     const SparseSDF& cam_sdf,
                                                     double cam_angle,
                                                     double cam_speed,
                                                     double follower_y,
                                                     double follower_vy,
                                                     double gradient_band) {
    std::vector<FieldSampleQuery> queries;
    queries.reserve(roller_graph.samples.size());
    SdfSampler sampler(cam_sdf.grid->tree(), cam_sdf.grid->transform());

    for (const auto& sample : roller_graph.samples) {
        ChVector3d world_pos = sample.local_pos + ChVector3d(0, follower_y, 0);
        ChVector3d world_vel(0, follower_vy, 0);
        ChVector3d cam_velocity = ChVector3d(0, 0, cam_speed).Cross(world_pos);
        ChVector3d local_pos = RotateZ(world_pos, -cam_angle);

        FieldSampleQuery query;
        query.world_pos = world_pos;
        query.world_vel = world_vel - cam_velocity;
        query.phi = SampleSDF(cam_sdf, sampler, local_pos);
        if (query.phi < gradient_band) {
            ChVector3d grad_local = SDFGradient(cam_sdf, sampler, local_pos);
            query.grad = RotateZ(grad_local, cam_angle);
        } else {
            query.grad = ChVector3d(0, 1, 0);
        }
        queries.push_back(query);
    }
    return queries;
}

static double MinCamPhi(const SurfaceGraph& roller_graph,
                        const SparseSDF& cam_sdf,
                        double cam_angle,
                        double follower_y) {
    SdfSampler sampler(cam_sdf.grid->tree(), cam_sdf.grid->transform());
    double min_phi = std::numeric_limits<double>::max();
    for (const auto& sample : roller_graph.samples) {
        ChVector3d world_pos = sample.local_pos + ChVector3d(0, follower_y, 0);
        ChVector3d local_pos = RotateZ(world_pos, -cam_angle);
        min_phi = std::min(min_phi, SampleSDF(cam_sdf, sampler, local_pos));
    }
    return min_phi;
}

static double RootFollowerY(const SurfaceGraph& roller_graph,
                            const SparseSDF& cam_sdf,
                            double cam_angle,
                            double reference_y) {
    double low = reference_y - 0.003;
    double high = reference_y + 0.003;
    double f_low = MinCamPhi(roller_graph, cam_sdf, cam_angle, low);
    double f_high = MinCamPhi(roller_graph, cam_sdf, cam_angle, high);
    for (int expand = 0; expand < 12 && !(f_low <= 0.0 && f_high >= 0.0); expand++) {
        low -= 0.002;
        high += 0.002;
        f_low = MinCamPhi(roller_graph, cam_sdf, cam_angle, low);
        f_high = MinCamPhi(roller_graph, cam_sdf, cam_angle, high);
    }
    if (!(f_low <= 0.0 && f_high >= 0.0)) {
        return reference_y;
    }

    for (int iter = 0; iter < 48; iter++) {
        double mid = 0.5 * (low + high);
        double f_mid = MinCamPhi(roller_graph, cam_sdf, cam_angle, mid);
        if (f_mid <= 0.0) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

static SummaryRow Summarize(const std::string& case_name,
                            const std::string& quantity,
                            const std::vector<double>& errors,
                            double tolerance) {
    double max_abs = 0.0;
    double sq = 0.0;
    for (double e : errors) {
        max_abs = std::max(max_abs, std::abs(e));
        sq += e * e;
    }
    double rms = errors.empty() ? 0.0 : std::sqrt(sq / static_cast<double>(errors.size()));
    return SummaryRow{case_name, quantity, max_abs, rms, tolerance, max_abs <= tolerance};
}

static std::vector<SummaryRow> RunCamCase(const std::string& project_root,
                                          std::ofstream& frames,
                                          const std::string& case_name,
                                          const std::filesystem::path& cam_obj,
                                          const std::filesystem::path& roller_obj,
                                          double cam_radius,
                                          double eccentricity,
                                          double roller_radius,
                                          double roller_thickness,
                                          double density,
                                          double gravity_y,
                                          double motor_speed,
                                          double dt,
                                          double total_time,
                                          double phase,
                                          double follower_initial_y,
                                          bool check_onset) {
    (void)project_root;
    Mesh cam_mesh = LoadObj(cam_obj);
    Mesh roller_mesh = LoadObj(roller_obj);
    SurfaceGraph roller_graph = MakeTriangleMeshSurfaceGraph(roller_mesh.vertices, roller_mesh.faces);
    SparseSDF cam_sdf = BuildSDF(cam_mesh, 1.0e-4, 24.0f);
    FieldContactPrimitiveTracker tracker;
    const int substeps_per_output = std::max(1, static_cast<int>(std::ceil(dt / 1.0e-4)));
    const double h = dt / static_cast<double>(substeps_per_output);
    FieldContactRuntimeSettings settings = MakeSettings(h, 2.0 * cam_sdf.voxel_size, 5.0e11, 8.0e6);
    const double mass = CylinderMass(roller_radius, roller_thickness, density);

    int steps = static_cast<int>(std::round(total_time / dt));
    std::vector<double> reference(steps + 1);
    std::vector<double> errors;
    errors.reserve(steps + 1);

    for (int i = 0; i <= steps; i++) {
        double time = static_cast<double>(i) * dt;
        double theta_reference = phase + motor_speed * time;
        reference[i] = CamEnvelopeY(eccentricity, cam_radius, roller_radius, theta_reference);
    }

    double follower_y = follower_initial_y;
    double follower_vy = 0.0;
    double observed_onset = std::numeric_limits<double>::quiet_NaN();
    double previous_phi = MinCamPhi(roller_graph, cam_sdf, 0.0, follower_y);
    double previous_phi_time = 0.0;

    for (int i = 0; i <= steps; i++) {
        double time = static_cast<double>(i) * dt;
        double theta_sdf = motor_speed * time;
        auto queries = BuildCamQueries(roller_graph,
                                       cam_sdf,
                                       theta_sdf,
                                       motor_speed,
                                       follower_y,
                                       follower_vy,
                                       settings.extraction.activation_band);
        FieldContactStepResult step =
            tracker.Evaluate(roller_graph, queries, ChVector3d(0, follower_y, 0), settings);
        double min_phi = MinCamPhi(roller_graph, cam_sdf, theta_sdf, follower_y);
        double err = follower_y - reference[i];
        errors.push_back(err);

        frames << case_name << "," << time << ",sparse_sdf_contact_force_dynamics,"
               << follower_y << "," << reference[i] << "," << err << ","
               << "0,0,0,0,0,0," << step.stats.patch_count << "," << min_phi << "\n";

        if (i == steps) {
            continue;
        }

        for (int sub = 0; sub < substeps_per_output; sub++) {
            const double t_sub = time + static_cast<double>(sub) * h;
            const double angle_sub = motor_speed * t_sub;
            auto sub_queries = BuildCamQueries(roller_graph,
                                               cam_sdf,
                                               angle_sub,
                                               motor_speed,
                                               follower_y,
                                               follower_vy,
                                               settings.extraction.activation_band);
            FieldContactStepResult sub_step =
                tracker.Evaluate(roller_graph, sub_queries, ChVector3d(0, follower_y, 0), settings);
            if (check_onset) {
                double phi = MinCamPhi(roller_graph, cam_sdf, angle_sub, follower_y);
                if (std::isnan(observed_onset) && previous_phi > 0.0 && phi <= 0.0) {
                    double alpha = previous_phi / (previous_phi - phi);
                    observed_onset = previous_phi_time + alpha * (t_sub - previous_phi_time);
                }
                previous_phi = phi;
                previous_phi_time = t_sub;
            }

            double force_y = sub_step.total_force.y() + mass * gravity_y;
            double ay = force_y / mass;
            follower_vy += ay * h;
            follower_y += follower_vy * h;
        }
    }

    std::vector<SummaryRow> out;
    if (!check_onset) {
        out.push_back(Summarize(case_name, "follower_y", errors, 6.0e-4));
    }
    if (check_onset) {
        double onset_error = observed_onset - 0.15;
        out.push_back(SummaryRow{case_name, "onset_time", std::abs(onset_error), std::abs(onset_error), dt,
                                 std::abs(onset_error) <= dt});
    }
    return out;
}

static double ElasticPostV1(double v1, double v2, double m1, double m2) {
    return ((m1 - m2) / (m1 + m2)) * v1 + (2.0 * m2 / (m1 + m2)) * v2;
}

static double ElasticPostV2(double v1, double v2, double m1, double m2) {
    return (2.0 * m1 / (m1 + m2)) * v1 + ((m2 - m1) / (m1 + m2)) * v2;
}

static double SphereMinPhi(const SurfaceGraph& graph,
                           const SparseSDF& sphere_sdf,
                           double ax,
                           double bx) {
    SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());
    double min_phi = std::numeric_limits<double>::max();
    for (const auto& sample : graph.samples) {
        ChVector3d local_to_b = sample.local_pos + ChVector3d(ax - bx, 0, 0);
        min_phi = std::min(min_phi, SampleSDF(sphere_sdf, sampler, local_to_b));
    }
    return min_phi;
}

static std::vector<SummaryRow> RunHeadonCase(std::ofstream& frames,
                                             const std::string& case_name,
                                             double radius,
                                             double rho_a,
                                             double rho_b,
                                             double ax0,
                                             double bx0,
                                             double av0,
                                             double bv0,
                                             double dt,
                                             double total_time) {
    Mesh sphere_mesh = BuildClosedSphereMesh(radius, 48, 96);
    SparseSDF sphere_sdf = BuildSDF(sphere_mesh, 5.0e-5, 64.0f);
    SurfaceGraph graph = MakeSphereSurfaceGraph(radius, 16, 32);
    FieldContactPrimitiveTracker tracker_a_on_b;
    FieldContactPrimitiveTracker tracker_b_on_a;
    FieldContactRuntimeSettings settings = MakeSettings(5.0e-7, 8.0 * sphere_sdf.voxel_size, 1.0e14, 0.0);
    const double sdf_contact_offset = SphereMinPhi(graph, sphere_sdf, 0.0, 2.0 * radius);

    int steps = static_cast<int>(std::round(total_time / dt));
    double volume = 4.0 * kPi * radius * radius * radius / 3.0;
    double ma = rho_a * volume;
    double mb = rho_b * volume;
    double av_post = ElasticPostV1(av0, bv0, ma, mb);
    double bv_post = ElasticPostV2(av0, bv0, ma, mb);

    double impact_time = (bx0 - ax0 - 2.0 * radius) / (av0 - bv0);
    double ax_hit = ax0 + av0 * impact_time;
    double bx_hit = bx0 + bv0 * impact_time;

    std::vector<double> ax_errors;
    std::vector<double> bx_errors;
    std::vector<double> av_errors;
    std::vector<double> bv_errors;

    double ax = ax0;
    double bx = bx0;
    double av = av0;
    double bv = bv0;
    int last_patch_count = 0;
    double last_min_phi = std::numeric_limits<double>::max();

    auto build_queries = [&](double surface_x, double target_x, double surface_v, double target_v) {
        std::vector<FieldSampleQuery> queries;
        queries.reserve(graph.samples.size());
        SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());
        for (const auto& sample : graph.samples) {
            ChVector3d local_to_target = sample.local_pos + ChVector3d(surface_x - target_x, 0, 0);
            FieldSampleQuery query;
            query.world_pos = sample.local_pos + ChVector3d(surface_x, 0, 0);
            query.world_vel = ChVector3d(surface_v - target_v, 0, 0);
            query.phi = SampleSDF(sphere_sdf, sampler, local_to_target) - sdf_contact_offset;
            query.grad = query.phi < settings.extraction.activation_band ?
                             SDFGradient(sphere_sdf, sampler, local_to_target) :
                             ChVector3d(surface_x < target_x ? -1.0 : 1.0, 0, 0);
            queries.push_back(query);
        }
        return queries;
    };

    auto evaluate_pair = [&]() {
        auto queries_a = build_queries(ax, bx, av, bv);
        auto queries_b = build_queries(bx, ax, bv, av);
        FieldContactStepResult a_on_b =
            tracker_a_on_b.Evaluate(graph, queries_a, ChVector3d(ax, 0, 0), settings);
        FieldContactStepResult b_on_a =
            tracker_b_on_a.Evaluate(graph, queries_b, ChVector3d(bx, 0, 0), settings);
        FieldContactPairResult pair =
            CombineBidirectionalFieldContactPair(a_on_b, ChVector3d(ax, 0, 0), b_on_a, ChVector3d(bx, 0, 0));
        last_patch_count = a_on_b.stats.patch_count + b_on_a.stats.patch_count;
        return pair;
    };

    double time_integrated = 0.0;
    for (int i = 0; i <= steps; i++) {
        double time = static_cast<double>(i) * dt;
        if (i > 0) {
            double target_time = time;
            while (time_integrated < target_time - 1.0e-14) {
                double remaining = target_time - time_integrated;
                double gap = bx - ax - 2.0 * radius;
                double rel_closing = std::max(0.0, av - bv);
                bool near_contact = gap < 0.0025;
                double h = near_contact ? std::min(remaining, 5.0e-7) : remaining;
                if (!near_contact && rel_closing > 1.0e-12) {
                    h = std::min(h, std::max(2.0e-5, 0.2 * gap / rel_closing));
                }

                ChVector3d force_a(0, 0, 0);
                ChVector3d force_b(0, 0, 0);
                if (gap < 0.0012) {
                    FieldContactPairResult pair = evaluate_pair();
                    force_a = pair.on_a.force;
                    force_b = pair.on_b.force;
                } else {
                    last_patch_count = 0;
                    last_min_phi = gap;
                }

                av += (force_a.x() / ma) * h;
                bv += (force_b.x() / mb) * h;
                ax += av * h;
                bx += bv * h;
                time_integrated += h;
            }
        } else {
            last_min_phi = SphereMinPhi(graph, sphere_sdf, ax, bx) - sdf_contact_offset;
        }

        double ref_ax = 0.0;
        double ref_bx = 0.0;
        double ref_av = 0.0;
        double ref_bv = 0.0;
        if (time <= impact_time) {
            ref_ax = ax0 + av0 * time;
            ref_bx = bx0 + bv0 * time;
            ref_av = av0;
            ref_bv = bv0;
        } else {
            ref_ax = ax_hit + av_post * (time - impact_time);
            ref_bx = bx_hit + bv_post * (time - impact_time);
            ref_av = av_post;
            ref_bv = bv_post;
        }
        ax_errors.push_back(ax - ref_ax);
        bx_errors.push_back(bx - ref_bx);
        if (std::abs(time - impact_time) > 0.25 * dt) {
            av_errors.push_back(av - ref_av);
            bv_errors.push_back(bv - ref_bv);
        }
        if (std::abs((bx - ax) - 2.0 * radius) < 0.004) {
            last_min_phi = std::min(SphereMinPhi(graph, sphere_sdf, ax, bx), SphereMinPhi(graph, sphere_sdf, bx, ax)) -
                           sdf_contact_offset;
        } else {
            last_min_phi = bx - ax - 2.0 * radius;
        }

        frames << case_name << "," << time << ",sparse_sdf_contact_force_dynamics,"
               << "0,0,0," << ax << "," << bx << "," << av << "," << bv << ","
               << ref_ax << "," << ref_bx << "," << last_patch_count << "," << last_min_phi << "\n";
    }

    std::vector<SummaryRow> out;
    out.push_back(Summarize(case_name, "sphere_a_x", ax_errors, 2.0e-3));
    out.push_back(Summarize(case_name, "sphere_b_x", bx_errors, 2.0e-3));
    out.push_back(Summarize(case_name, "sphere_a_vx", av_errors, 5.0e-2));
    out.push_back(Summarize(case_name, "sphere_b_vx", bv_errors, 5.0e-2));
    return out;
}

static void WriteSummary(const std::filesystem::path& path, const std::vector<SummaryRow>& rows) {
    std::ofstream out(path);
    out << "case,quantity,max_abs_error,rms_error,tolerance,passed\n";
    for (const auto& row : rows) {
        out << row.case_name << "," << row.quantity << "," << row.max_abs_error << ","
            << row.rms_error << "," << row.tolerance << "," << (row.passed ? "true" : "false") << "\n";
    }
}

}  // namespace

int main() {
    try {
        openvdb::initialize();
        const std::string root = GetProjectRoot();
        const auto cases = std::filesystem::path(root) / "paper_example" / "cases";
        const auto out_dir = std::filesystem::path(root) / "out" / "paper_example_dynamic_benchmarks";
        std::filesystem::create_directories(out_dir);

        std::ofstream frames(out_dir / "sparse_sdf_frames.csv");
        frames << std::setprecision(16);
        frames << "case,time,backend,backend_y,reference_y,y_error,sphere_a_x,sphere_b_x,"
                  "sphere_a_vx,sphere_b_vx,reference_sphere_a_x,reference_sphere_b_x,patch_count,min_phi\n";

        std::vector<SummaryRow> summary;
        auto append = [&](std::vector<SummaryRow> rows) {
            summary.insert(summary.end(), rows.begin(), rows.end());
        };

        append(RunCamCase(root,
                          frames,
                          "eccentric_roller",
                          cases / "eccentric_roller" / "models" / "eccentric_disk_cam.obj",
                          cases / "eccentric_roller" / "models" / "roller_follower.obj",
                          0.03,
                          0.006,
                          0.01,
                          0.018,
                          7800.0,
                          -9.81,
                          -2.0,
                          0.002,
                          1.5707963267948966,
                          0.0,
                          0.039547439866570375,
                          false));

        append(RunHeadonCase(frames,
                             "headon_spheres",
                             0.05,
                             1000.0,
                             1000.0,
                             -0.15,
                             0.15,
                             1.0,
                             0.0,
                             0.0005,
                             0.5));

        append(RunHeadonCase(frames,
                             "headon_spheres_mass_ratio",
                             0.05,
                             1000.0,
                             2000.0,
                             -0.15,
                             0.15,
                             0.0,
                             -1.0,
                             0.0005,
                             0.5));

        append(RunCamCase(root,
                          frames,
                          "onset_stress",
                          cases / "onset_stress" / "models" / "onset_cam.obj",
                          cases / "onset_stress" / "models" / "roller_follower.obj",
                          0.03,
                          0.006,
                          0.01,
                          0.018,
                          7800.0,
                          0.0,
                          -2.0,
                          0.001,
                          0.45,
                          kPi,
                          0.04136029035991934,
                          true));

        WriteSummary(out_dir / "comparison_summary.csv", summary);
        bool passed = std::all_of(summary.begin(), summary.end(), [](const SummaryRow& row) {
            return row.passed;
        });
        std::cout << "Wrote " << (out_dir / "sparse_sdf_frames.csv").string() << std::endl;
        std::cout << "Wrote " << (out_dir / "comparison_summary.csv").string() << std::endl;
        std::cout << "PASS=" << (passed ? "true" : "false") << std::endl;
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 2;
    }
}
