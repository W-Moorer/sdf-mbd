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
//   - simple_gear
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
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
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

struct TimingRow {
    std::string case_name;
    double elapsed_seconds = 0.0;
};

struct Mat3 {
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    ChVector3d operator*(const ChVector3d& v) const {
        return ChVector3d(m[0][0] * v.x() + m[0][1] * v.y() + m[0][2] * v.z(),
                          m[1][0] * v.x() + m[1][1] * v.y() + m[1][2] * v.z(),
                          m[2][0] * v.x() + m[2][1] * v.y() + m[2][2] * v.z());
    }
};

struct RmdBodyInfo {
    ChVector3d cm_marker_mm = ChVector3d(0, 0, 0);
    Mat3 part_rotation;
    ChVector3d surface_ref_marker_mm = ChVector3d(0, 0, 0);
    Mat3 surface_ref_rotation;
    double inertia_x_kg_m2 = 0.0;
    bool has_part = false;
    bool has_cm = false;
    bool has_surface_ref = false;
    bool has_inertia = false;
};

struct GearPose {
    ChVector3d center = ChVector3d(0, 0, 0);
    Mat3 initial_rotation;
};

struct GearReferenceRow {
    double time = 0.0;
    double omega_rx = 0.0;
    double alpha_rx = 0.0;
};

struct RecurDynContactSettings {
    double bpen = 1.0e-5;
    double max_pen = 6.0e-5;
    int korder = 2;
    double recurdyn_k = 100000.0;
    double recurdyn_c = 10.0;
    double pressure_at_bpen = 2.0e5;
    double damping_pressure = 1.0e8;
};

struct SimpleGearRmd {
    RmdBodyInfo gear21;
    RmdBodyInfo gear22;
    RecurDynContactSettings contact;
};

struct DirectionalContactResult {
    ChVector3d torque_on_surface = ChVector3d(0, 0, 0);
    ChVector3d torque_on_target = ChVector3d(0, 0, 0);
    ChVector3d force_on_surface = ChVector3d(0, 0, 0);
    ChVector3d force_on_target = ChVector3d(0, 0, 0);
    ChVector3d application_point_sum = ChVector3d(0, 0, 0);
    double application_point_weight = 0.0;
    ChVector3d max_penetration_point = ChVector3d(0, 0, 0);
    ChVector3d max_penetration_normal = ChVector3d(0, 0, 1);
    double active_area = 0.0;
    double max_pressure = 0.0;
    double min_phi = std::numeric_limits<double>::max();
    double max_effective_penetration = 0.0;
    int active_samples = 0;
    int patch_count = 0;
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

static Mat3 Multiply(const Mat3& a, const Mat3& b) {
    Mat3 out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out.m[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                out.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return out;
}

static Mat3 Transpose(const Mat3& a) {
    Mat3 out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out.m[i][j] = a.m[j][i];
        }
    }
    return out;
}

static Mat3 RotX(double angle) {
    Mat3 out;
    double c = std::cos(angle);
    double s = std::sin(angle);
    out.m[1][1] = c;
    out.m[1][2] = -s;
    out.m[2][1] = s;
    out.m[2][2] = c;
    return out;
}

static Mat3 RotY(double angle) {
    Mat3 out;
    double c = std::cos(angle);
    double s = std::sin(angle);
    out.m[0][0] = c;
    out.m[0][2] = s;
    out.m[2][0] = -s;
    out.m[2][2] = c;
    return out;
}

static Mat3 RotZ(double angle) {
    Mat3 out;
    double c = std::cos(angle);
    double s = std::sin(angle);
    out.m[0][0] = c;
    out.m[0][1] = -s;
    out.m[1][0] = s;
    out.m[1][1] = c;
    return out;
}

static Mat3 RecurDynEuler(double rx, double ry, double rz) {
    return Multiply(RotZ(rz), Multiply(RotY(ry), RotX(rx)));
}

static ChVector3d AngularVelocityXCross(double omega, const ChVector3d& r) {
    return ChVector3d(0, -omega * r.z(), omega * r.y());
}

static Mat3 BodyRotation(const GearPose& pose, double theta_rx) {
    return Multiply(RotX(theta_rx), pose.initial_rotation);
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

static std::string Trim(std::string text) {
    auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char c) {
        return std::isspace(c) != 0;
    });
    auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char c) {
        return std::isspace(c) != 0;
    }).base();
    if (first >= last) {
        return "";
    }
    return std::string(first, last);
}

static bool Contains(const std::string& text, const std::string& pattern) {
    return text.find(pattern) != std::string::npos;
}

static std::string ExtractQuotedName(const std::string& line) {
    size_t first = line.find('\'');
    if (first == std::string::npos) {
        return "";
    }
    size_t second = line.find('\'', first + 1);
    if (second == std::string::npos) {
        return "";
    }
    return line.substr(first + 1, second - first - 1);
}

static std::string CleanRmdNumberToken(std::string token) {
    token = Trim(token);
    while (!token.empty() && (token.back() == ',' || token.back() == ';')) {
        token.pop_back();
    }
    for (size_t i = 0; i < token.size();) {
        if (token[i] == 'D' || token[i] == 'd') {
            if (i + 1 < token.size() &&
                (token[i + 1] == '+' || token[i + 1] == '-' ||
                 std::isdigit(static_cast<unsigned char>(token[i + 1])) != 0)) {
                token[i] = 'E';
                i++;
            } else {
                token.erase(token.begin() + static_cast<std::ptrdiff_t>(i));
            }
        } else {
            i++;
        }
    }
    return token;
}

static std::vector<double> ParseNumbersAfterEquals(const std::string& line) {
    size_t eq = line.find('=');
    if (eq == std::string::npos) {
        return {};
    }
    std::string rhs = line.substr(eq + 1);
    std::replace(rhs.begin(), rhs.end(), ',', ' ');
    std::istringstream stream(rhs);
    std::vector<double> values;
    std::string token;
    while (stream >> token) {
        token = CleanRmdNumberToken(token);
        if (token.empty()) {
            continue;
        }
        try {
            values.push_back(std::stod(token));
        } catch (...) {
        }
    }
    return values;
}

static ChVector3d ParseRmdTriple(const std::string& line) {
    auto values = ParseNumbersAfterEquals(line);
    if (values.size() < 3) {
        throw std::runtime_error("Expected RMD triple in line: " + line);
    }
    return ChVector3d(values[0], values[1], values[2]);
}

static void RequireBodyInfo(const RmdBodyInfo& body, const std::string& name) {
    if (!body.has_part || !body.has_cm || !body.has_surface_ref || !body.has_inertia) {
        throw std::runtime_error("Incomplete RMD body data for " + name);
    }
}

static SimpleGearRmd LoadSimpleGearRmd(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open RMD: " + path.string());
    }

    enum class Block { None, Part, Marker, Contact };
    Block block = Block::None;
    SimpleGearRmd rmd;
    std::string current_part_name;
    std::string current_marker_name;

    std::string line;
    while (std::getline(in, line)) {
        std::string trimmed = Trim(line);
        if (trimmed.rfind("PART /", 0) == 0) {
            block = Block::Part;
            current_part_name.clear();
            current_marker_name.clear();
            continue;
        }
        if (trimmed.rfind("MARKER /", 0) == 0) {
            block = Block::Marker;
            current_part_name.clear();
            current_marker_name.clear();
            continue;
        }
        if (trimmed.rfind("GGEOMCONTACT /", 0) == 0) {
            block = Block::Contact;
            current_part_name.clear();
            current_marker_name.clear();
            continue;
        }

        if (Contains(trimmed, "NAME =")) {
            std::string name = ExtractQuotedName(trimmed);
            if (block == Block::Part) {
                current_part_name = name;
            } else if (block == Block::Marker) {
                current_marker_name = name;
            }
        }

        auto body_for_part = [&]() -> RmdBodyInfo* {
            if (current_part_name == "GEAR21") {
                return &rmd.gear21;
            }
            if (current_part_name == "GEAR22") {
                return &rmd.gear22;
            }
            return nullptr;
        };

        if (block == Block::Part) {
            if (Contains(trimmed, "IP =")) {
                if (auto* body = body_for_part()) {
                    auto values = ParseNumbersAfterEquals(trimmed);
                    if (!values.empty()) {
                        body->inertia_x_kg_m2 = values[0] * 1.0e-6;
                        body->has_inertia = true;
                    }
                }
            } else if (Contains(trimmed, "REULER =")) {
                if (auto* body = body_for_part()) {
                    ChVector3d euler = ParseRmdTriple(trimmed);
                    body->part_rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
                    body->has_part = true;
                }
            }
        } else if (block == Block::Marker) {
            RmdBodyInfo* cm_body = nullptr;
            RmdBodyInfo* ref_body = nullptr;
            if (current_marker_name == "GEAR21.CM") {
                cm_body = &rmd.gear21;
            } else if (current_marker_name == "GEAR22.CM") {
                cm_body = &rmd.gear22;
            } else if (Contains(current_marker_name, "GEAR21.BaseGSurfacePatchRefMarker")) {
                ref_body = &rmd.gear21;
            } else if (Contains(current_marker_name, "GEAR22.BaseGSurfacePatchRefMarker")) {
                ref_body = &rmd.gear22;
            }

            if (Contains(trimmed, "QP =")) {
                if (cm_body) {
                    cm_body->cm_marker_mm = ParseRmdTriple(trimmed);
                    cm_body->has_cm = true;
                } else if (ref_body) {
                    ref_body->surface_ref_marker_mm = ParseRmdTriple(trimmed);
                    ref_body->has_surface_ref = true;
                }
            } else if (Contains(trimmed, "REULER =") && ref_body) {
                ChVector3d euler = ParseRmdTriple(trimmed);
                ref_body->surface_ref_rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
            }
        } else if (block == Block::Contact) {
            auto values = ParseNumbersAfterEquals(trimmed);
            if (values.empty()) {
                continue;
            }
            if (Contains(trimmed, "BPEN =")) {
                rmd.contact.bpen = values[0] * 1.0e-3;
            } else if (Contains(trimmed, "MAXPEN =")) {
                rmd.contact.max_pen = values[0] * 1.0e-3;
            } else if (Contains(trimmed, "KORDER =")) {
                rmd.contact.korder = std::max(1, static_cast<int>(std::lround(values[0])));
            } else if (trimmed.rfind(", K =", 0) == 0 || trimmed.rfind("K =", 0) == 0) {
                rmd.contact.recurdyn_k = values[0];
            } else if (trimmed.rfind(", C =", 0) == 0 || trimmed.rfind("C =", 0) == 0) {
                rmd.contact.recurdyn_c = values[0];
            }
        }
    }

    RequireBodyInfo(rmd.gear21, "GEAR21");
    RequireBodyInfo(rmd.gear22, "GEAR22");
    return rmd;
}

static Mesh LoadObjAsBodyLocal(const std::filesystem::path& path,
                               const ChVector3d& surface_ref_marker_mm,
                               const Mat3& surface_ref_rotation,
                               const ChVector3d& body_cm_mm,
                               const Mat3& body_initial_rotation) {
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
            ChVector3d world_mm = surface_ref_marker_mm + surface_ref_rotation * ChVector3d(x, y, z);
            ChVector3d body_local_mm = Transpose(body_initial_rotation) * (world_mm - body_cm_mm);
            mesh.vertices.emplace_back(body_local_mm * 1.0e-3);
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

static std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::string token;
    std::istringstream stream(line);
    while (std::getline(stream, token, ',')) {
        out.push_back(token);
    }
    return out;
}

static std::vector<GearReferenceRow> LoadGear22Reference(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open reference CSV: " + path.string());
    }

    std::string line;
    std::getline(in, line);
    std::vector<GearReferenceRow> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        auto cols = SplitCsvLine(line);
        if (cols.size() < 10) {
            continue;
        }
        GearReferenceRow row;
        row.time = std::stod(cols[0]);
        row.omega_rx = std::stod(cols[6]);
        row.alpha_rx = std::stod(cols[9]);
        rows.push_back(row);
    }
    return rows;
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

static SparseSDF BuildAnalyticSphereSDF(double radius, double voxel_size, float half_width_voxels) {
    double half_width = static_cast<double>(half_width_voxels) * voxel_size;
    auto grid = openvdb::FloatGrid::create(static_cast<float>(half_width));
    auto transform = openvdb::math::Transform::createLinearTransform(voxel_size);
    grid->setTransform(transform);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    auto accessor = grid->getAccessor();
    int extent = static_cast<int>(std::ceil((radius + half_width) / voxel_size)) + 3;
    for (int i = -extent; i <= extent; i++) {
        for (int j = -extent; j <= extent; j++) {
            for (int k = -extent; k <= extent; k++) {
                openvdb::Vec3d p = transform->indexToWorld(openvdb::Vec3d(i, j, k));
                double phi = std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()) - radius;
                if (std::abs(phi) <= half_width) {
                    accessor.setValue(openvdb::Coord(i, j, k), static_cast<float>(phi));
                }
            }
        }
    }
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

static double RecurDynStylePressure(double effective_penetration,
                                    double closing_normal_velocity,
                                    const RecurDynContactSettings& settings) {
    if (effective_penetration <= 0.0) {
        return 0.0;
    }
    double x = std::min(effective_penetration, settings.max_pen) / settings.bpen;
    double elastic = settings.pressure_at_bpen * std::pow(x, static_cast<double>(settings.korder));
    double damping = settings.damping_pressure * x * std::max(-closing_normal_velocity, 0.0);
    return elastic + damping;
}

static DirectionalContactResult EvaluateDirectionalContact(const SurfaceGraph& surface_graph,
                                                           const SparseSDF& target_sdf,
                                                           const GearPose& surface_pose,
                                                           const GearPose& target_pose,
                                                           double theta_surface,
                                                           double theta_target,
                                                           double omega_surface,
                                                           double omega_target,
                                                           const RecurDynContactSettings& settings) {
    DirectionalContactResult result;
    std::vector<int> active_indices;
    active_indices.reserve(surface_graph.samples.size());
    SdfSampler sampler(target_sdf.grid->tree(), target_sdf.grid->transform());

    const Mat3 r_surface = BodyRotation(surface_pose, theta_surface);
    const Mat3 r_target = BodyRotation(target_pose, theta_target);
    const Mat3 r_target_t = Transpose(r_target);

    for (const auto& sample : surface_graph.samples) {
        ChVector3d r_surface_world = r_surface * sample.local_pos;
        ChVector3d world_pos = surface_pose.center + r_surface_world;
        ChVector3d r_target_world = world_pos - target_pose.center;
        ChVector3d target_local = r_target_t * r_target_world;
        double phi = SampleSDF(target_sdf, sampler, target_local);
        result.min_phi = std::min(result.min_phi, phi);

        double effective_penetration = std::min(std::max(settings.bpen - phi, 0.0), settings.max_pen);
        if (effective_penetration <= 0.0 || sample.area <= 0.0) {
            continue;
        }

        ChVector3d normal_world = SafeNormalize(r_target * SDFGradient(target_sdf, sampler, target_local),
                                                SafeNormalize(r_target_world, ChVector3d(0, 0, 1)));
        ChVector3d v_surface = AngularVelocityXCross(omega_surface, r_surface_world);
        ChVector3d v_target = AngularVelocityXCross(omega_target, r_target_world);
        double normal_velocity = (v_surface - v_target).Dot(normal_world);
        double pressure = RecurDynStylePressure(effective_penetration, normal_velocity, settings);
        ChVector3d dforce = normal_world * (pressure * sample.area);

        result.force_on_surface += dforce;
        result.force_on_target += dforce * -1.0;
        result.torque_on_surface += r_surface_world.Cross(dforce);
        result.torque_on_target += r_target_world.Cross(dforce * -1.0);
        result.active_area += sample.area;
        double point_weight = std::abs(pressure) * sample.area;
        result.application_point_sum += world_pos * point_weight;
        result.application_point_weight += point_weight;
        if (effective_penetration > result.max_effective_penetration) {
            result.max_effective_penetration = effective_penetration;
            result.max_penetration_point = world_pos;
            result.max_penetration_normal = normal_world;
            result.max_pressure = pressure;
        }
        result.active_samples++;
        active_indices.push_back(sample.id);
    }

    result.patch_count = static_cast<int>(surface_graph.FindConnectedComponents(active_indices).size());
    return result;
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
               << "0,0,0,0,0,0," << step.stats.patch_count << "," << min_phi << ",\n";

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
               << ref_ax << "," << ref_bx << "," << last_patch_count << "," << last_min_phi << ",\n";
    }

    std::vector<SummaryRow> out;
    out.push_back(Summarize(case_name, "sphere_a_x", ax_errors, 2.0e-3));
    out.push_back(Summarize(case_name, "sphere_b_x", bx_errors, 2.0e-3));
    out.push_back(Summarize(case_name, "sphere_a_vx", av_errors, 5.0e-2));
    out.push_back(Summarize(case_name, "sphere_b_vx", bv_errors, 5.0e-2));
    return out;
}

static double Rms(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum / static_cast<double>(values.size()));
}

static double MaxAbsValue(const std::vector<double>& values) {
    double out = 0.0;
    for (double value : values) {
        out = std::max(out, std::abs(value));
    }
    return out;
}

enum class HeadonContactMode {
    Bidirectional,
    AOnBOnly,
    BOnAOnly
};

static const char* HeadonModeName(HeadonContactMode mode) {
    switch (mode) {
        case HeadonContactMode::Bidirectional:
            return "bidirectional";
        case HeadonContactMode::AOnBOnly:
            return "a_on_b_only";
        case HeadonContactMode::BOnAOnly:
            return "b_on_a_only";
    }
    return "unknown";
}

struct HeadonValidationResult {
    std::string suite;
    std::string variant;
    std::string mode;
    double voxel_size = 0.0;
    double contact_substep = 0.0;
    int rings_a = 0;
    int sectors_a = 0;
    int rings_b = 0;
    int sectors_b = 0;
    double final_time = 0.0;
    double final_av = 0.0;
    double final_bv = 0.0;
    double ref_av = 0.0;
    double ref_bv = 0.0;
    double restitution = 0.0;
    double restitution_error = 0.0;
    double momentum_error = 0.0;
    double energy_rel_error = 0.0;
    double velocity_l2_error = 0.0;
    double max_penetration = 0.0;
    double contact_duration = 0.0;
    double impulse_error = 0.0;
    double elapsed_seconds = 0.0;
};

static HeadonValidationResult SimulateHeadonValidation(const std::string& suite,
                                                       const std::string& variant,
                                                       double voxel_size,
                                                       double contact_substep,
                                                       HeadonContactMode mode,
                                                       int rings_a,
                                                       int sectors_a,
                                                       int rings_b,
                                                       int sectors_b) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    const double radius = 0.05;
    const double rho_a = 1000.0;
    const double rho_b = 1000.0;
    const double ax0 = -0.15;
    const double bx0 = 0.15;
    const double av0 = 1.0;
    const double bv0 = 0.0;
    const double total_time = 0.23;

    Mesh sphere_mesh = BuildClosedSphereMesh(radius, 32, 64);
    SparseSDF sphere_sdf = BuildSDF(sphere_mesh, voxel_size, 32.0f);
    SurfaceGraph graph_a = MakeSphereSurfaceGraph(radius, rings_a, sectors_a);
    SurfaceGraph graph_b = MakeSphereSurfaceGraph(radius, rings_b, sectors_b);
    FieldContactPrimitiveTracker tracker_a_on_b;
    FieldContactPrimitiveTracker tracker_b_on_a;
    FieldContactRuntimeSettings settings = MakeSettings(contact_substep, 8.0 * sphere_sdf.voxel_size, 2.0e12, 0.0);

    const double offset_a = SphereMinPhi(graph_a, sphere_sdf, 0.0, 2.0 * radius);
    const double offset_b = SphereMinPhi(graph_b, sphere_sdf, 0.0, 2.0 * radius);

    const double volume = 4.0 * kPi * radius * radius * radius / 3.0;
    const double ma = rho_a * volume;
    const double mb = rho_b * volume;
    const double ref_av = ElasticPostV1(av0, bv0, ma, mb);
    const double ref_bv = ElasticPostV2(av0, bv0, ma, mb);
    const double initial_momentum = ma * av0 + mb * bv0;
    const double initial_energy = 0.5 * ma * av0 * av0 + 0.5 * mb * bv0 * bv0;
    const double analytic_impulse_a = ma * (ref_av - av0);

    double ax = ax0;
    double bx = bx0;
    double av = av0;
    double bv = bv0;
    double time = 0.0;
    double max_penetration = 0.0;
    double contact_duration = 0.0;
    double impulse_a = 0.0;

    auto build_queries = [&](const SurfaceGraph& graph,
                             double sdf_offset,
                             double surface_x,
                             double target_x,
                             double surface_v,
                             double target_v) {
        std::vector<FieldSampleQuery> queries;
        queries.reserve(graph.samples.size());
        SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());
        for (const auto& sample : graph.samples) {
            ChVector3d local_to_target = sample.local_pos + ChVector3d(surface_x - target_x, 0, 0);
            FieldSampleQuery query;
            query.world_pos = sample.local_pos + ChVector3d(surface_x, 0, 0);
            query.world_vel = ChVector3d(surface_v - target_v, 0, 0);
            query.phi = SampleSDF(sphere_sdf, sampler, local_to_target) - sdf_offset;
            query.grad = query.phi < settings.extraction.activation_band ?
                             SDFGradient(sphere_sdf, sampler, local_to_target) :
                             ChVector3d(surface_x < target_x ? -1.0 : 1.0, 0, 0);
            queries.push_back(query);
        }
        return queries;
    };

    auto evaluate_forces = [&]() {
        auto queries_a = build_queries(graph_a, offset_a, ax, bx, av, bv);
        auto queries_b = build_queries(graph_b, offset_b, bx, ax, bv, av);
        FieldContactStepResult a_on_b =
            tracker_a_on_b.Evaluate(graph_a, queries_a, ChVector3d(ax, 0, 0), settings);
        FieldContactStepResult b_on_a =
            tracker_b_on_a.Evaluate(graph_b, queries_b, ChVector3d(bx, 0, 0), settings);

        ChVector3d force_a(0, 0, 0);
        ChVector3d force_b(0, 0, 0);
        if (mode == HeadonContactMode::Bidirectional) {
            FieldContactPairResult pair =
                CombineBidirectionalFieldContactPair(a_on_b, ChVector3d(ax, 0, 0), b_on_a, ChVector3d(bx, 0, 0));
            force_a = pair.on_a.force;
            force_b = pair.on_b.force;
        } else if (mode == HeadonContactMode::AOnBOnly) {
            force_a = a_on_b.total_force;
            force_b = -a_on_b.total_force;
        } else {
            force_a = -b_on_a.total_force;
            force_b = b_on_a.total_force;
        }

        double min_phi = std::min(SphereMinPhi(graph_a, sphere_sdf, ax, bx) - offset_a,
                                  SphereMinPhi(graph_b, sphere_sdf, bx, ax) - offset_b);
        max_penetration = std::max(max_penetration, std::max(0.0, -min_phi));
        return std::pair<ChVector3d, ChVector3d>(force_a, force_b);
    };

    while (time < total_time - 1.0e-14) {
        double remaining = total_time - time;
        double gap = bx - ax - 2.0 * radius;
        double rel_closing = std::max(0.0, av - bv);
        bool near_contact = gap < 0.0025;
        double h = near_contact ? std::min(remaining, contact_substep) : remaining;
        if (!near_contact && rel_closing > 1.0e-12) {
            h = std::min(h, std::max(2.0e-5, 0.2 * gap / rel_closing));
        }

        ChVector3d force_a(0, 0, 0);
        ChVector3d force_b(0, 0, 0);
        if (gap < 0.0012) {
            auto forces = evaluate_forces();
            force_a = forces.first;
            force_b = forces.second;
            if (force_a.Length() > 1.0e-9 || force_b.Length() > 1.0e-9) {
                contact_duration += h;
            }
        }

        impulse_a += force_a.x() * h;
        av += (force_a.x() / ma) * h;
        bv += (force_b.x() / mb) * h;
        ax += av * h;
        bx += bv * h;
        time += h;
    }

    HeadonValidationResult result;
    result.suite = suite;
    result.variant = variant;
    result.mode = HeadonModeName(mode);
    result.voxel_size = voxel_size;
    result.contact_substep = contact_substep;
    result.rings_a = rings_a;
    result.sectors_a = sectors_a;
    result.rings_b = rings_b;
    result.sectors_b = sectors_b;
    result.final_time = total_time;
    result.final_av = av;
    result.final_bv = bv;
    result.ref_av = ref_av;
    result.ref_bv = ref_bv;
    result.restitution = (bv - av) / std::max(av0 - bv0, 1.0e-12);
    result.restitution_error = std::abs(result.restitution - 1.0);
    result.momentum_error = std::abs((ma * av + mb * bv) - initial_momentum);
    double final_energy = 0.5 * ma * av * av + 0.5 * mb * bv * bv;
    result.energy_rel_error = std::abs(final_energy - initial_energy) / std::max(initial_energy, 1.0e-12);
    result.velocity_l2_error = std::sqrt((av - ref_av) * (av - ref_av) + (bv - ref_bv) * (bv - ref_bv));
    result.max_penetration = max_penetration;
    result.contact_duration = contact_duration;
    result.impulse_error = std::abs(impulse_a - analytic_impulse_a);
    result.elapsed_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

static ChVector3d DirectionalApplicationPoint(const DirectionalContactResult& result, const ChVector3d& fallback) {
    return result.application_point_weight > 1.0e-16 ?
               result.application_point_sum / result.application_point_weight :
               fallback;
}

static ChVector3d DirectionalMaxPenetrationPoint(const DirectionalContactResult& result, const ChVector3d& fallback) {
    return result.active_samples > 0 ? result.max_penetration_point : fallback;
}

static double PearsonCorrelation(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.size() < 2) {
        return 0.0;
    }
    double mean_a = 0.0;
    double mean_b = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        mean_a += a[i];
        mean_b += b[i];
    }
    mean_a /= static_cast<double>(a.size());
    mean_b /= static_cast<double>(b.size());

    double cov = 0.0;
    double var_a = 0.0;
    double var_b = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double da = a[i] - mean_a;
        double db = b[i] - mean_b;
        cov += da * db;
        var_a += da * da;
        var_b += db * db;
    }
    double denom = std::sqrt(var_a * var_b);
    return denom > 1.0e-24 ? cov / denom : 0.0;
}

enum class GearAblationMode {
    BidirectionalDistributed,
    OneWayGear22Surface,
    OneWayGear21Surface,
    SymmetricSinglePoint,
    SymmetricMaxPenetrationPoint
};

static const char* GearModeName(GearAblationMode mode) {
    switch (mode) {
        case GearAblationMode::BidirectionalDistributed:
            return "bidirectional_distributed";
        case GearAblationMode::OneWayGear22Surface:
            return "one_way_gear22_surface";
        case GearAblationMode::OneWayGear21Surface:
            return "one_way_gear21_surface";
        case GearAblationMode::SymmetricSinglePoint:
            return "symmetric_single_point";
        case GearAblationMode::SymmetricMaxPenetrationPoint:
            return "sdf_max_penetration_point";
    }
    return "unknown";
}

struct GearAblationMetrics {
    std::string mode;
    double duration = 0.0;
    double rms_error = 0.0;
    double max_abs_error = 0.0;
    double max_omega_jump = 0.0;
    double max_alpha_jump = 0.0;
    double max_torque_jump = 0.0;
    double max_patch_event_torque_jump = 0.0;
    double patch_count_abs_error_correlation = 0.0;
    double mean_patch_count = 0.0;
    double mean_active_samples = 0.0;
    double max_effective_penetration = 0.0;
    double max_penetration_to_bpen = 0.0;
    double max_penetration_to_voxel = 0.0;
    double rms_torque = 0.0;
    double max_abs_torque = 0.0;
    double max_force_norm = 0.0;
    double net_contact_work = 0.0;
    double positive_contact_work = 0.0;
    double elapsed_seconds = 0.0;
};

struct GearSensitivityVariant {
    std::string group;
    std::string label;
    double voxel_size = 2.5e-5;
    double time_step = 1.0e-5;
    double pressure_at_bpen = 2.0e5;
    double damping_pressure = 1.0e8;
    double bpen = 1.0e-5;
};

struct GearSensitivityMetrics {
    GearSensitivityVariant variant;
    double duration = 0.0;
    double rms_error = 0.0;
    double max_abs_error = 0.0;
    double final_abs_error = 0.0;
    double max_omega_jump = 0.0;
    double max_alpha_jump = 0.0;
    double max_torque_jump = 0.0;
    double max_patch_event_omega_jump = 0.0;
    double patch_count_abs_error_correlation = 0.0;
    double patch_jump_omega_jump_correlation = 0.0;
    double mean_patch_count = 0.0;
    int max_patch_count = 0;
    double mean_active_samples = 0.0;
    double max_effective_penetration = 0.0;
    double elapsed_seconds = 0.0;
};

static void RunSimpleGearAblation(const std::filesystem::path& case_dir,
                                  const std::filesystem::path& out_dir) {
    SimpleGearRmd rmd = LoadSimpleGearRmd(case_dir / "simple gear.rmd");
    rmd.contact.pressure_at_bpen = 2.0e5;
    rmd.contact.damping_pressure = 1.0e8;

    GearPose pose21{rmd.gear21.cm_marker_mm * 1.0e-3, rmd.gear21.part_rotation};
    GearPose pose22{rmd.gear22.cm_marker_mm * 1.0e-3, rmd.gear22.part_rotation};
    Mesh gear21 = LoadObjAsBodyLocal(case_dir / "models" / "gear_21.obj",
                                     rmd.gear21.surface_ref_marker_mm,
                                     rmd.gear21.surface_ref_rotation,
                                     rmd.gear21.cm_marker_mm,
                                     rmd.gear21.part_rotation);
    Mesh gear22 = LoadObjAsBodyLocal(case_dir / "models" / "gear_22.obj",
                                     rmd.gear22.surface_ref_marker_mm,
                                     rmd.gear22.surface_ref_rotation,
                                     rmd.gear22.cm_marker_mm,
                                     rmd.gear22.part_rotation);
    SparseSDF gear21_sdf = BuildSDF(gear21, 2.5e-5, 20.0f);
    SparseSDF gear22_sdf = BuildSDF(gear22, 2.5e-5, 20.0f);
    SurfaceGraph gear21_graph = MakeTriangleMeshSurfaceGraph(gear21.vertices, gear21.faces);
    SurfaceGraph gear22_graph = MakeTriangleMeshSurfaceGraph(gear22.vertices, gear22.faces);
    auto reference = LoadGear22Reference(case_dir / "data" / "Gear22.csv");

    std::ofstream frames(out_dir / "gear_ablation_frames.csv");
    frames << std::setprecision(16);
    frames << "mode,time,omega22,analytic_reference,recurdyn_reference,analytic_error,recurdyn_error,"
              "torque,alpha22,force_norm,contact_power,patch_count,active_samples,min_phi,"
              "max_effective_penetration\n";

    std::vector<GearAblationMetrics> metrics;
    const double duration = 0.10;
    const double omega21 = 1.0;
    const double analytic_omega22 = -omega21;
    const double dt = 1.0e-5;
    const double startup_time = 0.02;
    std::vector<GearAblationMode> modes = {GearAblationMode::BidirectionalDistributed,
                                           GearAblationMode::OneWayGear22Surface,
                                           GearAblationMode::OneWayGear21Surface,
                                           GearAblationMode::SymmetricSinglePoint,
                                           GearAblationMode::SymmetricMaxPenetrationPoint};

    for (GearAblationMode mode : modes) {
        using Clock = std::chrono::steady_clock;
        const auto mode_start = Clock::now();
        double theta21 = 0.0;
        double theta22 = 0.0;
        double omega22 = 0.0;
        double time_integrated = 0.0;
        double last_torque = 0.0;
        double last_alpha = 0.0;
        double last_force_norm = 0.0;
        double last_min_phi = 0.0;
        double last_max_pen = 0.0;
        int last_patch_count = 0;
        int last_active_samples = 0;
        std::vector<double> errors;
        std::vector<double> abs_errors;
        std::vector<double> patch_counts_for_errors;
        std::vector<double> omega_jumps;
        std::vector<double> alpha_jumps;
        std::vector<double> torque_jumps;
        std::vector<double> patch_event_torque_jumps;
        std::vector<double> torques;
        double patch_sum = 0.0;
        double active_sample_sum = 0.0;
        double max_pen_observed = 0.0;
        double max_force_norm = 0.0;
        double net_contact_work = 0.0;
        double positive_contact_work = 0.0;
        int patch_samples = 0;
        bool have_previous_output = false;
        double previous_omega = 0.0;
        double previous_alpha = 0.0;
        double previous_torque = 0.0;
        int previous_patch_count = 0;

        for (const auto& ref : reference) {
            if (ref.time > duration + 1.0e-14) {
                break;
            }
            while (time_integrated < ref.time - 1.0e-14) {
                double step = std::min(dt, ref.time - time_integrated);
                DirectionalContactResult gear22_on_gear21 =
                    EvaluateDirectionalContact(gear22_graph,
                                               gear21_sdf,
                                               pose22,
                                               pose21,
                                               theta22,
                                               theta21,
                                               omega22,
                                               omega21,
                                               rmd.contact);
                DirectionalContactResult gear21_on_gear22 =
                    EvaluateDirectionalContact(gear21_graph,
                                               gear22_sdf,
                                               pose21,
                                               pose22,
                                               theta21,
                                               theta22,
                                               omega21,
                                               omega22,
                                               rmd.contact);

                ChVector3d selected_force(0, 0, 0);
                if (mode == GearAblationMode::BidirectionalDistributed) {
                    last_torque = 0.5 * (gear22_on_gear21.torque_on_surface.x() +
                                         gear21_on_gear22.torque_on_target.x());
                    selected_force = 0.5 * (gear22_on_gear21.force_on_surface +
                                            gear21_on_gear22.force_on_target);
                } else if (mode == GearAblationMode::OneWayGear22Surface) {
                    last_torque = gear22_on_gear21.torque_on_surface.x();
                    selected_force = gear22_on_gear21.force_on_surface;
                } else if (mode == GearAblationMode::OneWayGear21Surface) {
                    last_torque = gear21_on_gear22.torque_on_target.x();
                    selected_force = gear21_on_gear22.force_on_target;
                } else if (mode == GearAblationMode::SymmetricSinglePoint) {
                    bool has_22 = gear22_on_gear21.force_on_surface.Length() > 1.0e-14;
                    bool has_21 = gear21_on_gear22.force_on_target.Length() > 1.0e-14;
                    ChVector3d force(0, 0, 0);
                    ChVector3d point = pose22.center;
                    if (has_22 && has_21) {
                        force = 0.5 * (gear22_on_gear21.force_on_surface + gear21_on_gear22.force_on_target);
                        point = 0.5 * (DirectionalApplicationPoint(gear22_on_gear21, pose22.center) +
                                       DirectionalApplicationPoint(gear21_on_gear22, pose22.center));
                    } else if (has_22) {
                        force = gear22_on_gear21.force_on_surface;
                        point = DirectionalApplicationPoint(gear22_on_gear21, pose22.center);
                    } else if (has_21) {
                        force = gear21_on_gear22.force_on_target;
                        point = DirectionalApplicationPoint(gear21_on_gear22, pose22.center);
                    }
                    selected_force = force;
                    last_torque = (point - pose22.center).Cross(force).x();
                } else {
                    bool has_22 = gear22_on_gear21.force_on_surface.Length() > 1.0e-14;
                    bool has_21 = gear21_on_gear22.force_on_target.Length() > 1.0e-14;
                    ChVector3d force(0, 0, 0);
                    ChVector3d point = pose22.center;
                    if (has_22 && has_21) {
                        force = 0.5 * (gear22_on_gear21.force_on_surface + gear21_on_gear22.force_on_target);
                        point = 0.5 * (DirectionalMaxPenetrationPoint(gear22_on_gear21, pose22.center) +
                                       DirectionalMaxPenetrationPoint(gear21_on_gear22, pose22.center));
                    } else if (has_22) {
                        force = gear22_on_gear21.force_on_surface;
                        point = DirectionalMaxPenetrationPoint(gear22_on_gear21, pose22.center);
                    } else if (has_21) {
                        force = gear21_on_gear22.force_on_target;
                        point = DirectionalMaxPenetrationPoint(gear21_on_gear22, pose22.center);
                    }
                    selected_force = force;
                    last_torque = (point - pose22.center).Cross(force).x();
                }

                last_patch_count = gear22_on_gear21.patch_count + gear21_on_gear22.patch_count;
                last_active_samples = gear22_on_gear21.active_samples + gear21_on_gear22.active_samples;
                last_min_phi = std::min(gear22_on_gear21.min_phi, gear21_on_gear22.min_phi);
                last_max_pen = std::max(gear22_on_gear21.max_effective_penetration,
                                        gear21_on_gear22.max_effective_penetration);
                last_force_norm = selected_force.Length();
                max_pen_observed = std::max(max_pen_observed, last_max_pen);
                max_force_norm = std::max(max_force_norm, last_force_norm);
                last_alpha = last_torque / rmd.gear22.inertia_x_kg_m2;
                const double contact_power = last_torque * omega22;
                const double contact_work_increment = contact_power * step;
                net_contact_work += contact_work_increment;
                if (contact_work_increment > 0.0) {
                    positive_contact_work += contact_work_increment;
                }
                omega22 += last_alpha * step;
                theta22 += omega22 * step;
                theta21 += omega21 * step;
                time_integrated += step;
            }

            double error = omega22 - analytic_omega22;
            double recurdyn_error = omega22 - ref.omega_rx;
            if (ref.time >= startup_time) {
                errors.push_back(error);
                abs_errors.push_back(std::abs(error));
                patch_counts_for_errors.push_back(static_cast<double>(last_patch_count));
                torques.push_back(last_torque);
                if (have_previous_output) {
                    omega_jumps.push_back(std::abs(omega22 - previous_omega));
                    alpha_jumps.push_back(std::abs(last_alpha - previous_alpha));
                    torque_jumps.push_back(std::abs(last_torque - previous_torque));
                    if (last_patch_count != previous_patch_count) {
                        patch_event_torque_jumps.push_back(std::abs(last_torque - previous_torque));
                    }
                }
            }
            have_previous_output = true;
            previous_omega = omega22;
            previous_alpha = last_alpha;
            previous_torque = last_torque;
            previous_patch_count = last_patch_count;
            patch_sum += static_cast<double>(last_patch_count);
            active_sample_sum += static_cast<double>(last_active_samples);
            patch_samples++;

            frames << GearModeName(mode) << "," << ref.time << "," << omega22 << "," << analytic_omega22 << ","
                   << ref.omega_rx << "," << error << "," << recurdyn_error << "," << last_torque << ","
                   << last_alpha << "," << last_force_norm << "," << last_torque * omega22 << ","
                   << last_patch_count << "," << last_active_samples << ","
                   << last_min_phi << "," << last_max_pen << "\n";
        }

        GearAblationMetrics row;
        row.mode = GearModeName(mode);
        row.duration = duration;
        row.rms_error = Rms(errors);
        row.max_abs_error = MaxAbsValue(errors);
        row.max_omega_jump = MaxAbsValue(omega_jumps);
        row.max_alpha_jump = MaxAbsValue(alpha_jumps);
        row.max_torque_jump = MaxAbsValue(torque_jumps);
        row.max_patch_event_torque_jump = MaxAbsValue(patch_event_torque_jumps);
        row.patch_count_abs_error_correlation = PearsonCorrelation(patch_counts_for_errors, abs_errors);
        row.mean_patch_count = patch_samples > 0 ? patch_sum / static_cast<double>(patch_samples) : 0.0;
        row.mean_active_samples = patch_samples > 0 ? active_sample_sum / static_cast<double>(patch_samples) : 0.0;
        row.max_effective_penetration = max_pen_observed;
        row.max_penetration_to_bpen = rmd.contact.bpen > 0.0 ? max_pen_observed / rmd.contact.bpen : 0.0;
        row.max_penetration_to_voxel = gear21_sdf.voxel_size > 0.0 ? max_pen_observed / gear21_sdf.voxel_size : 0.0;
        row.rms_torque = Rms(torques);
        row.max_abs_torque = MaxAbsValue(torques);
        row.max_force_norm = max_force_norm;
        row.net_contact_work = net_contact_work;
        row.positive_contact_work = positive_contact_work;
        row.elapsed_seconds = std::chrono::duration<double>(Clock::now() - mode_start).count();
        metrics.push_back(row);
    }

    std::ofstream summary(out_dir / "gear_ablation_summary.csv");
    summary << std::setprecision(16);
    summary << "mode,duration,rms_error,max_abs_error,max_omega_jump,max_alpha_jump,max_torque_jump,"
               "max_patch_event_torque_jump,patch_count_abs_error_correlation,mean_patch_count,"
               "mean_active_samples,max_effective_penetration,max_penetration_to_bpen,max_penetration_to_voxel,"
               "rms_torque,max_abs_torque,max_force_norm,net_contact_work,positive_contact_work,"
               "elapsed_seconds\n";
    for (const auto& row : metrics) {
        summary << row.mode << "," << row.duration << "," << row.rms_error << "," << row.max_abs_error << ","
                << row.max_omega_jump << "," << row.max_alpha_jump << "," << row.max_torque_jump << ","
                << row.max_patch_event_torque_jump << "," << row.patch_count_abs_error_correlation << ","
                << row.mean_patch_count << "," << row.mean_active_samples << "," << row.max_effective_penetration
                << "," << row.max_penetration_to_bpen << "," << row.max_penetration_to_voxel << ","
                << row.rms_torque << "," << row.max_abs_torque << "," << row.max_force_norm << ","
                << row.net_contact_work << "," << row.positive_contact_work << "," << row.elapsed_seconds << "\n";
    }
}

static GearSensitivityMetrics SimulateSimpleGearSensitivity(const SimpleGearRmd& base_rmd,
                                                            const GearPose& pose21,
                                                            const GearPose& pose22,
                                                            const Mesh& gear21,
                                                            const Mesh& gear22,
                                                            const SurfaceGraph& gear21_graph,
                                                            const SurfaceGraph& gear22_graph,
                                                            const std::vector<GearReferenceRow>& reference,
                                                            const GearSensitivityVariant& variant) {
    using Clock = std::chrono::steady_clock;
    const auto variant_start = Clock::now();

    SimpleGearRmd rmd = base_rmd;
    rmd.contact.pressure_at_bpen = variant.pressure_at_bpen;
    rmd.contact.damping_pressure = variant.damping_pressure;
    rmd.contact.bpen = variant.bpen;

    SparseSDF gear21_sdf = BuildSDF(gear21, variant.voxel_size, 20.0f);
    SparseSDF gear22_sdf = BuildSDF(gear22, variant.voxel_size, 20.0f);

    const double duration = 0.10;
    const double omega21 = 1.0;
    const double analytic_omega22 = -omega21;
    const double startup_time = 0.02;

    double theta21 = 0.0;
    double theta22 = 0.0;
    double omega22 = 0.0;
    double time_integrated = 0.0;
    double last_torque = 0.0;
    double last_alpha = 0.0;
    double last_max_pen = 0.0;
    int last_patch_count = 0;
    int last_active_samples = 0;

    std::vector<double> errors;
    std::vector<double> abs_errors;
    std::vector<double> patch_counts_for_errors;
    std::vector<double> omega_jumps;
    std::vector<double> alpha_jumps;
    std::vector<double> torque_jumps;
    std::vector<double> patch_jumps;
    std::vector<double> patch_event_omega_jumps;
    double patch_sum = 0.0;
    double active_sample_sum = 0.0;
    double max_pen_observed = 0.0;
    int patch_samples = 0;
    int max_patch_count = 0;
    bool have_previous_output = false;
    double previous_omega = 0.0;
    double previous_alpha = 0.0;
    double previous_torque = 0.0;
    int previous_patch_count = 0;

    for (const auto& ref : reference) {
        if (ref.time > duration + 1.0e-14) {
            break;
        }
        while (time_integrated < ref.time - 1.0e-14) {
            double step = std::min(variant.time_step, ref.time - time_integrated);
            DirectionalContactResult gear22_on_gear21 =
                EvaluateDirectionalContact(gear22_graph,
                                           gear21_sdf,
                                           pose22,
                                           pose21,
                                           theta22,
                                           theta21,
                                           omega22,
                                           omega21,
                                           rmd.contact);
            DirectionalContactResult gear21_on_gear22 =
                EvaluateDirectionalContact(gear21_graph,
                                           gear22_sdf,
                                           pose21,
                                           pose22,
                                           theta21,
                                           theta22,
                                           omega21,
                                           omega22,
                                           rmd.contact);

            last_torque = 0.5 * (gear22_on_gear21.torque_on_surface.x() +
                                 gear21_on_gear22.torque_on_target.x());
            last_patch_count = gear22_on_gear21.patch_count + gear21_on_gear22.patch_count;
            last_active_samples = gear22_on_gear21.active_samples + gear21_on_gear22.active_samples;
            last_max_pen = std::max(gear22_on_gear21.max_effective_penetration,
                                    gear21_on_gear22.max_effective_penetration);
            max_pen_observed = std::max(max_pen_observed, last_max_pen);
            last_alpha = last_torque / rmd.gear22.inertia_x_kg_m2;

            omega22 += last_alpha * step;
            theta22 += omega22 * step;
            theta21 += omega21 * step;
            time_integrated += step;
        }

        double error = omega22 - analytic_omega22;
        if (ref.time >= startup_time) {
            errors.push_back(error);
            abs_errors.push_back(std::abs(error));
            patch_counts_for_errors.push_back(static_cast<double>(last_patch_count));
            if (have_previous_output) {
                double omega_jump = std::abs(omega22 - previous_omega);
                double patch_jump = std::abs(static_cast<double>(last_patch_count - previous_patch_count));
                omega_jumps.push_back(omega_jump);
                alpha_jumps.push_back(std::abs(last_alpha - previous_alpha));
                torque_jumps.push_back(std::abs(last_torque - previous_torque));
                patch_jumps.push_back(patch_jump);
                if (last_patch_count != previous_patch_count) {
                    patch_event_omega_jumps.push_back(omega_jump);
                }
            }
        }

        have_previous_output = true;
        previous_omega = omega22;
        previous_alpha = last_alpha;
        previous_torque = last_torque;
        previous_patch_count = last_patch_count;
        patch_sum += static_cast<double>(last_patch_count);
        active_sample_sum += static_cast<double>(last_active_samples);
        max_patch_count = std::max(max_patch_count, last_patch_count);
        patch_samples++;
    }

    GearSensitivityMetrics row;
    row.variant = variant;
    row.duration = duration;
    row.rms_error = Rms(errors);
    row.max_abs_error = MaxAbsValue(errors);
    row.final_abs_error = std::abs(omega22 - analytic_omega22);
    row.max_omega_jump = MaxAbsValue(omega_jumps);
    row.max_alpha_jump = MaxAbsValue(alpha_jumps);
    row.max_torque_jump = MaxAbsValue(torque_jumps);
    row.max_patch_event_omega_jump = MaxAbsValue(patch_event_omega_jumps);
    row.patch_count_abs_error_correlation = PearsonCorrelation(patch_counts_for_errors, abs_errors);
    row.patch_jump_omega_jump_correlation = PearsonCorrelation(patch_jumps, omega_jumps);
    row.mean_patch_count = patch_samples > 0 ? patch_sum / static_cast<double>(patch_samples) : 0.0;
    row.max_patch_count = max_patch_count;
    row.mean_active_samples = patch_samples > 0 ? active_sample_sum / static_cast<double>(patch_samples) : 0.0;
    row.max_effective_penetration = max_pen_observed;
    row.elapsed_seconds = std::chrono::duration<double>(Clock::now() - variant_start).count();
    return row;
}

static void RunSimpleGearSensitivity(const std::filesystem::path& case_dir,
                                     const std::filesystem::path& out_dir) {
    SimpleGearRmd rmd = LoadSimpleGearRmd(case_dir / "simple gear.rmd");
    rmd.contact.pressure_at_bpen = 2.0e5;
    rmd.contact.damping_pressure = 1.0e8;

    GearPose pose21{rmd.gear21.cm_marker_mm * 1.0e-3, rmd.gear21.part_rotation};
    GearPose pose22{rmd.gear22.cm_marker_mm * 1.0e-3, rmd.gear22.part_rotation};
    Mesh gear21 = LoadObjAsBodyLocal(case_dir / "models" / "gear_21.obj",
                                     rmd.gear21.surface_ref_marker_mm,
                                     rmd.gear21.surface_ref_rotation,
                                     rmd.gear21.cm_marker_mm,
                                     rmd.gear21.part_rotation);
    Mesh gear22 = LoadObjAsBodyLocal(case_dir / "models" / "gear_22.obj",
                                     rmd.gear22.surface_ref_marker_mm,
                                     rmd.gear22.surface_ref_rotation,
                                     rmd.gear22.cm_marker_mm,
                                     rmd.gear22.part_rotation);
    SurfaceGraph gear21_graph = MakeTriangleMeshSurfaceGraph(gear21.vertices, gear21.faces);
    SurfaceGraph gear22_graph = MakeTriangleMeshSurfaceGraph(gear22.vertices, gear22.faces);
    auto reference = LoadGear22Reference(case_dir / "data" / "Gear22.csv");
    if (reference.empty()) {
        throw std::runtime_error("Simple gear reference is empty");
    }

    const GearSensitivityVariant base{"baseline", "baseline", 2.5e-5, 1.0e-5, 2.0e5, 1.0e8, rmd.contact.bpen};
    std::vector<GearSensitivityVariant> variants = {
        base,
        GearSensitivityVariant{"voxel", "50um", 5.0e-5, base.time_step, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"voxel", "35um", 3.5e-5, base.time_step, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"voxel", "25um", 2.5e-5, base.time_step, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"time_step", "20us", base.voxel_size, 2.0e-5, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"time_step", "10us", base.voxel_size, 1.0e-5, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"time_step", "5us", base.voxel_size, 5.0e-6, base.pressure_at_bpen, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"pressure", "1e5Pa", base.voxel_size, base.time_step, 1.0e5, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"pressure", "2e5Pa", base.voxel_size, base.time_step, 2.0e5, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"pressure", "4e5Pa", base.voxel_size, base.time_step, 4.0e5, base.damping_pressure, base.bpen},
        GearSensitivityVariant{"damping", "5e7", base.voxel_size, base.time_step, base.pressure_at_bpen, 5.0e7, base.bpen},
        GearSensitivityVariant{"damping", "1e8", base.voxel_size, base.time_step, base.pressure_at_bpen, 1.0e8, base.bpen},
        GearSensitivityVariant{"damping", "2e8", base.voxel_size, base.time_step, base.pressure_at_bpen, 2.0e8, base.bpen},
        GearSensitivityVariant{"bpen", "5um", base.voxel_size, base.time_step, base.pressure_at_bpen, base.damping_pressure, 5.0e-6},
        GearSensitivityVariant{"bpen", "10um", base.voxel_size, base.time_step, base.pressure_at_bpen, base.damping_pressure, 1.0e-5},
        GearSensitivityVariant{"bpen", "20um", base.voxel_size, base.time_step, base.pressure_at_bpen, base.damping_pressure, 2.0e-5},
    };

    std::ofstream summary(out_dir / "gear_sensitivity_summary.csv");
    summary << std::setprecision(16);
    summary << "group,label,duration,voxel_size,time_step,pressure_at_bpen,damping_pressure,bpen,"
               "rms_error,max_abs_error,final_abs_error,max_omega_jump,max_alpha_jump,max_torque_jump,"
               "max_patch_event_omega_jump,patch_count_abs_error_correlation,patch_jump_omega_jump_correlation,"
               "mean_patch_count,max_patch_count,mean_active_samples,max_effective_penetration,elapsed_seconds\n";
    for (const auto& variant : variants) {
        GearSensitivityMetrics row =
            SimulateSimpleGearSensitivity(rmd, pose21, pose22, gear21, gear22, gear21_graph, gear22_graph, reference, variant);
        summary << row.variant.group << "," << row.variant.label << "," << row.duration << ","
                << row.variant.voxel_size << "," << row.variant.time_step << ","
                << row.variant.pressure_at_bpen << "," << row.variant.damping_pressure << ","
                << row.variant.bpen << "," << row.rms_error << "," << row.max_abs_error << ","
                << row.final_abs_error << "," << row.max_omega_jump << "," << row.max_alpha_jump << ","
                << row.max_torque_jump << "," << row.max_patch_event_omega_jump << ","
                << row.patch_count_abs_error_correlation << "," << row.patch_jump_omega_jump_correlation << ","
                << row.mean_patch_count << "," << row.max_patch_count << "," << row.mean_active_samples << ","
                << row.max_effective_penetration << "," << row.elapsed_seconds << "\n";
    }
}

struct SdfSpatialConvergenceRow {
    double voxel_size = 0.0;
    int sample_count = 0;
    double rms_distance_error = 0.0;
    double max_distance_error = 0.0;
    double rms_normal_angle = 0.0;
    double max_normal_angle = 0.0;
    double active_area_rel_error = 0.0;
    double force_rel_error = 0.0;
    double torque_rel_error = 0.0;
    double elapsed_seconds = 0.0;
};

static SdfSpatialConvergenceRow EvaluateSphereSdfSpatialConvergence(double voxel_size) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    const double radius = 0.005;
    SparseSDF sphere_sdf = BuildAnalyticSphereSDF(radius, voxel_size, 32.0f);
    SurfaceGraph graph = MakeSphereSurfaceGraph(radius, 32, 64);
    SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());

    std::vector<double> distance_errors;
    std::vector<double> normal_angles;
    for (const auto& sample : graph.samples) {
        ChVector3d n_exact = SafeNormalize(sample.local_pos, ChVector3d(1, 0, 0));
        for (double offset : {-1.5 * voxel_size, 0.0, 1.5 * voxel_size}) {
            ChVector3d p = sample.local_pos + n_exact * offset;
            double phi_exact = p.Length() - radius;
            double phi_sdf = SampleSDF(sphere_sdf, sampler, p);
            ChVector3d n_sdf = SDFGradient(sphere_sdf, sampler, p);
            distance_errors.push_back(phi_sdf - phi_exact);
            normal_angles.push_back(std::acos(ClampSigned(n_sdf.Dot(n_exact))));
        }
    }

    const double prescribed_penetration = 8.0e-5;
    const ChVector3d target_center(2.0 * radius - prescribed_penetration, 3.0e-5, 2.0e-5);
    const double pressure_stiffness = 2.0e12;
    ChVector3d force_sdf(0, 0, 0);
    ChVector3d force_exact(0, 0, 0);
    ChVector3d torque_sdf(0, 0, 0);
    ChVector3d torque_exact(0, 0, 0);
    double active_area_sdf = 0.0;
    double active_area_exact = 0.0;

    for (const auto& sample : graph.samples) {
        ChVector3d local_to_target = sample.local_pos - target_center;
        double phi_exact = local_to_target.Length() - radius;
        double phi_sdf = SampleSDF(sphere_sdf, sampler, local_to_target);
        if (phi_exact < 0.0) {
            ChVector3d n_exact = SafeNormalize(local_to_target, ChVector3d(-1, 0, 0));
            ChVector3d dforce = n_exact * (pressure_stiffness * (-phi_exact) * sample.area);
            force_exact += dforce;
            torque_exact += sample.local_pos.Cross(dforce);
            active_area_exact += sample.area;
        }
        if (phi_sdf < 0.0) {
            ChVector3d n_sdf = SDFGradient(sphere_sdf, sampler, local_to_target);
            ChVector3d dforce = n_sdf * (pressure_stiffness * (-phi_sdf) * sample.area);
            force_sdf += dforce;
            torque_sdf += sample.local_pos.Cross(dforce);
            active_area_sdf += sample.area;
        }
    }

    SdfSpatialConvergenceRow row;
    row.voxel_size = voxel_size;
    row.sample_count = static_cast<int>(distance_errors.size());
    row.rms_distance_error = Rms(distance_errors);
    row.max_distance_error = MaxAbsValue(distance_errors);
    row.rms_normal_angle = Rms(normal_angles);
    row.max_normal_angle = MaxAbsValue(normal_angles);
    row.active_area_rel_error =
        active_area_exact > 1.0e-16 ? std::abs(active_area_sdf - active_area_exact) / active_area_exact : 0.0;
    row.force_rel_error = force_exact.Length() > 1.0e-16 ? (force_sdf - force_exact).Length() / force_exact.Length() : 0.0;
    row.torque_rel_error =
        torque_exact.Length() > 1.0e-16 ? (torque_sdf - torque_exact).Length() / torque_exact.Length() : 0.0;
    row.elapsed_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return row;
}

static void WriteSdfSpatialConvergence(const std::filesystem::path& path,
                                       const std::vector<SdfSpatialConvergenceRow>& rows) {
    std::ofstream out(path);
    out << "voxel_size,sample_count,rms_distance_error,max_distance_error,rms_normal_angle,max_normal_angle,"
           "active_area_rel_error,force_rel_error,torque_rel_error,elapsed_seconds\n";
    for (const auto& row : rows) {
        out << row.voxel_size << "," << row.sample_count << "," << row.rms_distance_error << ","
            << row.max_distance_error << "," << row.rms_normal_angle << "," << row.max_normal_angle << ","
            << row.active_area_rel_error << "," << row.force_rel_error << "," << row.torque_rel_error << ","
            << row.elapsed_seconds << "\n";
    }
}

struct SpherePatchWrenchIntegral {
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
    ChVector3d center_sum = ChVector3d(0, 0, 0);
    double center_weight = 0.0;
    double active_area = 0.0;
    int active_samples = 0;

    ChVector3d Center() const {
        return center_weight > 1.0e-16 ? center_sum / center_weight : ChVector3d(0, 0, 0);
    }
};

static SpherePatchWrenchIntegral IntegrateSpherePatchExact(const SurfaceGraph& graph,
                                                           double radius,
                                                           const ChVector3d& target_center,
                                                           const ChVector3d& torque_reference,
                                                           double pressure_stiffness) {
    SpherePatchWrenchIntegral result;
    for (const auto& sample : graph.samples) {
        ChVector3d local_to_target = sample.local_pos - target_center;
        double phi = local_to_target.Length() - radius;
        if (phi >= 0.0 || sample.area <= 0.0) {
            continue;
        }
        ChVector3d normal = SafeNormalize(local_to_target, ChVector3d(-1, 0, 0));
        double scalar_weight = pressure_stiffness * (-phi) * sample.area;
        ChVector3d dforce = normal * scalar_weight;
        result.force += dforce;
        result.torque += (sample.local_pos - torque_reference).Cross(dforce);
        result.center_sum += sample.local_pos * scalar_weight;
        result.center_weight += scalar_weight;
        result.active_area += sample.area;
        result.active_samples++;
    }
    return result;
}

static SpherePatchWrenchIntegral IntegrateSpherePatchSdf(const SurfaceGraph& graph,
                                                         const SparseSDF& sphere_sdf,
                                                         SdfSampler& sampler,
                                                         const ChVector3d& target_center,
                                                         const ChVector3d& torque_reference,
                                                         double pressure_stiffness) {
    SpherePatchWrenchIntegral result;
    for (const auto& sample : graph.samples) {
        ChVector3d local_to_target = sample.local_pos - target_center;
        double phi = SampleSDF(sphere_sdf, sampler, local_to_target);
        if (phi >= 0.0 || sample.area <= 0.0) {
            continue;
        }
        ChVector3d normal = SDFGradient(sphere_sdf, sampler, local_to_target);
        double scalar_weight = pressure_stiffness * (-phi) * sample.area;
        ChVector3d dforce = normal * scalar_weight;
        result.force += dforce;
        result.torque += (sample.local_pos - torque_reference).Cross(dforce);
        result.center_sum += sample.local_pos * scalar_weight;
        result.center_weight += scalar_weight;
        result.active_area += sample.area;
        result.active_samples++;
    }
    return result;
}

static double RelativeVectorError(const ChVector3d& value, const ChVector3d& reference) {
    double ref_norm = reference.Length();
    return ref_norm > 1.0e-16 ? (value - reference).Length() / ref_norm : value.Length();
}

struct PatchWrenchConvergenceRow {
    std::string suite;
    double voxel_size = 0.0;
    int rings = 0;
    int sectors = 0;
    int sample_count = 0;
    int active_samples_exact = 0;
    int active_samples_sdf = 0;
    double active_area_rel_error = 0.0;
    double exact_quad_force_rel_error = 0.0;
    double exact_quad_torque_rel_error = 0.0;
    double sdf_force_rel_error = 0.0;
    double sdf_torque_rel_error = 0.0;
    double total_force_rel_error = 0.0;
    double total_torque_rel_error = 0.0;
    double center_error = 0.0;
    double elapsed_seconds = 0.0;
};

static PatchWrenchConvergenceRow EvaluatePatchWrenchConvergence(const std::string& suite,
                                                                double voxel_size,
                                                                int rings,
                                                                int sectors) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    const double radius = 0.005;
    const double prescribed_penetration = 4.0e-4;
    const ChVector3d target_center(2.0 * radius - prescribed_penetration, 4.0e-4, 3.0e-4);
    const ChVector3d torque_reference(-2.0e-3, 1.5e-3, -1.0e-3);
    const double pressure_stiffness = 2.0e12;

    SurfaceGraph reference_graph = MakeSphereSurfaceGraph(radius, 96, 192);
    SurfaceGraph graph = MakeSphereSurfaceGraph(radius, rings, sectors);
    SparseSDF sphere_sdf = BuildAnalyticSphereSDF(radius, voxel_size, 32.0f);
    SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());

    SpherePatchWrenchIntegral reference =
        IntegrateSpherePatchExact(reference_graph, radius, target_center, torque_reference, pressure_stiffness);
    SpherePatchWrenchIntegral exact =
        IntegrateSpherePatchExact(graph, radius, target_center, torque_reference, pressure_stiffness);
    SpherePatchWrenchIntegral sdf =
        IntegrateSpherePatchSdf(graph, sphere_sdf, sampler, target_center, torque_reference, pressure_stiffness);

    PatchWrenchConvergenceRow row;
    row.suite = suite;
    row.voxel_size = voxel_size;
    row.rings = rings;
    row.sectors = sectors;
    row.sample_count = static_cast<int>(graph.samples.size());
    row.active_samples_exact = exact.active_samples;
    row.active_samples_sdf = sdf.active_samples;
    row.active_area_rel_error =
        reference.active_area > 1.0e-16 ? std::abs(sdf.active_area - reference.active_area) / reference.active_area : 0.0;
    row.exact_quad_force_rel_error = RelativeVectorError(exact.force, reference.force);
    row.exact_quad_torque_rel_error = RelativeVectorError(exact.torque, reference.torque);
    row.sdf_force_rel_error = RelativeVectorError(sdf.force, exact.force);
    row.sdf_torque_rel_error = RelativeVectorError(sdf.torque, exact.torque);
    row.total_force_rel_error = RelativeVectorError(sdf.force, reference.force);
    row.total_torque_rel_error = RelativeVectorError(sdf.torque, reference.torque);
    row.center_error = (sdf.Center() - reference.Center()).Length();
    row.elapsed_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return row;
}

static void WritePatchWrenchConvergence(const std::filesystem::path& path,
                                        const std::vector<PatchWrenchConvergenceRow>& rows) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "suite,voxel_size,rings,sectors,sample_count,active_samples_exact,active_samples_sdf,"
           "active_area_rel_error,exact_quad_force_rel_error,exact_quad_torque_rel_error,"
           "sdf_force_rel_error,sdf_torque_rel_error,total_force_rel_error,total_torque_rel_error,"
           "center_error,elapsed_seconds\n";
    for (const auto& row : rows) {
        out << row.suite << "," << row.voxel_size << "," << row.rings << "," << row.sectors << ","
            << row.sample_count << "," << row.active_samples_exact << "," << row.active_samples_sdf << ","
            << row.active_area_rel_error << "," << row.exact_quad_force_rel_error << ","
            << row.exact_quad_torque_rel_error << "," << row.sdf_force_rel_error << ","
            << row.sdf_torque_rel_error << "," << row.total_force_rel_error << ","
            << row.total_torque_rel_error << "," << row.center_error << "," << row.elapsed_seconds << "\n";
    }
}

struct JcndProfileRow {
    std::string benchmark;
    int evaluations = 0;
    int surface_samples_per_body = 0;
    double voxel_size = 0.0;
    double mean_active_samples = 0.0;
    double mean_patch_count = 0.0;
    double mean_query_ms = 0.0;
    double mean_patch_history_wrench_ms = 0.0;
    double mean_pair_assembly_ms = 0.0;
    double mean_total_ms = 0.0;
};

static JcndProfileRow RunHeadonProfilingStudy() {
    using Clock = std::chrono::steady_clock;
    const double radius = 0.05;
    const double voxel_size = 5.0e-5;
    Mesh sphere_mesh = BuildClosedSphereMesh(radius, 32, 64);
    SparseSDF sphere_sdf = BuildSDF(sphere_mesh, voxel_size, 32.0f);
    SurfaceGraph graph = MakeSphereSurfaceGraph(radius, 12, 24);
    FieldContactPrimitiveTracker tracker_a_on_b;
    FieldContactPrimitiveTracker tracker_b_on_a;
    FieldContactRuntimeSettings settings = MakeSettings(1.0e-5, 8.0 * sphere_sdf.voxel_size, 2.0e12, 0.0);
    const double offset = SphereMinPhi(graph, sphere_sdf, 0.0, 2.0 * radius);

    auto build_queries = [&](double ax, double bx, double av, double bv) {
        std::vector<FieldSampleQuery> queries;
        queries.reserve(graph.samples.size());
        SdfSampler sampler(sphere_sdf.grid->tree(), sphere_sdf.grid->transform());
        for (const auto& sample : graph.samples) {
            ChVector3d local_to_target = sample.local_pos + ChVector3d(ax - bx, 0, 0);
            FieldSampleQuery query;
            query.world_pos = sample.local_pos + ChVector3d(ax, 0, 0);
            query.world_vel = ChVector3d(av - bv, 0, 0);
            query.phi = SampleSDF(sphere_sdf, sampler, local_to_target) - offset;
            query.grad = query.phi < settings.extraction.activation_band ?
                             SDFGradient(sphere_sdf, sampler, local_to_target) :
                             ChVector3d(ax < bx ? -1.0 : 1.0, 0, 0);
            queries.push_back(query);
        }
        return queries;
    };

    const int evaluations = 240;
    double query_seconds = 0.0;
    double patch_seconds = 0.0;
    double pair_seconds = 0.0;
    double active_sum = 0.0;
    double patch_sum = 0.0;

    for (int i = 0; i < evaluations; i++) {
        double phase = static_cast<double>(i) / static_cast<double>(evaluations);
        double penetration = 1.4e-4 + 4.0e-5 * std::sin(2.0 * kPi * phase);
        double ax = 0.0;
        double bx = 2.0 * radius - penetration;
        double av = 0.2;
        double bv = 0.0;

        const auto q0 = Clock::now();
        auto queries_a = build_queries(ax, bx, av, bv);
        auto queries_b = build_queries(bx, ax, bv, av);
        const auto q1 = Clock::now();

        const auto p0 = Clock::now();
        FieldContactStepResult a_on_b = tracker_a_on_b.Evaluate(graph, queries_a, ChVector3d(ax, 0, 0), settings);
        FieldContactStepResult b_on_a = tracker_b_on_a.Evaluate(graph, queries_b, ChVector3d(bx, 0, 0), settings);
        const auto p1 = Clock::now();

        const auto c0 = Clock::now();
        FieldContactPairResult pair =
            CombineBidirectionalFieldContactPair(a_on_b, ChVector3d(ax, 0, 0), b_on_a, ChVector3d(bx, 0, 0));
        const auto c1 = Clock::now();
        (void)pair;

        int active = 0;
        for (const auto& query : queries_a) {
            if (query.phi < settings.extraction.activation_band) {
                active++;
            }
        }
        for (const auto& query : queries_b) {
            if (query.phi < settings.extraction.activation_band) {
                active++;
            }
        }
        active_sum += static_cast<double>(active);
        patch_sum += static_cast<double>(a_on_b.stats.patch_count + b_on_a.stats.patch_count);
        query_seconds += std::chrono::duration<double>(q1 - q0).count();
        patch_seconds += std::chrono::duration<double>(p1 - p0).count();
        pair_seconds += std::chrono::duration<double>(c1 - c0).count();
    }

    JcndProfileRow row;
    row.benchmark = "headon_spheres_fixed_overlap";
    row.evaluations = evaluations;
    row.surface_samples_per_body = static_cast<int>(graph.samples.size());
    row.voxel_size = voxel_size;
    row.mean_active_samples = active_sum / static_cast<double>(evaluations);
    row.mean_patch_count = patch_sum / static_cast<double>(evaluations);
    row.mean_query_ms = 1000.0 * query_seconds / static_cast<double>(evaluations);
    row.mean_patch_history_wrench_ms = 1000.0 * patch_seconds / static_cast<double>(evaluations);
    row.mean_pair_assembly_ms = 1000.0 * pair_seconds / static_cast<double>(evaluations);
    row.mean_total_ms = row.mean_query_ms + row.mean_patch_history_wrench_ms + row.mean_pair_assembly_ms;
    return row;
}

static void WriteProfileSummary(const std::filesystem::path& path, const std::vector<JcndProfileRow>& rows) {
    std::ofstream out(path);
    out << "benchmark,evaluations,surface_samples_per_body,voxel_size,mean_active_samples,mean_patch_count,"
           "mean_query_ms,mean_patch_history_wrench_ms,mean_pair_assembly_ms,mean_total_ms\n";
    for (const auto& row : rows) {
        out << row.benchmark << "," << row.evaluations << "," << row.surface_samples_per_body << ","
            << row.voxel_size << "," << row.mean_active_samples << "," << row.mean_patch_count << ","
            << row.mean_query_ms << "," << row.mean_patch_history_wrench_ms << "," << row.mean_pair_assembly_ms
            << "," << row.mean_total_ms << "\n";
    }
}

static void WriteHeadonValidation(const std::filesystem::path& path,
                                  const std::vector<HeadonValidationResult>& rows) {
    std::ofstream out(path);
    out << "suite,variant,mode,voxel_size,contact_substep,rings_a,sectors_a,rings_b,sectors_b,final_time,"
           "final_av,final_bv,ref_av,ref_bv,restitution,restitution_error,momentum_error,energy_rel_error,"
           "velocity_l2_error,max_penetration,contact_duration,impulse_error,elapsed_seconds\n";
    for (const auto& row : rows) {
        out << row.suite << "," << row.variant << "," << row.mode << "," << row.voxel_size << ","
            << row.contact_substep << "," << row.rings_a << "," << row.sectors_a << "," << row.rings_b << ","
            << row.sectors_b << "," << row.final_time << "," << row.final_av << "," << row.final_bv << ","
            << row.ref_av << "," << row.ref_bv << "," << row.restitution << "," << row.restitution_error << ","
            << row.momentum_error << "," << row.energy_rel_error << "," << row.velocity_l2_error << ","
            << row.max_penetration << "," << row.contact_duration << "," << row.impulse_error << ","
            << row.elapsed_seconds << "\n";
    }
}

static void RunJcndValidation(const std::string& root) {
    const auto out_dir = std::filesystem::path(root) / "out" / "jcnd_validation";
    std::filesystem::create_directories(out_dir);

    std::vector<HeadonValidationResult> headon_rows;
    headon_rows.push_back(SimulateHeadonValidation("physics",
                                                   "baseline",
                                                   5.0e-5,
                                                   5.0e-6,
                                                   HeadonContactMode::Bidirectional,
                                                   12,
                                                   24,
                                                   12,
                                                   24));
    for (double voxel : {9.0e-5, 6.0e-5, 4.0e-5}) {
        headon_rows.push_back(SimulateHeadonValidation("sdf_resolution",
                                                       "voxel_sweep",
                                                       voxel,
                                                       1.0e-5,
                                                       HeadonContactMode::Bidirectional,
                                                       12,
                                                       24,
                                                       12,
                                                       24));
    }
    for (double h : {2.0e-5, 1.0e-5, 5.0e-6}) {
        headon_rows.push_back(SimulateHeadonValidation("time_step",
                                                       "substep_sweep",
                                                       5.0e-5,
                                                       h,
                                                       HeadonContactMode::Bidirectional,
                                                       12,
                                                       24,
                                                       12,
                                                       24));
    }
    headon_rows.push_back(SimulateHeadonValidation("bidirectional_ablation",
                                                   "asymmetric_surface_quadrature",
                                                   5.0e-5,
                                                   1.0e-5,
                                                   HeadonContactMode::Bidirectional,
                                                   6,
                                                   12,
                                                   12,
                                                   24));
    headon_rows.push_back(SimulateHeadonValidation("bidirectional_ablation",
                                                   "asymmetric_surface_quadrature",
                                                   5.0e-5,
                                                   1.0e-5,
                                                   HeadonContactMode::AOnBOnly,
                                                   6,
                                                   12,
                                                   12,
                                                   24));
    headon_rows.push_back(SimulateHeadonValidation("bidirectional_ablation",
                                                   "asymmetric_surface_quadrature",
                                                   5.0e-5,
                                                   1.0e-5,
                                                   HeadonContactMode::BOnAOnly,
                                                   6,
                                                   12,
                                                   12,
                                                   24));
    WriteHeadonValidation(out_dir / "headon_validation.csv", headon_rows);

    std::vector<SdfSpatialConvergenceRow> spatial_rows;
    for (double voxel : {1.0e-4, 7.0e-5, 5.0e-5, 3.5e-5}) {
        spatial_rows.push_back(EvaluateSphereSdfSpatialConvergence(voxel));
    }
    WriteSdfSpatialConvergence(out_dir / "sdf_spatial_convergence.csv", spatial_rows);

    std::vector<PatchWrenchConvergenceRow> patch_wrench_rows;
    for (const auto& resolution : {std::pair<int, int>{8, 16},
                                   std::pair<int, int>{12, 24},
                                   std::pair<int, int>{16, 32},
                                   std::pair<int, int>{24, 48},
                                   std::pair<int, int>{32, 64}}) {
        patch_wrench_rows.push_back(
            EvaluatePatchWrenchConvergence("surface_sampling", 3.5e-5, resolution.first, resolution.second));
    }
    for (double voxel : {1.0e-4, 7.0e-5, 5.0e-5, 3.5e-5}) {
        patch_wrench_rows.push_back(EvaluatePatchWrenchConvergence("sdf_voxel", voxel, 32, 64));
    }
    WritePatchWrenchConvergence(out_dir / "patch_wrench_convergence.csv", patch_wrench_rows);

    std::vector<JcndProfileRow> profile_rows;
    profile_rows.push_back(RunHeadonProfilingStudy());
    WriteProfileSummary(out_dir / "profile_summary.csv", profile_rows);

    const auto simple_gear_dir = std::filesystem::path(root) / "paper_example" / "cases" / "simple_gear";
    RunSimpleGearAblation(simple_gear_dir, out_dir);
    RunSimpleGearSensitivity(simple_gear_dir, out_dir);
}

static std::vector<SummaryRow> RunSimpleGearCase(std::ofstream& frames,
                                                 const std::string& case_name,
                                                 const std::filesystem::path& case_dir) {
    SimpleGearRmd rmd = LoadSimpleGearRmd(case_dir / "simple gear.rmd");
    rmd.contact.pressure_at_bpen = 2.0e5;
    rmd.contact.damping_pressure = 1.0e8;

    GearPose pose21{rmd.gear21.cm_marker_mm * 1.0e-3, rmd.gear21.part_rotation};
    GearPose pose22{rmd.gear22.cm_marker_mm * 1.0e-3, rmd.gear22.part_rotation};

    Mesh gear21 = LoadObjAsBodyLocal(case_dir / "models" / "gear_21.obj",
                                     rmd.gear21.surface_ref_marker_mm,
                                     rmd.gear21.surface_ref_rotation,
                                     rmd.gear21.cm_marker_mm,
                                     rmd.gear21.part_rotation);
    Mesh gear22 = LoadObjAsBodyLocal(case_dir / "models" / "gear_22.obj",
                                     rmd.gear22.surface_ref_marker_mm,
                                     rmd.gear22.surface_ref_rotation,
                                     rmd.gear22.cm_marker_mm,
                                     rmd.gear22.part_rotation);
    SparseSDF gear21_sdf = BuildSDF(gear21, 2.5e-5, 20.0f);
    SparseSDF gear22_sdf = BuildSDF(gear22, 2.5e-5, 20.0f);
    SurfaceGraph gear21_graph = MakeTriangleMeshSurfaceGraph(gear21.vertices, gear21.faces);
    SurfaceGraph gear22_graph = MakeTriangleMeshSurfaceGraph(gear22.vertices, gear22.faces);
    auto reference = LoadGear22Reference(case_dir / "data" / "Gear22.csv");
    if (reference.empty()) {
        throw std::runtime_error("Simple gear reference is empty");
    }

    const double omega21 = 1.0;
    const double analytic_omega22 = -omega21;
    const double dt = 5.0e-6;
    const double startup_time = 0.02;
    double theta21 = 0.0;
    double theta22 = 0.0;
    double omega22 = 0.0;
    double time_integrated = 0.0;
    double last_alpha = 0.0;
    double last_torque = 0.0;
    int last_patch_count = 0;
    double last_min_phi = 0.0;

    std::vector<double> post_startup_errors;
    std::vector<double> post_startup_jumps;
    double previous_output_omega = omega22;
    bool have_previous_output = false;

    for (const auto& ref : reference) {
        while (time_integrated < ref.time - 1.0e-14) {
            double step = std::min(dt, ref.time - time_integrated);
            DirectionalContactResult gear22_on_gear21 =
                EvaluateDirectionalContact(gear22_graph,
                                           gear21_sdf,
                                           pose22,
                                           pose21,
                                           theta22,
                                           theta21,
                                           omega22,
                                           omega21,
                                           rmd.contact);
            DirectionalContactResult gear21_on_gear22 =
                EvaluateDirectionalContact(gear21_graph,
                                           gear22_sdf,
                                           pose21,
                                           pose22,
                                           theta21,
                                           theta22,
                                           omega21,
                                           omega22,
                                           rmd.contact);

            double torque_from_gear22_surface = gear22_on_gear21.torque_on_surface.x();
            double torque_from_gear21_surface = gear21_on_gear22.torque_on_target.x();
            last_torque = 0.5 * (torque_from_gear22_surface + torque_from_gear21_surface);
            last_alpha = last_torque / rmd.gear22.inertia_x_kg_m2;
            last_patch_count = gear22_on_gear21.patch_count + gear21_on_gear22.patch_count;
            last_min_phi = std::min(gear22_on_gear21.min_phi, gear21_on_gear22.min_phi);

            omega22 += last_alpha * step;
            theta22 += omega22 * step;
            theta21 += omega21 * step;
            time_integrated += step;
        }

        double error = omega22 - analytic_omega22;
        if (ref.time >= startup_time) {
            post_startup_errors.push_back(error);
            if (have_previous_output) {
                post_startup_jumps.push_back(std::abs(omega22 - previous_output_omega));
            }
        }
        have_previous_output = true;
        previous_output_omega = omega22;
        frames << case_name << "," << ref.time << ",sparse_sdf_contact_force_dynamics,"
               << omega22 << "," << analytic_omega22 << "," << error << ","
               << "0,0,0,0,0,0," << last_patch_count << "," << last_min_phi << ","
               << ref.omega_rx << "\n";
        (void)last_alpha;
        (void)last_torque;
    }

    const double final_error = std::abs(omega22 - analytic_omega22);
    const double rms_after_startup = Rms(post_startup_errors);
    const double max_jump_after_startup =
        post_startup_jumps.empty() ? 0.0 : *std::max_element(post_startup_jumps.begin(), post_startup_jumps.end());

    std::vector<SummaryRow> out;
    out.push_back(SummaryRow{case_name, "gear22_rx_analytic_final", final_error, final_error, 9.0e-2,
                             final_error <= 9.0e-2});
    out.push_back(SummaryRow{case_name, "gear22_rx_analytic_rms_t_ge_0p02", rms_after_startup, rms_after_startup,
                             5.0e-2, rms_after_startup <= 5.0e-2});
    out.push_back(SummaryRow{case_name, "gear22_rx_max_jump_t_ge_0p02", max_jump_after_startup,
                             max_jump_after_startup, 1.5e-1, max_jump_after_startup <= 1.5e-1});
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

static void WriteTimingSummary(const std::filesystem::path& path,
                               const std::vector<TimingRow>& rows,
                               double total_elapsed_seconds) {
    std::ofstream out(path);
    out << "case,elapsed_seconds\n";
    for (const auto& row : rows) {
        out << row.case_name << "," << row.elapsed_seconds << "\n";
    }
    out << "total," << total_elapsed_seconds << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        using Clock = std::chrono::steady_clock;
        const auto total_start = Clock::now();
        openvdb::initialize();
        const std::string root = GetProjectRoot();
        for (int i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--jcnd-validation") {
                RunJcndValidation(root);
                std::cout << "JCND validation output: "
                          << (std::filesystem::path(root) / "out" / "jcnd_validation").string() << "\n";
                return 0;
            }
        }
        const auto cases = std::filesystem::path(root) / "paper_example" / "cases";
        const auto out_dir = std::filesystem::path(root) / "out" / "paper_example_dynamic_benchmarks";
        std::filesystem::create_directories(out_dir);

        std::ofstream frames(out_dir / "sparse_sdf_frames.csv");
        frames << std::setprecision(16);
        frames << "case,time,backend,backend_y,reference_y,y_error,sphere_a_x,sphere_b_x,"
                  "sphere_a_vx,sphere_b_vx,reference_sphere_a_x,reference_sphere_b_x,patch_count,min_phi,"
                  "recurdyn_y\n";

        std::vector<SummaryRow> summary;
        std::vector<TimingRow> timings;
        auto append = [&](std::vector<SummaryRow> rows) {
            summary.insert(summary.end(), rows.begin(), rows.end());
        };

        auto case_start = Clock::now();
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
        timings.push_back(
            TimingRow{"eccentric_roller",
                      std::chrono::duration<double>(Clock::now() - case_start).count()});

        case_start = Clock::now();
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
        timings.push_back(
            TimingRow{"headon_spheres",
                      std::chrono::duration<double>(Clock::now() - case_start).count()});

        case_start = Clock::now();
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
        timings.push_back(
            TimingRow{"headon_spheres_mass_ratio",
                      std::chrono::duration<double>(Clock::now() - case_start).count()});

        case_start = Clock::now();
        append(RunSimpleGearCase(frames,
                                 "simple_gear",
                                 cases / "simple_gear"));
        timings.push_back(
            TimingRow{"simple_gear",
                      std::chrono::duration<double>(Clock::now() - case_start).count()});

        WriteSummary(out_dir / "comparison_summary.csv", summary);
        const double total_elapsed_seconds = std::chrono::duration<double>(Clock::now() - total_start).count();
        WriteTimingSummary(out_dir / "timing_summary.csv", timings, total_elapsed_seconds);
        bool passed = std::all_of(summary.begin(), summary.end(), [](const SummaryRow& row) {
            return row.passed;
        });
        std::cout << "Wrote " << (out_dir / "sparse_sdf_frames.csv").string() << std::endl;
        std::cout << "Wrote " << (out_dir / "comparison_summary.csv").string() << std::endl;
        std::cout << "Wrote " << (out_dir / "timing_summary.csv").string() << std::endl;
        std::cout << "PASS=" << (passed ? "true" : "false") << std::endl;
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 2;
    }
}
