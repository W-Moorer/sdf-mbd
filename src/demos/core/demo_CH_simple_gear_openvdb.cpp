// =============================================================================
// Simple gear OpenVDB field-contact dynamics trial.
//
// Drives GEAR21 at +1 rad/s about the global X axis and integrates the RX
// angular speed of GEAR22 from sparse-SDF field-contact torque.
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
#include <stdexcept>
#include <string>
#include <vector>

#include "chrono/collision/ChFieldContactRuntime.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

using namespace chrono;
using namespace chrono::fieldcontact;

namespace {

constexpr double kMmToM = 1.0e-3;

struct Mat3 {
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    ChVector3d operator*(const ChVector3d& v) const {
        return ChVector3d(m[0][0] * v.x() + m[0][1] * v.y() + m[0][2] * v.z(),
                          m[1][0] * v.x() + m[1][1] * v.y() + m[1][2] * v.z(),
                          m[2][0] * v.x() + m[2][1] * v.y() + m[2][2] * v.z());
    }
};

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

struct Mesh {
    std::vector<ChVector3d> vertices;
    std::vector<TriangleFace> faces;
};

struct SparseSDF {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 2.5e-5;
};

struct ReferenceRow {
    double time = 0.0;
    double omega_rx = 0.0;
    double alpha_rx = 0.0;
};

struct GearPose {
    ChVector3d center = ChVector3d(0, 0, 0);
    Mat3 initial_rotation;
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

struct RecurDynContactSettings {
    double bpen = 1.0e-5;          // 0.01 mm
    double max_pen = 6.0e-5;       // 0.06 mm
    int korder = 2;
    double recurdyn_k = 100000.0;
    double recurdyn_c = 10.0;
    double pressure_at_bpen = 1.0e5;
    double damping_pressure = 2.0e7;
};

struct SimpleGearRmd {
    RmdBodyInfo gear21;
    RmdBodyInfo gear22;
    RecurDynContactSettings contact;
};

struct DirectionalContactResult {
    ChVector3d force_on_surface = ChVector3d(0, 0, 0);
    ChVector3d torque_on_surface = ChVector3d(0, 0, 0);
    ChVector3d force_on_target = ChVector3d(0, 0, 0);
    ChVector3d torque_on_target = ChVector3d(0, 0, 0);
    double min_phi = std::numeric_limits<double>::max();
    double max_effective_penetration = 0.0;
    int active_samples = 0;
    int patch_count = 0;
};

using SdfSampler = openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler>;

static std::string GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "assets")) {
            return path.string();
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path().string();
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

static std::string Trim(std::string text) {
    auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; });
    auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();
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
                        body->inertia_x_kg_m2 = values[0] * 1.0e-6;  // kg mm^2 -> kg m^2
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
            } else if (Contains(trimmed, "REULER =")) {
                if (ref_body) {
                    ChVector3d euler = ParseRmdTriple(trimmed);
                    ref_body->surface_ref_rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
                }
            }
        } else if (block == Block::Contact) {
            auto values = ParseNumbersAfterEquals(trimmed);
            if (values.empty()) {
                continue;
            }
            if (Contains(trimmed, "BPEN =")) {
                rmd.contact.bpen = values[0] * kMmToM;
            } else if (Contains(trimmed, "MAXPEN =")) {
                rmd.contact.max_pen = values[0] * kMmToM;
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

static std::vector<ReferenceRow> LoadGear22Reference(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open reference CSV: " + path.string());
    }

    std::string line;
    std::getline(in, line);
    std::vector<ReferenceRow> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        auto cols = SplitCsvLine(line);
        if (cols.size() < 10) {
            continue;
        }
        ReferenceRow row;
        row.time = std::stod(cols[0]);
        row.omega_rx = std::stod(cols[6]);
        row.alpha_rx = std::stod(cols[9]);
        rows.push_back(row);
    }
    return rows;
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
            mesh.vertices.emplace_back(body_local_mm * kMmToM);
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

static ChVector3d AngularVelocityXCross(double omega, const ChVector3d& r) {
    return ChVector3d(0, -omega * r.z(), omega * r.y());
}

static double SampleSDF(SdfSampler& sampler, const ChVector3d& p) {
    return static_cast<double>(sampler.wsSample(openvdb::Vec3d(p.x(), p.y(), p.z())));
}

static ChVector3d SDFGradient(const SparseSDF& sdf, SdfSampler& sampler, const ChVector3d& p) {
    const double h = 0.5 * sdf.voxel_size;
    ChVector3d grad((SampleSDF(sampler, p + ChVector3d(h, 0, 0)) -
                         SampleSDF(sampler, p - ChVector3d(h, 0, 0))) /
                        (2.0 * h),
                    (SampleSDF(sampler, p + ChVector3d(0, h, 0)) -
                         SampleSDF(sampler, p - ChVector3d(0, h, 0))) /
                        (2.0 * h),
                    (SampleSDF(sampler, p + ChVector3d(0, 0, h)) -
                         SampleSDF(sampler, p - ChVector3d(0, 0, h))) /
                        (2.0 * h));
    return SafeNormalize(grad, ChVector3d(0, 0, 1));
}

static Mat3 BodyRotation(const GearPose& pose, double theta_rx) {
    return Multiply(RotX(theta_rx), pose.initial_rotation);
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
        double phi = SampleSDF(sampler, target_local);
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
        result.torque_on_surface += r_surface_world.Cross(dforce);
        result.force_on_target += dforce * -1.0;
        result.torque_on_target += r_target_world.Cross(dforce * -1.0);
        result.max_effective_penetration = std::max(result.max_effective_penetration, effective_penetration);
        result.active_samples++;
        active_indices.push_back(sample.id);
    }

    result.patch_count = static_cast<int>(surface_graph.FindConnectedComponents(active_indices).size());
    return result;
}

static double Rms(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double v : values) {
        sum += v * v;
    }
    return std::sqrt(sum / static_cast<double>(values.size()));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto wall_start = std::chrono::steady_clock::now();
        openvdb::initialize();
        const std::string root = GetProjectRoot();
        const auto asset_dir = std::filesystem::path(root) / "assets" / "simple_gear";

        std::string out_name = "simple_gear_openvdb";
        double duration_limit = std::numeric_limits<double>::quiet_NaN();
        double h = 5.0e-6;
        double pressure_override = std::numeric_limits<double>::quiet_NaN();
        double damping_override = std::numeric_limits<double>::quiet_NaN();

        auto parse_arg_value = [](const std::string& arg, const std::string& key) -> std::string {
            const std::string prefix = key + "=";
            if (arg.rfind(prefix, 0) == 0) {
                return arg.substr(prefix.size());
            }
            return "";
        };

        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            std::string value = parse_arg_value(arg, "--out-name");
            if (!value.empty()) {
                out_name = value;
                continue;
            }
            value = parse_arg_value(arg, "--duration");
            if (!value.empty()) {
                duration_limit = std::stod(value);
                continue;
            }
            value = parse_arg_value(arg, "--dt");
            if (!value.empty()) {
                h = std::stod(value);
                continue;
            }
            value = parse_arg_value(arg, "--pressure-at-bpen");
            if (!value.empty()) {
                pressure_override = std::stod(value);
                continue;
            }
            value = parse_arg_value(arg, "--damping-pressure");
            if (!value.empty()) {
                damping_override = std::stod(value);
            }
        }

        const auto out_dir = std::filesystem::path(root) / "out" / out_name;
        std::filesystem::create_directories(out_dir);

        SimpleGearRmd rmd = LoadSimpleGearRmd(asset_dir / "simple gear.rmd");
        if (!std::isnan(pressure_override)) {
            rmd.contact.pressure_at_bpen = pressure_override;
        }
        if (!std::isnan(damping_override)) {
            rmd.contact.damping_pressure = damping_override;
        }

        GearPose pose21{rmd.gear21.cm_marker_mm * kMmToM, rmd.gear21.part_rotation};
        GearPose pose22{rmd.gear22.cm_marker_mm * kMmToM, rmd.gear22.part_rotation};

        Mesh gear21 = LoadObjAsBodyLocal(asset_dir / "gear_21.obj",
                                         rmd.gear21.surface_ref_marker_mm,
                                         rmd.gear21.surface_ref_rotation,
                                         rmd.gear21.cm_marker_mm,
                                         rmd.gear21.part_rotation);
        Mesh gear22 = LoadObjAsBodyLocal(asset_dir / "gear_22.obj",
                                         rmd.gear22.surface_ref_marker_mm,
                                         rmd.gear22.surface_ref_rotation,
                                         rmd.gear22.cm_marker_mm,
                                         rmd.gear22.part_rotation);
        SparseSDF gear21_sdf = BuildSDF(gear21, 2.5e-5, 20.0f);
        SparseSDF gear22_sdf = BuildSDF(gear22, 2.5e-5, 20.0f);
        SurfaceGraph gear21_graph = MakeTriangleMeshSurfaceGraph(gear21.vertices, gear21.faces);
        SurfaceGraph gear22_graph = MakeTriangleMeshSurfaceGraph(gear22.vertices, gear22.faces);
        auto reference = LoadGear22Reference(asset_dir / "data" / "Gear22.csv");
        if (reference.empty()) {
            throw std::runtime_error("Gear22 reference is empty");
        }

        const double omega21 = 1.0;
        double theta21 = 0.0;
        double theta22 = 0.0;
        double omega22 = 0.0;
        double time = 0.0;
        const double run_until = std::isnan(duration_limit) ? reference.back().time
                                                            : std::min(duration_limit, reference.back().time);

        std::ofstream frames(out_dir / "gear22_rx_velocity.csv");
        frames << std::setprecision(16);
        frames << "time,omega22_rx,reference_omega22_rx,omega_error,alpha22_rx,reference_alpha22_rx,"
                  "patch_count,active_samples,min_phi,max_effective_penetration,torque_x,"
                  "torque22_from_gear22_surface,torque22_from_gear21_surface\n";

        std::vector<double> omega_errors;
        std::vector<double> omega_errors_after_startup;
        double max_abs_error = 0.0;
        double max_abs_error_after_startup = 0.0;
        int max_patch_count = 0;
        int max_active_samples = 0;
        double min_phi_seen = std::numeric_limits<double>::max();
        double max_effective_penetration_seen = 0.0;
        double last_alpha = 0.0;
        double last_torque = 0.0;
        double last_torque_from_gear22_surface = 0.0;
        double last_torque_from_gear21_surface = 0.0;
        int last_patch_count = 0;
        int last_active_samples = 0;
        double last_min_phi = 0.0;
        double last_max_effective_penetration = 0.0;
        ReferenceRow last_reference = reference.front();
        size_t rows_written = 0;

        const auto integration_start = std::chrono::steady_clock::now();
        for (const auto& ref : reference) {
            if (ref.time > run_until + 1.0e-14) {
                break;
            }
            while (time < ref.time - 1.0e-14) {
                double step = std::min(h, ref.time - time);
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

                last_torque_from_gear22_surface = gear22_on_gear21.torque_on_surface.x();
                last_torque_from_gear21_surface = gear21_on_gear22.torque_on_target.x();
                last_torque = 0.5 * (last_torque_from_gear22_surface + last_torque_from_gear21_surface);
                last_alpha = last_torque / rmd.gear22.inertia_x_kg_m2;
                last_patch_count = gear22_on_gear21.patch_count + gear21_on_gear22.patch_count;
                last_active_samples = gear22_on_gear21.active_samples + gear21_on_gear22.active_samples;
                last_min_phi = std::min(gear22_on_gear21.min_phi, gear21_on_gear22.min_phi);
                last_max_effective_penetration = std::max(gear22_on_gear21.max_effective_penetration,
                                                          gear21_on_gear22.max_effective_penetration);
                max_patch_count = std::max(max_patch_count, last_patch_count);
                max_active_samples = std::max(max_active_samples, last_active_samples);
                min_phi_seen = std::min(min_phi_seen, last_min_phi);
                max_effective_penetration_seen =
                    std::max(max_effective_penetration_seen, last_max_effective_penetration);

                omega22 += last_alpha * step;
                theta22 += omega22 * step;
                theta21 += omega21 * step;
                time += step;
            }

            double error = omega22 - ref.omega_rx;
            omega_errors.push_back(error);
            max_abs_error = std::max(max_abs_error, std::abs(error));
            if (ref.time >= 0.02) {
                omega_errors_after_startup.push_back(error);
                max_abs_error_after_startup = std::max(max_abs_error_after_startup, std::abs(error));
            }
            frames << ref.time << "," << omega22 << "," << ref.omega_rx << "," << error << ","
                   << last_alpha << "," << ref.alpha_rx << "," << last_patch_count << ","
                   << last_active_samples << "," << last_min_phi << "," << last_max_effective_penetration << ","
                   << last_torque << "," << last_torque_from_gear22_surface << ","
                   << last_torque_from_gear21_surface << "\n";
            last_reference = ref;
            rows_written++;
        }
        const auto integration_end = std::chrono::steady_clock::now();
        if (rows_written == 0) {
            throw std::runtime_error("No reference rows written for requested duration");
        }
        const double integration_seconds =
            std::chrono::duration<double>(integration_end - integration_start).count();
        const double wall_seconds = std::chrono::duration<double>(integration_end - wall_start).count();

        std::ofstream summary(out_dir / "summary.csv");
        summary << "quantity,value\n";
        summary << "duration," << last_reference.time << "\n";
        summary << "dt," << h << "\n";
        summary << "rows_written," << rows_written << "\n";
        summary << "gear21_vertices," << gear21.vertices.size() << "\n";
        summary << "gear22_vertices," << gear22.vertices.size() << "\n";
        summary << "rmd_bpen_m," << rmd.contact.bpen << "\n";
        summary << "rmd_maxpen_m," << rmd.contact.max_pen << "\n";
        summary << "rmd_korder," << rmd.contact.korder << "\n";
        summary << "rmd_k," << rmd.contact.recurdyn_k << "\n";
        summary << "rmd_c," << rmd.contact.recurdyn_c << "\n";
        summary << "pressure_at_bpen," << rmd.contact.pressure_at_bpen << "\n";
        summary << "damping_pressure," << rmd.contact.damping_pressure << "\n";
        summary << "omega_rx_max_abs_error," << max_abs_error << "\n";
        summary << "omega_rx_rms_error," << Rms(omega_errors) << "\n";
        summary << "omega_rx_max_abs_error_t_ge_0p02," << max_abs_error_after_startup << "\n";
        summary << "omega_rx_rms_error_t_ge_0p02," << Rms(omega_errors_after_startup) << "\n";
        summary << "omega_rx_final," << omega22 << "\n";
        summary << "omega_rx_reference_final," << last_reference.omega_rx << "\n";
        summary << "max_patch_count," << max_patch_count << "\n";
        summary << "max_active_samples," << max_active_samples << "\n";
        summary << "min_phi_seen," << min_phi_seen << "\n";
        summary << "max_effective_penetration_seen," << max_effective_penetration_seen << "\n";
        summary << "integration_wall_time_seconds," << integration_seconds << "\n";
        summary << "wall_time_seconds_before_summary," << wall_seconds << "\n";

        std::cout << "Wrote " << (out_dir / "gear22_rx_velocity.csv").string() << std::endl;
        std::cout << "Wrote " << (out_dir / "summary.csv").string() << std::endl;
        std::cout << "Final omega22_rx=" << omega22
                  << " reference=" << last_reference.omega_rx
                  << " rms_error=" << Rms(omega_errors) << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 2;
    }
}
