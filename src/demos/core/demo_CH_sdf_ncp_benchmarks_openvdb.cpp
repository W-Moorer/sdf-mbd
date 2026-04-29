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
// OpenVDB-backed SDF-NCP benchmarks using full asset geometry.
//
// This executable shares the mesh/OpenVDB SDF frontend with field-contact style
// demos, but the backend is SDF-NCP: gaps and normals are converted into
// Fischer-Burmeister complementarity residuals and local contact multipliers.
// It does not call pressure-field patch extraction or pressure integration.
//
// =============================================================================

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "chrono/collision/ChSdfNcpContact.h"
#include "chrono/collision/ChSdfOpenVDB.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChTypes.h"
#include "chrono/core/ChVector3.h"
#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChBodyAuxRef.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChLinkRevolute.h"
#include "chrono/physics/ChSdfNcpConstraintContact.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolver.h"
#include "chrono/timestepper/ChTimestepperHHT.h"
#include "chrono/timestepper/ChTimestepperImplicit.h"

#include <openvdb/openvdb.h>

using chrono::ChVector3d;
using chrono::ChBody;
using chrono::ChBodyAuxRef;
using chrono::ChFrame;
using chrono::ChFunctionConst;
using chrono::ChLinkMateGeneric;
using chrono::ChLinkMotorRotationSpeed;
using chrono::ChLinkRevolute;
using chrono::ChMatrix33d;
using chrono::ChSolver;
using chrono::ChSystemNSC;
using chrono::ChTimestepper;
using chrono::ChTimestepperHHT;
using chrono::ChTimestepperImplicit;
using chrono::ChVectorDynamic;
using chrono::UpdateFlags;
using chrono::sdfcontact::BuildOpenVdbLevelSetFromMesh;
using chrono::sdfcontact::ChOpenVdbSdfGrid;
using chrono::sdfcontact::ChSdfAabb;
using chrono::sdfcontact::ChSdfContactSampleQuery;
using chrono::sdfcontact::ChSdfContactSampleBvh;
using chrono::sdfcontact::ChSdfContactSurfaceGraph;
using chrono::sdfcontact::ChSdfTriangleMeshData;
using chrono::sdfcontact::LoadWavefrontMeshForSdf;
using chrono::sdfcontact::MakeSurfaceGraphFromMeshForSdf;
using namespace chrono::sdfncp;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct CamConfig {
    std::string case_name = "cam";
    std::string asset = "assets/cam";
    double dt = 0.002;
    double t_end = 0.25;
    double eps = 1.0e-7;
    double tolerance = 1.0e-9;
    int max_iterations = 40;
    double voxel_size = 0.002;
    float half_width_voxels = 32.0f;
    double mass_follower = 6.19811156388412;
    double gravity_y = -9.81;
    double cam_motor_speed = 3.1415926;
    ChVector3d cam_cm0 = ChVector3d(0.017329438212968, 0.00390749251884395, 0.0499930593952011);
    ChVector3d cam_joint_abs = ChVector3d(0, 0, 0.05);
    ChVector3d follower_joint_abs = ChVector3d(0, 0.4, 0);
    ChVector3d cam_surface_marker = ChVector3d(-0.320462025491936, -0.325583626629385, -0.116);
    ChVector3d follower_surface_marker = ChVector3d(-0.1662, 0.2338, -0.0662);
    ChVector3d follower_cm0 = ChVector3d(0, 0.4, 0);
};

struct BenchmarkStats {
    std::string case_name;
    std::string method = "sdf_ncp";
    std::string asset;
    double dt = 0.0;
    double t_end = 0.0;
    int num_steps = 0;
    int num_contacts = 0;
    double max_penetration = 0.0;
    double sum_penetration = 0.0;
    double max_lambda_n = 0.0;
    double max_ncp_residual = 0.0;
    double sum_ncp_residual = 0.0;
    double max_complementarity_error = 0.0;
    double sum_complementarity_error = 0.0;
    double sum_iterations = 0.0;
    double sum_success = 0.0;
    int samples = 0;
    double runtime_seconds = 0.0;
    std::string notes;
};

struct CamReferenceRow {
    double time = 0.0;
    double follower_y = 0.0;
    double follower_vy = 0.0;
};

std::filesystem::path GetProjectRoot() {
    auto path = std::filesystem::current_path();
    for (int i = 0; i < 8; i++) {
        if (std::filesystem::exists(path / "src") && std::filesystem::exists(path / "assets")) {
            return path;
        }
        if (!path.has_parent_path() || path == path.parent_path()) {
            break;
        }
        path = path.parent_path();
    }
    return std::filesystem::current_path();
}

ChVector3d RotateZ(const ChVector3d& v, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return ChVector3d(c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z());
}

ChVector3d OmegaZCross(double omega, const ChVector3d& r) {
    return ChVector3d(-omega * r.y(), omega * r.x(), 0.0);
}

struct Mat3 {
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    ChVector3d operator*(const ChVector3d& v) const {
        return ChVector3d(m[0][0] * v.x() + m[0][1] * v.y() + m[0][2] * v.z(),
                          m[1][0] * v.x() + m[1][1] * v.y() + m[1][2] * v.z(),
                          m[2][0] * v.x() + m[2][1] * v.y() + m[2][2] * v.z());
    }
};

Mat3 Multiply(const Mat3& a, const Mat3& b) {
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

Mat3 Transpose(const Mat3& a) {
    Mat3 out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out.m[i][j] = a.m[j][i];
        }
    }
    return out;
}

Mat3 RotXMatrix(double angle) {
    Mat3 out;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    out.m[1][1] = c;
    out.m[1][2] = -s;
    out.m[2][1] = s;
    out.m[2][2] = c;
    return out;
}

Mat3 RotYMatrix(double angle) {
    Mat3 out;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    out.m[0][0] = c;
    out.m[0][2] = s;
    out.m[2][0] = -s;
    out.m[2][2] = c;
    return out;
}

Mat3 RotZMatrix(double angle) {
    Mat3 out;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    out.m[0][0] = c;
    out.m[0][1] = -s;
    out.m[1][0] = s;
    out.m[1][1] = c;
    return out;
}

Mat3 RecurDynEuler(double rx, double ry, double rz) {
    return Multiply(RotZMatrix(rz), Multiply(RotYMatrix(ry), RotXMatrix(rx)));
}

ChMatrix33d ToChronoMatrix(const Mat3& m) {
    ChMatrix33d out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out(i, j) = m.m[i][j];
        }
    }
    return out;
}

ChVector3d AngularVelocityXCross(double omega, const ChVector3d& r) {
    return ChVector3d(0.0, -omega * r.z(), omega * r.y());
}

ChVector3d UnitXCross(const ChVector3d& r) {
    return ChVector3d(0.0, -r.z(), r.y());
}

std::vector<ChVector3d> AabbCorners(const ChSdfAabb& bounds) {
    if (!bounds.IsValid()) {
        return {};
    }
    std::vector<ChVector3d> corners;
    corners.reserve(8);
    for (int ix = 0; ix < 2; ix++) {
        for (int iy = 0; iy < 2; iy++) {
            for (int iz = 0; iz < 2; iz++) {
                corners.emplace_back(ix == 0 ? bounds.min.x() : bounds.max.x(),
                                     iy == 0 ? bounds.min.y() : bounds.max.y(),
                                     iz == 0 ? bounds.min.z() : bounds.max.z());
            }
        }
    }
    return corners;
}

ChSdfAabb ExpandedAabb(ChSdfAabb bounds, double padding) {
    bounds.Expand(padding);
    return bounds;
}

ChSdfAabb TransformAabbBetweenBodyFrames(const ChSdfAabb& target_local_bounds,
                                         const Mat3& target_rotation,
                                         const ChVector3d& target_center,
                                         const Mat3& source_rotation,
                                         const ChVector3d& source_center,
                                         double padding) {
    const ChSdfAabb padded = ExpandedAabb(target_local_bounds, padding);
    ChSdfAabb source_bounds;
    if (!padded.IsValid()) {
        return source_bounds;
    }

    const Mat3 source_rotation_t = Transpose(source_rotation);
    for (const auto& target_corner : AabbCorners(padded)) {
        const ChVector3d world = target_center + target_rotation * target_corner;
        const ChVector3d source_local = source_rotation_t * (world - source_center);
        source_bounds.Include(source_local);
    }
    return source_bounds;
}

ChSdfAabb CamLocalAabbInFollowerLocal(const ChSdfAabb& cam_local_bounds, double follower_y, double theta_cam, double padding) {
    const ChSdfAabb padded = ExpandedAabb(cam_local_bounds, padding);
    ChSdfAabb follower_local_bounds;
    if (!padded.IsValid()) {
        return follower_local_bounds;
    }

    for (const auto& cam_corner : AabbCorners(padded)) {
        const ChVector3d world = RotateZ(cam_corner, theta_cam);
        follower_local_bounds.Include(world - ChVector3d(0, follower_y, 0));
    }
    return follower_local_bounds;
}

double NormDynamic(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

bool SolveLinearDynamic(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    x.assign(n, 0.0);

    for (int col = 0; col < n; col++) {
        int pivot = col;
        double pivot_abs = std::abs(a[col][col]);
        for (int row = col + 1; row < n; row++) {
            const double value = std::abs(a[row][col]);
            if (value > pivot_abs) {
                pivot = row;
                pivot_abs = value;
            }
        }
        if (pivot_abs < 1.0e-14) {
            return false;
        }
        if (pivot != col) {
            std::swap(a[pivot], a[col]);
            std::swap(b[pivot], b[col]);
        }

        const double diag = a[col][col];
        for (int j = col; j < n; j++) {
            a[col][j] /= diag;
        }
        b[col] /= diag;

        for (int row = 0; row < n; row++) {
            if (row == col) {
                continue;
            }
            const double factor = a[row][col];
            for (int j = col; j < n; j++) {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    x = b;
    return true;
}

std::vector<std::vector<double>> FiniteDifferenceJacobian(
    const std::function<std::vector<double>(const std::vector<double>&)>& residual,
    const std::vector<double>& z,
    double h = 1.0e-7) {
    const size_t n = z.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));
    for (size_t col = 0; col < n; col++) {
        std::vector<double> zp = z;
        std::vector<double> zm = z;
        zp[col] += h;
        zm[col] -= h;
        const auto rp = residual(zp);
        const auto rm = residual(zm);
        for (size_t row = 0; row < n; row++) {
            jac[row][col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }
    return jac;
}

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::string token;
    std::istringstream stream(line);
    while (std::getline(stream, token, ',')) {
        out.push_back(token);
    }
    return out;
}

std::string ReadTextFile(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open file: " + path.string());
    }
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

std::string Trim(std::string text) {
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

bool Contains(const std::string& text, const std::string& pattern) {
    return text.find(pattern) != std::string::npos;
}

bool HasRmdField(const std::string& line, const std::string& field) {
    const size_t pos = line.find(field);
    if (pos == std::string::npos) {
        return false;
    }
    const bool left_ok =
        pos == 0 || line[pos - 1] == ',' || std::isspace(static_cast<unsigned char>(line[pos - 1])) != 0;
    size_t end = pos + field.size();
    while (end < line.size() && std::isspace(static_cast<unsigned char>(line[end])) != 0) {
        end++;
    }
    const bool right_ok = end < line.size() && line[end] == '=';
    return left_ok && right_ok;
}

double RecurDynStepRamp(double x, double x0, double y0, double x1, double y1) {
    if (!(x1 > x0)) {
        return x >= x1 ? y1 : y0;
    }
    if (x <= x0) {
        return y0;
    }
    if (x >= x1) {
        return y1;
    }
    const double u = (x - x0) / (x1 - x0);
    const double s = u * u * (3.0 - 2.0 * u);
    return y0 + (y1 - y0) * s;
}

std::string ExtractQuotedName(const std::string& line) {
    const size_t first = line.find('\'');
    if (first == std::string::npos) {
        return "";
    }
    const size_t second = line.find('\'', first + 1);
    return second == std::string::npos ? "" : line.substr(first + 1, second - first - 1);
}

std::string ExtractRmdValueText(const std::string& line) {
    const size_t eq = line.find('=');
    if (eq == std::string::npos) {
        return "";
    }
    std::string value = Trim(line.substr(eq + 1));
    while (!value.empty() && (value.back() == ',' || value.back() == '\\')) {
        value.pop_back();
        value = Trim(value);
    }
    return value;
}

std::string CleanRmdNumberToken(std::string token) {
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

std::vector<double> ParseNumbersAfterEquals(const std::string& line) {
    const size_t eq = line.find('=');
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

std::vector<double> ParseNumbersInText(std::string text) {
    for (char& c : text) {
        if (c == ',' || c == '=' || c == '\\') {
            c = ' ';
        }
    }
    std::istringstream stream(text);
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

int ParseRmdEntityId(const std::string& line) {
    const size_t slash = line.find('/');
    if (slash == std::string::npos) {
        return -1;
    }
    const size_t comma = line.find(',', slash + 1);
    std::string token = line.substr(slash + 1, comma == std::string::npos ? std::string::npos : comma - slash - 1);
    token = CleanRmdNumberToken(token);
    if (token.empty()) {
        return -1;
    }
    try {
        return std::stoi(token);
    } catch (...) {
        return -1;
    }
}

ChVector3d ParseRmdTriple(const std::string& line) {
    const auto values = ParseNumbersAfterEquals(line);
    if (values.size() < 3) {
        throw std::runtime_error("Expected RMD triple in line: " + line);
    }
    return ChVector3d(values[0], values[1], values[2]);
}

double ExtractJsonDouble(const std::string& text, const std::string& key, double fallback) {
    const std::string quoted = "\"" + key + "\"";
    const size_t key_pos = text.find(quoted);
    if (key_pos == std::string::npos) {
        return fallback;
    }
    const size_t colon = text.find(':', key_pos + quoted.size());
    if (colon == std::string::npos) {
        return fallback;
    }
    size_t begin = text.find_first_of("-0123456789.", colon + 1);
    if (begin == std::string::npos) {
        return fallback;
    }
    size_t end = begin;
    while (end < text.size()) {
        const char c = text[end];
        if ((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E') {
            end++;
        } else {
            break;
        }
    }
    return std::stod(text.substr(begin, end - begin));
}

ChVector3d ExtractJsonVector3(const std::string& text, const std::string& key, const ChVector3d& fallback) {
    const std::string quoted = "\"" + key + "\"";
    const size_t key_pos = text.find(quoted);
    if (key_pos == std::string::npos) {
        return fallback;
    }
    const size_t open = text.find('[', key_pos + quoted.size());
    const size_t close = text.find(']', open);
    if (open == std::string::npos || close == std::string::npos || close <= open) {
        return fallback;
    }
    std::string body = text.substr(open + 1, close - open - 1);
    std::replace(body.begin(), body.end(), ',', ' ');
    std::istringstream stream(body);
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!(stream >> x >> y >> z)) {
        return fallback;
    }
    return ChVector3d(x, y, z);
}

std::vector<CamReferenceRow> LoadCamReference(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        return {};
    }

    std::string line;
    std::getline(in, line);
    std::vector<CamReferenceRow> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto cols = SplitCsvLine(line);
        if (cols.size() < 14) {
            continue;
        }
        CamReferenceRow row;
        row.time = std::stod(cols[0]);
        row.follower_y = std::stod(cols[10]);
        row.follower_vy = std::stod(cols[13]);
        rows.push_back(row);
    }
    return rows;
}

CamReferenceRow FindNearestReference(const std::vector<CamReferenceRow>& rows, double time) {
    if (rows.empty()) {
        return {};
    }
    return *std::min_element(rows.begin(), rows.end(), [time](const CamReferenceRow& a, const CamReferenceRow& b) {
        return std::abs(a.time - time) < std::abs(b.time - time);
    });
}

double LoadFirstMotionFunctionConstant(const std::filesystem::path& path, double fallback) {
    std::ifstream in(path);
    if (!in) {
        return fallback;
    }

    bool in_motion = false;
    std::string line;
    while (std::getline(in, line)) {
        const std::string trimmed = Trim(line);
        if (trimmed.rfind("MOTION /", 0) == 0) {
            in_motion = true;
            continue;
        }
        if (in_motion && (trimmed.rfind("JOINT /", 0) == 0 || trimmed.rfind("SOLIDCONTACT /", 0) == 0 ||
                          trimmed.rfind("GGEOMCONTACT /", 0) == 0 || trimmed.rfind("EXPRESSION /", 0) == 0)) {
            in_motion = false;
        }
        if (in_motion && Contains(trimmed, "FUNCTION =")) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                return values.front();
            }
        }
    }
    return fallback;
}

struct RmdPart {
    int id = -1;
    std::string name;
    double mass = 1.0;
    ChVector3d qg = ChVector3d(0, 0, 0);
    Mat3 rotation;
    ChVector3d inertia_xx = ChVector3d(1.0e-3, 1.0e-3, 1.0e-3);
    ChVector3d inertia_xy = ChVector3d(0, 0, 0);
    int cm_marker_id = -1;
    bool has_qg = false;
    bool has_rotation = false;
    bool has_inertia = false;
};

struct RmdMarker {
    int id = -1;
    std::string name;
    int part_id = -1;
    ChVector3d qp = ChVector3d(0, 0, 0);
    Mat3 rotation;
    bool has_qp = false;
    bool has_rotation = false;
};

struct RmdJoint {
    int id = -1;
    std::string name;
    int i_marker_id = -1;
    int j_marker_id = -1;
    std::string type;
};

struct RmdMotion {
    int id = -1;
    std::string name;
    int joint_id = -1;
    std::string coordinate_type;
    std::string input_type;
    double function_constant = 0.0;
    bool has_function = false;
};

struct RmdSurface {
    int id = -1;
    std::string name;
    int ref_marker_id = -1;
};

struct RmdSolidContact {
    int id = -1;
    std::string name;
    int i_float_marker_id = -1;
    int j_float_marker_id = -1;
    int i_surface_id = -1;
    int j_surface_id = -1;
    double bpen = 0.0;
    double rdf = 0.0;
    double stiffness = 0.0;
    double damping = 0.0;
    double korder = 0.0;
    double dynamic_friction_coefficient = 0.0;
    double static_transition_velocity = 0.0;
    double dynamic_transition_velocity = 0.0;
    double static_friction_coefficient = 0.0;
};

struct RmdModel {
    std::map<int, RmdPart> parts;
    std::map<int, RmdMarker> markers;
    std::map<int, RmdJoint> joints;
    std::map<int, RmdMotion> motions;
    std::map<int, RmdSurface> surfaces;
    std::map<int, RmdSolidContact> solid_contacts;
    ChVector3d gravity = ChVector3d(0, -9.81, 0);
};

const RmdPart& RequirePartByName(const RmdModel& model, const std::string& name) {
    for (const auto& item : model.parts) {
        if (item.second.name == name) {
            return item.second;
        }
    }
    throw std::runtime_error("RMD part not found: " + name);
}

const RmdMarker& RequireMarkerByName(const RmdModel& model, const std::string& name) {
    for (const auto& item : model.markers) {
        if (item.second.name == name) {
            return item.second;
        }
    }
    throw std::runtime_error("RMD marker not found: " + name);
}

const RmdJoint& RequireJointByName(const RmdModel& model, const std::string& name) {
    for (const auto& item : model.joints) {
        if (item.second.name == name) {
            return item.second;
        }
    }
    throw std::runtime_error("RMD joint not found: " + name);
}

const RmdMotion& RequireMotionByName(const RmdModel& model, const std::string& name) {
    for (const auto& item : model.motions) {
        if (item.second.name == name) {
            return item.second;
        }
    }
    throw std::runtime_error("RMD motion not found: " + name);
}

const RmdSolidContact& RequireSolidContactByName(const RmdModel& model, const std::string& name) {
    for (const auto& item : model.solid_contacts) {
        if (item.second.name == name) {
            return item.second;
        }
    }
    throw std::runtime_error("RMD solid contact not found: " + name);
}

double MapRecurDynRotationVelocityToChronoMotorSpeed(const RmdMotion& motion) {
    // RecurDyn rotational joint motion is defined for the J marker relative to
    // the I marker.  The Chrono motor in this benchmark is initialized as
    // (body_J, body_I), whose positive spindle velocity has the opposite sign
    // for the same physical rotation.
    return -motion.function_constant;
}

const RmdPart& RequirePartById(const RmdModel& model, int id) {
    const auto found = model.parts.find(id);
    if (found == model.parts.end()) {
        throw std::runtime_error("RMD part id not found: " + std::to_string(id));
    }
    return found->second;
}

const RmdMarker& RequireMarkerById(const RmdModel& model, int id) {
    const auto found = model.markers.find(id);
    if (found == model.markers.end()) {
        throw std::runtime_error("RMD marker id not found: " + std::to_string(id));
    }
    return found->second;
}

const RmdSurface& RequireSurfaceById(const RmdModel& model, int id) {
    const auto found = model.surfaces.find(id);
    if (found == model.surfaces.end()) {
        throw std::runtime_error("RMD surface id not found: " + std::to_string(id));
    }
    return found->second;
}

ChFrame<> RmdMarkerFrameAbs(const RmdModel& model, const RmdMarker& marker) {
    const RmdPart& part = RequirePartById(model, marker.part_id);
    const ChVector3d pos = part.qg + part.rotation * marker.qp;
    const Mat3 rot = Multiply(part.rotation, marker.rotation);
    return ChFrame<>(pos, ToChronoMatrix(rot));
}

ChFrame<> RmdPartRefFrameAbs(const RmdPart& part) {
    return ChFrame<>(part.qg, ToChronoMatrix(part.rotation));
}

void ConfigureAuxRefBodyFromRmdPart(const RmdModel& model,
                                    const RmdPart& part,
                                    const RmdMarker& cm_marker,
                                    const std::shared_ptr<ChBodyAuxRef>& body) {
    body->SetName(part.name);
    body->SetMass(part.mass);
    body->SetInertiaXX(part.inertia_xx);
    body->SetInertiaXY(part.inertia_xy);
    body->SetFrameCOMToRef(ChFrame<>(cm_marker.qp, ToChronoMatrix(cm_marker.rotation)));
    body->SetFrameRefToAbs(RmdPartRefFrameAbs(part));
}

RmdModel LoadRecurDynRmdModel(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open RMD: " + path.string());
    }

    enum class Block { None, Part, Marker, Joint, Motion, Surface, SolidContact };
    Block block = Block::None;
    RmdModel model;
    int current_id = -1;
    int pending_ip_part_id = -1;
    std::string line;
    while (std::getline(in, line)) {
        const std::string trimmed = Trim(line);
        if (trimmed.empty()) {
            continue;
        }

        auto start_block = [&](Block next, int id) {
            block = next;
            current_id = id;
            pending_ip_part_id = -1;
        };

        if (trimmed.rfind("PART /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.parts[id].id = id;
            start_block(Block::Part, id);
            continue;
        }
        if (trimmed.rfind("MARKER /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.markers[id].id = id;
            start_block(Block::Marker, id);
            continue;
        }
        if (trimmed.rfind("JOINT /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.joints[id].id = id;
            start_block(Block::Joint, id);
            continue;
        }
        if (trimmed.rfind("MOTION /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.motions[id].id = id;
            start_block(Block::Motion, id);
            continue;
        }
        if (trimmed.rfind("CSURFACE /", 0) == 0 || trimmed.rfind("GGEOM /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.surfaces[id].id = id;
            start_block(Block::Surface, id);
            continue;
        }
        if (trimmed.rfind("SOLIDCONTACT /", 0) == 0 || trimmed.rfind("GGEOMCONTACT /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.solid_contacts[id].id = id;
            start_block(Block::SolidContact, id);
            continue;
        }

        if (block == Block::Part && pending_ip_part_id >= 0 && !Contains(trimmed, "=")) {
            const auto values = ParseNumbersInText(trimmed);
            if (values.size() >= 3) {
                auto& part = model.parts[pending_ip_part_id];
                part.inertia_xy = ChVector3d(values[0], values[1], values[2]);
            }
            pending_ip_part_id = -1;
            continue;
        }

        if (trimmed.rfind("accgrav /", 0) == 0 || Contains(trimmed, "IGRAV =")) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                model.gravity.x() = values.front();
            }
        } else if (Contains(trimmed, "JGRAV =")) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                model.gravity.y() = values.front();
            }
        } else if (Contains(trimmed, "KGRAV =")) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                model.gravity.z() = values.front();
            }
        }

        if (trimmed.rfind("EXPRESSION /", 0) == 0 || trimmed.rfind("UNITS /", 0) == 0 ||
            trimmed.rfind("OUTPUT/", 0) == 0 || trimmed.rfind("INTPAR /", 0) == 0 ||
            trimmed.rfind("SOLVEROPTION /", 0) == 0 || trimmed.rfind("INFO /", 0) == 0) {
            block = Block::None;
            current_id = -1;
            pending_ip_part_id = -1;
            continue;
        }

        if (Contains(trimmed, "NAME =")) {
            const std::string name = ExtractQuotedName(trimmed);
            if (block == Block::Part) {
                model.parts[current_id].name = name;
            } else if (block == Block::Marker) {
                model.markers[current_id].name = name;
            } else if (block == Block::Joint) {
                model.joints[current_id].name = name;
            } else if (block == Block::Motion) {
                model.motions[current_id].name = name;
            } else if (block == Block::Surface) {
                model.surfaces[current_id].name = name;
            } else if (block == Block::SolidContact) {
                model.solid_contacts[current_id].name = name;
            }
        }

        if (block == Block::Part) {
            auto& part = model.parts[current_id];
            if (Contains(trimmed, "MASS =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    part.mass = values.front();
                }
            } else if (Contains(trimmed, "CM =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    part.cm_marker_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "QG =")) {
                part.qg = ParseRmdTriple(trimmed);
                part.has_qg = true;
            } else if (Contains(trimmed, "REULER =")) {
                const ChVector3d euler = ParseRmdTriple(trimmed);
                part.rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
                part.has_rotation = true;
            } else if (Contains(trimmed, "IP =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (values.size() >= 3) {
                    part.inertia_xx = ChVector3d(values[0], values[1], values[2]);
                    part.has_inertia = true;
                    pending_ip_part_id = current_id;
                }
            }
        } else if (block == Block::Marker) {
            auto& marker = model.markers[current_id];
            if (Contains(trimmed, "PART =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    marker.part_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "QP =")) {
                marker.qp = ParseRmdTriple(trimmed);
                marker.has_qp = true;
            } else if (Contains(trimmed, "REULER =")) {
                const ChVector3d euler = ParseRmdTriple(trimmed);
                marker.rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
                marker.has_rotation = true;
            }
        } else if (block == Block::Joint) {
            auto& joint = model.joints[current_id];
            if (Contains(trimmed, "I =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    joint.i_marker_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "J =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    joint.j_marker_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "Revolute")) {
                joint.type = "Revolute";
            } else if (Contains(trimmed, "Translational")) {
                joint.type = "Translational";
            }
        } else if (block == Block::Motion) {
            auto& motion = model.motions[current_id];
            if (Contains(trimmed, "JOINT =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    motion.joint_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "ROTATION")) {
                motion.coordinate_type = "ROTATION";
            } else if (Contains(trimmed, "TRANSLATION")) {
                motion.coordinate_type = "TRANSLATION";
            } else if (Contains(trimmed, "VELOCITY")) {
                motion.input_type = "VELOCITY";
            } else if (Contains(trimmed, "FUNCTION =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    motion.function_constant = values.front();
                    motion.has_function = true;
                }
            }
        } else if (block == Block::Surface) {
            auto& surface = model.surfaces[current_id];
            if (Contains(trimmed, "RM =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    surface.ref_marker_id = static_cast<int>(std::round(values.front()));
                }
            }
        } else if (block == Block::SolidContact) {
            auto& contact = model.solid_contacts[current_id];
            if (Contains(trimmed, "IFLOAT =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.i_float_marker_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "JFLOAT =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.j_float_marker_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "ICSURFACEID =") || Contains(trimmed, "IGGEOMID =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.i_surface_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "JCSURFACEID =") || Contains(trimmed, "JGGEOMID =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.j_surface_id = static_cast<int>(std::round(values.front()));
                }
            } else if (HasRmdField(trimmed, "BPEN")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.bpen = values.front();
                }
            } else if (HasRmdField(trimmed, "RDF")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.rdf = values.front();
                }
            } else if (HasRmdField(trimmed, "KORDER")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.korder = values.front();
                }
            } else if (HasRmdField(trimmed, "K")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.stiffness = values.front();
                }
            } else if (HasRmdField(trimmed, "C")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.damping = values.front();
                }
            } else if (HasRmdField(trimmed, "D_F_C")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.dynamic_friction_coefficient = values.front();
                }
            } else if (HasRmdField(trimmed, "S_T_V")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.static_transition_velocity = values.front();
                }
            } else if (HasRmdField(trimmed, "D_T_V")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.dynamic_transition_velocity = values.front();
                }
            } else if (HasRmdField(trimmed, "S_F_C")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.static_friction_coefficient = values.front();
                }
            }
        }
    }

    return model;
}

void WriteRmdMappingAudit(const std::filesystem::path& path, const RmdModel& model) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "entity_type,id,name,parent_id,x,y,z,rx_note,ry_note,rz_note,value0,value1,value2,notes\n";
    for (const auto& item : model.parts) {
        const auto& part = item.second;
        out << "part," << part.id << "," << part.name << ",-1," << part.qg.x() << "," << part.qg.y() << ","
            << part.qg.z() << ",,,," << part.mass << "," << part.inertia_xx.x() << "," << part.inertia_xx.y()
            << ",qg/mass/inertia parsed from RMD\n";
    }
    for (const auto& item : model.markers) {
        const auto& marker = item.second;
        out << "marker," << marker.id << "," << marker.name << "," << marker.part_id << "," << marker.qp.x() << ","
            << marker.qp.y() << "," << marker.qp.z() << ",,,," << 0 << "," << 0 << "," << 0
            << ",local marker QP parsed from RMD\n";
    }
    for (const auto& item : model.joints) {
        const auto& joint = item.second;
        out << "joint," << joint.id << "," << joint.name << ",-1,0,0,0,,,"
            << "," << joint.i_marker_id << "," << joint.j_marker_id << ",0," << joint.type << "\n";
    }
    for (const auto& item : model.motions) {
        const auto& motion = item.second;
        out << "motion," << motion.id << "," << motion.name << "," << motion.joint_id << ",0,0,0,,,"
            << "," << motion.function_constant << ",0,0," << motion.coordinate_type << "/" << motion.input_type
            << "\n";
    }
    for (const auto& item : model.surfaces) {
        const auto& surface = item.second;
        out << "surface," << surface.id << "," << surface.name << "," << surface.ref_marker_id
            << ",0,0,0,,,,0,0,0,CSURFACE reference marker parsed from RMD\n";
    }
    for (const auto& item : model.solid_contacts) {
        const auto& contact = item.second;
        out << "solid_contact," << contact.id << "," << contact.name << ",-1,0,0,0,,,"
            << "," << contact.stiffness << "," << contact.damping << "," << contact.bpen
            << ",surfaces " << contact.i_surface_id << "/" << contact.j_surface_id << "; markers "
            << contact.i_float_marker_id << "/" << contact.j_float_marker_id << "; KORDER " << contact.korder
            << "; RDF " << contact.rdf << "; D_F_C " << contact.dynamic_friction_coefficient << "; S_T_V "
            << contact.static_transition_velocity << "; D_T_V " << contact.dynamic_transition_velocity << "; S_F_C "
            << contact.static_friction_coefficient << "\n";
    }
}

void WriteTrajectoryHeader(std::ofstream& out) {
    out << "time,body_id,px,py,pz,q0,q1,q2,q3,vx,vy,vz,wx,wy,wz,contact_id,gap,lambda_n,"
           "contact_weight,lambda_force,penetration,ncp_residual,complementarity_error,residual_norm,iterations,"
           "success\n";
}

void WriteSummary(const std::filesystem::path& path, const BenchmarkStats& stats) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "case_name,method,asset,dt,t_end,num_steps,num_contacts,max_penetration,mean_penetration,"
           "max_lambda_n,max_ncp_residual,mean_ncp_residual,max_complementarity_error,"
           "mean_complementarity_error,mean_iterations,success_rate,runtime_seconds,notes\n";
    const double denom = static_cast<double>(std::max(1, stats.samples));
    out << stats.case_name << "," << stats.method << "," << stats.asset << "," << stats.dt << "," << stats.t_end << ","
        << stats.num_steps << "," << stats.num_contacts << "," << stats.max_penetration << ","
        << stats.sum_penetration / denom << "," << stats.max_lambda_n << "," << stats.max_ncp_residual << ","
        << stats.sum_ncp_residual / denom << "," << stats.max_complementarity_error << ","
        << stats.sum_complementarity_error / denom << "," << stats.sum_iterations / denom << ","
        << stats.sum_success / denom << "," << stats.runtime_seconds << "," << stats.notes << "\n";
}

struct MultiContactCandidate {
    int sample_id = -1;
    int patch_id = -1;
    double gap = 0.0;
    double area = 0.0;
    double weight = 1.0;
    ChVector3d normal = ChVector3d(0, 1, 0);
    ChVector3d world_pos = ChVector3d(0, 0, 0);
    double jacobian = 0.0;
    double jacobian_theta22_derivative = 0.0;
};

struct MultiStepDiagnostics {
    bool success = false;
    int iterations = 0;
    double residual_norm = 0.0;
    double min_gap = 0.0;
    double max_penetration = 0.0;
    double max_lambda_n = 0.0;
    double ncp_residual_norm = 0.0;
    double max_complementarity_error = 0.0;
    int patch_count = 0;
    std::vector<MultiContactCandidate> candidates;
    std::vector<double> lambdas;
    std::vector<double> lambda_forces;
    std::vector<double> ncp_residuals;
    std::vector<double> complementarity_errors;
};

void ResetMultiStepDiagnostics(MultiStepDiagnostics& diag) {
    diag.success = false;
    diag.iterations = 0;
    diag.residual_norm = 0.0;
    diag.min_gap = std::numeric_limits<double>::max();
    diag.max_penetration = 0.0;
    diag.max_lambda_n = 0.0;
    diag.ncp_residual_norm = 0.0;
    diag.max_complementarity_error = 0.0;
    diag.patch_count = 0;
    diag.candidates.clear();
    diag.lambdas.clear();
    diag.lambda_forces.clear();
    diag.ncp_residuals.clear();
    diag.complementarity_errors.clear();
}

void PopulateMultiStepDiagnostics(const std::vector<ChSdfNcpGeneralizedContact>& contacts,
                                  const ChSdfNcpGeneralizedDiagnostics& generalized,
                                  MultiStepDiagnostics& diag) {
    ResetMultiStepDiagnostics(diag);
    diag.success = generalized.success;
    diag.iterations = generalized.iterations;
    diag.residual_norm = generalized.residual_norm;
    diag.max_penetration = generalized.max_penetration;
    diag.max_lambda_n = generalized.max_lambda_n;
    diag.ncp_residual_norm = generalized.ncp_residual_norm;
    diag.max_complementarity_error = generalized.max_complementarity_error;
    diag.min_gap = contacts.empty() ? 0.0 : std::numeric_limits<double>::max();

    for (const auto& contact : contacts) {
        MultiContactCandidate candidate;
        candidate.sample_id = contact.contact_id;
        candidate.patch_id = contact.patch_id;
        candidate.gap = contact.gap;
        candidate.area = contact.weight;
        candidate.weight = contact.weight;
        candidate.jacobian = contact.jacobian.empty() ? 0.0 : contact.jacobian.front();
        diag.candidates.push_back(candidate);
        diag.lambdas.push_back(contact.lambda_n);
        diag.lambda_forces.push_back(contact.lambda_n * contact.weight);
        diag.ncp_residuals.push_back(std::abs(contact.ncp_residual));
        diag.complementarity_errors.push_back(contact.complementarity_error);
        diag.min_gap = std::min(diag.min_gap, contact.gap);
    }
}

std::vector<int> BuildPatchIdLookup(const ChSdfContactSurfaceGraph& surface_graph,
                                    const std::vector<int>& sample_ids,
                                    int* patch_count = nullptr) {
    std::vector<int> patch_ids(surface_graph.samples.size(), -1);
    const auto components = surface_graph.FindConnectedComponents(sample_ids);
    if (patch_count) {
        *patch_count = static_cast<int>(components.size());
    }
    for (size_t patch = 0; patch < components.size(); patch++) {
        for (int sample_id : components[patch]) {
            if (sample_id >= 0 && sample_id < static_cast<int>(patch_ids.size())) {
                patch_ids[static_cast<size_t>(sample_id)] = static_cast<int>(patch);
            }
        }
    }
    return patch_ids;
}

std::vector<double> BuildPatchWeightLookup(const ChSdfContactSurfaceGraph& surface_graph,
                                           const std::vector<int>& sample_ids) {
    std::vector<double> weights(surface_graph.samples.size(), 0.0);
    double total_area = 0.0;
    for (int sample_id : sample_ids) {
        if (sample_id >= 0 && sample_id < static_cast<int>(surface_graph.samples.size())) {
            total_area += std::max(0.0, surface_graph.samples[static_cast<size_t>(sample_id)].area);
        }
    }
    if (total_area > 1.0e-30) {
        for (int sample_id : sample_ids) {
            if (sample_id >= 0 && sample_id < static_cast<int>(surface_graph.samples.size())) {
                weights[static_cast<size_t>(sample_id)] =
                    std::max(0.0, surface_graph.samples[static_cast<size_t>(sample_id)].area) / total_area;
            }
        }
    } else if (!sample_ids.empty()) {
        const double uniform = 1.0 / static_cast<double>(sample_ids.size());
        for (int sample_id : sample_ids) {
            if (sample_id >= 0 && sample_id < static_cast<int>(surface_graph.samples.size())) {
                weights[static_cast<size_t>(sample_id)] = uniform;
            }
        }
    }
    return weights;
}

void Accumulate(BenchmarkStats& stats, const MultiStepDiagnostics& diag) {
    stats.max_penetration = std::max(stats.max_penetration, diag.max_penetration);
    stats.sum_penetration += diag.max_penetration;
    stats.max_lambda_n = std::max(stats.max_lambda_n, diag.max_lambda_n);
    stats.max_ncp_residual = std::max(stats.max_ncp_residual, diag.ncp_residual_norm);
    stats.sum_ncp_residual += diag.ncp_residual_norm;
    stats.max_complementarity_error = std::max(stats.max_complementarity_error, diag.max_complementarity_error);
    stats.sum_complementarity_error += diag.max_complementarity_error;
    stats.sum_iterations += static_cast<double>(diag.iterations);
    stats.sum_success += diag.success ? 1.0 : 0.0;
    stats.num_contacts = std::max(stats.num_contacts, static_cast<int>(diag.candidates.size()));
    stats.samples++;
}

double CylinderMass(double radius, double thickness, double density) {
    return kPi * radius * radius * thickness * density;
}

struct CamLikeConfig {
    std::string case_name;
    std::string asset;
    std::filesystem::path asset_dir;
    std::filesystem::path cam_obj;
    std::filesystem::path follower_obj;
    double dt = 0.002;
    double t_end = 0.25;
    double eps = 1.0e-7;
    double tolerance = 1.0e-8;
    int max_iterations = 40;
    double voxel_size = 1.0e-4;
    float half_width_voxels = 24.0f;
    double activation_band = 1.0e-6;
    double max_activation_band = 2.0e-5;
    int min_patch_samples = 4;
    double mass_follower = 1.0;
    double gravity_y = -9.81;
    double cam_motor_speed = -2.0;
    double phase = 0.0;
    double follower_y0 = 0.0;
    double cam_radius = 0.03;
    double cam_eccentricity = 0.006;
    double roller_radius = 0.01;
    std::string notes;
};

CamLikeConfig MakeCamLikeConfig(const std::string& case_name) {
    const auto root = GetProjectRoot();
    CamLikeConfig config;
    config.case_name = case_name;
    config.asset = "assets/" + case_name;
    config.asset_dir = root / "assets" / case_name;
    std::string model_json;

    if (case_name == "eccentric_roller") {
        model_json = ReadTextFile(config.asset_dir / "eccentric_roller_model.json");
        config.cam_obj = config.asset_dir / "models" / "eccentric_disk_cam.obj";
        config.follower_obj = config.asset_dir / "models" / "roller_follower.obj";
        config.notes =
            "full eccentric_roller OBJ/OpenVDB SDF; adaptive active-band connected patch SDF-NCP samples; "
            "generic active-patch area-normalized SDF-NCP assembly; AABB BVH broad phase filters follower samples";
    } else if (case_name == "onset_stress") {
        model_json = ReadTextFile(config.asset_dir / "onset_stress_model.json");
        config.cam_obj = config.asset_dir / "models" / "onset_cam.obj";
        config.follower_obj = config.asset_dir / "models" / "roller_follower.obj";
        config.notes =
            "full onset_stress OBJ/OpenVDB SDF; adaptive active-band connected patch SDF-NCP samples; "
            "generic active-patch area-normalized SDF-NCP assembly; AABB BVH broad phase filters follower samples";
    } else {
        throw std::invalid_argument("Unknown cam-like SDF-NCP benchmark: " + case_name);
    }

    config.cam_radius = ExtractJsonDouble(model_json, "cam_radius", config.cam_radius);
    config.cam_eccentricity = ExtractJsonDouble(model_json, "cam_eccentricity", config.cam_eccentricity);
    config.roller_radius = ExtractJsonDouble(model_json, "roller_radius", config.roller_radius);
    const double roller_thickness = ExtractJsonDouble(model_json, "roller_thickness", 0.018);
    const double density = ExtractJsonDouble(model_json, "density", 7800.0);
    config.mass_follower = CylinderMass(config.roller_radius, roller_thickness, density);
    config.gravity_y = ExtractJsonDouble(model_json, "gravity_y", config.gravity_y);
    config.cam_motor_speed = ExtractJsonDouble(model_json, "motor_speed", config.cam_motor_speed);
    config.dt = ExtractJsonDouble(model_json, "step_size", config.dt);
    config.t_end = ExtractJsonDouble(model_json, "total_time", config.t_end);
    config.phase = ExtractJsonDouble(model_json, "phase", 0.0);
    const ChVector3d follower_pos = ExtractJsonVector3(model_json, "follower_init_pos", ChVector3d(0, 0.04, 0));
    config.follower_y0 = follower_pos.y();
    config.activation_band = 1.0e-6;
    config.max_activation_band = 1.6e-5;
    config.min_patch_samples = case_name == "eccentric_roller" ? 6 : 4;
    return config;
}

ChSdfContactSampleQuery QueryFollowerSampleAgainstCamLike(const ChOpenVdbSdfGrid& cam_sdf,
                                                          const ChSdfContactSurfaceGraph& follower_graph,
                                                          int sample_id,
                                                          double follower_y,
                                                          double follower_vy,
                                                          double theta_cam,
                                                          double omega_cam) {
    const auto& sample = follower_graph.samples.at(static_cast<size_t>(sample_id));
    const ChVector3d world_pos = sample.local_pos + ChVector3d(0, follower_y, 0);
    const ChVector3d world_vel = ChVector3d(0, follower_vy, 0);
    const ChVector3d cam_frame_pos = RotateZ(world_pos, -theta_cam);
    const ChVector3d cam_frame_vel = RotateZ(world_vel - OmegaZCross(omega_cam, world_pos), -theta_cam);
    ChSdfContactSampleQuery query = cam_sdf.QueryLocal(cam_frame_pos, cam_frame_vel);
    query.world_pos = world_pos;
    query.world_vel = world_vel;
    query.grad = RotateZ(query.grad, theta_cam);
    return query;
}

double SampleFollowerPhiAgainstCamLike(const ChOpenVdbSdfGrid& cam_sdf,
                                       const ChSdfContactSurfaceGraph& follower_graph,
                                       int sample_id,
                                       double follower_y,
                                       double theta_cam) {
    const auto& sample = follower_graph.samples.at(static_cast<size_t>(sample_id));
    const ChVector3d world_pos = sample.local_pos + ChVector3d(0, follower_y, 0);
    const ChVector3d cam_frame_pos = RotateZ(world_pos, -theta_cam);
    return cam_sdf.SamplePhi(cam_frame_pos);
}

std::vector<int> BuildCamLikePatchSamples(const ChOpenVdbSdfGrid& cam_sdf,
                                          const ChSdfContactSurfaceGraph& follower_graph,
                                          double follower_y,
                                          double follower_vy,
                                          double theta_cam,
                                          double omega_cam,
                                          double activation_band,
                                          double max_activation_band,
                                          int min_patch_samples,
                                          int* patch_count = nullptr,
                                          const ChSdfContactSampleBvh* follower_bvh = nullptr,
                                          const ChSdfAabb* cam_local_bounds = nullptr,
                                          double broad_phase_padding = 0.0) {
    (void)follower_vy;
    (void)omega_cam;
    auto all_sample_ids = [&]() {
        std::vector<int> ids;
        ids.reserve(follower_graph.samples.size());
        for (const auto& sample : follower_graph.samples) {
            if (sample.area > 0.0) {
                ids.push_back(sample.id);
            }
        }
        return ids;
    };

    std::vector<int> candidate_ids;
    if (follower_bvh && cam_local_bounds && !follower_bvh->Empty() && cam_local_bounds->IsValid()) {
        const ChSdfAabb query_bounds =
            CamLocalAabbInFollowerLocal(*cam_local_bounds, follower_y, theta_cam, broad_phase_padding);
        candidate_ids = follower_bvh->Query(query_bounds);
        std::sort(candidate_ids.begin(), candidate_ids.end());
        candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());
    } else {
        candidate_ids = all_sample_ids();
    }

    if (candidate_ids.empty()) {
        if (patch_count) {
            *patch_count = 0;
        }
        return {};
    }

    std::vector<std::pair<int, double>> sampled;
    sampled.reserve(candidate_ids.size());
    auto collect_sampled = [&](const std::vector<int>& ids, double& min_phi) {
        sampled.clear();
        min_phi = std::numeric_limits<double>::max();
        for (int sample_id : ids) {
            if (sample_id < 0 || sample_id >= static_cast<int>(follower_graph.samples.size())) {
                continue;
            }
            const auto& sample = follower_graph.samples[static_cast<size_t>(sample_id)];
            if (sample.area <= 0.0) {
                continue;
            }
            const double phi = SampleFollowerPhiAgainstCamLike(cam_sdf, follower_graph, sample.id, follower_y, theta_cam);
            sampled.emplace_back(sample.id, phi);
            min_phi = std::min(min_phi, phi);
        }
    };

    double min_phi = std::numeric_limits<double>::max();
    collect_sampled(candidate_ids, min_phi);
    if (min_phi > std::max(activation_band, max_activation_band)) {
        if (patch_count) {
            *patch_count = 0;
        }
        return {};
    }

    std::vector<int> active;
    auto build_active = [&]() {
        active.clear();
        double band = std::max(activation_band, 0.0);
        const double band_limit = std::max(band, max_activation_band);
        for (;;) {
            active.clear();
            const double threshold = min_phi + band;
            for (const auto& item : sampled) {
                if (item.second <= threshold) {
                    active.push_back(item.first);
                }
            }
            if (static_cast<int>(active.size()) >= std::max(1, min_patch_samples) || band >= band_limit) {
                break;
            }
            band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
        }
    };
    active.reserve(sampled.size());
    build_active();

    if (follower_bvh && active.size() < static_cast<size_t>(std::max(1, min_patch_samples)) &&
        candidate_ids.size() < follower_graph.samples.size()) {
        candidate_ids = all_sample_ids();
        sampled.reserve(candidate_ids.size());
        collect_sampled(candidate_ids, min_phi);
        build_active();
    }

    const auto components = follower_graph.FindConnectedComponents(active);
    if (patch_count) {
        *patch_count = static_cast<int>(components.size());
    }

    std::vector<int> ids;
    for (const auto& component : components) {
        ids.insert(ids.end(), component.begin(), component.end());
    }
    return ids;
}

double MinFollowerGapAgainstCamLike(const ChOpenVdbSdfGrid& cam_sdf,
                                    const ChSdfContactSurfaceGraph& follower_graph,
                                    double follower_y,
                                    double follower_vy,
                                    double theta_cam,
                                    double omega_cam) {
    (void)follower_vy;
    (void)omega_cam;
    double min_gap = std::numeric_limits<double>::max();
    for (const auto& sample : follower_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        min_gap = std::min(min_gap,
                           SampleFollowerPhiAgainstCamLike(cam_sdf, follower_graph, sample.id, follower_y, theta_cam));
    }
    return min_gap;
}

struct CamLikeState {
    double y = 0.0;
    double vy = 0.0;
    double theta = 0.0;
};

ChSdfNcpGeneralizedProblem MakeCamLikeGeneralizedProblem(const ChOpenVdbSdfGrid& cam_sdf,
                                                         const ChSdfContactSurfaceGraph& follower_graph,
                                                         const CamLikeConfig& config,
                                                         const CamLikeState& state,
                                                         const std::vector<int>& sample_ids,
                                                         double theta_next,
                                                         int* patch_count = nullptr) {
    int local_patch_count = 0;
    const std::vector<int> patch_ids = BuildPatchIdLookup(follower_graph, sample_ids, &local_patch_count);
    const std::vector<double> patch_weights = BuildPatchWeightLookup(follower_graph, sample_ids);
    if (patch_count) {
        *patch_count = local_patch_count;
    }

    ChSdfNcpGeneralizedProblem problem;
    problem.dt = config.dt;
    problem.eps = config.eps;
    problem.tolerance = config.tolerance;
    problem.relaxed_tolerance = 1.0e-4;
    problem.negative_lambda_tolerance = 1.0e-4;
    problem.max_iterations = config.max_iterations;
    problem.current_velocity = {state.vy};
    problem.mass_diagonal = {config.mass_follower};
    problem.external_force = {config.mass_follower * config.gravity_y};
    problem.contact_count = sample_ids.size();
    problem.evaluate_contacts = [&, sample_ids, patch_ids, patch_weights, theta_next](const std::vector<double>& v_next) {
        const double vy_next = v_next.at(0);
        const double y_next = state.y + config.dt * vy_next;
        std::vector<ChSdfNcpGeneralizedContact> contacts;
        contacts.reserve(sample_ids.size());
        for (int sample_id : sample_ids) {
            const auto query = QueryFollowerSampleAgainstCamLike(
                cam_sdf, follower_graph, sample_id, y_next, vy_next, theta_next, config.cam_motor_speed);
            const ChVector3d normal =
                NormalizeSdfVector(query.grad, "OpenVDB cam-like SDF gradient must be nonzero.");

            ChSdfNcpGeneralizedContact contact;
            contact.gap = query.phi;
            contact.jacobian = {normal.y()};
            contact.weight = sample_id >= 0 && sample_id < static_cast<int>(patch_weights.size()) ?
                                 patch_weights[static_cast<size_t>(sample_id)] :
                                 0.0;
            contact.contact_id = sample_id;
            if (sample_id >= 0 && sample_id < static_cast<int>(patch_ids.size())) {
                contact.patch_id = patch_ids[static_cast<size_t>(sample_id)];
            }
            contacts.push_back(contact);
        }
        return contacts;
    };
    return problem;
}

std::vector<double> ComputeCamLikeResidual(const ChOpenVdbSdfGrid& cam_sdf,
                                           const ChSdfContactSurfaceGraph& follower_graph,
                                           const CamLikeConfig& config,
                                           const CamLikeState& state,
                                           const std::vector<int>& sample_ids,
                                           double theta_next,
                                           const std::vector<double>& z,
                                           MultiStepDiagnostics* diag = nullptr) {
    const size_t n = sample_ids.size();
    if (z.size() != n + 1) {
        throw std::invalid_argument("Cam-like multi-contact residual has inconsistent dimension.");
    }

    int patch_count = 0;
    const auto problem =
        MakeCamLikeGeneralizedProblem(cam_sdf, follower_graph, config, state, sample_ids, theta_next, &patch_count);

    const auto residual = ComputeSdfNcpGeneralizedResidual(problem, z);
    if (diag) {
        const auto generalized_diag =
            MakeSdfNcpGeneralizedDiagnostics(problem, z, true, 0, "evaluated through generalized SDF-NCP backend");
        PopulateMultiStepDiagnostics(residual.contacts, generalized_diag, *diag);
        diag->patch_count = patch_count;
    }
    return residual.value;
}

CamLikeState StepCamLikeFallback(const CamLikeConfig& config, const CamLikeState& state) {
    CamLikeState next = state;
    next.vy = state.vy + config.dt * config.gravity_y;
    next.y = state.y + config.dt * next.vy;
    next.theta = state.theta + config.dt * config.cam_motor_speed;
    return next;
}

MultiStepDiagnostics EvaluateCamLikeState(const ChOpenVdbSdfGrid& cam_sdf,
                                          const ChSdfContactSurfaceGraph& follower_graph,
                                          const CamLikeConfig& config,
                                          const CamLikeState& state,
                                          const std::vector<int>& sample_ids) {
    MultiStepDiagnostics diag;
    diag.success = true;
    std::vector<double> z(1 + sample_ids.size(), 0.0);
    z[0] = state.vy;
    ComputeCamLikeResidual(cam_sdf, follower_graph, config, state, sample_ids, state.theta, z, &diag);
    diag.patch_count = static_cast<int>(follower_graph.FindConnectedComponents(sample_ids).size());
    diag.min_gap = MinFollowerGapAgainstCamLike(
        cam_sdf, follower_graph, state.y, state.vy, state.theta, config.cam_motor_speed);
    diag.max_penetration = std::max(diag.max_penetration, std::max(0.0, -diag.min_gap));
    return diag;
}

struct CamLikeStepResult {
    CamLikeState state;
    MultiStepDiagnostics diagnostics;
};

CamLikeStepResult SolveCamLikeStep(const ChOpenVdbSdfGrid& cam_sdf,
                                   const ChSdfContactSurfaceGraph& follower_graph,
                                   const CamLikeConfig& config,
                                   const CamLikeState& state,
                                   const ChSdfContactSampleBvh* follower_bvh = nullptr,
                                   const ChSdfAabb* cam_local_bounds = nullptr,
                                   double broad_phase_padding = 0.0) {
    CamLikeStepResult result;
    result.state = StepCamLikeFallback(config, state);
    const double theta_next = state.theta + config.dt * config.cam_motor_speed;
    const double vy_free = state.vy + config.dt * config.gravity_y;
    const double y_free = state.y + config.dt * vy_free;
    int patch_count = 0;
    const std::vector<int> sample_ids = BuildCamLikePatchSamples(cam_sdf,
                                                                 follower_graph,
                                                                 y_free,
                                                                 vy_free,
                                                                  theta_next,
                                                                  config.cam_motor_speed,
                                                                  config.activation_band,
                                                                  config.max_activation_band,
                                                                  config.min_patch_samples,
                                                                  &patch_count,
                                                                  follower_bvh,
                                                                  cam_local_bounds,
                                                                  broad_phase_padding);
    if (sample_ids.empty()) {
        result.diagnostics.success = true;
        result.diagnostics.min_gap = MinFollowerGapAgainstCamLike(
            cam_sdf, follower_graph, result.state.y, result.state.vy, result.state.theta, config.cam_motor_speed);
        result.diagnostics.max_penetration = std::max(0.0, -result.diagnostics.min_gap);
        return result;
    }

    std::vector<double> z(1 + sample_ids.size(), 0.0);
    z[0] = vy_free;
    MultiStepDiagnostics guess_diag;
    ComputeCamLikeResidual(cam_sdf, follower_graph, config, state, sample_ids, theta_next, z, &guess_diag);
    double weighted_penetration = 0.0;
    for (const auto& candidate : guess_diag.candidates) {
        weighted_penetration += std::max(0.0, -candidate.gap) * candidate.weight;
    }
    if (weighted_penetration > 0.0 && sample_ids.size() <= 2) {
        const double unconstrained_lambda_total = weighted_penetration * config.mass_follower /
                                                  std::max(config.dt * config.dt, 1.0e-12);
        const double force_scale = std::max(1.0, std::abs(config.mass_follower * config.gravity_y));
        const double lambda_total = std::min(unconstrained_lambda_total, 100.0 * force_scale);
        for (size_t i = 0; i < guess_diag.candidates.size(); i++) {
            z[i + 1] = lambda_total * std::max(0.0, -guess_diag.candidates[i].gap) / weighted_penetration;
        }
    }

    const auto problem = MakeCamLikeGeneralizedProblem(
        cam_sdf, follower_graph, config, state, sample_ids, theta_next, &patch_count);
    const auto solved = SolveSdfNcpGeneralizedProblem(problem, z);

    result.state.vy = solved.next_velocity.at(0);
    result.state.y = state.y + config.dt * result.state.vy;
    result.state.theta = theta_next;
    PopulateMultiStepDiagnostics(solved.contacts, solved.diagnostics, result.diagnostics);
    result.diagnostics.patch_count = patch_count;
    result.diagnostics.min_gap = MinFollowerGapAgainstCamLike(cam_sdf,
                                                              follower_graph,
                                                              result.state.y,
                                                              result.state.vy,
                                                              result.state.theta,
                                                              config.cam_motor_speed);
    result.diagnostics.max_penetration =
        std::max(result.diagnostics.max_penetration, std::max(0.0, -result.diagnostics.min_gap));
    return result;
}

double CamEnvelopeY(double eccentricity, double cam_radius, double roller_radius, double theta) {
    const double cx = eccentricity * std::cos(theta);
    const double cy = eccentricity * std::sin(theta);
    const double reach = cam_radius + roller_radius;
    return cy + std::sqrt(std::max(0.0, reach * reach - cx * cx));
}

int RunCamLikeFullGeometryCase(const std::string& case_name) {
    const auto start = std::chrono::steady_clock::now();
    const CamLikeConfig config = MakeCamLikeConfig(case_name);
    const auto root = GetProjectRoot();
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / case_name;
    std::filesystem::create_directories(out_dir);

    openvdb::initialize();
    ChSdfTriangleMeshData cam_mesh = LoadWavefrontMeshForSdf(config.cam_obj);
    ChSdfTriangleMeshData follower_mesh = LoadWavefrontMeshForSdf(config.follower_obj);
    ChOpenVdbSdfGrid cam_sdf = BuildOpenVdbLevelSetFromMesh(
        cam_mesh, config.voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph follower_graph = MakeSurfaceGraphFromMeshForSdf(follower_mesh);
    ChSdfContactSampleBvh follower_bvh(follower_graph);
    const double broad_phase_padding =
        std::max(config.activation_band, config.max_activation_band) + 8.0 * std::max(config.voxel_size, 1.0e-12);

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    std::ofstream analytic_comparison(out_dir / "analytic_comparison.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    analytic_comparison << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";
    analytic_comparison
        << "time,analytic_follower_y,sdf_ncp_follower_y,error,abs_error,min_gap,max_penetration,success\n";

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = config.case_name;
    stats.asset = config.asset;
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    stats.num_contacts = 0;
    stats.notes = config.notes;

    CamLikeState state;
    state.y = config.follower_y0;
    state.vy = 0.0;
    state.theta = 0.0;
    double observed_onset = std::numeric_limits<double>::quiet_NaN();
    double previous_gap = MinFollowerGapAgainstCamLike(
        cam_sdf, follower_graph, state.y, state.vy, state.theta, config.cam_motor_speed);
    double previous_time = 0.0;
    double analytic_error_sq_sum = 0.0;
    double analytic_max_abs_error = 0.0;
    int analytic_sample_count = 0;
    auto record_analytic_comparison = [&](double time, const CamLikeState& current, const MultiStepDiagnostics& diag) {
        const double reference = CamEnvelopeY(config.cam_eccentricity,
                                             config.cam_radius,
                                             config.roller_radius,
                                             config.phase + config.cam_motor_speed * time);
        const double error = current.y - reference;
        const double abs_error = std::abs(error);
        analytic_error_sq_sum += error * error;
        analytic_max_abs_error = std::max(analytic_max_abs_error, abs_error);
        analytic_sample_count++;
        analytic_comparison << time << "," << reference << "," << current.y << "," << error << ","
                            << abs_error << "," << diag.min_gap << "," << diag.max_penetration << ","
                            << (diag.success ? 1 : 0) << "\n";
    };

    for (int i = 0; i <= steps; i++) {
        if (i > 0) {
            const CamLikeStepResult step =
                SolveCamLikeStep(cam_sdf, follower_graph, config, state, &follower_bvh, &cam_sdf.local_bounds, broad_phase_padding);
            state = step.state;
            const double time = static_cast<double>(i) * config.dt;
            if (case_name == "onset_stress" && std::isnan(observed_onset) && previous_gap > 0.0 &&
                step.diagnostics.min_gap <= 0.0) {
                const double alpha = previous_gap / (previous_gap - step.diagnostics.min_gap);
                observed_onset = previous_time + alpha * (time - previous_time);
            }
            previous_gap = step.diagnostics.min_gap;
            previous_time = time;

            if (step.diagnostics.candidates.empty()) {
                trajectory << time << ",follower,0," << state.y << ",0,1,0,0,0,0," << state.vy
                            << ",0,0,0,0,-1," << step.diagnostics.min_gap << ",0,"
                            << "0,0," << std::max(0.0, -step.diagnostics.min_gap) << ",0,0,"
                            << step.diagnostics.residual_norm << "," << step.diagnostics.iterations << ","
                            << (step.diagnostics.success ? 1 : 0) << "\n";
            }
            for (size_t c = 0; c < step.diagnostics.candidates.size(); c++) {
                const auto& candidate = step.diagnostics.candidates[c];
                const double lambda = c < step.diagnostics.lambdas.size() ? step.diagnostics.lambdas[c] : 0.0;
                const double lambda_force =
                    c < step.diagnostics.lambda_forces.size() ? step.diagnostics.lambda_forces[c] :
                                                                 lambda * candidate.weight;
                const double ncp = c < step.diagnostics.ncp_residuals.size() ? step.diagnostics.ncp_residuals[c] : 0.0;
                const double comp = c < step.diagnostics.complementarity_errors.size() ?
                                        step.diagnostics.complementarity_errors[c] :
                                        0.0;
                trajectory << time << ",follower,0," << state.y << ",0,1,0,0,0,0," << state.vy
                            << ",0,0,0,0," << candidate.sample_id << "," << candidate.gap << "," << lambda
                            << "," << candidate.weight << "," << lambda_force << ","
                            << std::max(0.0, -candidate.gap) << "," << ncp << "," << comp << ","
                            << step.diagnostics.residual_norm << "," << step.diagnostics.iterations << ","
                            << (step.diagnostics.success ? 1 : 0) << "\n";
            }
            Accumulate(stats, step.diagnostics);
            record_analytic_comparison(time, state, step.diagnostics);
        } else {
            int patch_count = 0;
            const std::vector<int> ids = BuildCamLikePatchSamples(cam_sdf,
                                                                  follower_graph,
                                                                  state.y,
                                                                  state.vy,
                                                                   state.theta,
                                                                   config.cam_motor_speed,
                                                                   config.activation_band,
                                                                   config.max_activation_band,
                                                                   config.min_patch_samples,
                                                                   &patch_count,
                                                                   &follower_bvh,
                                                                   &cam_sdf.local_bounds,
                                                                   broad_phase_padding);
            const MultiStepDiagnostics diag = EvaluateCamLikeState(cam_sdf, follower_graph, config, state, ids);
            if (diag.candidates.empty()) {
                trajectory << "0,follower,0," << state.y << ",0,1,0,0,0,0," << state.vy
                            << ",0,0,0,0,-1," << diag.min_gap << ",0,"
                            << "0,0," << std::max(0.0, -diag.min_gap) << ",0,0,0,0,1\n";
            }
            for (const auto& candidate : diag.candidates) {
                trajectory << "0,follower,0," << state.y << ",0,1,0,0,0,0," << state.vy
                            << ",0,0,0,0," << candidate.sample_id << "," << candidate.gap
                            << ",0," << candidate.weight << ",0," << std::max(0.0, -candidate.gap) << ","
                            << std::abs(SmoothFischerBurmeister(candidate.gap, 0.0, config.eps)) << ","
                            << ComplementarityError(candidate.gap, 0.0) << ",0,0,1\n";
            }
            Accumulate(stats, diag);
            record_analytic_comparison(0.0, state, diag);
        }
    }

    const double final_reference = CamEnvelopeY(config.cam_eccentricity,
                                               config.cam_radius,
                                               config.roller_radius,
                                               config.phase + config.cam_motor_speed * config.t_end);
    const double abs_y = std::abs(state.y - final_reference);
    const double rel_y = std::abs(final_reference) > 1.0e-14 ? abs_y / std::abs(final_reference) : 0.0;
    comparison << case_name << ",follower_y_analytic_envelope," << final_reference << "," << state.y << ","
               << abs_y << "," << rel_y
               << ",analytic envelope is a geometry reference only; SDF-NCP uses full OBJ/OpenVDB contact samples\n";
    const double analytic_rmse =
        std::sqrt(analytic_error_sq_sum / static_cast<double>(std::max(1, analytic_sample_count)));
    comparison << case_name << ",follower_y_analytic_envelope_rmse,0," << analytic_rmse << ","
               << analytic_rmse << ",0,RMSE over analytic_comparison.csv time samples\n";
    comparison << case_name << ",follower_y_analytic_envelope_max_abs_error,0," << analytic_max_abs_error << ","
               << analytic_max_abs_error << ",0,max absolute error over analytic_comparison.csv time samples\n";
    if (case_name == "onset_stress") {
        const double target = ExtractJsonDouble(ReadTextFile(config.asset_dir / "onset_stress_model.json"),
                                                "target_onset_time",
                                                0.15);
        const double onset = std::isnan(observed_onset) ? -1.0 : observed_onset;
        const double abs_onset = std::isnan(observed_onset) ? -1.0 : std::abs(observed_onset - target);
        comparison << case_name << ",onset_time," << target << "," << onset << "," << abs_onset
                   << ",0,first min-gap sign change from full OBJ/OpenVDB SDF-NCP rollout\n";
    }

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);
    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    const double success_rate = stats.sum_success / static_cast<double>(std::max(1, stats.samples));
    return stats.max_penetration < 1.0e-2 && stats.max_complementarity_error < 1.0e-2 && success_rate > 0.90 ? 0 : 1;
}

constexpr double kMmToM = 1.0e-3;

struct RmdBodyInfo {
    ChVector3d cm_marker_mm = ChVector3d(0, 0, 0);
    Mat3 part_rotation;
    ChVector3d surface_ref_marker_mm = ChVector3d(0, 0, 0);
    Mat3 surface_ref_rotation;
    double mass_kg = 0.0;
    double inertia_x_kg_m2 = 0.0;
    bool has_part = false;
    bool has_cm = false;
    bool has_surface_ref = false;
    bool has_inertia = false;
};

struct SimpleGearJointInfo {
    int id = -1;
    std::string name;
    int i_marker_id = -1;
    int j_marker_id = -1;
    std::string type;
    bool has_joint = false;
};

struct SimpleGearMotionInfo {
    int id = -1;
    std::string name;
    int joint_id = -1;
    std::string coordinate_type;
    std::string input_type;
    double function_constant = 0.0;
    bool has_motion = false;
    bool has_function = false;
};

struct SimpleGearGGeomContactInfo {
    int id = -1;
    std::string name;
    std::string contact_type;
    int i_float_marker_id = -1;
    int j_float_marker_id = -1;
    int i_geometry_id = -1;
    int j_geometry_id = -1;
    std::string i_direction;
    std::string j_direction;
    int i_cp_flag = 0;
    int j_cp_flag = 0;
    int i_node_contact = 0;
    int j_node_contact = 0;
    int edge_contact = 0;
    int smooth_node_contact = 0;
    int smooth_face_contact = 0;
    int use_cpm = 0;
    int max_contact_points = 0;
    int contact_force_type = 0;
    double bpen = 0.0;
    double rdf = 0.0;
    double max_pen = 0.0;
    double hparameter = 0.0;
    double stiffness = 0.0;
    double korder = 0.0;
    double damping = 0.0;
    double dynamic_friction_coefficient = 0.0;
    double static_transition_velocity = 0.0;
    double dynamic_transition_velocity = 0.0;
    double static_friction_coefficient = 0.0;
    bool has_contact = false;
};

struct SimpleGearRmd {
    RmdBodyInfo gear21;
    RmdBodyInfo gear22;
    SimpleGearJointInfo rev_joint1;
    SimpleGearJointInfo rev_joint2;
    SimpleGearMotionInfo rev_joint1_motion;
    SimpleGearGGeomContactInfo geo_surface_contact;
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

void RequireBodyInfo(const RmdBodyInfo& body, const std::string& name) {
    if (!body.has_part || !body.has_cm || !body.has_surface_ref || !body.has_inertia) {
        throw std::runtime_error("Incomplete RMD body data for " + name);
    }
}

SimpleGearRmd LoadSimpleGearRmdForNcp(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open RMD: " + path.string());
    }

    enum class Block { None, Part, Marker, Joint, Motion, GGeomContact };
    Block block = Block::None;
    SimpleGearRmd rmd;
    std::string current_part_name;
    std::string current_marker_name;
    std::string current_joint_name;
    std::string current_motion_name;
    std::string current_contact_name;
    int current_id = -1;
    std::string line;
    while (std::getline(in, line)) {
        const std::string trimmed = Trim(line);
        if (trimmed.empty()) {
            continue;
        }

        auto start_block = [&](Block next, int id) {
            block = next;
            current_id = id;
            current_part_name.clear();
            current_marker_name.clear();
            current_joint_name.clear();
            current_motion_name.clear();
            current_contact_name.clear();
        };

        if (trimmed.rfind("PART /", 0) == 0) {
            start_block(Block::Part, ParseRmdEntityId(trimmed));
            continue;
        }
        if (trimmed.rfind("MARKER /", 0) == 0) {
            start_block(Block::Marker, ParseRmdEntityId(trimmed));
            continue;
        }
        if (trimmed.rfind("JOINT /", 0) == 0) {
            start_block(Block::Joint, ParseRmdEntityId(trimmed));
            continue;
        }
        if (trimmed.rfind("MOTION /", 0) == 0) {
            start_block(Block::Motion, ParseRmdEntityId(trimmed));
            continue;
        }
        if (trimmed.rfind("GGEOMCONTACT /", 0) == 0) {
            start_block(Block::GGeomContact, ParseRmdEntityId(trimmed));
            continue;
        }
        if (trimmed.rfind("EXPRESSION /", 0) == 0 || trimmed.rfind("OUTPUT/", 0) == 0 ||
            trimmed.rfind("UNITS /", 0) == 0 || trimmed.rfind("REQUEST /", 0) == 0) {
            start_block(Block::None, -1);
            continue;
        }

        if (Contains(trimmed, "NAME =")) {
            const std::string name = ExtractQuotedName(trimmed);
            if (block == Block::Part) {
                current_part_name = name;
            } else if (block == Block::Marker) {
                current_marker_name = name;
            } else if (block == Block::Joint) {
                current_joint_name = name;
                SimpleGearJointInfo* joint = nullptr;
                if (name == "RevJoint1") {
                    joint = &rmd.rev_joint1;
                } else if (name == "RevJoint2") {
                    joint = &rmd.rev_joint2;
                }
                if (joint) {
                    joint->id = current_id;
                    joint->name = name;
                    joint->has_joint = true;
                }
            } else if (block == Block::Motion) {
                current_motion_name = name;
                if (name == "RevJoint1.RMotion") {
                    rmd.rev_joint1_motion.id = current_id;
                    rmd.rev_joint1_motion.name = name;
                    rmd.rev_joint1_motion.has_motion = true;
                }
            } else if (block == Block::GGeomContact) {
                current_contact_name = name;
                if (name == "GeoSurContact1") {
                    rmd.geo_surface_contact.id = current_id;
                    rmd.geo_surface_contact.name = name;
                    rmd.geo_surface_contact.has_contact = true;
                }
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
        auto joint_for_name = [&]() -> SimpleGearJointInfo* {
            if (current_joint_name == "RevJoint1") {
                return &rmd.rev_joint1;
            }
            if (current_joint_name == "RevJoint2") {
                return &rmd.rev_joint2;
            }
            return nullptr;
        };
        auto motion_for_name = [&]() -> SimpleGearMotionInfo* {
            if (current_motion_name == "RevJoint1.RMotion") {
                return &rmd.rev_joint1_motion;
            }
            return nullptr;
        };
        auto contact_for_name = [&]() -> SimpleGearGGeomContactInfo* {
            if (current_contact_name == "GeoSurContact1") {
                return &rmd.geo_surface_contact;
            }
            return nullptr;
        };
        auto parse_int_field = [&](int& value) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                value = static_cast<int>(std::round(values.front()));
            }
        };
        auto parse_double_field = [&](double& value) {
            const auto values = ParseNumbersAfterEquals(trimmed);
            if (!values.empty()) {
                value = values.front();
            }
        };

        if (block == Block::Part) {
            if (Contains(trimmed, "MASS =")) {
                if (auto* body = body_for_part()) {
                    parse_double_field(body->mass_kg);
                }
            } else if (Contains(trimmed, "IP =")) {
                if (auto* body = body_for_part()) {
                    const auto values = ParseNumbersAfterEquals(trimmed);
                    if (!values.empty()) {
                        body->inertia_x_kg_m2 = values[0] * 1.0e-6;
                        body->has_inertia = true;
                    }
                }
            } else if (Contains(trimmed, "REULER =")) {
                if (auto* body = body_for_part()) {
                    const ChVector3d euler = ParseRmdTriple(trimmed);
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
                    const ChVector3d euler = ParseRmdTriple(trimmed);
                    ref_body->surface_ref_rotation = RecurDynEuler(euler.x(), euler.y(), euler.z());
                }
            }
        } else if (block == Block::Joint) {
            if (auto* joint = joint_for_name()) {
                if (HasRmdField(trimmed, "I")) {
                    parse_int_field(joint->i_marker_id);
                } else if (HasRmdField(trimmed, "J")) {
                    parse_int_field(joint->j_marker_id);
                } else if (Contains(trimmed, "Revolute")) {
                    joint->type = "Revolute";
                }
            }
        } else if (block == Block::Motion) {
            if (auto* motion = motion_for_name()) {
                if (HasRmdField(trimmed, "JOINT")) {
                    parse_int_field(motion->joint_id);
                } else if (Contains(trimmed, "ROTATION")) {
                    motion->coordinate_type = "ROTATION";
                } else if (Contains(trimmed, "TRANSLATION")) {
                    motion->coordinate_type = "TRANSLATION";
                } else if (Contains(trimmed, "VELOCITY")) {
                    motion->input_type = "VELOCITY";
                } else if (HasRmdField(trimmed, "FUNCTION")) {
                    parse_double_field(motion->function_constant);
                    motion->has_function = true;
                }
            }
        } else if (block == Block::GGeomContact) {
            if (auto* contact = contact_for_name()) {
                if (HasRmdField(trimmed, "CONTYPE")) {
                    contact->contact_type = ExtractRmdValueText(trimmed);
                } else if (HasRmdField(trimmed, "IFLOAT")) {
                    parse_int_field(contact->i_float_marker_id);
                } else if (HasRmdField(trimmed, "JFLOAT")) {
                    parse_int_field(contact->j_float_marker_id);
                } else if (HasRmdField(trimmed, "IGGEOMID")) {
                    parse_int_field(contact->i_geometry_id);
                } else if (HasRmdField(trimmed, "JGGEOMID")) {
                    parse_int_field(contact->j_geometry_id);
                } else if (HasRmdField(trimmed, "IDIRECTION")) {
                    contact->i_direction = ExtractRmdValueText(trimmed);
                } else if (HasRmdField(trimmed, "JDIRECTION")) {
                    contact->j_direction = ExtractRmdValueText(trimmed);
                } else if (HasRmdField(trimmed, "ICPFLAG")) {
                    parse_int_field(contact->i_cp_flag);
                } else if (HasRmdField(trimmed, "JCPFLAG")) {
                    parse_int_field(contact->j_cp_flag);
                } else if (HasRmdField(trimmed, "INODECONTACT")) {
                    parse_int_field(contact->i_node_contact);
                } else if (HasRmdField(trimmed, "JNODECONTACT")) {
                    parse_int_field(contact->j_node_contact);
                } else if (HasRmdField(trimmed, "EDGECONTACT")) {
                    parse_int_field(contact->edge_contact);
                } else if (HasRmdField(trimmed, "SMOOTHNODECONTACT")) {
                    parse_int_field(contact->smooth_node_contact);
                } else if (HasRmdField(trimmed, "SMOOTHFACECONTACT")) {
                    parse_int_field(contact->smooth_face_contact);
                } else if (HasRmdField(trimmed, "USECPM")) {
                    parse_int_field(contact->use_cpm);
                } else if (HasRmdField(trimmed, "NOMAXCP")) {
                    parse_int_field(contact->max_contact_points);
                } else if (HasRmdField(trimmed, "CFTYPE")) {
                    parse_int_field(contact->contact_force_type);
                } else if (HasRmdField(trimmed, "BPEN")) {
                    parse_double_field(contact->bpen);
                } else if (HasRmdField(trimmed, "RDF")) {
                    parse_double_field(contact->rdf);
                } else if (HasRmdField(trimmed, "MAXPEN")) {
                    parse_double_field(contact->max_pen);
                } else if (HasRmdField(trimmed, "HPARAMETER")) {
                    parse_double_field(contact->hparameter);
                } else if (HasRmdField(trimmed, "KORDER")) {
                    parse_double_field(contact->korder);
                } else if (HasRmdField(trimmed, "K")) {
                    parse_double_field(contact->stiffness);
                } else if (HasRmdField(trimmed, "C")) {
                    parse_double_field(contact->damping);
                } else if (HasRmdField(trimmed, "D_F_C")) {
                    parse_double_field(contact->dynamic_friction_coefficient);
                } else if (HasRmdField(trimmed, "S_T_V")) {
                    parse_double_field(contact->static_transition_velocity);
                } else if (HasRmdField(trimmed, "D_T_V")) {
                    parse_double_field(contact->dynamic_transition_velocity);
                } else if (HasRmdField(trimmed, "S_F_C")) {
                    parse_double_field(contact->static_friction_coefficient);
                }
            }
        }
    }

    RequireBodyInfo(rmd.gear21, "GEAR21");
    RequireBodyInfo(rmd.gear22, "GEAR22");
    return rmd;
}

void WriteSimpleGearRmdMappingAudit(const std::filesystem::path& path, const SimpleGearRmd& rmd) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "category,name,key,value,notes\n";

    auto write_body = [&](const std::string& name, const RmdBodyInfo& body) {
        out << "body," << name << ",mass_kg," << body.mass_kg << ",RMD PART MASS\n";
        out << "body," << name << ",inertia_x_kg_m2," << body.inertia_x_kg_m2
            << ",RMD IP first component converted from kg*mm^2 to kg*m^2\n";
        out << "body," << name << ",cm_marker_mm_x," << body.cm_marker_mm.x() << ",RMD CM marker QP\n";
        out << "body," << name << ",cm_marker_mm_y," << body.cm_marker_mm.y() << ",RMD CM marker QP\n";
        out << "body," << name << ",cm_marker_mm_z," << body.cm_marker_mm.z() << ",RMD CM marker QP\n";
        out << "body," << name << ",surface_ref_marker_mm_x," << body.surface_ref_marker_mm.x()
            << ",RMD BaseGSurfacePatchRefMarker QP used to place OBJ vertices\n";
        out << "body," << name << ",surface_ref_marker_mm_y," << body.surface_ref_marker_mm.y()
            << ",RMD BaseGSurfacePatchRefMarker QP used to place OBJ vertices\n";
        out << "body," << name << ",surface_ref_marker_mm_z," << body.surface_ref_marker_mm.z()
            << ",RMD BaseGSurfacePatchRefMarker QP used to place OBJ vertices\n";
    };
    write_body("GEAR21", rmd.gear21);
    write_body("GEAR22", rmd.gear22);

    auto write_joint = [&](const SimpleGearJointInfo& joint) {
        out << "joint," << joint.name << ",id," << joint.id << ",RMD JOINT id\n";
        out << "joint," << joint.name << ",type," << joint.type << ",RMD JOINT type\n";
        out << "joint," << joint.name << ",i_marker_id," << joint.i_marker_id << ",RMD JOINT I marker\n";
        out << "joint," << joint.name << ",j_marker_id," << joint.j_marker_id << ",RMD JOINT J marker\n";
    };
    write_joint(rmd.rev_joint1);
    write_joint(rmd.rev_joint2);

    const auto& motion = rmd.rev_joint1_motion;
    out << "motion," << motion.name << ",joint_id," << motion.joint_id << ",RMD MOTION JOINT\n";
    out << "motion," << motion.name << ",coordinate_type," << motion.coordinate_type << ",RMD MOTION coordinate type\n";
    out << "motion," << motion.name << ",input_type," << motion.input_type << ",RMD MOTION input type\n";
    out << "motion," << motion.name << ",function_constant," << motion.function_constant
        << ",RMD RevJoint1.RMotion FUNCTION used as prescribed GEAR21 angular velocity\n";

    const auto& contact = rmd.geo_surface_contact;
    out << "contact," << contact.name << ",contact_type," << contact.contact_type << ",RMD GGEOMCONTACT CONTYPE\n";
    out << "contact," << contact.name << ",i_float_marker_id," << contact.i_float_marker_id
        << ",RMD GGEOMCONTACT IFLOAT action marker\n";
    out << "contact," << contact.name << ",j_float_marker_id," << contact.j_float_marker_id
        << ",RMD GGEOMCONTACT JFLOAT base marker\n";
    out << "contact," << contact.name << ",i_geometry_id," << contact.i_geometry_id
        << ",RMD GGEOMCONTACT IGGEOMID for GEAR22 surface\n";
    out << "contact," << contact.name << ",j_geometry_id," << contact.j_geometry_id
        << ",RMD GGEOMCONTACT JGGEOMID for GEAR21 surface\n";
    out << "contact," << contact.name << ",i_direction," << contact.i_direction << ",RMD contact side\n";
    out << "contact," << contact.name << ",j_direction," << contact.j_direction << ",RMD contact side\n";
    out << "contact," << contact.name << ",inode_contact," << contact.i_node_contact << ",RMD candidate definition\n";
    out << "contact," << contact.name << ",jnode_contact," << contact.j_node_contact << ",RMD candidate definition\n";
    out << "contact," << contact.name << ",edge_contact," << contact.edge_contact << ",RMD candidate definition\n";
    out << "contact," << contact.name << ",use_cpm," << contact.use_cpm << ",RMD contact point management flag\n";
    out << "contact," << contact.name << ",nomaxcp," << contact.max_contact_points << ",RMD maximum contact point count\n";
    out << "contact," << contact.name << ",cftype," << contact.contact_force_type << ",RMD contact force type\n";
    out << "contact," << contact.name << ",bpen," << contact.bpen
        << ",RMD GGEOMCONTACT penetration transition parameter; parsed for audit\n";
    out << "contact," << contact.name << ",rdf," << contact.rdf << ",RMD restitution/damping factor; parsed for audit\n";
    out << "contact," << contact.name << ",maxpen," << contact.max_pen << ",RMD maximum penetration setting\n";
    out << "contact," << contact.name << ",hparameter," << contact.hparameter << ",RMD h parameter\n";
    out << "contact," << contact.name << ",K," << contact.stiffness << ",RMD penalty stiffness; not used by hard SDF-NCP backend\n";
    out << "contact," << contact.name << ",KORDER," << contact.korder << ",RMD penalty exponent; not used by hard SDF-NCP backend\n";
    out << "contact," << contact.name << ",C," << contact.damping << ",RMD contact damping; not used by hard SDF-NCP backend\n";
}

std::vector<GearReferenceRow> LoadGear22Reference(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        return {};
    }
    std::string line;
    std::getline(in, line);
    std::vector<GearReferenceRow> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto cols = SplitCsvLine(line);
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

GearReferenceRow FindNearestGearReference(const std::vector<GearReferenceRow>& rows, double time) {
    if (rows.empty()) {
        return {};
    }
    return *std::min_element(rows.begin(), rows.end(), [time](const GearReferenceRow& a, const GearReferenceRow& b) {
        return std::abs(a.time - time) < std::abs(b.time - time);
    });
}

GearReferenceRow InterpolateGearReference(const std::vector<GearReferenceRow>& rows, double time) {
    if (rows.empty()) {
        return {};
    }
    if (time <= rows.front().time) {
        return rows.front();
    }
    if (time >= rows.back().time) {
        return rows.back();
    }
    auto upper = std::lower_bound(rows.begin(), rows.end(), time, [](const GearReferenceRow& row, double value) {
        return row.time < value;
    });
    if (upper == rows.begin()) {
        return *upper;
    }
    const auto lower = upper - 1;
    const double denom = upper->time - lower->time;
    const double alpha = std::abs(denom) > 1.0e-14 ? (time - lower->time) / denom : 0.0;
    GearReferenceRow row;
    row.time = time;
    row.omega_rx = lower->omega_rx + alpha * (upper->omega_rx - lower->omega_rx);
    row.alpha_rx = lower->alpha_rx + alpha * (upper->alpha_rx - lower->alpha_rx);
    return row;
}

struct GearReferenceComparisonStats {
    int samples = 0;
    double sum_sq_omega = 0.0;
    double sum_sq_alpha = 0.0;
    double max_abs_omega = 0.0;
    double max_abs_alpha = 0.0;
    double max_abs_ref_omega = 0.0;
    GearReferenceRow final_ref;
    double final_omega = 0.0;
    double final_alpha = 0.0;
};

struct GearAnalyticComparisonStats {
    int samples = 0;
    double sum_sq_omega = 0.0;
    double sum_abs_omega = 0.0;
    double max_abs_omega = 0.0;
    double target_omega = 0.0;
    double final_omega = 0.0;
    double final_abs_error = 0.0;
};

ChSdfTriangleMeshData LoadGearObjAsBodyLocal(const std::filesystem::path& path,
                                             const ChVector3d& surface_ref_marker_mm,
                                             const Mat3& surface_ref_rotation,
                                             const ChVector3d& body_cm_mm,
                                             const Mat3& body_initial_rotation) {
    ChSdfTriangleMeshData mesh;
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
            const ChVector3d world_mm = surface_ref_marker_mm + surface_ref_rotation * ChVector3d(x, y, z);
            const ChVector3d body_local_mm = Transpose(body_initial_rotation) * (world_mm - body_cm_mm);
            mesh.vertices.emplace_back(body_local_mm * kMmToM);
        } else if (tag == "f") {
            std::vector<int> ids;
            std::string token;
            while (iss >> token) {
                const size_t slash = token.find('/');
                std::string v = slash == std::string::npos ? token : token.substr(0, slash);
                int idx = std::stoi(v);
                if (idx < 0) {
                    idx = static_cast<int>(mesh.vertices.size()) + idx + 1;
                }
                ids.push_back(idx - 1);
            }
            for (size_t i = 1; i + 1 < ids.size(); i++) {
                mesh.faces.push_back(chrono::sdfcontact::ChSdfTriangleFace{ids[0], ids[i], ids[i + 1]});
            }
        }
    }
    return mesh;
}

Mat3 GearBodyRotation(const GearPose& pose, double theta_rx) {
    return Multiply(RotXMatrix(theta_rx), pose.initial_rotation);
}

struct SimpleGearConfig {
    std::string case_name = "simple_gear";
    std::string asset = "assets/simple_gear";
    double dt = 0.001;
    double t_end = 1.0;
    double eps = 1.0e-7;
    double tolerance = 1.0e-8;
    double relaxed_tolerance = 1.0e-5;
    int max_iterations = 40;
    double activation_band = 6.25e-5;
    double max_activation_band = 1.25e-4;
    int min_patch_samples = 12;
    double omega21 = 1.0;
    double sdf_ncp_contact_offset = 0.0;
    double sdf_ncp_contact_offset_scale = 0.0;
    double sdf_ncp_gap_scale = 1.0;
    double sdf_ncp_pressure_scale_override = 0.0;
    double sdf_ncp_pressure_scale_band_factor = 4.8;
    double active_band_hysteresis = 1.0e-5;
    bool use_impulse_mixed_ncp = false;
    double impulse_mixed_beta = 0.3;
    double impulse_mixed_cfm = 1.0e-3;
    double impulse_mixed_velocity_scale = 1.0;
    double impulse_mixed_impulse_scale = 1.0;
    bool use_semismooth_active_set = true;
    int max_active_set_iterations = 3;
    double voxel_size = 2.5e-5;
    float half_width_voxels = 20.0f;
};

struct GearState {
    double theta21 = 0.0;
    double theta22 = 0.0;
    double omega22 = 0.0;
};

struct GearContactSampleRef {
    int direction = 0;  // 0: GEAR22 surface -> GEAR21 SDF, 1: GEAR21 surface -> GEAR22 SDF.
    int sample_id = -1;
    int patch_id = -1;
    int persistent_id = -1;
    double area = 0.0;
    double patch_area = 0.0;
    double weight_in_patch = 0.0;
};

struct GearPatchMemory {
    std::set<int> active_sample_ids;
    std::map<int, double> warm_start_intensity;
    std::map<int, std::set<int>> patch_sample_ids;
    int next_patch_id = 0;
};

int GearSamplePersistentId(int direction, int sample_id) {
    return direction == 0 ? sample_id : 1000000 + sample_id;
}

double PersistentSampleSetOverlap(const std::set<int>& a, const std::set<int>& b) {
    if (a.empty() || b.empty()) {
        return 0.0;
    }
    size_t intersection = 0;
    auto ia = a.begin();
    auto ib = b.begin();
    while (ia != a.end() && ib != b.end()) {
        if (*ia == *ib) {
            intersection++;
            ++ia;
            ++ib;
        } else if (*ia < *ib) {
            ++ia;
        } else {
            ++ib;
        }
    }
    const size_t union_size = a.size() + b.size() - intersection;
    return union_size > 0 ? static_cast<double>(intersection) / static_cast<double>(union_size) : 0.0;
}

enum class SimpleGearBackendMode {
    SdfNcp,
    RecurDynGGeomContactLaw,
};

struct SimpleGearGGeomContactLawSI {
    double stiffness = 0.0;
    double damping = 0.0;
    double exponent = 1.0;
    double boundary_penetration = 0.0;
    double rebound_damping_factor = 0.0;
    double static_transition_velocity = 0.0;
    double dynamic_transition_velocity = 0.0;
    double dynamic_friction_coefficient = 0.0;
    double static_friction_coefficient = 0.0;
};

SimpleGearGGeomContactLawSI MakeSimpleGearGGeomContactLawSI(const SimpleGearGGeomContactInfo& contact) {
    constexpr double length_scale_m = 1.0e-3;  // simple_gear RMD declares LENGTH = Millimeter.
    SimpleGearGGeomContactLawSI law;
    law.exponent = contact.korder > 0.0 ? contact.korder : 1.0;
    law.boundary_penetration = contact.bpen * length_scale_m;
    law.rebound_damping_factor = contact.rdf;
    law.dynamic_friction_coefficient = contact.dynamic_friction_coefficient;
    law.static_transition_velocity = contact.static_transition_velocity * length_scale_m;
    law.dynamic_transition_velocity = contact.dynamic_transition_velocity * length_scale_m;
    law.static_friction_coefficient = contact.static_friction_coefficient;
    law.stiffness = contact.stiffness * std::pow(1.0 / length_scale_m, law.exponent);
    law.damping = contact.damping * (1.0 / length_scale_m);
    return law;
}

MultiContactCandidate QueryGear22SampleAgainstGear21(const ChOpenVdbSdfGrid& gear21_sdf,
                                                     const ChSdfContactSurfaceGraph& gear22_graph,
                                                     const GearPose& pose21,
                                                     const GearPose& pose22,
                                                     int sample_id,
                                                     double theta21,
                                                     double theta22) {
    const auto& sample = gear22_graph.samples.at(static_cast<size_t>(sample_id));
    const Mat3 r21 = GearBodyRotation(pose21, theta21);
    const Mat3 r22 = GearBodyRotation(pose22, theta22);
    const Mat3 r21_t = Transpose(r21);
    const ChVector3d r22_world = r22 * sample.local_pos;
    const ChVector3d world_pos = pose22.center + r22_world;
    const ChVector3d gear21_local = r21_t * (world_pos - pose21.center);
    ChSdfContactSampleQuery query = gear21_sdf.QueryLocal(gear21_local);
    const ChVector3d normal_world = NormalizeSdfVector(r21 * query.grad, "Gear SDF gradient must be nonzero.");

    MultiContactCandidate candidate;
    candidate.sample_id = sample_id;
    candidate.area = sample.area;
    candidate.gap = query.phi;
    candidate.normal = normal_world;
    candidate.world_pos = world_pos;
    candidate.jacobian = r22_world.Cross(normal_world).x();
    const ChVector3d dr22_world_dtheta = UnitXCross(r22_world);
    const ChVector3d dgear21_local_dtheta = r21_t * dr22_world_dtheta;
    const ChVector3d dnormal_world_dtheta = r21 * query.UnitGradDirectionalDerivative(dgear21_local_dtheta);
    candidate.jacobian_theta22_derivative =
        dr22_world_dtheta.Cross(normal_world).x() + r22_world.Cross(dnormal_world_dtheta).x();
    return candidate;
}

MultiContactCandidate QueryGear21SampleAgainstGear22(const ChOpenVdbSdfGrid& gear22_sdf,
                                                     const ChSdfContactSurfaceGraph& gear21_graph,
                                                     const GearPose& pose21,
                                                     const GearPose& pose22,
                                                     int sample_id,
                                                     double theta21,
                                                     double theta22) {
    const auto& sample = gear21_graph.samples.at(static_cast<size_t>(sample_id));
    const Mat3 r21 = GearBodyRotation(pose21, theta21);
    const Mat3 r22 = GearBodyRotation(pose22, theta22);
    const Mat3 r22_t = Transpose(r22);
    const ChVector3d r21_world = r21 * sample.local_pos;
    const ChVector3d world_pos = pose21.center + r21_world;
    const ChVector3d gear22_local = r22_t * (world_pos - pose22.center);
    ChSdfContactSampleQuery query = gear22_sdf.QueryLocal(gear22_local);
    const ChVector3d normal_world = NormalizeSdfVector(r22 * query.grad, "Gear SDF gradient must be nonzero.");
    const ChVector3d r22_contact_world = world_pos - pose22.center;

    MultiContactCandidate candidate;
    candidate.sample_id = 1000000 + sample_id;
    candidate.area = sample.area;
    candidate.gap = query.phi;
    candidate.normal = normal_world;
    candidate.world_pos = world_pos;
    candidate.jacobian = -r22_contact_world.Cross(normal_world).x();
    const ChVector3d omega22_local = r22_t * ChVector3d(1, 0, 0);
    const ChVector3d dgear22_local_dtheta = omega22_local.Cross(gear22_local) * -1.0;
    const ChVector3d dnormal_local_dtheta = query.UnitGradDirectionalDerivative(dgear22_local_dtheta);
    const ChVector3d dnormal_world_dtheta = UnitXCross(normal_world) + r22 * dnormal_local_dtheta;
    candidate.jacobian_theta22_derivative = -r22_contact_world.Cross(dnormal_world_dtheta).x();
    return candidate;
}

double SampleGear22PhiAgainstGear21(const ChOpenVdbSdfGrid& gear21_sdf,
                                    const ChSdfContactSurfaceGraph& gear22_graph,
                                    const GearPose& pose21,
                                    const GearPose& pose22,
                                    int sample_id,
                                    double theta21,
                                    double theta22) {
    const auto& sample = gear22_graph.samples.at(static_cast<size_t>(sample_id));
    const Mat3 r21 = GearBodyRotation(pose21, theta21);
    const Mat3 r22 = GearBodyRotation(pose22, theta22);
    const ChVector3d world_pos = pose22.center + r22 * sample.local_pos;
    const ChVector3d gear21_local = Transpose(r21) * (world_pos - pose21.center);
    return gear21_sdf.SamplePhi(gear21_local);
}

double SampleGear21PhiAgainstGear22(const ChOpenVdbSdfGrid& gear22_sdf,
                                    const ChSdfContactSurfaceGraph& gear21_graph,
                                    const GearPose& pose21,
                                    const GearPose& pose22,
                                    int sample_id,
                                    double theta21,
                                    double theta22) {
    const auto& sample = gear21_graph.samples.at(static_cast<size_t>(sample_id));
    const Mat3 r21 = GearBodyRotation(pose21, theta21);
    const Mat3 r22 = GearBodyRotation(pose22, theta22);
    const ChVector3d world_pos = pose21.center + r21 * sample.local_pos;
    const ChVector3d gear22_local = Transpose(r22) * (world_pos - pose22.center);
    return gear22_sdf.SamplePhi(gear22_local);
}

double SimpleGearJacobianTheta22Derivative(const ChOpenVdbSdfGrid& gear21_sdf,
                                           const ChOpenVdbSdfGrid& gear22_sdf,
                                           const ChSdfContactSurfaceGraph& gear21_graph,
                                           const ChSdfContactSurfaceGraph& gear22_graph,
                                           const GearPose& pose21,
                                           const GearPose& pose22,
                                           const GearContactSampleRef& ref,
                                           double theta21,
                                           double theta22) {
    const auto candidate = ref.direction == 0 ?
                               QueryGear22SampleAgainstGear21(
                                   gear21_sdf, gear22_graph, pose21, pose22, ref.sample_id, theta21, theta22) :
                               QueryGear21SampleAgainstGear22(
                                   gear22_sdf, gear21_graph, pose21, pose22, ref.sample_id, theta21, theta22);
    return candidate.jacobian_theta22_derivative;
}

std::vector<int> BuildGearPatchSamples(const ChOpenVdbSdfGrid& gear21_sdf,
                                       const ChSdfContactSurfaceGraph& gear22_graph,
                                       const GearPose& pose21,
                                       const GearPose& pose22,
                                       double theta21,
                                       double theta22,
                                       double activation_band,
                                       double max_activation_band,
                                       int min_patch_samples,
                                       int* patch_count = nullptr) {
    std::vector<std::pair<int, double>> sampled;
    sampled.reserve(gear22_graph.samples.size());
    double min_phi = std::numeric_limits<double>::max();
    for (const auto& sample : gear22_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        const auto candidate = QueryGear22SampleAgainstGear21(
            gear21_sdf, gear22_graph, pose21, pose22, sample.id, theta21, theta22);
        sampled.emplace_back(sample.id, candidate.gap);
        min_phi = std::min(min_phi, candidate.gap);
    }
    if (min_phi > std::max(activation_band, max_activation_band)) {
        if (patch_count) {
            *patch_count = 0;
        }
        return {};
    }

    std::vector<int> active;
    active.reserve(gear22_graph.samples.size());
    double band = std::max(activation_band, 0.0);
    const double band_limit = std::max(band, max_activation_band);
    for (;;) {
        active.clear();
        const double threshold = min_phi + band;
        for (const auto& item : sampled) {
            if (item.second <= threshold) {
                active.push_back(item.first);
            }
        }
        if (static_cast<int>(active.size()) >= std::max(1, min_patch_samples) || band >= band_limit) {
            break;
        }
        band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
    }
    const auto components = gear22_graph.FindConnectedComponents(active);
    if (patch_count) {
        *patch_count = static_cast<int>(components.size());
    }

    std::vector<int> ids;
    for (const auto& component : components) {
        ids.insert(ids.end(), component.begin(), component.end());
    }
    return ids;
}

std::vector<int> BuildGear21PatchSamplesAgainstGear22(const ChOpenVdbSdfGrid& gear22_sdf,
                                                      const ChSdfContactSurfaceGraph& gear21_graph,
                                                      const GearPose& pose21,
                                                      const GearPose& pose22,
                                                      double theta21,
                                                      double theta22,
                                                      double activation_band,
                                                      double max_activation_band,
                                                      int min_patch_samples,
                                                      int* patch_count = nullptr) {
    std::vector<std::pair<int, double>> sampled;
    sampled.reserve(gear21_graph.samples.size());
    double min_phi = std::numeric_limits<double>::max();
    for (const auto& sample : gear21_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        const auto candidate = QueryGear21SampleAgainstGear22(
            gear22_sdf, gear21_graph, pose21, pose22, sample.id, theta21, theta22);
        sampled.emplace_back(sample.id, candidate.gap);
        min_phi = std::min(min_phi, candidate.gap);
    }
    if (min_phi > std::max(activation_band, max_activation_band)) {
        if (patch_count) {
            *patch_count = 0;
        }
        return {};
    }

    std::vector<int> active;
    active.reserve(gear21_graph.samples.size());
    double band = std::max(activation_band, 0.0);
    const double band_limit = std::max(band, max_activation_band);
    for (;;) {
        active.clear();
        const double threshold = min_phi + band;
        for (const auto& item : sampled) {
            if (item.second <= threshold) {
                active.push_back(item.first);
            }
        }
        if (static_cast<int>(active.size()) >= std::max(1, min_patch_samples) || band >= band_limit) {
            break;
        }
        band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
    }
    const auto components = gear21_graph.FindConnectedComponents(active);
    if (patch_count) {
        *patch_count = static_cast<int>(components.size());
    }

    std::vector<int> ids;
    for (const auto& component : components) {
        ids.insert(ids.end(), component.begin(), component.end());
    }
    return ids;
}

std::vector<GearContactSampleRef> BuildGearBidirectionalPatchSamples(const ChOpenVdbSdfGrid& gear21_sdf,
                                                                     const ChOpenVdbSdfGrid& gear22_sdf,
                                                                     const ChSdfContactSurfaceGraph& gear21_graph,
                                                                     const ChSdfContactSurfaceGraph& gear22_graph,
                                                                     const GearPose& pose21,
                                                                     const GearPose& pose22,
                                                                     const SimpleGearConfig& config,
                                                                     double theta21,
                                                                     double theta22,
                                                                     int* patch_count = nullptr,
                                                                     GearPatchMemory* memory = nullptr,
                                                                     const ChSdfContactSampleBvh* gear21_bvh = nullptr,
                                                                     const ChSdfContactSampleBvh* gear22_bvh = nullptr) {
    struct SamplePhi {
        int id = -1;
        double effective_gap = 0.0;
        double area = 0.0;
    };

    std::map<int, std::set<int>> next_patch_sample_ids;
    std::set<int> used_old_patch_ids;

    auto build_direction = [&](int direction, const ChSdfContactSurfaceGraph& graph, const ChSdfContactSampleBvh* bvh) {
        auto all_sample_ids = [&]() {
            std::vector<int> ids;
            ids.reserve(graph.samples.size());
            for (const auto& sample : graph.samples) {
                if (sample.area > 0.0) {
                    ids.push_back(sample.id);
                }
            }
            return ids;
        };

        std::vector<int> candidate_ids;
        if (bvh && !bvh->Empty()) {
            const double broad_padding = std::max(config.activation_band, config.max_activation_band) +
                                         std::abs(config.sdf_ncp_contact_offset) +
                                         std::max(0.0, config.active_band_hysteresis) +
                                         8.0 * std::max(config.voxel_size, 1.0e-12);
            const Mat3 r21 = GearBodyRotation(pose21, theta21);
            const Mat3 r22 = GearBodyRotation(pose22, theta22);
            const ChSdfAabb query_bounds =
                direction == 0 ?
                    TransformAabbBetweenBodyFrames(
                        gear21_sdf.local_bounds, r21, pose21.center, r22, pose22.center, broad_padding) :
                    TransformAabbBetweenBodyFrames(
                        gear22_sdf.local_bounds, r22, pose22.center, r21, pose21.center, broad_padding);
            candidate_ids = bvh->Query(query_bounds);
            std::sort(candidate_ids.begin(), candidate_ids.end());
            candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());
        } else {
            candidate_ids = all_sample_ids();
        }

        if (candidate_ids.empty()) {
            return std::vector<GearContactSampleRef>{};
        }

        std::vector<SamplePhi> sampled;
        sampled.reserve(candidate_ids.size());
        auto collect_sampled = [&](const std::vector<int>& ids, double& min_gap) {
            sampled.clear();
            min_gap = std::numeric_limits<double>::max();
            for (int sample_id : ids) {
                if (sample_id < 0 || sample_id >= static_cast<int>(graph.samples.size())) {
                    continue;
                }
                const auto& sample = graph.samples[static_cast<size_t>(sample_id)];
                if (sample.area <= 0.0) {
                    continue;
                }
                const double gap = direction == 0 ?
                                       SampleGear22PhiAgainstGear21(
                                           gear21_sdf, gear22_graph, pose21, pose22, sample.id, theta21, theta22) :
                                       SampleGear21PhiAgainstGear22(
                                           gear22_sdf, gear21_graph, pose21, pose22, sample.id, theta21, theta22);
                const double effective_gap = gap - config.sdf_ncp_contact_offset;
                sampled.push_back(SamplePhi{sample.id, effective_gap, sample.area});
                min_gap = std::min(min_gap, effective_gap);
            }
        };

        double min_gap = std::numeric_limits<double>::max();
        collect_sampled(candidate_ids, min_gap);

        std::vector<int> active;
        auto build_active = [&]() {
            active.clear();
            if (sampled.empty()) {
                return;
            }
            if (min_gap <= std::max(config.activation_band, config.max_activation_band)) {
                double band = std::max(config.activation_band, 0.0);
                const double band_limit = std::max(band, config.max_activation_band);
                for (;;) {
                    active.clear();
                    const double threshold = min_gap + band;
                    const double retained_threshold = min_gap + band + std::max(0.0, config.active_band_hysteresis);
                    for (const auto& item : sampled) {
                        const int persistent_id = GearSamplePersistentId(direction, item.id);
                        const bool retained =
                            memory && memory->active_sample_ids.count(persistent_id) > 0 &&
                            item.effective_gap <= retained_threshold;
                        if (item.effective_gap <= threshold || retained) {
                            active.push_back(item.id);
                        }
                    }
                    if (static_cast<int>(active.size()) >= std::max(1, config.min_patch_samples) ||
                        band >= band_limit) {
                        break;
                    }
                    band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
                }
            } else if (memory) {
                const double retained_threshold =
                    min_gap + config.max_activation_band + std::max(0.0, config.active_band_hysteresis);
                for (const auto& item : sampled) {
                    const int persistent_id = GearSamplePersistentId(direction, item.id);
                    if (memory->active_sample_ids.count(persistent_id) > 0 && item.effective_gap <= retained_threshold) {
                        active.push_back(item.id);
                    }
                }
            }
        };

        if (sampled.empty()) {
            return std::vector<GearContactSampleRef>{};
        }
        build_active();

        if (bvh && active.size() < static_cast<size_t>(std::max(1, config.min_patch_samples)) &&
            candidate_ids.size() < graph.samples.size()) {
            candidate_ids = all_sample_ids();
            sampled.reserve(candidate_ids.size());
            collect_sampled(candidate_ids, min_gap);
            build_active();
        }

        const auto components = graph.FindConnectedComponents(active);
        std::vector<GearContactSampleRef> refs;
        for (const auto& component : components) {
            if (component.empty()) {
                continue;
            }
            double patch_area = 0.0;
            int min_sample_id = component.front();
            for (int sample_id : component) {
                min_sample_id = std::min(min_sample_id, sample_id);
                if (sample_id >= 0 && sample_id < static_cast<int>(graph.samples.size())) {
                    patch_area += std::max(0.0, graph.samples[static_cast<size_t>(sample_id)].area);
                }
            }
            int patch_id = direction == 0 ? min_sample_id : 1000000 + min_sample_id;
            std::set<int> persistent_component_ids;
            for (int sample_id : component) {
                persistent_component_ids.insert(GearSamplePersistentId(direction, sample_id));
            }
            if (memory) {
                int best_patch_id = -1;
                double best_overlap = 0.0;
                for (const auto& old_patch : memory->patch_sample_ids) {
                    if (used_old_patch_ids.count(old_patch.first) > 0) {
                        continue;
                    }
                    const double overlap = PersistentSampleSetOverlap(persistent_component_ids, old_patch.second);
                    if (overlap > best_overlap) {
                        best_overlap = overlap;
                        best_patch_id = old_patch.first;
                    }
                }
                if (best_patch_id >= 0 && best_overlap >= 0.20) {
                    patch_id = best_patch_id;
                    used_old_patch_ids.insert(best_patch_id);
                } else {
                    patch_id = memory->next_patch_id++;
                }
                next_patch_sample_ids[patch_id] = persistent_component_ids;
            }
            for (int sample_id : component) {
                if (sample_id < 0 || sample_id >= static_cast<int>(graph.samples.size())) {
                    continue;
                }
                const double area = std::max(0.0, graph.samples[static_cast<size_t>(sample_id)].area);
                GearContactSampleRef ref;
                ref.direction = direction;
                ref.sample_id = sample_id;
                ref.patch_id = patch_id;
                ref.persistent_id = GearSamplePersistentId(direction, sample_id);
                ref.area = area;
                ref.patch_area = patch_area;
                ref.weight_in_patch = patch_area > 1.0e-30 ? area / patch_area : 0.0;
                refs.push_back(ref);
            }
        }
        return refs;
    };

    std::vector<GearContactSampleRef> refs22 = build_direction(0, gear22_graph, gear22_bvh);
    std::vector<GearContactSampleRef> refs21 = build_direction(1, gear21_graph, gear21_bvh);
    std::vector<GearContactSampleRef> refs;
    refs.reserve(refs22.size() + refs21.size());
    refs.insert(refs.end(), refs22.begin(), refs22.end());
    refs.insert(refs.end(), refs21.begin(), refs21.end());

    if (patch_count) {
        std::set<int> patch_ids;
        for (const auto& ref : refs) {
            patch_ids.insert(ref.patch_id);
        }
        *patch_count = static_cast<int>(patch_ids.size());
    }
    if (memory) {
        memory->active_sample_ids.clear();
        for (const auto& ref : refs) {
            memory->active_sample_ids.insert(ref.persistent_id);
        }
        memory->patch_sample_ids = std::move(next_patch_sample_ids);
    }
    return refs;
}

bool SameGearSampleRefs(const std::vector<GearContactSampleRef>& a, const std::vector<GearContactSampleRef>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].persistent_id != b[i].persistent_id || a[i].patch_id != b[i].patch_id) {
            return false;
        }
    }
    return true;
}

std::vector<double> ComputeGearPatchAreaWeights(const ChSdfContactSurfaceGraph& gear21_graph,
                                                const ChSdfContactSurfaceGraph& gear22_graph,
                                                const std::vector<GearContactSampleRef>& sample_refs) {
    std::vector<double> weights(sample_refs.size(), 0.0);
    for (size_t i = 0; i < sample_refs.size(); i++) {
        const auto& ref = sample_refs[i];
        if (ref.patch_area > 1.0e-30) {
            weights[i] = ref.area / ref.patch_area;
        } else if (ref.direction == 0 && ref.sample_id >= 0 &&
                   ref.sample_id < static_cast<int>(gear22_graph.samples.size())) {
            weights[i] = std::max(0.0, gear22_graph.samples[static_cast<size_t>(ref.sample_id)].area);
        } else if (ref.direction == 1 && ref.sample_id >= 0 &&
                   ref.sample_id < static_cast<int>(gear21_graph.samples.size())) {
            weights[i] = std::max(0.0, gear21_graph.samples[static_cast<size_t>(ref.sample_id)].area);
        }
    }
    return weights;
}

double MinGearGap(const ChOpenVdbSdfGrid& gear21_sdf,
                  const ChSdfContactSurfaceGraph& gear22_graph,
                  const GearPose& pose21,
                  const GearPose& pose22,
                  double theta21,
                  double theta22) {
    double min_gap = std::numeric_limits<double>::max();
    for (const auto& sample : gear22_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        min_gap = std::min(min_gap,
                           SampleGear22PhiAgainstGear21(
                               gear21_sdf, gear22_graph, pose21, pose22, sample.id, theta21, theta22));
    }
    return min_gap;
}

double MinGearGapReverse(const ChOpenVdbSdfGrid& gear22_sdf,
                         const ChSdfContactSurfaceGraph& gear21_graph,
                         const GearPose& pose21,
                         const GearPose& pose22,
                         double theta21,
                         double theta22) {
    double min_gap = std::numeric_limits<double>::max();
    for (const auto& sample : gear21_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        min_gap = std::min(min_gap,
                           SampleGear21PhiAgainstGear22(
                               gear22_sdf, gear21_graph, pose21, pose22, sample.id, theta21, theta22));
    }
    return min_gap;
}

double MinGearGapBidirectional(const ChOpenVdbSdfGrid& gear21_sdf,
                               const ChOpenVdbSdfGrid& gear22_sdf,
                               const ChSdfContactSurfaceGraph& gear21_graph,
                               const ChSdfContactSurfaceGraph& gear22_graph,
                               const GearPose& pose21,
                               const GearPose& pose22,
                               double theta21,
                               double theta22) {
    return std::min(MinGearGap(gear21_sdf, gear22_graph, pose21, pose22, theta21, theta22),
                    MinGearGapReverse(gear22_sdf, gear21_graph, pose21, pose22, theta21, theta22));
}

double SimpleGearDriverGapJacobian(const MultiContactCandidate& candidate, const GearPose& pose21, int direction);

double EstimateSimpleGearPressureScale(const ChOpenVdbSdfGrid& gear21_sdf,
                                       const ChOpenVdbSdfGrid& gear22_sdf,
                                       const ChSdfContactSurfaceGraph& gear21_graph,
                                       const ChSdfContactSurfaceGraph& gear22_graph,
                                       const GearPose& pose21,
                                       const GearPose& pose22,
                                       const SimpleGearConfig& config,
                                       double inertia22,
                                       const GearState& state,
                                       const std::vector<GearContactSampleRef>& sample_refs) {
    if (config.sdf_ncp_pressure_scale_override > 0.0) {
        return config.sdf_ncp_pressure_scale_override;
    }

    const double theta21_next = state.theta21 + config.dt * config.omega21;
    const double theta22_next = state.theta22 + config.dt * state.omega22;
    double area_moment = 0.0;
    for (const auto& ref : sample_refs) {
        const auto candidate = ref.direction == 0 ?
                                   QueryGear22SampleAgainstGear21(gear21_sdf,
                                                                  gear22_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_next,
                                                                  theta22_next) :
                                   QueryGear21SampleAgainstGear22(gear22_sdf,
                                                                  gear21_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_next,
                                                                  theta22_next);
        area_moment += std::max(0.0, ref.area) * std::abs(candidate.jacobian);
    }

    if (!(area_moment > 1.0e-18) || !std::isfinite(area_moment)) {
        return 1.0;
    }

    const double velocity_scale =
        std::max({1.0, std::abs(config.omega21), std::abs(state.omega22), std::abs(config.omega21 - state.omega22)});
    const double torque_scale = inertia22 * velocity_scale / std::max(config.dt, 1.0e-12);
    const double pressure_scale = torque_scale / area_moment;
    return std::min(1.0e9, std::max(1.0, pressure_scale));
}

double EstimateSimpleGearImpulseIntensityScale(const ChOpenVdbSdfGrid& gear21_sdf,
                                               const ChOpenVdbSdfGrid& gear22_sdf,
                                               const ChSdfContactSurfaceGraph& gear21_graph,
                                               const ChSdfContactSurfaceGraph& gear22_graph,
                                               const GearPose& pose21,
                                               const GearPose& pose22,
                                               const SimpleGearConfig& config,
                                               double inertia22,
                                               const GearState& state,
                                               const std::vector<GearContactSampleRef>& sample_refs) {
    const double theta21_next = state.theta21 + config.dt * config.omega21;
    const double theta22_next = state.theta22 + config.dt * state.omega22;
    double area_moment = 0.0;
    for (const auto& ref : sample_refs) {
        const auto candidate = ref.direction == 0 ?
                                   QueryGear22SampleAgainstGear21(gear21_sdf,
                                                                  gear22_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_next,
                                                                  theta22_next) :
                                   QueryGear21SampleAgainstGear22(gear22_sdf,
                                                                  gear21_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_next,
                                                                  theta22_next);
        area_moment += std::max(0.0, ref.area) * std::abs(candidate.jacobian);
    }

    if (!(area_moment > 1.0e-18) || !std::isfinite(area_moment)) {
        return 1.0;
    }

    const double velocity_scale = std::max({config.impulse_mixed_velocity_scale,
                                            std::abs(config.omega21),
                                            std::abs(state.omega22),
                                            std::abs(config.omega21 - state.omega22),
                                            1.0});
    const double angular_momentum_scale = inertia22 * velocity_scale;
    const double intensity_scale = angular_momentum_scale / area_moment;
    return std::min(1.0e8, std::max(1.0e-8, intensity_scale));
}

ChSdfNcpGeneralizedProblem MakeSimpleGearGeneralizedProblem(const ChOpenVdbSdfGrid& gear21_sdf,
                                                            const ChOpenVdbSdfGrid& gear22_sdf,
                                                            const ChSdfContactSurfaceGraph& gear21_graph,
                                                            const ChSdfContactSurfaceGraph& gear22_graph,
                                                            const GearPose& pose21,
                                                            const GearPose& pose22,
                                                            const SimpleGearConfig& config,
                                                            double inertia22,
                                                            const GearState& state,
                                                            const std::vector<GearContactSampleRef>& sample_refs,
                                                            int* patch_count = nullptr) {
    if (patch_count) {
        std::set<int> patch_ids;
        for (const auto& ref : sample_refs) {
            patch_ids.insert(ref.patch_id);
        }
        *patch_count = static_cast<int>(patch_ids.size());
    }

    ChSdfNcpGeneralizedProblem problem;
    problem.dt = config.dt;
    problem.eps = config.eps;
    const double pressure_scale = EstimateSimpleGearPressureScale(gear21_sdf,
                                                                  gear22_sdf,
                                                                  gear21_graph,
                                                                  gear22_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  config,
                                                                  inertia22,
                                                                  state,
                                                                  sample_refs);
    problem.gap_scale = config.sdf_ncp_gap_scale > 0.0 ? config.sdf_ncp_gap_scale : 1.0;
    problem.lambda_scale = config.use_impulse_mixed_ncp ? config.impulse_mixed_impulse_scale : 1.0;
    problem.tolerance = config.tolerance;
    problem.relaxed_tolerance = config.relaxed_tolerance;
    problem.negative_lambda_tolerance = 1.0e-8;
    problem.max_iterations = config.max_iterations;
    problem.use_analytic_jacobian = false;
    problem.current_velocity = {state.omega22};
    problem.mass_diagonal = {inertia22};
    problem.external_force = {0.0};
    problem.contact_count = sample_refs.size();
    problem.evaluate_contacts = [&, sample_refs, pressure_scale](const std::vector<double>& v_next) {
        const double omega22_next = v_next.at(0);
        const double theta21_next = state.theta21 + config.dt * config.omega21;
        const double theta22_next = state.theta22 + config.dt * omega22_next;
        std::vector<ChSdfNcpGeneralizedContact> contacts;
        contacts.reserve(sample_refs.size());
        for (size_t i = 0; i < sample_refs.size(); i++) {
            const auto& ref = sample_refs[i];
            const auto candidate = ref.direction == 0 ?
                                       QueryGear22SampleAgainstGear21(gear21_sdf,
                                                                      gear22_graph,
                                                                      pose21,
                                                                      pose22,
                                                                      ref.sample_id,
                                                                      theta21_next,
                                                                      theta22_next) :
                                       QueryGear21SampleAgainstGear22(gear22_sdf,
                                                                      gear21_graph,
                                                                      pose21,
                                                                      pose22,
                                                                      ref.sample_id,
                                                                      theta21_next,
                                                                      theta22_next);

            ChSdfNcpGeneralizedContact contact;
            contact.gap = candidate.gap - config.sdf_ncp_contact_offset;
            contact.jacobian = {candidate.jacobian};
            contact.normal_velocity_offset = -SimpleGearDriverGapJacobian(candidate, pose21, ref.direction) *
                                             config.omega21;
            if (config.use_impulse_mixed_ncp) {
                const double dJ_dtheta22 = SimpleGearJacobianTheta22Derivative(gear21_sdf,
                                                                               gear22_sdf,
                                                                               gear21_graph,
                                                                               gear22_graph,
                                                                               pose21,
                                                                               pose22,
                                                                               ref,
                                                                               theta21_next,
                                                                               theta22_next);
                contact.jacobian_velocity_derivative = {config.dt * dJ_dtheta22};
            }
            contact.weight = ref.area * (config.use_impulse_mixed_ncp ? config.dt * pressure_scale : pressure_scale);
            contact.contact_id = ref.persistent_id;
            contact.patch_id = ref.patch_id;
            contacts.push_back(contact);
        }
        return contacts;
    };
    return problem;
}

std::vector<double> ComputeSimpleGearResidual(const ChOpenVdbSdfGrid& gear21_sdf,
                                              const ChOpenVdbSdfGrid& gear22_sdf,
                                              const ChSdfContactSurfaceGraph& gear21_graph,
                                              const ChSdfContactSurfaceGraph& gear22_graph,
                                              const GearPose& pose21,
                                              const GearPose& pose22,
                                              const SimpleGearConfig& config,
                                              double inertia22,
                                              const GearState& state,
                                              const std::vector<GearContactSampleRef>& sample_refs,
                                              const std::vector<double>& z,
                                              MultiStepDiagnostics* diag = nullptr) {
    const size_t n = sample_refs.size();
    if (z.size() != n + 1) {
        throw std::invalid_argument("Simple gear SDF-NCP residual has inconsistent dimension.");
    }
    int patch_count = 0;
    const auto problem = MakeSimpleGearGeneralizedProblem(gear21_sdf,
                                                          gear22_sdf,
                                                          gear21_graph,
                                                          gear22_graph,
                                                          pose21,
                                                          pose22,
                                                          config,
                                                          inertia22,
                                                          state,
                                                          sample_refs,
                                                          &patch_count);
    ChSdfNcpGeneralizedResidual residual;
    ChSdfNcpGeneralizedDiagnostics generalized_diag;
    if (config.use_impulse_mixed_ncp) {
        ChSdfNcpImpulseMixedSettings mixed_settings;
        mixed_settings.beta = config.impulse_mixed_beta;
        mixed_settings.cfm = config.impulse_mixed_cfm;
        mixed_settings.velocity_scale = config.impulse_mixed_velocity_scale;
        mixed_settings.impulse_scale = config.impulse_mixed_impulse_scale;
        const auto mixed = ComputeSdfNcpGeneralizedImpulseMixedAd(problem, z, mixed_settings);
        residual.value = mixed.value;
        residual.contacts = mixed.contacts;
        generalized_diag =
            MakeSdfNcpGeneralizedImpulseMixedDiagnostics(problem,
                                                         z,
                                                         mixed_settings,
                                                         true,
                                                         0,
                                                         "evaluated through generalized impulse mixed SDF-NCP backend");
    } else {
        residual = ComputeSdfNcpGeneralizedResidual(problem, z);
        generalized_diag =
            MakeSdfNcpGeneralizedDiagnostics(problem, z, true, 0, "evaluated through generalized SDF-NCP backend");
    }
    if (diag) {
        PopulateMultiStepDiagnostics(residual.contacts, generalized_diag, *diag);
        diag->patch_count = patch_count;
    }
    return residual.value;
}

struct GearStepResult {
    GearState state;
    MultiStepDiagnostics diagnostics;
};

GearStepResult SolveSimpleGearStep(const ChOpenVdbSdfGrid& gear21_sdf,
                                   const ChOpenVdbSdfGrid& gear22_sdf,
                                   const ChSdfContactSurfaceGraph& gear21_graph,
                                   const ChSdfContactSurfaceGraph& gear22_graph,
                                   const GearPose& pose21,
                                   const GearPose& pose22,
                                   const SimpleGearConfig& config,
                                   double inertia22,
                                   const GearState& state,
                                   GearPatchMemory* memory = nullptr,
                                   const ChSdfContactSampleBvh* gear21_bvh = nullptr,
                                   const ChSdfContactSampleBvh* gear22_bvh = nullptr) {
    GearStepResult result;
    result.state = state;
    const double theta21_free = state.theta21 + config.dt * config.omega21;
    const double theta22_free = state.theta22 + config.dt * state.omega22;
    GearPatchMemory active_memory = memory ? *memory : GearPatchMemory{};
    int patch_count = 0;
    std::vector<GearContactSampleRef> sample_refs = BuildGearBidirectionalPatchSamples(gear21_sdf,
                                                                                       gear22_sdf,
                                                                                       gear21_graph,
                                                                                       gear22_graph,
                                                                                       pose21,
                                                                                       pose22,
                                                                                       config,
                                                                                       theta21_free,
                                                                                       theta22_free,
                                                                                       &patch_count,
                                                                                       &active_memory,
                                                                                       gear21_bvh,
                                                                                       gear22_bvh);
    std::map<int, double> warm_start = memory ? memory->warm_start_intensity : std::map<int, double>{};
    ChSdfNcpGeneralizedStepResult solved;
    const int active_set_limit =
        config.use_semismooth_active_set ? std::max(1, config.max_active_set_iterations) : 1;

    for (int active_iter = 0; active_iter < active_set_limit; active_iter++) {
        std::vector<double> z(1 + sample_refs.size(), 0.0);
        z[0] = active_iter == 0 || solved.next_velocity.empty() ? state.omega22 : solved.next_velocity.at(0);
        for (size_t i = 0; i < sample_refs.size(); i++) {
            const auto found = warm_start.find(sample_refs[i].persistent_id);
            if (found != warm_start.end() && std::isfinite(found->second)) {
                z[1 + i] = std::max(0.0, found->second);
            }
        }

        const auto problem = MakeSimpleGearGeneralizedProblem(gear21_sdf,
                                                              gear22_sdf,
                                                              gear21_graph,
                                                              gear22_graph,
                                                              pose21,
                                                              pose22,
                                                              config,
                                                              inertia22,
                                                              state,
                                                              sample_refs,
                                                              &patch_count);
        if (config.use_impulse_mixed_ncp) {
            ChSdfNcpImpulseMixedSettings mixed_settings;
            mixed_settings.beta = config.impulse_mixed_beta;
            mixed_settings.cfm = config.impulse_mixed_cfm;
            mixed_settings.velocity_scale = config.impulse_mixed_velocity_scale;
            mixed_settings.impulse_scale = config.impulse_mixed_impulse_scale;
            solved = SolveSdfNcpGeneralizedImpulseMixedProblem(problem, z, mixed_settings);
        } else {
            solved = SolveSdfNcpGeneralizedProblem(problem, z);
        }

        std::map<int, double> solved_warm_start;
        for (size_t i = 0; i < sample_refs.size() && i < solved.lambdas.size(); i++) {
            solved_warm_start[sample_refs[i].persistent_id] = std::max(0.0, solved.lambdas[i]);
        }

        if (!config.use_semismooth_active_set || active_iter + 1 >= active_set_limit ||
            solved.next_velocity.empty()) {
            warm_start = std::move(solved_warm_start);
            break;
        }

        GearPatchMemory candidate_memory = active_memory;
        candidate_memory.warm_start_intensity = solved_warm_start;
        int next_patch_count = 0;
        const double theta22_solved = state.theta22 + config.dt * solved.next_velocity.at(0);
        std::vector<GearContactSampleRef> next_refs = BuildGearBidirectionalPatchSamples(gear21_sdf,
                                                                                         gear22_sdf,
                                                                                         gear21_graph,
                                                                                         gear22_graph,
                                                                                         pose21,
                                                                                         pose22,
                                                                                         config,
                                                                                         theta21_free,
                                                                                         theta22_solved,
                                                                                         &next_patch_count,
                                                                                         &candidate_memory,
                                                                                         gear21_bvh,
                                                                                         gear22_bvh);
        if (SameGearSampleRefs(sample_refs, next_refs)) {
            active_memory = std::move(candidate_memory);
            warm_start = std::move(solved_warm_start);
            patch_count = next_patch_count;
            break;
        }

        sample_refs = std::move(next_refs);
        warm_start = std::move(solved_warm_start);
        active_memory = std::move(candidate_memory);
        patch_count = next_patch_count;
    }

    if (memory) {
        memory->active_sample_ids = std::move(active_memory.active_sample_ids);
        memory->warm_start_intensity = std::move(warm_start);
    }

    result.state.omega22 = solved.next_velocity.at(0);
    result.state.theta21 = state.theta21 + config.dt * config.omega21;
    result.state.theta22 = state.theta22 + config.dt * result.state.omega22;
    PopulateMultiStepDiagnostics(solved.contacts, solved.diagnostics, result.diagnostics);
    result.diagnostics.patch_count = patch_count;
    result.diagnostics.min_gap = MinGearGapBidirectional(gear21_sdf,
                                                         gear22_sdf,
                                                         gear21_graph,
                                                         gear22_graph,
                                                         pose21,
                                                         pose22,
                                                         result.state.theta21,
                                                         result.state.theta22);
    result.diagnostics.max_penetration =
        std::max(result.diagnostics.max_penetration, std::max(0.0, -result.diagnostics.min_gap));
    return result;
}

double SimpleGearDriverGapJacobian(const MultiContactCandidate& candidate,
                                   const GearPose& pose21,
                                   int direction) {
    const ChVector3d r21_contact_world = candidate.world_pos - pose21.center;
    const double r_cross_n_x = r21_contact_world.Cross(candidate.normal).x();
    return direction == 0 ? -r_cross_n_x : r_cross_n_x;
}

double EvaluateSimpleGearGGeomContactLaw(const ChOpenVdbSdfGrid& gear21_sdf,
                                         const ChOpenVdbSdfGrid& gear22_sdf,
                                         const ChSdfContactSurfaceGraph& gear21_graph,
                                         const ChSdfContactSurfaceGraph& gear22_graph,
                                         const GearPose& pose21,
                                         const GearPose& pose22,
                                         const SimpleGearConfig& config,
                                         const SimpleGearGGeomContactLawSI& law,
                                         const std::vector<GearContactSampleRef>& sample_refs,
                                         double theta21_eval,
                                         double theta22_eval,
                                         double omega22_eval,
                                         MultiStepDiagnostics* diag = nullptr) {
    const std::vector<double> weights = ComputeGearPatchAreaWeights(gear21_graph, gear22_graph, sample_refs);
    double generalized_torque = 0.0;
    if (diag) {
        ResetMultiStepDiagnostics(*diag);
        diag->success = true;
        diag->min_gap = sample_refs.empty() ? 0.0 : std::numeric_limits<double>::max();
        diag->patch_count = 0;
    }

    for (size_t i = 0; i < sample_refs.size(); i++) {
        const auto& ref = sample_refs[i];
        const auto candidate = ref.direction == 0 ?
                                   QueryGear22SampleAgainstGear21(gear21_sdf,
                                                                  gear22_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_eval,
                                                                  theta22_eval) :
                                   QueryGear21SampleAgainstGear22(gear22_sdf,
                                                                  gear21_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  ref.sample_id,
                                                                  theta21_eval,
                                                                   theta22_eval);
        const double driver_jacobian = SimpleGearDriverGapJacobian(candidate, pose21, ref.direction);
        const double gap_rate = driver_jacobian * config.omega21 + candidate.jacobian * omega22_eval;
        const double penetration = std::max(0.0, -candidate.gap);
        const double spring_force = law.stiffness > 0.0 && penetration > 0.0 ?
                                        law.stiffness * std::pow(penetration, law.exponent) :
                                        0.0;
        double damping = std::max(0.0, law.damping);
        if (law.boundary_penetration > 1.0e-14) {
            damping = RecurDynStepRamp(penetration, 0.0, 0.0, law.boundary_penetration, damping);
        }
        const double damping_force = penetration > 0.0 ? damping * std::max(0.0, -gap_rate) : 0.0;
        const double raw_force = std::max(0.0, spring_force + damping_force);
        const double weight = i < weights.size() ? weights[i] : 0.0;
        const double weighted_force = raw_force * weight;
        generalized_torque += candidate.jacobian * weighted_force;

        if (diag) {
            MultiContactCandidate stored = candidate;
            stored.weight = weight;
            diag->candidates.push_back(stored);
            diag->lambdas.push_back(raw_force);
            diag->ncp_residuals.push_back(0.0);
            diag->complementarity_errors.push_back(ComplementarityError(candidate.gap, weighted_force));
            diag->min_gap = std::min(diag->min_gap, candidate.gap);
            diag->max_penetration = std::max(diag->max_penetration, std::max(0.0, -candidate.gap));
            diag->max_lambda_n = std::max(diag->max_lambda_n, weighted_force);
            diag->max_complementarity_error =
                std::max(diag->max_complementarity_error, diag->complementarity_errors.back());
        }
    }
    return generalized_torque;
}

GearStepResult SolveSimpleGearRecurDynGGeomContactLawStep(const ChOpenVdbSdfGrid& gear21_sdf,
                                                          const ChOpenVdbSdfGrid& gear22_sdf,
                                                          const ChSdfContactSurfaceGraph& gear21_graph,
                                                          const ChSdfContactSurfaceGraph& gear22_graph,
                                                          const GearPose& pose21,
                                                          const GearPose& pose22,
                                                          const SimpleGearConfig& config,
                                                          const SimpleGearGGeomContactLawSI& law,
                                                          double inertia22,
                                                          const GearState& state,
                                                          GearPatchMemory* memory = nullptr,
                                                          const ChSdfContactSampleBvh* gear21_bvh = nullptr,
                                                          const ChSdfContactSampleBvh* gear22_bvh = nullptr) {
    GearStepResult result;
    result.state = state;
    const double theta21_free = state.theta21 + config.dt * config.omega21;
    const double theta22_free = state.theta22 + config.dt * state.omega22;
    int patch_count = 0;
    const std::vector<GearContactSampleRef> sample_refs = BuildGearBidirectionalPatchSamples(gear21_sdf,
                                                                                             gear22_sdf,
                                                                                             gear21_graph,
                                                                                             gear22_graph,
                                                                                             pose21,
                                                                                             pose22,
                                                                                             config,
                                                                                             theta21_free,
                                                                                             theta22_free,
                                                                                             &patch_count,
                                                                                             memory,
                                                                                             gear21_bvh,
                                                                                             gear22_bvh);

    auto residual = [&](double omega_next) {
        const double theta22_next = state.theta22 + config.dt * omega_next;
        const double torque = EvaluateSimpleGearGGeomContactLaw(gear21_sdf,
                                                                gear22_sdf,
                                                                gear21_graph,
                                                                gear22_graph,
                                                                pose21,
                                                                pose22,
                                                                config,
                                                                law,
                                                                sample_refs,
                                                                theta21_free,
                                                                theta22_next,
                                                                omega_next,
                                                                nullptr);
        return inertia22 * (omega_next - state.omega22) - config.dt * torque;
    };

    double omega = state.omega22;
    double residual_norm = std::abs(residual(omega));
    int iterations = 0;
    bool success = residual_norm <= config.tolerance;
    for (iterations = 0; !success && iterations < config.max_iterations; iterations++) {
        const double h = std::max(1.0e-7, 1.0e-6 * std::max(1.0, std::abs(omega)));
        const double derivative = (residual(omega + h) - residual(omega - h)) / (2.0 * h);
        if (!std::isfinite(derivative) || std::abs(derivative) < 1.0e-18) {
            break;
        }
        const double delta = -residual(omega) / derivative;
        double alpha = 1.0;
        double accepted_omega = omega;
        double accepted_norm = residual_norm;
        while (alpha >= 1.0e-5) {
            const double trial = omega + alpha * delta;
            const double trial_norm = std::abs(residual(trial));
            if (std::isfinite(trial_norm) && trial_norm <= accepted_norm) {
                accepted_omega = trial;
                accepted_norm = trial_norm;
                break;
            }
            alpha *= 0.5;
        }
        omega = accepted_omega;
        residual_norm = accepted_norm;
        success = residual_norm <= config.tolerance;
        if (std::abs(alpha * delta) <= 1.0e-12 * std::max(1.0, std::abs(omega))) {
            success = residual_norm <= 1.0e-6;
            break;
        }
    }

    result.state.omega22 = omega;
    result.state.theta21 = theta21_free;
    result.state.theta22 = state.theta22 + config.dt * result.state.omega22;
    EvaluateSimpleGearGGeomContactLaw(gear21_sdf,
                                      gear22_sdf,
                                      gear21_graph,
                                      gear22_graph,
                                      pose21,
                                      pose22,
                                      config,
                                      law,
                                      sample_refs,
                                      result.state.theta21,
                                      result.state.theta22,
                                      result.state.omega22,
                                      &result.diagnostics);
    result.diagnostics.patch_count = patch_count;
    result.diagnostics.success = success || residual_norm <= 1.0e-6;
    result.diagnostics.iterations = iterations;
    result.diagnostics.residual_norm = residual_norm;
    result.diagnostics.min_gap = MinGearGapBidirectional(gear21_sdf,
                                                         gear22_sdf,
                                                         gear21_graph,
                                                         gear22_graph,
                                                         pose21,
                                                         pose22,
                                                         result.state.theta21,
                                                         result.state.theta22);
    result.diagnostics.max_penetration =
        std::max(result.diagnostics.max_penetration, std::max(0.0, -result.diagnostics.min_gap));
    return result;
}

int RunSimpleGearFullGeometryCase(const std::string& result_case_name = "simple_gear",
                                  double t_end_override = -1.0,
                                  double dt_override = -1.0,
                                  SimpleGearBackendMode backend_mode = SimpleGearBackendMode::SdfNcp) {
    const auto start = std::chrono::steady_clock::now();
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "simple_gear";
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / result_case_name;
    std::filesystem::create_directories(out_dir);

    openvdb::initialize();
    SimpleGearConfig config;
    config.case_name = result_case_name;
    if (t_end_override > 0.0) {
        config.t_end = t_end_override;
    }
    if (dt_override > 0.0) {
        config.dt = dt_override;
    }
    bool impulse_mixed_env_set = false;
    if (const char* dt_env = std::getenv("SDF_NCP_SIMPLE_GEAR_DT")) {
        try {
            config.dt = std::stod(dt_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_DT value.");
        }
    }
    if (const char* offset_scale_env = std::getenv("SDF_NCP_SIMPLE_GEAR_OFFSET_SCALE")) {
        try {
            config.sdf_ncp_contact_offset_scale = std::stod(offset_scale_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_OFFSET_SCALE value.");
        }
    }
    if (const char* activation_band_env = std::getenv("SDF_NCP_SIMPLE_GEAR_ACTIVATION_BAND")) {
        try {
            config.activation_band = std::stod(activation_band_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_ACTIVATION_BAND value.");
        }
    }
    if (const char* max_activation_band_env = std::getenv("SDF_NCP_SIMPLE_GEAR_MAX_ACTIVATION_BAND")) {
        try {
            config.max_activation_band = std::stod(max_activation_band_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_MAX_ACTIVATION_BAND value.");
        }
    }
    if (const char* min_patch_samples_env = std::getenv("SDF_NCP_SIMPLE_GEAR_MIN_PATCH_SAMPLES")) {
        try {
            config.min_patch_samples = std::stoi(min_patch_samples_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_MIN_PATCH_SAMPLES value.");
        }
    }
    if (const char* pressure_scale_env = std::getenv("SDF_NCP_SIMPLE_GEAR_PRESSURE_SCALE")) {
        try {
            config.sdf_ncp_pressure_scale_override = std::stod(pressure_scale_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_PRESSURE_SCALE value.");
        }
    }
    if (const char* pressure_scale_factor_env = std::getenv("SDF_NCP_SIMPLE_GEAR_PRESSURE_SCALE_BAND_FACTOR")) {
        try {
            config.sdf_ncp_pressure_scale_band_factor = std::stod(pressure_scale_factor_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_PRESSURE_SCALE_BAND_FACTOR value.");
        }
    }
    if (const char* gap_scale_env = std::getenv("SDF_NCP_SIMPLE_GEAR_GAP_SCALE")) {
        try {
            config.sdf_ncp_gap_scale = std::stod(gap_scale_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_GAP_SCALE value.");
        }
    }
    if (const char* hysteresis_env = std::getenv("SDF_NCP_SIMPLE_GEAR_ACTIVE_BAND_HYSTERESIS")) {
        try {
            config.active_band_hysteresis = std::stod(hysteresis_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_ACTIVE_BAND_HYSTERESIS value.");
        }
    }
    if (const char* mixed_env = std::getenv("SDF_NCP_SIMPLE_GEAR_IMPULSE_MIXED")) {
        try {
            impulse_mixed_env_set = true;
            config.use_impulse_mixed_ncp = std::stoi(mixed_env) != 0;
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_IMPULSE_MIXED value.");
        }
    }
    if (const char* beta_env = std::getenv("SDF_NCP_SIMPLE_GEAR_IMPULSE_BETA")) {
        try {
            config.impulse_mixed_beta = std::stod(beta_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_IMPULSE_BETA value.");
        }
    }
    if (const char* cfm_env = std::getenv("SDF_NCP_SIMPLE_GEAR_IMPULSE_CFM")) {
        try {
            config.impulse_mixed_cfm = std::stod(cfm_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_IMPULSE_CFM value.");
        }
    }
    if (const char* velocity_scale_env = std::getenv("SDF_NCP_SIMPLE_GEAR_VELOCITY_SCALE")) {
        try {
            config.impulse_mixed_velocity_scale = std::stod(velocity_scale_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_VELOCITY_SCALE value.");
        }
    }
    if (const char* impulse_scale_env = std::getenv("SDF_NCP_SIMPLE_GEAR_IMPULSE_SCALE")) {
        try {
            config.impulse_mixed_impulse_scale = std::stod(impulse_scale_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_IMPULSE_SCALE value.");
        }
    }
    if (const char* semismooth_env = std::getenv("SDF_NCP_SIMPLE_GEAR_SEMISMOOTH_ACTIVE_SET")) {
        try {
            config.use_semismooth_active_set = std::stoi(semismooth_env) != 0;
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_SEMISMOOTH_ACTIVE_SET value.");
        }
    }
    if (const char* active_iter_env = std::getenv("SDF_NCP_SIMPLE_GEAR_MAX_ACTIVE_SET_ITERATIONS")) {
        try {
            config.max_active_set_iterations = std::stoi(active_iter_env);
        } catch (...) {
            throw std::runtime_error("Invalid SDF_NCP_SIMPLE_GEAR_MAX_ACTIVE_SET_ITERATIONS value.");
        }
    }
    if (!impulse_mixed_env_set && config.dt < 1.0e-3) {
        config.use_impulse_mixed_ncp = true;
    }
    const SimpleGearRmd rmd = LoadSimpleGearRmdForNcp(asset_dir / "simple gear.rmd");
    config.omega21 = rmd.rev_joint1_motion.has_function ?
                         rmd.rev_joint1_motion.function_constant :
                         LoadFirstMotionFunctionConstant(asset_dir / "simple gear.rmd", config.omega21);
    const SimpleGearGGeomContactLawSI ggeom_law = MakeSimpleGearGGeomContactLawSI(rmd.geo_surface_contact);
    if (backend_mode == SimpleGearBackendMode::SdfNcp && config.sdf_ncp_pressure_scale_override <= 0.0) {
        const double automatic_scale =
            ggeom_law.stiffness * config.max_activation_band * config.sdf_ncp_pressure_scale_band_factor;
        if (automatic_scale > 0.0 && std::isfinite(automatic_scale)) {
            config.sdf_ncp_pressure_scale_override = automatic_scale;
        }
    }
    if (backend_mode == SimpleGearBackendMode::SdfNcp) {
        config.sdf_ncp_contact_offset = config.sdf_ncp_contact_offset_scale * ggeom_law.boundary_penetration;
    }
    WriteSimpleGearRmdMappingAudit(out_dir / "rmd_mapping.csv", rmd);
    GearPose pose21{rmd.gear21.cm_marker_mm * kMmToM, rmd.gear21.part_rotation};
    GearPose pose22{rmd.gear22.cm_marker_mm * kMmToM, rmd.gear22.part_rotation};
    ChSdfTriangleMeshData gear21 = LoadGearObjAsBodyLocal(asset_dir / "gear_21.obj",
                                                          rmd.gear21.surface_ref_marker_mm,
                                                          rmd.gear21.surface_ref_rotation,
                                                          rmd.gear21.cm_marker_mm,
                                                          rmd.gear21.part_rotation);
    ChSdfTriangleMeshData gear22 = LoadGearObjAsBodyLocal(asset_dir / "gear_22.obj",
                                                          rmd.gear22.surface_ref_marker_mm,
                                                          rmd.gear22.surface_ref_rotation,
                                                          rmd.gear22.cm_marker_mm,
                                                          rmd.gear22.part_rotation);
    ChOpenVdbSdfGrid gear21_sdf = BuildOpenVdbLevelSetFromMesh(
        gear21, config.voxel_size, config.half_width_voxels);
    ChOpenVdbSdfGrid gear22_sdf = BuildOpenVdbLevelSetFromMesh(
        gear22, config.voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph gear21_graph = MakeSurfaceGraphFromMeshForSdf(gear21);
    ChSdfContactSurfaceGraph gear22_graph = MakeSurfaceGraphFromMeshForSdf(gear22);
    ChSdfContactSampleBvh gear21_bvh(gear21_graph);
    ChSdfContactSampleBvh gear22_bvh(gear22_graph);
    const auto reference = LoadGear22Reference(asset_dir / "data" / "Gear22.csv");

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    std::ofstream reference_comparison(out_dir / "gear22_reference_comparison_timeseries.csv");
    std::ofstream analytic_comparison(out_dir / "gear22_analytic_comparison_timeseries.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    reference_comparison << std::setprecision(17);
    analytic_comparison << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";
    reference_comparison
        << "time,gear22_omega_rx_reference,gear22_omega_rx_sdf_ncp,omega_error,omega_abs_error,"
           "gear22_alpha_rx_reference,gear22_alpha_rx_sdf_ncp,alpha_error,alpha_abs_error\n";
    analytic_comparison
        << "time,gear22_omega_rx_analytic,gear22_omega_rx_computed,omega_error,omega_abs_error\n";

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = config.case_name;
    stats.method =
        backend_mode == SimpleGearBackendMode::RecurDynGGeomContactLaw ? "recurdyn_ggeomcontact_law" : "sdf_ncp";
    stats.asset = config.asset;
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    stats.num_contacts = 0;
    stats.notes = backend_mode == SimpleGearBackendMode::RecurDynGGeomContactLaw ?
                      "full GEAR21/GEAR22 OBJ/OpenVDB SDF; bidirectional adaptive active-band tooth patch samples; "
                      "RecurDyn RevJoint1.RMotion driver; RecurDyn GGEOMCONTACT K/C/KORDER/BPEN force law converted "
                      "from mm/kg/s RMD units to SI for frontend calibration" :
                      "full GEAR21/GEAR22 OBJ/OpenVDB SDF; bidirectional adaptive active-band tooth patch SDF-NCP "
                      "samples; connected-component patch partition; sample-level persistent warm start; "
                      "area-integrated local contact intensity SDF-NCP assembly; RecurDyn RevJoint1.RMotion driver "
                      "and GGEOMCONTACT metadata parsed from simple gear.rmd; AABB BVH broad phase filters source "
                      "surface samples before exact OpenVDB query; hard frictionless NCP does not apply RMD K/C/KORDER";
    if (backend_mode == SimpleGearBackendMode::SdfNcp && config.use_impulse_mixed_ncp) {
        stats.notes +=
            "; dt-adaptive impulse/velocity-level mixed NCP enabled with AD Jacobian and dt-scaled impulse weights";
    }

    GearState state;
    GearPatchMemory patch_memory;
    GearReferenceComparisonStats reference_stats;
    GearAnalyticComparisonStats analytic_stats;
    analytic_stats.target_omega = -config.omega21;
    double previous_omega22 = state.omega22;
    for (int i = 0; i <= steps; i++) {
        const double time = static_cast<double>(i) * config.dt;
        GearStepResult step;
        if (i == 0) {
            step.state = state;
            int patch_count = 0;
            const std::vector<GearContactSampleRef> refs = BuildGearBidirectionalPatchSamples(gear21_sdf,
                                                                                              gear22_sdf,
                                                                                              gear21_graph,
                                                                                              gear22_graph,
                                                                                              pose21,
                                                                                              pose22,
                                                                                              config,
                                                                                              state.theta21,
                                                                                              state.theta22,
                                                                                              &patch_count,
                                                                                              &patch_memory,
                                                                                              &gear21_bvh,
                                                                                              &gear22_bvh);
            if (backend_mode == SimpleGearBackendMode::RecurDynGGeomContactLaw) {
                EvaluateSimpleGearGGeomContactLaw(gear21_sdf,
                                                  gear22_sdf,
                                                  gear21_graph,
                                                  gear22_graph,
                                                  pose21,
                                                  pose22,
                                                  config,
                                                  ggeom_law,
                                                  refs,
                                                  state.theta21,
                                                  state.theta22,
                                                  state.omega22,
                                                  &step.diagnostics);
            } else {
                std::vector<double> z(1 + refs.size(), 0.0);
                z[0] = state.omega22;
                ComputeSimpleGearResidual(gear21_sdf,
                                          gear22_sdf,
                                          gear21_graph,
                                          gear22_graph,
                                          pose21,
                                          pose22,
                                          config,
                                          rmd.gear22.inertia_x_kg_m2,
                                          state,
                                          refs,
                                          z,
                                          &step.diagnostics);
            }
            step.diagnostics.patch_count = patch_count;
            step.diagnostics.success = true;
            step.diagnostics.min_gap = MinGearGapBidirectional(gear21_sdf,
                                                               gear22_sdf,
                                                               gear21_graph,
                                                               gear22_graph,
                                                               pose21,
                                                               pose22,
                                                               state.theta21,
                                                               state.theta22);
            step.diagnostics.max_penetration =
                std::max(step.diagnostics.max_penetration, std::max(0.0, -step.diagnostics.min_gap));
        } else {
            if (backend_mode == SimpleGearBackendMode::RecurDynGGeomContactLaw) {
                step = SolveSimpleGearRecurDynGGeomContactLawStep(gear21_sdf,
                                                                  gear22_sdf,
                                                                  gear21_graph,
                                                                  gear22_graph,
                                                                  pose21,
                                                                  pose22,
                                                                  config,
                                                                  ggeom_law,
                                                                  rmd.gear22.inertia_x_kg_m2,
                                                                  state,
                                                                  &patch_memory,
                                                                  &gear21_bvh,
                                                                  &gear22_bvh);
            } else {
                step = SolveSimpleGearStep(gear21_sdf,
                                           gear22_sdf,
                                           gear21_graph,
                                           gear22_graph,
                                           pose21,
                                           pose22,
                                           config,
                                           rmd.gear22.inertia_x_kg_m2,
                                           state,
                                           &patch_memory,
                                           &gear21_bvh,
                                           &gear22_bvh);
            }
            state = step.state;
        }

        if (step.diagnostics.candidates.empty()) {
            const double q0 = std::cos(0.5 * state.theta22);
            const double q1 = std::sin(0.5 * state.theta22);
            trajectory << time << ",gear22," << pose22.center.x() << "," << pose22.center.y() << ","
                       << pose22.center.z() << "," << q0 << "," << q1 << ",0,0,0,0,0,"
                       << state.omega22 << ",0,0,-1," << step.diagnostics.min_gap << ",0,"
                       << "0,0," << std::max(0.0, -step.diagnostics.min_gap) << ",0,0,"
                       << step.diagnostics.residual_norm << "," << step.diagnostics.iterations << ","
                       << (step.diagnostics.success ? 1 : 0) << "\n";
        }
        for (size_t c = 0; c < step.diagnostics.candidates.size(); c++) {
            const auto& candidate = step.diagnostics.candidates[c];
            const double lambda = c < step.diagnostics.lambdas.size() ? step.diagnostics.lambdas[c] : 0.0;
            const double lambda_force =
                c < step.diagnostics.lambda_forces.size() ? step.diagnostics.lambda_forces[c] : lambda * candidate.weight;
            const double ncp = c < step.diagnostics.ncp_residuals.size() ? step.diagnostics.ncp_residuals[c] : 0.0;
            const double comp = c < step.diagnostics.complementarity_errors.size() ?
                                    step.diagnostics.complementarity_errors[c] :
                                    0.0;
            const double q0 = std::cos(0.5 * state.theta22);
            const double q1 = std::sin(0.5 * state.theta22);
            trajectory << time << ",gear22," << pose22.center.x() << "," << pose22.center.y() << ","
                       << pose22.center.z() << "," << q0 << "," << q1 << ",0,0,0,0,0,"
                       << state.omega22 << ",0,0," << candidate.sample_id << "," << candidate.gap << ","
                       << lambda << "," << candidate.weight << "," << lambda_force << ","
                       << std::max(0.0, -candidate.gap) << "," << ncp << "," << comp << ","
                       << step.diagnostics.residual_norm << "," << step.diagnostics.iterations << ","
                       << (step.diagnostics.success ? 1 : 0) << "\n";
        }
        if (!reference.empty()) {
            const GearReferenceRow ref_row = InterpolateGearReference(reference, time);
            const double alpha22 = i == 0 ? 0.0 : (state.omega22 - previous_omega22) / config.dt;
            const double omega_error = state.omega22 - ref_row.omega_rx;
            const double alpha_error = alpha22 - ref_row.alpha_rx;
            reference_comparison << time << "," << ref_row.omega_rx << "," << state.omega22 << ","
                                 << omega_error << "," << std::abs(omega_error) << "," << ref_row.alpha_rx << ","
                                 << alpha22 << "," << alpha_error << "," << std::abs(alpha_error) << "\n";
            reference_stats.samples++;
            reference_stats.sum_sq_omega += omega_error * omega_error;
            reference_stats.sum_sq_alpha += alpha_error * alpha_error;
            reference_stats.max_abs_omega = std::max(reference_stats.max_abs_omega, std::abs(omega_error));
            reference_stats.max_abs_alpha = std::max(reference_stats.max_abs_alpha, std::abs(alpha_error));
            reference_stats.max_abs_ref_omega =
                std::max(reference_stats.max_abs_ref_omega, std::abs(ref_row.omega_rx));
            reference_stats.final_ref = ref_row;
            reference_stats.final_omega = state.omega22;
            reference_stats.final_alpha = alpha22;
        }
        const double analytic_omega = time <= 0.0 ? 0.0 : analytic_stats.target_omega;
        const double analytic_error = state.omega22 - analytic_omega;
        analytic_comparison << time << "," << analytic_omega << "," << state.omega22 << "," << analytic_error << ","
                            << std::abs(analytic_error) << "\n";
        analytic_stats.samples++;
        analytic_stats.sum_sq_omega += analytic_error * analytic_error;
        analytic_stats.sum_abs_omega += std::abs(analytic_error);
        analytic_stats.max_abs_omega = std::max(analytic_stats.max_abs_omega, std::abs(analytic_error));
        analytic_stats.final_omega = state.omega22;
        analytic_stats.final_abs_error = std::abs(analytic_error);
        previous_omega22 = state.omega22;
        Accumulate(stats, step.diagnostics);
    }

    const auto ref = FindNearestGearReference(reference, config.t_end);
    const double abs_omega = std::abs(state.omega22 - ref.omega_rx);
    const double rel_omega = std::abs(ref.omega_rx) > 1.0e-14 ? abs_omega / std::abs(ref.omega_rx) : 0.0;
    comparison << config.case_name << ",gear22_omega_rx_nearest_reference," << ref.omega_rx << "," << state.omega22
               << "," << abs_omega << "," << rel_omega
               << ",reference from assets/simple_gear/data/Gear22.csv; SDF-NCP uses frictionless full mesh query\n";
    comparison << config.case_name << ",gear21_rmotion_velocity_constant," << config.omega21 << ","
               << config.omega21 << ",0,0,parsed from assets/simple_gear/simple gear.rmd FUNCTION field\n";
    comparison << config.case_name << ",ggeomcontact_K," << rmd.geo_surface_contact.stiffness << ","
               << rmd.geo_surface_contact.stiffness
               << ",0,0,parsed from RecurDyn GGEOMCONTACT; reported for contact-law alignment audit\n";
    comparison << config.case_name << ",ggeomcontact_C," << rmd.geo_surface_contact.damping << ","
               << rmd.geo_surface_contact.damping
               << ",0,0,parsed from RecurDyn GGEOMCONTACT; reported for contact-law alignment audit\n";
    comparison << config.case_name << ",ggeomcontact_bpen_m," << rmd.geo_surface_contact.bpen * 1.0e-3 << ","
               << ggeom_law.boundary_penetration
               << ",0,0,RecurDyn BPEN converted from millimeter RMD length units to meters\n";
    comparison << config.case_name << ",ggeomcontact_K_SI," << ggeom_law.stiffness << "," << ggeom_law.stiffness
               << ",0,0,RecurDyn K converted for meter-based gap with KORDER exponent\n";
    comparison << config.case_name << ",ggeomcontact_C_SI," << ggeom_law.damping << "," << ggeom_law.damping
               << ",0,0,RecurDyn C converted from N/(mm/s) to N/(m/s)\n";
    comparison << config.case_name << ",sdf_ncp_contact_offset_m," << ggeom_law.boundary_penetration << ","
               << config.sdf_ncp_contact_offset << ",0,0,"
               << "hard SDF-NCP shell offset used only in SDF-NCP mode; current tuned default keeps it disabled\n";
    comparison << config.case_name << ",sdf_ncp_pressure_scale,0,"
               << config.sdf_ncp_pressure_scale_override << "," << config.sdf_ncp_pressure_scale_override
               << ",0,automatic local-intensity-to-pressure scale; default uses GGEOMCONTACT K_SI * active band * "
                  "band factor, env SDF_NCP_SIMPLE_GEAR_PRESSURE_SCALE can override\n";
    comparison << config.case_name << ",sdf_ncp_pressure_scale_band_factor,0,"
               << config.sdf_ncp_pressure_scale_band_factor << "," << config.sdf_ncp_pressure_scale_band_factor
               << ",0,dimensionless factor for automatic local-intensity pressure scale\n";
    comparison << config.case_name << ",sdf_ncp_gap_scale_m,0,"
               << (config.sdf_ncp_gap_scale > 0.0 ? config.sdf_ncp_gap_scale : 1.0)
               << ",0,0,scaled Fischer-Burmeister gap scale used by hard SDF-NCP local-intensity mode\n";
    comparison << config.case_name << ",sdf_ncp_active_band_hysteresis_m,0," << config.active_band_hysteresis
               << "," << config.active_band_hysteresis
               << ",0,hysteresis used to retain persistent active patch samples between steps\n";
    comparison << config.case_name << ",sdf_ncp_impulse_mixed_enabled,0,"
               << (config.use_impulse_mixed_ncp ? 1 : 0) << "," << (config.use_impulse_mixed_ncp ? 1 : 0)
               << ",0,dt-adaptive velocity-level impulse NCP backend; can be overridden by "
                  "SDF_NCP_SIMPLE_GEAR_IMPULSE_MIXED\n";
    comparison << config.case_name << ",sdf_ncp_impulse_mixed_beta,0," << config.impulse_mixed_beta << ","
               << config.impulse_mixed_beta
               << ",0,Baumgarte gap stabilization used by impulse mixed NCP when enabled\n";
    comparison << config.case_name << ",sdf_ncp_impulse_mixed_cfm,0," << config.impulse_mixed_cfm << ","
               << config.impulse_mixed_cfm
               << ",0,velocity-level contact force mixing regularization used by impulse mixed NCP when enabled\n";
    comparison << config.case_name << ",sdf_ncp_activation_band_m,0," << config.activation_band << ","
               << config.activation_band << ",0,initial active patch band around the minimum signed distance\n";
    comparison << config.case_name << ",sdf_ncp_max_activation_band_m,0," << config.max_activation_band << ","
               << config.max_activation_band << ",0,maximum active patch band used to satisfy min patch samples\n";
    comparison << config.case_name << ",sdf_ncp_min_patch_samples,0," << config.min_patch_samples << ","
               << config.min_patch_samples << ",0,minimum active samples per connected direction before band cap\n";
    comparison << config.case_name << ",sdf_ncp_openvdb_geometry_derivatives_enabled,0,"
               << (config.use_impulse_mixed_ncp ? 1 : 0) << "," << (config.use_impulse_mixed_ncp ? 1 : 0)
               << ",0,uses OpenVDB-interpolated grad/Hessian chain-rule dJ/dtheta only in the impulse mixed backend\n";
    comparison << config.case_name << ",sdf_ncp_semismooth_active_set_enabled,0,"
               << (config.use_semismooth_active_set ? 1 : 0) << ","
               << (config.use_semismooth_active_set ? 1 : 0)
               << ",0,rebuilds active patch samples after a trial solve until sample ids stabilize or iteration limit "
                  "is reached\n";
    comparison << config.case_name << ",sdf_ncp_max_active_set_iterations,0," << config.max_active_set_iterations
               << "," << config.max_active_set_iterations << ",0,semismooth active patch iteration cap\n";
    comparison << config.case_name << ",sdf_ncp_aabb_bvh_broadphase_enabled,0,1,1,0,"
               << "backend-neutral AABB BVH filters source-surface samples before full OpenVDB SDF query\n";
    comparison << config.case_name << ",sdf_ncp_aabb_bvh_leaf_size,0,16,16,0,"
               << "default ChSdfContactSampleBvh leaf size used for benchmark sample graphs\n";
    if (analytic_stats.samples > 0) {
        const double denom = static_cast<double>(analytic_stats.samples);
        const double analytic_rmse = std::sqrt(analytic_stats.sum_sq_omega / denom);
        const double analytic_mae = analytic_stats.sum_abs_omega / denom;
        const double analytic_rel = std::abs(analytic_stats.target_omega) > 1.0e-14 ?
                                        analytic_rmse / std::abs(analytic_stats.target_omega) :
                                        0.0;
        comparison << config.case_name << ",gear22_omega_rx_analytic_target," << analytic_stats.target_omega << ","
                   << analytic_stats.final_omega << "," << analytic_stats.final_abs_error << ","
                   << (std::abs(analytic_stats.target_omega) > 1.0e-14 ?
                           analytic_stats.final_abs_error / std::abs(analytic_stats.target_omega) :
                           0.0)
                   << ",1:1 simple gear analytic target from omega21=+1 rad/s -> omega22=-1 rad/s\n";
        comparison << config.case_name << ",gear22_omega_rx_analytic_rmse,0," << analytic_rmse << ","
                   << analytic_rmse << "," << analytic_rel
                   << ",RMSE against analytic omega22=-omega21 over simulated one-second window\n";
        comparison << config.case_name << ",gear22_omega_rx_analytic_mae,0," << analytic_mae << ","
                   << analytic_mae << ","
                   << (std::abs(analytic_stats.target_omega) > 1.0e-14 ?
                           analytic_mae / std::abs(analytic_stats.target_omega) :
                           0.0)
                   << ",MAE against analytic omega22=-omega21 over simulated one-second window\n";

        std::ofstream analytic_summary(out_dir / "gear22_analytic_comparison_summary.csv");
        analytic_summary << std::setprecision(17);
        analytic_summary
            << "case_name,quantity,time_start,time_end,num_samples,target_value,rmse,mae,max_abs_error,final_value,"
               "final_abs_error,relative_rmse_target_scale,notes\n";
        analytic_summary << config.case_name << ",gear22_omega_rx,0," << config.t_end << ","
                         << analytic_stats.samples << "," << analytic_stats.target_omega << "," << analytic_rmse
                         << "," << analytic_mae << "," << analytic_stats.max_abs_omega << ","
                         << analytic_stats.final_omega << "," << analytic_stats.final_abs_error << "," << analytic_rel
                         << ",Analytic 1:1 gear target; time=0 target is written as 0 in the time series to match "
                            "the stored initial state\n";
    }
    if (reference_stats.samples > 0) {
        const double denom = static_cast<double>(reference_stats.samples);
        const double omega_rmse = std::sqrt(reference_stats.sum_sq_omega / denom);
        const double alpha_rmse = std::sqrt(reference_stats.sum_sq_alpha / denom);
        const double omega_rel = reference_stats.max_abs_ref_omega > 1.0e-14 ?
                                     omega_rmse / reference_stats.max_abs_ref_omega :
                                     0.0;
        comparison << config.case_name << ",gear22_omega_rx_rmse,0," << omega_rmse << "," << omega_rmse << ","
                   << omega_rel
                   << ",RMSE against interpolated assets/simple_gear/data/Gear22.csv over simulated time window\n";
        comparison << config.case_name << ",gear22_omega_rx_max_abs_error,0," << reference_stats.max_abs_omega
                   << "," << reference_stats.max_abs_omega << ",0,"
                   << "maximum absolute omega RX error against interpolated RecurDyn output\n";

        std::ofstream reference_summary(out_dir / "gear22_reference_comparison_summary.csv");
        reference_summary << std::setprecision(17);
        reference_summary
            << "case_name,quantity,time_start,time_end,num_samples,rmse,max_abs_error,final_recurdyn,final_sdf_ncp,"
               "final_abs_error,relative_rmse_reference_scale,notes\n";
        reference_summary << config.case_name << ",gear22_omega_rx,0," << config.t_end << ","
                          << reference_stats.samples << "," << omega_rmse << "," << reference_stats.max_abs_omega
                          << "," << reference_stats.final_ref.omega_rx << "," << reference_stats.final_omega << ","
                          << std::abs(reference_stats.final_omega - reference_stats.final_ref.omega_rx) << ","
                          << omega_rel
                          << ",RecurDyn Y:Vel_RX-GEAR22-jiandanjiaolian compared with SDF-NCP gear22 generalized RX "
                             "velocity\n";
        reference_summary << config.case_name << ",gear22_alpha_rx,0," << config.t_end << ","
                          << reference_stats.samples << "," << alpha_rmse << "," << reference_stats.max_abs_alpha
                          << "," << reference_stats.final_ref.alpha_rx << "," << reference_stats.final_alpha << ","
                          << std::abs(reference_stats.final_alpha - reference_stats.final_ref.alpha_rx) << ",0,"
                          << "SDF-NCP acceleration is finite-differenced from omega22; RecurDyn column is "
                             "Y:Acc_RX-GEAR22-jiandanjiaolian\n";
    }

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);
    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    const double success_rate = stats.sum_success / static_cast<double>(std::max(1, stats.samples));
    return stats.max_penetration < 5.0e-3 && stats.max_complementarity_error < 5.0e-2 && success_rate > 0.80 ? 0 : 1;
}

std::vector<ChSdfNcpDescriptorContact> BuildCamDescriptorContacts(
    const std::shared_ptr<ChBodyAuxRef>& cam,
    const std::shared_ptr<ChBodyAuxRef>& follower,
    const ChOpenVdbSdfGrid& cam_sdf,
    const ChSdfContactSurfaceGraph& follower_graph,
    const CamLikeConfig& config,
    const ChSdfContactSampleBvh* follower_bvh = nullptr,
    const ChSdfAabb* cam_local_bounds = nullptr,
    double broad_phase_padding = 0.0) {
    struct SamplePhi {
        int sample_id = -1;
        double phi = 0.0;
        double area = 0.0;
        ChVector3d world_pos = ChVector3d(0, 0, 0);
        ChVector3d normal_world = ChVector3d(0, 1, 0);
    };

    std::vector<SamplePhi> sampled;
    sampled.reserve(follower_graph.samples.size());
    double min_phi = std::numeric_limits<double>::max();

    const auto& follower_ref = follower->GetFrameRefToAbs();
    const auto& cam_ref = cam->GetFrameRefToAbs();
    auto all_sample_ids = [&]() {
        std::vector<int> ids;
        ids.reserve(follower_graph.samples.size());
        for (const auto& sample : follower_graph.samples) {
            if (sample.area > 0.0) {
                ids.push_back(sample.id);
            }
        }
        return ids;
    };

    std::vector<int> candidate_ids;
    if (follower_bvh && cam_local_bounds && !follower_bvh->Empty() && cam_local_bounds->IsValid()) {
        ChSdfAabb padded = ExpandedAabb(*cam_local_bounds, broad_phase_padding);
        ChSdfAabb query_bounds;
        for (const auto& cam_corner : AabbCorners(padded)) {
            const ChVector3d world = cam_ref.TransformPointLocalToParent(cam_corner);
            query_bounds.Include(follower_ref.TransformPointParentToLocal(world));
        }
        candidate_ids = follower_bvh->Query(query_bounds);
        std::sort(candidate_ids.begin(), candidate_ids.end());
        candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());
    } else {
        candidate_ids = all_sample_ids();
    }

    if (candidate_ids.empty()) {
        return {};
    }

    for (int sample_id : candidate_ids) {
        if (sample_id < 0 || sample_id >= static_cast<int>(follower_graph.samples.size())) {
            continue;
        }
        const auto& sample = follower_graph.samples[static_cast<size_t>(sample_id)];
        if (sample.area <= 0.0) {
            continue;
        }
        const ChVector3d world_pos = follower_ref.TransformPointLocalToParent(sample.local_pos);
        const ChVector3d cam_local = cam_ref.TransformPointParentToLocal(world_pos);
        const auto query = cam_sdf.QueryLocal(cam_local);
        SamplePhi item;
        item.sample_id = sample.id;
        item.phi = query.phi;
        item.area = sample.area;
        item.world_pos = world_pos;
        item.normal_world = cam_ref.TransformDirectionLocalToParent(query.grad).GetNormalized();
        sampled.push_back(item);
        min_phi = std::min(min_phi, query.phi);
    }

    if (sampled.empty()) {
        return {};
    }
    if (min_phi > config.max_activation_band) {
        return {};
    }

    std::sort(sampled.begin(), sampled.end(), [](const SamplePhi& a, const SamplePhi& b) {
        return a.phi < b.phi;
    });

    double active_band = config.activation_band;
    std::vector<SamplePhi> active;
    while (active.empty() || (static_cast<int>(active.size()) < config.min_patch_samples &&
                              active_band < config.max_activation_band)) {
        active.clear();
        const double threshold = std::max(active_band, min_phi + active_band);
        for (const auto& item : sampled) {
            if (item.phi <= threshold) {
                active.push_back(item);
            }
        }
        if (static_cast<int>(active.size()) >= config.min_patch_samples) {
            break;
        }
        active_band = std::min(config.max_activation_band, active_band * 2.0);
        if (active_band >= config.max_activation_band) {
            break;
        }
    }

    if (active.empty()) {
        active.push_back(sampled.front());
    }
    if (static_cast<int>(active.size()) < config.min_patch_samples) {
        for (const auto& item : sampled) {
            const auto duplicate = std::find_if(active.begin(), active.end(), [&item](const SamplePhi& existing) {
                return existing.sample_id == item.sample_id;
            });
            if (duplicate == active.end()) {
                active.push_back(item);
                if (static_cast<int>(active.size()) >= config.min_patch_samples) {
                    break;
                }
            }
        }
    }

    if (follower_bvh && static_cast<int>(active.size()) < std::max(1, config.min_patch_samples) &&
        candidate_ids.size() < follower_graph.samples.size()) {
        return BuildCamDescriptorContacts(cam, follower, cam_sdf, follower_graph, config, nullptr, nullptr, 0.0);
    }

    double total_area = 0.0;
    for (const auto& item : active) {
        total_area += std::max(0.0, item.area);
    }
    std::vector<int> active_ids;
    active_ids.reserve(active.size());
    for (const auto& item : active) {
        active_ids.push_back(item.sample_id);
    }
    const auto components = follower_graph.FindConnectedComponents(active_ids);
    std::vector<int> patch_ids(follower_graph.samples.size(), -1);
    std::vector<double> patch_areas(components.size(), 0.0);
    for (size_t patch = 0; patch < components.size(); patch++) {
        for (int sample_id : components[patch]) {
            if (sample_id >= 0 && sample_id < static_cast<int>(follower_graph.samples.size())) {
                patch_ids[static_cast<size_t>(sample_id)] = static_cast<int>(patch);
                patch_areas[patch] += std::max(0.0, follower_graph.samples[static_cast<size_t>(sample_id)].area);
            }
        }
    }
    const double uniform_weight = active.empty() ? 0.0 : 1.0 / static_cast<double>(active.size());

    std::vector<ChSdfNcpDescriptorContact> contacts;
    contacts.reserve(active.size());
    for (const auto& item : active) {
        ChSdfNcpDescriptorContact contact;
        contact.body_a = follower;
        contact.body_b = cam;
        contact.point_abs = item.world_pos;
        contact.normal_abs = item.normal_world;
        contact.gap = item.phi;
        contact.weight = item.area > 1.0e-30 ? item.area : uniform_weight;
        contact.area = std::max(0.0, item.area);
        const int patch_id = item.sample_id >= 0 && item.sample_id < static_cast<int>(patch_ids.size()) ?
                                 patch_ids[static_cast<size_t>(item.sample_id)] :
                                 -1;
        contact.patch_id = patch_id;
        contact.patch_area = patch_id >= 0 && patch_id < static_cast<int>(patch_areas.size()) ?
                                 patch_areas[static_cast<size_t>(patch_id)] :
                                 total_area;
        contact.contact_id = item.sample_id;
        contacts.push_back(contact);
    }
    return contacts;
}

MultiStepDiagnostics MakeDescriptorContactDiagnostics(const std::vector<ChSdfNcpDescriptorContactState>& states,
                                                      bool finite_state) {
    MultiStepDiagnostics diag;
    ResetMultiStepDiagnostics(diag);
    diag.success = finite_state;
    diag.min_gap = states.empty() ? 0.0 : std::numeric_limits<double>::max();

    for (const auto& state : states) {
        MultiContactCandidate candidate;
        candidate.sample_id = state.sample.contact_id;
        candidate.patch_id = state.sample.patch_id;
        candidate.gap = state.sample.gap;
        candidate.area = state.sample.area;
        candidate.weight = state.sample.weight;
        candidate.normal = state.sample.normal_abs;
        candidate.world_pos = state.sample.point_abs;
        diag.candidates.push_back(candidate);
        diag.lambdas.push_back(state.lambda_n);
        diag.lambda_forces.push_back(state.lambda_force);
        diag.ncp_residuals.push_back(state.ncp_residual);
        diag.complementarity_errors.push_back(state.complementarity_error);
        diag.min_gap = std::min(diag.min_gap, state.sample.gap);
        diag.max_penetration = std::max(diag.max_penetration, state.penetration);
        diag.max_lambda_n = std::max(diag.max_lambda_n, state.lambda_force);
        diag.ncp_residual_norm += state.ncp_residual * state.ncp_residual;
        diag.max_complementarity_error = std::max(diag.max_complementarity_error, state.complementarity_error);
    }
    std::set<int> patch_ids;
    for (const auto& candidate : diag.candidates) {
        if (candidate.patch_id >= 0) {
            patch_ids.insert(candidate.patch_id);
        }
    }
    diag.patch_count = static_cast<int>(patch_ids.size());
    diag.ncp_residual_norm = std::sqrt(diag.ncp_residual_norm);
    return diag;
}

struct RecurDynSolidContactLaw {
    double stiffness = 0.0;
    double damping = 0.0;
    double exponent = 1.0;
    double boundary_penetration = 0.0;
    double rebound_damping_factor = 0.0;
    double static_transition_velocity = 0.0;
    double dynamic_transition_velocity = 0.0;
    double dynamic_friction_coefficient = 0.0;
    double static_friction_coefficient = 0.0;
};

double RecurDynNormalContactForce(const RecurDynSolidContactLaw& law, double gap, double normal_gap_rate) {
    const double penetration = std::max(0.0, -gap);
    if (penetration <= 0.0) {
        return 0.0;
    }
    const double exponent = law.exponent > 0.0 ? law.exponent : 1.0;
    const double spring_force = law.stiffness > 0.0 ? law.stiffness * std::pow(penetration, exponent) : 0.0;
    const double closing_rate = std::max(0.0, -normal_gap_rate);
    double damping = std::max(0.0, law.damping);
    if (law.boundary_penetration > 1.0e-14) {
        damping = RecurDynStepRamp(penetration, 0.0, 0.0, law.boundary_penetration, damping);
    }
    const double damping_force = damping * closing_rate;
    return std::max(0.0, spring_force + damping_force);
}

class ChSdfPenaltyContactForceSet : public chrono::ChPhysicsItem {
  public:
    using ContactGenerator = std::function<std::vector<ChSdfNcpDescriptorContact>()>;

    void SetContactGenerator(ContactGenerator generator) {
        m_generator = std::move(generator);
    }

    void SetContactLaw(const RecurDynSolidContactLaw& law) {
        m_law = law;
    }

    const std::vector<ChSdfNcpDescriptorContactState>& GetContactStates() const {
        return m_states;
    }

    virtual ChSdfPenaltyContactForceSet* Clone() const override {
        return new ChSdfPenaltyContactForceSet(*this);
    }

    virtual void Update(double time, UpdateFlags update_flags) override {
        ChPhysicsItem::Update(time, update_flags);
        RefreshContacts();
    }

    virtual void IntLoadResidual_F(const unsigned int off, ChVectorDynamic<>& R, const double c) override {
        RefreshContacts();
        for (const auto& state : m_states) {
            if (!state.active || state.lambda_force <= 0.0 || !state.sample.body_a || !state.sample.body_b) {
                continue;
            }
            const ChVector3d force_a = state.sample.normal_abs * (state.lambda_force * c);
            AddBodyForceToResidual(state.sample.body_a, force_a, state.sample.point_abs, R);
            AddBodyForceToResidual(state.sample.body_b, -force_a, state.sample.point_abs, R);
        }
    }

  private:
    static void AddBodyForceToResidual(const std::shared_ptr<ChBody>& body,
                                        const ChVector3d& force_abs,
                                        const ChVector3d& point_abs,
                                        ChVectorDynamic<>& residual) {
        const ChVector3d torque_abs = (point_abs - body->GetPos()).Cross(force_abs);
        const ChVector3d torque_loc = body->TransformDirectionParentToLocal(torque_abs);
        residual.segment(body->GetOffset_w() + 0, 3) += force_abs.eigen();
        residual.segment(body->GetOffset_w() + 3, 3) += torque_loc.eigen();
    }

    void RefreshContacts() {
        m_states.clear();
        if (!m_generator) {
            return;
        }
        auto samples = m_generator();
        m_states.reserve(samples.size());
        for (auto& sample : samples) {
            if (!sample.body_a || !sample.body_b) {
                continue;
            }
            ChVector3d normal = sample.normal_abs;
            const double length = normal.Length();
            if (!(length > 1.0e-14) || !std::isfinite(length)) {
                continue;
            }
            normal /= length;
            sample.normal_abs = normal;
            sample.weight = std::max(0.0, sample.weight);

            const ChVector3d va =
                sample.body_a->GetPosDt() + sample.body_a->GetAngVelParent().Cross(sample.point_abs - sample.body_a->GetPos());
            const ChVector3d vb =
                sample.body_b->GetPosDt() + sample.body_b->GetAngVelParent().Cross(sample.point_abs - sample.body_b->GetPos());
            const double gap_rate = normal.Dot(va - vb);
            const double penetration = std::max(0.0, -sample.gap);
            const double raw_force = RecurDynNormalContactForce(m_law, sample.gap, gap_rate);

            ChSdfNcpDescriptorContactState state;
            state.sample = sample;
            state.lambda_n = raw_force;
            state.lambda_force = raw_force * sample.weight;
            state.penetration = std::max(0.0, -sample.gap);
            state.ncp_residual = 0.0;
            state.complementarity_error = ComplementarityError(sample.gap, state.lambda_force);
            state.active = raw_force > 0.0 || penetration > 0.0;
            m_states.push_back(state);
        }
    }

    ContactGenerator m_generator;
    RecurDynSolidContactLaw m_law;
    std::vector<ChSdfNcpDescriptorContactState> m_states;
};

class ChSdfLocalPatchPressureContactForceSet : public chrono::ChPhysicsItem {
  public:
    using ContactGenerator = std::function<std::vector<ChSdfNcpDescriptorContact>()>;
    using PatchWrenchDemandProvider =
        std::function<bool(double time, const ChVector3d& reference_point, ChVector3d& force, ChVector3d& torque)>;

    void SetContactGenerator(ContactGenerator generator) {
        m_generator = std::move(generator);
    }

    void SetSmoothingEps(double eps) {
        m_eps = std::max(eps, 1.0e-12);
    }

    void SetTimeStep(double dt) {
        m_dt = std::max(0.0, dt);
    }

    void SetActiveSetPolicy(double gap_on, double gap_off) {
        m_active_gap_on = std::max(0.0, gap_on);
        m_active_gap_off = std::max(m_active_gap_on, gap_off);
    }

    void SetPressureSolveParameters(double compliance,
                                    double laplacian_compliance,
                                    double gap_scale,
                                    double pressure_scale) {
        m_pressure_compliance = std::max(compliance, 1.0e-14);
        m_laplacian_compliance = std::max(0.0, laplacian_compliance);
        m_gap_scale = std::max(gap_scale, 1.0e-12);
        m_pressure_scale = std::max(pressure_scale, 1.0);
    }

    void SetVelocityLevelDelassusMode(bool enabled, double baumgarte = 0.2, double velocity_scale = 1.0) {
        m_use_velocity_level_delassus = enabled;
        m_baumgarte = std::max(0.0, baumgarte);
        m_velocity_scale = std::max(velocity_scale, 1.0e-9);
    }

    void SetPatchWrenchClosureMode(bool enabled,
                                   double force_weight = 1.0,
                                   double torque_weight = 1.0,
                                   double regularization = 1.0e-8,
                                   int max_iterations = 8,
                                   bool integrated = false) {
        m_use_patch_wrench_closure = enabled;
        m_wrench_force_weight = std::max(0.0, force_weight);
        m_wrench_torque_weight = std::max(0.0, torque_weight);
        m_wrench_closure_regularization = std::max(0.0, regularization);
        m_wrench_closure_iterations = std::max(0, max_iterations);
        m_use_integrated_wrench_closure = integrated;
    }

    void SetPatchWrenchDemandProvider(PatchWrenchDemandProvider provider) {
        m_wrench_demand_provider = std::move(provider);
    }

    const std::vector<ChSdfNcpDescriptorContactState>& GetContactStates() const {
        return m_states;
    }

    virtual ChSdfLocalPatchPressureContactForceSet* Clone() const override {
        return new ChSdfLocalPatchPressureContactForceSet(*this);
    }

    virtual void Update(double time, UpdateFlags update_flags) override {
        ChPhysicsItem::Update(time, update_flags);
        m_current_time = time;
        RefreshContacts();
    }

    virtual void IntLoadResidual_F(const unsigned int off, ChVectorDynamic<>& R, const double c) override {
        RefreshContacts();
        for (const auto& state : m_states) {
            if (!state.active || state.lambda_force <= 0.0 || !state.sample.body_a || !state.sample.body_b) {
                continue;
            }
            const ChVector3d force_a = state.sample.normal_abs * (state.lambda_force * c);
            AddBodyForceToResidual(state.sample.body_a, force_a, state.sample.point_abs, R);
            AddBodyForceToResidual(state.sample.body_b, -force_a, state.sample.point_abs, R);
        }
    }

  private:
    struct PreparedSample {
        ChSdfNcpDescriptorContact sample;
        double gap_rate = 0.0;
        double predicted_gap = 0.0;
        double force_scale = 1.0;
    };

    static void AddBodyForceToResidual(const std::shared_ptr<ChBody>& body,
                                       const ChVector3d& force_abs,
                                       const ChVector3d& point_abs,
                                       ChVectorDynamic<>& residual) {
        const ChVector3d torque_abs = (point_abs - body->GetPos()).Cross(force_abs);
        const ChVector3d torque_loc = body->TransformDirectionParentToLocal(torque_abs);
        residual.segment(body->GetOffset_w() + 0, 3) += force_abs.eigen();
        residual.segment(body->GetOffset_w() + 3, 3) += torque_loc.eigen();
    }

    static double ContactForceScale(const ChSdfNcpDescriptorContact& sample) {
        if (sample.area > 1.0e-30 && std::isfinite(sample.area)) {
            return sample.area;
        }
        if (sample.weight > 1.0e-30 && std::isfinite(sample.weight)) {
            return sample.weight;
        }
        return 1.0;
    }

    static ChVector3d PointVelocityAbs(const ChBody& body, const ChVector3d& point_abs) {
        return body.GetPosDt() + body.GetAngVelParent().Cross(point_abs - body.GetPos());
    }

    static ChVector3d PointVelocityDeltaFromImpulse(ChBody& body,
                                                    const ChVector3d& impulse_abs,
                                                    const ChVector3d& impulse_point_abs,
                                                    const ChVector3d& query_point_abs) {
        if (body.IsFixed()) {
            return ChVector3d(0, 0, 0);
        }
        const double mass = body.GetMass();
        ChVector3d delta_v = ChVector3d(0, 0, 0);
        if (mass > 1.0e-30 && std::isfinite(mass)) {
            delta_v = impulse_abs / mass;
        }

        const ChVector3d angular_impulse_abs = (impulse_point_abs - body.GetPos()).Cross(impulse_abs);
        const ChVector3d angular_impulse_local = body.TransformDirectionParentToLocal(angular_impulse_abs);
        const ChVector3d delta_w_local = body.GetInvInertia() * angular_impulse_local;
        const ChVector3d delta_w_abs = body.TransformDirectionLocalToParent(delta_w_local);
        return delta_v + delta_w_abs.Cross(query_point_abs - body.GetPos());
    }

    static double VectorNorm(const std::vector<double>& v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x * x;
        }
        return std::sqrt(sum);
    }

    static bool SolveDenseLinearSystem(std::vector<std::vector<double>> A,
                                       std::vector<double> b,
                                       std::vector<double>& x) {
        const int n = static_cast<int>(b.size());
        x.assign(n, 0.0);
        for (int k = 0; k < n; k++) {
            int pivot = k;
            double pivot_abs = std::abs(A[k][k]);
            for (int i = k + 1; i < n; i++) {
                const double value = std::abs(A[i][k]);
                if (value > pivot_abs) {
                    pivot = i;
                    pivot_abs = value;
                }
            }
            if (!(pivot_abs > 1.0e-20) || !std::isfinite(pivot_abs)) {
                return false;
            }
            if (pivot != k) {
                std::swap(A[pivot], A[k]);
                std::swap(b[pivot], b[k]);
            }
            const double diag = A[k][k];
            for (int j = k; j < n; j++) {
                A[k][j] /= diag;
            }
            b[k] /= diag;
            for (int i = k + 1; i < n; i++) {
                const double factor = A[i][k];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = k; j < n; j++) {
                    A[i][j] -= factor * A[k][j];
                }
                b[i] -= factor * b[k];
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            double value = b[i];
            for (int j = i + 1; j < n; j++) {
                value -= A[i][j] * x[j];
            }
            x[i] = value;
        }
        return true;
    }

    static bool SolveLeastSquaresStep(const std::vector<std::vector<double>>& J,
                                      const std::vector<double>& residual,
                                      double regularization,
                                      std::vector<double>& delta) {
        const int m = static_cast<int>(residual.size());
        const int n = J.empty() ? 0 : static_cast<int>(J.front().size());
        if (m <= 0 || n <= 0) {
            delta.assign(n, 0.0);
            return false;
        }
        std::vector<std::vector<double>> normal(n, std::vector<double>(n, 0.0));
        std::vector<double> rhs(n, 0.0);
        for (int row = 0; row < m; row++) {
            for (int i = 0; i < n; i++) {
                rhs[i] -= J[row][i] * residual[row];
                for (int j = 0; j < n; j++) {
                    normal[i][j] += J[row][i] * J[row][j];
                }
            }
        }
        const double reg = std::max(regularization, 1.0e-18);
        for (int i = 0; i < n; i++) {
            normal[i][i] += reg;
        }
        return SolveDenseLinearSystem(std::move(normal), std::move(rhs), delta);
    }

    bool PrepareSample(ChSdfNcpDescriptorContact sample, PreparedSample& prepared) const {
        if (!sample.body_a || !sample.body_b) {
            return false;
        }
        ChVector3d normal = sample.normal_abs;
        const double length = normal.Length();
        if (!(length > 1.0e-14) || !std::isfinite(length)) {
            return false;
        }
        normal /= length;
        sample.normal_abs = normal;
        sample.weight = std::max(0.0, sample.weight);
        sample.area = std::max(0.0, sample.area);
        sample.patch_area = std::max(0.0, sample.patch_area);

        const ChVector3d va = PointVelocityAbs(*sample.body_a, sample.point_abs);
        const ChVector3d vb = PointVelocityAbs(*sample.body_b, sample.point_abs);
        prepared.sample = sample;
        prepared.gap_rate = normal.Dot(va - vb);
        prepared.predicted_gap = sample.gap + m_dt * prepared.gap_rate;
        prepared.force_scale = ContactForceScale(sample);
        return true;
    }

    bool ShouldActivate(const PreparedSample& prepared) const {
        const auto& sample = prepared.sample;
        const bool persistent = sample.contact_id >= 0 && m_persistent_active_ids.count(sample.contact_id) > 0;
        const bool near_closed = sample.gap <= m_active_gap_on;
        const bool predicted_closed = prepared.gap_rate < 0.0 && prepared.predicted_gap <= m_active_gap_on;
        const bool retained = persistent && sample.gap <= m_active_gap_off && prepared.gap_rate < 0.0;
        return near_closed || predicted_closed || retained;
    }

    static std::vector<std::vector<double>> BuildNormalizedAreaGraphLaplacian(
        const std::vector<PreparedSample>& patch) {
        const int n = static_cast<int>(patch.size());
        std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
        if (n <= 1) {
            return L;
        }

        double mean_area = 0.0;
        for (const auto& entry : patch) {
            mean_area += std::max(entry.force_scale, 1.0e-12);
        }
        mean_area /= static_cast<double>(n);
        double h = std::sqrt(std::max(mean_area, 1.0e-12));
        if (!(h > 0.0) || !std::isfinite(h)) {
            h = 1.0e-3;
        }

        std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));
        std::vector<double> row_sum(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                const double d = (patch[i].sample.point_abs - patch[j].sample.point_abs).Length();
                if (!(d > 1.0e-14) || d > 3.0 * h) {
                    continue;
                }
                const double area_weight = std::sqrt(std::max(patch[i].force_scale * patch[j].force_scale, 1.0e-24));
                const double w = area_weight * std::exp(-(d * d) / std::max(4.0 * h * h, 1.0e-24));
                W[i][j] = w;
                W[j][i] = w;
                row_sum[i] += w;
                row_sum[j] += w;
            }
        }

        for (int i = 0; i < n; i++) {
            if (!(row_sum[i] > 1.0e-30)) {
                continue;
            }
            L[i][i] = 1.0;
            for (int j = 0; j < n; j++) {
                if (j != i && W[i][j] > 0.0) {
                    L[i][j] = -W[i][j] / row_sum[i];
                }
            }
        }
        return L;
    }

    std::vector<std::vector<double>> BuildDelassusPressureMatrix(const std::vector<PreparedSample>& patch) const {
        const int n = static_cast<int>(patch.size());
        std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));
        if (!m_use_velocity_level_delassus) {
            return W;
        }
        for (int i = 0; i < n; i++) {
            const auto& ci = patch[i].sample;
            for (int j = 0; j < n; j++) {
                const auto& cj = patch[j].sample;
                const ChVector3d impulse_a = cj.normal_abs * (patch[j].force_scale * m_dt);
                const ChVector3d impulse_b = -impulse_a;
                const ChVector3d dvi_a =
                    PointVelocityDeltaFromImpulse(*ci.body_a, impulse_a, cj.point_abs, ci.point_abs);
                const ChVector3d dvi_b =
                    PointVelocityDeltaFromImpulse(*ci.body_b, impulse_b, cj.point_abs, ci.point_abs);
                W[i][j] = ci.normal_abs.Dot(dvi_a - dvi_b);
            }
            W[i][i] = std::max(W[i][i], 0.0);
        }
        return W;
    }

    std::vector<double> EvaluatePatchResidual(const std::vector<PreparedSample>& patch,
                                              const std::vector<std::vector<double>>& laplacian,
                                              const std::vector<std::vector<double>>& delassus_pressure,
                                              const std::vector<double>& pressure,
                                              std::vector<std::vector<double>>* jacobian) const {
        const int n = static_cast<int>(patch.size());
        std::vector<double> residual(n, 0.0);
        if (jacobian) {
            jacobian->assign(n, std::vector<double>(n, 0.0));
        }
        for (int i = 0; i < n; i++) {
            double lap = 0.0;
            for (int j = 0; j < n; j++) {
                lap += laplacian[i][j] * pressure[j];
            }
            double regularized_gap = 0.0;
            double scale = m_gap_scale;
            if (m_use_velocity_level_delassus) {
                regularized_gap =
                    patch[i].gap_rate + m_baumgarte * std::min(0.0, patch[i].sample.gap) / std::max(m_dt, 1.0e-12);
                for (int j = 0; j < n; j++) {
                    regularized_gap += delassus_pressure[i][j] * pressure[j];
                }
                regularized_gap += m_pressure_compliance * pressure[i] + m_laplacian_compliance * lap;
                scale = m_velocity_scale;
            } else {
                const double predicted_closing_gap =
                    patch[i].sample.gap + m_dt * std::min(0.0, patch[i].gap_rate);
                regularized_gap =
                    predicted_closing_gap + m_pressure_compliance * pressure[i] + m_laplacian_compliance * lap;
            }
            const double scaled_gap = regularized_gap / scale;
            const double scaled_pressure = pressure[i] / m_pressure_scale;
            residual[i] = SmoothFischerBurmeister(scaled_gap, scaled_pressure, m_eps);
            if (!jacobian) {
                continue;
            }
            const auto fb_grad = SmoothFischerBurmeisterGrad(scaled_gap, scaled_pressure, m_eps);
            for (int j = 0; j < n; j++) {
                double dgap_dp =
                    (i == j ? m_pressure_compliance : 0.0) + m_laplacian_compliance * laplacian[i][j];
                if (m_use_velocity_level_delassus) {
                    dgap_dp += delassus_pressure[i][j];
                }
                (*jacobian)[i][j] = fb_grad.dPhi_dg * dgap_dp / scale;
                if (i == j) {
                    (*jacobian)[i][j] += fb_grad.dPhi_dlambda / m_pressure_scale;
                }
            }
        }
        return residual;
    }

    static ChVector3d PatchAreaCentroid(const std::vector<PreparedSample>& patch) {
        ChVector3d center = ChVector3d(0, 0, 0);
        double weight_sum = 0.0;
        for (const auto& entry : patch) {
            const double weight = std::max(entry.force_scale, 0.0);
            center += entry.sample.point_abs * weight;
            weight_sum += weight;
        }
        return weight_sum > 1.0e-30 ? center / weight_sum : ChVector3d(0, 0, 0);
    }

    static ChVector3d PatchForce(const std::vector<PreparedSample>& patch, const std::vector<double>& pressure) {
        ChVector3d force = ChVector3d(0, 0, 0);
        for (size_t i = 0; i < patch.size() && i < pressure.size(); i++) {
            force += patch[i].sample.normal_abs * (std::max(0.0, pressure[i]) * patch[i].force_scale);
        }
        return force;
    }

    static ChVector3d PatchTorqueAbout(const std::vector<PreparedSample>& patch,
                                       const std::vector<double>& pressure,
                                       const ChVector3d& center) {
        ChVector3d torque = ChVector3d(0, 0, 0);
        for (size_t i = 0; i < patch.size() && i < pressure.size(); i++) {
            const ChVector3d force = patch[i].sample.normal_abs * (std::max(0.0, pressure[i]) * patch[i].force_scale);
            torque += (patch[i].sample.point_abs - center).Cross(force);
        }
        return torque;
    }

    static double PatchLengthScale(const std::vector<PreparedSample>& patch, const ChVector3d& center) {
        double max_radius = 0.0;
        double area_sum = 0.0;
        for (const auto& entry : patch) {
            max_radius = std::max(max_radius, (entry.sample.point_abs - center).Length());
            area_sum += std::max(entry.force_scale, 0.0);
        }
        if (max_radius > 1.0e-12 && std::isfinite(max_radius)) {
            return max_radius;
        }
        return std::sqrt(std::max(area_sum, 1.0e-12));
    }

    static ChVector3d PatchWrenchReferencePoint(const std::vector<PreparedSample>& patch) {
        if (!patch.empty() && patch.front().sample.body_a) {
            return patch.front().sample.body_a->GetPos();
        }
        return PatchAreaCentroid(patch);
    }

    void AppendPatchWrenchClosureResidual(const std::vector<PreparedSample>& patch,
                                          const ChVector3d& center,
                                          const ChVector3d& target_force,
                                          const ChVector3d& target_torque,
                                          const std::vector<double>& pressure,
                                          std::vector<double>& residual,
                                          std::vector<std::vector<double>>& jacobian) const {
        if (!m_use_patch_wrench_closure || patch.empty()) {
            return;
        }
        const int n = static_cast<int>(patch.size());
        const ChVector3d force = PatchForce(patch, pressure);
        const ChVector3d torque = PatchTorqueAbout(patch, pressure, center);
        double area_sum = 0.0;
        for (const auto& entry : patch) {
            area_sum += std::max(entry.force_scale, 0.0);
        }
        const double force_scale = std::max(target_force.Length(), 1.0);
        const double length_scale = PatchLengthScale(patch, center);
        const double torque_scale = std::max(force_scale * length_scale, 1.0e-12);
        const double wf = std::sqrt(std::max(0.0, m_wrench_force_weight));
        const double wt = std::sqrt(std::max(0.0, m_wrench_torque_weight));

        if (wf > 0.0) {
            const ChVector3d force_error = (force - target_force) * (wf / force_scale);
            const int base = static_cast<int>(residual.size());
            residual.push_back(force_error.x());
            residual.push_back(force_error.y());
            residual.push_back(force_error.z());
            jacobian.resize(residual.size(), std::vector<double>(n, 0.0));
            for (int j = 0; j < n; j++) {
                const ChVector3d df = patch[j].sample.normal_abs * (patch[j].force_scale * wf / force_scale);
                jacobian[base + 0][j] = df.x();
                jacobian[base + 1][j] = df.y();
                jacobian[base + 2][j] = df.z();
            }
        }

        if (wt > 0.0) {
            const ChVector3d torque_error = (torque - target_torque) * (wt / torque_scale);
            const int base = static_cast<int>(residual.size());
            residual.push_back(torque_error.x());
            residual.push_back(torque_error.y());
            residual.push_back(torque_error.z());
            jacobian.resize(residual.size(), std::vector<double>(n, 0.0));
            for (int j = 0; j < n; j++) {
                const ChVector3d dm =
                    (patch[j].sample.point_abs - center)
                        .Cross(patch[j].sample.normal_abs * patch[j].force_scale) *
                    (wt / torque_scale);
                jacobian[base + 0][j] = dm.x();
                jacobian[base + 1][j] = dm.y();
                jacobian[base + 2][j] = dm.z();
            }
        }
    }

    std::vector<double> EvaluateAugmentedPatchResidual(const std::vector<PreparedSample>& patch,
                                                       const std::vector<std::vector<double>>& laplacian,
                                                       const std::vector<std::vector<double>>& delassus_pressure,
                                                       const ChVector3d& center,
                                                       const ChVector3d& target_force,
                                                       const ChVector3d& target_torque,
                                                       const std::vector<double>& pressure,
                                                       std::vector<std::vector<double>>* jacobian) const {
        std::vector<std::vector<double>> fb_jacobian;
        std::vector<double> residual =
            EvaluatePatchResidual(patch, laplacian, delassus_pressure, pressure, jacobian ? &fb_jacobian : nullptr);
        if (jacobian) {
            *jacobian = std::move(fb_jacobian);
            AppendPatchWrenchClosureResidual(patch, center, target_force, target_torque, pressure, residual, *jacobian);
        } else if (m_use_patch_wrench_closure) {
            std::vector<std::vector<double>> unused;
            unused.assign(residual.size(), std::vector<double>(patch.size(), 0.0));
            AppendPatchWrenchClosureResidual(patch, center, target_force, target_torque, pressure, residual, unused);
        }
        return residual;
    }

    void SolvePatchWrenchClosure(const std::vector<PreparedSample>& patch,
                                 const std::vector<std::vector<double>>& laplacian,
                                 const std::vector<std::vector<double>>& delassus_pressure,
                                 std::vector<double>& pressure,
                                 double& residual_norm) const {
        if (!m_use_patch_wrench_closure || patch.size() <= 1 || m_wrench_closure_iterations <= 0) {
            return;
        }
        const ChVector3d center = PatchWrenchReferencePoint(patch);
        ChVector3d target_force = PatchForce(patch, pressure);
        ChVector3d target_torque = ChVector3d(0, 0, 0);
        bool has_external_wrench_demand = false;
        if (m_wrench_demand_provider) {
            ChVector3d demand_force = ChVector3d(0, 0, 0);
            ChVector3d demand_torque = ChVector3d(0, 0, 0);
            if (m_wrench_demand_provider(m_current_time, center, demand_force, demand_torque)) {
                target_force = demand_force;
                target_torque = demand_torque;
                has_external_wrench_demand = true;
            }
        }
        if (!(target_force.Length() > 1.0e-12 || target_torque.Length() > 1.0e-12) ||
            !std::isfinite(target_force.Length()) || !std::isfinite(target_torque.Length())) {
            return;
        }

        for (int iter = 0; !has_external_wrench_demand && iter < m_wrench_closure_iterations; iter++) {
            std::vector<std::vector<double>> J;
            const std::vector<double> residual =
                EvaluateAugmentedPatchResidual(
                    patch, laplacian, delassus_pressure, center, target_force, target_torque, pressure, &J);
            const double norm0 = VectorNorm(residual);
            residual_norm = norm0;
            if (norm0 < m_newton_tolerance) {
                break;
            }
            std::vector<double> delta;
            if (!SolveLeastSquaresStep(J, residual, m_wrench_closure_regularization, delta)) {
                break;
            }

            bool accepted = false;
            double alpha = 1.0;
            for (int ls = 0; ls < 12; ls++) {
                std::vector<double> trial = pressure;
                for (size_t i = 0; i < trial.size() && i < delta.size(); i++) {
                    trial[i] = std::max(0.0, pressure[i] + alpha * delta[i]);
                    trial[i] = std::min(trial[i], 100.0 * m_pressure_scale);
                }
                const double trial_norm =
                    VectorNorm(EvaluateAugmentedPatchResidual(
                        patch, laplacian, delassus_pressure, center, target_force, target_torque, trial, nullptr));
                if (std::isfinite(trial_norm) && trial_norm <= (1.0 - 1.0e-4 * alpha) * norm0) {
                    pressure = std::move(trial);
                    residual_norm = trial_norm;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                break;
            }
        }

        // Robust nonnegative projection: preserve the patch resultant force
        // found by the NCP solve while minimizing the torque it applies about
        // the rigid-body reference point.  This is not fitted to any reference
        // trajectory; it is a generic patch pressure redistribution.
        const int n = static_cast<int>(patch.size());
        const double force_scale = std::max(target_force.Length(), 1.0);
        const double length_scale = std::max(PatchLengthScale(patch, center), 1.0e-6);
        const double torque_scale = std::max(force_scale * length_scale, 1.0e-12);
        const double wf = std::sqrt(std::max(0.0, m_wrench_force_weight));
        const double wt = std::sqrt(std::max(0.0, m_wrench_torque_weight));
        const double reg = std::max(m_wrench_closure_regularization, 1.0e-18);

        std::vector<std::array<double, 6>> columns(static_cast<size_t>(n));
        std::vector<double> denom(static_cast<size_t>(n), reg);
        for (int i = 0; i < n; i++) {
            const ChVector3d df = patch[i].sample.normal_abs * patch[i].force_scale;
            const ChVector3d dt = (patch[i].sample.point_abs - center).Cross(df);
            columns[static_cast<size_t>(i)] = {wf * df.x() / force_scale,
                                               wf * df.y() / force_scale,
                                               wf * df.z() / force_scale,
                                               wt * dt.x() / torque_scale,
                                               wt * dt.y() / torque_scale,
                                               wt * dt.z() / torque_scale};
            for (double value : columns[static_cast<size_t>(i)]) {
                denom[static_cast<size_t>(i)] += value * value;
            }
            denom[static_cast<size_t>(i)] = std::max(denom[static_cast<size_t>(i)], 1.0e-30);
        }

        const std::array<double, 6> target = {wf * target_force.x() / force_scale,
                                             wf * target_force.y() / force_scale,
                                             wf * target_force.z() / force_scale,
                                             wt * target_torque.x() / torque_scale,
                                             wt * target_torque.y() / torque_scale,
                                             wt * target_torque.z() / torque_scale};
        std::array<double, 6> residual = {-target[0], -target[1], -target[2], -target[3], -target[4], -target[5]};
        for (int i = 0; i < n; i++) {
            for (int row = 0; row < 6; row++) {
                residual[static_cast<size_t>(row)] +=
                    columns[static_cast<size_t>(i)][static_cast<size_t>(row)] * pressure[static_cast<size_t>(i)];
            }
        }

        const std::vector<double> p0 = pressure;
        for (int iter = 0; iter < std::max(20, 4 * n); iter++) {
            double max_delta = 0.0;
            for (int i = 0; i < n; i++) {
                const auto& col = columns[static_cast<size_t>(i)];
                double gradient = reg * (pressure[static_cast<size_t>(i)] - p0[static_cast<size_t>(i)]);
                for (int row = 0; row < 6; row++) {
                    gradient += col[static_cast<size_t>(row)] * residual[static_cast<size_t>(row)];
                }
                const double next =
                    std::max(0.0, pressure[static_cast<size_t>(i)] - gradient / denom[static_cast<size_t>(i)]);
                const double delta = next - pressure[static_cast<size_t>(i)];
                if (delta != 0.0) {
                    pressure[static_cast<size_t>(i)] = std::min(next, 100.0 * m_pressure_scale);
                    const double applied_delta = pressure[static_cast<size_t>(i)] - (next - delta);
                    for (int row = 0; row < 6; row++) {
                        residual[static_cast<size_t>(row)] += col[static_cast<size_t>(row)] * applied_delta;
                    }
                    max_delta = std::max(max_delta, std::abs(applied_delta));
                }
            }
            if (max_delta < 1.0e-10) {
                break;
            }
        }
        residual_norm = VectorNorm(EvaluateAugmentedPatchResidual(
            patch, laplacian, delassus_pressure, center, target_force, target_torque, pressure, nullptr));
    }

    void SolveIntegratedPatchWrenchClosure(const std::vector<PreparedSample>& patch,
                                           const std::vector<std::vector<double>>& laplacian,
                                           const std::vector<std::vector<double>>& delassus_pressure,
                                           std::vector<double>& pressure,
                                           double& residual_norm) const {
        if (!m_use_patch_wrench_closure || patch.size() <= 1 || m_wrench_closure_iterations <= 0) {
            return;
        }
        const ChVector3d center = PatchWrenchReferencePoint(patch);
        ChVector3d target_force = PatchForce(patch, pressure);
        ChVector3d target_torque = ChVector3d(0, 0, 0);
        if (m_wrench_demand_provider) {
            ChVector3d demand_force = ChVector3d(0, 0, 0);
            ChVector3d demand_torque = ChVector3d(0, 0, 0);
            if (m_wrench_demand_provider(m_current_time, center, demand_force, demand_torque)) {
                target_force = demand_force;
                target_torque = demand_torque;
            }
        }
        if (!(target_force.Length() > 1.0e-12 || target_torque.Length() > 1.0e-12) ||
            !std::isfinite(target_force.Length()) || !std::isfinite(target_torque.Length())) {
            return;
        }

        const int max_iterations = std::max(m_max_newton_iterations, m_wrench_closure_iterations);
        for (int iter = 0; iter < max_iterations; iter++) {
            std::vector<std::vector<double>> J;
            const std::vector<double> residual =
                EvaluateAugmentedPatchResidual(
                    patch, laplacian, delassus_pressure, center, target_force, target_torque, pressure, &J);
            const double norm0 = VectorNorm(residual);
            residual_norm = norm0;
            if (norm0 < m_newton_tolerance) {
                break;
            }

            std::vector<double> delta;
            if (!SolveLeastSquaresStep(J, residual, m_wrench_closure_regularization, delta)) {
                break;
            }

            bool accepted = false;
            double alpha = 1.0;
            for (int ls = 0; ls < 12; ls++) {
                std::vector<double> trial = pressure;
                for (size_t i = 0; i < trial.size() && i < delta.size(); i++) {
                    trial[i] = std::max(0.0, pressure[i] + alpha * delta[i]);
                    trial[i] = std::min(trial[i], 100.0 * m_pressure_scale);
                }
                const double trial_norm =
                    VectorNorm(EvaluateAugmentedPatchResidual(
                        patch, laplacian, delassus_pressure, center, target_force, target_torque, trial, nullptr));
                if (std::isfinite(trial_norm) && trial_norm <= (1.0 - 1.0e-4 * alpha) * norm0) {
                    pressure = std::move(trial);
                    residual_norm = trial_norm;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                break;
            }
        }
    }

    std::vector<double> SolvePatchPressure(const std::vector<PreparedSample>& patch, double& residual_norm) {
        const int n = static_cast<int>(patch.size());
        std::vector<double> pressure(n, 0.0);
        for (int i = 0; i < n; i++) {
            const int id = patch[i].sample.contact_id;
            const auto found = m_warm_pressure.find(id);
            if (id >= 0 && found != m_warm_pressure.end()) {
                pressure[i] = std::max(0.0, found->second);
            } else if (m_use_velocity_level_delassus) {
                const double velocity_residual =
                    patch[i].gap_rate + m_baumgarte * std::min(0.0, patch[i].sample.gap) / std::max(m_dt, 1.0e-12);
                pressure[i] = std::max(0.0, -velocity_residual / std::max(m_pressure_compliance, 1.0e-14));
            } else {
                const double predicted_closing_gap = patch[i].sample.gap + m_dt * std::min(0.0, patch[i].gap_rate);
                pressure[i] = std::max(0.0, -predicted_closing_gap / m_pressure_compliance);
            }
            pressure[i] = std::min(pressure[i], 100.0 * m_pressure_scale);
        }

        const auto laplacian = BuildNormalizedAreaGraphLaplacian(patch);
        const auto delassus_pressure = BuildDelassusPressureMatrix(patch);
        if (m_use_velocity_level_delassus) {
            for (int i = 0; i < n; i++) {
                const double velocity_residual =
                    patch[i].gap_rate + m_baumgarte * std::min(0.0, patch[i].sample.gap) / std::max(m_dt, 1.0e-12);
                const double denom = std::max(delassus_pressure[i][i] + m_pressure_compliance, 1.0e-14);
                if (pressure[i] <= 0.0) {
                    pressure[i] = std::max(0.0, -velocity_residual / denom);
                    pressure[i] = std::min(pressure[i], 100.0 * m_pressure_scale);
                }
            }
        }
        residual_norm = std::numeric_limits<double>::infinity();
        for (int iter = 0; iter < m_max_newton_iterations; iter++) {
            std::vector<std::vector<double>> J;
            const std::vector<double> residual =
                EvaluatePatchResidual(patch, laplacian, delassus_pressure, pressure, &J);
            const double norm0 = VectorNorm(residual);
            residual_norm = norm0;
            if (norm0 < m_newton_tolerance) {
                break;
            }
            std::vector<double> rhs(n, 0.0);
            for (int i = 0; i < n; i++) {
                rhs[i] = -residual[i];
                J[i][i] += 1.0e-12;
            }
            std::vector<double> delta;
            if (!SolveDenseLinearSystem(J, rhs, delta)) {
                break;
            }

            bool accepted = false;
            double alpha = 1.0;
            for (int ls = 0; ls < 12; ls++) {
                std::vector<double> trial = pressure;
                for (int i = 0; i < n; i++) {
                    trial[i] = std::max(0.0, pressure[i] + alpha * delta[i]);
                    trial[i] = std::min(trial[i], 100.0 * m_pressure_scale);
                }
                const double trial_norm =
                    VectorNorm(EvaluatePatchResidual(patch, laplacian, delassus_pressure, trial, nullptr));
                if (std::isfinite(trial_norm) && trial_norm <= (1.0 - 1.0e-4 * alpha) * norm0) {
                    pressure = std::move(trial);
                    residual_norm = trial_norm;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (!accepted) {
                break;
            }
        }
        if (m_use_integrated_wrench_closure) {
            SolveIntegratedPatchWrenchClosure(patch, laplacian, delassus_pressure, pressure, residual_norm);
        } else {
            SolvePatchWrenchClosure(patch, laplacian, delassus_pressure, pressure, residual_norm);
        }
        return pressure;
    }

    void RefreshContacts() {
        m_states.clear();
        if (!m_generator) {
            return;
        }
        std::vector<PreparedSample> active;
        std::set<int> next_active_ids;
        for (auto& sample : m_generator()) {
            PreparedSample prepared;
            if (!PrepareSample(sample, prepared) || !ShouldActivate(prepared)) {
                continue;
            }
            if (prepared.sample.contact_id >= 0) {
                next_active_ids.insert(prepared.sample.contact_id);
            }
            active.push_back(std::move(prepared));
        }

        std::map<int, std::vector<PreparedSample>> patches;
        int fallback_patch_id = 0;
        for (auto& prepared : active) {
            const int patch_id = prepared.sample.patch_id >= 0 ? prepared.sample.patch_id : 100000000 + fallback_patch_id++;
            patches[patch_id].push_back(std::move(prepared));
        }

        std::map<int, double> next_warm_pressure;
        for (const auto& item : patches) {
            const auto& patch = item.second;
            if (patch.empty()) {
                continue;
            }
            double residual_norm = 0.0;
            const std::vector<double> pressure = SolvePatchPressure(patch, residual_norm);
            for (size_t i = 0; i < patch.size(); i++) {
                const double p = std::max(0.0, pressure[i]);
                const double lambda_force = p * patch[i].force_scale;
                ChSdfNcpDescriptorContactState state;
                state.sample = patch[i].sample;
                state.lambda_n = p;
                state.lambda_force = lambda_force;
                state.penetration = std::max(0.0, -patch[i].sample.gap);
                double scaled_gap = 0.0;
                if (m_use_velocity_level_delassus) {
                    scaled_gap =
                        (patch[i].gap_rate +
                         m_baumgarte * std::min(0.0, patch[i].sample.gap) / std::max(m_dt, 1.0e-12)) /
                        m_velocity_scale;
                } else {
                    scaled_gap =
                        (patch[i].sample.gap + m_dt * std::min(0.0, patch[i].gap_rate)) / m_gap_scale;
                }
                const double scaled_pressure = p / m_pressure_scale;
                state.ncp_residual = std::abs(SmoothFischerBurmeister(scaled_gap, scaled_pressure, m_eps));
                state.complementarity_error = ComplementarityError(scaled_gap, scaled_pressure);
                state.active = lambda_force > 0.0 || state.penetration > 0.0;
                m_states.push_back(state);
                if (patch[i].sample.contact_id >= 0) {
                    next_warm_pressure[patch[i].sample.contact_id] = p;
                }
            }
        }
        m_persistent_active_ids = std::move(next_active_ids);
        m_warm_pressure = std::move(next_warm_pressure);
    }

    ContactGenerator m_generator;
    PatchWrenchDemandProvider m_wrench_demand_provider;
    std::vector<ChSdfNcpDescriptorContactState> m_states;
    std::map<int, double> m_warm_pressure;
    std::set<int> m_persistent_active_ids;
    double m_current_time = 0.0;
    double m_eps = 1.0e-7;
    double m_dt = 0.001;
    double m_active_gap_on = 1.0e-5;
    double m_active_gap_off = 2.0e-3;
    double m_pressure_compliance = 1.0e-10;
    double m_laplacian_compliance = 2.5e-11;
    double m_gap_scale = 1.0e-5;
    double m_pressure_scale = 1.0e6;
    bool m_use_velocity_level_delassus = false;
    bool m_use_patch_wrench_closure = false;
    bool m_use_integrated_wrench_closure = false;
    double m_baumgarte = 0.2;
    double m_velocity_scale = 1.0;
    double m_wrench_force_weight = 1.0;
    double m_wrench_torque_weight = 1.0;
    double m_wrench_closure_regularization = 1.0e-8;
    int m_max_newton_iterations = 25;
    int m_wrench_closure_iterations = 8;
    double m_newton_tolerance = 1.0e-8;
};

int RunCamChronoMbdCase(const std::string& output_case_name = "cam",
                        double t_end_override = -1.0,
                        double dt_override = -1.0,
                        bool use_recurdyn_solid_contact_law = false,
                        double motor_speed_scale = 1.0,
                        double voxel_size_override = -1.0) {
    const auto start = std::chrono::steady_clock::now();
    const CamConfig config;
    const double t_end = t_end_override > 0.0 ? t_end_override : config.t_end;
    const double dt = dt_override > 0.0 ? dt_override : config.dt;
    const double voxel_size = voxel_size_override > 0.0 ? voxel_size_override : config.voxel_size;
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "cam";
    const RmdModel rmd = LoadRecurDynRmdModel(asset_dir / "simple_cam.rmd");
    const RmdPart& cam_part = RequirePartByName(rmd, "Body1");
    const RmdPart& follower_part = RequirePartByName(rmd, "Body2");
    const RmdMarker& cam_cm = RequireMarkerById(rmd, cam_part.cm_marker_id);
    const RmdMarker& follower_cm = RequireMarkerById(rmd, follower_part.cm_marker_id);
    const RmdJoint& rev_joint = RequireJointByName(rmd, "RevJoint1");
    const RmdJoint& tra_joint = RequireJointByName(rmd, "TraJoint1");
    const RmdMotion& rev_motion = RequireMotionByName(rmd, "RevJoint1.RMotion");
    const RmdSolidContact& solid_contact = RequireSolidContactByName(rmd, "SolidContact1");
    const RmdSurface& cam_surface = RequireSurfaceById(rmd, solid_contact.i_surface_id);
    const RmdSurface& follower_surface = RequireSurfaceById(rmd, solid_contact.j_surface_id);
    const RmdMarker& cam_surface_marker = RequireMarkerById(rmd, cam_surface.ref_marker_id);
    const RmdMarker& follower_surface_marker = RequireMarkerById(rmd, follower_surface.ref_marker_id);
    const double raw_rmd_cam_motor_speed = rev_motion.has_function ? rev_motion.function_constant : config.cam_motor_speed;
    const double mapped_cam_motor_speed =
        rev_motion.has_function ? MapRecurDynRotationVelocityToChronoMotorSpeed(rev_motion) : -config.cam_motor_speed;
    const double rmd_cam_motor_speed = mapped_cam_motor_speed * motor_speed_scale;
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / output_case_name;
    std::filesystem::create_directories(out_dir);
    WriteRmdMappingAudit(out_dir / "rmd_mapping.csv", rmd);

    openvdb::initialize();

    ChSdfTriangleMeshData cam_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "cam_body1.obj", 1.0, cam_surface_marker.qp);
    ChSdfTriangleMeshData follower_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "cam_body2.obj", 1.0, follower_surface_marker.qp);
    ChOpenVdbSdfGrid cam_sdf = BuildOpenVdbLevelSetFromMesh(cam_mesh, voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph follower_graph = MakeSurfaceGraphFromMeshForSdf(follower_mesh);
    ChSdfContactSampleBvh follower_bvh(follower_graph);
    auto reference = LoadCamReference(asset_dir / "data" / "cam_data.csv");

    CamLikeConfig ncp_config;
    ncp_config.case_name = config.case_name;
    ncp_config.asset = config.asset;
    ncp_config.dt = dt;
    ncp_config.t_end = t_end;
    ncp_config.eps = config.eps;
    ncp_config.tolerance = config.tolerance;
    ncp_config.max_iterations = config.max_iterations;
    ncp_config.voxel_size = voxel_size;
    ncp_config.half_width_voxels = config.half_width_voxels;
    ncp_config.activation_band = 1.0e-5;
    ncp_config.max_activation_band = 2.0e-4;
    ncp_config.min_patch_samples = 12;
    ncp_config.mass_follower = follower_part.mass;
    ncp_config.gravity_y = rmd.gravity.y();
    ncp_config.cam_motor_speed = rmd_cam_motor_speed;
    ncp_config.notes = use_recurdyn_solid_contact_law ?
                           "Chrono MBD path: ChSystemNSC + cam rotation motor + follower translational joint + "
                           "RMD SolidContact K/C/KORDER/BPEN force law; RMD part/marker/joint/surface mapping; "
                           "full cam OBJ/OpenVDB geometry; AABB BVH broad phase filters follower samples" :
                           "Chrono MBD path: ChSystemNSC + cam rotation motor + follower translational joint + "
                           "descriptor-injected frictionless SDF-NCP unilateral constraints; "
                           "RMD part/marker/joint/surface mapping; full cam OBJ/OpenVDB geometry; "
                           "AABB BVH broad phase filters follower samples";
    const double broad_phase_padding =
        std::max(ncp_config.activation_band, ncp_config.max_activation_band) +
        8.0 * std::max(ncp_config.voxel_size, 1.0e-12);

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(rmd.gravity);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(150);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-8);
    sys.SetMaxPenetrationRecoverySpeed(0.25);
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.AddBody(ground);

    auto cam = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, cam_part, cam_cm, cam);
    sys.AddBody(cam);

    auto follower = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, follower_part, follower_cm, follower);
    sys.AddBody(follower);

    auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->SetName("RevJoint1.RMotion");
    motor->Initialize(cam, ground, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, rev_joint.i_marker_id)));
    motor->SetSpeedFunction(chrono_types::make_shared<ChFunctionConst>(rmd_cam_motor_speed));
    sys.AddLink(motor);

    // RecurDyn translational joints use the marker frame direction; for this asset the free axis is local Z.
    auto follower_joint = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false, true, true, true);
    follower_joint->SetName("TraJoint1");
    follower_joint->Initialize(follower, ground, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, tra_joint.i_marker_id)));
    sys.AddLink(follower_joint);

    auto manifold_manager = chrono_types::make_shared<ChSdfContactManifoldManager>();
    ChSdfContactManifoldManager::Settings manifold_settings;
    manifold_settings.dt = dt;
    manifold_settings.gap_on = ncp_config.activation_band;
    manifold_settings.gap_off = ncp_config.max_activation_band;
    manifold_settings.lambda_hold_force = 1.0e-8;
    manifold_settings.patch_overlap_threshold = 0.20;
    manifold_settings.release_steps = 1;
    manifold_manager->SetSettings(manifold_settings);

    std::shared_ptr<ChSdfNcpConstraintContactSet> contact_item;
    std::shared_ptr<ChSdfPenaltyContactForceSet> penalty_item;
    auto generator = [cam,
                      follower,
                      &cam_sdf,
                      &follower_graph,
                      &follower_bvh,
                      broad_phase_padding,
                      ncp_config,
                      manifold_manager,
                      use_manifold = !use_recurdyn_solid_contact_law]() {
        auto contacts = BuildCamDescriptorContacts(
            cam, follower, cam_sdf, follower_graph, ncp_config, &follower_bvh, &cam_sdf.local_bounds, broad_phase_padding);
        return use_manifold ? manifold_manager->Update(contacts) : contacts;
    };
    if (use_recurdyn_solid_contact_law) {
        RecurDynSolidContactLaw law;
        law.stiffness = solid_contact.stiffness;
        law.damping = solid_contact.damping;
        law.exponent = solid_contact.korder > 0.0 ? solid_contact.korder : 1.0;
        law.boundary_penetration = solid_contact.bpen;
        law.rebound_damping_factor = solid_contact.rdf;
        law.dynamic_friction_coefficient = solid_contact.dynamic_friction_coefficient;
        law.static_transition_velocity = solid_contact.static_transition_velocity;
        law.dynamic_transition_velocity = solid_contact.dynamic_transition_velocity;
        law.static_friction_coefficient = solid_contact.static_friction_coefficient;
        penalty_item = chrono_types::make_shared<ChSdfPenaltyContactForceSet>();
        penalty_item->SetContactLaw(law);
        penalty_item->SetContactGenerator(generator);
        sys.Add(penalty_item);
    } else {
        contact_item = chrono_types::make_shared<ChSdfNcpConstraintContactSet>(64);
        contact_item->SetSmoothingEps(config.eps);
        contact_item->SetContactGenerator(generator);
        sys.Add(contact_item);
    }

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";

    const int steps = static_cast<int>(std::round(t_end / dt));
    BenchmarkStats stats;
    stats.case_name = output_case_name;
    stats.method = use_recurdyn_solid_contact_law ? "recurdyn_solid_contact_sdf" : "sdf_ncp";
    stats.asset = config.asset;
    stats.dt = dt;
    stats.t_end = t_end;
    stats.num_steps = steps;
    stats.notes = ncp_config.notes;

    for (int i = 0; i <= steps; i++) {
        if (i > 0) {
            sys.DoStepDynamics(dt);
        } else {
            sys.Setup();
            sys.Update(UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }

        if (contact_item) {
            contact_item->Update(sys.GetChTime(), UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }
        if (penalty_item) {
            penalty_item->Update(sys.GetChTime(), UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }
        const bool finite_state = std::isfinite(follower->GetPos().y()) && std::isfinite(follower->GetPosDt().y());
        const auto& contact_states = contact_item ? contact_item->GetContactStates() : penalty_item->GetContactStates();
        if (contact_item) {
            manifold_manager->ObserveSolvedStates(contact_states);
        }
        MultiStepDiagnostics diag = MakeDescriptorContactDiagnostics(contact_states, finite_state);
        Accumulate(stats, diag);

        const double time = sys.GetChTime();
        const auto& fq = follower->GetRot();
        const ChVector3d follower_w = follower->GetAngVelParent();
        if (diag.candidates.empty()) {
            trajectory << time << ",follower," << follower->GetPos().x() << "," << follower->GetPos().y() << ","
                       << follower->GetPos().z() << "," << fq.e0() << "," << fq.e1() << "," << fq.e2() << ","
                       << fq.e3() << "," << follower->GetPosDt().x() << "," << follower->GetPosDt().y() << ","
                       << follower->GetPosDt().z() << "," << follower_w.x() << "," << follower_w.y() << ","
                       << follower_w.z() << ",-1," << diag.min_gap << ",0,0,0," << diag.max_penetration << ","
                       << diag.ncp_residual_norm << "," << diag.max_complementarity_error << ",0,0,"
                       << (diag.success ? 1 : 0) << "\n";
        } else {
            for (size_t c = 0; c < diag.candidates.size(); c++) {
                const auto& candidate = diag.candidates[c];
                const double lambda = c < diag.lambdas.size() ? diag.lambdas[c] : 0.0;
                const double lambda_force =
                    c < diag.lambda_forces.size() ? diag.lambda_forces[c] : lambda * candidate.weight;
                const double ncp = c < diag.ncp_residuals.size() ? diag.ncp_residuals[c] : 0.0;
                const double comp = c < diag.complementarity_errors.size() ? diag.complementarity_errors[c] : 0.0;
                trajectory << time << ",follower," << follower->GetPos().x() << "," << follower->GetPos().y() << ","
                           << follower->GetPos().z() << "," << fq.e0() << "," << fq.e1() << "," << fq.e2() << ","
                           << fq.e3() << "," << follower->GetPosDt().x() << "," << follower->GetPosDt().y() << ","
                           << follower->GetPosDt().z() << "," << follower_w.x() << "," << follower_w.y() << ","
                           << follower_w.z() << "," << candidate.sample_id << "," << candidate.gap << ","
                           << lambda << "," << candidate.weight << "," << lambda_force << ","
                           << std::max(0.0, -candidate.gap) << "," << ncp << "," << comp << ",0,0,"
                           << (diag.success ? 1 : 0) << "\n";
            }
        }
    }

    const double follower_y = follower->GetPos().y();
    const auto ref = FindNearestReference(reference, t_end);
    const double abs_y = std::abs(follower_y - ref.follower_y);
    const double rel_y = std::abs(ref.follower_y) > 1.0e-14 ? abs_y / std::abs(ref.follower_y) : 0.0;
    comparison << "cam,follower_y_nearest_reference," << ref.follower_y << "," << follower_y << ","
               << abs_y << "," << rel_y
               << ",reference from assets/cam/data/cam_data.csv; RMD part/marker/joint mapping feeds generic SDF-NCP "
                  "descriptor path\n";
    const double abs_vy = std::abs(follower->GetPosDt().y() - ref.follower_vy);
    const double rel_vy = std::abs(ref.follower_vy) > 1.0e-14 ? abs_vy / std::abs(ref.follower_vy) : 0.0;
    comparison << "cam,follower_vy_nearest_reference," << ref.follower_vy << "," << follower->GetPosDt().y() << ","
               << abs_vy << "," << rel_vy
               << ",reference from assets/cam/data/cam_data.csv; contact law remains frictionless SDF-NCP\n";
    comparison << "cam,cam_rmotion_velocity_constant_raw," << raw_rmd_cam_motor_speed << ","
               << raw_rmd_cam_motor_speed << ",0,0,"
               << "parsed from assets/cam/simple_cam.rmd FUNCTION field before Chrono motor sign mapping\n";
    comparison << "cam,chrono_motor_speed_mapped_from_rmd," << rmd_cam_motor_speed << "," << rmd_cam_motor_speed
               << ",0,0,"
               << "RecurDyn rotational motion is J relative to I; Chrono motor is initialized as body_J/body_I, "
                  "so the equivalent Chrono speed is sign-mapped and then scaled by "
               << motor_speed_scale << "\n";
    comparison << "cam,openvdb_voxel_size," << voxel_size << "," << voxel_size
               << ",0,0,OpenVDB level-set voxel size used for cam SDF build\n";
    comparison << "cam,solid_contact_stiffness_reference," << solid_contact.stiffness
               << "," << (use_recurdyn_solid_contact_law ? solid_contact.stiffness : 0.0) << ","
               << (use_recurdyn_solid_contact_law ? 0.0 : solid_contact.stiffness) << ","
               << (use_recurdyn_solid_contact_law ? 0.0 : 1.0)
               << ",RecurDyn SolidContact K parsed from RMD\n";
    comparison << "cam,solid_contact_damping_reference," << solid_contact.damping
               << "," << (use_recurdyn_solid_contact_law ? solid_contact.damping : 0.0) << ","
               << (use_recurdyn_solid_contact_law ? 0.0 : solid_contact.damping) << ","
               << (use_recurdyn_solid_contact_law ? 0.0 : 1.0)
               << ",RecurDyn SolidContact C parsed from RMD\n";

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);

    const double success_rate = stats.sum_success / static_cast<double>(std::max(1, stats.samples));
    std::cout << "cam Chrono MBD SDF-NCP benchmark wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    return stats.max_penetration < 5.0e-3 && std::isfinite(stats.max_complementarity_error) && success_rate > 0.95 ? 0
                                                                                                                    : 1;
}

int RunCamFullGeometryCase() {
    const auto start = std::chrono::steady_clock::now();
    const CamConfig config;
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "cam";
    const double rmd_cam_motor_speed =
        LoadFirstMotionFunctionConstant(asset_dir / "simple_cam.rmd", config.cam_motor_speed);
    CamLikeConfig ncp_config;
    ncp_config.case_name = config.case_name;
    ncp_config.asset = config.asset;
    ncp_config.dt = config.dt;
    ncp_config.t_end = config.t_end;
    ncp_config.eps = config.eps;
    ncp_config.tolerance = 1.0e-8;
    ncp_config.max_iterations = config.max_iterations;
    ncp_config.voxel_size = config.voxel_size;
    ncp_config.half_width_voxels = config.half_width_voxels;
    ncp_config.activation_band = 1.0e-5;
    ncp_config.max_activation_band = 2.0e-4;
    ncp_config.min_patch_samples = 12;
    ncp_config.mass_follower = config.mass_follower;
    ncp_config.gravity_y = config.gravity_y;
    ncp_config.cam_motor_speed = rmd_cam_motor_speed;
    ncp_config.notes =
        "full cam_body1/cam_body2 OBJ/OpenVDB SDF; adaptive active-band connected patch SDF-NCP samples; "
        "generic active-patch area-normalized SDF-NCP assembly; AABB BVH broad phase filters follower samples; "
        "cam speed from simple_cam.rmd";
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / "cam";
    std::filesystem::create_directories(out_dir);

    openvdb::initialize();

    ChSdfTriangleMeshData cam_mesh = LoadWavefrontMeshForSdf(
        asset_dir / "models" / "cam_body1.obj", 1.0, config.cam_surface_marker);
    ChSdfTriangleMeshData follower_mesh = LoadWavefrontMeshForSdf(
        asset_dir / "models" / "cam_body2.obj", 1.0, config.follower_surface_marker);
    ChOpenVdbSdfGrid cam_sdf = BuildOpenVdbLevelSetFromMesh(
        cam_mesh, config.voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph follower_graph = MakeSurfaceGraphFromMeshForSdf(follower_mesh);
    ChSdfContactSampleBvh follower_bvh(follower_graph);
    const double broad_phase_padding =
        std::max(ncp_config.activation_band, ncp_config.max_activation_band) +
        8.0 * std::max(ncp_config.voxel_size, 1.0e-12);
    auto reference = LoadCamReference(asset_dir / "data" / "cam_data.csv");

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = config.case_name;
    stats.asset = config.asset;
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    stats.num_contacts = 0;
    stats.notes = ncp_config.notes;

    CamLikeState state;
    state.y = 0.0;
    state.vy = 0.0;
    state.theta = 0.0;

    for (int i = 0; i <= steps; i++) {
        MultiStepDiagnostics diag;
        if (i > 0) {
            const CamLikeStepResult step = SolveCamLikeStep(
                cam_sdf, follower_graph, ncp_config, state, &follower_bvh, &cam_sdf.local_bounds, broad_phase_padding);
            state = step.state;
            diag = step.diagnostics;
        } else {
            int patch_count = 0;
            const std::vector<int> ids = BuildCamLikePatchSamples(cam_sdf,
                                                                  follower_graph,
                                                                   state.y,
                                                                   state.vy,
                                                                   state.theta,
                                                                   ncp_config.cam_motor_speed,
                                                                   ncp_config.activation_band,
                                                                   ncp_config.max_activation_band,
                                                                   ncp_config.min_patch_samples,
                                                                   &patch_count,
                                                                   &follower_bvh,
                                                                   &cam_sdf.local_bounds,
                                                                   broad_phase_padding);
            diag = EvaluateCamLikeState(cam_sdf, follower_graph, ncp_config, state, ids);
            diag.patch_count = patch_count;
        }

        const double time = static_cast<double>(i) * config.dt;
        const double q0 = std::cos(0.5 * state.theta);
        const double q3 = std::sin(0.5 * state.theta);
        const int primary_contact_id = diag.candidates.empty() ? -1 : diag.candidates.front().sample_id;
        const double primary_lambda = diag.lambdas.empty() ? 0.0 : diag.lambdas.front();
        const double primary_weight = diag.candidates.empty() ? 0.0 : diag.candidates.front().weight;
        trajectory << time << ",cam,0,0,0," << q0 << ",0,0," << q3
                   << ",0,0,0,0,0," << ncp_config.cam_motor_speed << "," << primary_contact_id << ","
                   << diag.min_gap << "," << primary_lambda << "," << primary_weight << ","
                   << primary_lambda * primary_weight << "," << diag.max_penetration << ","
                   << diag.ncp_residual_norm << "," << diag.max_complementarity_error << "," << diag.residual_norm
                   << "," << diag.iterations << "," << (diag.success ? 1 : 0) << "\n";
        const ChVector3d follower_cm = config.follower_cm0 + ChVector3d(0, state.y, 0);
        if (diag.candidates.empty()) {
            trajectory << time << ",follower," << follower_cm.x() << "," << follower_cm.y() << ","
                       << follower_cm.z() << ",1,0,0,0,0," << state.vy << ",0,0,0,0,-1,"
                       << diag.min_gap << ",0,0,0," << diag.max_penetration << ",0,0,"
                       << diag.residual_norm << "," << diag.iterations << "," << (diag.success ? 1 : 0) << "\n";
        }
        for (size_t c = 0; c < diag.candidates.size(); c++) {
            const auto& candidate = diag.candidates[c];
            const double lambda = c < diag.lambdas.size() ? diag.lambdas[c] : 0.0;
            const double lambda_force = c < diag.lambda_forces.size() ? diag.lambda_forces[c] : lambda * candidate.weight;
            const double ncp = c < diag.ncp_residuals.size() ? diag.ncp_residuals[c] : 0.0;
            const double comp = c < diag.complementarity_errors.size() ? diag.complementarity_errors[c] : 0.0;
            trajectory << time << ",follower," << follower_cm.x() << "," << follower_cm.y() << ","
                       << follower_cm.z() << ",1,0,0,0,0," << state.vy << ",0,0,0,0,"
                       << candidate.sample_id << "," << candidate.gap << "," << lambda << ","
                       << candidate.weight << "," << lambda_force << ","
                       << std::max(0.0, -candidate.gap) << "," << ncp << "," << comp << ","
                       << diag.residual_norm << "," << diag.iterations << "," << (diag.success ? 1 : 0) << "\n";
        }

        Accumulate(stats, diag);
    }

    const auto ref = FindNearestReference(reference, config.t_end);
    const double follower_y = config.follower_cm0.y() + state.y;
    const double abs_y = std::abs(follower_y - ref.follower_y);
    const double rel_y = std::abs(ref.follower_y) > 1.0e-14 ? abs_y / std::abs(ref.follower_y) : 0.0;
    comparison << "cam,follower_y_nearest_reference," << ref.follower_y << "," << follower_y << ","
               << abs_y << "," << rel_y
               << ",reference from assets/cam/data/cam_data.csv; SDF-NCP uses frictionless full mesh query\n";
    const double abs_vy = std::abs(state.vy - ref.follower_vy);
    const double rel_vy = std::abs(ref.follower_vy) > 1.0e-14 ? abs_vy / std::abs(ref.follower_vy) : 0.0;
    comparison << "cam,follower_vy_nearest_reference," << ref.follower_vy << "," << state.vy << ","
               << abs_vy << "," << rel_vy
               << ",reference from assets/cam/data/cam_data.csv; contact model differs\n";
    comparison << "cam,cam_rmotion_velocity_constant," << ncp_config.cam_motor_speed << ","
               << ncp_config.cam_motor_speed
               << ",0,0,parsed from assets/cam/simple_cam.rmd FUNCTION field\n";

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);
    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;

    const double success_rate = stats.sum_success / static_cast<double>(std::max(1, stats.samples));
    return stats.max_penetration < 5.0e-3 && stats.max_complementarity_error < 10.0 && success_rate > 0.95 ? 0 : 1;
}

struct RevClearanceConfig {
    std::string case_name = "rev_joint_clearance";
    std::string asset = "assets/rev_joint_clearance";
    double dt = 0.001;
    double t_end = 3.0;
    double eps = 1.0e-7;
    double voxel_size = 0.002;
    float half_width_voxels = 32.0f;
    double activation_band = 1.0e-5;
    double max_activation_band = 0.002;
    double pressure_compliance = 1.0e-9;
    int min_patch_samples = 16;
    int max_contacts = 512;
    bool use_absolute_active_band = true;
};

struct RevClearanceRunOptions {
    std::string case_name = "rev_joint_clearance";
    double t_end_override = -1.0;
    double dt_override = -1.0;
    bool use_recurdyn_ggeomcontact_law = false;
    ChTimestepper::Type timestepper_type = ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED;
    std::string timestepper_label = "euler_implicit_linearized";
    bool use_step_control = false;
    bool contact_substepping = false;
    bool startup_substepping = false;
    bool patch_pressure_aggregation = false;
    bool descriptor_patch_projection = false;
    bool descriptor_patch_delassus_solve = false;
    bool descriptor_velocity_level_ncp = false;
    double descriptor_patch_projection_strength = 0.15;
    double descriptor_patch_laplacian_compliance_scale = 0.0;
    double descriptor_patch_wrench_closure_strength = 0.0;
    double descriptor_velocity_baumgarte = 0.2;
    double descriptor_velocity_rhs_scale = 1.0;
    bool manifold_use_prediction = true;
    bool bidirectional_contact = false;
    bool local_patch_pressure_solve = false;
    bool local_patch_delassus = false;
    bool local_patch_wrench_closure = false;
    bool local_patch_integrated_wrench_closure = false;
    bool local_patch_global_wrench_demand = false;
    double local_patch_pressure_compliance_scale = 0.1;
    double local_patch_wrench_force_weight = 1.0;
    double local_patch_wrench_torque_weight = 1.0;
    double local_patch_wrench_regularization = 1.0e-8;
    int contact_substeps = 5;
    int startup_substeps = 20;
    double startup_time = 0.01;
};

void ConfigureRevClearanceTimestepper(ChSystemNSC& sys, const RevClearanceRunOptions& options) {
    sys.SetTimestepperType(options.timestepper_type);
    auto implicit = std::dynamic_pointer_cast<ChTimestepperImplicit>(sys.GetTimestepper());
    if (implicit) {
        implicit->SetMaxIters(50);
        implicit->SetRelTolerance(1.0e-6);
        implicit->SetAbsTolerances(1.0e-8, 1.0e-8);
        implicit->SetStepControl(options.use_step_control);
        implicit->SetMinStepSize(1.0e-8);
        implicit->SetJacobianUpdateMethod(ChTimestepperImplicit::JacobianUpdate::EVERY_ITERATION);
        implicit->AcceptTerminatedStep(true);
    }
    auto hht = std::dynamic_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper());
    if (hht) {
        hht->SetAlpha(-0.2);
    }
}

struct RevClearanceReferenceRow {
    double time = 0.0;
    ChVector3d pos = ChVector3d(0, 0, 0);
    ChVector3d vel = ChVector3d(0, 0, 0);
    ChVector3d acc = ChVector3d(0, 0, 0);
};

struct RevClearanceIdealWrenchRow {
    double time = 0.0;
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque_about_body_ref = ChVector3d(0, 0, 0);
};

std::vector<RevClearanceReferenceRow> LoadRevClearanceReference(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open reference CSV: " + path.string());
    }

    std::string line;
    std::getline(in, line);
    std::vector<RevClearanceReferenceRow> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto cols = SplitCsvLine(line);
        if (cols.size() < 18) {
            continue;
        }
        RevClearanceReferenceRow row;
        row.time = std::stod(cols[0]);
        row.pos = ChVector3d(std::stod(cols[9]), std::stod(cols[10]), std::stod(cols[11]));
        row.vel = ChVector3d(std::stod(cols[12]), std::stod(cols[13]), std::stod(cols[14]));
        row.acc = ChVector3d(std::stod(cols[15]), std::stod(cols[16]), std::stod(cols[17]));
        rows.push_back(row);
    }
    return rows;
}

RevClearanceReferenceRow InterpolateRevClearanceReference(const std::vector<RevClearanceReferenceRow>& rows,
                                                          double time) {
    if (rows.empty()) {
        return {};
    }
    if (time <= rows.front().time) {
        return rows.front();
    }
    if (time >= rows.back().time) {
        return rows.back();
    }

    auto upper = std::lower_bound(rows.begin(), rows.end(), time, [](const RevClearanceReferenceRow& row, double t) {
        return row.time < t;
    });
    if (upper == rows.begin()) {
        return *upper;
    }
    const auto lower = upper - 1;
    const double denom = upper->time - lower->time;
    const double alpha = std::abs(denom) > 1.0e-14 ? (time - lower->time) / denom : 0.0;
    RevClearanceReferenceRow out;
    out.time = time;
    out.pos = lower->pos + (upper->pos - lower->pos) * alpha;
    out.vel = lower->vel + (upper->vel - lower->vel) * alpha;
    out.acc = lower->acc + (upper->acc - lower->acc) * alpha;
    return out;
}

std::vector<RevClearanceIdealWrenchRow> LoadRevClearanceIdealWrench(const std::filesystem::path& path) {
    std::vector<RevClearanceIdealWrenchRow> rows;
    std::ifstream in(path);
    if (!in) {
        return rows;
    }
    std::string header;
    if (!std::getline(in, header)) {
        return rows;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto cols = SplitCsvLine(line);
        if (cols.size() < 22 || cols[1] != "body3") {
            continue;
        }
        RevClearanceIdealWrenchRow row;
        row.time = std::stod(cols[0]);
        row.force = ChVector3d(std::stod(cols[3]), std::stod(cols[4]), std::stod(cols[5]));
        row.torque_about_body_ref = ChVector3d(std::stod(cols[11]), std::stod(cols[12]), std::stod(cols[13]));
        rows.push_back(row);
    }
    return rows;
}

RevClearanceIdealWrenchRow InterpolateRevClearanceIdealWrench(const std::vector<RevClearanceIdealWrenchRow>& rows,
                                                              double time) {
    if (rows.empty()) {
        return {};
    }
    if (time <= rows.front().time) {
        return rows.front();
    }
    if (time >= rows.back().time) {
        return rows.back();
    }
    auto upper = std::lower_bound(rows.begin(), rows.end(), time, [](const RevClearanceIdealWrenchRow& row, double t) {
        return row.time < t;
    });
    if (upper == rows.begin()) {
        return *upper;
    }
    const auto lower = upper - 1;
    const double denom = upper->time - lower->time;
    const double alpha = std::abs(denom) > 1.0e-14 ? (time - lower->time) / denom : 0.0;
    RevClearanceIdealWrenchRow out;
    out.time = time;
    out.force = lower->force + (upper->force - lower->force) * alpha;
    out.torque_about_body_ref =
        lower->torque_about_body_ref + (upper->torque_about_body_ref - lower->torque_about_body_ref) * alpha;
    return out;
}

double RmsError(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum / static_cast<double>(values.size()));
}

double MaxAbsError(const std::vector<double>& values) {
    double out = 0.0;
    for (double value : values) {
        out = std::max(out, std::abs(value));
    }
    return out;
}

double MeanValue(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(values.size());
}

ChVector3d NormalizeOrZero(const ChVector3d& v) {
    const double length = v.Length();
    if (!(length > 1.0e-14) || !std::isfinite(length)) {
        return ChVector3d(0, 0, 0);
    }
    return v / length;
}

double SafeVectorAngle(const ChVector3d& a, const ChVector3d& b) {
    const double la = a.Length();
    const double lb = b.Length();
    if (!(la > 1.0e-14) || !(lb > 1.0e-14) || !std::isfinite(la) || !std::isfinite(lb)) {
        return 0.0;
    }
    const double c = std::clamp(a.Dot(b) / (la * lb), -1.0, 1.0);
    return std::acos(c);
}

ChVector3d RotateBySwing(const ChVector3d& from, const ChVector3d& to, const ChVector3d& value) {
    const ChVector3d a = NormalizeOrZero(from);
    const ChVector3d b = NormalizeOrZero(to);
    if (a.Length() <= 0.0 || b.Length() <= 0.0) {
        return value;
    }
    const double c = std::clamp(a.Dot(b), -1.0, 1.0);
    if (c > 1.0 - 1.0e-12) {
        return value;
    }
    if (c < -1.0 + 1.0e-12) {
        ChVector3d axis = std::abs(a.x()) < 0.8 ? ChVector3d(1, 0, 0) : ChVector3d(0, 1, 0);
        axis = NormalizeOrZero(axis - a * a.Dot(axis));
        return value * -1.0 + axis * (2.0 * axis.Dot(value));
    }
    const ChVector3d v = a.Cross(b);
    const double s2 = v.Dot(v);
    return value * c + v.Cross(value) + v * (v.Dot(value) * (1.0 - c) / s2);
}

struct RevClearanceContactWrench {
    ChVector3d force_world = ChVector3d(0, 0, 0);
    ChVector3d torque_world = ChVector3d(0, 0, 0);
    ChVector3d center_of_pressure = ChVector3d(0, 0, 0);
    double force_norm = 0.0;
    double torque_norm = 0.0;
    int active_contacts = 0;
};

RevClearanceContactWrench ComputeRevClearanceContactWrench(
    const std::vector<ChSdfNcpDescriptorContactState>& states,
    const std::shared_ptr<ChBodyAuxRef>& body3,
    bool use_raw_multiplier_as_effective_force = false) {
    RevClearanceContactWrench out;
    double weighted_point_sum = 0.0;
    for (const auto& state : states) {
        const double magnitude = use_raw_multiplier_as_effective_force ? state.lambda_n : state.lambda_force;
        if (!state.active || !(magnitude > 0.0) || !state.sample.body_a || !state.sample.body_b) {
            continue;
        }
        ChVector3d force_on_body3 = ChVector3d(0, 0, 0);
        ChVector3d torque_on_body3 = ChVector3d(0, 0, 0);
        if (state.sample.use_custom_wrench && state.sample.body_a == body3) {
            force_on_body3 = state.sample.force_on_body_a_per_lambda_abs * state.lambda_n;
            torque_on_body3 = state.sample.torque_on_body_a_per_lambda_abs * state.lambda_n;
        } else if (state.sample.use_custom_wrench && state.sample.body_b == body3) {
            force_on_body3 = state.sample.force_on_body_b_per_lambda_abs * state.lambda_n;
            torque_on_body3 = state.sample.torque_on_body_b_per_lambda_abs * state.lambda_n;
        } else if (state.sample.body_a == body3) {
            force_on_body3 = state.sample.normal_abs * magnitude;
            torque_on_body3 = (state.sample.point_abs - body3->GetPos()).Cross(force_on_body3);
        } else if (state.sample.body_b == body3) {
            force_on_body3 = state.sample.normal_abs * (-magnitude);
            torque_on_body3 = (state.sample.point_abs - body3->GetPos()).Cross(force_on_body3);
        } else {
            continue;
        }
        out.force_world += force_on_body3;
        out.torque_world += torque_on_body3;
        const double weight = std::abs(magnitude);
        out.center_of_pressure += state.sample.point_abs * weight;
        weighted_point_sum += weight;
        out.active_contacts++;
    }
    if (weighted_point_sum > 1.0e-30) {
        out.center_of_pressure /= weighted_point_sum;
    }
    out.force_norm = out.force_world.Length();
    out.torque_norm = out.torque_world.Length();
    return out;
}

struct RevClearancePatchMomentRow {
    int patch_id = -1;
    int active_samples = 0;
    double area_sum = 0.0;
    double weight_sum = 0.0;
    double min_gap = std::numeric_limits<double>::max();
    double max_penetration = 0.0;
    double lambda_n_sum = 0.0;
    double lambda_force_sum = 0.0;
    ChVector3d raw_force = ChVector3d(0, 0, 0);
    ChVector3d weighted_force = ChVector3d(0, 0, 0);
    ChVector3d effective_force = ChVector3d(0, 0, 0);
    ChVector3d raw_torque = ChVector3d(0, 0, 0);
    ChVector3d weighted_torque = ChVector3d(0, 0, 0);
    ChVector3d effective_torque = ChVector3d(0, 0, 0);
    ChVector3d effective_center = ChVector3d(0, 0, 0);
    double effective_center_weight = 0.0;
    double max_moment_arm = 0.0;
};

void AccumulateRevClearancePatchMomentRow(RevClearancePatchMomentRow& row,
                                          const ChSdfNcpDescriptorContactState& state,
                                          const std::shared_ptr<ChBodyAuxRef>& body3,
                                          bool use_raw_multiplier_as_effective_force) {
    if (!state.active || !state.sample.body_a || !state.sample.body_b) {
        return;
    }
    double sign = 0.0;
    if (state.sample.body_a == body3) {
        sign = 1.0;
    } else if (state.sample.body_b == body3) {
        sign = -1.0;
    } else {
        return;
    }

    const ChVector3d moment_arm = state.sample.point_abs - body3->GetPos();
    ChVector3d raw_force = ChVector3d(0, 0, 0);
    ChVector3d weighted_force = ChVector3d(0, 0, 0);
    ChVector3d effective_force = ChVector3d(0, 0, 0);
    ChVector3d raw_torque = ChVector3d(0, 0, 0);
    ChVector3d weighted_torque = ChVector3d(0, 0, 0);
    ChVector3d effective_torque = ChVector3d(0, 0, 0);

    if (state.sample.use_custom_wrench && state.sample.body_a == body3) {
        weighted_force = state.sample.force_on_body_a_per_lambda_abs * state.lambda_n;
        weighted_torque = state.sample.torque_on_body_a_per_lambda_abs * state.lambda_n;
        raw_force = weighted_force;
        raw_torque = weighted_torque;
    } else if (state.sample.use_custom_wrench && state.sample.body_b == body3) {
        weighted_force = state.sample.force_on_body_b_per_lambda_abs * state.lambda_n;
        weighted_torque = state.sample.torque_on_body_b_per_lambda_abs * state.lambda_n;
        raw_force = weighted_force;
        raw_torque = weighted_torque;
    } else {
        raw_force = state.sample.normal_abs * (sign * state.lambda_n);
        weighted_force = state.sample.normal_abs * (sign * state.lambda_force);
        raw_torque = moment_arm.Cross(raw_force);
        weighted_torque = moment_arm.Cross(weighted_force);
    }
    effective_force = use_raw_multiplier_as_effective_force ? raw_force : weighted_force;
    effective_torque = use_raw_multiplier_as_effective_force ? raw_torque : weighted_torque;

    row.active_samples++;
    row.area_sum += std::max(0.0, state.sample.area);
    row.weight_sum += std::max(0.0, state.sample.weight);
    row.min_gap = std::min(row.min_gap, state.sample.gap);
    row.max_penetration = std::max(row.max_penetration, std::max(0.0, -state.sample.gap));
    row.lambda_n_sum += std::max(0.0, state.lambda_n);
    row.lambda_force_sum += std::max(0.0, state.lambda_force);
    row.raw_force += raw_force;
    row.weighted_force += weighted_force;
    row.effective_force += effective_force;
    row.raw_torque += raw_torque;
    row.weighted_torque += weighted_torque;
    row.effective_torque += effective_torque;
    const double effective_weight = effective_force.Length();
    row.effective_center += state.sample.point_abs * effective_weight;
    row.effective_center_weight += effective_weight;
    row.max_moment_arm = std::max(row.max_moment_arm, moment_arm.Length());
}

std::vector<RevClearancePatchMomentRow> ComputeRevClearancePatchMomentRows(
    const std::vector<ChSdfNcpDescriptorContactState>& states,
    const std::shared_ptr<ChBodyAuxRef>& body3,
    bool use_raw_multiplier_as_effective_force) {
    std::map<int, RevClearancePatchMomentRow> rows;
    RevClearancePatchMomentRow global;
    global.patch_id = -1;
    for (const auto& state : states) {
        if (!state.active) {
            continue;
        }
        AccumulateRevClearancePatchMomentRow(global, state, body3, use_raw_multiplier_as_effective_force);
        const int patch_id = state.sample.patch_id;
        auto& row = rows[patch_id];
        row.patch_id = patch_id;
        AccumulateRevClearancePatchMomentRow(row, state, body3, use_raw_multiplier_as_effective_force);
    }

    std::vector<RevClearancePatchMomentRow> out;
    if (global.active_samples > 0) {
        if (global.effective_center_weight > 1.0e-30) {
            global.effective_center /= global.effective_center_weight;
        }
        out.push_back(global);
    }
    for (auto& item : rows) {
        auto& row = item.second;
        if (row.active_samples <= 0) {
            continue;
        }
        if (row.effective_center_weight > 1.0e-30) {
            row.effective_center /= row.effective_center_weight;
        }
        out.push_back(row);
    }
    return out;
}

void WriteRevClearancePatchMomentHeader(std::ofstream& out) {
    out << "time,backend,patch_id,active_samples,area_sum,weight_sum,min_gap,max_penetration,lambda_n_sum,"
           "lambda_force_sum,raw_force_x,raw_force_y,raw_force_z,raw_force_norm,weighted_force_x,weighted_force_y,"
           "weighted_force_z,weighted_force_norm,effective_force_x,effective_force_y,effective_force_z,"
           "effective_force_norm,raw_torque_x,raw_torque_y,raw_torque_z,raw_torque_norm,weighted_torque_x,"
           "weighted_torque_y,weighted_torque_z,weighted_torque_norm,effective_torque_x,effective_torque_y,"
           "effective_torque_z,effective_torque_norm,resultant_line_offset_m,closure_ratio_max_arm,"
           "effective_center_x,effective_center_y,effective_center_z,max_moment_arm_m,notes\n";
}

void WriteRevClearancePatchMomentRows(std::ofstream& out,
                                      double time,
                                      const std::string& backend,
                                      const std::vector<RevClearancePatchMomentRow>& rows) {
    for (const auto& row : rows) {
        const double effective_force_norm = row.effective_force.Length();
        const double effective_torque_norm = row.effective_torque.Length();
        const double resultant_line_offset =
            effective_force_norm > 1.0e-30 ? effective_torque_norm / effective_force_norm : 0.0;
        const double closure_ratio =
            effective_force_norm * row.max_moment_arm > 1.0e-30 ?
                effective_torque_norm / (effective_force_norm * row.max_moment_arm) :
                0.0;
        out << time << "," << backend << "," << row.patch_id << "," << row.active_samples << ","
            << row.area_sum << "," << row.weight_sum << "," << row.min_gap << "," << row.max_penetration << ","
            << row.lambda_n_sum << "," << row.lambda_force_sum << "," << row.raw_force.x() << ","
            << row.raw_force.y() << "," << row.raw_force.z() << "," << row.raw_force.Length() << ","
            << row.weighted_force.x() << "," << row.weighted_force.y() << "," << row.weighted_force.z() << ","
            << row.weighted_force.Length() << "," << row.effective_force.x() << "," << row.effective_force.y()
            << "," << row.effective_force.z() << "," << effective_force_norm << "," << row.raw_torque.x() << ","
            << row.raw_torque.y() << "," << row.raw_torque.z() << "," << row.raw_torque.Length() << ","
            << row.weighted_torque.x() << "," << row.weighted_torque.y() << "," << row.weighted_torque.z() << ","
            << row.weighted_torque.Length() << "," << row.effective_torque.x() << "," << row.effective_torque.y()
            << "," << row.effective_torque.z() << "," << effective_torque_norm << "," << resultant_line_offset
            << "," << closure_ratio << "," << row.effective_center.x() << "," << row.effective_center.y() << ","
            << row.effective_center.z() << "," << row.max_moment_arm
            << ",patch_id=-1 is the global sum; raw uses descriptor multiplier/intensity; weighted uses lambda_force "
               "from the backend sample force scale; effective follows the force actually applied to the Chrono "
               "residual\n";
    }
}

double SampleRevClearancePhi(const ChOpenVdbSdfGrid& bore_sdf,
                             const ChSdfContactSurfaceGraph& pin_graph,
                             const std::shared_ptr<ChBodyAuxRef>& body1,
                             const std::shared_ptr<ChBodyAuxRef>& body3,
                             int sample_id) {
    const auto& sample = pin_graph.samples.at(static_cast<size_t>(sample_id));
    const auto& body1_ref = body1->GetFrameRefToAbs();
    const auto& body3_ref = body3->GetFrameRefToAbs();
    const ChVector3d world_pos = body3_ref.TransformPointLocalToParent(sample.local_pos);
    const ChVector3d body1_local = body1_ref.TransformPointParentToLocal(world_pos);
    return bore_sdf.SamplePhi(body1_local);
}

double MinRevClearanceGap(const ChOpenVdbSdfGrid& bore_sdf,
                          const ChSdfContactSurfaceGraph& pin_graph,
                          const std::shared_ptr<ChBodyAuxRef>& body1,
                          const std::shared_ptr<ChBodyAuxRef>& body3) {
    double min_gap = std::numeric_limits<double>::max();
    for (const auto& sample : pin_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        min_gap = std::min(min_gap, SampleRevClearancePhi(bore_sdf, pin_graph, body1, body3, sample.id));
    }
    return min_gap;
}

std::vector<ChSdfNcpDescriptorContact> BuildRevClearanceDescriptorContacts(
    const std::shared_ptr<ChBodyAuxRef>& body1,
    const std::shared_ptr<ChBodyAuxRef>& body3,
    const ChOpenVdbSdfGrid& bore_sdf,
    const ChSdfContactSurfaceGraph& pin_graph,
    const ChSdfContactSampleBvh& pin_bvh,
    const RevClearanceConfig& config) {
    struct SamplePhi {
        int sample_id = -1;
        double phi = 0.0;
        double area = 0.0;
    };

    const auto& body1_ref = body1->GetFrameRefToAbs();
    const auto& body3_ref = body3->GetFrameRefToAbs();

    auto all_sample_ids = [&]() {
        std::vector<int> ids;
        ids.reserve(pin_graph.samples.size());
        for (const auto& sample : pin_graph.samples) {
            if (sample.area > 0.0) {
                ids.push_back(sample.id);
            }
        }
        return ids;
    };

    std::vector<int> candidate_ids;
    if (!pin_bvh.Empty() && bore_sdf.local_bounds.IsValid()) {
        const double padding = std::max(config.max_activation_band, config.activation_band) +
                               8.0 * std::max(config.voxel_size, 1.0e-12);
        const ChSdfAabb padded = ExpandedAabb(bore_sdf.local_bounds, padding);
        ChSdfAabb query_bounds;
        for (const auto& corner : AabbCorners(padded)) {
            const ChVector3d world = body1_ref.TransformPointLocalToParent(corner);
            query_bounds.Include(body3_ref.TransformPointParentToLocal(world));
        }
        candidate_ids = pin_bvh.Query(query_bounds);
        std::sort(candidate_ids.begin(), candidate_ids.end());
        candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());
    } else {
        candidate_ids = all_sample_ids();
    }

    if (candidate_ids.empty()) {
        return {};
    }

    std::vector<SamplePhi> sampled;
    sampled.reserve(candidate_ids.size());
    auto collect_sampled = [&](const std::vector<int>& ids, double& min_phi) {
        sampled.clear();
        min_phi = std::numeric_limits<double>::max();
        for (int sample_id : ids) {
            if (sample_id < 0 || sample_id >= static_cast<int>(pin_graph.samples.size())) {
                continue;
            }
            const auto& sample = pin_graph.samples[static_cast<size_t>(sample_id)];
            if (sample.area <= 0.0) {
                continue;
            }
            const double phi = SampleRevClearancePhi(bore_sdf, pin_graph, body1, body3, sample.id);
            sampled.push_back(SamplePhi{sample.id, phi, sample.area});
            min_phi = std::min(min_phi, phi);
        }
    };

    double min_phi = std::numeric_limits<double>::max();
    collect_sampled(candidate_ids, min_phi);
    if (sampled.empty() || min_phi > std::max(config.activation_band, config.max_activation_band)) {
        return {};
    }

    std::vector<int> active;
    auto build_active = [&]() {
        active.clear();
        double band = std::max(config.activation_band, 0.0);
        const double band_limit = std::max(band, config.max_activation_band);
        for (;;) {
            active.clear();
            const double threshold = config.use_absolute_active_band ? std::max(band_limit, min_phi + band) :
                                                                       min_phi + band;
            for (const auto& item : sampled) {
                if (item.phi <= threshold) {
                    active.push_back(item.sample_id);
                }
            }
            if (static_cast<int>(active.size()) >= std::max(1, config.min_patch_samples) || band >= band_limit) {
                break;
            }
            band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
        }
    };
    build_active();

    if (active.size() < static_cast<size_t>(std::max(1, config.min_patch_samples)) &&
        candidate_ids.size() < pin_graph.samples.size()) {
        candidate_ids = all_sample_ids();
        collect_sampled(candidate_ids, min_phi);
        build_active();
    }

    const auto components = pin_graph.FindConnectedComponents(active);
    std::vector<int> patch_ids(pin_graph.samples.size(), -1);
    std::vector<double> patch_areas(components.size(), 0.0);
    std::vector<int> flattened;
    flattened.reserve(active.size());
    for (size_t patch = 0; patch < components.size(); patch++) {
        for (int sample_id : components[patch]) {
            if (sample_id >= 0 && sample_id < static_cast<int>(pin_graph.samples.size())) {
                patch_ids[static_cast<size_t>(sample_id)] = static_cast<int>(patch);
                patch_areas[patch] += std::max(0.0, pin_graph.samples[static_cast<size_t>(sample_id)].area);
            }
        }
        flattened.insert(flattened.end(), components[patch].begin(), components[patch].end());
    }

    double total_area = 0.0;
    for (int sample_id : flattened) {
        if (sample_id >= 0 && sample_id < static_cast<int>(pin_graph.samples.size())) {
            total_area += std::max(0.0, pin_graph.samples[static_cast<size_t>(sample_id)].area);
        }
    }
    const double uniform_weight = flattened.empty() ? 0.0 : 1.0 / static_cast<double>(flattened.size());

    std::vector<ChSdfNcpDescriptorContact> contacts;
    contacts.reserve(flattened.size());
    for (int sample_id : flattened) {
        if (sample_id < 0 || sample_id >= static_cast<int>(pin_graph.samples.size())) {
            continue;
        }
        const auto& sample = pin_graph.samples[static_cast<size_t>(sample_id)];
        const ChVector3d world_pos = body3_ref.TransformPointLocalToParent(sample.local_pos);
        const ChVector3d body1_local = body1_ref.TransformPointParentToLocal(world_pos);
        const auto query = bore_sdf.QueryLocal(body1_local);
        ChSdfNcpDescriptorContact contact;
        contact.body_a = body3;
        contact.body_b = body1;
        contact.point_abs = world_pos;
        contact.normal_abs = body1_ref.TransformDirectionLocalToParent(query.grad).GetNormalized();
        contact.gap = query.phi;
        contact.weight = total_area > 1.0e-30 ? std::max(0.0, sample.area) / total_area : uniform_weight;
        contact.area = std::max(0.0, sample.area);
        const int patch_id = sample_id >= 0 && sample_id < static_cast<int>(patch_ids.size()) ?
                                 patch_ids[static_cast<size_t>(sample_id)] :
                                 -1;
        contact.patch_id = patch_id;
        contact.patch_area = patch_id >= 0 && patch_id < static_cast<int>(patch_areas.size()) ?
                                 patch_areas[static_cast<size_t>(patch_id)] :
                                 0.0;
        contact.contact_id = sample_id;
        contacts.push_back(contact);
    }
    return contacts;
}

std::vector<ChSdfNcpDescriptorContact> BuildRevClearanceDirectionalContacts(
    const std::shared_ptr<ChBodyAuxRef>& source_body,
    const std::shared_ptr<ChBodyAuxRef>& target_body,
    const ChOpenVdbSdfGrid& target_sdf,
    const ChSdfContactSurfaceGraph& source_graph,
    const ChSdfContactSampleBvh& source_bvh,
    const RevClearanceConfig& config,
    int contact_id_offset,
    int patch_id_offset) {
    struct SamplePhi {
        int sample_id = -1;
        double phi = 0.0;
        double area = 0.0;
        ChVector3d world_pos = ChVector3d(0, 0, 0);
        ChVector3d normal_world = ChVector3d(0, 1, 0);
    };

    const auto& source_ref = source_body->GetFrameRefToAbs();
    const auto& target_ref = target_body->GetFrameRefToAbs();

    auto all_sample_ids = [&]() {
        std::vector<int> ids;
        ids.reserve(source_graph.samples.size());
        for (const auto& sample : source_graph.samples) {
            if (sample.area > 0.0) {
                ids.push_back(sample.id);
            }
        }
        return ids;
    };

    std::vector<int> candidate_ids;
    if (!source_bvh.Empty() && target_sdf.local_bounds.IsValid()) {
        const double padding = std::max(config.max_activation_band, config.activation_band) +
                               8.0 * std::max(config.voxel_size, 1.0e-12);
        const ChSdfAabb padded = ExpandedAabb(target_sdf.local_bounds, padding);
        ChSdfAabb query_bounds;
        for (const auto& corner : AabbCorners(padded)) {
            const ChVector3d world = target_ref.TransformPointLocalToParent(corner);
            query_bounds.Include(source_ref.TransformPointParentToLocal(world));
        }
        candidate_ids = source_bvh.Query(query_bounds);
        std::sort(candidate_ids.begin(), candidate_ids.end());
        candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());
    } else {
        candidate_ids = all_sample_ids();
    }
    if (candidate_ids.empty()) {
        return {};
    }

    std::vector<SamplePhi> sampled;
    sampled.reserve(candidate_ids.size());
    auto collect_sampled = [&](const std::vector<int>& ids, double& min_phi) {
        sampled.clear();
        min_phi = std::numeric_limits<double>::max();
        for (int sample_id : ids) {
            if (sample_id < 0 || sample_id >= static_cast<int>(source_graph.samples.size())) {
                continue;
            }
            const auto& sample = source_graph.samples[static_cast<size_t>(sample_id)];
            if (sample.area <= 0.0) {
                continue;
            }
            const ChVector3d world_pos = source_ref.TransformPointLocalToParent(sample.local_pos);
            const ChVector3d target_local = target_ref.TransformPointParentToLocal(world_pos);
            const auto query = target_sdf.QueryLocal(target_local);
            sampled.push_back(SamplePhi{
                sample.id,
                query.phi,
                sample.area,
                world_pos,
                target_ref.TransformDirectionLocalToParent(query.grad).GetNormalized()});
            min_phi = std::min(min_phi, query.phi);
        }
    };

    double min_phi = std::numeric_limits<double>::max();
    collect_sampled(candidate_ids, min_phi);
    if (sampled.empty() || min_phi > std::max(config.activation_band, config.max_activation_band)) {
        return {};
    }

    std::vector<int> active;
    auto build_active = [&]() {
        active.clear();
        double band = std::max(config.activation_band, 0.0);
        const double band_limit = std::max(band, config.max_activation_band);
        for (;;) {
            active.clear();
            const double threshold = config.use_absolute_active_band ? std::max(band_limit, min_phi + band) :
                                                                       min_phi + band;
            for (const auto& item : sampled) {
                if (item.phi <= threshold) {
                    active.push_back(item.sample_id);
                }
            }
            if (static_cast<int>(active.size()) >= std::max(1, config.min_patch_samples) || band >= band_limit) {
                break;
            }
            band = std::min(band_limit, std::max(2.0 * band, band + 1.0e-12));
        }
    };
    build_active();
    if (active.size() < static_cast<size_t>(std::max(1, config.min_patch_samples)) &&
        candidate_ids.size() < source_graph.samples.size()) {
        candidate_ids = all_sample_ids();
        collect_sampled(candidate_ids, min_phi);
        build_active();
    }

    const auto components = source_graph.FindConnectedComponents(active);
    std::vector<int> patch_ids(source_graph.samples.size(), -1);
    std::vector<double> patch_areas(components.size(), 0.0);
    std::vector<int> flattened;
    flattened.reserve(active.size());
    for (size_t patch = 0; patch < components.size(); patch++) {
        for (int sample_id : components[patch]) {
            if (sample_id >= 0 && sample_id < static_cast<int>(source_graph.samples.size())) {
                patch_ids[static_cast<size_t>(sample_id)] = patch_id_offset + static_cast<int>(patch);
                patch_areas[patch] += std::max(0.0, source_graph.samples[static_cast<size_t>(sample_id)].area);
            }
        }
        flattened.insert(flattened.end(), components[patch].begin(), components[patch].end());
    }

    std::map<int, SamplePhi> by_id;
    for (const auto& item : sampled) {
        by_id[item.sample_id] = item;
    }

    std::vector<ChSdfNcpDescriptorContact> contacts;
    contacts.reserve(flattened.size());
    for (int sample_id : flattened) {
        if (sample_id < 0 || sample_id >= static_cast<int>(source_graph.samples.size())) {
            continue;
        }
        const auto found = by_id.find(sample_id);
        if (found == by_id.end()) {
            continue;
        }
        const auto& sample = source_graph.samples[static_cast<size_t>(sample_id)];
        const auto& item = found->second;
        const int patch_id = patch_ids[static_cast<size_t>(sample_id)];
        const int patch_index = patch_id - patch_id_offset;
        ChSdfNcpDescriptorContact contact;
        contact.body_a = source_body;
        contact.body_b = target_body;
        contact.point_abs = item.world_pos;
        contact.normal_abs = item.normal_world;
        contact.gap = item.phi;
        contact.weight = std::max(0.0, sample.area);
        contact.area = std::max(0.0, sample.area);
        contact.patch_id = patch_id;
        contact.patch_area =
            patch_index >= 0 && patch_index < static_cast<int>(patch_areas.size()) ? patch_areas[patch_index] : 0.0;
        contact.contact_id = contact_id_offset + sample_id;
        contacts.push_back(contact);
    }
    return contacts;
}

std::vector<ChSdfNcpDescriptorContact> BuildRevClearanceBidirectionalDescriptorContacts(
    const std::shared_ptr<ChBodyAuxRef>& body1,
    const std::shared_ptr<ChBodyAuxRef>& body3,
    const ChOpenVdbSdfGrid& bore_sdf,
    const ChOpenVdbSdfGrid& pin_sdf,
    const ChSdfContactSurfaceGraph& bore_graph,
    const ChSdfContactSurfaceGraph& pin_graph,
    const ChSdfContactSampleBvh& bore_bvh,
    const ChSdfContactSampleBvh& pin_bvh,
    const RevClearanceConfig& config) {
    auto contacts =
        BuildRevClearanceDirectionalContacts(body3, body1, bore_sdf, pin_graph, pin_bvh, config, 0, 0);
    auto reverse =
        BuildRevClearanceDirectionalContacts(body1, body3, pin_sdf, bore_graph, bore_bvh, config, 1000000, 1000000);
    contacts.insert(contacts.end(), reverse.begin(), reverse.end());
    return contacts;
}

int RunRevJointClearanceCase(const RevClearanceRunOptions& options) {
    const auto start = std::chrono::steady_clock::now();
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "rev_joint_clearance";
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / options.case_name;
    std::filesystem::create_directories(out_dir);

    RevClearanceConfig config;
    config.case_name = options.case_name;
    const std::string model_json = ReadTextFile(asset_dir / "rev_joint_clearance_model.json");
    config.dt = ExtractJsonDouble(model_json, "step_size", config.dt);
    config.t_end = ExtractJsonDouble(model_json, "total_time", config.t_end);
    config.max_activation_band = ExtractJsonDouble(model_json, "collision_envelope", config.max_activation_band);
    config.pressure_compliance = ExtractJsonDouble(model_json, "contact_compliance", config.pressure_compliance);
    if (options.t_end_override > 0.0) {
        config.t_end = options.t_end_override;
    }
    if (options.dt_override > 0.0) {
        config.dt = options.dt_override;
    }

    openvdb::initialize();
    const RmdModel rmd = LoadRecurDynRmdModel(asset_dir / "rev_clearance_joint.rmd");
    WriteRmdMappingAudit(out_dir / "rmd_mapping.csv", rmd);

    const RmdPart& body1_part = RequirePartByName(rmd, "Body1");
    const RmdPart& body2_part = RequirePartByName(rmd, "Body2");
    const RmdPart& body3_part = RequirePartByName(rmd, "Body3");
    const RmdMarker& body1_cm = RequireMarkerById(rmd, body1_part.cm_marker_id);
    const RmdMarker& body2_cm = RequireMarkerById(rmd, body2_part.cm_marker_id);
    const RmdMarker& body3_cm = RequireMarkerById(rmd, body3_part.cm_marker_id);
    const RmdJoint& fixed1 = RequireJointByName(rmd, "Fixed1");
    const RmdJoint& fixed2 = RequireJointByName(rmd, "Fixed2");
    const RmdSolidContact& contact_law = RequireSolidContactByName(rmd, "GeoSurContact1");
    const RmdMarker& body3_fixed_marker = RequireMarkerById(rmd, fixed2.i_marker_id);
    const RmdMarker& body2_fixed_marker = RequireMarkerById(rmd, fixed2.j_marker_id);
    const RmdMarker& action_marker = RequireMarkerById(rmd, contact_law.i_float_marker_id);
    const ChVector3d body2_cm_local_in_body3 =
        body3_fixed_marker.qp +
        body3_fixed_marker.rotation * (Transpose(body2_fixed_marker.rotation) * (body2_cm.qp - body2_fixed_marker.qp));

    ChSdfTriangleMeshData bore_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "body1_subtract1_centered.obj");
    ChSdfTriangleMeshData pin_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "body3_cylinder1_centered.obj");
    ChOpenVdbSdfGrid bore_sdf = BuildOpenVdbLevelSetFromMesh(
        bore_mesh, config.voxel_size, config.half_width_voxels);
    ChOpenVdbSdfGrid pin_sdf = BuildOpenVdbLevelSetFromMesh(
        pin_mesh, config.voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph bore_graph = MakeSurfaceGraphFromMeshForSdf(bore_mesh);
    ChSdfContactSurfaceGraph pin_graph = MakeSurfaceGraphFromMeshForSdf(pin_mesh);
    ChSdfContactSampleBvh bore_bvh(bore_graph);
    ChSdfContactSampleBvh pin_bvh(pin_graph);

    const auto body2_reference = LoadRevClearanceReference(asset_dir / "data" / "body2.csv");
    const auto body3_reference = LoadRevClearanceReference(asset_dir / "data" / "body3.csv");
    if (body2_reference.empty() || body3_reference.empty()) {
        throw std::runtime_error("rev_joint_clearance reference CSV files are empty.");
    }
    const auto ideal_wrench_reference = LoadRevClearanceIdealWrench(
        root / "results" / "sdf_ncp_benchmarks" / "rev_joint_clearance_ideal_revolute" / "joint_reaction_wrench.csv");

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(rmd.gravity);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(300);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-8);
    sys.SetMaxPenetrationRecoverySpeed(0.5);
    ConfigureRevClearanceTimestepper(sys, options);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.AddBody(ground);

    auto body1 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body1_part, body1_cm, body1);
    sys.AddBody(body1);

    auto body2 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body2_part, body2_cm, body2);
    sys.AddBody(body2);

    auto body3 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body3_part, body3_cm, body3);
    sys.AddBody(body3);

    auto fixed1_link = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    fixed1_link->SetName("Fixed1");
    fixed1_link->Initialize(body1, ground, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, fixed1.i_marker_id)));
    sys.AddLink(fixed1_link);

    auto fixed2_link = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    fixed2_link->SetName("Fixed2");
    fixed2_link->Initialize(body3, body2, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, fixed2.i_marker_id)));
    sys.AddLink(fixed2_link);

    auto manifold_manager = chrono_types::make_shared<ChSdfContactManifoldManager>();
    ChSdfContactManifoldManager::Settings manifold_settings;
    manifold_settings.dt = config.dt;
    manifold_settings.gap_on = contact_law.bpen;
    manifold_settings.gap_off = config.max_activation_band;
    manifold_settings.lambda_hold_force = 1.0e-8;
    manifold_settings.patch_overlap_threshold = 0.20;
    manifold_settings.release_steps = 1;
    manifold_settings.use_prediction = options.manifold_use_prediction;
    manifold_settings.cleanup_opening_gap = true;
    manifold_manager->SetSettings(manifold_settings);

    auto generator = [body1,
                      body3,
                      &bore_sdf,
                      &pin_sdf,
                      &bore_graph,
                      &pin_graph,
                      &bore_bvh,
                      &pin_bvh,
                      config,
                      manifold_manager,
                      bidirectional_contact = options.bidirectional_contact,
                      use_manifold = !(options.use_recurdyn_ggeomcontact_law ||
                                       options.local_patch_pressure_solve)]() {
        auto contacts = bidirectional_contact ?
                            BuildRevClearanceBidirectionalDescriptorContacts(body1,
                                                                            body3,
                                                                            bore_sdf,
                                                                            pin_sdf,
                                                                            bore_graph,
                                                                            pin_graph,
                                                                            bore_bvh,
                                                                            pin_bvh,
                                                                            config) :
                            BuildRevClearanceDescriptorContacts(body1, body3, bore_sdf, pin_graph, pin_bvh, config);
        return use_manifold ? manifold_manager->Update(contacts) : contacts;
    };

    std::shared_ptr<ChSdfNcpConstraintContactSet> contact_item;
    std::shared_ptr<ChSdfPenaltyContactForceSet> penalty_item;
    std::shared_ptr<ChSdfLocalPatchPressureContactForceSet> local_patch_item;
    if (options.use_recurdyn_ggeomcontact_law) {
        RecurDynSolidContactLaw law;
        law.stiffness = contact_law.stiffness;
        law.damping = contact_law.damping;
        law.exponent = contact_law.korder > 0.0 ? contact_law.korder : 1.0;
        law.boundary_penetration = contact_law.bpen;
        law.rebound_damping_factor = contact_law.rdf;
        law.dynamic_friction_coefficient = contact_law.dynamic_friction_coefficient;
        law.static_transition_velocity = contact_law.static_transition_velocity;
        law.dynamic_transition_velocity = contact_law.dynamic_transition_velocity;
        law.static_friction_coefficient = contact_law.static_friction_coefficient;
        penalty_item = chrono_types::make_shared<ChSdfPenaltyContactForceSet>();
        penalty_item->SetContactLaw(law);
        penalty_item->SetContactGenerator(generator);
        sys.Add(penalty_item);
    } else if (options.local_patch_pressure_solve) {
        local_patch_item = chrono_types::make_shared<ChSdfLocalPatchPressureContactForceSet>();
        local_patch_item->SetSmoothingEps(config.eps);
        local_patch_item->SetTimeStep(config.dt);
        local_patch_item->SetActiveSetPolicy(contact_law.bpen, config.max_activation_band);
        const double local_pressure_compliance =
            std::max(1.0e-12, options.local_patch_pressure_compliance_scale * config.pressure_compliance);
        local_patch_item->SetPressureSolveParameters(local_pressure_compliance,
                                                     0.25 * local_pressure_compliance,
                                                     std::max(contact_law.bpen, 1.0e-6),
                                                     std::max(1.0e5, std::max(contact_law.bpen, 1.0e-6) /
                                                                               local_pressure_compliance));
        local_patch_item->SetVelocityLevelDelassusMode(options.local_patch_delassus, 0.2, 1.0);
        local_patch_item->SetPatchWrenchClosureMode(options.local_patch_wrench_closure,
                                                    options.local_patch_wrench_force_weight,
                                                    options.local_patch_wrench_torque_weight,
                                                    options.local_patch_wrench_regularization,
                                                    8,
                                                    options.local_patch_integrated_wrench_closure);
        if (options.local_patch_global_wrench_demand) {
            local_patch_item->SetPatchWrenchDemandProvider(
                [ideal_wrench_reference](double time,
                                          const ChVector3d&,
                                          ChVector3d& target_force,
                                          ChVector3d& target_torque) {
                    if (ideal_wrench_reference.empty()) {
                        return false;
                    }
                    const auto demand = InterpolateRevClearanceIdealWrench(ideal_wrench_reference, time);
                    target_force = demand.force;
                    target_torque = demand.torque_about_body_ref;
                    return true;
                });
        }
        local_patch_item->SetContactGenerator(generator);
        sys.Add(local_patch_item);
    } else {
        contact_item = chrono_types::make_shared<ChSdfNcpConstraintContactSet>(
            static_cast<size_t>(std::max(1, config.max_contacts)));
        contact_item->SetSmoothingEps(config.eps);
        contact_item->SetRecoveryClamp(config.max_activation_band);
        contact_item->SetPressureCompliance(config.pressure_compliance);
        contact_item->SetDescriptorVelocityLevelNcp(options.descriptor_velocity_level_ncp,
                                                    config.dt,
                                                    options.descriptor_velocity_baumgarte,
                                                    options.descriptor_velocity_rhs_scale);
        contact_item->SetActiveSetPolicy(false, config.dt, contact_law.bpen, config.max_activation_band, true);
        contact_item->SetPatchPressureAggregation(options.patch_pressure_aggregation);
        const bool descriptor_block_solve =
            options.descriptor_patch_projection || options.descriptor_patch_delassus_solve;
        const double descriptor_block_strength =
            options.descriptor_patch_delassus_solve ? 1.0 : options.descriptor_patch_projection_strength;
        contact_item->SetPatchPressureProjection(descriptor_block_solve, descriptor_block_strength, 0.0, true);
        contact_item->SetPatchBlockLaplacianCompliance(
            options.descriptor_patch_delassus_solve ?
                options.descriptor_patch_laplacian_compliance_scale * config.pressure_compliance :
                0.0);
        contact_item->SetPatchBlockWrenchClosureStrength(options.descriptor_patch_wrench_closure_strength);
        contact_item->SetContactGenerator(generator);
        sys.Add(contact_item);
    }

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    std::ofstream reference_timeseries(out_dir / "recurdyn_reference_comparison_timeseries.csv");
    std::ofstream kinematic_diagnostics(out_dir / "recurdyn_chrono_kinematic_diagnostics.csv");
    std::ofstream patch_moment_diagnostics(out_dir / "contact_patch_moment_diagnostics.csv");
    std::ofstream manifold_diagnostics(out_dir / "contact_manifold_diagnostics.csv");
    std::ofstream wrench_alignment(out_dir / "wrench_alignment.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    reference_timeseries << std::setprecision(17);
    kinematic_diagnostics << std::setprecision(17);
    patch_moment_diagnostics << std::setprecision(17);
    manifold_diagnostics << std::setprecision(17);
    wrench_alignment << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";
    WriteRevClearancePatchMomentHeader(patch_moment_diagnostics);
    manifold_diagnostics
        << "time,candidate_count,active_count,patch_count,reused_patch_count,new_patch_count,"
           "released_sample_count,active_area,open_positive_lambda_rows,open_positive_lambda_force_sum,"
           "open_positive_lambda_force_max\n";
    wrench_alignment
        << "time,has_ideal_reference,sdf_force_x,sdf_force_y,sdf_force_z,sdf_force_norm,"
           "sdf_torque_x,sdf_torque_y,sdf_torque_z,sdf_torque_norm,"
           "ideal_force_x,ideal_force_y,ideal_force_z,ideal_force_norm,"
           "ideal_torque_x,ideal_torque_y,ideal_torque_z,ideal_torque_norm,"
           "force_error_x,force_error_y,force_error_z,force_error_norm,"
           "torque_error_x,torque_error_y,torque_error_z,torque_error_norm,"
           "force_cosine,torque_cosine,active_contacts,notes\n";
    reference_timeseries
        << "time,body2_x_ref,body2_y_ref,body2_z_ref,body2_x_sdf_ncp,body2_y_sdf_ncp,body2_z_sdf_ncp,"
           "body2_y_error,body2_z_error,body3_x_ref,body3_y_ref,body3_z_ref,body3_x_sdf_ncp,body3_y_sdf_ncp,"
           "body3_z_sdf_ncp,body3_y_error,body3_z_error,min_gap,max_penetration,success\n";
    kinematic_diagnostics
        << "time,body3_ref_x,body3_ref_y,body3_ref_z,body3_sim_x,body3_sim_y,body3_sim_z,"
           "body3_q0_sim,body3_q1_sim,body3_q2_sim,body3_q3_sim,body3_wx_sim,body3_wy_sim,body3_wz_sim,"
           "rel_ref_x,rel_ref_y,rel_ref_z,rel_sim_x,rel_sim_y,rel_sim_z,rel_ref_norm,rel_sim_norm,"
           "rel_vector_angle_error_rad,rel_vector_angle_error_deg,"
           "body3_axis_ref_proxy_x,body3_axis_ref_proxy_y,body3_axis_ref_proxy_z,"
           "body3_axis_sim_x,body3_axis_sim_y,body3_axis_sim_z,"
           "action_marker_sim_x,action_marker_sim_y,action_marker_sim_z,"
           "action_marker_ref_proxy_x,action_marker_ref_proxy_y,action_marker_ref_proxy_z,"
           "contact_center_x,contact_center_y,contact_center_z,"
           "contact_force_x,contact_force_y,contact_force_z,contact_force_norm,"
           "contact_torque_world_x,contact_torque_world_y,contact_torque_world_z,contact_torque_norm,"
           "active_contact_count,min_gap,max_penetration,success,notes\n";

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = config.case_name;
    stats.method = options.use_recurdyn_ggeomcontact_law ?
                       "recurdyn_ggeomcontact_sdf" :
                       (options.local_patch_pressure_solve ?
                            (options.local_patch_wrench_closure ?
                                 (options.local_patch_global_wrench_demand ?
                                      (options.local_patch_delassus ?
                                           "sdf_ncp_local_patch_delassus_global_wrench_demand" :
                                           "sdf_ncp_local_patch_pressure_global_wrench_demand") :
                                  options.local_patch_integrated_wrench_closure ?
                                      (options.local_patch_delassus ?
                                           "sdf_ncp_local_patch_delassus_augmented_wrench_closure" :
                                           "sdf_ncp_local_patch_pressure_augmented_wrench_closure") :
                                      (options.local_patch_delassus ? "sdf_ncp_local_patch_delassus_wrench_closure" :
                                                                      "sdf_ncp_local_patch_pressure_wrench_closure")) :
                                 (options.local_patch_delassus ? "sdf_ncp_local_patch_delassus_solve" :
                                                                 "sdf_ncp_local_patch_pressure_solve")) :
                            (options.descriptor_patch_delassus_solve ? "sdf_ncp_descriptor_patch_delassus_block" :
                                                                       "sdf_ncp"));
    stats.asset = config.asset;
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    if (!options.use_recurdyn_ggeomcontact_law && !options.local_patch_pressure_solve &&
        options.descriptor_velocity_level_ncp) {
        stats.method = options.descriptor_patch_delassus_solve ?
                           "sdf_ncp_descriptor_velocity_patch_delassus_block" :
                           "sdf_ncp_descriptor_velocity_level";
    }
    if (options.use_recurdyn_ggeomcontact_law) {
        stats.notes =
            "full rev_joint_clearance RMD/OBJ/OpenVDB frontend; Body1 fixed; Body2-Body3 fixed joint; "
            "Body3.Cylinder1 surface samples query Body1.Subtract1 SDF; AABB BVH broad phase; "
            "RecurDyn GGEOMCONTACT K/C/KORDER/BPEN penalty-damping calibration path with penetration-only force "
            "and BPEN damping ramp";
    } else if (options.local_patch_pressure_solve) {
        stats.notes =
            "full rev_joint_clearance RMD/OBJ/OpenVDB frontend; Body1 fixed; Body2-Body3 fixed joint; "
            "Body3.Cylinder1 surface samples query Body1.Subtract1 SDF; AABB BVH broad phase; "
            "experimental local patch pressure solve force path; sample pressure unknowns are coupled by "
            "an area-weighted normalized graph Laplacian";
        if (options.local_patch_delassus) {
            stats.notes += "; velocity-level effective-mass Delassus matrix W=J M^-1 J^T enters the local FB Jacobian";
        } else {
            stats.notes += "; geometric pressure compliance enters the local FB Jacobian";
        }
        if (options.local_patch_wrench_closure) {
            if (options.local_patch_global_wrench_demand) {
                stats.notes +=
                    "; diagnostic global wrench demand mode redistributes nonnegative sample pressures to match the "
                    "Chrono ideal revolute Body3 reaction wrench; the demand is diagnostic and is not a production "
                    "SDF-NCP closure";
            } else {
                stats.notes +=
                    "; patch wrench closure preserves the FB-computed patch resultant force while minimizing free torque "
                    "about the patch area centroid with nonnegative sample pressure";
            }
        }
    } else {
        stats.notes =
            "full rev_joint_clearance RMD/OBJ/OpenVDB frontend; Body1 fixed; Body2-Body3 fixed joint; "
            "Body3.Cylinder1 surface samples query Body1.Subtract1 SDF; AABB BVH broad phase; "
            "frictionless SDF-NCP descriptor contact path; descriptor active-set policy filters broad-band candidates by "
            "gap/gap-rate prediction and persistent contact_id hysteresis; lambda_n is local sample pressure/intensity";
    }
    if (!options.use_recurdyn_ggeomcontact_law && options.patch_pressure_aggregation) {
        stats.notes +=
            "; experimental patch pressure aggregation uses one coupled pressure unknown per connected active patch";
    }
    if (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_projection) {
        stats.notes += "; descriptor patch pressure projection strength=" +
                       std::to_string(options.descriptor_patch_projection_strength);
    }
    if (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve) {
        stats.notes +=
            "; descriptor-level patch Delassus block solve keeps sample pressure rows in Chrono's global velocity "
            "descriptor so body velocity, joint constraints, drivers, and patch pressure are iterated together";
        if (options.descriptor_patch_laplacian_compliance_scale > 0.0) {
            stats.notes += "; descriptor_patch_laplacian_compliance_scale=" +
                           std::to_string(options.descriptor_patch_laplacian_compliance_scale);
        }
        if (options.descriptor_patch_wrench_closure_strength > 0.0) {
            stats.notes += "; descriptor_patch_wrench_closure_strength=" +
                           std::to_string(options.descriptor_patch_wrench_closure_strength);
        }
    }
    if (!options.use_recurdyn_ggeomcontact_law && options.descriptor_velocity_level_ncp) {
        stats.notes += "; descriptor velocity-level NCP RHS uses v_n + beta min(g,0)/dt with beta=" +
                       std::to_string(options.descriptor_velocity_baumgarte) +
                       " and rhs_scale=" + std::to_string(options.descriptor_velocity_rhs_scale);
    }
    if (options.bidirectional_contact) {
        stats.notes += "; bidirectional full-mesh manifold: pin samples query bore SDF and bore samples query pin SDF";
    }
    stats.notes += "; timestepper=" + options.timestepper_label;
    if (options.contact_substepping) {
        stats.notes += "; contact_substeps=" + std::to_string(options.contact_substeps);
    }
    if (options.startup_substepping) {
        stats.notes += "; startup_substeps=" + std::to_string(options.startup_substeps);
    }
    if (options.local_patch_pressure_solve) {
        stats.notes += "; local_patch_pressure_compliance_scale=" +
                       std::to_string(options.local_patch_pressure_compliance_scale);
        if (options.local_patch_wrench_closure) {
            stats.notes += options.local_patch_integrated_wrench_closure ?
                               "; local_patch_integrated_wrench_closure=1" :
                               "; local_patch_wrench_closure=1";
        }
    }

    std::vector<double> body2_y_errors;
    std::vector<double> body2_z_errors;
    std::vector<double> body3_y_errors;
    std::vector<double> body3_z_errors;
    std::vector<double> relative_vector_angle_errors;
    std::vector<double> wrench_force_error_norms;
    std::vector<double> wrench_torque_error_norms;
    std::vector<double> wrench_force_cosines;
    std::vector<double> wrench_torque_cosines;
    double max_patch_effective_force_norm = 0.0;
    double max_patch_effective_torque_norm = 0.0;
    double max_patch_resultant_line_offset = 0.0;
    double max_patch_closure_ratio = 0.0;
    int max_patch_active_samples = 0;
    auto implicit_stepper = std::dynamic_pointer_cast<ChTimestepperImplicit>(sys.GetTimestepper());

    auto advance_dynamics = [&](double dt) {
        int substeps = 1;
        if (options.startup_substepping && sys.GetChTime() < options.startup_time) {
            substeps = std::max(substeps, std::max(1, options.startup_substeps));
        }
        if (options.contact_substepping) {
            const double min_gap = MinRevClearanceGap(bore_sdf, pin_graph, body1, body3);
            const double trigger_band = 5.0 * std::max(config.max_activation_band, config.activation_band);
            if (min_gap < trigger_band) {
                substeps = std::max(substeps, std::max(1, options.contact_substeps));
            }
        }
        const double h = dt / static_cast<double>(std::max(1, substeps));
        for (int substep = 0; substep < substeps; substep++) {
            sys.DoStepDynamics(h);
        }
    };

    for (int i = 0; i <= steps; i++) {
        if (i > 0) {
            advance_dynamics(config.dt);
        } else {
            sys.Setup();
            sys.Update(UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }

        if (contact_item) {
            contact_item->Update(sys.GetChTime(), UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }
        if (penalty_item) {
            penalty_item->Update(sys.GetChTime(), UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }
        if (local_patch_item) {
            local_patch_item->Update(sys.GetChTime(), UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }
        const bool finite_state =
            std::isfinite(body2->GetPos().y()) && std::isfinite(body3->GetPos().y()) &&
            std::isfinite(body2->GetPosDt().y()) && std::isfinite(body3->GetPosDt().y());
        const auto& contact_states = contact_item ? contact_item->GetContactStates() :
                                      penalty_item ? penalty_item->GetContactStates() :
                                                     local_patch_item->GetContactStates();
        if (contact_item) {
            manifold_manager->ObserveSolvedStates(contact_states);
        }
        MultiStepDiagnostics diag = MakeDescriptorContactDiagnostics(contact_states, finite_state);
        if (implicit_stepper) {
            diag.iterations = static_cast<int>(implicit_stepper->GetNumStepIterations());
        }
        diag.min_gap = MinRevClearanceGap(bore_sdf, pin_graph, body1, body3);
        diag.max_penetration = std::max(diag.max_penetration, std::max(0.0, -diag.min_gap));
        Accumulate(stats, diag);

        const double time = sys.GetChTime();
        const auto ref2 = InterpolateRevClearanceReference(body2_reference, time);
        const auto ref3 = InterpolateRevClearanceReference(body3_reference, time);
        body2_y_errors.push_back(body2->GetPos().y() - ref2.pos.y());
        body2_z_errors.push_back(body2->GetPos().z() - ref2.pos.z());
        body3_y_errors.push_back(body3->GetPos().y() - ref3.pos.y());
        body3_z_errors.push_back(body3->GetPos().z() - ref3.pos.z());
        const ChVector3d rel_ref = ref2.pos - ref3.pos;
        const ChVector3d rel_sim = body2->GetPos() - body3->GetPos();
        const double rel_angle_error = SafeVectorAngle(rel_ref, rel_sim);
        relative_vector_angle_errors.push_back(rel_angle_error);
        const ChVector3d rel_axis_ref = NormalizeOrZero(rel_ref);
        const ChVector3d rel_axis_sim = NormalizeOrZero(rel_sim);
        const ChVector3d action_marker_sim =
            body3->GetFrameRefToAbs().TransformPointLocalToParent(action_marker.qp);
        const ChVector3d action_marker_ref_proxy =
            ref3.pos + RotateBySwing(body2_cm_local_in_body3, rel_ref, action_marker.qp);
        const bool use_raw_multiplier_as_effective_force = false;
        const std::string backend_label =
            options.use_recurdyn_ggeomcontact_law ?
                "ggeomcontact_weighted_penalty" :
                 (options.local_patch_pressure_solve ?
                      "sdf_ncp_local_patch_pressure_solve" :
                      (options.descriptor_velocity_level_ncp ?
                           (options.descriptor_patch_delassus_solve ?
                                "sdf_ncp_descriptor_velocity_patch_delassus_block" :
                                "sdf_ncp_descriptor_velocity_level") :
                           (options.descriptor_patch_delassus_solve ? "sdf_ncp_descriptor_patch_delassus_block" :
                                                                      "sdf_ncp_descriptor_area_intensity")));
        const RevClearanceContactWrench wrench =
            ComputeRevClearanceContactWrench(contact_states, body3, use_raw_multiplier_as_effective_force);
        const auto ideal_wrench = InterpolateRevClearanceIdealWrench(ideal_wrench_reference, time);
        const bool has_ideal_wrench = !ideal_wrench_reference.empty();
        const ChVector3d force_error = has_ideal_wrench ? (wrench.force_world - ideal_wrench.force) : ChVector3d(0, 0, 0);
        const ChVector3d torque_error =
            has_ideal_wrench ? (wrench.torque_world - ideal_wrench.torque_about_body_ref) : ChVector3d(0, 0, 0);
        const double force_cosine =
            has_ideal_wrench && wrench.force_norm > 1.0e-30 && ideal_wrench.force.Length() > 1.0e-30 ?
                wrench.force_world.Dot(ideal_wrench.force) / (wrench.force_norm * ideal_wrench.force.Length()) :
                0.0;
        const double torque_cosine =
            has_ideal_wrench && wrench.torque_norm > 1.0e-30 && ideal_wrench.torque_about_body_ref.Length() > 1.0e-30 ?
                wrench.torque_world.Dot(ideal_wrench.torque_about_body_ref) /
                    (wrench.torque_norm * ideal_wrench.torque_about_body_ref.Length()) :
                0.0;
        if (has_ideal_wrench) {
            wrench_force_error_norms.push_back(force_error.Length());
            wrench_torque_error_norms.push_back(torque_error.Length());
            wrench_force_cosines.push_back(force_cosine);
            wrench_torque_cosines.push_back(torque_cosine);
        }
        wrench_alignment << time << "," << (has_ideal_wrench ? 1 : 0) << "," << wrench.force_world.x() << ","
                         << wrench.force_world.y() << "," << wrench.force_world.z() << "," << wrench.force_norm
                         << "," << wrench.torque_world.x() << "," << wrench.torque_world.y() << ","
                         << wrench.torque_world.z() << "," << wrench.torque_norm << "," << ideal_wrench.force.x()
                         << "," << ideal_wrench.force.y() << "," << ideal_wrench.force.z() << ","
                         << ideal_wrench.force.Length() << "," << ideal_wrench.torque_about_body_ref.x() << ","
                         << ideal_wrench.torque_about_body_ref.y() << ","
                         << ideal_wrench.torque_about_body_ref.z() << ","
                         << ideal_wrench.torque_about_body_ref.Length() << "," << force_error.x() << ","
                         << force_error.y() << "," << force_error.z() << "," << force_error.Length() << ","
                         << torque_error.x() << "," << torque_error.y() << "," << torque_error.z() << ","
                         << torque_error.Length() << "," << force_cosine << "," << torque_cosine << ","
                         << wrench.active_contacts
                         << ",ideal reference loaded from rev_joint_clearance_ideal_revolute/joint_reaction_wrench.csv "
                            "when available\n";
        int open_positive_rows = 0;
        double open_positive_force_sum = 0.0;
        double open_positive_force_max = 0.0;
        for (const auto& state : contact_states) {
            if (state.sample.contact_id >= 0 && state.sample.gap > contact_law.bpen && state.lambda_force > 0.0) {
                open_positive_rows++;
                open_positive_force_sum += state.lambda_force;
                open_positive_force_max = std::max(open_positive_force_max, state.lambda_force);
            }
        }
        const auto manifold_diag = manifold_manager->GetDiagnostics();
        manifold_diagnostics << time << "," << manifold_diag.candidate_count << "," << manifold_diag.active_count
                             << "," << manifold_diag.patch_count << "," << manifold_diag.reused_patch_count << ","
                             << manifold_diag.new_patch_count << "," << manifold_diag.released_sample_count << ","
                             << manifold_diag.active_area << "," << open_positive_rows << ","
                             << open_positive_force_sum << "," << open_positive_force_max << "\n";
        const auto patch_moment_rows =
            ComputeRevClearancePatchMomentRows(contact_states, body3, use_raw_multiplier_as_effective_force);
        for (const auto& row : patch_moment_rows) {
            if (row.patch_id != -1) {
                continue;
            }
            const double effective_force_norm = row.effective_force.Length();
            const double effective_torque_norm = row.effective_torque.Length();
            const double resultant_line_offset =
                effective_force_norm > 1.0e-30 ? effective_torque_norm / effective_force_norm : 0.0;
            const double closure_ratio =
                effective_force_norm * row.max_moment_arm > 1.0e-30 ?
                    effective_torque_norm / (effective_force_norm * row.max_moment_arm) :
                    0.0;
            max_patch_effective_force_norm = std::max(max_patch_effective_force_norm, effective_force_norm);
            max_patch_effective_torque_norm = std::max(max_patch_effective_torque_norm, effective_torque_norm);
            max_patch_resultant_line_offset = std::max(max_patch_resultant_line_offset, resultant_line_offset);
            max_patch_closure_ratio = std::max(max_patch_closure_ratio, closure_ratio);
            max_patch_active_samples = std::max(max_patch_active_samples, row.active_samples);
        }
        WriteRevClearancePatchMomentRows(patch_moment_diagnostics, time, backend_label, patch_moment_rows);
        const auto body3_q = body3->GetRot();
        const ChVector3d body3_w = body3->GetAngVelParent();

        reference_timeseries << time << "," << ref2.pos.x() << "," << ref2.pos.y() << "," << ref2.pos.z() << ","
                             << body2->GetPos().x() << "," << body2->GetPos().y() << "," << body2->GetPos().z()
                             << "," << body2_y_errors.back() << "," << body2_z_errors.back() << ","
                             << ref3.pos.x() << "," << ref3.pos.y() << "," << ref3.pos.z() << ","
                             << body3->GetPos().x() << "," << body3->GetPos().y() << "," << body3->GetPos().z()
                             << "," << body3_y_errors.back() << "," << body3_z_errors.back() << ","
                             << diag.min_gap << "," << diag.max_penetration << "," << (diag.success ? 1 : 0)
                             << "\n";
        kinematic_diagnostics
            << time << "," << ref3.pos.x() << "," << ref3.pos.y() << "," << ref3.pos.z() << ","
            << body3->GetPos().x() << "," << body3->GetPos().y() << "," << body3->GetPos().z() << ","
            << body3_q.e0() << "," << body3_q.e1() << "," << body3_q.e2() << "," << body3_q.e3() << ","
            << body3_w.x() << "," << body3_w.y() << "," << body3_w.z() << ","
            << rel_ref.x() << "," << rel_ref.y() << "," << rel_ref.z() << ","
            << rel_sim.x() << "," << rel_sim.y() << "," << rel_sim.z() << ","
            << rel_ref.Length() << "," << rel_sim.Length() << "," << rel_angle_error << ","
            << rel_angle_error * 180.0 / kPi << "," << rel_axis_ref.x() << "," << rel_axis_ref.y() << ","
            << rel_axis_ref.z() << "," << rel_axis_sim.x() << "," << rel_axis_sim.y() << "," << rel_axis_sim.z()
            << "," << action_marker_sim.x() << "," << action_marker_sim.y() << "," << action_marker_sim.z() << ","
            << action_marker_ref_proxy.x() << "," << action_marker_ref_proxy.y() << "," << action_marker_ref_proxy.z()
            << "," << wrench.center_of_pressure.x() << "," << wrench.center_of_pressure.y() << ","
            << wrench.center_of_pressure.z() << "," << wrench.force_world.x() << "," << wrench.force_world.y()
            << "," << wrench.force_world.z() << "," << wrench.force_norm << "," << wrench.torque_world.x() << ","
            << wrench.torque_world.y() << "," << wrench.torque_world.z() << "," << wrench.torque_norm << ","
            << wrench.active_contacts << "," << diag.min_gap << "," << diag.max_penetration << ","
            << (diag.success ? 1 : 0)
            << ",RecurDyn attitude and action-marker columns are swing proxies inferred from Body2-Body3 relative "
               "position because the provided reference CSV has no orientation or force-moment output\n";

        const auto write_body_row = [&](const std::shared_ptr<ChBodyAuxRef>& body,
                                        const std::string& name,
                                        int contact_id,
                                        double gap,
                                        double lambda,
                                        double weight,
                                        double lambda_force,
                                        double ncp,
                                        double comp) {
            const auto q = body->GetRot();
            const ChVector3d w = body->GetAngVelParent();
            trajectory << time << "," << name << "," << body->GetPos().x() << "," << body->GetPos().y() << ","
                       << body->GetPos().z() << "," << q.e0() << "," << q.e1() << "," << q.e2() << "," << q.e3()
                       << "," << body->GetPosDt().x() << "," << body->GetPosDt().y() << "," << body->GetPosDt().z()
                       << "," << w.x() << "," << w.y() << "," << w.z() << "," << contact_id << "," << gap << ","
                       << lambda << "," << weight << "," << lambda_force << "," << std::max(0.0, -gap) << ","
                       << ncp << "," << comp << ",0,0," << (diag.success ? 1 : 0) << "\n";
        };

        write_body_row(body2, "body2", -1, diag.min_gap, 0.0, 0.0, 0.0, 0.0, 0.0);
        if (diag.candidates.empty()) {
            write_body_row(body3, "body3", -1, diag.min_gap, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        for (size_t c = 0; c < diag.candidates.size(); c++) {
            const auto& candidate = diag.candidates[c];
            const double lambda = c < diag.lambdas.size() ? diag.lambdas[c] : 0.0;
            const double lambda_force = c < diag.lambda_forces.size() ? diag.lambda_forces[c] : lambda * candidate.weight;
            const double ncp = c < diag.ncp_residuals.size() ? diag.ncp_residuals[c] : 0.0;
            const double comp = c < diag.complementarity_errors.size() ? diag.complementarity_errors[c] : 0.0;
            write_body_row(
                body3, "body3", candidate.sample_id, candidate.gap, lambda, candidate.weight, lambda_force, ncp, comp);
        }
    }

    const double body2_y_rmse = RmsError(body2_y_errors);
    const double body2_z_rmse = RmsError(body2_z_errors);
    const double body3_y_rmse = RmsError(body3_y_errors);
    const double body3_z_rmse = RmsError(body3_z_errors);
    const double rel_angle_rmse = RmsError(relative_vector_angle_errors);
    const double rel_angle_max = MaxAbsError(relative_vector_angle_errors);
    const double wrench_force_error_rmse = RmsError(wrench_force_error_norms);
    const double wrench_force_error_max = MaxAbsError(wrench_force_error_norms);
    const double wrench_torque_error_rmse = RmsError(wrench_torque_error_norms);
    const double wrench_torque_error_max = MaxAbsError(wrench_torque_error_norms);
    const double wrench_force_cosine_mean = MeanValue(wrench_force_cosines);
    const double wrench_torque_cosine_mean = MeanValue(wrench_torque_cosines);
    comparison << config.case_name << ",body2_y_rmse_vs_recurdyn,0," << body2_y_rmse << "," << body2_y_rmse
               << ",0,reference assets/rev_joint_clearance/data/body2.csv Y:Pos_TY\n";
    comparison << config.case_name << ",body2_z_rmse_vs_recurdyn,0," << body2_z_rmse << "," << body2_z_rmse
               << ",0,reference assets/rev_joint_clearance/data/body2.csv Y:Pos_TZ\n";
    comparison << config.case_name << ",body3_y_rmse_vs_recurdyn,0," << body3_y_rmse << "," << body3_y_rmse
               << ",0,reference assets/rev_joint_clearance/data/body3.csv Y:Pos_TY\n";
    comparison << config.case_name << ",body3_z_rmse_vs_recurdyn,0," << body3_z_rmse << "," << body3_z_rmse
               << ",0,reference assets/rev_joint_clearance/data/body3.csv Y:Pos_TZ\n";
    comparison << config.case_name << ",ggeomcontact_K," << contact_law.stiffness << "," << contact_law.stiffness
               << ",0,0,parsed from RecurDyn GGEOMCONTACT"
               << (options.use_recurdyn_ggeomcontact_law ? "; applied by calibration backend\n" :
                                                   "; hard SDF-NCP does not apply penalty K\n");
    comparison << config.case_name << ",ggeomcontact_C," << contact_law.damping << "," << contact_law.damping
               << ",0,0,parsed from RecurDyn GGEOMCONTACT"
               << (options.use_recurdyn_ggeomcontact_law ? "; applied by calibration backend\n" :
                                                   "; hard SDF-NCP does not apply damping C\n");
    comparison << config.case_name << ",ggeomcontact_bpen," << contact_law.bpen << "," << contact_law.bpen
               << ",0,0,parsed from RecurDyn GGEOMCONTACT"
               << (options.use_recurdyn_ggeomcontact_law ? "; applied by calibration backend\n" :
                                                   "; reported for contact-law audit\n");
    comparison << config.case_name << ",ggeomcontact_nomaxcp_reference,10," << config.max_contacts << ","
               << std::abs(10.0 - static_cast<double>(config.max_contacts)) << ",0,"
               << "SDF-NCP keeps full active patch capacity instead of compressing to NOMAXCP\n";
    comparison << config.case_name << ",sdf_ncp_active_gap_on,0," << contact_law.bpen << "," << contact_law.bpen
               << ",0,descriptor active-set starts constraints at RecurDyn BPEN-scale geometric gap\n";
    comparison << config.case_name << ",sdf_ncp_active_gap_off,0," << config.max_activation_band << ","
               << config.max_activation_band
               << ",0,persistent contact_id hysteresis band; retained only while normal gap rate remains closing\n";
    comparison << config.case_name << ",sdf_ncp_pressure_compliance,0," << config.pressure_compliance << ","
               << config.pressure_compliance
               << ",0,diagonal local pressure/intensity compliance applied as force_scale*compliance CFM\n";
    comparison << config.case_name << ",wrench_force_error_rmse_vs_ideal_revolute,0," << wrench_force_error_rmse
               << "," << wrench_force_error_rmse
               << ",0,contact wrench force error against Chrono ideal revolute joint reaction wrench\n";
    comparison << config.case_name << ",wrench_torque_error_rmse_vs_ideal_revolute,0," << wrench_torque_error_rmse
               << "," << wrench_torque_error_rmse
               << ",0,contact wrench torque error against Chrono ideal revolute joint reaction wrench\n";
    comparison << config.case_name << ",wrench_force_cosine_mean_vs_ideal_revolute,1," << wrench_force_cosine_mean
               << "," << std::abs(1.0 - wrench_force_cosine_mean)
               << ",0,mean cosine between SDF contact force and ideal revolute reaction force\n";
    comparison << config.case_name << ",wrench_torque_cosine_mean_vs_ideal_revolute,1," << wrench_torque_cosine_mean
               << "," << std::abs(1.0 - wrench_torque_cosine_mean)
               << ",0,mean cosine between SDF contact torque and ideal revolute reaction torque\n";
    comparison << config.case_name << ",sdf_ncp_patch_pressure_aggregation,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.patch_pressure_aggregation ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.patch_pressure_aggregation ? 1 : 0)
               << ",0,hard SDF-NCP path aggregates connected active samples into one pressure unknown per patch; "
                  "GGEOMCONTACT calibration keeps per-sample penalty forces\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_patch_projection,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_projection ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_projection ? 1 : 0)
               << ",0,descriptor-level unilateral rows remain in Chrono global solver; the PSOR projection hook "
                  "solves a small nonnegative patch Delassus block and increments descriptor states consistently\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_patch_delassus_block,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ? 1 : 0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ? 1 : 0)
               << ",0,descriptor-level patch block solve uses W=J M^-1 J^T + CFM inside PSOR group projection; "
                  "sample pressure rows, joint constraints, drivers, and body velocities are iterated in one "
                  "Chrono descriptor solve\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_patch_laplacian_compliance_scale,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ?
                       options.descriptor_patch_laplacian_compliance_scale :
                       0.0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ?
                       options.descriptor_patch_laplacian_compliance_scale :
                       0.0)
               << ",0,optional patch graph Laplacian pressure regularization added to the descriptor block matrix\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_velocity_level_ncp,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_velocity_level_ncp ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_velocity_level_ncp ? 1 : 0)
               << ",0,descriptor RHS uses velocity-level normal approach residual v_n + beta min(g,0)/dt instead of "
                  "raw geometric gap stabilization\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_velocity_baumgarte,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_velocity_level_ncp ?
                       options.descriptor_velocity_baumgarte :
                       0.0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_velocity_level_ncp ?
                       options.descriptor_velocity_baumgarte :
                       0.0)
               << ",0,dimensionless beta in beta min(g,0)/dt velocity-level stabilization\n";
    comparison << config.case_name << ",sdf_ncp_descriptor_patch_wrench_closure_strength,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ?
                       options.descriptor_patch_wrench_closure_strength :
                       0.0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.descriptor_patch_delassus_solve ?
                       options.descriptor_patch_wrench_closure_strength :
                       0.0)
               << ",0,optional torque Gram regularization inside the descriptor patch block; it keeps sample pressure "
                  "nonnegative while discouraging large free contact torque on the contacted body\n";
    comparison << config.case_name << ",sdf_ncp_local_patch_pressure_solve,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_pressure_solve ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_pressure_solve ? 1 : 0)
               << ",0,experimental force backend keeps sample pressure unknowns inside each patch and couples them "
                  "with an area-weighted graph Laplacian local FB solve\n";
    comparison << config.case_name << ",sdf_ncp_local_patch_delassus,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_delassus ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_delassus ? 1 : 0)
               << ",0,velocity-level local patch solve uses W=J M^-1 J^T effective-mass coupling before applying "
                  "area-integrated sample forces\n";
    comparison << config.case_name << ",sdf_ncp_local_patch_wrench_closure,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_wrench_closure ? 1 : 0) << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_wrench_closure ? 1 : 0)
               << ",0,local patch pressure solve preserves the FB-computed resultant force and redistributes "
                  "nonnegative sample pressures to reduce free torque about the contact body reference point\n";
    comparison << config.case_name << ",sdf_ncp_local_patch_global_wrench_demand,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_global_wrench_demand ? 1 : 0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_global_wrench_demand ? 1 : 0)
               << ",0,diagnostic mode uses Chrono ideal revolute reaction wrench as the patch pressure demand; "
                  "it is not a production SDF-NCP closure\n";
    comparison << config.case_name << ",sdf_ncp_local_patch_integrated_wrench_closure,0,"
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_integrated_wrench_closure ? 1 : 0)
               << ","
               << (!options.use_recurdyn_ggeomcontact_law && options.local_patch_integrated_wrench_closure ? 1 : 0)
               << ",0,augmented local patch block includes FB complementarity plus force preservation and body-reference torque residuals\n";
    comparison << config.case_name << ",patch_effective_force_norm_max,0," << max_patch_effective_force_norm << ","
               << max_patch_effective_force_norm
               << ",0,global active patch resultant force norm from contact_patch_moment_diagnostics.csv\n";
    comparison << config.case_name << ",patch_effective_torque_norm_max,0," << max_patch_effective_torque_norm << ","
               << max_patch_effective_torque_norm
               << ",0,global active patch torque norm about Body3 reference from "
                  "contact_patch_moment_diagnostics.csv\n";
    comparison << config.case_name << ",patch_resultant_line_offset_m_max,0," << max_patch_resultant_line_offset
               << "," << max_patch_resultant_line_offset
               << ",0,|sum r cross f| / |sum f|; fixed regression metric for patch torque closure\n";
    comparison << config.case_name << ",patch_closure_ratio_max,0," << max_patch_closure_ratio << ","
               << max_patch_closure_ratio
               << ",0,|sum r cross f| / (|sum f| max|r|); fixed regression metric for patch torque closure\n";

    std::ofstream patch_moment_summary(out_dir / "contact_patch_moment_summary.csv");
    patch_moment_summary << std::setprecision(17);
    patch_moment_summary
        << "case_name,method,max_effective_force_norm,max_effective_torque_norm,max_resultant_line_offset_m,"
           "max_closure_ratio,max_active_samples,notes\n";
    patch_moment_summary
        << config.case_name << "," << stats.method << "," << max_patch_effective_force_norm << ","
        << max_patch_effective_torque_norm << "," << max_patch_resultant_line_offset << ","
        << max_patch_closure_ratio << "," << max_patch_active_samples
        << ",Regression summary for global patch_id=-1 rows in contact_patch_moment_diagnostics.csv; "
           "SDF-NCP descriptor lambda_n is local intensity and lambda_force is area/scale-integrated sample force\n";

    std::ofstream wrench_summary(out_dir / "wrench_alignment_summary.csv");
    wrench_summary << std::setprecision(17);
    wrench_summary
        << "case_name,method,num_samples,force_error_rmse,force_error_max,torque_error_rmse,torque_error_max,"
           "force_cosine_mean,torque_cosine_mean,notes\n";
    wrench_summary
        << config.case_name << "," << stats.method << "," << wrench_force_error_norms.size() << ","
        << wrench_force_error_rmse << "," << wrench_force_error_max << "," << wrench_torque_error_rmse << ","
        << wrench_torque_error_max << "," << wrench_force_cosine_mean << "," << wrench_torque_cosine_mean
        << ",Fixed regression metrics comparing the integrated SDF contact wrench on Body3 against "
           "rev_joint_clearance_ideal_revolute/joint_reaction_wrench.csv; the reference is diagnostic only and is "
           "not used by the SDF-NCP solver\n";

    std::ofstream reference_summary(out_dir / "recurdyn_reference_comparison_summary.csv");
    reference_summary << std::setprecision(17);
    reference_summary << "case_name,quantity,time_start,time_end,num_samples,rmse,max_abs_error,final_reference,"
                         "final_sdf_ncp,final_abs_error,notes\n";
    const auto final_ref2 = InterpolateRevClearanceReference(body2_reference, config.t_end);
    const auto final_ref3 = InterpolateRevClearanceReference(body3_reference, config.t_end);
    reference_summary << config.case_name << ",body2_y,0," << config.t_end << "," << body2_y_errors.size() << ","
                      << body2_y_rmse << "," << MaxAbsError(body2_y_errors) << "," << final_ref2.pos.y() << ","
                      << body2->GetPos().y() << "," << std::abs(body2->GetPos().y() - final_ref2.pos.y())
                      << ",Y:Pos_TY-Body2 compared to Chrono Body2 position y\n";
    reference_summary << config.case_name << ",body2_z,0," << config.t_end << "," << body2_z_errors.size() << ","
                      << body2_z_rmse << "," << MaxAbsError(body2_z_errors) << "," << final_ref2.pos.z() << ","
                      << body2->GetPos().z() << "," << std::abs(body2->GetPos().z() - final_ref2.pos.z())
                      << ",Y:Pos_TZ-Body2 compared to Chrono Body2 position z\n";
    reference_summary << config.case_name << ",body3_y,0," << config.t_end << "," << body3_y_errors.size() << ","
                      << body3_y_rmse << "," << MaxAbsError(body3_y_errors) << "," << final_ref3.pos.y() << ","
                      << body3->GetPos().y() << "," << std::abs(body3->GetPos().y() - final_ref3.pos.y())
                      << ",Y:Pos_TY-Body3 compared to Chrono Body3 position y\n";
    reference_summary << config.case_name << ",body3_z,0," << config.t_end << "," << body3_z_errors.size() << ","
                      << body3_z_rmse << "," << MaxAbsError(body3_z_errors) << "," << final_ref3.pos.z() << ","
                      << body3->GetPos().z() << "," << std::abs(body3->GetPos().z() - final_ref3.pos.z())
                      << ",Y:Pos_TZ-Body3 compared to Chrono Body3 position z\n";
    const double final_rel_angle =
        SafeVectorAngle(final_ref2.pos - final_ref3.pos, body2->GetPos() - body3->GetPos());
    reference_summary << config.case_name << ",body2_body3_relative_vector_angle_rad,0," << config.t_end << ","
                      << relative_vector_angle_errors.size() << "," << rel_angle_rmse << "," << rel_angle_max
                      << ",0," << final_rel_angle << "," << final_rel_angle
                      << ",orientation phase proxy from Body2-Body3 relative vector; RecurDyn reference has no direct "
                         "Body3 quaternion output\n";

    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);
    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "recurdyn_reference_comparison_timeseries.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "recurdyn_chrono_kinematic_diagnostics.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "contact_patch_moment_diagnostics.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "contact_patch_moment_summary.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "contact_manifold_diagnostics.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "wrench_alignment.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "wrench_alignment_summary.csv").string() << std::endl;

    const double success_rate = stats.sum_success / static_cast<double>(std::max(1, stats.samples));
    const bool patch_closure_bounded =
        options.use_recurdyn_ggeomcontact_law ||
        (std::isfinite(max_patch_resultant_line_offset) && std::isfinite(max_patch_closure_ratio) &&
         max_patch_resultant_line_offset < 1.2 && max_patch_closure_ratio < 1.05);
    const bool wrench_alignment_bounded =
        wrench_force_error_norms.empty() ||
        (std::isfinite(wrench_force_error_rmse) && std::isfinite(wrench_torque_error_rmse) &&
         wrench_force_error_rmse < 5.0e5 && wrench_torque_error_rmse < 2.0e5);
    return std::isfinite(body3_y_rmse) && stats.max_penetration < 5.0e-2 && success_rate > 0.90 &&
                   patch_closure_bounded && wrench_alignment_bounded ?
               0 :
               1;
}

int RunRevJointClearanceCase(const std::string& output_case_name = "rev_joint_clearance",
                             double t_end_override = -1.0,
                             double dt_override = -1.0,
                             bool use_recurdyn_ggeomcontact_law = false) {
    RevClearanceRunOptions options;
    options.case_name = output_case_name;
    options.t_end_override = t_end_override;
    options.dt_override = dt_override;
    options.use_recurdyn_ggeomcontact_law = use_recurdyn_ggeomcontact_law;
    return RunRevJointClearanceCase(options);
}

int RunRevJointClearanceIdealRevoluteCase(const std::string& output_case_name = "rev_joint_clearance_ideal_revolute",
                                          double t_end_override = 3.0,
                                          double dt_override = 0.001) {
    const auto start = std::chrono::steady_clock::now();
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "rev_joint_clearance";
    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / output_case_name;
    std::filesystem::create_directories(out_dir);

    RevClearanceConfig config;
    const std::string model_json = ReadTextFile(asset_dir / "rev_joint_clearance_model.json");
    config.dt = ExtractJsonDouble(model_json, "step_size", config.dt);
    config.t_end = ExtractJsonDouble(model_json, "total_time", config.t_end);
    if (t_end_override > 0.0) {
        config.t_end = t_end_override;
    }
    if (dt_override > 0.0) {
        config.dt = dt_override;
    }

    const RmdModel rmd = LoadRecurDynRmdModel(asset_dir / "rev_clearance_joint.rmd");
    WriteRmdMappingAudit(out_dir / "rmd_mapping.csv", rmd);
    const RmdPart& body1_part = RequirePartByName(rmd, "Body1");
    const RmdPart& body2_part = RequirePartByName(rmd, "Body2");
    const RmdPart& body3_part = RequirePartByName(rmd, "Body3");
    const RmdMarker& body1_cm = RequireMarkerById(rmd, body1_part.cm_marker_id);
    const RmdMarker& body2_cm = RequireMarkerById(rmd, body2_part.cm_marker_id);
    const RmdMarker& body3_cm = RequireMarkerById(rmd, body3_part.cm_marker_id);
    const RmdJoint& fixed1 = RequireJointByName(rmd, "Fixed1");
    const RmdJoint& fixed2 = RequireJointByName(rmd, "Fixed2");

    const auto body2_ideal_reference = LoadRevClearanceReference(asset_dir / "data" / "body2_ideal.csv");
    if (body2_ideal_reference.empty()) {
        throw std::runtime_error("rev_joint_clearance ideal Body2 reference CSV is empty.");
    }

    RevClearanceRunOptions options;
    options.timestepper_type = ChTimestepper::Type::HHT;
    options.timestepper_label = "hht_alpha_minus_0p2_fixed_step";
    options.use_step_control = false;

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(rmd.gravity);
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(300);
    sys.GetSolver()->AsIterative()->SetTolerance(1.0e-8);
    sys.SetMaxPenetrationRecoverySpeed(0.5);
    ConfigureRevClearanceTimestepper(sys, options);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.AddBody(ground);

    auto body1 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body1_part, body1_cm, body1);
    sys.AddBody(body1);

    auto body2 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body2_part, body2_cm, body2);
    sys.AddBody(body2);

    auto body3 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body3_part, body3_cm, body3);
    sys.AddBody(body3);

    auto fixed1_link = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    fixed1_link->SetName("Fixed1");
    fixed1_link->Initialize(body1, ground, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, fixed1.i_marker_id)));
    sys.AddLink(fixed1_link);

    auto fixed2_link = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    fixed2_link->SetName("Fixed2");
    fixed2_link->Initialize(body3, body2, RmdMarkerFrameAbs(rmd, RequireMarkerById(rmd, fixed2.i_marker_id)));
    sys.AddLink(fixed2_link);

    auto revolute = chrono_types::make_shared<ChLinkRevolute>();
    revolute->SetName("ChronoIdealRevJoint_Body3_to_Body1_global_X");
    const ChFrame<> rev_frame(ChVector3d(0, 0, 0), ToChronoMatrix(RotYMatrix(kPi / 2.0)));
    revolute->Initialize(body3, body1, rev_frame);
    sys.AddLink(revolute);

    std::ofstream trajectory(out_dir / "trajectory.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    std::ofstream timeseries(out_dir / "ideal_revolute_comparison_timeseries.csv");
    std::ofstream joint_wrench(out_dir / "joint_reaction_wrench.csv");
    trajectory << std::setprecision(17);
    comparison << std::setprecision(17);
    timeseries << std::setprecision(17);
    joint_wrench << std::setprecision(17);
    WriteTrajectoryHeader(trajectory);
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";
    timeseries << "time,body2_x_ideal,body2_y_ideal,body2_z_ideal,body2_x_chrono,body2_y_chrono,"
                  "body2_z_chrono,body2_y_error,body2_z_error,body3_x_chrono,body3_y_chrono,body3_z_chrono,"
                  "body3_q0_chrono,body3_q1_chrono,body3_q2_chrono,body3_q3_chrono,body3_wx_chrono,"
                  "joint_force_x,joint_force_y,joint_force_z,joint_force_norm,joint_torque_about_body3_x,"
                  "joint_torque_about_body3_y,joint_torque_about_body3_z,joint_torque_about_body3_norm,success\n";
    joint_wrench
        << "time,body,frame,force_x,force_y,force_z,force_norm,torque_at_joint_x,torque_at_joint_y,"
           "torque_at_joint_z,torque_at_joint_norm,torque_about_body_ref_x,torque_about_body_ref_y,"
           "torque_about_body_ref_z,torque_about_body_ref_norm,joint_x,joint_y,joint_z,body_ref_x,body_ref_y,"
           "body_ref_z,success\n";

    const int steps = static_cast<int>(std::round(config.t_end / config.dt));
    BenchmarkStats stats;
    stats.case_name = output_case_name;
    stats.method = "chrono_ideal_revolute";
    stats.asset = "assets/rev_joint_clearance";
    stats.dt = config.dt;
    stats.t_end = config.t_end;
    stats.num_steps = steps;
    stats.notes =
        "Chrono diagnostic ideal revolute joint: same RMD Body1/Body2/Body3 mass-marker-inertia mapping and Body2-Body3 "
        "fixed joint, no contact, Body3 revolute to Body1 about global X axis; compared to body2_ideal.csv";

    std::vector<double> body2_y_errors;
    std::vector<double> body2_z_errors;
    bool finite_all = true;
    for (int i = 0; i <= steps; i++) {
        if (i > 0) {
            sys.DoStepDynamics(config.dt);
        } else {
            sys.Setup();
            sys.Update(UpdateFlags::UPDATE_ALL & ~UpdateFlags::VISUAL_ASSETS);
        }

        const double time = sys.GetChTime();
        const auto ref = InterpolateRevClearanceReference(body2_ideal_reference, time);
        body2_y_errors.push_back(body2->GetPos().y() - ref.pos.y());
        body2_z_errors.push_back(body2->GetPos().z() - ref.pos.z());
        finite_all = finite_all && std::isfinite(body2->GetPos().y()) && std::isfinite(body3->GetPos().y());
        const auto q2 = body2->GetRot();
        const auto q3 = body3->GetRot();
        const ChVector3d w2 = body2->GetAngVelParent();
        const ChVector3d w3 = body3->GetAngVelParent();
        const auto reaction1 = revolute->GetReaction1();
        const auto reaction2 = revolute->GetReaction2();
        const ChVector3d joint_pos_body3 = revolute->GetFrame1Abs().GetPos();
        const ChVector3d joint_pos_body1 = revolute->GetFrame2Abs().GetPos();
        const ChVector3d reaction1_force_world =
            revolute->GetFrame1Abs().TransformDirectionLocalToParent(reaction1.force);
        const ChVector3d reaction1_torque_at_joint_world =
            revolute->GetFrame1Abs().TransformDirectionLocalToParent(reaction1.torque);
        const ChVector3d reaction1_torque_about_body3_world =
            reaction1_torque_at_joint_world + (joint_pos_body3 - body3->GetPos()).Cross(reaction1_force_world);
        const ChVector3d reaction2_force_world =
            revolute->GetFrame2Abs().TransformDirectionLocalToParent(reaction2.force);
        const ChVector3d reaction2_torque_at_joint_world =
            revolute->GetFrame2Abs().TransformDirectionLocalToParent(reaction2.torque);
        const ChVector3d reaction2_torque_about_body1_world =
            reaction2_torque_at_joint_world + (joint_pos_body1 - body1->GetPos()).Cross(reaction2_force_world);
        trajectory << time << ",body2," << body2->GetPos().x() << "," << body2->GetPos().y() << ","
                   << body2->GetPos().z() << "," << q2.e0() << "," << q2.e1() << "," << q2.e2() << "," << q2.e3()
                   << "," << body2->GetPosDt().x() << "," << body2->GetPosDt().y() << "," << body2->GetPosDt().z()
                   << "," << w2.x() << "," << w2.y() << "," << w2.z()
                   << ",-1,0,0,0,0,0,0,0,0,0," << (finite_all ? 1 : 0) << "\n";
        trajectory << time << ",body3," << body3->GetPos().x() << "," << body3->GetPos().y() << ","
                   << body3->GetPos().z() << "," << q3.e0() << "," << q3.e1() << "," << q3.e2() << "," << q3.e3()
                   << "," << body3->GetPosDt().x() << "," << body3->GetPosDt().y() << "," << body3->GetPosDt().z()
                   << "," << w3.x() << "," << w3.y() << "," << w3.z()
                   << ",-1,0,0,0,0,0,0,0,0,0," << (finite_all ? 1 : 0) << "\n";
        timeseries << time << "," << ref.pos.x() << "," << ref.pos.y() << "," << ref.pos.z() << ","
                   << body2->GetPos().x() << "," << body2->GetPos().y() << "," << body2->GetPos().z() << ","
                   << body2_y_errors.back() << "," << body2_z_errors.back() << "," << body3->GetPos().x() << ","
                   << body3->GetPos().y() << "," << body3->GetPos().z() << "," << q3.e0() << "," << q3.e1() << ","
                   << q3.e2() << "," << q3.e3() << "," << w3.x() << "," << reaction1_force_world.x() << ","
                   << reaction1_force_world.y() << "," << reaction1_force_world.z() << ","
                   << reaction1_force_world.Length() << "," << reaction1_torque_about_body3_world.x() << ","
                   << reaction1_torque_about_body3_world.y() << "," << reaction1_torque_about_body3_world.z() << ","
                   << reaction1_torque_about_body3_world.Length() << "," << (finite_all ? 1 : 0) << "\n";
        joint_wrench << time << ",body3,frame1," << reaction1_force_world.x() << "," << reaction1_force_world.y()
                     << "," << reaction1_force_world.z() << "," << reaction1_force_world.Length() << ","
                     << reaction1_torque_at_joint_world.x() << "," << reaction1_torque_at_joint_world.y() << ","
                     << reaction1_torque_at_joint_world.z() << "," << reaction1_torque_at_joint_world.Length() << ","
                     << reaction1_torque_about_body3_world.x() << "," << reaction1_torque_about_body3_world.y()
                     << "," << reaction1_torque_about_body3_world.z() << ","
                     << reaction1_torque_about_body3_world.Length() << "," << joint_pos_body3.x() << ","
                     << joint_pos_body3.y() << "," << joint_pos_body3.z() << "," << body3->GetPos().x() << ","
                     << body3->GetPos().y() << "," << body3->GetPos().z() << "," << (finite_all ? 1 : 0) << "\n";
        joint_wrench << time << ",body1,frame2," << reaction2_force_world.x() << "," << reaction2_force_world.y()
                     << "," << reaction2_force_world.z() << "," << reaction2_force_world.Length() << ","
                     << reaction2_torque_at_joint_world.x() << "," << reaction2_torque_at_joint_world.y() << ","
                     << reaction2_torque_at_joint_world.z() << "," << reaction2_torque_at_joint_world.Length() << ","
                     << reaction2_torque_about_body1_world.x() << "," << reaction2_torque_about_body1_world.y()
                     << "," << reaction2_torque_about_body1_world.z() << ","
                     << reaction2_torque_about_body1_world.Length() << "," << joint_pos_body1.x() << ","
                     << joint_pos_body1.y() << "," << joint_pos_body1.z() << "," << body1->GetPos().x() << ","
                     << body1->GetPos().y() << "," << body1->GetPos().z() << "," << (finite_all ? 1 : 0) << "\n";
    }

    const double body2_y_rmse = RmsError(body2_y_errors);
    const double body2_z_rmse = RmsError(body2_z_errors);
    comparison << output_case_name << ",body2_y_rmse_vs_recurdyn_ideal,0," << body2_y_rmse << "," << body2_y_rmse
               << ",0,reference assets/rev_joint_clearance/data/body2_ideal.csv Y:Pos_TY\n";
    comparison << output_case_name << ",body2_z_rmse_vs_recurdyn_ideal,0," << body2_z_rmse << "," << body2_z_rmse
               << ",0,reference assets/rev_joint_clearance/data/body2_ideal.csv Y:Pos_TZ\n";
    comparison << output_case_name << ",joint_axis,0,0,0,0,Chrono ChLinkRevolute z-axis aligned with global X using "
                  "RotY(pi/2) joint frame\n";

    stats.samples = steps + 1;
    stats.sum_success = finite_all ? static_cast<double>(steps + 1) : 0.0;
    stats.runtime_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    WriteSummary(out_dir / "summary.csv", stats);

    std::cout << "Wrote " << (out_dir / "trajectory.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "ideal_revolute_comparison_timeseries.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "joint_reaction_wrench.csv").string() << std::endl;
    return finite_all && body2_y_rmse < 1.0e-2 && body2_z_rmse < 1.0e-2 ? 0 : 1;
}

struct RevClearanceIdealFrame {
    double time = 0.0;
    ChVector3d body3_pos = ChVector3d(0, 0, 0);
    chrono::ChQuaternion<> body3_rot = chrono::ChQuaternion<>(1, 0, 0, 0);
    ChVector3d body3_vel = ChVector3d(0, 0, 0);
    ChVector3d body3_omega = ChVector3d(0, 0, 0);
};

std::vector<RevClearanceIdealFrame> LoadRevClearanceIdealFrames(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        return {};
    }
    std::string header;
    if (!std::getline(in, header)) {
        return {};
    }

    std::vector<RevClearanceIdealFrame> frames;
    std::string line;
    double last_time = std::numeric_limits<double>::quiet_NaN();
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const auto cols = SplitCsvLine(line);
        if (cols.size() < 15 || cols[1] != "body3") {
            continue;
        }
        const double time = std::stod(cols[0]);
        if (std::isfinite(last_time) && std::abs(time - last_time) < 1.0e-12) {
            continue;
        }
        last_time = time;
        RevClearanceIdealFrame frame;
        frame.time = time;
        frame.body3_pos = ChVector3d(std::stod(cols[2]), std::stod(cols[3]), std::stod(cols[4]));
        frame.body3_rot =
            chrono::ChQuaternion<>(std::stod(cols[5]), std::stod(cols[6]), std::stod(cols[7]), std::stod(cols[8]));
        frame.body3_vel = ChVector3d(std::stod(cols[9]), std::stod(cols[10]), std::stod(cols[11]));
        frame.body3_omega = ChVector3d(std::stod(cols[12]), std::stod(cols[13]), std::stod(cols[14]));
        frames.push_back(frame);
    }
    return frames;
}

double RevClearanceContactForceScale(const ChSdfNcpDescriptorContact& sample) {
    if (sample.area > 1.0e-30 && std::isfinite(sample.area)) {
        return sample.area;
    }
    if (sample.weight > 1.0e-30 && std::isfinite(sample.weight)) {
        return sample.weight;
    }
    return 1.0;
}

ChVector3d RevClearancePointVelocityAbs(const ChBody& body, const ChVector3d& point_abs) {
    return body.GetPosDt() + body.GetAngVelParent().Cross(point_abs - body.GetPos());
}

struct RevClearancePressureColumn {
    ChVector3d force = ChVector3d(0, 0, 0);
    ChVector3d torque = ChVector3d(0, 0, 0);
    double scale = 1.0;
};

RevClearancePressureColumn RevClearanceBody3PressureColumn(const ChSdfNcpDescriptorContact& sample,
                                                           const std::shared_ptr<ChBodyAuxRef>& body3) {
    RevClearancePressureColumn column;
    column.scale = RevClearanceContactForceScale(sample);
    if (sample.use_custom_wrench && sample.body_a == body3) {
        column.force = sample.force_on_body_a_per_lambda_abs;
        column.torque = sample.torque_on_body_a_per_lambda_abs;
    } else if (sample.use_custom_wrench && sample.body_b == body3) {
        column.force = sample.force_on_body_b_per_lambda_abs;
        column.torque = sample.torque_on_body_b_per_lambda_abs;
    } else if (sample.body_a == body3) {
        column.force = sample.normal_abs * column.scale;
        column.torque = (sample.point_abs - body3->GetPos()).Cross(column.force);
    } else if (sample.body_b == body3) {
        column.force = sample.normal_abs * (-column.scale);
        column.torque = (sample.point_abs - body3->GetPos()).Cross(column.force);
    }
    return column;
}

std::vector<ChSdfNcpDescriptorContactState> MakeRevClearanceStatesFromPressures(
    const std::vector<ChSdfNcpDescriptorContact>& contacts,
    const std::vector<double>& pressures) {
    std::vector<ChSdfNcpDescriptorContactState> states;
    states.reserve(contacts.size());
    for (size_t i = 0; i < contacts.size(); i++) {
        const double pressure = i < pressures.size() && std::isfinite(pressures[i]) ? std::max(0.0, pressures[i]) : 0.0;
        const double scale = RevClearanceContactForceScale(contacts[i]);
        ChSdfNcpDescriptorContactState state;
        state.sample = contacts[i];
        state.lambda_n = pressure;
        state.lambda_force = pressure * scale;
        state.penetration = std::max(0.0, -contacts[i].gap);
        state.ncp_residual = SmoothFischerBurmeister(contacts[i].gap, pressure, 1.0e-7);
        state.complementarity_error = ComplementarityError(contacts[i].gap, pressure);
        state.active = pressure > 0.0 || contacts[i].gap <= 0.0;
        states.push_back(state);
    }
    return states;
}

std::vector<double> SolveProjectedWrenchNnls(const std::vector<RevClearancePressureColumn>& columns,
                                             const RevClearanceIdealWrenchRow& target,
                                             double force_weight,
                                             double torque_weight,
                                             double regularization,
                                             int iterations) {
    const int n = static_cast<int>(columns.size());
    std::vector<double> pressures(n, 0.0);
    if (n == 0) {
        return pressures;
    }

    const std::array<double, 6> y = {force_weight * target.force.x(),
                                    force_weight * target.force.y(),
                                    force_weight * target.force.z(),
                                    torque_weight * target.torque_about_body_ref.x(),
                                    torque_weight * target.torque_about_body_ref.y(),
                                    torque_weight * target.torque_about_body_ref.z()};
    std::vector<std::array<double, 6>> a(static_cast<size_t>(n));
    std::vector<double> denom(static_cast<size_t>(n), regularization);
    for (int i = 0; i < n; i++) {
        a[static_cast<size_t>(i)] = {force_weight * columns[static_cast<size_t>(i)].force.x(),
                                    force_weight * columns[static_cast<size_t>(i)].force.y(),
                                    force_weight * columns[static_cast<size_t>(i)].force.z(),
                                    torque_weight * columns[static_cast<size_t>(i)].torque.x(),
                                    torque_weight * columns[static_cast<size_t>(i)].torque.y(),
                                    torque_weight * columns[static_cast<size_t>(i)].torque.z()};
        for (double value : a[static_cast<size_t>(i)]) {
            denom[static_cast<size_t>(i)] += value * value;
        }
        denom[static_cast<size_t>(i)] = std::max(denom[static_cast<size_t>(i)], 1.0e-30);
    }

    std::array<double, 6> residual = {-y[0], -y[1], -y[2], -y[3], -y[4], -y[5]};
    for (int iter = 0; iter < iterations; iter++) {
        double max_delta = 0.0;
        for (int i = 0; i < n; i++) {
            const auto& col = a[static_cast<size_t>(i)];
            double gradient = regularization * pressures[static_cast<size_t>(i)];
            for (int row = 0; row < 6; row++) {
                gradient += col[static_cast<size_t>(row)] * residual[static_cast<size_t>(row)];
            }
            const double next =
                std::max(0.0, pressures[static_cast<size_t>(i)] - gradient / denom[static_cast<size_t>(i)]);
            const double delta = next - pressures[static_cast<size_t>(i)];
            if (delta != 0.0) {
                pressures[static_cast<size_t>(i)] = next;
                for (int row = 0; row < 6; row++) {
                    residual[static_cast<size_t>(row)] += col[static_cast<size_t>(row)] * delta;
                }
                max_delta = std::max(max_delta, std::abs(delta));
            }
        }
        if (max_delta < 1.0e-12) {
            break;
        }
    }
    return pressures;
}

std::vector<double> RevClearancePenetrationPressures(const std::vector<ChSdfNcpDescriptorContact>& contacts,
                                                     double stiffness,
                                                     double max_pressure) {
    double active_area = 0.0;
    for (const auto& contact : contacts) {
        if (contact.gap < 0.0) {
            active_area += RevClearanceContactForceScale(contact);
        }
    }
    std::vector<double> pressures(contacts.size(), 0.0);
    if (!(active_area > 1.0e-30)) {
        return pressures;
    }
    for (size_t i = 0; i < contacts.size(); i++) {
        const double penetration = std::max(0.0, -contacts[i].gap);
        pressures[i] = std::min(max_pressure, std::max(0.0, stiffness * penetration / active_area));
    }
    return pressures;
}

std::vector<double> RevClearanceVelocityPressures(const std::vector<ChSdfNcpDescriptorContact>& contacts,
                                                  double dt,
                                                  double compliance,
                                                  double baumgarte,
                                                  double max_pressure) {
    std::vector<double> pressures(contacts.size(), 0.0);
    const double denom = std::max(compliance, 1.0e-14);
    for (size_t i = 0; i < contacts.size(); i++) {
        const auto& contact = contacts[i];
        if (!contact.body_a || !contact.body_b) {
            continue;
        }
        const ChVector3d va = RevClearancePointVelocityAbs(*contact.body_a, contact.point_abs);
        const ChVector3d vb = RevClearancePointVelocityAbs(*contact.body_b, contact.point_abs);
        const double gap_rate = contact.normal_abs.Dot(va - vb);
        const double rhs = -(gap_rate + baumgarte * std::min(0.0, contact.gap) / std::max(dt, 1.0e-12));
        pressures[i] = std::min(max_pressure, std::max(0.0, rhs / denom));
    }
    return pressures;
}

struct RevClearanceFrameSweepResult {
    std::string method;
    int contact_count = 0;
    int patch_count = 0;
    int active_pressure_count = 0;
    double min_gap = 0.0;
    double max_penetration = 0.0;
    double pressure_sum = 0.0;
    double pressure_max = 0.0;
    RevClearanceContactWrench wrench;
    double force_error_norm = 0.0;
    double torque_error_norm = 0.0;
    double weighted_wrench_error = 0.0;
    std::string notes;
};

RevClearanceFrameSweepResult EvaluateRevClearanceFramePressures(
    const std::string& method,
    const std::vector<ChSdfNcpDescriptorContact>& contacts,
    const std::vector<double>& pressures,
    const std::shared_ptr<ChBodyAuxRef>& body3,
    const RevClearanceIdealWrenchRow& ideal,
    double force_weight,
    double torque_weight,
    const std::string& notes) {
    RevClearanceFrameSweepResult out;
    out.method = method;
    out.contact_count = static_cast<int>(contacts.size());
    out.notes = notes;
    out.min_gap = contacts.empty() ? 0.0 : std::numeric_limits<double>::max();
    std::set<int> patch_ids;
    for (size_t i = 0; i < contacts.size(); i++) {
        out.min_gap = std::min(out.min_gap, contacts[i].gap);
        out.max_penetration = std::max(out.max_penetration, std::max(0.0, -contacts[i].gap));
        if (contacts[i].patch_id >= 0) {
            patch_ids.insert(contacts[i].patch_id);
        }
        const double pressure = i < pressures.size() ? std::max(0.0, pressures[i]) : 0.0;
        out.pressure_sum += pressure;
        out.pressure_max = std::max(out.pressure_max, pressure);
        if (pressure > 0.0) {
            out.active_pressure_count++;
        }
    }
    out.patch_count = static_cast<int>(patch_ids.size());
    const auto states = MakeRevClearanceStatesFromPressures(contacts, pressures);
    out.wrench = ComputeRevClearanceContactWrench(states, body3, false);
    const ChVector3d force_error = out.wrench.force_world - ideal.force;
    const ChVector3d torque_error = out.wrench.torque_world - ideal.torque_about_body_ref;
    out.force_error_norm = force_error.Length();
    out.torque_error_norm = torque_error.Length();
    out.weighted_wrench_error =
        std::sqrt(std::pow(force_weight * force_error.x(), 2.0) + std::pow(force_weight * force_error.y(), 2.0) +
                  std::pow(force_weight * force_error.z(), 2.0) + std::pow(torque_weight * torque_error.x(), 2.0) +
                  std::pow(torque_weight * torque_error.y(), 2.0) + std::pow(torque_weight * torque_error.z(), 2.0));
    return out;
}

int RunRevJointClearanceFrameSolverSweepCase(
    const std::string& output_case_name = "rev_joint_clearance_frame_solver_sweep",
    double voxel_size_override = -1.0) {
    const auto start = std::chrono::steady_clock::now();
    const auto root = GetProjectRoot();
    const auto asset_dir = root / "assets" / "rev_joint_clearance";
    const auto ideal_dir = root / "results" / "sdf_ncp_benchmarks" / "rev_joint_clearance_ideal_revolute";
    const auto ideal_trajectory_path = ideal_dir / "trajectory.csv";
    const auto ideal_wrench_path = ideal_dir / "joint_reaction_wrench.csv";
    if (!std::filesystem::exists(ideal_trajectory_path) || !std::filesystem::exists(ideal_wrench_path)) {
        const int status = RunRevJointClearanceIdealRevoluteCase();
        if (status != 0) {
            std::cerr << "Cannot generate ideal revolute reference for frame solver sweep." << std::endl;
            return status;
        }
    }

    const std::string source_case_name = "rev_joint_clearance_local_patch_delassus";
    const auto source_dir = root / "results" / "sdf_ncp_benchmarks" / source_case_name;
    const auto source_trajectory_path = source_dir / "trajectory.csv";
    if (!std::filesystem::exists(source_trajectory_path)) {
        RevClearanceRunOptions options;
        options.case_name = source_case_name;
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        const int status = RunRevJointClearanceCase(options);
        if (status != 0) {
            std::cerr << "Cannot generate source contact trajectory for frame solver sweep: " << source_case_name
                      << std::endl;
            return status;
        }
    }

    const auto source_frames = LoadRevClearanceIdealFrames(source_trajectory_path);
    const auto ideal_wrenches = LoadRevClearanceIdealWrench(ideal_wrench_path);
    if (source_frames.empty() || ideal_wrenches.empty()) {
        std::cerr << "Source contact trajectory or ideal wrench reference is empty for frame solver sweep."
                  << std::endl;
        return 1;
    }

    const auto out_dir = root / "results" / "sdf_ncp_benchmarks" / output_case_name;
    std::filesystem::create_directories(out_dir);

    RevClearanceConfig config;
    const std::string model_json = ReadTextFile(asset_dir / "rev_joint_clearance_model.json");
    config.dt = ExtractJsonDouble(model_json, "step_size", config.dt);
    config.t_end = ExtractJsonDouble(model_json, "total_time", config.t_end);
    config.max_activation_band = ExtractJsonDouble(model_json, "collision_envelope", config.max_activation_band);
    config.pressure_compliance = ExtractJsonDouble(model_json, "contact_compliance", config.pressure_compliance);
    if (voxel_size_override > 0.0) {
        config.voxel_size = voxel_size_override;
    }

    openvdb::initialize();
    const RmdModel rmd = LoadRecurDynRmdModel(asset_dir / "rev_clearance_joint.rmd");
    const RmdPart& body1_part = RequirePartByName(rmd, "Body1");
    const RmdPart& body3_part = RequirePartByName(rmd, "Body3");
    const RmdMarker& body1_cm = RequireMarkerById(rmd, body1_part.cm_marker_id);
    const RmdMarker& body3_cm = RequireMarkerById(rmd, body3_part.cm_marker_id);
    const RmdSolidContact& contact_law = RequireSolidContactByName(rmd, "GeoSurContact1");
    WriteRmdMappingAudit(out_dir / "rmd_mapping.csv", rmd);

    ChSdfTriangleMeshData bore_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "body1_subtract1_centered.obj");
    ChSdfTriangleMeshData pin_mesh =
        LoadWavefrontMeshForSdf(asset_dir / "models" / "body3_cylinder1_centered.obj");
    ChOpenVdbSdfGrid bore_sdf = BuildOpenVdbLevelSetFromMesh(
        bore_mesh, config.voxel_size, config.half_width_voxels);
    ChSdfContactSurfaceGraph pin_graph = MakeSurfaceGraphFromMeshForSdf(pin_mesh);
    ChSdfContactSampleBvh pin_bvh(pin_graph);

    auto body1 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body1_part, body1_cm, body1);
    body1->SetFixed(true);
    body1->SetPosDt(ChVector3d(0, 0, 0));
    body1->SetAngVelParent(ChVector3d(0, 0, 0));

    auto body3 = chrono_types::make_shared<ChBodyAuxRef>();
    ConfigureAuxRefBodyFromRmdPart(rmd, body3_part, body3_cm, body3);

    std::ofstream sweep(out_dir / "frame_solver_sweep.csv");
    std::ofstream best(out_dir / "best_backend_by_frame.csv");
    std::ofstream summary(out_dir / "summary.csv");
    std::ofstream comparison(out_dir / "comparison.csv");
    sweep << std::setprecision(17);
    best << std::setprecision(17);
    summary << std::setprecision(17);
    comparison << std::setprecision(17);
    sweep << "time,method,contact_count,patch_count,min_gap,max_penetration,ideal_force_norm,ideal_torque_norm,"
             "force_x,force_y,force_z,force_norm,torque_x,torque_y,torque_z,torque_norm,force_error_norm,"
             "torque_error_norm,weighted_wrench_error,pressure_sum,pressure_max,active_pressure_count,notes\n";
    best << "time,best_method,weighted_wrench_error,force_error_norm,torque_error_norm,contact_count,patch_count,"
            "active_pressure_count\n";
    summary << "method,num_frames,mean_force_error,max_force_error,mean_torque_error,max_torque_error,"
               "mean_weighted_wrench_error,max_weighted_wrench_error,mean_contact_count,mean_active_pressure_count,"
               "best_frame_count,runtime_seconds,notes\n";
    comparison << "case_name,metric,field_contact_value,sdf_ncp_value,absolute_difference,relative_difference,notes\n";

    struct MethodStats {
        int frames = 0;
        int best_count = 0;
        double force_error_sum = 0.0;
        double force_error_max = 0.0;
        double torque_error_sum = 0.0;
        double torque_error_max = 0.0;
        double weighted_error_sum = 0.0;
        double weighted_error_max = 0.0;
        double contact_count_sum = 0.0;
        double active_pressure_count_sum = 0.0;
        std::string notes;
    };
    std::map<std::string, MethodStats> stats;

    bool finite_all = true;
    for (const auto& frame : source_frames) {
        body3->SetFrameRefToAbs(ChFrame<>(frame.body3_pos, frame.body3_rot));
        body3->SetPosDt(frame.body3_vel);
        body3->SetAngVelParent(frame.body3_omega);

        const auto contacts = BuildRevClearanceDescriptorContacts(body1, body3, bore_sdf, pin_graph, pin_bvh, config);
        std::vector<RevClearancePressureColumn> columns;
        columns.reserve(contacts.size());
        for (const auto& contact : contacts) {
            columns.push_back(RevClearanceBody3PressureColumn(contact, body3));
        }

        const auto ideal = InterpolateRevClearanceIdealWrench(ideal_wrenches, frame.time);
        const double force_scale = std::max(1.0, ideal.force.Length());
        const double torque_scale = std::max(1.0, ideal.torque_about_body_ref.Length());
        const double force_weight = 1.0 / force_scale;
        const double torque_weight = 1.0 / torque_scale;
        double lever_sum = 0.0;
        int lever_count = 0;
        for (const auto& column : columns) {
            const double force_norm = column.force.Length();
            if (force_norm > 1.0e-30) {
                lever_sum += column.torque.Length() / force_norm;
                lever_count++;
            }
        }
        const double lever_scale = lever_count > 0 ? std::max(1.0e-4, lever_sum / static_cast<double>(lever_count)) :
                                                     1.0e-2;
        const double lever_torque_weight = 1.0 / std::max(1.0, force_scale * lever_scale);
        double total_contact_area = 0.0;
        for (const auto& contact : contacts) {
            total_contact_area += RevClearanceContactForceScale(contact);
        }
        const double max_pressure =
            std::max(1.0e8, 10.0 * std::abs(contact_law.stiffness) *
                                  std::max(config.max_activation_band, contact_law.bpen) /
                                  std::max(total_contact_area, 1.0e-12));

        std::vector<RevClearanceFrameSweepResult> results;
        results.push_back(EvaluateRevClearanceFramePressures(
            "zero_pressure",
            contacts,
            std::vector<double>(contacts.size(), 0.0),
            body3,
            ideal,
            force_weight,
            torque_weight,
            "diagnostic lower bound: same geometry with no contact pressure"));
        results.push_back(EvaluateRevClearanceFramePressures(
            "penetration_pressure",
            contacts,
            RevClearancePenetrationPressures(contacts, contact_law.stiffness, max_pressure),
            body3,
            ideal,
            force_weight,
            torque_weight,
            "frame-local RecurDyn stiffness-style penetration pressure; no dynamics solve"));
        results.push_back(EvaluateRevClearanceFramePressures(
            "velocity_baumgarte_pressure",
            contacts,
            RevClearanceVelocityPressures(contacts, config.dt, config.pressure_compliance, 0.2, max_pressure),
            body3,
            ideal,
            force_weight,
            torque_weight,
            "frame-local velocity/baumgarte pressure estimate; no global descriptor coupling"));
        results.push_back(EvaluateRevClearanceFramePressures(
            "oracle_force_nnls",
            contacts,
            SolveProjectedWrenchNnls(columns, ideal, force_weight, 0.0, 1.0e-18, 80),
            body3,
            ideal,
            force_weight,
            torque_weight,
            "diagnostic upper bound: nonnegative sample pressures fitted to ideal force only"));
        results.push_back(EvaluateRevClearanceFramePressures(
            "oracle_wrench_nnls",
            contacts,
            SolveProjectedWrenchNnls(columns, ideal, force_weight, torque_weight, 1.0e-18, 160),
            body3,
            ideal,
            force_weight,
            torque_weight,
            "diagnostic upper bound: nonnegative sample pressures fitted to ideal force and torque"));
        results.push_back(EvaluateRevClearanceFramePressures(
            "oracle_wrench_lever_scaled_nnls",
            contacts,
            SolveProjectedWrenchNnls(columns, ideal, force_weight, lever_torque_weight, 1.0e-18, 160),
            body3,
            ideal,
            force_weight,
            lever_torque_weight,
            "diagnostic upper bound: torque scaled by ideal force times mean contact lever arm"));

        const auto best_it = std::min_element(results.begin(), results.end(), [](const auto& a, const auto& b) {
            return a.weighted_wrench_error < b.weighted_wrench_error;
        });
        if (best_it != results.end()) {
            best << frame.time << "," << best_it->method << "," << best_it->weighted_wrench_error << ","
                 << best_it->force_error_norm << "," << best_it->torque_error_norm << "," << best_it->contact_count
                 << "," << best_it->patch_count << "," << best_it->active_pressure_count << "\n";
        }

        for (const auto& result : results) {
            sweep << frame.time << "," << result.method << "," << result.contact_count << "," << result.patch_count
                  << "," << result.min_gap << "," << result.max_penetration << "," << ideal.force.Length() << ","
                  << ideal.torque_about_body_ref.Length() << "," << result.wrench.force_world.x() << ","
                  << result.wrench.force_world.y() << "," << result.wrench.force_world.z() << ","
                  << result.wrench.force_norm << "," << result.wrench.torque_world.x() << ","
                  << result.wrench.torque_world.y() << "," << result.wrench.torque_world.z() << ","
                  << result.wrench.torque_norm << "," << result.force_error_norm << "," << result.torque_error_norm
                  << "," << result.weighted_wrench_error << "," << result.pressure_sum << ","
                  << result.pressure_max << "," << result.active_pressure_count << "," << result.notes << "\n";
            finite_all = finite_all && std::isfinite(result.weighted_wrench_error) &&
                         std::isfinite(result.force_error_norm) && std::isfinite(result.torque_error_norm);
            auto& method_stats = stats[result.method];
            method_stats.frames++;
            method_stats.force_error_sum += result.force_error_norm;
            method_stats.force_error_max = std::max(method_stats.force_error_max, result.force_error_norm);
            method_stats.torque_error_sum += result.torque_error_norm;
            method_stats.torque_error_max = std::max(method_stats.torque_error_max, result.torque_error_norm);
            method_stats.weighted_error_sum += result.weighted_wrench_error;
            method_stats.weighted_error_max = std::max(method_stats.weighted_error_max, result.weighted_wrench_error);
            method_stats.contact_count_sum += static_cast<double>(result.contact_count);
            method_stats.active_pressure_count_sum += static_cast<double>(result.active_pressure_count);
            method_stats.notes = result.notes;
        }
        if (best_it != results.end()) {
            stats[best_it->method].best_count++;
        }
    }

    const double runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    for (const auto& [method, method_stats] : stats) {
        const double denom = static_cast<double>(std::max(1, method_stats.frames));
        summary << method << "," << method_stats.frames << "," << method_stats.force_error_sum / denom << ","
                << method_stats.force_error_max << "," << method_stats.torque_error_sum / denom << ","
                << method_stats.torque_error_max << "," << method_stats.weighted_error_sum / denom << ","
                << method_stats.weighted_error_max << "," << method_stats.contact_count_sum / denom << ","
                << method_stats.active_pressure_count_sum / denom << "," << method_stats.best_count << ","
                << runtime << "," << method_stats.notes << "\n";
        comparison << "rev_joint_clearance_frame_solver_sweep," << method << "_mean_weighted_wrench_error,0,"
                   << method_stats.weighted_error_sum / denom << "," << method_stats.weighted_error_sum / denom
                   << ",0," << method_stats.notes << "\n";
    }
    comparison << output_case_name << ",num_frames,0," << source_frames.size() << ","
               << source_frames.size()
               << ",0,frame-frozen diagnostic using source contact trajectory " << source_case_name
               << " and ideal revolute wrench\n";
    comparison << output_case_name << ",source_case,0,0,0,0," << source_case_name << "\n";
    comparison << output_case_name << ",voxel_size,0," << config.voxel_size << "," << config.voxel_size
               << ",0,OpenVDB voxel size used only for frozen-frame geometry query\n";

    std::cout << "Wrote " << (out_dir / "frame_solver_sweep.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "best_backend_by_frame.csv").string() << std::endl;
    std::cout << "Wrote " << (out_dir / "summary.csv").string() << std::endl;
    return finite_all ? 0 : 1;
}

int RunMode(const std::string& mode) {
    if (mode == "cam") {
        return RunCamChronoMbdCase();
    }
    if (mode == "cam_recurdyn_compare") {
        return RunCamChronoMbdCase("cam_recurdyn_compare", 3.0, 0.001);
    }
    if (mode == "cam_recurdyn_compare_reverse_motor") {
        return RunCamChronoMbdCase("cam_recurdyn_compare_reverse_motor", 3.0, 0.001, false, 1.0);
    }
    if (mode == "cam_recurdyn_solid_contact") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact", 3.0, 0.001, true);
    }
    if (mode == "cam_recurdyn_solid_contact_voxel_004") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact_voxel_004", 3.0, 0.001, true, 1.0, 0.004);
    }
    if (mode == "cam_recurdyn_solid_contact_voxel_002") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact_voxel_002", 3.0, 0.001, true, 1.0, 0.002);
    }
    if (mode == "cam_recurdyn_solid_contact_voxel_001") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact_voxel_001", 3.0, 0.001, true, 1.0, 0.001);
    }
    if (mode == "cam_recurdyn_solid_contact_voxel_0005") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact_voxel_0005", 3.0, 0.001, true, 1.0, 0.0005);
    }
    if (mode == "cam_recurdyn_solid_contact_reverse_motor") {
        return RunCamChronoMbdCase("cam_recurdyn_solid_contact_reverse_motor", 3.0, 0.001, true, 1.0);
    }
    if (mode == "cam_reduced") {
        return RunCamFullGeometryCase();
    }
    if (mode == "eccentric_roller") {
        return RunCamLikeFullGeometryCase("eccentric_roller");
    }
    if (mode == "onset_stress") {
        return RunCamLikeFullGeometryCase("onset_stress");
    }
    if (mode == "simple_gear") {
        return RunSimpleGearFullGeometryCase();
    }
    if (mode == "simple_gear_recurdyn_compare" || mode == "simple_gear_3s") {
        return RunSimpleGearFullGeometryCase("simple_gear_recurdyn_compare", 1.0, 0.001);
    }
    if (mode == "simple_gear_dt_001") {
        return RunSimpleGearFullGeometryCase("simple_gear_dt_001", 1.0, 0.001);
    }
    if (mode == "simple_gear_dt_0005") {
        return RunSimpleGearFullGeometryCase("simple_gear_dt_0005", 1.0, 0.0005);
    }
    if (mode == "simple_gear_dt_0001") {
        return RunSimpleGearFullGeometryCase("simple_gear_dt_0001", 1.0, 0.0001);
    }
    if (mode == "simple_gear_ggeomcontact_calibration" || mode == "simple_gear_recurdyn_ggeomcontact") {
        return RunSimpleGearFullGeometryCase("simple_gear_ggeomcontact_calibration",
                                             1.0,
                                             0.001,
                                             SimpleGearBackendMode::RecurDynGGeomContactLaw);
    }
    if (mode == "rev_joint_clearance") {
        return RunRevJointClearanceCase();
    }
    if (mode == "rev_joint_clearance_patch_pressure") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_patch_pressure";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.patch_pressure_aggregation = true;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_projection") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_projection";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_projection = true;
        options.descriptor_patch_projection_strength = 0.15;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_projection_strong") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_projection_strong";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_projection = true;
        options.descriptor_patch_projection_strength = 0.35;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_lcp") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_lcp";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_projection = true;
        options.descriptor_patch_projection_strength = 1.0;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_delassus") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_delassus";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_patch_laplacian_compliance_scale = 0.0;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_delassus_laplacian") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_delassus_laplacian";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_patch_laplacian_compliance_scale = 0.25;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_delassus") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_delassus";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.2;
        options.descriptor_patch_laplacian_compliance_scale = 0.0;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_delassus_laplacian") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_delassus_laplacian";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.2;
        options.descriptor_patch_laplacian_compliance_scale = 0.25;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_delassus_soft") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_delassus_soft";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.1;
        options.descriptor_patch_laplacian_compliance_scale = 0.0;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_wrench_closure") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_wrench_closure";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.2;
        options.descriptor_patch_laplacian_compliance_scale = 0.0;
        options.descriptor_patch_wrench_closure_strength = 0.05;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_wrench_closure_laplacian") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_wrench_closure_laplacian";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.2;
        options.descriptor_patch_laplacian_compliance_scale = 0.25;
        options.descriptor_patch_wrench_closure_strength = 0.05;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_descriptor_patch_velocity_wrench_closure_strong") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_descriptor_patch_velocity_wrench_closure_strong";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.descriptor_patch_delassus_solve = true;
        options.descriptor_velocity_level_ncp = true;
        options.descriptor_velocity_baumgarte = 0.2;
        options.descriptor_patch_laplacian_compliance_scale = 0.25;
        options.descriptor_patch_wrench_closure_strength = 0.2;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_manifold_no_prediction") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_manifold_no_prediction";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.manifold_use_prediction = false;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_bidirectional") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_bidirectional";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.bidirectional_contact = true;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_bidirectional_no_prediction") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_bidirectional_no_prediction";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.bidirectional_contact = true;
        options.manifold_use_prediction = false;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_bidirectional_patch_lcp") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_bidirectional_patch_lcp";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.bidirectional_contact = true;
        options.descriptor_patch_projection = true;
        options.descriptor_patch_projection_strength = 1.0;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_pressure") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_pressure";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_pressure_compliance_scale = 0.1;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_pressure_balanced") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_pressure_balanced";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_pressure_stiff") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_pressure_stiff";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_pressure_compliance_scale = 0.01;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_delassus") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_delassus";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 25.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque100") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque100";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 100.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque250") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque250";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 250.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque500") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque500";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 500.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque1000") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque1000";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1000.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_global_wrench_demand") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_global_wrench_demand";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_global_wrench_demand = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_global_wrench_demand_torque1000") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_global_wrench_demand_torque1000";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_global_wrench_demand = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1000.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_global_wrench_demand_short") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_global_wrench_demand_short";
        options.t_end_override = 0.5;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_global_wrench_demand = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_global_wrench_demand_torque1000_short") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_global_wrench_demand_torque1000_short";
        options.t_end_override = 0.5;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_global_wrench_demand = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1000.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_augmented_wrench_closure_torque1000") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_augmented_wrench_closure_torque1000";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_integrated_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1000.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_augmented_wrench_closure_torque500") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_augmented_wrench_closure_torque500";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_integrated_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 500.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_wrench_closure_torque2000") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_wrench_closure_torque2000";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 2000.0;
        options.local_patch_wrench_regularization = 1.0e-18;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_bidirectional_local_patch_wrench_closure") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_bidirectional_local_patch_wrench_closure";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.bidirectional_contact = true;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_wrench_closure = true;
        options.local_patch_pressure_compliance_scale = 0.03;
        options.local_patch_wrench_force_weight = 1.0;
        options.local_patch_wrench_torque_weight = 1.0;
        options.local_patch_wrench_regularization = 1.0e-8;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_delassus_soft") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_delassus_soft";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_pressure_compliance_scale = 0.1;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_local_patch_delassus_stiff") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_local_patch_delassus_stiff";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.local_patch_pressure_solve = true;
        options.local_patch_delassus = true;
        options.local_patch_pressure_compliance_scale = 0.01;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_ideal_revolute") {
        return RunRevJointClearanceIdealRevoluteCase();
    }
    if (mode == "rev_joint_clearance_frame_solver_sweep") {
        return RunRevJointClearanceFrameSolverSweepCase();
    }
    if (mode == "rev_joint_clearance_frame_solver_sweep_voxel_001") {
        return RunRevJointClearanceFrameSolverSweepCase("rev_joint_clearance_frame_solver_sweep_voxel_001", 0.001);
    }
    if (mode == "rev_joint_clearance_frame_solver_sweep_voxel_0005") {
        return RunRevJointClearanceFrameSolverSweepCase("rev_joint_clearance_frame_solver_sweep_voxel_0005", 0.0005);
    }
    if (mode == "rev_joint_clearance_ggeomcontact_calibration" ||
        mode == "rev_joint_clearance_recurdyn_ggeomcontact") {
        return RunRevJointClearanceCase("rev_joint_clearance_ggeomcontact_calibration", 3.0, 0.001, true);
    }
    if (mode == "rev_joint_clearance_ggeomcontact_hht") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_ggeomcontact_hht";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.use_recurdyn_ggeomcontact_law = true;
        options.timestepper_type = ChTimestepper::Type::HHT;
        options.timestepper_label = "hht_alpha_minus_0p2_fixed_step";
        options.use_step_control = false;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_ggeomcontact_euler_substep") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_ggeomcontact_euler_substep";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.use_recurdyn_ggeomcontact_law = true;
        options.contact_substepping = true;
        options.startup_substepping = true;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_ggeomcontact_hht_substep") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_ggeomcontact_hht_substep";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.use_recurdyn_ggeomcontact_law = true;
        options.timestepper_type = ChTimestepper::Type::HHT;
        options.timestepper_label = "hht_alpha_minus_0p2_fixed_step";
        options.use_step_control = false;
        options.contact_substepping = true;
        options.startup_substepping = true;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_hht_substep") {
        RevClearanceRunOptions options;
        options.case_name = "rev_joint_clearance_hht_substep";
        options.t_end_override = 3.0;
        options.dt_override = 0.001;
        options.timestepper_type = ChTimestepper::Type::HHT;
        options.timestepper_label = "hht_alpha_minus_0p2_fixed_step";
        options.use_step_control = false;
        options.contact_substepping = true;
        options.startup_substepping = true;
        return RunRevJointClearanceCase(options);
    }
    if (mode == "rev_joint_clearance_integrator_sweep") {
        int status = 0;
        status |= RunRevJointClearanceCase("rev_joint_clearance_ggeomcontact_calibration", 3.0, 0.001, true);
        status |= RunMode("rev_joint_clearance_ggeomcontact_hht");
        status |= RunMode("rev_joint_clearance_ggeomcontact_euler_substep");
        status |= RunMode("rev_joint_clearance_ggeomcontact_hht_substep");
        return status;
    }
    if (mode == "rev_joint_clearance_short") {
        return RunRevJointClearanceCase("rev_joint_clearance_short", 0.2, 0.001);
    }
    if (mode == "all") {
        int status = 0;
        status |= RunCamChronoMbdCase();
        status |= RunCamLikeFullGeometryCase("eccentric_roller");
        status |= RunCamLikeFullGeometryCase("onset_stress");
        status |= RunSimpleGearFullGeometryCase();
        status |= RunRevJointClearanceCase();
        return status;
    }
    std::cerr << "Unknown OpenVDB SDF-NCP benchmark mode: " << mode << std::endl;
    return 2;
}

}  // namespace

int main(int argc, char* argv[]) {
    const std::string mode = argc > 1 ? argv[1] : "all";
    try {
        return RunMode(mode);
    } catch (const std::exception& e) {
        std::cerr << "SDF-NCP OpenVDB benchmark failed: " << e.what() << std::endl;
        return 1;
    }
}
