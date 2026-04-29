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
#include "chrono/physics/ChSdfNcpConstraintContact.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolver.h"

#include <openvdb/openvdb.h>

using chrono::ChVector3d;
using chrono::ChBody;
using chrono::ChBodyAuxRef;
using chrono::ChFrame;
using chrono::ChFunctionConst;
using chrono::ChLinkMateGeneric;
using chrono::ChLinkMotorRotationSpeed;
using chrono::ChMatrix33d;
using chrono::ChSolver;
using chrono::ChSystemNSC;
using chrono::ChTimestepper;
using chrono::ChVectorDynamic;
using chrono::UpdateFlags;
using chrono::sdfcontact::BuildOpenVdbLevelSetFromMesh;
using chrono::sdfcontact::ChOpenVdbSdfGrid;
using chrono::sdfcontact::ChSdfContactSampleQuery;
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
        if (trimmed.rfind("CSURFACE /", 0) == 0) {
            const int id = ParseRmdEntityId(trimmed);
            model.surfaces[id].id = id;
            start_block(Block::Surface, id);
            continue;
        }
        if (trimmed.rfind("SOLIDCONTACT /", 0) == 0) {
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
            } else if (Contains(trimmed, "ICSURFACEID =")) {
                const auto values = ParseNumbersAfterEquals(trimmed);
                if (!values.empty()) {
                    contact.i_surface_id = static_cast<int>(std::round(values.front()));
                }
            } else if (Contains(trimmed, "JCSURFACEID =")) {
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
            << "; RDF " << contact.rdf << "\n";
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
            "generic active-patch area-normalized SDF-NCP assembly";
    } else if (case_name == "onset_stress") {
        model_json = ReadTextFile(config.asset_dir / "onset_stress_model.json");
        config.cam_obj = config.asset_dir / "models" / "onset_cam.obj";
        config.follower_obj = config.asset_dir / "models" / "roller_follower.obj";
        config.notes =
            "full onset_stress OBJ/OpenVDB SDF; adaptive active-band connected patch SDF-NCP samples; "
            "generic active-patch area-normalized SDF-NCP assembly";
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

std::vector<int> BuildCamLikePatchSamples(const ChOpenVdbSdfGrid& cam_sdf,
                                          const ChSdfContactSurfaceGraph& follower_graph,
                                          double follower_y,
                                          double follower_vy,
                                          double theta_cam,
                                          double omega_cam,
                                          double activation_band,
                                          double max_activation_band,
                                          int min_patch_samples,
                                          int* patch_count = nullptr) {
    std::vector<std::pair<int, double>> sampled;
    sampled.reserve(follower_graph.samples.size());
    double min_phi = std::numeric_limits<double>::max();
    for (const auto& sample : follower_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        const auto query = QueryFollowerSampleAgainstCamLike(
            cam_sdf, follower_graph, sample.id, follower_y, follower_vy, theta_cam, omega_cam);
        sampled.emplace_back(sample.id, query.phi);
        min_phi = std::min(min_phi, query.phi);
    }
    if (min_phi > std::max(activation_band, max_activation_band)) {
        if (patch_count) {
            *patch_count = 0;
        }
        return {};
    }

    std::vector<int> active;
    active.reserve(follower_graph.samples.size());
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
    double min_gap = std::numeric_limits<double>::max();
    for (const auto& sample : follower_graph.samples) {
        if (sample.area <= 0.0) {
            continue;
        }
        const auto query = QueryFollowerSampleAgainstCamLike(
            cam_sdf, follower_graph, sample.id, follower_y, follower_vy, theta_cam, omega_cam);
        min_gap = std::min(min_gap, query.phi);
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
                                   const CamLikeState& state) {
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
                                                                  &patch_count);
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
            const CamLikeStepResult step = SolveCamLikeStep(cam_sdf, follower_graph, config, state);
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
                const double ncp = c < step.diagnostics.ncp_residuals.size() ? step.diagnostics.ncp_residuals[c] : 0.0;
                const double comp = c < step.diagnostics.complementarity_errors.size() ?
                                        step.diagnostics.complementarity_errors[c] :
                                        0.0;
                trajectory << time << ",follower,0," << state.y << ",0,1,0,0,0,0," << state.vy
                            << ",0,0,0,0," << candidate.sample_id << "," << candidate.gap << "," << lambda
                            << "," << candidate.weight << "," << lambda * candidate.weight << ","
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
                                                                   &patch_count);
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
};

int GearSamplePersistentId(int direction, int sample_id) {
    return direction == 0 ? sample_id : 1000000 + sample_id;
}

enum class SimpleGearBackendMode {
    SdfNcp,
    RecurDynGGeomContactLaw,
};

struct SimpleGearGGeomContactLawSI {
    double stiffness = 0.0;
    double damping = 0.0;
    double exponent = 1.0;
    double buffer_penetration = 0.0;
};

SimpleGearGGeomContactLawSI MakeSimpleGearGGeomContactLawSI(const SimpleGearGGeomContactInfo& contact) {
    constexpr double length_scale_m = 1.0e-3;  // simple_gear RMD declares LENGTH = Millimeter.
    SimpleGearGGeomContactLawSI law;
    law.exponent = contact.korder > 0.0 ? contact.korder : 1.0;
    law.buffer_penetration = contact.bpen * length_scale_m;
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
                                                                     GearPatchMemory* memory = nullptr) {
    struct SamplePhi {
        int id = -1;
        double effective_gap = 0.0;
        double area = 0.0;
    };

    auto build_direction = [&](int direction, const ChSdfContactSurfaceGraph& graph) {
        std::vector<SamplePhi> sampled;
        sampled.reserve(graph.samples.size());
        double min_gap = std::numeric_limits<double>::max();
        for (const auto& sample : graph.samples) {
            if (sample.area <= 0.0) {
                continue;
            }
            const auto candidate = direction == 0 ?
                                       QueryGear22SampleAgainstGear21(gear21_sdf,
                                                                      gear22_graph,
                                                                      pose21,
                                                                      pose22,
                                                                      sample.id,
                                                                      theta21,
                                                                      theta22) :
                                       QueryGear21SampleAgainstGear22(gear22_sdf,
                                                                      gear21_graph,
                                                                      pose21,
                                                                      pose22,
                                                                      sample.id,
                                                                      theta21,
                                                                      theta22);
            const double effective_gap = candidate.gap - config.sdf_ncp_contact_offset;
            sampled.push_back(SamplePhi{sample.id, effective_gap, sample.area});
            min_gap = std::min(min_gap, effective_gap);
        }

        std::vector<int> active;
        if (sampled.empty()) {
            return std::vector<GearContactSampleRef>{};
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
                if (static_cast<int>(active.size()) >= std::max(1, config.min_patch_samples) || band >= band_limit) {
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
            const int patch_id = direction == 0 ? min_sample_id : 1000000 + min_sample_id;
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

    std::vector<GearContactSampleRef> refs22 = build_direction(0, gear22_graph);
    std::vector<GearContactSampleRef> refs21 = build_direction(1, gear21_graph);
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
        const auto candidate = QueryGear22SampleAgainstGear21(
            gear21_sdf, gear22_graph, pose21, pose22, sample.id, theta21, theta22);
        min_gap = std::min(min_gap, candidate.gap);
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
        const auto candidate = QueryGear21SampleAgainstGear22(
            gear22_sdf, gear21_graph, pose21, pose22, sample.id, theta21, theta22);
        min_gap = std::min(min_gap, candidate.gap);
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
                                   GearPatchMemory* memory = nullptr) {
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
                                                                                       &active_memory);
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
                                                                                         &candidate_memory);
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
        const double law_penetration = std::max(0.0, law.buffer_penetration - candidate.gap);
        const double spring_force = law.stiffness > 0.0 && law_penetration > 0.0 ?
                                        law.stiffness * std::pow(law_penetration, law.exponent) :
                                        0.0;
        const double damping_force = law_penetration > 0.0 ? law.damping * std::max(0.0, -gap_rate) : 0.0;
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
                                                          GearPatchMemory* memory = nullptr) {
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
                                                                                             memory);

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
        config.sdf_ncp_contact_offset = config.sdf_ncp_contact_offset_scale * ggeom_law.buffer_penetration;
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
                      "and GGEOMCONTACT metadata parsed from simple gear.rmd; hard frictionless NCP does not apply "
                      "RMD K/C/KORDER";
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
                                                                                              &patch_memory);
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
                                                                  &patch_memory);
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
                                           &patch_memory);
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
            const double ncp = c < step.diagnostics.ncp_residuals.size() ? step.diagnostics.ncp_residuals[c] : 0.0;
            const double comp = c < step.diagnostics.complementarity_errors.size() ?
                                    step.diagnostics.complementarity_errors[c] :
                                    0.0;
            const double q0 = std::cos(0.5 * state.theta22);
            const double q1 = std::sin(0.5 * state.theta22);
            trajectory << time << ",gear22," << pose22.center.x() << "," << pose22.center.y() << ","
                       << pose22.center.z() << "," << q0 << "," << q1 << ",0,0,0,0,0,"
                       << state.omega22 << ",0,0," << candidate.sample_id << "," << candidate.gap << ","
                       << lambda << "," << candidate.weight << "," << lambda * candidate.weight << ","
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
               << ggeom_law.buffer_penetration
               << ",0,0,RecurDyn BPEN converted from millimeter RMD length units to meters\n";
    comparison << config.case_name << ",ggeomcontact_K_SI," << ggeom_law.stiffness << "," << ggeom_law.stiffness
               << ",0,0,RecurDyn K converted for meter-based gap with KORDER exponent\n";
    comparison << config.case_name << ",ggeomcontact_C_SI," << ggeom_law.damping << "," << ggeom_law.damping
               << ",0,0,RecurDyn C converted from N/(mm/s) to N/(m/s)\n";
    comparison << config.case_name << ",sdf_ncp_contact_offset_m," << ggeom_law.buffer_penetration << ","
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
    const CamLikeConfig& config) {
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
    for (const auto& sample : follower_graph.samples) {
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

    double total_area = 0.0;
    for (const auto& item : active) {
        total_area += std::max(0.0, item.area);
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
        contact.weight = total_area > 1.0e-30 ? std::max(0.0, item.area) / total_area : uniform_weight;
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
        candidate.patch_id = 0;
        candidate.gap = state.sample.gap;
        candidate.area = state.sample.weight;
        candidate.weight = state.sample.weight;
        candidate.normal = state.sample.normal_abs;
        candidate.world_pos = state.sample.point_abs;
        diag.candidates.push_back(candidate);
        diag.lambdas.push_back(state.lambda_n);
        diag.ncp_residuals.push_back(state.ncp_residual);
        diag.complementarity_errors.push_back(state.complementarity_error);
        diag.min_gap = std::min(diag.min_gap, state.sample.gap);
        diag.max_penetration = std::max(diag.max_penetration, state.penetration);
        diag.max_lambda_n = std::max(diag.max_lambda_n, state.lambda_force);
        diag.ncp_residual_norm += state.ncp_residual * state.ncp_residual;
        diag.max_complementarity_error = std::max(diag.max_complementarity_error, state.complementarity_error);
    }
    diag.ncp_residual_norm = std::sqrt(diag.ncp_residual_norm);
    return diag;
}

struct RecurDynSolidContactLaw {
    double stiffness = 0.0;
    double damping = 0.0;
    double exponent = 1.0;
    double buffer_penetration = 0.0;
};

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
        const ChVector3d point_loc = body->TransformPointParentToLocal(point_abs);
        const ChVector3d force_loc = body->TransformDirectionParentToLocal(force_abs);
        const ChVector3d torque_loc = point_loc.Cross(force_loc);
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
            const double penetration = std::max(0.0, m_law.buffer_penetration - sample.gap);
            const double spring_force =
                m_law.stiffness > 0.0 && penetration > 0.0 ? m_law.stiffness * std::pow(penetration, m_law.exponent)
                                                            : 0.0;
            const double damping_force = penetration > 0.0 ? m_law.damping * std::max(0.0, -gap_rate) : 0.0;
            const double raw_force = std::max(0.0, spring_force + damping_force);

            ChSdfNcpDescriptorContactState state;
            state.sample = sample;
            state.lambda_n = raw_force;
            state.lambda_force = raw_force * sample.weight;
            state.penetration = std::max(0.0, -sample.gap);
            state.ncp_residual = 0.0;
            state.complementarity_error = ComplementarityError(sample.gap, state.lambda_force);
            state.active = raw_force > 0.0 || sample.gap <= m_law.buffer_penetration;
            m_states.push_back(state);
        }
    }

    ContactGenerator m_generator;
    RecurDynSolidContactLaw m_law;
    std::vector<ChSdfNcpDescriptorContactState> m_states;
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
                           "full cam OBJ/OpenVDB geometry" :
                           "Chrono MBD path: ChSystemNSC + cam rotation motor + follower translational joint + "
                           "descriptor-injected frictionless SDF-NCP unilateral constraints; "
                           "RMD part/marker/joint/surface mapping; full cam OBJ/OpenVDB geometry";

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

    std::shared_ptr<ChSdfNcpConstraintContactSet> contact_item;
    std::shared_ptr<ChSdfPenaltyContactForceSet> penalty_item;
    auto generator = [cam, follower, &cam_sdf, &follower_graph, ncp_config]() {
        return BuildCamDescriptorContacts(cam, follower, cam_sdf, follower_graph, ncp_config);
    };
    if (use_recurdyn_solid_contact_law) {
        RecurDynSolidContactLaw law;
        law.stiffness = solid_contact.stiffness;
        law.damping = solid_contact.damping;
        law.exponent = solid_contact.korder > 0.0 ? solid_contact.korder : 1.0;
        law.buffer_penetration = solid_contact.bpen;
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
                const double ncp = c < diag.ncp_residuals.size() ? diag.ncp_residuals[c] : 0.0;
                const double comp = c < diag.complementarity_errors.size() ? diag.complementarity_errors[c] : 0.0;
                trajectory << time << ",follower," << follower->GetPos().x() << "," << follower->GetPos().y() << ","
                           << follower->GetPos().z() << "," << fq.e0() << "," << fq.e1() << "," << fq.e2() << ","
                           << fq.e3() << "," << follower->GetPosDt().x() << "," << follower->GetPosDt().y() << ","
                           << follower->GetPosDt().z() << "," << follower_w.x() << "," << follower_w.y() << ","
                           << follower_w.z() << "," << candidate.sample_id << "," << candidate.gap << ","
                           << lambda << "," << candidate.weight << "," << lambda * candidate.weight << ","
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
        "generic active-patch area-normalized SDF-NCP assembly; cam speed from simple_cam.rmd";
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
            const CamLikeStepResult step = SolveCamLikeStep(cam_sdf, follower_graph, ncp_config, state);
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
                                                                   &patch_count);
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
            const double ncp = c < diag.ncp_residuals.size() ? diag.ncp_residuals[c] : 0.0;
            const double comp = c < diag.complementarity_errors.size() ? diag.complementarity_errors[c] : 0.0;
            trajectory << time << ",follower," << follower_cm.x() << "," << follower_cm.y() << ","
                       << follower_cm.z() << ",1,0,0,0,0," << state.vy << ",0,0,0,0,"
                       << candidate.sample_id << "," << candidate.gap << "," << lambda << ","
                       << candidate.weight << "," << lambda * candidate.weight << ","
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
    if (mode == "all") {
        int status = 0;
        status |= RunCamChronoMbdCase();
        status |= RunCamLikeFullGeometryCase("eccentric_roller");
        status |= RunCamLikeFullGeometryCase("onset_stress");
        status |= RunSimpleGearFullGeometryCase();
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
