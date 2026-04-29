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
// Optional OpenVDB SDF geometry frontend.
//
// This header contains only mesh/SDF loading and query utilities. It deliberately
// does not contain field-contact pressure laws, patch extraction, NCP residuals,
// contact multipliers, or solver logic. Pressure-field and SDF-NCP demos can use
// this frontend to obtain ChSdfContactSampleQuery values from the same sparse
// SDF representation.
//
// =============================================================================

#ifndef CH_SDF_OPENVDB_H
#define CH_SDF_OPENVDB_H

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "chrono/collision/ChSdfContactGeometry.h"
#include "chrono/core/ChVector3.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

namespace chrono {
namespace sdfcontact {

struct ChSdfTriangleMeshData {
    std::vector<ChVector3d> vertices;
    std::vector<ChSdfTriangleFace> faces;
};

struct ChOpenVdbSdfGrid {
    openvdb::FloatGrid::Ptr grid;
    double voxel_size = 0.0;

    bool Empty() const {
        return !grid;
    }

    double SamplePhi(const ChVector3d& local_pos) const {
        if (!grid) {
            throw std::runtime_error("ChOpenVdbSdfGrid has no grid.");
        }

        openvdb::tools::GridSampler<openvdb::FloatGrid::TreeType, openvdb::tools::BoxSampler> sampler(
            grid->tree(), grid->transform());
        return static_cast<double>(sampler.wsSample(openvdb::Vec3d(local_pos.x(), local_pos.y(), local_pos.z())));
    }

    ChVector3d SampleRawGradient(const ChVector3d& local_pos) const {
        if (voxel_size <= 0.0) {
            throw std::runtime_error("ChOpenVdbSdfGrid voxel_size must be positive.");
        }

        const double h = 0.5 * voxel_size;
        ChVector3d grad((SamplePhi(local_pos + ChVector3d(h, 0, 0)) -
                             SamplePhi(local_pos - ChVector3d(h, 0, 0))) /
                            (2.0 * h),
                        (SamplePhi(local_pos + ChVector3d(0, h, 0)) -
                             SamplePhi(local_pos - ChVector3d(0, h, 0))) /
                            (2.0 * h),
                        (SamplePhi(local_pos + ChVector3d(0, 0, h)) -
                             SamplePhi(local_pos - ChVector3d(0, 0, h))) /
                            (2.0 * h));
        return grad;
    }

    ChVector3d SampleGradient(const ChVector3d& local_pos) const {
        const ChVector3d grad = SampleRawGradient(local_pos);
        return NormalizeOrFallback(grad, ChVector3d(0, 1, 0));
    }

    void SampleHessian(const ChVector3d& local_pos,
                       ChVector3d& hessian_x,
                       ChVector3d& hessian_y,
                       ChVector3d& hessian_z) const {
        if (voxel_size <= 0.0) {
            throw std::runtime_error("ChOpenVdbSdfGrid voxel_size must be positive.");
        }

        const double h = voxel_size;
        const double inv_h2 = 1.0 / (h * h);
        const double inv_4h2 = 0.25 * inv_h2;
        const double c = SamplePhi(local_pos);

        const ChVector3d ex(h, 0, 0);
        const ChVector3d ey(0, h, 0);
        const ChVector3d ez(0, 0, h);

        const double xp = SamplePhi(local_pos + ex);
        const double xm = SamplePhi(local_pos - ex);
        const double yp = SamplePhi(local_pos + ey);
        const double ym = SamplePhi(local_pos - ey);
        const double zp = SamplePhi(local_pos + ez);
        const double zm = SamplePhi(local_pos - ez);

        const double hxx = (xp - 2.0 * c + xm) * inv_h2;
        const double hyy = (yp - 2.0 * c + ym) * inv_h2;
        const double hzz = (zp - 2.0 * c + zm) * inv_h2;
        const double hxy = (SamplePhi(local_pos + ex + ey) - SamplePhi(local_pos + ex - ey) -
                            SamplePhi(local_pos - ex + ey) + SamplePhi(local_pos - ex - ey)) *
                           inv_4h2;
        const double hxz = (SamplePhi(local_pos + ex + ez) - SamplePhi(local_pos + ex - ez) -
                            SamplePhi(local_pos - ex + ez) + SamplePhi(local_pos - ex - ez)) *
                           inv_4h2;
        const double hyz = (SamplePhi(local_pos + ey + ez) - SamplePhi(local_pos + ey - ez) -
                            SamplePhi(local_pos - ey + ez) + SamplePhi(local_pos - ey - ez)) *
                           inv_4h2;

        hessian_x = ChVector3d(hxx, hxy, hxz);
        hessian_y = ChVector3d(hxy, hyy, hyz);
        hessian_z = ChVector3d(hxz, hyz, hzz);
    }

    ChSdfContactSampleQuery QueryLocal(const ChVector3d& local_pos,
                                       const ChVector3d& local_vel = ChVector3d(0, 0, 0)) const {
        ChSdfContactSampleQuery query;
        query.phi = SamplePhi(local_pos);
        query.raw_grad = SampleRawGradient(local_pos);
        query.raw_grad_norm = query.raw_grad.Length();
        query.grad = NormalizeOrFallback(query.raw_grad, ChVector3d(0, 1, 0));
        SampleHessian(local_pos, query.hessian_x, query.hessian_y, query.hessian_z);
        query.has_hessian = true;
        query.world_pos = local_pos;
        query.world_vel = local_vel;
        return query;
    }
};

inline ChSdfTriangleMeshData LoadWavefrontMeshForSdf(const std::filesystem::path& path,
                                                     double scale = 1.0,
                                                     const ChVector3d& translation = ChVector3d(0, 0, 0)) {
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
            mesh.vertices.emplace_back(translation + ChVector3d(x, y, z) * scale);
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
                mesh.faces.push_back(ChSdfTriangleFace{ids[0], ids[i], ids[i + 1]});
            }
        }
    }

    return mesh;
}

inline ChSdfTriangleMeshData TranslateMeshForSdf(const ChSdfTriangleMeshData& mesh, const ChVector3d& translation) {
    ChSdfTriangleMeshData out = mesh;
    for (auto& vertex : out.vertices) {
        vertex += translation;
    }
    return out;
}

inline ChOpenVdbSdfGrid BuildOpenVdbLevelSetFromMesh(const ChSdfTriangleMeshData& mesh,
                                                     double voxel_size,
                                                     float half_width_voxels) {
    if (voxel_size <= 0.0) {
        throw std::invalid_argument("OpenVDB SDF voxel_size must be positive.");
    }
    if (mesh.vertices.empty() || mesh.faces.empty()) {
        throw std::invalid_argument("OpenVDB SDF mesh must contain vertices and faces.");
    }

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
    auto grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
        *transform, points, triangles, half_width_voxels);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    ChOpenVdbSdfGrid out;
    out.grid = grid;
    out.voxel_size = voxel_size;
    return out;
}

inline ChSdfContactSurfaceGraph MakeSurfaceGraphFromMeshForSdf(const ChSdfTriangleMeshData& mesh) {
    return MakeTriangleMeshSurfaceGraph(mesh.vertices, mesh.faces);
}

}  // namespace sdfcontact
}  // namespace chrono

#endif
