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
// Milestone 2A: OpenVDB SDF asset + body binding + probe query verification
//
// This demo verifies:
//   1. OpenVDB can be built and linked in the current project
//   2. A Chrono body can be bound to an SDF asset
//   3. World -> body local -> SDF local coordinate transformation works
//   4. Phi (signed distance) and gradient queries are correct
//   5. Results match analytical sphere distance as a sanity check
//
// Output: out/milestone_02/sdf_probe_output.csv
// =============================================================================

// -- Chrono includes --
#include "chrono/physics/ChBody.h"
#include "chrono/core/ChFrame.h"

// -- OpenVDB includes (from vcpkg) --
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/math/Stencils.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace chrono;

// =============================================================================
// Minimal SDF asset object
// =============================================================================

/// A lightweight SDF asset that wraps an OpenVDB FloatGrid
/// and tracks the local frame in which the SDF is defined.
struct SDFAsset {
    openvdb::FloatGrid::Ptr grid;
    std::string name;
    double voxel_size;
    double half_width;

    SDFAsset() : voxel_size(0.0), half_width(0.0) {}
};

// =============================================================================
// Body-SDF binding
// =============================================================================

/// Represents a Chrono body bound to an SDF asset.
/// The SDF is defined in the body's local coordinate frame.
struct BodySDFBinding {
    ChBody* body;
    SDFAsset* asset;

    BodySDFBinding() : body(nullptr), asset(nullptr) {}

    /// Transform a world-space point into the SDF's local coordinate space.
    /// Uses Chrono's ChFrame API: body frame transforms world -> body local.
    openvdb::Vec3d WorldToSDFLocal(const ChVector3d& world_pt) const {
        ChFrame<double> body_frame(body->GetPos(), body->GetRot());
        ChVector3d local = body_frame.TransformPointParentToLocal(world_pt);
        return openvdb::Vec3d(local.x(), local.y(), local.z());
    }

    /// Query signed distance (phi) at a world-space point.
    /// Uses trilinear interpolation in OpenVDB index space.
    double ProbePhi(const ChVector3d& world_pt) const {
        openvdb::Vec3d sdf_local = WorldToSDFLocal(world_pt);
        openvdb::tools::GridSampler<
            openvdb::FloatGrid::TreeType,
            openvdb::tools::BoxSampler
        > sampler(
            asset->grid->tree(),
            asset->grid->transform()
        );
        return sampler.wsSample(sdf_local);
    }

    /// Query gradient at a world-space point using central finite differences.
    ChVector3d ProbeGradient(const ChVector3d& world_pt) const {
        openvdb::Vec3d sdf_local = WorldToSDFLocal(world_pt);
        openvdb::Coord ijk = asset->grid->transform().worldToIndexCellCentered(sdf_local);
        openvdb::math::GradStencil<openvdb::FloatGrid> stencil(*asset->grid);
        stencil.moveTo(ijk);
        openvdb::Vec3d grad = stencil.gradient();
        return ChVector3d(grad.x(), grad.y(), grad.z());
    }
};

// =============================================================================
// Analytical sphere signed distance (for validation)
// =============================================================================

double AnalyticSpherePhi(const ChVector3d& pt, const ChVector3d& center, double radius) {
    return (pt - center).Length() - radius;
}

ChVector3d AnalyticSphereNormal(const ChVector3d& pt, const ChVector3d& center) {
    ChVector3d d = pt - center;
    double len = d.Length();
    if (len < 1e-12) return ChVector3d(1, 0, 0);
    return d / len;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "=== Milestone 2A: OpenVDB SDF Probe Demo ===" << std::endl;

    // Initialize OpenVDB
    openvdb::initialize();

    // -- SDF geometry parameters (sphere) --
    double sphere_radius = 1.0;
    double voxel_size = 0.05;
    float half_width = 3.0f;
    double band_width = half_width * voxel_size;  // narrow band half width in world units

    // Create an OpenVDB sphere level set
    std::cout << "Creating OpenVDB sphere level set..." << std::endl;
    std::cout << "  Radius: " << sphere_radius << " m" << std::endl;
    std::cout << "  Voxel size: " << voxel_size << " m" << std::endl;
    std::cout << "  Half width: " << half_width << " voxels" << std::endl;
    std::cout << "  Band width: +/- " << band_width << " m (narrow band limit)" << std::endl;

    auto sphere_grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        static_cast<float>(sphere_radius),
        openvdb::Vec3f(0.0f, 0.0f, 0.0f),
        static_cast<float>(voxel_size),
        half_width
    );
    sphere_grid->setGridClass(openvdb::GRID_LEVEL_SET);

    // -- Construct SDF asset --
    SDFAsset sdf_asset;
    sdf_asset.grid = sphere_grid;
    sdf_asset.name = "sphere_sdf";
    sdf_asset.voxel_size = voxel_size;
    sdf_asset.half_width = half_width;

    // -- Construct a Chrono body --
    auto body = chrono_types::make_shared<ChBody>();
    body->SetName("SDF_Body");
    body->SetMass(1.0);
    body->SetFixed(true);

    // Place the body at a non-trivial pose to test coordinate transforms
    ChVector3d body_pos(2.0, 1.0, 0.5);
    ChQuaterniond body_rot = QuatFromAngleAxis(0.3, ChVector3d(0, 1, 0).GetNormalized());
    body->SetPos(body_pos);
    body->SetRot(body_rot);

    std::cout << "  Body position: (" << body_pos.x() << ", " << body_pos.y() << ", " << body_pos.z() << ")" << std::endl;
    std::cout << "  Body rotation: angle=0.3 rad around Y axis" << std::endl;

    // -- Create body-SDF binding --
    BodySDFBinding binding;
    binding.body = body.get();
    binding.asset = &sdf_asset;

    // -- Define probe points in world space --
    ChVector3d sphere_center_world = body_pos;

    struct ProbePoint {
        int id;
        ChVector3d world_pt;
    };

    std::vector<ProbePoint> probe_points;

    // Probe along X-axis through sphere center (within narrow band)
    // The sphere center in world space is body_pos.
    // Surface points in world space: body_pos + body_rot * (radius * direction)
    // We probe along the local X axis of the body (which is rotated in world space)
    ChVector3d body_local_x = body->GetRot().Rotate(ChVector3d(1, 0, 0));
    ChVector3d body_local_y = body->GetRot().Rotate(ChVector3d(0, 1, 0));
    ChVector3d body_local_z = body->GetRot().Rotate(ChVector3d(0, 0, 1));

    for (int i = -4; i <= 4; i++) {
        double offset = sphere_radius + i * (band_width / 4.0);
        ProbePoint p;
        p.id = static_cast<int>(probe_points.size());
        p.world_pt = sphere_center_world + body_local_x * offset;
        probe_points.push_back(p);
    }

    // Probe along Y direction (within narrow band)
    for (int i = -4; i <= 4; i++) {
        double offset = sphere_radius + i * (band_width / 4.0);
        ProbePoint p;
        p.id = static_cast<int>(probe_points.size());
        p.world_pt = sphere_center_world + body_local_y * offset;
        probe_points.push_back(p);
    }

    // Probe along Z direction (within narrow band)
    for (int i = -4; i <= 4; i++) {
        double offset = sphere_radius + i * (band_width / 4.0);
        ProbePoint p;
        p.id = static_cast<int>(probe_points.size());
        p.world_pt = sphere_center_world + body_local_z * offset;
        probe_points.push_back(p);
    }

    // Off-axis surface points (within narrow band)
    double d1 = sphere_radius + band_width * 0.5;
    double d2 = sphere_radius - band_width * 0.5;
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_x * d1});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_y * d1});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_z * d1});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_x * d2});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_y * d2});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_z * d2});
    // Exact surface points
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_x * sphere_radius});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_y * sphere_radius});
    probe_points.push_back({static_cast<int>(probe_points.size()), sphere_center_world + body_local_z * sphere_radius});

    std::cout << "\nProbing " << probe_points.size() << " points..." << std::endl;

    // -- Set up CSV output --
    std::string csv_path = "out/milestone_02/sdf_probe_output.csv";
    if (argc > 1) {
        csv_path = argv[1];
    }

    std::ofstream csv_file(csv_path);
    if (!csv_file.is_open()) {
        std::cerr << "Error: Cannot open output file: " << csv_path << std::endl;
        return 1;
    }

    csv_file << "point_id,world_x,world_y,world_z,local_x,local_y,local_z,"
             << "sdf_phi,analytic_phi,abs_error,in_band,"
             << "grad_x,grad_y,grad_z,"
             << "analytic_nx,analytic_ny,analytic_nz,"
             << "grad_dot_normal" << std::endl;

    // -- Run probe queries --
    double max_abs_error = 0.0;
    double sum_abs_error = 0.0;
    int in_band_count = 0;
    int grad_valid_count = 0;
    double avg_grad_normal_dot = 0.0;

    for (const auto& pp : probe_points) {
        // Get local coordinates
        openvdb::Vec3d local = binding.WorldToSDFLocal(pp.world_pt);

        // Query SDF
        double phi_sdf = binding.ProbePhi(pp.world_pt);

        // Gradient via finite difference (uses world-space coordinate)
        ChVector3d grad;
        double dx = voxel_size * 0.5;
        ChVector3d px = pp.world_pt + ChVector3d(dx, 0, 0);
        ChVector3d mx = pp.world_pt - ChVector3d(dx, 0, 0);
        ChVector3d py = pp.world_pt + ChVector3d(0, dx, 0);
        ChVector3d my = pp.world_pt - ChVector3d(0, dx, 0);
        ChVector3d pz = pp.world_pt + ChVector3d(0, 0, dx);
        ChVector3d mz = pp.world_pt - ChVector3d(0, 0, dx);
        grad = ChVector3d(
            (binding.ProbePhi(px) - binding.ProbePhi(mx)) / (2.0 * dx),
            (binding.ProbePhi(py) - binding.ProbePhi(my)) / (2.0 * dx),
            (binding.ProbePhi(pz) - binding.ProbePhi(mz)) / (2.0 * dx)
        );

        // Analytical signed distance in body local frame (sphere centered at origin)
        ChVector3d local_vec(local.x(), local.y(), local.z());
        double phi_analytic = AnalyticSpherePhi(local_vec, ChVector3d(0, 0, 0), sphere_radius);
        ChVector3d normal_analytic = AnalyticSphereNormal(local_vec, ChVector3d(0, 0, 0));

        // Check if point is within the narrow band
        bool in_band = (std::abs(phi_analytic) <= band_width);

        // Error: only compare when within band (outside band, SDF is clipped)
        double abs_error = 0.0;
        if (in_band) {
            abs_error = std::abs(phi_sdf - phi_analytic);
            sum_abs_error += abs_error;
            in_band_count++;
            if (abs_error > max_abs_error) {
                max_abs_error = abs_error;
            }
        } else {
            // Outside band, SDF returns +/- band_width; compare with clipped value
            double phi_clipped = (phi_analytic > 0) ? band_width : -band_width;
            abs_error = std::abs(phi_sdf - phi_clipped);
        }

        // Gradient vs normal dot product (only meaningful near surface)
        double grad_len = grad.Length();
        double dot_product = 0.0;
        if (grad_len > 1e-12) {
            ChVector3d grad_unit = grad / grad_len;
            dot_product = grad_unit.Dot(normal_analytic);
            if (std::abs(phi_analytic) < band_width * 0.8) {
                avg_grad_normal_dot += std::abs(dot_product);
                grad_valid_count++;
            }
        }

        // Console output for first/last few points
        if (pp.id < 3 || pp.id >= static_cast<int>(probe_points.size()) - 2) {
            std::cout << "  pt[" << pp.id << "] world=(" << pp.world_pt.x() << ","
                      << pp.world_pt.y() << "," << pp.world_pt.z() << ") "
                      << "phi_sdf=" << phi_sdf << " phi_analytic=" << phi_analytic
                      << " err=" << abs_error << " in_band=" << (in_band ? "Y" : "N") << std::endl;
        }

        // CSV output
        csv_file << pp.id << ","
                 << std::fixed
                 << pp.world_pt.x() << "," << pp.world_pt.y() << "," << pp.world_pt.z() << ","
                 << local.x() << "," << local.y() << "," << local.z() << ","
                 << phi_sdf << "," << phi_analytic << "," << abs_error << "," << (in_band ? 1 : 0) << ","
                 << grad.x() << "," << grad.y() << "," << grad.z() << ","
                 << normal_analytic.x() << "," << normal_analytic.y() << "," << normal_analytic.z() << ","
                 << dot_product
                 << std::endl;
    }

    csv_file.close();

    double mean_abs_error = (in_band_count > 0) ? (sum_abs_error / in_band_count) : 0.0;
    double avg_dot = (grad_valid_count > 0) ? (avg_grad_normal_dot / grad_valid_count) : 0.0;

    // -- Summary --
    std::cout << "\n=== Probe Results ===" << std::endl;
    std::cout << "  Points probed: " << probe_points.size() << std::endl;
    std::cout << "  Points in band: " << in_band_count << std::endl;
    std::cout << "  Max |error| (in band): " << max_abs_error << std::endl;
    std::cout << "  Mean |error| (in band): " << mean_abs_error << std::endl;
    std::cout << "  Avg |grad dot normal|: " << avg_dot << std::endl;
    std::cout << "  Output file: " << csv_path << std::endl;

    // Verification: SDF values within band should be close to analytical,
    // and gradient direction should align with analytical normal
    bool pass = (max_abs_error < 0.02 && mean_abs_error < 0.005 && avg_dot > 0.95);
    std::cout << "  Verification: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass ? 0 : 1;
}
