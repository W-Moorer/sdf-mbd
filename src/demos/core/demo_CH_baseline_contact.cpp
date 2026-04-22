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
// Baseline headless contact demo: sphere dropping onto a static ground plate
// using Chrono native SMC contact (no SDF, no OpenVDB, no custom collision).
//
// Output: CSV file with timestep, position, velocity, contact count, contact force
// =============================================================================

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChContactContainerSMC.h"
#include "chrono/physics/ChBodyEasy.h"

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace chrono;

int main(int argc, char* argv[]) {
    std::cout << "=== Chrono Baseline Contact Demo ===" << std::endl;
    std::cout << "Chrono version: " << CHRONO_VERSION << std::endl;
    std::cout << "Contact method: SMC (Smooth Contact)" << std::endl;
    std::cout << "Collision system: Bullet (native)" << std::endl;

    // Simulation parameters
    double time_step = 5e-4;
    double total_time = 2.0;
    double output_interval = 0.01;

    // Physical parameters
    double gravity = -9.81;
    double sphere_radius = 0.5;
    double sphere_density = 1000.0;
    double drop_height = 3.0;
    double friction = 0.3f;
    double restitution = 0.1f;
    double young_modulus = 1e7f;

    // Create the SMC physical system
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, gravity, 0));

    // Use Bullet collision system (native Chrono capability)
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Solver settings for stability
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(100);
    sys.GetSolver()->AsIterative()->SetTolerance(1e-6);

    // Contact material
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetFriction(friction);
    mat->SetRestitution(restitution);
    mat->SetYoungModulus(young_modulus);
    mat->SetPoissonRatio(0.3f);

    // Set default effective curvature radius for SMC
    ChCollisionInfo::SetDefaultEffectiveCurvatureRadius(sphere_radius);

    // Create static ground plate
    auto ground = chrono_types::make_shared<ChBodyEasyBox>(
        10.0, 0.2, 10.0,    // x, y, z size
        1000.0,             // density
        false,              // no visualization
        true,               // enable collision
        mat                 // contact material
    );
    ground->SetPos(ChVector3d(0, -0.1, 0));
    ground->SetFixed(true);
    ground->EnableCollision(true);
    sys.AddBody(ground);

    // Create falling sphere
    auto sphere = chrono_types::make_shared<ChBodyEasySphere>(
        sphere_radius,       // radius
        sphere_density,      // density
        false,               // no visualization
        true,                // enable collision
        mat                  // contact material
    );
    sphere->SetPos(ChVector3d(0, drop_height, 0));
    sphere->SetFixed(false);
    sphere->EnableCollision(true);
    sys.AddBody(sphere);

    // Print simulation info
    std::cout << "\nSimulation setup:" << std::endl;
    std::cout << "  Sphere radius: " << sphere_radius << " m" << std::endl;
    std::cout << "  Sphere mass: " << sphere->GetMass() << " kg" << std::endl;
    std::cout << "  Drop height: " << drop_height << " m" << std::endl;
    std::cout << "  Time step: " << time_step << " s" << std::endl;
    std::cout << "  Total time: " << total_time << " s" << std::endl;
    std::cout << "  Output interval: " << output_interval << " s" << std::endl;

    // Open CSV output file
    std::string csv_filename = "out/milestone_01/baseline_contact_output.csv";

    // Allow overriding output path via command line argument
    if (argc > 1) {
        csv_filename = argv[1];
    }

    std::ofstream csv_file(csv_filename);
    if (!csv_file.is_open()) {
        std::cerr << "Error: Cannot open output file: " << csv_filename << std::endl;
        return 1;
    }

    // CSV header
    csv_file << "time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,contact_count,contact_force_x,contact_force_y,contact_force_z,contact_torque_x,contact_torque_y,contact_torque_z" << std::endl;

    // Simulation loop
    double next_output_time = 0.0;
    int step_count = 0;
    int output_count = 0;
    bool simulation_ok = true;

    std::cout << "\nRunning simulation..." << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;

    while (sys.GetChTime() < total_time) {
        sys.DoStepDynamics(time_step);
        step_count++;

        double time = sys.GetChTime();

        if (time >= next_output_time - 1e-10) {
            next_output_time += output_interval;
            output_count++;

            // Get sphere state
            ChVector3d pos = sphere->GetPos();
            ChVector3d vel = sphere->GetPosDt();
            ChVector3d cforce = sphere->GetContactForce();
            ChVector3d ctorque = sphere->GetContactTorque();

            // Get contact count from the system
            int contact_count = static_cast<int>(sys.GetNumContacts());

            // Console output (summary)
            if (output_count % 10 == 1 || time >= total_time - 1e-10) {
                std::cout << std::fixed;
                std::cout.precision(4);
                std::cout << "t=" << time
                          << " y=" << pos.y()
                          << " vy=" << vel.y()
                          << " contacts=" << contact_count
                          << " Fy=" << cforce.y()
                          << std::endl;
            }

            // CSV output (full data)
            csv_file << std::fixed;
            csv_file.precision(8);
            csv_file << time << ","
                     << pos.x() << "," << pos.y() << "," << pos.z() << ","
                     << vel.x() << "," << vel.y() << "," << vel.z() << ","
                     << contact_count << ","
                     << cforce.x() << "," << cforce.y() << "," << cforce.z() << ","
                     << ctorque.x() << "," << ctorque.y() << "," << ctorque.z()
                     << std::endl;

            // Check for divergence
            if (std::abs(pos.y()) > 100.0 || std::isnan(pos.y()) || std::isinf(pos.y())) {
                std::cerr << "Warning: Simulation may have diverged at t=" << time << std::endl;
                simulation_ok = false;
                break;
            }
        }
    }

    csv_file.close();

    // Final summary
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "Simulation completed." << std::endl;
    std::cout << "  Total steps: " << step_count << std::endl;
    std::cout << "  Outputs written: " << output_count << std::endl;
    std::cout << "  Final time: " << sys.GetChTime() << " s" << std::endl;
    std::cout << "  Final sphere Y position: " << sphere->GetPos().y() << " m" << std::endl;
    std::cout << "  Final sphere Y velocity: " << sphere->GetPosDt().y() << " m/s" << std::endl;
    std::cout << "  Output file: " << csv_filename << std::endl;

    if (simulation_ok) {
        std::cout << "  Status: SUCCESS" << std::endl;
    } else {
        std::cout << "  Status: WARNING - Possible divergence detected" << std::endl;
    }

    return simulation_ok ? 0 : 1;
}
