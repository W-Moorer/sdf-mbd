# Milestone 01 - Chrono Baseline Contact Demo

## Overview

This milestone establishes a stable, reproducible Chrono baseline simulation using
native SMC (Smooth Contact) capabilities only. No SDF, OpenVDB, or custom collision
logic is included.

## Scenario

A sphere (radius=0.5m, mass=523.6kg) drops from 3m height onto a static ground plate.
Contact is modeled using Chrono's built-in SMC penalty method with Bullet collision
detection.

## Build

```powershell
cmake -B build -S . -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target demo_CH_baseline_contact -j 8
```

## Run

```powershell
# Default output to out/milestone_01/baseline_contact_output.csv
build\bin\Release\demo_CH_baseline_contact.exe

# Custom output path
build\bin\Release\demo_CH_baseline_contact.exe custom/path/output.csv
```

## Output

CSV file at `out/milestone_01/baseline_contact_output.csv` containing:
- time, position (x/y/z), velocity (x/y/z)
- contact count, contact force (x/y/z), contact torque (x/y/z)

## Key Results

- Final position: y = 0.4901 m (sphere resting on plate)
- Final contact force: Fy = 5136.5 N (matches gravity: 523.6 * 9.81)
- Reproducible: multiple runs produce identical output
