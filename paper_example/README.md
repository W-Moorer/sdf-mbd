# Paper Examples

This folder locks the paper benchmark inputs and regression checks for the current sparse-SDF contact backend.

The examples are self-contained here:

- `cases/eccentric_roller`
- `cases/headon_spheres`
- `cases/headon_spheres_mass_ratio`
- `cases/simple_gear`

The executable source is `paper_examples_openvdb.cpp`. It uses the project backend (`chrono::fieldcontact` plus OpenVDB sparse narrow-band SDF) and writes:

- `out/paper_example_dynamic_benchmarks/sparse_sdf_frames.csv`
- `out/paper_example_dynamic_benchmarks/comparison_summary.csv`
- `out/paper_example_dynamic_benchmarks/timing_summary.csv`
- `paper_example/figures/*.pdf`

Run from the repository root:

```powershell
cmake --build build --config Release --target demo_CH_paper_examples_openvdb
build\bin\Release\demo_CH_paper_examples_openvdb.exe
python paper_example\plot_paper_examples.py --project-root .
python paper_example\check_paper_examples.py --project-root .
```

When field-contact regression tests are enabled in CMake, CTest runs the executable and then validates the summary against `manifest.json`.

The `simple_gear` case locks the RecurDyn RMD, GEAR21/GEAR22 OBJ surfaces, and the GEAR22 RX reference curve. Its regression checks verify the final driven-wheel angular speed, post-startup RMS error, and post-startup adjacent-sample jump for the bidirectional sparse-SDF contact response.
