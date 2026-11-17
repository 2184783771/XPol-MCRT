# XPol-MCRT

X-ray Polarization Monte Carlo Radiative Transfer, based on STOKES v1.2 by
Rene Goosmann and Frederic Marin. This folder contains the current simulation
core, its configuration, and a cross-platform CMake build file. Analysis,
paper-generation scripts, historical inputs, and previous outputs are excluded.

Build with CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Run from this directory so the program finds `config.ini`. On Linux the
executable is `build/XPol-MCRT`; with a Visual Studio generator on Windows it
is usually `build/Release/XPol-MCRT.exe`. The supplied configuration runs
10,000,000 photons, so lower `PhotonNum` for a quick startup test. Results are
written under `outputs/<OutputFile>/` for an unqualified output name.

For Visual Studio 2022, open this source folder as a CMake project.
Visual Studio regenerates its project from `CMakeLists.txt` and uses the
configured debugger working directory and stack setting.

Supported corona types are `slab`, `toroidal`, and `lamp-post`. The legacy
standalone `Sphere` and `Torus` scattering-region options and the misleading
`CoronaType sphere` alias have been removed. A finite lamp-post emitter still
uses spherical geometry, and the toroidal corona still uses torus intersections.

## Source layout

- `main.cpp`: minimal entry point.
- `run-simulation.cpp` and `run-state.h`: photon transport and run counters.
- `run-report.cpp`: console diagnostics and metadata CSV output.
- `model.cpp`: model defaults, validation, and accessors.
- `model-parse-basic.cpp` and `model-parse-xrb.cpp`: grouped `config.ini` options.
- `corona-lamppost.cpp` and `corona-toroidal.cpp`:
  model-specific source and boundary geometry.
- `photon.cpp` and `observe.cpp`: photon physics and detectors.
