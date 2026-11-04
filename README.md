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

For Visual Studio 2022, open `vs2022/XPolMCRT.sln`. This solution is generated
from `CMakeLists.txt` and includes the required debugger working directory and
stack setting. The root `XPol-MCRT.sln` was created by the existing-code wizard
and does not inherit these settings. `vs2022/` is a local build directory and
should not be committed; after cloning, open the source folder as a CMake
project or regenerate the solution.

## Source layout

- `main-stokes-v1.2.cpp`: minimal entry point.
- `run-simulation.cpp` and `run-state.h`: photon transport and run counters.
- `run-report.cpp`: console diagnostics and metadata CSV output.
- `model-stokes-v1.2.cpp`: model defaults, validation, and accessors.
- `model-parse-basic.cpp` and `model-parse-xrb.cpp`: grouped `config.ini` options.
- `photon-stokes-v1.2.cpp` and `observe-stokes-v1.2.cpp`: photon physics and detectors.
