# AGENTS.md

This file defines the default operating guide for coding agents working in this repository.

## Project Summary
- Language: C++20
- Build system: CMake
- Main target: `StructureEstimationC`
- Domain: structural estimation / simulation with large CSV inputs

## Repository Map
- `main.cpp`: top-level entry file in repo root
- `archive/`: current location of most project source files (`*.cpp`, `*.h`)
- `StructuralData_InflationAdj_NormalizeAll/`: large input datasets (CSV)
- `ThirdParty/`: vendored dependencies (Eigen, Boost headers, pcg, libigl, etc.)

## Build And Run
- Configure from repo root:
  - `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
- Build:
  - `cmake --build build -j`
- Run:
  - `./build/StructureEstimationC`

## Dependency Notes
- `ipopt` is required through `pkg-config` (`pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)`).
- Expected include paths:
  - Eigen: `$EIGEN3_PATH` or `ThirdParty/eigen-3.4.0`
  - Boost: `$boost_PATH` or `ThirdParty/boost_1_78_0`
  - PCG: `$PCG_PATH/include` or `ThirdParty/pcg-cpp-0.98/include`
- `PKG_CONFIG_PATH` may need to include `$HOME/.local/lib/pkgconfig`.

## Important Caveat (Current Layout)
- `CMakeLists.txt` currently lists many source files as if they were in repo root.
- Most of those files are currently in `archive/`.
- If a build fails with missing source files, update `CMakeLists.txt` paths (or restore source locations) before other changes.

## Editing Rules For Agents
- Do not modify files in `ThirdParty/` unless explicitly requested.
- Do not rewrite, reformat, or mass-edit files in `StructuralData_InflationAdj_NormalizeAll/`.
- Keep changes minimal and scoped to requested behavior.
- Preserve existing numerical logic and parameter ordering unless the task explicitly asks for model changes.

## Validation Expectations
- For code changes, run at least:
  - `cmake --build build -j`
- If runtime behavior is touched, run:
  - `./build/StructureEstimationC`
- Report any dependency or environment blockers clearly (especially `ipopt`).

## Git Hygiene
- Keep commits focused (one logical change per commit when possible).
- Avoid committing generated build artifacts (`build/`).
- Avoid committing IDE state changes unless explicitly requested.
