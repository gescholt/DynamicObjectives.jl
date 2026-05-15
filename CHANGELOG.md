# Changelog

All notable changes to Dynamic_objectives.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-04-30

Initial public release.

### Added

- ODE-based objective functions for global optimization, built on ModelingToolkit + SciML.
- **Glued objectives** — combined trajectory-tracking objectives across multiple state variables and time points.
- **Model catalog** — Lotka-Volterra (2D, 3D, 4D), FitzHugh-Nagumo (3D), daisy world (4D), and other parametrized ODE systems.
- **Dynamic objective construction** — symbolic-to-numeric pipeline driven by `ModelingToolkit.System` definitions.
- **Solver dispatch** — Tsit5 (default), Rosenbrock23, KenCarp4, Vern9 via `OrdinaryDiffEq*` per-method packages.
- **NaN penalty handling** — robust trajectory failure modes via configurable penalty values.
- **Globtim integration** (weakdep) — `DynamicObjectivesGlobtimExt` activates when `Globtim` + `GlobtimPostProcessing` + `Optim` are loaded; provides screening helpers, post-processing pipelines, and refinement adapters.

### Notes

- `Globtim`, `GlobtimPostProcessing`, and `Optim` are weak dependencies. Load them alongside `Dynamic_objectives` to activate the integration extension.
- Julia 1.12+ required.
