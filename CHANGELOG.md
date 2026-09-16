# Changelog

All notable changes to Dynamic_objectives.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-10

### Added

- **Driven FitzHugh--Nagumo family** — `define_fhn_driven_auto`,
  `define_fhn_driven_A015`, `define_fhn_driven_A030`, `define_fhn_driven3_auto`.
  Spike-count aliasing creases the loss into tens of local minima, which makes
  the family a completeness test rather than a recovery test.
- **Coupled FHN and coupled LV families** across a coupling ladder
  (`define_fhn_coupled_A030_k*`, `define_lv2d_coupled_k*`,
  `define_lv2d_sciml_coupled_k*`, and 24 randomized-wiring
  `define_lv2d_sciml_rc*` instances). At zero coupling the summed-square loss
  splits exactly, so the four-dimensional critical set is the product of the two
  factor sets and is known without a four-dimensional grid.
- **Nested DAISY and Goodwin ladders** — `define_daisy_ex3_model_5D` through
  `_8D`, `define_goodwin_oscillator_5D`/`_6D`, and the Hill-exponent dial
  `define_goodwin_oscillator_4D_hill2`/`_hill4`/`_hill6`. Each rung is an exact
  slice of the next, so cost can be measured against dimension with the
  landscape held fixed.
- **New benchmark systems**: Brusselator 2D, Michaelis--Menten 2D, MM-chain 3D,
  SIR 2D, SEIR 3D, Hindmarsh--Rose 3D, two-compartment PK 3D, and a
  locally-identifiable 1D model.
- **Grid-based interestingness scoring** — `interestingness_score`,
  `score_top_candidates`, `print_grid_score`, and `build_experiment_objective`.

### Fixed

- **`PrettyTables` pinned to 2.** Version 3 changed the `backend` keyword from
  a `Symbol` to a `Val`, which broke every table-printing path
  (`expected Symbol, got Val{:text}`). This is the failure the package CI has
  been hitting.
- Canonical initial-condition binding in the last 15 call sites; the package
  tests were not harmless about this.
- Four numerics-audit fixes covering `Dual`/`BigFloat` inputs reaching generic
  signatures.

### Changed

- **Relicensed from GPL-3.0 to MIT.** The published `LICENSE` was already MIT.
- Comments and docstrings no longer cite internal issue ids or monorepo-only
  script paths.

## [0.1.1] - 2026-08-05

### Fixed

- Pin the SciML stack to the tested majors — `SciMLBase = "2"`, `ModelingToolkit = "10"`, `OrdinaryDiffEq{Tsit5,Verner,Rosenbrock,SDIRK} = "1"`. The looser bounds let a fresh resolve pull the OrdinaryDiffEq-7-era stack (DiffEqBase 7), which **removed `verbose::Bool`** — the ODE-solving code passes `verbose = false`, so CI errored (`ArgumentError: Passing a Bool for verbose is no longer supported`) while local tests passed on the cached Manifest. Verified via fresh-resolve (DiffEqBase 6.214.1, tests pass).

### Changed

- Source-comment cleanup ahead of public registration (removed internal issue-tracker prefixes from a couple of code comments).

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
