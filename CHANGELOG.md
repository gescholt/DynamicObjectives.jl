# Changelog

All notable changes to Dynamic_objectives.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `define_fitzhugh_nagumo_ramsay_model`: FitzHugh-Nagumo in the form of Ramsay et
  al. (2007), `dR/dt = -(1/g)(V - a + bR)`, voltage observed.
- Package extensions `DynamicObjectivesPrettyTablesExt`, `DynamicObjectivesTermExt`,
  `DynamicObjectivesUnicodePlotsExt` and `DynamicObjectivesSDIRKExt` (see Changed).

### Changed

- **Breaking.** Terminal display needs the package that draws it, loaded next to
  DynamicObjectives; PrettyTables, Term and UnicodePlots are no longer
  dependencies. Without the package the function raises an error that names it
  (`display_results requires PrettyTables.jl. Add `using PrettyTables` ...`).
  - `using PrettyTables`: `display_parameters`, `display_error_metrics`,
    `display_results`, `display_gradient_analysis`, `display_quality_summary`,
    `display_degree_comparison`, `display_top_critical_points`, and the preview
    table of `display_time_series`.
  - `using Term`: `display_section`, `display_subsection`, `display_step`,
    `display_model_summary`, `display_optimization_result`,
    `display_optimization_progress`, and drawing of progress bars and spinners
    (`progress_bar`, `update_progress!`, `finish_progress!`, `with_progress`,
    `spinner`, `spin!`, `stop_spinner!`, `with_spinner`). The `ProgressBar` and
    `Spinner` types stay in the core package.
  - `using UnicodePlots`: `display_comparison` and the plots of
    `display_time_series` (`show_plot = false` needs PrettyTables only).
  - `DisplayConfig`, `set_display_config!`, `format_duration` and
    `format_critical_points_markdown` need none of them.
  - So that each function needs one package, the panel headers above tables and
    plots are gone: `display_error_metrics` and `display_time_series` print their
    heading as the table title, `display_comparison` prints a plain-text heading.
- **Breaking.** The solver names `"TRBDF2"` and `"KenCarp4"` need
  `using OrdinaryDiffEqSDIRK`; OrdinaryDiffEqSDIRK is no longer a dependency.
  Without it `_resolve_solver` (and every TOML config naming those solvers) raises
  `solver "TRBDF2" requires OrdinaryDiffEqSDIRK.jl. Add `using OrdinaryDiffEqSDIRK` ...`.
  Tsit5, AutoTsit5, Vern7, Vern9, Rosenbrock23 and Rodas5 are unchanged.
- The dependency closure drops from 204 to 171 registered packages.

### Removed

- **Breaking.** `TRBDF2` and `KenCarp4` are no longer exported; take them from
  `OrdinaryDiffEqSDIRK`.

### Documentation

- `define_fitzhugh_nagumo_3D_model`, `define_fitzhugh_nagumo_3D_model_two_outputs`
  and `define_fhn_3d_locally_id_model` have the recovery sign reversed relative to
  FitzHugh-Nagumo (`dR/dt = +(1/g)(V - a + bR)`). Their dynamics are unchanged, so
  catalogue entries and published results stay reproducible, but their docstrings
  now say that they are cubic test systems, not neuron models (at `FHN3D_1` the
  voltage decreases monotonically to about -72 by `t = 20`). The driven and
  coupled FitzHugh-Nagumo families use the standard form and are not affected.

## [0.2.0] - 2026-11-19

First registered release.

### Removed

- The campaign scripts in `examples/` (per-model screening runs, sweeps,
  diagnostics, integration templates) and `scripts/run_screening.jl`. They read
  result directories outside the package. `examples/` keeps a few self-contained
  demos (see `examples/README.md`).

### Changed

- `examples/configs/lv4d_basic.toml` reads the catalogue shipped in
  `paper/catalogue/`, and `examples/exp_candidates/lv4d.toml` writes its catalogue
  next to the config.
- The error for a `globtim_results/...` path with no `globtim_results/` directory
  above the config now says what to pass instead.
- `GlobtimPostProcessing` weak-dependency compat widened to `0.1, 0.2`.

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
