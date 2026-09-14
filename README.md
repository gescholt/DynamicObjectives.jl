# Dynamic\_objectives

Dynamical system objective functions for global optimization and parameter estimation experiments.

## Overview

This package provides **29 ODE benchmark models** for testing parameter estimation algorithms. Each model creates objective functions that measure the error between simulated and reference ODE trajectories.

## Installation

Requires **Julia 1.12** or newer.

DynamicObjectives is **not yet registered in Julia General**, so it installs by URL.
Press `]` at the Julia prompt to enter Pkg mode (backspace exits):

```julia-repl
pkg> add https://github.com/gescholt/DynamicObjectives.jl
```

That is enough for **Pipeline 1** (screening / candidate discovery). **Pipeline 2** (CP
recovery) lives in a package extension and additionally needs Globtim, which is
registered:

```julia-repl
pkg> add Globtim GlobtimPostProcessing
```

Or, for development against a local checkout:

```julia-repl
pkg> develop /path/to/DynamicObjectives
```

## Quick Start

```julia
using DynamicObjectives

# Define model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Create objective function
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

error_func = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    return_inf_on_error = true
)

# Test
@assert error_func(p_true) < 1e-6
```

## TOML-Driven Pipelines

This package provides two TOML-driven pipelines for ODE parameter estimation:

```
Pipeline 1: Candidate Discovery (screening)
  screening TOML  -->  run_screening_from_config()  -->  catalogue JSONL

Pipeline 2: CP Recovery (Globtim experiment)
  experiment TOML  -->  run_experiment_from_config()  -->  results
```

### Pipeline 1: Candidate Discovery

For a given ODE model, screening sweeps the parameter space to find viable
`p_true` values where the trajectory is bounded and the objective landscape
is well-behaved. Top candidates are saved to a catalogue JSONL file.

Screening configs live in [`examples/exp_candidates/`](examples/exp_candidates/).

```julia
using DynamicObjectives
result = run_screening_from_config("DynamicObjectives/examples/exp_candidates/lv4d.toml")
```

Available screening configs:
- `lv4d.toml` — Constrained Lotka-Volterra 4D
- `fhn3d.toml` — FitzHugh-Nagumo 3D (stiff-aware solver)
- `daisy4d.toml` — DAISY Ex3 4D

### Pipeline 2: CP Recovery

Once screening has produced a catalogue, the experiment pipeline runs Globtim
polynomial CP recovery on a selected catalogue entry.

Experiment configs live in [`examples/configs/`](examples/configs/).

```julia
using DynamicObjectives, Globtim, GlobtimPostProcessing
result = run_experiment_from_config("DynamicObjectives/examples/configs/lv4d_basic.toml")
```

Available experiment configs:
- `lv4d_basic.toml` — LV4D at radius=0.005
- `fhn3d_basic.toml` — FHN3D with AutoTsit5
- `fhn3d_tight.toml` — FHN3D with tighter tolerances
- `daisy4d_sweep.toml` — DAISY4D at radius=0.005

## Running Tests

```bash
# All tests
./run_tests.sh

# Specific test suite
./run_tests.sh basic      # Basic tests
./run_tests.sh validate   # Validate all models
```

## Documentation

- **[API Reference](https://gescholt.github.io/DynamicObjectives.jl/api.html)** — docstrings for every exported model, objective builder, and helper

The narrative documentation lives in the [Globtim docs](https://gescholt.github.io/Globtim.jl/dev/ode_parameter_estimation/), next to the optimizer that consumes these objectives:

- **[Model Catalog](https://gescholt.github.io/Globtim.jl/dev/model_catalog/)** — all 29 models with configurations and recommended testing sequences
- **[Glued Objectives](https://gescholt.github.io/Globtim.jl/dev/glued_objectives/)** — higher-dimensional problems with known critical-point sets

## License

MIT License
