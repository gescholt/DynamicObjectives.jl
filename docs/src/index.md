# Dynamic\_objectives Documentation

**Dynamical system objective functions for global optimization and parameter estimation experiments.**

## The Problem

Parameter estimation in dynamical systems is fundamentally hard. Given observed time series data from an ODE system, how do you find the parameters that best explain the data? Standard optimization algorithms can get stuck in local minima, and the objective function landscape for ODE fitting is often highly multimodal.

## The Approach

Dynamic\_objectives provides **13 benchmark ODE models** specifically designed for testing global optimization algorithms on parameter estimation tasks. Each model:

1. **Defines an ODE system** — Using ModelingToolkit.jl for symbolic representation
2. **Generates synthetic data** — Time series from true parameters with optional noise
3. **Creates objective functions** — Returns error between simulated and reference trajectories
4. **Handles failures gracefully** — Returns `Inf` on solver failures (optimization-friendly)

The result: a systematic way to benchmark parameter estimation algorithms across problems of varying difficulty.

## Key Features

### 13 ODE Benchmark Models

| Difficulty | Models | Parameters |
|------------|--------|------------|
| Easy | Lotka-Volterra 2D (v1, v2, v3) | 2 params |
| Medium | Lotka-Volterra 3D, DAISY benchmarks | 3-4 params |
| Hard | FitzHugh-Nagumo, Identifiability tests | 2-3 params |
| Very Hard | Generalized Lotka-Volterra 4D | 20 params |

See [Model Catalog](model_catalog.md) for complete details.

### Robust Error Handling

- **Returns `Inf` on failures** — No exceptions during optimization
- **Timeout mechanism** — Prevents ODE solver from stalling
- **Multiple distance metrics** — L1, L2, log-L2 norms

### Flexible Data Generation

- **Synthetic time series** with configurable sampling
- **Optional noise injection** for robustness testing
- **Full/partial observability** — 1-4 measured outputs

### Display Infrastructure

Pretty, readable output for models, parameters, and results:

```julia
display_model_summary(model)
display_parameters(param_names, values)
display_time_series(data)
```

## Installation

This package is a local development package:

```julia
using Pkg
Pkg.develop(path="/path/to/Dynamic_objectives")
```

Or add to your `Project.toml`:

```toml
[deps]
Dynamic_objectives = "a7860e71-583f-4874-99b9-9df52d7a8e5c"

[sources]
Dynamic_objectives = {path = "../Dynamic_objectives"}
```

## Quick Start

```julia
using Dynamic_objectives

# 1. Define an ODE system
model, params, states, outputs = define_lotka_volterra_2D_model_v3()

# 2. Set up the parameter estimation problem
p_true = [1.0, 1.0]
ic = [1.0, 1.0]
time_interval = [0.0, 10.0]
numpoints = 25

# 3. Create objective function
error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, numpoints,
    L2_norm, first, nothing;
    return_inf_on_error = true,
    eval_timeout = 10.0
)

# 4. Test at true parameters
@assert error_func(p_true) < 1e-6  # Should be ~0
```

For more details, see [Getting Started](getting_started.md).

## Contents

```@contents
Pages = ["getting_started.md", "architecture.md", "model_catalog.md", "testing_guide.md", "integration.md", "postprocessing_requirements.md", "api_reference.md", "performance_report.md"]
Depth = 2
```
