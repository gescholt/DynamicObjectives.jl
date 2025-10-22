# Dynamic_objectives

Dynamical system objective functions for global optimization experiments.

## Overview

This package provides objective functions based on dynamical systems for parameter estimation and global optimization. It was extracted from `globtimcore` to enable faster compilation and better modularity.

## Features

- **ODE System Definitions**: Lotka-Volterra models, DAISY benchmark models, FitzHugh-Nagumo, and test models for identifiability studies
- **Error Metrics**: L1, L2, log-L2 distance functions
- **Data Generation**: Synthetic time series generation from ODE systems
- **Flexible Objective Functions**: Configurable error functions with timeout and noise injection support

## Included Models

### Lotka-Volterra Systems
- 2D variants (v1, v2, v3)
- 3D variants (v1, v2)
- 4D generalized (20 parameters)
- 4D constrained with skew-symmetric perturbations (4 parameters)

### DAISY Benchmark Models
- DAISY Example 3 (4D with/without input)

### Other Systems
- FitzHugh-Nagumo (neuronal excitability)
- Simple 1D/2D locally identifiable test models

## Usage

```julia
using Dynamic_objectives

# Define a 4D Lotka-Volterra system
model, params, states, outputs = define_daisy_ex3_model_4D()

# Create an objective function
p_true = [0.1, 0.2, 0.3, 0.4]
ic = [1.0, 2.0, 1.0, 1.0]
time_interval = [0.0, 10.0]

error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, 25
)

# Evaluate at test parameters
p_test = [0.15, 0.25, 0.35, 0.45]
error = error_func(p_test)
```

## Installation

This package is currently a local development package. To use it:

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

## Testing

```julia
using Pkg
Pkg.test("Dynamic_objectives")
```

## Migration from globtimcore

If you were previously using:
```julia
include(joinpath(PROJECT_ROOT, "Examples", "systems", "DynamicalSystems.jl"))
```

Replace with:
```julia
using Dynamic_objectives
```

All functions remain the same, but now:
- ✅ Precompiled (faster startup)
- ✅ No world age warnings
- ✅ Better organized
- ✅ Independently testable

## License

MIT License
