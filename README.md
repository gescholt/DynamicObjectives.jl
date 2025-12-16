# Dynamic\_objectives

Dynamical system objective functions for global optimization and parameter estimation experiments.

## Overview

This package provides **13 ODE benchmark models** for testing parameter estimation algorithms. Each model creates objective functions that measure the error between simulated and reference ODE trajectories.

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/Dynamic_objectives")
```

## Quick Start

```julia
using Dynamic_objectives

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

## Running Tests

```bash
# All tests
./run_tests.sh

# Specific test suite
./run_tests.sh basic      # Basic tests
./run_tests.sh validate   # Validate all models
```

## Documentation

Full documentation is available in the `docs/` folder:

- **[Getting Started](docs/src/getting_started.md)** - Installation and setup
- **[Model Catalog](docs/src/model_catalog.md)** - All 13 models with configurations
- **[Testing Guide](docs/src/testing_guide.md)** - Testing strategies
- **[Integration](docs/src/integration.md)** - Globtim integration
- **[API Reference](docs/src/api_reference.md)** - Function reference

Build the documentation:

```bash
cd docs
julia --project=. make.jl
```

## License

MIT License
