# Dynamic_objectives

Dynamical system objective functions for global optimization and parameter estimation experiments.

## Overview

This package provides **13 time parameter estimation benchmark problems** based on ordinary differential equation (ODE) systems. It was extracted from `globtimcore` to enable faster compilation and better modularity, and is specifically designed for testing global optimization algorithms on parameter estimation tasks.

## Key Features

- **13 ODE Benchmark Models**: Ranging from 1 to 20 parameters
- **Time Parameter Estimation**: All problems estimate ODE parameters from time series data
- **Robust Error Handling**: Returns `Inf` on failures instead of throwing errors (optimization-friendly)
- **Timeout Mechanism**: Per-evaluation timeout prevents ODE solver from stalling during global optimization
- **Multiple Distance Metrics**: L1, L2, log-L2 norms with support for custom functions
- **Flexible Data Generation**: Synthetic time series with optional noise injection
- **Full/Partial Observability**: Models with 1-4 measured outputs
- **Display Infrastructure**: Pretty, readable output for models, parameters, time series, and optimization results (NEW!)

## Included Models

### Lotka-Volterra Systems (8 models)
Predator-prey and competition models with varying complexity:
- **2D variants** (v1, v2, v3): 2 parameters, easy benchmarks
- **3D variants** (v1, v2): 3 parameters, medium difficulty
- **4D generalized**: 20 parameters, very hard (full interaction matrix)
- **4D constrained**: 4 parameters, medium-hard (skew-symmetric perturbations)

### DAISY Benchmark Models (2 models)
From the DAISY (Differential Algebra for Identifiability of SYstems) suite:
- DAISY Example 3 (4D with input): 4 parameters, medium difficulty
- DAISY Example 3 (4D no input): 4 parameters, medium difficulty

### Other Systems (3 models)
- **FitzHugh-Nagumo**: 3 parameters, hard (neuronal excitability, oscillatory)
- **Simple 1D/2D locally identifiable models**: 1-2 parameters, hard (test identifiability issues)

## Quick Start

### Basic Parameter Estimation

```julia
using Dynamic_objectives

# 1. Define an ODE system
model, params, states, outputs = define_daisy_ex3_model_4D()

# 2. Set up the parameter estimation problem
p_true = [0.1, 0.2, 0.3, 0.4]        # True parameters to recover
ic = [1.0, 2.0, 1.0, 1.0]            # Initial conditions
time_interval = [0.0, 10.0]          # Simulation time
numpoints = 25                        # Number of time points

# 3. Create objective function for optimization
error_func = make_error_distance(
    model,                            # ODE system
    outputs,                          # Measured variables
    ic,                               # Initial conditions
    p_true,                           # True parameters (generates reference data)
    time_interval,                    # Time interval
    numpoints,                        # Number of sampling points
    L2_norm,                          # Distance metric
    first,                            # Aggregation (for multi-output)
    nothing;                          # Noise function (none)
    return_inf_on_error = true,       # Return Inf on ODE solver failure
    eval_timeout = 10.0               # Timeout per evaluation (seconds)
)

# 4. Test the objective function
@assert error_func(p_true) < 1e-6    # Should be ~0 at true parameters
p_test = p_true .+ 0.1
@assert error_func(p_test) > 0        # Should be >0 at perturbed parameters

# 5. Use with global optimizer (e.g., globtim)
# The optimizer will minimize error_func to recover p_true
```

## Ready-to-Run Examples

See **[examples/](examples/)** directory for complete testing campaign scripts:

- **[quick_start.jl](examples/quick_start.jl)** ⭐ - Your first test (start here!)
- **[single_test_template.jl](examples/single_test_template.jl)** - Template for testing any model
- **[batch_test_campaign.jl](examples/batch_test_campaign.jl)** - Run multiple models systematically
- **[test_configs.jl](examples/test_configs.jl)** - Pre-configured test cases for all 13 models

```bash
# Run your first test
julia --project=. examples/quick_start.jl
```

See **[examples/README.md](examples/README.md)** for complete usage guide.

## Timeout Mechanism

A key feature for global optimization: the `eval_timeout` parameter prevents the ODE solver from stalling on difficult parameter regions.

```julia
error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, numpoints;
    eval_timeout = 5.0  # Kill evaluation after 5 seconds
)
```

**When to use timeout:**
- High-dimensional problems (e.g., 20D Lotka-Volterra)
- Oscillatory systems (e.g., FitzHugh-Nagumo)
- Stiff ODEs or unbounded parameters
- Any problem where optimizer explores extreme parameter values

**How it works:**
- Each objective function evaluation runs in an async task
- If evaluation exceeds `eval_timeout`, returns `Inf` (guides optimizer away)
- Prevents entire optimization from hanging on bad parameter sets

## Distance Functions

Three built-in distance metrics:

```julia
# L2 norm (Euclidean distance) - DEFAULT
error_L2 = make_error_distance(..., L2_norm, ...)

# L1 norm (Manhattan distance, scaled by 100)
error_L1 = make_error_distance(..., L1_norm, ...)

# log-L2 norm (log-scale, good for wide-range errors)
error_log = make_error_distance(..., log_L2_norm, ...)

# Custom distance function
my_distance(y_true, y_test) = maximum(abs.(y_true - y_test))
error_custom = make_error_distance(..., my_distance, ...)
```

### Pretty Display Infrastructure

```julia
# Display model information
display_model_summary(model)

# Display parameters in a formatted table
param_names = [:α, :β, :γ, :δ]
display_parameters(param_names, [p_true p_test], labels=["True", "Test"])

# Generate and display time series data with plots
problem = ODEProblem(...)
data = sample_data(problem, model, outputs, time_interval, p_true, ic, 25)
display_time_series(data, show_plot=true)

# Compare different parameter sets
data_test = sample_data(problem, model, outputs, time_interval, p_test, ic, 25)
display_comparison(data, data_test, labels=["True", "Test"])

# Display error metrics
errors = Dict("L1" => 0.123, "L2" => 0.045, "log_L2" => -3.45)
display_error_metrics(errors)

# Display optimization results
result = (params=p_optimal, error=0.001, iterations=150, converged=true, time_elapsed=2.5)
display_optimization_result(result)
```

See `examples/display_demo.jl` for a comprehensive demonstration.

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

### Quick Test Script (Recommended)

```bash
# Run all tests (basic + validation)
./run_tests.sh

# Or run specific test suite
./run_tests.sh basic      # Basic tests only
./run_tests.sh validate   # Validate all 13 models
```

### Manual Testing

```julia
# Basic test suite
using Pkg
Pkg.test("Dynamic_objectives")

# Comprehensive validation of all 13 models
include("test/validate_all_models.jl")
```

## Testing with globtim

For comprehensive testing campaigns with your local globtim optimizer, see the **[TESTING_GUIDE.md](TESTING_GUIDE.md)** which includes:

- **Complete model catalog** with all 13 models
- **Recommended parameter bounds** for each model
- **Initial conditions** and time intervals
- **Difficulty ratings** (Easy/Medium/Hard/Very Hard)
- **Testing protocol** for systematic evaluation
- **Performance metrics** to track
- **Four-phase testing campaign** (Easy → Medium → Hard → Very Hard)
- **Expected results** for identifiable vs. non-identifiable systems
- **Debugging tips** and troubleshooting

### Quick Testing Recommendations

**Start with (Phase 1 - Easy):**
1. `define_lotka_volterra_2D_model_v3_two_outputs()` - 2 params, 2 outputs
2. `define_lotka_volterra_2D_model()` - 2 params, 1 output

**Progress to (Phase 2 - Medium):**
3. `define_lotka_volterra_3D_model_v2()` - 3 params
4. `define_daisy_ex3_model_4D_no_input()` - 4 params

**Challenge with (Phase 3 - Hard):**
5. `define_fitzhugh_nagumo_3D_model()` - 3 params, oscillatory
6. `define_simple_2D_model_locally_identifiable_square()` - 2 params, non-identifiable

**Stress test (Phase 4 - Very Hard):**
7. `define_generalized_lotka_volterra_4D()` - 20 params!

See [TESTING_GUIDE.md](TESTING_GUIDE.md) for complete details.

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
