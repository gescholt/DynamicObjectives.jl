# globtim Integration Guide

Integration of Dynamic_objectives parameter estimation benchmarks with the globtimcore global optimization engine.

## Overview

This integration allows you to use globtim's polynomial approximation-based optimization to solve parameter estimation problems from Dynamic_objectives. The integration provides:

- **Signature compatibility**: Adapts Dynamic_objectives error functions to globtim's interface
- **Complete workflow**: Helper functions for running experiments and analyzing results
- **Comprehensive testing**: Formal test suite with validation criteria
- **Batch processing**: Scripts for testing multiple models systematically

## Architecture

### Package Dependencies

```
Dynamic_objectives (this package)
    ├── Defines ODE models and error metrics
    ├── Depends on: ModelingToolkit, OrdinaryDiffEq
    └── Integration module: src/globtim_integration.jl
            ├── Wrapper functions
            └── Depends on: Globtim (globtimcore)

globtimcore (Globtim package)
    ├── Polynomial approximation
    ├── Critical point solving (HomotopyContinuation)
    └── BFGS refinement
```

### Key Signature Difference

**Dynamic_objectives** creates single-argument error functions:
```julia
error_func(p::Vector{Float64}) -> Float64
```

**globtimcore** expects two-argument objective functions:
```julia
objective(point::Vector{Float64}, params) -> Float64
```

The integration wrapper bridges this gap using closures that capture model configuration.

## Quick Start

### 1. Basic Usage - Single Model

```julia
using Dynamic_objectives
using Globtim

# Define a model (simplest: 2D Lotka-Volterra with 2 outputs)
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Configuration
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Create globtim-compatible objective
objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0],  # time_interval
    30            # numpoints
)

# Run optimization
result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 8,
    degree_range = 4:8,
    output_dir = "results/lv2d_test"
)

# Check results
if result[:success_rate] > 0.0
    println("Best parameters: ", result[:best_params])
    println("Recovery error: ", result[:recovery_error])
end
```

### 2. Using Example Scripts

Run a single model test:
```bash
julia --project=. examples/globtim_single_test.jl
```

Run batch tests:
```bash
# Test EASY models only (4 models, 2D)
julia --project=. examples/globtim_batch_test.jl easy

# Test MEDIUM models only (5 models, 3-4D)
julia --project=. examples/globtim_batch_test.jl medium

# Test all (default)
julia --project=. examples/globtim_batch_test.jl all
```

Analyze results:
```bash
julia --project=. examples/analyze_globtim_results.jl
```

### 3. Running Tests

```bash
# Run full test suite (including globtim integration tests)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run only integration tests
julia --project=. test/test_globtim_integration.jl
```

## API Reference

### `create_globtim_objective`

Create a globtim-compatible objective function from a Dynamic_objectives model.

**Signature:**
```julia
create_globtim_objective(
    model, outputs, ic, p_true, time_interval, numpoints;
    distance = L2_norm,
    aggregate = first,
    eval_timeout = nothing,
    return_inf_on_error = true
) -> Function
```

**Arguments:**
- `model::ODESystem`: ModelingToolkit ODE system
- `outputs::Vector{Equation}`: Measured quantities
- `ic::Vector{Float64}`: Initial conditions
- `p_true::Vector{Float64}`: True parameters (generates reference data)
- `time_interval`: `[t_start, t_end]`
- `numpoints::Int`: Number of time samples

**Keyword Arguments:**
- `distance`: Distance function (`L2_norm`, `L1_norm`, `log_L2_norm`)
- `aggregate`: How to combine multi-output errors (default: `first`)
- `eval_timeout`: Timeout in seconds (prevents ODE hanging)
- `return_inf_on_error`: Return `Inf` on ODE failure (for optimization)

**Returns:**
- Function with signature `f(point::Vector{Float64}, params) -> Float64`

**Example:**
```julia
model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

objective = create_globtim_objective(
    model, outputs, [1.0, 0.5], [1.0, 0.5],
    [0.0, 20.0], 30;
    eval_timeout = 10.0  # Use for oscillatory models
)

# Use with globtim
result = run_standard_experiment(
    objective_function = objective,
    problem_params = (;),
    domain_bounds = [(0.0, 3.0), (0.0, 2.0)],
    ...
)
```

### `run_globtim_optimization`

Complete workflow for optimizing a Dynamic_objectives model with globtim.

**Signature:**
```julia
run_globtim_optimization(
    objective, bounds, p_true;
    GN = 8,
    degree_range = 4:8,
    output_dir = "globtim_results",
    model_name = "model",
    max_time = 3600,
    basis = :chebyshev
) -> Dict
```

**Arguments:**
- `objective::Function`: globtim-compatible objective `f(point, params) -> Float64`
- `bounds::Vector{Tuple{Float64, Float64}}`: Parameter bounds
- `p_true::Vector{Float64}`: True parameters (for recovery error)

**Keyword Arguments:**
- `GN::Int`: Grid points per dimension (total: `GN^dim`)
- `degree_range`: Range of polynomial degrees to test
- `output_dir::String`: Directory for output files
- `model_name::String`: Name for metadata
- `max_time::Int`: Maximum time in seconds
- `basis::Symbol`: Polynomial basis (`:chebyshev` or `:legendre`)

**Returns:**
Dictionary with:
- `:success_rate`: Fraction of degrees that succeeded
- `:degree_results`: Array of results for each degree
- `:total_time`: Total optimization time
- `:best_params`: Best parameters across all degrees
- `:best_objective`: Best objective value
- `:recovery_error`: `||p_best - p_true||`
- `:total_critical_points`: Total critical points found

## Configuration Guidelines

### Grid Size (`GN`)

Grid size determines polynomial approximation quality and computational cost.

| Dimension | Conservative | Standard | Thorough |
|-----------|--------------|----------|----------|
| **2D** | GN=6 (36 evals) | GN=8 (64 evals) | GN=12 (144 evals) |
| **3D** | GN=5 (125 evals) | GN=7 (343 evals) | GN=10 (1000 evals) |
| **4D** | GN=5 (625 evals) | GN=8 (4096 evals) | GN=12 (20736 evals) |
| **20D** | GN=3 (3.5M evals) | GN=4 (1.1e12 evals) | Not feasible |

**Recommendations:**
- Start with `GN=8` for 2D-4D problems
- Use `GN=4-5` for high-dimensional problems (>10D)
- Increase GN if solution quality is poor

### Polynomial Degree (`degree_range`)

Controls polynomial complexity. globtim tests multiple degrees automatically.

**Recommendations:**
- **Standard**: `4:8` (good balance for most problems)
- **Fast iteration**: `4:6` (fewer degrees, faster)
- **High accuracy**: `6:10` (more degrees, better fit)

Higher degrees:
- ✅ Better capture complex landscapes
- ❌ More critical points to solve (slower)
- ❌ Potential for overfitting on noisy data

### Timeout (`eval_timeout`)

Prevents ODE solver from hanging on difficult parameter regions.

**When to use:**
- **Always use** for oscillatory models (FitzHugh-Nagumo)
- **Optional** for well-behaved models (Lotka-Volterra)
- **Not needed** for linear models (DAISY benchmarks)

**Values:**
- `nothing`: No timeout (default for EASY models)
- `5.0`: Tight timeout for fast iteration
- `10.0`: Standard timeout for MEDIUM models
- `15.0`: Loose timeout for HARD models

**Warning:** Too tight a timeout causes many `Inf` returns, degrading polynomial approximation.

### Distance Metrics

Choose based on problem characteristics:

| Metric | Formula | Use Case |
|--------|---------|----------|
| **L2_norm** | `sqrt(sum((y_true - y_test).^2))` | Default, standard least-squares |
| **L1_norm** | `sum(abs.(y_true - y_test))` | Robust to outliers |
| **log_L2_norm** | `log(1 + L2_norm(...))` | Wide dynamic range (e.g., exponential growth) |

## Model Catalog

### EASY Models (2D, >75% expected success)

| Model | Function | Parameters | Outputs | Recommended Bounds |
|-------|----------|------------|---------|-------------------|
| LV 2D v1 | `define_lotka_volterra_2D_model()` | 2 | 1 | `[(-1, 2), (-1, 0)]` |
| LV 2D v2 | `define_lotka_volterra_2D_model_v2()` | 2 | 1 | `[(-1, 2), (-1, 0)]` |
| LV 2D v3 | `define_lotka_volterra_2D_model_v3()` | 2 | 1 | `[(0, 3), (0, 2)]` |
| LV 2D v3 (2 out) | `define_lotka_volterra_2D_model_v3_two_outputs()` | 2 | 2 | `[(0, 3), (0, 2)]` |

**Start here**: LV 2D v3 with 2 outputs (full observability, easiest)

### MEDIUM Models (3-4D, >50% expected success)

| Model | Function | Parameters | Recommended Config |
|-------|----------|------------|-------------------|
| LV 3D v2 | `define_lotka_volterra_3D_model_v2()` | 3 | GN=8, deg=4:8 |
| DAISY Ex3 (no input) | `define_daisy_ex3_model_4D_no_input()` | 4 | GN=8, deg=4:8 |
| LV 4D Constrained | `define_constrained_lotka_volterra_4D()` | 4 | GN=6-8, deg=4:8 |

### HARD Models (Challenging)

| Model | Function | Parameters | Challenge | Config |
|-------|----------|------------|-----------|--------|
| FitzHugh-Nagumo | `define_fitzhugh_nagumo_3D_model()` | 3 | Oscillatory | timeout=10.0, GN=8 |
| Simple 2D square | `define_simple_2D_model_locally_identifiable_square()` | 2 | Non-identifiable | deg=6:10 |

### VERY HARD (Stress Test)

| Model | Function | Parameters | Challenge | Config |
|-------|----------|------------|-----------|--------|
| LV 4D Generalized | `define_generalized_lotka_volterra_4D()` | 20 | High-dimensional | GN=4, deg=4:6 |

## Troubleshooting

### Problem: Many grid points return `Inf`

**Symptoms:**
- Polynomial approximation fails
- No critical points found
- All degrees fail

**Causes:**
1. Bounds too wide (ODE solver fails outside valid region)
2. Timeout too strict (returns `Inf` too often)
3. Model is stiff or unstable

**Solutions:**
- Tighten bounds based on prior knowledge
- Increase or disable timeout
- Try different ODE solver tolerances
- Check model definition for errors

### Problem: No solution recovers true parameters

**Symptoms:**
- globtim finds critical points
- Best solution has high recovery error (>0.1)
- Objective value is non-zero

**Causes:**
1. Non-identifiable problem (multiple parameter sets give same output)
2. Local minima trap
3. Polynomial approximation quality insufficient

**Solutions:**
- Check if problem is identifiable (read model description)
- Increase `GN` for better approximation
- Try more polynomial degrees (`degree_range = 6:10`)
- Check if true parameters are within bounds

### Problem: Very slow for 4D problems

**Symptoms:**
- 4D models take >30 minutes
- Grid evaluation takes most of the time

**Causes:**
- `GN=8` → 4096 evaluations
- ODE integration is expensive
- HomotopyContinuation solving many critical points

**Solutions:**
- Reduce `GN` to 6 (1296 evals) or 5 (625 evals)
- Reduce `degree_range` to 4:6
- Enable timeout to prevent slow ODE solves
- Use fewer time points (`numpoints = 20` instead of 50)

### Problem: Tests fail with package not found

**Symptoms:**
```julia
ERROR: ArgumentError: Package Globtim not found in current path
```

**Solution:**
The integration requires globtimcore to be available. Add it to the environment:

```julia
using Pkg
# If globtimcore is a local package:
Pkg.develop(path="../globtimcore")

# If it's registered:
Pkg.add("Globtim")
```

## Performance Expectations

### EASY Models (2D)
- **Success rate**: >75%
- **Time**: <60 seconds per model
- **Recovery error**: <0.01 for successful runs
- **Function evaluations**: <500

### MEDIUM Models (3-4D)
- **Success rate**: >50%
- **Time**: 1-5 minutes per model
- **Recovery error**: <0.1 for successful runs
- **Function evaluations**: <2000

### HARD Models
- **Success rate**: >30% (expected)
- **Time**: 5-30 minutes
- **May not fully converge** (local minima, non-identifiability)

## Output Files

globtim creates the following outputs in the specified `output_dir`:

### `results_summary.json`
JSON file with complete experiment metadata:
- Model configuration
- globtim settings
- Results for each polynomial degree
- Best solution found
- Recovery error (if `p_true` provided)

### `critical_points_deg_N.csv`
CSV files (one per degree tested) containing all critical points found:
- Parameter values
- Objective function value
- Status (local min, saddle, etc.)

### Example Analysis

```julia
using JSON3
using CSV
using DataFrames

# Load summary
summary = JSON3.read(read("results/lv2d/results_summary.json", String))

# Load critical points for degree 6
cp = CSV.read("results/lv2d/critical_points_deg_6.csv", DataFrame)

# Find best point
best_idx = argmin(cp.objective)
best_params = Vector(cp[best_idx, 1:end-1])  # Exclude objective column
```

## Integration with globtimpostprocessing

For advanced analysis, use the globtimpostprocessing package:

```julia
using GlobtimPostProcessing

# Load experiment results
experiment = load_experiment_results("results/lv2d")

# Compute statistics
stats = compute_statistics(experiment)

# Analyze critical points
analysis = analyze_critical_points(experiment)

# Generate report
report = generate_report(experiment, stats)
```

See `examples/analyze_globtim_results.jl` for integration example.

## References

- **Dynamic_objectives**: See `README.md` for model details and testing guide
- **globtimcore**: See globtimcore documentation for optimization algorithm details
- **INTEGRATION_OBJECTIVES.md**: Original integration goals and success criteria
- **test_configs.jl**: Pre-configured test cases for all 13 models

## Contributing

When adding new models or improving integration:

1. Add model to appropriate difficulty phase in `test_configs.jl`
2. Add test case to `test/test_globtim_integration.jl`
3. Run tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
4. Update this guide with configuration recommendations
5. Document any special requirements (timeout, bounds, etc.)

## License

See LICENSE file in repository root.

---

*Last updated: 2025-11-20*
*Dynamic_objectives version: 0.1.0*
*globtimcore integration: TDD implementation*
