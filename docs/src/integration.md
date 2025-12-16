# Globtim Integration

**Last Updated**: 2025-11-23 (Phase 3)

## Overview

This guide explains how to integrate Dynamic\_objectives with the globtim ecosystem for critical point finding and parameter recovery through a **2-stage pipeline**:

1. **Stage 1 (globtimcore)**: Find raw critical points using polynomial approximation + HomotopyContinuation
2. **Stage 2 (globtimpostprocessing)**: Refine critical points using local optimization on original ODE objective

**Key Benefit**: After Phase 2, no wrapper functions needed! Dynamic\_objectives' 1-argument functions work directly with globtimcore.

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                    Dynamic_objectives                        │
│  Define ODE models, create objective functions               │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│ Stage 1: globtimcore                                         │
│  Polynomial approximation + HomotopyContinuation             │
│  Output: critical_points_raw_deg_X.csv                       │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│ Stage 2: globtimpostprocessing                               │
│  Local BFGS refinement                                       │
│  Output: critical_points_refined_deg_X.csv                   │
└──────────────────────────────────────────────────────────────┘
```

### Standalone Integration

**Design Principle**: Dynamic\_objectives remains **standalone** with no formal package dependencies on globtim packages.

**Integration Method**: Environment activation + local dev versions

**Setup** (one-time):

**Option 1: Quick setup** (recommended):
```bash
cd /path/to/Dynamic_objectives
./setup_dev_packages.jl  # Automated setup script
```

**Option 2: Manual setup**:
```bash
cd /path/to/Dynamic_objectives
julia --project=. -e 'using Pkg; Pkg.develop(path="../globtimcore"); Pkg.develop(path="../globtimpostprocessing")'
```

After this setup, you can use all three packages together in your scripts.

**Verification**:
```bash
# Verify packages are available
julia --project=. -e 'using Globtim, GlobtimPostProcessing; println("Setup successful")'
```

## Phase 2 Migration (November 2025)

### What Changed

**Before Phase 2** (1-stage):
- `run_standard_experiment()` found AND refined critical points
- CSV files: `critical_points_deg_X.csv` (already refined)
- Return dict had `:critical_points` key

**After Phase 2** (2-stage):
- `run_standard_experiment()` finds only RAW critical points
- CSV files: `critical_points_raw_deg_X.csv` (not refined)
- Refinement moved to `refine_experiment_results()` in globtimpostprocessing
- Return dict has `:degree_results` key (array of DegreeResult)
- Refined CSV files: `critical_points_refined_deg_X.csv`

### Why This Matters

- **Old code will break**: Scripts expecting refined points from Stage 1 will get raw points
- **New workflow required**: Must explicitly call Stage 2 for refinement
- **Better separation**: Core algorithms (Stage 1) independent of refinement (Stage 2)
- **More flexible**: Can skip refinement, use custom refinement, or refine later

## Quick Start

### Recommended Workflow

```bash
# Step 1: Verify model works (fast, uses rich display)
julia --project=. examples/verify_model.jl

# Step 2: Quick integration test (GN=10, fast)
julia --project=. examples/test_simple_workflow.jl

# Step 3: Full integration test (GN=50, production quality)
julia --project=. examples/test_globtim_integration.jl
```

### Minimal Example

```julia
using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra

# Include ExperimentCLI for config
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

# Setup model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Create objective (1-arg function)
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

# Configure experiment
config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,
    degree_range = 4:6,
    max_time = 3600.0,
    basis = :chebyshev
)

# Stage 1: Find raw critical points
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = "test_results",
    true_params = p_true
)

# Stage 2: Refine
refinement_config = ode_refinement_config(max_time_per_point = 30.0, show_progress = false)
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)

# Verify
best_params = refined.refined_points[refined.best_refined_idx]
recovery_error = norm(best_params .- p_true) / norm(p_true)
println("Recovery error: $(round(100*recovery_error, digits=2))%")
```

## Integration Patterns

Choose based on your needs:

| Pattern | Use Case | Complexity | Control | Template File |
|---------|----------|------------|---------|---------------|
| **Simple** | Quick tests, prototypes | Low | Low | `templates/integration_template_simple.jl` |
| **Pipeline** | Production workflows | Medium | Medium | `templates/integration_template_pipeline.jl` |
| **Advanced** | HPC, custom experiments | High | High | `templates/integration_template_advanced.jl` |

### Pattern 2: Pipeline (Recommended)

**When to use:**
- Production workflows
- Want both raw and refined results
- Need configurable parameters

**See**: `examples/test_simple_workflow.jl`

**Pros:**
- Full polynomial approximation + HomotopyContinuation
- Returns both raw and refined results
- Configurable via ExperimentParams
- Reasonable code complexity

**Cons:**
- Requires manual ExperimentCLI include
- More boilerplate than Simple pattern

### Pattern 3: Advanced Manual

**When to use:**
- Maximum control over all parameters
- Running on HPC cluster
- Custom experiment configurations
- Need to inspect intermediate results

**See**: `examples/test_globtim_integration.jl`, `templates/integration_template_advanced.jl`

## Function Signature Compatibility

**After Phase 2**, globtimcore automatically detects function signatures.

### 1-Argument Functions (Dynamic\_objectives Pattern)

```julia
# Dynamic_objectives creates 1-arg functions
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing
)

# Works directly with globtimcore (no wrapper!)
result = run_standard_experiment(
    objective_function = objective,  # Detected as 1-arg, used as-is
    problem_params = nothing,        # Signal that no params needed
    domain_bounds = bounds,
    experiment_config = config
)
```

**Why it works**: globtimcore checks function signature and adapts accordingly. If `problem_params = nothing`, it uses the function directly.

## Configuration

### ExperimentParams (globtimcore)

```julia
config = ExperimentParams(
    domain_size = 1.5,              # Search domain size
    GN = 50,                        # Grid size (GN^2 points for 2D)
    degree_range = 4:8,             # Polynomial degrees to try
    max_time = 3600.0,              # Max time in seconds
    basis = :chebyshev,             # Polynomial basis
)
```

**Grid Size Guidelines**:

| GN | Points (2D) | Speed | Use Case |
|----|-------------|-------|----------|
| 10 | 100 | Fast (~seconds) | Quick test |
| 20 | 400 | Reasonable (~minutes) | Testing |
| 50 | 2500 | Slow (~10-30 min) | Production |
| 100 | 10000 | Very slow | HPC batch |

### RefinementConfig (globtimpostprocessing)

```julia
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,      # Timeout per point
    show_progress = false,          # Suppress progress output
)
```

## Output Files

After running both stages:

```
output_dir/
├── experiment_config.json                  # Experiment configuration
├── metadata.json                           # Experiment metadata
├── critical_points_raw_deg_4.csv           # Raw points from Stage 1
├── critical_points_raw_deg_5.csv
├── critical_points_raw_deg_6.csv
├── critical_points_refined_deg_4.csv       # Refined points from Stage 2
├── critical_points_refined_deg_5.csv
├── critical_points_refined_deg_6.csv
├── refinement_comparison_deg_4.csv         # Side-by-side comparison
├── refinement_comparison_deg_5.csv
├── refinement_comparison_deg_6.csv
└── refinement_summary.json                 # Overall statistics
```

## Common Pitfalls & Troubleshooting

### Anti-Pattern 1: Importing `StandardExperimentConfig`

```julia
# WRONG - doesn't exist after Phase 2
using Globtim: StandardExperimentConfig

config = StandardExperimentConfig(max_degree=8, grid_size=50)  # Error!
```

**Fix**: Use `ExperimentParams` from ExperimentCLI.

### Anti-Pattern 2: Forgetting `problem_params = nothing`

```julia
# WRONG - will fail with 1-arg functions
result = run_standard_experiment(
    objective_function = my_1arg_function,
    domain_bounds = bounds,
    experiment_config = config
)
```

**Fix**: Always specify `problem_params`:
```julia
result = run_standard_experiment(
    objective_function = my_1arg_function,
    problem_params = nothing,  # Signals 1-arg function
    domain_bounds = bounds,
    experiment_config = config
)
```

### Anti-Pattern 3: Expecting refined points from Stage 1

```julia
# WRONG - Stage 1 returns RAW points only
result = run_standard_experiment(...)
best_params = result[:critical_points][1]  # These are RAW, not refined!
```

**Fix**: Always refine in Stage 2:
```julia
# Stage 1: Raw points
result = run_standard_experiment(...)

# Stage 2: Refine
refined = refine_experiment_results(result[:output_dir], objective, config)
best_params = refined.refined_points[refined.best_refined_idx]  # Refined
```

### Slow Grid Evaluation

**Solutions**:
1. **Reduce GN**: Start with GN=10 for testing
2. **Reduce degree\_range**: Use 4:6 instead of 4:12
3. **Increase eval\_timeout**: Prevent timeouts on expensive ODE solves
4. **Run verify\_model.jl first**: Get performance estimates

### ODE Integration Failures

**Symptoms**:
```
Warning: Objective returned NaN at point [...]
Warning: ODE integration failed
```

**Solutions**:
1. **Check initial conditions**: Ensure IC is in valid state space region
2. **Check parameter bounds**: May be exploring unphysical parameter regions
3. **Increase eval\_timeout**: ODE solver may need more time for stiff systems
4. **Narrow domain\_size**: Reduce search space to reasonable parameter regions

### No Critical Points Found

**Solutions**:
1. **Check domain bounds**: May be excluding regions with critical points
2. **Increase GN**: More grid points = better polynomial approximation
3. **Expand degree\_range**: Try higher degrees
4. **Check objective function**: Ensure it's well-behaved (smooth, finite everywhere)

## Performance Tips

### Before Running Expensive Experiments

```bash
# Run verification script to estimate performance
julia --project=. examples/verify_model.jl
```

This script shows:
- Model summary and parameter values
- ODE response visualization
- Objective function tests
- Time estimates for different grid sizes

### Grid Size Selection

Use the estimates from `verify_model.jl`:
- If 1 eval takes 10ms, GN=50 (2500 points) ≈ 25 seconds
- If 1 eval takes 100ms, GN=50 ≈ 4 minutes
- If 1 eval takes 1s, GN=50 ≈ 42 minutes (consider reducing to GN=20)

## Additional Documentation

- [Architecture](architecture.md) - Standalone architecture principles
- [Testing Guide](testing_guide.md) - Comprehensive testing strategy
- [Model Catalog](model_catalog.md) - Quick reference for all 13 models
- [Getting Started](getting_started.md) - Development environment setup

## Questions?

See example scripts in `examples/` or templates in `examples/templates/` for working code.
