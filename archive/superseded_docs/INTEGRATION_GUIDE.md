# Integration Guide: Dynamic_objectives + globtimcore + globtimpostprocessing

**Last Updated:** 2025-11-23 (Phase 2 Migration)

## Overview

This guide explains how to integrate Dynamic_objectives with the globtim ecosystem for critical point finding and parameter recovery.

### Architecture

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
│  Output: raw_critical_points_deg_X.csv                       │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│ Stage 2: globtimpostprocessing                               │
│  Local BFGS refinement                                       │
│  Output: refined_critical_points_deg_X.csv                   │
└──────────────────────────────────────────────────────────────┘
```

## Phase 2 Migration Changes (November 2025)

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

- **Old code will break** - scripts expecting refined points from Stage 1 will get raw points
- **New workflow required** - must explicitly call Stage 2 for refinement
- **Better separation** - core algorithms (Stage 1) independent of refinement (Stage 2)
- **More flexible** - can skip refinement, use custom refinement, or refine later

## Integration Patterns

Dynamic_objectives supports **three integration patterns**. Choose based on your needs:

| Pattern | Use Case | Complexity | Control | Files |
|---------|----------|------------|---------|-------|
| **Simple** | Quick tests, prototypes | ⭐ Low | ⭐ Low | `integration_template_simple.jl` |
| **Pipeline** | Production workflows | ⭐⭐ Medium | ⭐⭐ Medium | `integration_template_pipeline.jl` |
| **Advanced** | HPC, custom experiments | ⭐⭐⭐ High | ⭐⭐⭐ High | `integration_template_advanced.jl` |

### Pattern 1: Simple Wrapper

**When to use:**
- Quick testing or exploration
- Don't need fine control over parameters
- Okay with sensible defaults

**Example:**
```julia
using Dynamic_objectives

# Create objective function
objective = make_error_distance(model, outputs, ic, p_true, tspan, n_points, L2_norm, first, nothing)

# Run optimization (simplified wrapper)
result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 50,
    degree_range = 4:8
)

# Check results
println("Recovery error: $(round(100*result[:recovery_error], digits=2))%")
```

**Pros:**
- ✅ Minimal code
- ✅ No config structs needed
- ✅ Built-in progress monitoring

**Cons:**
- ❌ Uses simplified grid-based approach (not full HomotopyContinuation)
- ❌ Limited customization
- ❌ Not suitable for HPC

**Template:** `examples/templates/integration_template_simple.jl`

### Pattern 2: Pipeline Wrapper

**When to use:**
- Production workflows
- Want both raw and refined results
- Need some customization but not maximum control

**Example:**
```julia
using Dynamic_objectives
using Globtim: run_standard_experiment
using GlobtimPostProcessing: refine_experiment_results, ode_refinement_config

# Include ExperimentCLI for config
include(joinpath(dirname(@__DIR__), "..", "globtimcore", "src", "ExperimentCLI.jl"))
using .ExperimentCLI

# Create config
config = ExperimentParams(
    domain_size = 1.5,
    GN = 50,
    degree_range = 4:8,
    basis = :chebyshev
)

# Stage 1: Find raw critical points
raw_result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    true_params = p_true
)

# Stage 2: Refine
refinement_config = ode_refinement_config(max_time_per_point = 30.0)
refined_result = refine_experiment_results(
    raw_result[:output_dir],
    objective,
    refinement_config
)
```

**Pros:**
- ✅ Full polynomial approximation + HomotopyContinuation
- ✅ Returns both raw and refined results
- ✅ Configurable via ExperimentParams
- ✅ Reasonable code complexity

**Cons:**
- ❌ Requires manual ExperimentCLI include
- ❌ More boilerplate than Simple pattern

**Template:** `examples/templates/integration_template_pipeline.jl`

### Pattern 3: Advanced Manual

**When to use:**
- Maximum control over all parameters
- Running on HPC cluster
- Custom experiment configurations
- Need to inspect intermediate results

**Example:**
```julia
using Dynamic_objectives
using Globtim: run_standard_experiment
using GlobtimPostProcessing: refine_experiment_results, ode_refinement_config

# Include ExperimentCLI
include(joinpath(dirname(@__DIR__), "..", "globtimcore", "src", "ExperimentCLI.jl"))
using .ExperimentCLI

# Create detailed config
config = ExperimentParams(
    domain_size = 1.5,
    GN = 50,
    degree_range = 4:12,
    max_time = 3600.0,
    basis = :chebyshev,
    optim_f_tol = 1e-6,
    optim_x_tol = 1e-6,
    max_iterations = 300,
    enable_gradient_computation = false,
    enable_hessian_computation = false,
    enable_bfgs_refinement = false
)

# Stage 1: Full control
raw_result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict("experiment_type" => "custom", "notes" => "Special config"),
    true_params = p_true
)

# Inspect per-degree results
for deg_result in raw_result[:degree_results]
    println("Degree $(deg_result.degree): $(deg_result.n_critical_points) points")
    println("  Best: $(deg_result.best_objective)")
    println("  Status: $(deg_result.status)")
end

# Stage 2: Custom refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 60.0,
    f_tol = 1e-12,
    x_tol = 1e-12,
    max_iterations = 2000,
    verbose = true
)

refined_result = refine_experiment_results(
    raw_result[:output_dir],
    objective,
    refinement_config
)
```

**Pros:**
- ✅ Complete control over all parameters
- ✅ Access to all result details
- ✅ Suitable for HPC batch jobs
- ✅ Can inspect intermediate results

**Cons:**
- ❌ More boilerplate code
- ❌ Requires understanding of all config options

**Template:** `examples/templates/integration_template_advanced.jl`

## Common Pitfalls & Anti-Patterns

### ❌ Anti-Pattern 1: Importing `StandardExperimentConfig`

```julia
# WRONG - doesn't exist
using Globtim: StandardExperimentConfig

config = StandardExperimentConfig(max_degree=8, grid_size=50)  # Error!
```

**Fix:**
```julia
# CORRECT - use ExperimentParams
include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
using .ExperimentCLI

config = ExperimentParams(
    degree_range = 4:8,
    GN = 50
)
```

### ❌ Anti-Pattern 2: Using Old Result Dict Keys

```julia
# WRONG - Phase 1 schema
result = run_standard_experiment(...)
n_points = result[:n_critical_points]  # Key doesn't exist!
critical_points = result[:critical_points]  # Key doesn't exist!
```

**Fix:**
```julia
# CORRECT - Phase 2 schema
result = run_standard_experiment(...)
n_points = result[:total_critical_points]  # ✓
degree_results = result[:degree_results]  # Array of DegreeResult
critical_points = vcat([dr.critical_points for dr in degree_results]...)  # ✓
```

### ❌ Anti-Pattern 3: Expecting Refined Points from Stage 1

```julia
# WRONG - Stage 1 now returns RAW points only
result = run_standard_experiment(...)
# These are NOT refined!
best_params = result[:degree_results][1].best_estimate
```

**Fix:**
```julia
# CORRECT - Must refine explicitly
raw_result = run_standard_experiment(...)
refined_result = refine_experiment_results(raw_result[:output_dir], objective, config)
best_params = refined_result[:refined_points][refined_result[:best_refined_idx]]
```

### ❌ Anti-Pattern 4: Using Old CSV Filenames

```julia
# WRONG - Phase 1 filename
csv_path = joinpath(output_dir, "critical_points_deg_8.csv")
```

**Fix:**
```julia
# CORRECT - Phase 2 filenames
csv_raw = joinpath(output_dir, "critical_points_raw_deg_8.csv")      # Stage 1 output
csv_refined = joinpath(output_dir, "critical_points_refined_deg_8.csv")  # Stage 2 output
```

### ❌ Anti-Pattern 5: Wrong `run_standard_experiment` Signature

```julia
# WRONG - old positional args
result = run_standard_experiment(objective, bounds, config)
```

**Fix:**
```julia
# CORRECT - keyword arguments
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict(),
    true_params = p_true
)
```

## Configuration Reference

### ExperimentParams Fields

```julia
ExperimentParams(
    # Domain configuration
    domain_size::Float64 = 0.1          # Search domain size

    # Grid configuration
    GN::Int = 5                          # Grid points per dimension

    # Polynomial approximation
    degree_range::AbstractRange{Int} = 4:4   # Degrees to try
    basis::Symbol = :chebyshev           # :chebyshev or :legendre

    # Time limits
    max_time::Float64 = 45.0             # Max seconds per degree

    # Optimization tolerances (grid refinement)
    optim_f_tol::Float64 = 1e-6
    optim_x_tol::Float64 = 1e-6
    max_iterations::Int = 300

    # Feature toggles (Phase 2: should be false)
    enable_gradient_computation::Bool = false
    enable_hessian_computation::Bool = false
    enable_bfgs_refinement::Bool = false
)
```

### Refinement Configuration

```julia
# For ODE problems (typical Dynamic_objectives use case)
refinement_config = ode_refinement_config(
    max_time_per_point::Float64 = 30.0   # Max seconds per point
    f_tol::Float64 = 1e-6                # Function tolerance
    x_tol::Float64 = 1e-6                # Parameter tolerance
    max_iterations::Int = 1000           # Max BFGS iterations
    verbose::Bool = false                # Print detailed output
)

# For general problems
refinement_config = general_refinement_config(...)
```

## Result Structures

### Stage 1 Output (run_standard_experiment)

```julia
result = Dict(
    :experiment_id => "experiment_20251123_142030",
    :total_critical_points => 45,
    :degrees_processed => 5,
    :success_rate => 1.0,
    :degree_results => [DegreeResult, ...],  # Array of per-degree results
    :output_dir => "/path/to/results",
    :total_time => 123.45,
    :timestamp => "20251123_142030",
    :schema_version => "2.0.0",
    :params_dict => @dict(GN, degree_range, domain_size_param, max_time)
)

# Access per-degree results
for deg_result in result[:degree_results]
    deg_result.degree             # Int: polynomial degree
    deg_result.status             # String: "success" or "failed"
    deg_result.n_critical_points  # Int: number of critical points found
    deg_result.critical_points    # Vector{Vector{Float64}}: the points
    deg_result.objective_values   # Vector{Float64}: objective at each point
    deg_result.best_objective     # Float64: best objective value
    deg_result.best_estimate      # Vector{Float64}: best parameter estimate
    deg_result.total_computation_time  # Float64: seconds
    deg_result.output_dir         # String: where CSV was saved
end
```

### Stage 2 Output (refine_experiment_results)

```julia
refined_result = Dict(
    :n_raw => 45,                           # Number of raw points
    :n_converged => 38,                     # Number that converged
    :refined_points => Vector{Vector{Float64}},  # Refined parameter vectors
    :refined_values => Vector{Float64},     # Objective values at refined points
    :best_refined_value => 1.234e-8,        # Best objective value
    :best_refined_idx => 12,                # Index of best point
    :mean_improvement => 15.3,              # Average improvement factor
    :refinement_time => 45.67               # Total refinement time
)
```

## File Outputs

### Directory Structure

```
output_dir/
├── experiment_summary.json                   # Experiment metadata
├── critical_points_raw_deg_4.csv            # Raw points from degree 4
├── critical_points_raw_deg_5.csv
├── critical_points_raw_deg_6.csv
├── ...
├── critical_points_refined_deg_4.csv        # Refined points from degree 4
├── critical_points_refined_deg_5.csv
└── ...
```

### CSV Format

**Raw critical points** (`critical_points_raw_deg_X.csv`):
```csv
param_1,param_2,...,param_N,objective_value
1.0234,0.5123,...,0.2345,0.000123
...
```

**Refined critical points** (`critical_points_refined_deg_X.csv`):
```csv
param_1,param_2,...,param_N,objective_value,converged
1.0000,0.5000,...,0.2500,0.00000001,true
...
```

## Testing Integration

### Quick Test

```bash
cd Dynamic_objectives
julia --project=. examples/test_simple_workflow.jl
```

Expected output:
```
Testing simple 2-stage workflow
================================================================================

Stage 1: run_standard_experiment...
✓ Found critical points across 5 degrees
  Total critical points: 23
  Best raw value: 0.001234

Stage 2: refine_experiment_results...
✓ Refined 21/23 points
  Best refined value: 0.00000001
  Improvement: 123.4x

Recovery error: 0.12%
✅ SUCCESS!
```

### Verify Templates

```bash
# Test simple template
julia --project=. examples/templates/integration_template_simple.jl

# Test pipeline template
julia --project=. examples/templates/integration_template_pipeline.jl

# Test advanced template
julia --project=. examples/templates/integration_template_advanced.jl
```

## Troubleshooting

### Problem: `StandardExperimentConfig` not found

**Cause:** Trying to import a type that doesn't exist

**Solution:** Use `ExperimentParams` and include `ExperimentCLI.jl` manually

### Problem: `result[:n_critical_points]` key error

**Cause:** Using Phase 1 schema with Phase 2 code

**Solution:** Use `result[:total_critical_points]` instead

### Problem: Raw points are not refined

**Cause:** Expecting Stage 1 to refine points (Phase 1 behavior)

**Solution:** Explicitly call `refine_experiment_results()` for Stage 2

### Problem: CSV file not found

**Cause:** Looking for old filename `critical_points_deg_X.csv`

**Solution:** Use new filename `critical_points_raw_deg_X.csv` or `critical_points_refined_deg_X.csv`

### Problem: `run_globtim_pipeline()` error

**Cause:** Function uses old API (needs Phase 2 update)

**Solution:** Use manual 2-stage pattern (Pipeline or Advanced template)

## Migration Checklist

If you have existing code using the old API, follow these steps:

- [ ] Replace `StandardExperimentConfig` with `ExperimentParams`
- [ ] Include `ExperimentCLI.jl` manually
- [ ] Update `run_standard_experiment()` calls to use keyword arguments
- [ ] Update result dict key access (`:n_critical_points` → `:total_critical_points`)
- [ ] Update result dict key access (`:critical_points` → `:degree_results`)
- [ ] Add explicit Stage 2 refinement call
- [ ] Update CSV filename references
- [ ] Update any custom post-processing to handle new schema
- [ ] Test with `examples/test_simple_workflow.jl`

## Further Reading

- `globtimcore/.claude/CLAUDE.md` - Core package documentation
- `globtimpostprocessing/.claude/CLAUDE.md` - Post-processing documentation
- `Dynamic_objectives/REFINEMENT_INTEGRATION_GUIDE.md` - Additional refinement details
- `.claude/CLAUDE.md` - Overall architecture and package responsibilities

## Questions?

If you encounter issues not covered here:

1. Check the template files in `examples/templates/`
2. Review `examples/test_simple_workflow.jl` for working example
3. Consult package-specific documentation
4. Check GitLab issues for known problems
