# Integration Guide: Dynamic_objectives + globtimcore + globtimpostprocessing

**Last Updated**: 2025-11-23 (Phase 3)
**Status**: ✅ Current

## Overview

This guide explains how to integrate Dynamic_objectives with the globtim ecosystem for critical point finding and parameter recovery through a **2-stage pipeline**:

1. **Stage 1 (globtimcore)**: Find raw critical points using polynomial approximation + HomotopyContinuation
2. **Stage 2 (globtimpostprocessing)**: Refine critical points using local optimization on original ODE objective

**Key Benefit**: After Phase 2, no wrapper functions needed! Dynamic_objectives' 1-argument functions work directly with globtimcore.

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

**Design Principle**: Dynamic_objectives remains **standalone** with no formal package dependencies on globtim packages.

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
julia --project=. -e 'using Globtim, GlobtimPostProcessing; println("✓ Setup successful")'
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
best_params = refined[:refined_points][refined[:best_refined_idx]]
recovery_error = norm(best_params .- p_true) / norm(p_true)
println("Recovery error: $(round(100*recovery_error, digits=2))%")
```

## Integration Patterns

Choose based on your needs:

| Pattern | Use Case | Complexity | Control | Template File |
|---------|----------|------------|---------|---------------|
| **Simple** | Quick tests, prototypes | ⭐ Low | ⭐ Low | `templates/integration_template_simple.jl` |
| **Pipeline** | Production workflows | ⭐⭐ Medium | ⭐⭐ Medium | `templates/integration_template_pipeline.jl` |
| **Advanced** | HPC, custom experiments | ⭐⭐⭐ High | ⭐⭐⭐ High | `templates/integration_template_advanced.jl` |

### Pattern 1: Simple Wrapper (deprecated)

**Note**: The `run_globtim_optimization()` wrapper is simplified and uses grid-based approach (not full HomotopyContinuation). For production use, prefer Pattern 2 or 3.

### Pattern 2: Pipeline (Recommended)

**When to use:**
- Production workflows
- Want both raw and refined results
- Need configurable parameters

**See**: `examples/test_simple_workflow.jl`

**Pros:**
- ✅ Full polynomial approximation + HomotopyContinuation
- ✅ Returns both raw and refined results
- ✅ Configurable via ExperimentParams
- ✅ Reasonable code complexity

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

**Pros:**
- ✅ Complete control over all parameters
- ✅ Access to all result details
- ✅ Suitable for HPC batch jobs
- ✅ Can inspect intermediate results

**Cons:**
- More boilerplate code
- Requires understanding of all config options

## Function Signature Compatibility

**After Phase 2**, globtimcore automatically detects function signatures.

### 1-Argument Functions (Dynamic_objectives Pattern)

```julia
# Dynamic_objectives creates 1-arg functions
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing
)

# Works directly with globtimcore (no wrapper!)
result = run_standard_experiment(
    objective_function = objective,  # ✅ Detected as 1-arg, used as-is
    problem_params = nothing,        # ✅ Signal that no params needed
    domain_bounds = bounds,
    experiment_config = config
)
```

**Why it works**: globtimcore checks function signature and adapts accordingly. If `problem_params = nothing`, it uses the function directly.

### 2-Argument Functions (globtim Native)

```julia
# Native globtim format
function objective(point::Vector{Float64}, params)
    # Use params for problem configuration
    return compute_cost(point, params)
end

result = run_standard_experiment(
    objective_function = objective,
    problem_params = my_config,  # Passed to objective
    domain_bounds = bounds,
    experiment_config = config
)
```

## Configuration

### ExperimentParams (globtimcore)

```julia
config = ExperimentParams(
    domain_size = 1.5,              # Search domain size
    GN = 50,                        # Grid size (GN^2 points for 2D)
    degree_range = 4:8,             # Polynomial degrees to try
    max_time = 3600.0,              # Max time in seconds
    basis = :chebyshev,             # Polynomial basis
    # Optional advanced parameters:
    # optim_f_tol = 1e-6,
    # optim_x_tol = 1e-6,
    # max_iterations = 300,
)
```

**Grid Size Guidelines**:
- GN = 10: Quick test (100 points for 2D) - ⭐ Fast (~seconds)
- GN = 20: Medium (400 points for 2D) - ⭐⭐ Reasonable (~minutes)
- GN = 50: Production (2500 points for 2D) - ⭐⭐⭐ Slow but accurate (~10-30 min)

**Performance Tuning**:
```julia
# Development/debugging: Fast iteration
config = ExperimentParams(domain_size=1.5, GN=10, degree_range=4:6)

# Testing: Reasonable quality, manageable time
config = ExperimentParams(domain_size=1.5, GN=20, degree_range=4:8)

# Production: High quality, long runtime
config = ExperimentParams(domain_size=1.5, GN=50, degree_range=4:12)

# HPC batch jobs: Maximum quality
config = ExperimentParams(domain_size=1.5, GN=100, degree_range=4:16)
```

**Balancing Quality vs. Speed**:
- **Grid size (GN)**: More points = better initial approximation, longer Stage 1
- **Degree range**: Higher degrees = more critical points found, longer HC solving
- **Domain size**: Larger domain = more exploration, more critical points
- **Refinement timeout**: Longer timeout = more refined points converge

**Rule of thumb**: Start with GN=10 and degree_range=4:6. If recovery is good, you're done. If not, incrementally increase GN or expand degree_range.

### RefinementConfig (globtimpostprocessing)

```julia
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,      # Timeout per point
    show_progress = false,          # Suppress progress output
    # Optional:
    # f_tol = 1e-12,
    # x_tol = 1e-12,
    # max_iterations = 2000,
)
```

**Custom Refinement Configurations**:

```julia
# Fast refinement (quick tests)
quick_config = ode_refinement_config(
    max_time_per_point = 10.0,
    f_tol = 1e-8,
    max_iterations = 500
)

# Aggressive refinement (high-precision)
precise_config = ode_refinement_config(
    max_time_per_point = 120.0,
    f_tol = 1e-14,
    x_tol = 1e-14,
    max_iterations = 5000
)

# Parallel-friendly (no progress bar, suitable for batch jobs)
batch_config = ode_refinement_config(
    max_time_per_point = 60.0,
    show_progress = false,
    f_tol = 1e-12
)

# Stiff ODE-optimized (slower timeouts for expensive evaluations)
stiff_config = ode_refinement_config(
    max_time_per_point = 300.0,    # 5 minutes per point
    f_tol = 1e-10,
    max_iterations = 1000
)
```

**Choosing refinement settings**:
- **Quick tests**: Use fast refinement with relaxed tolerances
- **Production runs**: Use default or precise config
- **Expensive ODE models**: Increase `max_time_per_point` significantly
- **HPC batch jobs**: Disable `show_progress` to avoid log spam

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

### ❌ Anti-Pattern 1: Importing `StandardExperimentConfig`

```julia
# WRONG - doesn't exist after Phase 2
using Globtim: StandardExperimentConfig

config = StandardExperimentConfig(max_degree=8, grid_size=50)  # Error!
```

**Fix**: Use `ExperimentParams` from ExperimentCLI:
```julia
include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
using .ExperimentCLI

config = ExperimentParams(domain_size=1.5, GN=50, degree_range=4:8)
```

### ❌ Anti-Pattern 2: Forgetting `problem_params = nothing`

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
    problem_params = nothing,  # ✅ Signals 1-arg function
    domain_bounds = bounds,
    experiment_config = config
)
```

### ❌ Anti-Pattern 3: Expecting refined points from Stage 1

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
best_params = refined[:refined_points][refined[:best_refined_idx]]  # ✅ Refined
```

### Slow Grid Evaluation

**Problem**: Grid evaluation taking too long

**Solutions**:
1. **Reduce GN**: Start with GN=10 for testing
2. **Reduce degree_range**: Use 4:6 instead of 4:12
3. **Increase eval_timeout**: Prevent timeouts on expensive ODE solves
4. **Run verify_model.jl first**: Get performance estimates before expensive runs

### ODE Integration Failures

**Problem**: Many NaN or Inf values during grid evaluation

**Symptoms**:
```
Warning: Objective returned NaN at point [...]
Warning: ODE integration failed
```

**Solutions**:
1. **Check initial conditions**: Ensure IC is in valid state space region
2. **Check parameter bounds**: May be exploring unphysical parameter regions
3. **Increase eval_timeout**: ODE solver may need more time for stiff systems
4. **Use tighter ODE tolerances**: Adjust in `make_error_distance()` if available
5. **Narrow domain_size**: Reduce search space to more reasonable parameter regions

**Example fix**:
```julia
# Before: Too wide, includes unphysical regions
bounds = [(0.0, 10.0), (0.0, 10.0)]

# After: Narrower, focused on realistic parameters
bounds = [(0.1, 3.0), (0.1, 2.0)]  # Exclude zero and extreme values
```

### No Critical Points Found

**Problem**: HomotopyContinuation finds no critical points for some degrees

**Solutions**:
1. **Check domain bounds**: May be excluding regions with critical points
2. **Increase GN**: More grid points = better polynomial approximation
3. **Expand degree_range**: Try higher degrees
4. **Check objective function**: Ensure it's well-behaved (smooth, finite everywhere)

**Diagnostic**:
```bash
# Check CSV files - if empty, no critical points found
ls -lh output_dir/critical_points_raw_deg_*.csv
```

### Poor Parameter Recovery

**Problem**: Best refined point is far from p_true

**Diagnostic steps**:
1. **Check raw results**: Are raw critical points close to p_true?
   - If YES: Refinement problem → increase refinement timeout or adjust tolerances
   - If NO: Stage 1 problem → increase GN or expand degree_range

2. **Check convergence**: Did refinement converge?
   ```julia
   # Look at refinement comparison file
   df = CSV.read("output_dir/refinement_comparison_deg_6.csv", DataFrame)
   # Check "converged" column - how many refined points converged?
   ```

3. **Visualize**: Are there multiple local minima?

**Solutions**:
- **Increase GN**: Better initial polynomial approximation
- **Expand degree_range**: Higher degrees capture more complexity
- **Increase refinement timeout**: Give BFGS more time to converge
- **Check objective function**: May have many similar local minima (identifiability issue)

### MethodError or Type Issues

**Problem**: MethodError when calling globtimcore functions

**Common causes**:
1. **Outdated globtimcore**: Pull latest changes
   ```bash
   cd ../globtimcore && git pull
   julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.precompile()'
   ```

2. **Wrong function signature**: Ensure you're using 1-argument objective
   ```julia
   # WRONG: Returns (value, gradient)
   objective(p) = (compute_cost(p), compute_grad(p))

   # CORRECT: Returns scalar value only
   objective(p) = compute_cost(p)
   ```

3. **Missing problem_params**: Always specify `problem_params = nothing` for 1-arg functions

## Performance Tips

### Before Running Expensive Experiments

```bash
# Run verification script to estimate performance
julia --project=. examples/verify_model.jl
```

This script:
- Shows model summary and parameter values
- Visualizes ODE response with terminal plots
- Tests objective function on sample points
- Estimates time for different grid sizes
- Helps you choose appropriate GN value

### Grid Size Selection

Use the estimates from `verify_model.jl`:
- If 1 eval takes 10ms, GN=50 (2500 points) ≈ 25 seconds grid evaluation
- If 1 eval takes 100ms, GN=50 ≈ 4 minutes
- If 1 eval takes 1s, GN=50 ≈ 42 minutes (consider reducing to GN=20)

## Example Scripts

### Quick Verification
```bash
julia --project=. examples/verify_model.jl
```
Fast verification with rich terminal display showing:
- Model structure
- Parameter values
- ODE response plots
- Performance estimates

### Simple Workflow
```bash
julia --project=. examples/test_simple_workflow.jl
```
Complete 2-stage pipeline with:
- GN=10 (fast, 100 grid points)
- Degree range 4:6
- Rich progress display
- Parameter recovery verification

### Full Integration Test
```bash
julia --project=. examples/test_globtim_integration.jl
```
Production-quality test with:
- GN=50 (2500 grid points)
- Degree range 4:8
- Complete verification
- Timing analysis

## Next Steps

1. **Setup**: Run the one-time `Pkg.develop()` setup (see "Standalone Integration" above)
2. **Verify**: Run `verify_model.jl` to understand your model's performance
3. **Test**: Run `test_simple_workflow.jl` with small grid (GN=10)
4. **Iterate**: Adjust configuration based on results
5. **Production**: Run `test_globtim_integration.jl` with larger grid (GN=50)

## Additional Documentation

- **ARCHITECTURE.md**: Standalone architecture principles
- **TESTING_GUIDE.md**: Comprehensive testing strategy
- **MODEL_CATALOG.md**: Quick reference for all 13 models
- **SETUP.md**: Development environment setup

## Questions?

See example scripts in `examples/` or templates in `examples/templates/` for working code.

---

**Last Updated**: 2025-11-23 (Phase 3)
**Package Version**: Dynamic_objectives v0.1.0
