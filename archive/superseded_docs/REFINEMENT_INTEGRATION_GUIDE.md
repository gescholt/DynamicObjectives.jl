# Critical Point Refinement - Integration Guide for Dynamic_objectives

**Status**: 📘 **INTEGRATION GUIDE**

**Updated**: 2025-11-22

## Testing Setup - REQUIRED FIRST STEP

**Dynamic_objectives is intentionally standalone** (no package dependencies on globtimcore or globtimpostprocessing). To run integration tests or use the refinement pipeline, you must first set up local dev versions:

### For Running Integration Tests

```bash
cd /Users/ghscholt/GlobalOptim/Dynamic_objectives

# Activate Dynamic_objectives environment
julia --project=. -e '
using Pkg

# Dev the local versions of globtimcore and globtimpostprocessing
Pkg.develop(path="../globtimcore")
Pkg.develop(path="../globtimpostprocessing")

# Verify they are available
using Globtim
using GlobtimPostProcessing
println("✅ Integration packages loaded successfully!")
'
```

**After this one-time setup**, the integration tests will work:
```bash
./run_tests.sh
# or
julia --project=. -e 'using Pkg; Pkg.test()'
```

### For Using Refinement Pipeline in Scripts

Once dev versions are set up (above), you can use the simple pattern:

```julia
using Pkg
Pkg.activate(".")  # Dynamic_objectives environment
using Dynamic_objectives
using Globtim  # Now available via dev
using GlobtimPostProcessing  # Now available via dev

# Your integration code here...
```

**IMPORTANT**: The `Pkg.develop()` commands only need to be run ONCE. After that, Dynamic_objectives will always find the local dev versions.

---

## Overview

This guide shows how to integrate Dynamic_objectives with the new 2-stage critical point refinement pipeline:

1. **Stage 1 (globtimcore)**: Find raw critical points using polynomial approximation + HomotopyContinuation
2. **Stage 2 (globtimpostprocessing)**: Refine critical points using local optimization on original ODE objective

**Key Benefit**: No wrapper functions needed! Dynamic_objectives' 1-argument functions work directly with globtimcore after Phase 2.

## Architecture: Standalone Integration

Dynamic_objectives remains **independent** with no direct package dependencies:

```
┌─────────────────────────────────────────────────────┐
│           Dynamic_objectives                        │
│  (ODE models, parameter estimation problems)        │
└─────────────────────────────────────────────────────┘
         │                            │
         │ calls                      │ calls
         ▼                            ▼
    ┌──────────┐              ┌──────────────────┐
    │globtim   │              │globtimpost       │
    │  core    │─────────────▶│  processing      │
    └──────────┘   produces   └──────────────────┘
                   CSV files
```

**Design Pattern**: Environment activation + standalone scripts

**No Changes to Project.toml**: Keep Dynamic_objectives lightweight

## Integration Pattern

### Pattern 1: Full Pipeline (Recommended)

**Use Case**: ODE parameter estimation where you need refined critical points

**File**: `examples/globtim_integration/full_pipeline.jl`

```julia
#!/usr/bin/env julia

# Full 2-stage pipeline: globtimcore → globtimpostprocessing

using Pkg

# Step 1: Load Dynamic_objectives to get objective function
Pkg.activate(".")
using DynamicObjectives

# Create ODE-based objective (1-argument function)
objective = create_lotka_volterra_objective(
    observed_data,
    param_bounds
)

# Step 2: Activate globtimcore environment and find raw critical points
println("\n" * "="^80)
println("STAGE 1: Finding raw critical points (globtimcore)")
println("="^80)

Pkg.activate("/Users/ghscholt/GlobalOptim/globtimcore")
using Globtim

result_raw = run_standard_experiment(
    objective,  # 1-arg function works directly!
    param_bounds,
    StandardExperimentConfig(
        max_degree = 18,
        grid_size = 100
    )
)

println("\nRaw critical points saved to: ", result_raw.output_dir)
println("Found $(result_raw.n_critical_points) raw critical points")

# Step 3: Activate globtimpostprocessing and refine
println("\n" * "="^80)
println("STAGE 2: Refining critical points (globtimpostprocessing)")
println("="^80)

Pkg.activate("/Users/ghscholt/GlobalOptim/globtimpostprocessing")
using GlobtimPostProcessing

result_refined = refine_experiment_results(
    result_raw.output_dir,
    objective,  # Same function, no wrapper needed
    ode_refinement_config()  # 60s timeout, gradient-free
)

println("\n" * "="^80)
println("RESULTS SUMMARY")
println("="^80)
println("Raw critical points:      ", result_refined.n_raw)
println("Converged after refine:   ", result_refined.n_converged)
println("Mean improvement:         ", result_refined.mean_improvement)
println("Best raw value:           ", result_refined.best_raw_value)
println("Best refined value:       ", result_refined.best_refined_value)
println("\nBest parameter estimate:")
println(result_refined.refined_points[result_refined.best_refined_idx])
```

**Output Files**:
```
globtim_results/lv4d_exp_TIMESTAMP/
├── critical_points_raw_deg_18.csv          # From globtimcore
├── critical_points_refined_deg_18.csv      # From globtimpostprocessing
├── refinement_comparison_deg_18.csv        # Side-by-side comparison
└── refinement_summary.json                 # Statistics
```

---

### Pattern 2: Raw Points Only (No Refinement)

**Use Case**: Quick exploration or when refinement is too slow

```julia
using Pkg
Pkg.activate("/Users/ghscholt/GlobalOptim/globtimcore")
using Globtim

# Load Dynamic_objectives function
include("../src/models/lotka_volterra.jl")
objective = create_lotka_volterra_objective(data, bounds)

# Find raw points only
result = run_standard_experiment(
    objective,
    bounds,
    StandardExperimentConfig(max_degree=12)  # Lower degree = faster
)

# Use raw points directly
using CSV, DataFrames
csv_file = joinpath(result.output_dir, "critical_points_raw_deg_$(result.degree).csv")
critical_points = CSV.read(csv_file, DataFrame)

# Find best raw point
best_idx = argmin(critical_points.objective)
best_params = [critical_points[best_idx, Symbol("p$i")] for i in 1:length(bounds)]
println("Best raw parameter estimate: ", best_params)
```

---

### Pattern 3: Refinement Only (Pre-existing Raw Points)

**Use Case**: You already have raw critical points from a previous run

```julia
using Pkg
Pkg.activate("/Users/ghscholt/GlobalOptim/globtimpostprocessing")
using GlobtimPostProcessing

# Load Dynamic_objectives function
include("../src/models/lotka_volterra.jl")
objective = create_lotka_volterra_objective(data, bounds)

# Refine existing experiment
result = refine_experiment_results(
    "/path/to/globtim_results/lv4d_exp_20251122_143022",
    objective,
    ode_refinement_config()
)

# Access best refined point
best_params = result.refined_points[result.best_refined_idx]
```

---

## Function Signature Compatibility

**After Phase 2**, globtimcore automatically detects function signatures:

### 1-Argument Functions (Dynamic_objectives Pattern)

```julia
# Dynamic_objectives creates 1-arg functions
function my_ode_objective(params::Vector{Float64})
    # Solve ODE with these parameters
    solution = solve_ode(params)
    # Compute cost
    return compute_cost(solution, observed_data)
end

# Works directly with globtimcore (no wrapper!)
result = run_standard_experiment(
    my_ode_objective,  # ✅ Detected as 1-arg, used as-is
    bounds,
    config
)
```

### 2-Argument Functions (Legacy Pattern)

```julia
# Old-style 2-arg function
function legacy_objective(params::Vector{Float64}, problem_data)
    return compute_cost(params, problem_data)
end

# Still works with problem_params
result = run_standard_experiment(
    legacy_objective,
    bounds,
    config;
    problem_params = my_data  # ✅ Detected as 2-arg, wrapped automatically
)
```

**Key Point**: You DON'T need to write wrapper functions anymore in `src/globtim_integration.jl`!

---

## Updating globtim_integration.jl (Optional)

**Current Code** (`src/globtim_integration.jl`):

The current adapter with `create_globtim_objective()` that wraps 1-arg → 2-arg is NO LONGER NEEDED after Phase 2.

**Recommended Update**:

```julia
"""
    globtim_integration.jl

Integration helpers for using Dynamic_objectives with globtimcore.

After Phase 2: 1-arg functions work directly, no wrapper needed!
"""

# BEFORE PHASE 2 (Old approach):
# function create_globtim_objective(dynam_obj_func, data)
#     return (params, _) -> dynam_obj_func(params)  # ❌ Unnecessary wrapper
# end

# AFTER PHASE 2 (New approach):
# Just pass your function directly to globtimcore!
# No helper needed - see examples/globtim_integration/full_pipeline.jl

"""
    run_globtim_pipeline(objective_func, bounds; max_degree=18, refine=true)

Convenience function to run the full 2-stage pipeline.

# Arguments
- `objective_func`: 1-argument function `f(p::Vector{Float64}) -> Float64`
- `bounds`: Vector of (min, max) tuples for each parameter
- `max_degree`: Polynomial degree for approximation (default: 18)
- `refine`: Whether to refine critical points (default: true)

# Returns
Tuple `(raw_result, refined_result)` or just `raw_result` if `refine=false`

# Example
```julia
using DynamicObjectives

objective = create_lotka_volterra_objective(data, bounds)
raw_result, refined_result = run_globtim_pipeline(
    objective,
    bounds;
    max_degree = 18,
    refine = true
)
```
"""
function run_globtim_pipeline(
    objective_func,
    bounds::Vector{Tuple{Float64, Float64}};
    max_degree::Int = 18,
    grid_size::Int = 100,
    refine::Bool = true
)
    using Pkg

    # Stage 1: globtimcore
    Pkg.activate("/Users/ghscholt/GlobalOptim/globtimcore")
    @eval using Globtim

    println("Stage 1: Finding raw critical points...")
    raw_result = @eval run_standard_experiment(
        $objective_func,
        $bounds,
        StandardExperimentConfig(
            max_degree = $max_degree,
            grid_size = $grid_size
        )
    )

    if !refine
        return raw_result
    end

    # Stage 2: globtimpostprocessing
    Pkg.activate("/Users/ghscholt/GlobalOptim/globtimpostprocessing")
    @eval using GlobtimPostProcessing

    println("Stage 2: Refining critical points...")
    refined_result = @eval refine_experiment_results(
        $(raw_result.output_dir),
        $objective_func,
        ode_refinement_config()
    )

    return (raw_result, refined_result)
end

# Export for convenience
export run_globtim_pipeline
```

---

## Example: Lotka-Volterra 4D Parameter Estimation

**File**: `examples/globtim_integration/lotka_volterra_4d.jl`

```julia
#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__ * "/../..")  # Dynamic_objectives

using DynamicObjectives
using DifferentialEquations
using CSV, DataFrames

# Load observed data
data = CSV.read("data/lotka_volterra_observed.csv", DataFrame)

# Define parameter bounds
bounds = [
    (0.5, 2.0),   # α (prey growth)
    (0.5, 2.0),   # β (predation rate)
    (0.5, 2.0),   # γ (predator death)
    (0.5, 2.0)    # δ (predator efficiency)
]

# Create objective function
function lotka_volterra_objective(params::Vector{Float64})
    α, β, γ, δ = params

    # ODE system
    function lv_system!(du, u, p, t)
        prey, predator = u
        du[1] = α * prey - β * prey * predator      # dprey/dt
        du[2] = δ * prey * predator - γ * predator  # dpredator/dt
    end

    # Solve ODE
    u0 = [1.0, 1.0]  # Initial populations
    tspan = (0.0, 10.0)
    prob = ODEProblem(lv_system!, u0, tspan)

    try
        sol = solve(prob, Tsit5(); saveat=data.time)

        # Compute cost (sum of squared errors)
        cost = sum((sol[1, :] .- data.prey).^2) +
               sum((sol[2, :] .- data.predator).^2)

        return cost
    catch e
        # If ODE solve fails (stiff/unstable), return large cost
        return Inf
    end
end

# Run full pipeline
println("Starting Lotka-Volterra 4D parameter estimation...")
println("Bounds: ", bounds)

# Option 1: Manual pipeline (explicit control)
Pkg.activate("/Users/ghscholt/GlobalOptim/globtimcore")
using Globtim

raw_result = run_standard_experiment(
    lotka_volterra_objective,
    bounds,
    StandardExperimentConfig(max_degree=18)
)

Pkg.activate("/Users/ghscholt/GlobalOptim/globtimpostprocessing")
using GlobtimPostProcessing

refined_result = refine_experiment_results(
    raw_result.output_dir,
    lotka_volterra_objective,
    ode_refinement_config()
)

# Display results
println("\n" * "="^80)
println("RESULTS")
println("="^80)
println("Best refined parameters:")
best_params = refined_result.refined_points[refined_result.best_refined_idx]
println("  α = ", best_params[1])
println("  β = ", best_params[2])
println("  γ = ", best_params[3])
println("  δ = ", best_params[4])
println("\nObjective value: ", refined_result.best_refined_value)
println("\nResults saved to: ", raw_result.output_dir)
```

---

## Troubleshooting

### Issue 1: "No method matching objective_function(::Vector{Float64})"

**Cause**: Function signature mismatch

**Solution**: Ensure your function accepts exactly 1 argument:
```julia
# ✅ CORRECT
function my_objective(params::Vector{Float64})
    return sum(params.^2)
end

# ❌ WRONG
function my_objective(params)  # Missing type annotation
    return sum(params.^2)
end
```

---

### Issue 2: ODE Solver Hangs During Refinement

**Cause**: Stiff ODE or unstable parameters

**Solution**: Use `ode_refinement_config()` with timeout:
```julia
refined = refine_experiment_results(
    result_raw.output_dir,
    objective,
    ode_refinement_config(max_time_per_point=60.0)  # 60 second timeout
)
```

The timeout prevents hanging and marks timed-out points as failed.

---

### Issue 3: "Package globtimcore not found"

**Cause**: Incorrect activation path

**Solution**: Use absolute path:
```julia
Pkg.activate("/Users/ghscholt/GlobalOptim/globtimcore")
```

Or check your current directory:
```julia
pwd()  # Should be /Users/ghscholt/GlobalOptim/Dynamic_objectives
```

---

### Issue 4: Refinement Returns `Inf` for All Points

**Cause**: ODE solver failing for all parameter combinations

**Solutions**:
1. Check parameter bounds (are they physically meaningful?)
2. Check ODE initial conditions
3. Try different ODE solver algorithm
4. Add error handling in objective function:

```julia
function robust_ode_objective(params::Vector{Float64})
    try
        sol = solve(prob, Tsit5())
        return compute_cost(sol)
    catch e
        @warn "ODE solve failed for params=$params: $e"
        return Inf
    end
end
```

---

## Performance Tips

### 1. Start with Lower Degree

```julia
# Fast exploration (minutes)
config_fast = StandardExperimentConfig(max_degree=12, grid_size=50)

# High accuracy (hours)
config_accurate = StandardExperimentConfig(max_degree=18, grid_size=100)
```

### 2. Parallelize Refinement (Future)

Currently sequential, but `refine_critical_points_batch()` supports:
```julia
config = RefinementConfig(parallel=true)  # Not yet implemented
```

### 3. Adjust Refinement Timeout

For fast objectives:
```julia
config = ode_refinement_config(max_time_per_point=10.0)  # 10 seconds
```

For expensive ODE solves:
```julia
config = ode_refinement_config(max_time_per_point=120.0)  # 2 minutes
```

---

## Files in This Repository

**No changes needed** to Dynamic_objectives for Phase 2 integration!

**Optional updates**:
- `src/globtim_integration.jl` - Simplify or remove wrapper functions
- `examples/globtim_integration/` - Add new example scripts using 2-stage pipeline

**Recommended new examples**:
- `examples/globtim_integration/full_pipeline.jl` - Complete 2-stage example
- `examples/globtim_integration/lotka_volterra_4d.jl` - LV parameter estimation
- `examples/globtim_integration/fitzhugh_nagumo.jl` - FHN parameter estimation

---

## Next Steps

### After Phase 2 (globtimcore) Completes:

1. **Test integration**:
   ```bash
   cd /Users/ghscholt/GlobalOptim/Dynamic_objectives
   julia examples/globtim_integration/full_pipeline.jl
   ```

2. **Create new examples** using 2-stage pipeline

3. **Update documentation** to remove references to wrapper functions

4. **Optional**: Add `run_globtim_pipeline()` convenience function to `src/globtim_integration.jl`

---

## Summary

**Key Benefits**:
- ✅ No wrapper functions needed (1-arg functions work directly)
- ✅ Explicit 2-stage pipeline (clear separation of concerns)
- ✅ Timeout support for stiff ODE problems
- ✅ Standalone integration (no package dependencies)
- ✅ Full control over each stage

**Integration Pattern**:
```julia
# Stage 1: Find raw critical points
Pkg.activate("globtimcore")
raw = run_standard_experiment(objective, bounds, config)

# Stage 2: Refine critical points
Pkg.activate("globtimpostprocessing")
refined = refine_experiment_results(raw.output_dir, objective, ode_config)
```

**Contact**:
- See `docs/API_DESIGN_REFINEMENT.md` for full specification
- See `docs/REFINEMENT_MIGRATION_COORDINATION.md` for coordination across repos
- See `globtimpostprocessing/REFINEMENT_PHASE1_STATUS.md` for refinement API details
