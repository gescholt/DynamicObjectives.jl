# PostProcessing Requirements

This document specifies the requirements for integrating Dynamic\_objectives with the globtimpostprocessing package for local refinement of grid search results.

## Overview

After `run_globtim_optimization()` finds promising parameter regions via grid search, we need to refine these results using local optimization methods provided by globtimpostprocessing.

## Current State

**What we have:**
- Grid-based global optimization that finds approximate solutions
- Real-time progress monitoring via `display_optimization_progress()`
- Results stored in structured format (Dict with :best\_params, :best\_objective, etc.)
- Optional file output (JSON summary + CSV critical points)

**What we need:**
- Local refinement starting from grid search results
- Improved parameter accuracy (20% error → <1% error)
- Progress monitoring during refinement
- Seamless integration with existing workflow

---

## Requirements

### 1. Input Format

**From Dynamic\_objectives to globtimpostprocessing:**

```julia
# What we can provide:
refinement_input = Dict(
    # Starting point(s) for local optimization
    :initial_points => Vector{Vector{Float64}},  # Best grid point(s)

    # Objective function (same as used in grid search)
    :objective => Function,  # f(point::Vector{Float64}, params) -> Float64

    # Search bounds
    :bounds => Vector{Tuple{Float64, Float64}},

    # Optional: Grid search context
    :grid_objective_values => Vector{Float64},  # All grid evaluations
    :grid_points => Vector{Vector{Float64}},    # All grid points

    # Optional: True parameters (for validation)
    :p_true => Vector{Float64}
)
```

**Alternative: File-based input**

If globtimpostprocessing expects files, we can provide:
- `critical_points_deg_*.csv` (already generated)
- `results_summary.json` (already generated)
- Additional metadata as needed

### 2. Expected Interface

**Desired function signature:**

```julia
using GlobtimPostprocessing

# Option A: Direct function call
result = refine_solution(
    objective::Function,
    initial_point::Vector{Float64},
    bounds::Vector{Tuple{Float64, Float64}};
    method = :nelder_mead,  # or :bfgs, :lbfgs, etc.
    max_iterations = 1000,
    tolerance = 1e-8,
    callback = nothing,  # For progress monitoring
    verbose = true
)

# Option B: Multi-start refinement
result = refine_multistart(
    objective::Function,
    initial_points::Vector{Vector{Float64}},
    bounds::Vector{Tuple{Float64, Float64}};
    n_best = 5,  # Refine top 5 grid points
    method = :nelder_mead,
    parallel = false,
    callback = nothing
)

# Option C: File-based workflow
result = refine_from_files(
    results_dir::String,  # Directory with critical_points*.csv
    objective::Function;
    method = :nelder_mead,
    output_dir = "refined_results"
)
```

### 3. Output Format

**From globtimpostprocessing back to Dynamic\_objectives:**

```julia
# Expected return structure:
refinement_result = Dict(
    # Best refined solution
    :best_params => Vector{Float64},
    :best_objective => Float64,

    # Convergence info
    :converged => Bool,
    :iterations => Int,
    :function_evaluations => Int,
    :time_elapsed => Float64,

    # Improvement metrics
    :initial_objective => Float64,
    :final_objective => Float64,
    :improvement_factor => Float64,

    # Optional: Optimization trajectory
    :parameter_history => Vector{Vector{Float64}},  # For visualization
    :objective_history => Vector{Float64},

    # Multi-start specific (if applicable)
    :all_solutions => Vector{Dict},  # All local minima found
    :n_unique_solutions => Int
)
```

### 4. Progress Monitoring Integration

**Callback mechanism for real-time updates:**

```julia
# We can provide a callback that uses our display infrastructure:
function refinement_callback(iteration, params, objective_value)
    display_optimization_progress(iteration, params, objective_value)
    return false  # Return true to stop optimization
end

# Usage:
result = refine_solution(
    objective,
    initial_point,
    bounds;
    callback = refinement_callback
)
```

**Requirements for callback:**
- Called at each iteration (or every N iterations)
- Provides: iteration number, current parameters, current objective value
- Return `false` to continue, `true` to stop early

### 5. Method Requirements

**Supported optimization methods:**

Must support at least one of:
- **Derivative-free methods** (preferred for ODE-based objectives):
  - Nelder-Mead simplex
  - Powell's method
  - Coordinate descent

- **Gradient-based methods** (if numerical gradients available):
  - BFGS
  - L-BFGS
  - Gradient descent

**Constraints handling:**
- Box constraints (parameter bounds) - **REQUIRED**
- General nonlinear constraints - optional

### 6. Error Handling

**Expected behavior:**

```julia
# If refinement fails, should not crash but return status
result = refine_solution(...)

if result[:converged]
    # Use refined solution
    println("Refined to: ", result[:best_params])
else
    # Fall back to grid search result
    @warn "Refinement failed: $(result[:termination_reason])"
    # Still have grid search result to use
end
```

### 7. Integration Workflow

**Proposed end-to-end workflow:**

```julia
using Dynamic_objectives
using GlobtimPostprocessing

# Step 1: Grid search (global optimization)
grid_result = run_globtim_optimization(
    objective, bounds, p_true;
    GN = 8,
    show_progress = true,
    save_results = true,
    output_dir = "grid_results"
)

# Step 2: Local refinement (via globtimpostprocessing)
refined_result = refine_solution(
    objective,
    grid_result[:best_params],  # Start from best grid point
    bounds;
    method = :nelder_mead,
    callback = (iter, p, obj) -> display_optimization_progress(iter, p, obj)
)

# Step 3: Display final results
println("\nGrid search result:")
println("  Params: ", grid_result[:best_params])
println("  Error: ", grid_result[:recovery_error])

println("\nRefined result:")
println("  Params: ", refined_result[:best_params])
println("  Error: ", norm(refined_result[:best_params] - p_true) / norm(p_true))
println("  Improvement: $(refined_result[:improvement_factor])x better")
```

---

## Questions for globtimpostprocessing

To finalize integration, we need to know:

1. **Package structure:**
   - What is the exact package name? `GlobtimPostprocessing` or different?
   - Is it registered or do we need to add from URL?
   - What are the main exported functions?

2. **Function signatures:**
   - What is the exact signature for local refinement?
   - Does it accept callbacks for progress monitoring?
   - Does it return structured results or raw output?

3. **File-based vs programmatic:**
   - Does it expect to read from files (CSV/JSON)?
   - Can it work purely in-memory (pass objective function directly)?
   - What directory structure does it expect?

4. **Objective function format:**
   - Does it expect `f(x)` or `f(x, params)`?
   - Does it accept `AbstractVector` or require `Vector{Float64}`?
   - Any special return value conventions?

5. **Dependencies:**
   - What optimization backend does it use? (Optim.jl, NLopt.jl, etc.)
   - Any special requirements or limitations?

---

## Integration Points

### Option A: Wrapper Function (Recommended)

Create a wrapper in `globtim_integration.jl`:

```julia
"""
    refine_grid_result(grid_result, objective, bounds; kwargs...)

Refine grid search results using local optimization from globtimpostprocessing.
"""
function refine_grid_result(
    grid_result::Dict,
    objective::Function,
    bounds::Vector{<:Tuple};
    method = :nelder_mead,
    show_progress = true,
    kwargs...
)
    # Extract starting point from grid result
    initial_point = grid_result[:best_params]

    # Setup progress callback if requested
    callback = if show_progress
        (iter, p, obj) -> display_optimization_progress(iter, p, obj)
    else
        nothing
    end

    # Call globtimpostprocessing
    refined = GlobtimPostprocessing.refine_solution(
        objective, initial_point, bounds;
        method = method,
        callback = callback,
        kwargs...
    )

    return refined
end
```

### Option B: Extended run\_globtim\_optimization

Add refinement as optional step:

```julia
result = run_globtim_optimization(
    objective, bounds, p_true;
    GN = 8,
    enable_refinement = true,  # NEW
    refinement_method = :nelder_mead,  # NEW
    show_progress = true
)

# Result includes both grid and refined solutions
result[:grid_params]      # Grid search result
result[:refined_params]   # Refined result (if enable_refinement=true)
result[:best_params]      # Best overall (refined if available, else grid)
```

---

## Success Criteria

Integration is successful when:

1. Can refine grid search results with <1% error (currently ~20%)
2. Progress monitoring works during refinement
3. Graceful fallback if refinement fails
4. Minimal additional dependencies
5. Clear documentation and examples
6. Works with all existing models

---

## Next Steps

1. **Investigate globtimpostprocessing:**
   - Locate package source code
   - Read documentation
   - Identify main API functions

2. **Create integration wrapper:**
   - Implement `refine_grid_result()` function
   - Handle different output formats
   - Add progress monitoring

3. **Test integration:**
   - Test with 2D model (should get <1% error)
   - Test with 3D/4D models
   - Verify all models benefit from refinement

4. **Update documentation:**
   - Add refinement examples
   - Update test\_integration.jl to show refinement
   - Document best practices
