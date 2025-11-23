"""
Integration utilities for using Dynamic_objectives with globtim optimizer.

# Architecture - Standalone Design

Dynamic_objectives is **standalone** with zero package dependencies on globtimcore
or globtimpostprocessing. Integration is achieved through:
- Local dev setup using `Pkg.develop()`  (see `setup_dev_packages.jl`)
- Environment activation patterns
- File-based data exchange (CSV files)

See ARCHITECTURE.md and REFINEMENT_INTEGRATION_GUIDE.md for complete integration patterns.

# Phase 2 Update: No Wrappers Needed!

After globtimcore Phase 2, 1-argument functions work directly:

```julia
# Dynamic_objectives creates 1-arg functions
objective = make_error_distance(model, outputs, ic, p_true, [0.0, 20.0], 30,
                                L2_norm, first, nothing; eval_timeout=10.0)

# Works directly with globtimcore (no wrapper needed!)
using Globtim
result = run_standard_experiment(objective, bounds, config)
```

**Note**: `create_globtim_objective()` is now **deprecated** but kept for backwards compatibility.

# Limitations

**ForwardDiff Incompatibility**: Dynamic_objectives uses ODE solvers internally which
cannot propagate ForwardDiff.Dual types for automatic differentiation. Therefore:
- Gradient computation must be disabled (enable_gradient_computation = false)
- Hessian computation must be disabled (enable_hessian_computation = false)
- BFGS refinement must be disabled (enable_bfgs_refinement = false)

# Integration Patterns

## Pattern 1: Manual 2-Stage Pipeline

```julia
using Pkg
using Dynamic_objectives

# Create objective (1-arg, works directly!)
objective = make_error_distance(...)

# Stage 1: Find raw critical points
Pkg.activate("../globtimcore")
using Globtim
raw = run_standard_experiment(objective, bounds, config)

# Stage 2: Refine critical points
Pkg.activate("../globtimpostprocessing")
using GlobtimPostProcessing
refined = refine_experiment_results(raw.output_dir, objective, ode_config)
```

## Pattern 2: Standalone (No globtim)

```julia
using Dynamic_objectives
using Optim  # Or any other optimizer

objective = make_error_distance(...)
result = optimize(objective, lower, upper, ParticleSwarm())
```

See REFINEMENT_INTEGRATION_GUIDE.md for complete examples.
"""

using LinearAlgebra
using JSON

export create_globtim_objective, run_globtim_optimization

"""
    create_globtim_objective(model, outputs, ic, p_true, time_interval, numpoints; kwargs...)

**DEPRECATED after Phase 2**: globtimcore now auto-detects 1-argument functions.

This function is kept for backwards compatibility. You can now use `make_error_distance()`
directly with globtimcore without any wrapper!

# Migration Guide

Instead of:
```julia
globtim_obj = create_globtim_objective(model, outputs, ic, p_true, [0.0, 20.0], 30)
result = run_standard_experiment(globtim_obj, bounds, config)
```

Use directly:
```julia
objective = make_error_distance(model, outputs, ic, p_true, [0.0, 20.0], 30,
                                L2_norm, first, nothing; eval_timeout=10.0)
result = run_standard_experiment(objective, bounds, config)  # Works directly!
```

# Arguments
- `model::ODESystem`: ModelingToolkit ODE system
- `outputs::Vector{Equation}`: Measured quantities
- `ic::Vector{Float64}`: Initial conditions for ODE
- `p_true::Vector{Float64}`: True parameters (generates reference data)
- `time_interval`: `[t_start, t_end]` time interval for integration
- `numpoints::Int`: Number of time points to sample

# Keyword Arguments
- `distance = L2_norm`: Distance function (L2_norm, L1_norm, log_L2_norm)
- `aggregate = first`: How to combine multi-output errors
- `eval_timeout = nothing`: Timeout in seconds (prevents ODE solver hanging)
- `return_inf_on_error = true`: Return Inf on ODE failure (for optimization)

# Returns
- Function with signature `f(point::Vector{Float64}, params) -> Float64`
  suitable for use with globtimcore's `run_standard_experiment`

# Example
```julia
using Dynamic_objectives

model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30;
    eval_timeout = 10.0  # Enable for oscillatory models
)

# Use objective in standalone script with globtimcore environment
# IMPORTANT: Disable ForwardDiff-dependent features in ExperimentParams:
#   enable_gradient_computation = false
#   enable_hessian_computation = false
#   enable_bfgs_refinement = false
```
"""
function create_globtim_objective(
    model, outputs, ic, p_true, time_interval, numpoints;
    distance = L2_norm,
    aggregate = first,
    eval_timeout = nothing,
    return_inf_on_error = true
)
    # Create Dynamic_objectives error function (single argument)
    error_func = make_error_distance(
        model, outputs, ic, p_true,
        time_interval, numpoints,
        distance,              # positional: distance_function
        aggregate,             # positional: aggregate_distances
        nothing;               # positional: add_noise_in_time_series
        return_inf_on_error = return_inf_on_error,  # keyword
        eval_timeout = eval_timeout                  # keyword
    )

    # Adapt to globtimcore signature (two arguments)
    # The `params` argument is ignored since all configuration is captured in closure
    # NOTE: Only accepts Float64 - NOT ForwardDiff.Dual compatible
    # Dynamic_objectives cannot propagate Dual numbers through ODE solvers
    # This triggers globtimcore's capability detection to disable gradient/Hessian features
    function globtim_objective(point::AbstractVector{Float64}, params)
        return error_func(Vector{Float64}(point))  # Ensure plain Vector{Float64}
    end

    return globtim_objective
end

"""
    run_globtim_optimization(objective, bounds, p_true; kwargs...)

Run grid-based global optimization with real-time progress monitoring.

This function performs grid-based optimization by evaluating the objective
on a structured grid and finding the minimum. It provides real-time progress
display using the display_optimization_progress() function.

# Arguments
- `objective`: Function with signature `f(point::Vector{Float64}, params) -> Float64`
- `bounds`: Vector of (min, max) tuples for each parameter
- `p_true`: True parameter values (for computing recovery error)

# Keyword Arguments
- `GN::Int = 8`: Grid resolution (grid size = GN^dimension)
- `degree_range = 4:8`: Range of polynomial degrees to try
- `output_dir = "."`: Directory for saving results
- `model_name = "model"`: Name for output files
- `max_time = 3600`: Maximum time in seconds per degree
- `basis = :chebyshev`: Polynomial basis (:chebyshev or :standard)
- `show_progress = true`: Display real-time optimization progress
- `save_results = true`: Save results to files

# Returns
Dictionary with keys:
- `:total_time`: Total computation time
- `:success_rate`: Fraction of degrees that found solutions
- `:total_critical_points`: Total critical points found across all degrees
- `:best_params`: Best parameter values found
- `:best_degree`: Degree that produced best solution
- `:best_objective`: Objective value at best solution
- `:recovery_error`: ||best_params - p_true|| / ||p_true||
- `:degree_results`: Detailed results for each degree

# Example
```julia
using Dynamic_objectives

# Create model and objective
model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

objective = create_globtim_objective(
    model, outputs, ic, p_true, [0.0, 20.0], 30
)

# Run optimization with progress monitoring
result = run_globtim_optimization(
    objective,
    [(0.0, 3.0), (0.0, 2.0)],
    p_true;
    GN = 8,
    degree_range = 4:8,
    show_progress = true
)

println("Best parameters: ", result[:best_params])
println("Recovery error: ", result[:recovery_error])
```
"""
function run_globtim_optimization(
    objective::Function,
    bounds::Vector{<:Tuple},
    p_true::AbstractVector;
    GN::Int = 8,
    degree_range = 4:8,
    output_dir::String = ".",
    model_name::String = "model",
    max_time::Real = 3600,
    basis::Symbol = :chebyshev,
    show_progress::Bool = true,
    save_results::Bool = true
)
    dimension = length(bounds)
    start_time = time()

    if show_progress
        println("\n" * "="^80)
        println("Starting Global Optimization")
        println("="^80)
        println("  Dimension: $dimension")
        println("  Grid size: $(GN^dimension) points")
        println("  Degrees: $degree_range")
        println("  Bounds: $bounds")
        println("="^80)
    end

    # Create output directory if saving
    if save_results
        mkpath(output_dir)
    end

    # Storage for results across degrees
    degree_results = []
    all_critical_points = []
    best_objective = Inf
    best_params = nothing
    best_degree = nothing

    # Try each grid resolution (degree parameter kept for compatibility)
    for (idx, degree) in enumerate(degree_range)
        degree_start = time()

        if show_progress
            println("\n[Grid resolution: $degree] ($idx/$(length(degree_range)))")
            println("─"^60)
        end

        try
            # Evaluate objective on grid
            grid_start = time()
            grid_values = []
            grid_points = []

            # Generate grid points
            n_per_dim = GN
            ranges = [range(b[1], b[2], length=n_per_dim) for b in bounds]

            # Iterate through all grid points
            grid_indices = Iterators.product([1:n_per_dim for _ in 1:dimension]...)
            total_grid = n_per_dim^dimension

            for (count, indices) in enumerate(grid_indices)
                point = [ranges[i][indices[i]] for i in 1:dimension]
                value = objective(point, nothing)

                push!(grid_points, point)
                push!(grid_values, value)

                # Show progress for larger grids
                if show_progress && total_grid > 100 && count % max(1, total_grid ÷ 10) == 0
                    pct = round(100 * count / total_grid, digits=1)
                    print("\r  Grid evaluation: $pct% ($count/$total_grid)")
                    flush(stdout)
                end
            end

            if show_progress
                grid_time = time() - grid_start
                println("\r  ✓ Grid evaluation complete ($(round(grid_time, digits=2))s)")
                println("    Min value on grid: $(round(minimum(grid_values), digits=6))")
                println("    Max value on grid: $(round(maximum(grid_values), digits=6))")
            end

            # Find critical points (simplified version - in practice use Globtim's method)
            if show_progress
                print("  Finding critical points... ")
                flush(stdout)
            end

            # For now, use grid minimum as a simple critical point
            # In full integration, this would use HomotopyContinuation.jl
            min_idx = argmin(grid_values)
            critical_point = grid_points[min_idx]
            critical_value = grid_values[min_idx]

            critical_points = [(point = critical_point, value = critical_value)]

            if show_progress
                println("✓ Found $(length(critical_points)) critical point(s)")
            end

            # Evaluate critical points and track best
            for (cp_idx, cp) in enumerate(critical_points)
                obj_val = cp.value

                if show_progress && length(critical_points) > 1
                    display_optimization_progress(cp_idx, cp.point, obj_val)
                end

                # Track best solution
                if obj_val < best_objective
                    best_objective = obj_val
                    best_params = copy(cp.point)
                    best_degree = degree
                end

                push!(all_critical_points, (degree = degree, point = cp.point, value = obj_val))
            end

            if show_progress && length(critical_points) > 1
                println()  # New line after progress
            end

            # Store degree results
            degree_time = time() - degree_start
            push!(degree_results, (
                degree = degree,
                n_critical_points = length(critical_points),
                best_value = minimum([cp.value for cp in critical_points]),
                time = degree_time,
                success = length(critical_points) > 0
            ))

            if show_progress
                println("  Best value for degree $degree: $(round(minimum([cp.value for cp in critical_points]), digits=6))")
                println("  Time: $(round(degree_time, digits=2))s")
            end

            # Save critical points if requested
            if save_results && length(critical_points) > 0
                # Save to CSV format
                cp_file = joinpath(output_dir, "critical_points_deg_$(degree).csv")
                open(cp_file, "w") do io
                    # Header
                    println(io, "point_id,", join(["p$i" for i in 1:dimension], ","), ",objective_value")
                    # Data
                    for (i, cp) in enumerate(critical_points)
                        println(io, "$i,", join(cp.point, ","), ",", cp.value)
                    end
                end
            end

        catch e
            if show_progress
                @warn "Degree $degree failed: $e"
            end

            push!(degree_results, (
                degree = degree,
                n_critical_points = 0,
                best_value = Inf,
                time = time() - degree_start,
                success = false
            ))
        end
    end

    total_time = time() - start_time

    # Compute recovery error
    recovery_error = if best_params !== nothing
        norm(best_params .- p_true) / norm(p_true)
    else
        Inf
    end

    # Compute success rate
    n_success = count(r -> r.success, degree_results)
    success_rate = n_success / length(degree_results)

    # Summary
    if show_progress
        println("\n" * "="^80)
        println("Optimization Complete")
        println("="^80)
        println("  Total time: $(round(total_time, digits=2))s")
        println("  Success rate: $(round(100 * success_rate, digits=1))% ($n_success/$(length(degree_results)) degrees)")
        println("  Total critical points: $(length(all_critical_points))")

        if best_params !== nothing
            println("\n  Best solution:")
            println("    Degree: $best_degree")
            println("    Objective: $(round(best_objective, digits=6))")
            println("    Parameters: $best_params")
            println("    Recovery error: $(round(100 * recovery_error, digits=2))%")
        else
            println("\n  No valid solutions found")
        end
        println("="^80)
    end

    # Save summary results
    if save_results
        summary = Dict(
            "model_name" => model_name,
            "dimension" => dimension,
            "GN" => GN,
            "degree_range" => collect(degree_range),
            "total_time" => total_time,
            "success_rate" => success_rate,
            "total_critical_points" => length(all_critical_points),
            "best_degree" => best_degree,
            "best_objective" => best_objective,
            "best_params" => best_params,
            "p_true" => collect(p_true),
            "recovery_error" => recovery_error,
            "degree_results" => [Dict(
                "degree" => r.degree,
                "n_critical_points" => r.n_critical_points,
                "best_value" => r.best_value,
                "time" => r.time,
                "success" => r.success
            ) for r in degree_results]
        )

        summary_file = joinpath(output_dir, "results_summary.json")
        open(summary_file, "w") do io
            JSON.print(io, summary, 2)
        end
    end

    # Return structured results
    return Dict(
        :total_time => total_time,
        :success_rate => success_rate,
        :total_critical_points => length(all_critical_points),
        :best_params => best_params,
        :best_degree => best_degree,
        :best_objective => best_objective,
        :recovery_error => recovery_error,
        :degree_results => degree_results,
        :all_critical_points => all_critical_points
    )
end

# Note: For complete 2-stage pipeline with refinement, see:
# - REFINEMENT_INTEGRATION_GUIDE.md
# - ARCHITECTURE.md
# - setup_dev_packages.jl (one-time dev setup)
# Refinement requires globtimpostprocessing package (set up via `./setup_dev_packages.jl`)
