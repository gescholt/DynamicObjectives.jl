"""
Integration utilities for using Dynamic_objectives with globtim optimizer.

This module provides wrapper functions that adapt Dynamic_objectives parameter estimation
benchmarks to work with the globtimcore optimization framework.

# Main Functions

- `create_globtim_objective`: Create globtim-compatible objective from Dynamic_objectives model
- `run_globtim_optimization`: Complete workflow for optimizing a model with globtim

# Usage

```julia
using Dynamic_objectives
using Globtim

# Define model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Create globtim-compatible objective
objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30
)

# Run optimization
result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 8,
    degree_range = 4:8,
    output_dir = "results/lv2d"
)
```

# Architecture

The key challenge is signature compatibility:
- Dynamic_objectives: `error_func(p::Vector{Float64}) -> Float64`
- globtimcore expects: `objective(point::Vector{Float64}, params) -> Float64`

The wrapper functions bridge this gap by creating closures that capture the model
configuration and adapt the function signature.
"""

using Globtim
using LinearAlgebra

export create_globtim_objective, run_globtim_optimization

"""
    create_globtim_objective(model, outputs, ic, p_true, time_interval, numpoints; kwargs...)

Create a globtim-compatible objective function from a Dynamic_objectives model.

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
  suitable for use with `run_standard_experiment` from globtimcore

# Example
```julia
model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30;
    eval_timeout = 10.0  # Enable for oscillatory models
)

# Use with globtim
result = run_standard_experiment(
    objective_function = objective,
    problem_params = (;),
    domain_bounds = [(0.0, 3.0), (0.0, 2.0)],
    ...
)
```
"""
function create_globtim_objective(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints;
    distance = L2_norm,
    aggregate = first,
    eval_timeout = nothing,
    return_inf_on_error = true
)
    # Create Dynamic_objectives error function
    # This captures all problem configuration in the closure
    error_func = make_error_distance(
        model, outputs, ic, p_true,
        time_interval, numpoints,
        distance, aggregate, nothing;
        return_inf_on_error = return_inf_on_error,
        eval_timeout = eval_timeout
    )

    # Wrapper for globtim's two-argument interface
    # The second argument (params) is ignored - all info is in the closure
    function globtim_objective(point::Vector{Float64}, params)
        return error_func(point)
    end

    return globtim_objective
end


"""
    run_globtim_optimization(objective, bounds, p_true; kwargs...)

Complete workflow for optimizing a Dynamic_objectives model with globtim.

# Arguments
- `objective::Function`: globtim-compatible objective function `f(point, params) -> Float64`
- `bounds::Vector{Tuple{Float64, Float64}}`: Parameter bounds `[(lb₁, ub₁), ..., (lbₙ, ubₙ)]`
- `p_true::Vector{Float64}`: True parameters (for computing recovery error)

# Keyword Arguments
- `GN::Int = 8`: Grid points per dimension (total points = GN^dim)
- `degree_range = 4:8`: Range of polynomial degrees to test
- `output_dir::String = "globtim_results"`: Directory for output files
- `model_name::String = "model"`: Name for this model (used in metadata)
- `max_time::Int = 3600`: Maximum time for optimization (seconds)
- `basis::Symbol = :chebyshev`: Polynomial basis (:chebyshev or :legendre)

# Returns
- `result::Dict`: Dictionary containing:
  - `:success_rate`: Fraction of degrees that succeeded
  - `:degree_results`: Array of results for each degree
  - `:total_time`: Total optimization time
  - `:best_params`: Best parameters found across all degrees
  - `:best_objective`: Best objective value
  - `:recovery_error`: `||p_best - p_true||` if p_true provided

# Example
```julia
# Create objective (see create_globtim_objective)
objective = create_globtim_objective(...)

# Run optimization
result = run_globtim_optimization(
    objective,
    [(0.0, 3.0), (0.0, 2.0)],  # bounds
    [1.0, 0.5];                 # p_true
    GN = 8,
    degree_range = 4:8,
    output_dir = "results/lv2d",
    model_name = "LV_2D_v3"
)

# Check results
if result[:success_rate] > 0.0
    println("Best params: ", result[:best_params])
    println("Recovery error: ", result[:recovery_error])
end
```
"""
function run_globtim_optimization(
    objective::Function,
    bounds::Vector{Tuple{Float64, Float64}},
    p_true::Vector{Float64};
    GN::Int = 8,
    degree_range = 4:8,
    output_dir::String = "globtim_results",
    model_name::String = "model",
    max_time::Int = 3600,
    basis::Symbol = :chebyshev
)
    # Create output directory
    mkpath(output_dir)

    # Create experiment configuration
    experiment_config = ExperimentParams(
        GN = GN,
        degree_range = degree_range,
        domain_size = 0.5,  # Not used (domain_bounds takes precedence)
        max_time = max_time,
        basis = basis
    )

    # Metadata for tracking
    metadata = Dict{String, Any}(
        "model_name" => model_name,
        "true_params" => p_true,
        "parameter_dimension" => length(p_true),
        "bounds" => bounds,
        "GN" => GN,
        "degree_range" => string(degree_range),
        "basis" => string(basis)
    )

    # Run optimization using globtimcore StandardExperiment
    result = run_standard_experiment(
        objective_function = objective,
        problem_params = (;),  # Empty named tuple (all info in closure)
        domain_bounds = bounds,
        experiment_config = experiment_config,
        output_dir = output_dir,
        metadata = metadata,
        true_params = p_true  # Enables recovery_error calculation
    )

    # Extract best solution across all degrees
    best_objective = Inf
    best_params = nothing
    best_degree = nothing

    for deg_result in result[:degree_results]
        if deg_result.status == "success" && deg_result.best_objective < best_objective
            best_objective = deg_result.best_objective
            best_params = deg_result.best_estimate
            best_degree = deg_result.degree
        end
    end

    # Calculate recovery error if solution found
    recovery_error = if best_params !== nothing
        norm(best_params - p_true)
    else
        Inf
    end

    # Add summary fields to result
    result[:best_params] = best_params
    result[:best_objective] = best_objective
    result[:best_degree] = best_degree
    result[:recovery_error] = recovery_error

    return result
end
