"""
Integration utilities for using Dynamic_objectives with globtim optimizer.

This module provides a helper function to create globtim-compatible objectives from
Dynamic_objectives parameter estimation models.

# Architecture

The key challenge is signature compatibility:
- Dynamic_objectives: `error_func(p::Vector{Float64}) -> Float64`
- globtimcore expects: `objective(point::Vector{Float64}, params) -> Float64`

The helper function bridges this gap by creating a closure that captures the model
configuration and adapts the function signature.

# Usage

This module is lightweight and only provides the objective adapter. To run full
globtim experiments, use the standalone scripts in `examples/globtim_integration/`:

```julia
# In Dynamic_objectives environment
using Dynamic_objectives

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

# Create globtim-compatible objective
objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30
)

# Use objective in globtimcore environment (see examples/globtim_integration/)
```

See `examples/globtim_integration/run_single_model.jl` for complete workflow.
"""

using LinearAlgebra

export create_globtim_objective

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
        time_interval, numpoints;
        distance = distance,
        aggregate = aggregate,
        eval_timeout = eval_timeout,
        return_inf_on_error = return_inf_on_error
    )

    # Adapt to globtimcore signature (two arguments)
    # The `params` argument is ignored since all configuration is captured in closure
    function globtim_objective(point::Vector{Float64}, params)
        return error_func(point)
    end

    return globtim_objective
end

# Note: For complete workflow including optimization, see standalone scripts:
# - examples/globtim_integration/run_single_model.jl
# - examples/globtim_integration/run_batch_models.jl
# These scripts use globtimcore's StandardExperiment infrastructure directly.
