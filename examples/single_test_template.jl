#!/usr/bin/env julia
# Single Model Test Template
# Usage: julia --project=. examples/single_test_template.jl

using DynamicObjectives
using LinearAlgebra
using Printf

println("="^80)
println("DynamicObjectives - Single Model Test Template")
println("="^80)

# ==============================================================================
# CONFIGURATION - MODIFY THIS SECTION
# ==============================================================================

# Choose your model (see docs/src/model_catalog.md for all options)
model_name = "LV 2D v3 (2 outputs)"
model_fn = define_lotka_volterra_2D_model_v3_two_outputs

# True parameters (what the optimizer should find)
p_true = [1.0, 0.5]

# Initial conditions
ic = [1.0, 0.5]

# Parameter bounds for optimization
bounds = [
    (0.0, 3.0),  # Parameter 1: a
    (0.0, 2.0),   # Parameter 2: b
]

# Time series configuration
time_interval = [0.0, 20.0]
numpoints = 30

# Distance function: L2_norm, L1_norm, or log_L2_norm
distance_function = L2_norm

# Timeout per evaluation (seconds, nothing = no timeout)
eval_timeout = nothing  # Use 10.0 for hard problems

# ==============================================================================
# SETUP - NO NEED TO MODIFY
# ==============================================================================

println("\nSetting up model: $model_name")
println("  True parameters: $p_true")
println("  Bounds: $bounds")
println("  Time interval: $time_interval")
println("  Sample points: $numpoints")

# Define the model
model, params, states, outputs = model_fn()
println("  ✓ Model defined ($(length(params)) parameters, $(length(outputs)) outputs)")

# Create objective function
error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints,
    distance_function;
    return_inf_on_error = true,
    eval_timeout = eval_timeout,
)
println("  ✓ Objective function created")

# Verify setup
println("\nVerifying setup...")
error_at_true = error_func(p_true)
@printf "  Error at true params: %.6e\n" error_at_true

if error_at_true > 1e-6
    @warn "Error at true parameters is not near zero! Check your setup."
else
    println("  ✓ Setup verified (error ≈ 0 at true params)")
end

# Test at perturbed parameters
p_test = p_true .+ 0.1
error_at_test = error_func(p_test)
@printf "  Error at perturbed params: %.6e\n" error_at_test

if error_at_test <= 0
    @warn "Error at perturbed parameters is not > 0! Check your setup."
else
    println("  ✓ Objective function responds to parameter changes")
end

# ==============================================================================
# GLOBTIM INTEGRATION PLACEHOLDER
# ==============================================================================

println("\n" * "="^80)
println("Ready for optimization!")
println("="^80)

println("\nTo integrate with globtim, use:")
println("""
# Pseudo-code for globtim integration:
using Globtim  # Your optimizer

result = globtim_optimize(
    objective = error_func,
    bounds = bounds,
    # ... your globtim options ...
)

p_best = result.minimizer
error_best = result.minimum

println("Best parameters found: \$p_best")
println("Best error: \$error_best")
println("True parameters: \$p_true")
println("Parameter error: \$(norm(p_best - p_true))")
""")

# ==============================================================================
# EXAMPLE: Grid Search (for testing without globtim)
# ==============================================================================

println("\n" * "="^80)
println("Example: Simple Grid Search")
println("="^80)

function simple_grid_search(error_func, bounds, n_samples = 5)
    println("\nRunning grid search with $n_samples points per dimension...")

    ndim = length(bounds)
    grids = [range(b[1], b[2], length = n_samples) for b in bounds]

    best_p = p_true
    best_error = Inf
    n_evals = 0

    for p in Iterators.product(grids...)
        p_vec = collect(p)
        error = error_func(p_vec)
        n_evals += 1

        if error < best_error
            best_error = error
            best_p = p_vec
        end
    end

    return best_p, best_error, n_evals
end

# Run simple grid search
p_best, error_best, n_evals = simple_grid_search(error_func, bounds, 5)

println("\nGrid search results:")
@printf "  Evaluations: %d\n" n_evals
@printf "  Best error: %.6e\n" error_best
println("  Best params: $p_best")
println("  True params: $p_true")
@printf "  Parameter error: %.6e\n" norm(p_best - p_true)

if error_best < 1e-4
    println("\n✓ Grid search found good solution!")
else
    println("\n⚠ Grid search did not find optimal solution (expected for coarse grid)")
end

println("\n" * "="^80)
println("Next steps:")
println("  1. Verify the objective function works with your setup")
println("  2. Replace grid search with your globtim optimizer")
println("  3. See examples/test_configs.jl for pre-configured test cases")
println("="^80)
