#!/usr/bin/env julia
"""
Standalone globtim Integration Test - Single Model

Tests globtimcore optimization on a single Dynamic_objectives model.
This script uses globtimcore's infrastructure as intended.

Usage:
    julia run_single_model.jl

Environment:
    - Activates globtimcore project (uses its full environment)
    - Loads both Dynamic_objectives and Globtim packages
    - Uses globtimcore's StandardExperiment infrastructure

Example:
    cd Dynamic_objectives/examples/globtim_integration
    julia run_single_model.jl
"""

using Pkg

# Activate globtimcore environment (required for StandardExperiment)
GLOBTIMCORE_PATH = joinpath(@__DIR__, "..", "..", "..", "globtimcore")
println("Activating globtimcore environment: $GLOBTIMCORE_PATH")
Pkg.activate(GLOBTIMCORE_PATH)

# Load packages from both environments
using Glob

tim
using Dynamic_objectives

# Include globtimcore modules (standard pattern from their examples)
GLOBTIM_SRC = joinpath(dirname(pathof(Globtim)), "..", "src")

if !isdefined(Main, :ExperimentCLI)
    include(joinpath(GLOBTIM_SRC, "ExperimentCLI.jl"))
end
using .ExperimentCLI

if !isdefined(Main, :StandardExperiment)
    include(joinpath(GLOBTIM_SRC, "StandardExperiment.jl"))
end
using .StandardExperiment

println("="^80)
println("globtim Integration Test - Single Model (2D LV v3)")
println("="^80)

# ==============================================================================
# Model Definition
# ==============================================================================

println("\n[1/5] Defining model...")
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]
time_interval = [0.0, 20.0]
numpoints = 30

println("  Model: Lotka-Volterra 2D v3 (2 outputs)")
println("  True parameters: $p_true")
println("  Parameter bounds: $bounds")
println("  Time interval: $time_interval")
println("  Sample points: $numpoints")

# ==============================================================================
# Create globtim-Compatible Objective
# ==============================================================================

println("\n[2/5] Creating objective function...")
objective = create_globtim_objective(
    model, outputs, ic, p_true,
    time_interval, numpoints;
    distance = L2_norm,
    aggregate = first,
    eval_timeout = nothing,  # Not needed for Lotka-Volterra
    return_inf_on_error = true
)

# Test objective at true parameters (should be ≈0)
test_error = objective(p_true, (;))
println("  ✓ Objective function created")
@printf "  Error at true params: %.6e (should be ≈0)\n" test_error

# ==============================================================================
# Run globtim Optimization
# ==============================================================================

println("\n[3/5] Running globtim optimization...")
GN = 8
degree_range = 4:8
output_dir = joinpath(@__DIR__, "..", "..", "test_results", "globtim_single", "LV_2D_v3")

println("  Configuration:")
println("    Grid points (GN): $GN")
println("    Grid size: $(GN)^2 = $(GN^2) evaluations")
println("    Polynomial degrees: $degree_range")
println("    Basis: Chebyshev")

# Create experiment configuration
experiment_config = ExperimentParams(
    GN = GN,
    degree_range = degree_range,
    domain_size = 0.5,  # Not used (domain_bounds takes precedence)
    max_time = 3600,
    basis = :chebyshev
)

# Metadata
metadata = Dict{String, Any}(
    "model_name" => "LV_2D_v3_two_outputs",
    "true_params" => p_true,
    "parameter_dimension" => length(p_true),
    "bounds" => bounds,
    "GN" => GN,
    "degree_range" => string(degree_range),
    "basis" => "chebyshev"
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

println("  ✓ Optimization complete")

# ==============================================================================
# Extract and Validate Results
# ==============================================================================

println("\n[4/5] Analyzing results...")

# Extract best solution across all degrees
best_objective = Inf
best_params = nothing
best_degree = nothing
n_success = 0

for deg_result in result[:degree_results]
    if deg_result.status == "success"
        n_success += 1
        if deg_result.best_objective < best_objective
            best_objective = deg_result.best_objective
            best_params = deg_result.best_estimate
            best_degree = deg_result.degree
        end
    end
end

success_rate = n_success / length(result[:degree_results])
recovery_error = best_params !== nothing ? norm(best_params - p_true) : Inf

println("  Success rate: $(n_success)/$(length(result[:degree_results])) ($(round(success_rate * 100, digits=1))%)")
println("  Total time: $(round(result[:total_time], digits=2)) seconds")
println("  Total critical points: $(result[:total_critical_points])")

if best_params !== nothing
    println("\n  Best solution:")
    println("    Degree: $best_degree")
    println("    Parameters: $best_params")
    println("    True params: $p_true")
    println("    Objective: $(round(best_objective, sigdigits=6))")
    println("    Recovery error: $(round(recovery_error, sigdigits=6))")
end

# ==============================================================================
# Validation
# ==============================================================================

println("\n[5/5] Validation...")

# Success criteria for EASY models (2D)
EASY_SUCCESS_THRESHOLD = 0.75  # >75% success rate
EASY_RECOVERY_THRESHOLD = 0.01  # <0.01 recovery error

passed = true

if success_rate < EASY_SUCCESS_THRESHOLD
    println("  ✗ FAIL: Success rate $(round(success_rate * 100, digits=1))% < $(EASY_SUCCESS_THRESHOLD * 100)%")
    passed = false
else
    println("  ✓ PASS: Success rate $(round(success_rate * 100, digits=1))% ≥ $(EASY_SUCCESS_THRESHOLD * 100)%")
end

if best_params !== nothing
    if recovery_error < EASY_RECOVERY_THRESHOLD
        println("  ✓ PASS: Recovery error $(round(recovery_error, sigdigits=4)) < $EASY_RECOVERY_THRESHOLD")
    else
        println("  ✗ FAIL: Recovery error $(round(recovery_error, sigdigits=4)) ≥ $EASY_RECOVERY_THRESHOLD")
        passed = false
    end
else
    println("  ✗ FAIL: No solution found")
    passed = false
end

# ==============================================================================
# Summary
# ==============================================================================

println("\n" * "="^80)
if passed
    println("✅ INTEGRATION TEST PASSED")
    println("\nThe globtim integration successfully:")
    println("  - Created a compatible objective function")
    println("  - Ran optimization using globtimcore infrastructure")
    println("  - Recovered true parameters within tolerance")
else
    println("❌ INTEGRATION TEST FAILED")
    println("\nSee results above for details.")
    println("  Output directory: $output_dir")
    println("  Check logs and critical points for debugging.")
end
println("="^80)

# Exit with appropriate code
exit(passed ? 0 : 1)
