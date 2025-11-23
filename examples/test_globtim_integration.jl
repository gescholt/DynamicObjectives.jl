#!/usr/bin/env julia
"""
Test integration with globtimcore and globtimpostprocessing (Phase 2 compatible)

This script tests the 2-stage workflow:
1. Use Globtim to find raw critical points
2. Use GlobtimPostProcessing to refine them

Requirements:
- Globtim and GlobtimPostProcessing must be available
- Run ./setup_dev_packages.jl first to set up local dev versions

Phase 2 Updates:
- Uses ExperimentParams instead of StandardExperimentConfig
- Uses keyword arguments for run_standard_experiment
- Compatible with new result schema
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra

# Include ExperimentCLI module for config
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

println("="^80)
println("Testing Globtim Integration Workflow (Phase 2)")
println("="^80)
println()

# ==============================================================================
# Step 1: Create test objective from Dynamic_objectives
# ==============================================================================

println("Step 1: Creating test objective...")
println("-"^80)

# Use simple 2D Lotka-Volterra model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# True parameters
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Create objective function (1-argument format)
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    return_inf_on_error = true,
    eval_timeout = 10.0
)

# Verify objective works
obj_at_true = objective(p_true)
println("  ✓ Objective created")
println("  Objective at true params: $(round(obj_at_true, digits=8))")
@assert obj_at_true < 1e-6 "Objective should be near zero at true parameters"
println()

# ==============================================================================
# Step 2: Run globtimcore to get raw critical points
# ==============================================================================

println("Step 2: Finding raw critical points with Globtim...")
println("-"^80)

# Configure experiment (Phase 2: use ExperimentParams)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 50,                # Grid size
    degree_range = 4:8,     # Polynomial degrees
    max_time = 3600.0,
    basis = :chebyshev
)

# Create output directory
output_dir = mkpath(joinpath(dirname(@__DIR__), "test_results", "integration_test"))

# Run standard experiment (Phase 2: keyword arguments)
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict("experiment_type" => "integration_test"),
    true_params = p_true
)

println("  ✓ Globtim experiment complete")
println("  Output directory: $(result[:output_dir])")
println("  Degrees processed: $(result[:degrees_processed])")
println("  Critical points found: $(result[:total_critical_points])")

# Get best raw value across all degrees
best_raw_value = minimum([dr.best_objective for dr in result[:degree_results]])
println("  Best raw value: $(round(best_raw_value, digits=6))")
println()

# ==============================================================================
# Step 3: Refine critical points with globtimpostprocessing
# ==============================================================================

println("Step 3: Refining critical points with GlobtimPostProcessing...")
println("-"^80)

# Configure refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    verbose = false
)

# Refine results
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)

println()
println("  ✓ Refinement complete")
println("  Converged: $(refined[:n_converged])/$(refined[:n_raw])")
println("  Mean improvement: $(round(refined[:mean_improvement], digits=2))x")
println("  Best refined value: $(round(refined[:best_refined_value], digits=8))")
println()

# ==============================================================================
# Step 4: Verify results
# ==============================================================================

println("Step 4: Verifying results...")
println("-"^80)

# Get best refined parameters
best_params = refined[:refined_points][refined[:best_refined_idx]]
recovery_error = norm(best_params .- p_true) / norm(p_true)

println("  True parameters:    $p_true")
println("  Refined parameters: $(round.(best_params, digits=4))")
println("  Recovery error:     $(round(100*recovery_error, digits=2))%")
println()

# Success criteria
if recovery_error < 0.01
    println("  ✅ SUCCESS: Recovery error < 1%")
elseif recovery_error < 0.05
    println("  ⚠️  PARTIAL: Recovery error < 5%")
else
    println("  ❌ FAIL: Recovery error > 5%")
end

println()
println("="^80)
println("Integration Test Complete!")
println("="^80)
