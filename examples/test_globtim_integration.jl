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
println("Globtim Integration Test (Phase 3)")
println("="^80)
println()

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

println("Step 1: Model Setup and Verification")
println("-"^80)

# Use simple 2D Lotka-Volterra model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# True parameters
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Display model information
display_model_summary(model)
println()

# Display true parameters
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)
println()

# Create objective function (1-argument format)
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    return_inf_on_error = true,
    eval_timeout = 10.0
)

# Verify objective works at true parameters
obj_at_true = objective(p_true)
println("✓ Objective created and verified")
println("  Objective at true params: $(round(obj_at_true, digits=8))")
@assert obj_at_true < 1e-6 "Objective should be near zero at true parameters"
println()

# ==============================================================================
# Step 2: Configure and Run globtimcore
# ==============================================================================

println("Step 2: Configure Globtim Experiment")
println("-"^80)

# Configure experiment (Phase 2: use ExperimentParams)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 50,                # Grid size
    degree_range = 4:8,     # Polynomial degrees
    max_time = 3600.0,
    basis = :chebyshev
)

println("Configuration:")
println("  Grid size (GN): $(config.GN) ($(config.GN^2) points for 2D)")
println("  Degree range: $(config.degree_range)")
println("  Domain size: $(config.domain_size)")
println("  Basis: $(config.basis)")
println()

# Create output directory
output_dir = mkpath(joinpath(dirname(@__DIR__), "test_results", "integration_test"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

println("Step 3: Stage 1 - Find Raw Critical Points (globtimcore)")
println("-"^80)

# Run standard experiment (Phase 2: keyword arguments)
t_start = time()
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict("experiment_type" => "integration_test"),
    true_params = p_true
)
stage1_time = time() - t_start

# Get best raw value across all degrees
best_raw_value = minimum([dr.best_objective for dr in result[:degree_results]])

println()
println("✓ Stage 1 Complete")
println("  Output directory: $(result[:output_dir])")
println("  Degrees processed: $(result[:degrees_processed])")
println("  Total critical points: $(result[:total_critical_points])")
println("  Best raw value: $(round(best_raw_value, digits=6))")
println("  Time elapsed: $(round(stage1_time, digits=1))s")
println()

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

println("Step 4: Stage 2 - Refine Critical Points (globtimpostprocessing)")
println("-"^80)

# Configure refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    show_progress = false
)

println("Refinement configuration:")
println("  Max time per point: $(refinement_config.max_time_per_point)s")
println("  Method: BFGS local optimization")
println()

# Refine results
t_start = time()
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
stage2_time = time() - t_start

println()
println("✓ Stage 2 Complete")
println("  Raw points: $(refined.n_raw)")
println("  Converged: $(refined.n_converged)")
println("  Success rate: $(round(100*refined.n_converged/refined.n_raw, digits=1))%")
println("  Mean improvement: $(round(refined.mean_improvement, digits=2))x")
println("  Best refined value: $(round(refined.best_refined_value, digits=8))")
println("  Time elapsed: $(round(stage2_time, digits=1))s")
println()

# ==============================================================================
# Step 5: Parameter Recovery Verification
# ==============================================================================

println("Step 5: Parameter Recovery Verification")
println("-"^80)

# Check if any points converged
if refined.n_converged > 0
    # Get best refined parameters
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    # Display parameter comparison
    display_parameters(
        [Symbol("α"), Symbol("β")],
        hcat(p_true, best_params),
        labels=["True", "Recovered"]
    )
    println()

    println("Recovery metrics:")
    println("  Relative error: $(round(100*recovery_error, digits=2))%")
    println("  Objective value: $(round(refined.best_refined_value, digits=8))")
    println("  Total time: $(round(stage1_time + stage2_time, digits=1))s")
    println("  Stage 1/Stage 2 ratio: $(round(stage1_time/stage2_time, digits=2))")
    println()

    # Success criteria
    if recovery_error < 0.01
        println("✅ SUCCESS: Excellent recovery (< 1% error)")
    elseif recovery_error < 0.05
        println("✅ SUCCESS: Good recovery (< 5% error)")
    else
        println("⚠️  Needs improvement (> 5% error)")
    end
else
    # No points converged - refinement needs improvement
    # TODO: Investigate postprocessing step - understand optimization methods,
    #       convergence criteria, and alternative refinement approaches
    println("⚠️  Refinement did not converge for any of the $(refined.n_raw) raw critical points")
    println("  Total time: $(round(stage1_time + stage2_time, digits=1))s")
end
println()

println("="^80)
println("Integration Test Complete!")
println("="^80)
