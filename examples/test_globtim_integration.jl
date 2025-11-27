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

display_section("Globtim Integration Test (Phase 3)")

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

display_section("Step 1: Model Setup and Verification")

# Use simple 2D Lotka-Volterra model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# True parameters
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.5, 1.5), (0.25, 0.75)]  # Tighter domain around true params

# Display model information
display_model_summary(model)

# Display true parameters
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)

# Create objective function (1-argument format)
# Note: eval_timeout removed - it adds 30x overhead per call due to @async/timedwait
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 300,
    L2_norm, first, nothing;
    return_inf_on_error = true
)

# Verify objective works at true parameters
obj_at_true = objective(p_true)
println("Objective created and verified")
println("  Objective at true params: $(round(obj_at_true, digits=8))")
@assert obj_at_true < 1e-6 "Objective should be near zero at true parameters"

# ==============================================================================
# Step 2: Configure and Run globtimcore
# ==============================================================================

display_section("Step 2: Configure Globtim Experiment")

# Configure experiment (Phase 2: use ExperimentParams)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 20,                # Grid size (20^2 = 400 points for 2D)
    degree_range = 8:12,    # Polynomial degrees (fewer degrees)
    max_time = 600.0,       # 10 minute max
    basis = :chebyshev
)

display_results([
    "Grid size (GN)" => "$(config.GN) ($(config.GN^2) points for 2D)",
    "Degree range" => "$(config.degree_range)",
    "Domain size" => config.domain_size,
    "Basis" => "$(config.basis)"
], title="Configuration")

# Create output directory
output_dir = mkpath(joinpath(dirname(@__DIR__), "test_results", "integration_test"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

display_section("Step 3: Stage 1 - Find Raw Critical Points (globtimcore)")

# Run standard experiment (Phase 2: keyword arguments)
t_start = time()
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String, Any}("experiment_type" => "integration_test"),
    true_params = p_true
)
stage1_time = time() - t_start

# Get best raw value across all degrees
best_raw_value = minimum([dr.best_objective for dr in result[:degree_results]])

display_results([
    "Output directory" => result[:output_dir],
    "Degrees processed" => result[:degrees_processed],
    "Total critical points" => result[:total_critical_points],
    "Best raw value" => best_raw_value,
    "Time elapsed (s)" => stage1_time
], title="Stage 1 Complete")

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

display_section("Step 4: Stage 2 - Refine Critical Points (globtimpostprocessing)")

# Configure refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 10.0,  # Reduced from 30s
    show_progress = true        # Show progress so we can see what's happening
)

display_results([
    "Max time per point" => "$(refinement_config.max_time_per_point)s",
    "Method" => "BFGS local optimization"
], title="Refinement Configuration")

# Refine results
t_start = time()
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
stage2_time = time() - t_start

display_results([
    "Raw points" => refined.n_raw,
    "Converged" => refined.n_converged,
    "Success rate" => "$(round(100*refined.n_converged/refined.n_raw, digits=1))%",
    "Mean improvement" => "$(round(refined.mean_improvement, digits=2))x",
    "Best refined value" => refined.best_refined_value,
    "Time elapsed (s)" => stage2_time
], title="Stage 2 Complete")

# ==============================================================================
# Step 5: Parameter Recovery Verification
# ==============================================================================

display_section("Step 5: Parameter Recovery Verification")

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

    display_results([
        "Relative error" => "$(round(100*recovery_error, digits=2))%",
        "Objective value" => refined.best_refined_value,
        "Total time" => "$(round(stage1_time + stage2_time, digits=1))s",
        "Stage 1/Stage 2 ratio" => round(stage1_time/stage2_time, digits=2)
    ], title="Recovery Metrics")

    # Success criteria
    println()
    if recovery_error < 0.01
        println("SUCCESS: Excellent recovery (< 1% error)")
    elseif recovery_error < 0.05
        println("SUCCESS: Good recovery (< 5% error)")
    else
        println("Needs improvement (> 5% error)")
    end
else
    println("Refinement did not converge for any of the $(refined.n_raw) raw critical points")
    println("  Total time: $(round(stage1_time + stage2_time, digits=1))s")
end

display_section("Integration Test Complete!")
