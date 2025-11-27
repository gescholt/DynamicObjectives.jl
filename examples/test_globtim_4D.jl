#!/usr/bin/env julia
"""
Test 4D parameter estimation using constrained Lotka-Volterra model

This script tests the 2-stage workflow on a 4-dimensional problem:
1. Use Globtim to find raw critical points
2. Use GlobtimPostProcessing to refine them
3. Validate with gradient analysis

Model: Constrained Lotka-Volterra 4D
- 4 parameters to estimate
- 4 state variables
- Higher computational complexity than 2D

Requirements:
- Globtim and GlobtimPostProcessing must be available
- Run ./setup_dev_packages.jl first to set up local dev versions
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using ForwardDiff

# Include ExperimentCLI module for config
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

display_section("4D Parameter Estimation Test", subtitle="""
Testing globtim pipeline on 4-dimensional parameter estimation problem
Model: Constrained Lotka-Volterra 4D
Dimensions: 4 parameters, 4 states""")

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

display_section("Step 1: Model Setup and Verification")

# Use constrained Lotka-Volterra 4D model
model, params, states, outputs = define_constrained_lotka_volterra_4D()

# True parameters (4D)
p_true = [1.0, 0.5, 1.0, 0.5]
ic = [1.0, 0.5, 0.5, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 3.0), (0.0, 2.0)]

# Display model information
display_model_summary(model)

# Display true parameters
param_names = [Symbol("p$i") for i in 1:4]
display_parameters(
    param_names,
    hcat(p_true),
    labels=["True Values"]
)

# Create objective function
# Note: Use more time points for 4D to get better signal
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 120,
    L2_norm, first, nothing;
    return_inf_on_error = true
)

# Verify objective works at true parameters
obj_at_true = objective(p_true)
println("Objective created and verified")
println("  Objective at true params: $(round(obj_at_true, digits=8))")
@assert obj_at_true < 1e-6 "Objective should be near zero at true parameters"

# ==============================================================================
# Step 2: Configure Globtim Experiment
# ==============================================================================

display_section("Step 2: Configure Globtim Experiment")

# Configure experiment for 4D
# Note: GN=10 gives 10^4 = 10,000 grid points (much denser than 2D's 400)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,                # 10^4 = 10,000 grid points for 4D
    degree_range = 6:10,    # Lower degrees for 4D (computational complexity)
    max_time = 1800.0,      # 30 minutes max
    basis = :chebyshev
)

display_results([
    "Grid size (GN)" => "$(config.GN) ($(config.GN^4) points for 4D)",
    "Degree range" => "$(config.degree_range)",
    "Domain size" => config.domain_size,
    "Basis" => "$(config.basis)",
    "Max time" => "$(config.max_time)s (30 min)"
], title="4D Configuration")

# Create output directory
output_dir = mkpath(joinpath(dirname(@__DIR__), "test_results", "4D_test"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

display_section("Step 3: Stage 1 - Find Raw Critical Points (4D)")
display_subsection("Running globtimcore on 4D problem (this may take a while)...")

# Run standard experiment
t_start = time()
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String, Any}("experiment_type" => "4D_test", "model" => "constrained_LV_4D"),
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
    "Time elapsed (s)" => stage1_time,
    "Time elapsed (min)" => round(stage1_time / 60, digits=2)
], title="Stage 1 Complete")

# Show degree comparison
degree_comparison = [(
    degree = dr.degree,
    n_critical_points = dr.n_critical_points,
    best_objective = dr.best_objective
) for dr in result[:degree_results]]
display_degree_comparison(degree_comparison)

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

display_section("Step 4: Stage 2 - Refine Critical Points (4D)")
display_subsection("Running globtimpostprocessing for local optimization...")

# Configure refinement (more time for 4D)
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,  # More time for 4D optimization
    show_progress = true        # Show progress
)

display_results([
    "Max time per point" => "$(refinement_config.max_time_per_point)s",
    "Method" => "BFGS local optimization",
    "Objective" => "L2 distance from true data"
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
    "Time elapsed (s)" => stage2_time,
    "Time elapsed (min)" => round(stage2_time / 60, digits=2)
], title="Stage 2 Complete")

# Display quality summary
display_quality_summary(refined, title="Critical Point Quality (4D)")

# ==============================================================================
# Step 4.5: Gradient Validation
# ==============================================================================

display_section("Step 4.5: Gradient Validation (4D)")

# Compute gradient norms for all converged points
if refined.n_converged > 0
    println("Computing gradient norms for $(refined.n_converged) converged points...")

    grad_norms = Float64[]
    for i in 1:refined.n_converged
        params = refined.refined_points[i]
        grad = ForwardDiff.gradient(objective, params)
        grad_norm = norm(grad)
        push!(grad_norms, grad_norm)
    end

    # Display gradient analysis
    display_gradient_analysis(grad_norms, tolerance=1e-6)

    # Show best point gradient
    best_grad = ForwardDiff.gradient(objective, refined.refined_points[refined.best_refined_idx])
    println("Best point gradient: ||∇f|| = $(norm(best_grad))")
else
    println("No converged points to validate")
end

# ==============================================================================
# Step 5: Parameter Recovery Verification (4D)
# ==============================================================================

display_section("Step 5: Parameter Recovery Verification (4D)")

# Check if any points converged
if refined.n_converged > 0
    # Get best refined parameters
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    # Display parameter comparison
    display_parameters(
        param_names,
        hcat(p_true, best_params),
        labels=["True", "Recovered"]
    )

    # Compute per-parameter errors
    param_errors = abs.(best_params .- p_true) ./ abs.(p_true)

    display_results([
        "Relative error (overall)" => "$(round(100*recovery_error, digits=2))%",
        "Parameter 1 error" => "$(round(100*param_errors[1], digits=2))%",
        "Parameter 2 error" => "$(round(100*param_errors[2], digits=2))%",
        "Parameter 3 error" => "$(round(100*param_errors[3], digits=2))%",
        "Parameter 4 error" => "$(round(100*param_errors[4], digits=2))%",
        "Objective value" => refined.best_refined_value,
        "Total time (min)" => round((stage1_time + stage2_time) / 60, digits=2),
        "Stage 1/Stage 2 ratio" => round(stage1_time/stage2_time, digits=2)
    ], title="4D Recovery Metrics")

    # Success criteria
    println()
    if recovery_error < 0.01
        println("SUCCESS: Excellent 4D recovery (< 1% error)")
    elseif recovery_error < 0.05
        println("SUCCESS: Good 4D recovery (< 5% error)")
    elseif recovery_error < 0.10
        println("ACCEPTABLE: Reasonable 4D recovery (< 10% error)")
        println("Note: 4D problems are significantly harder than 2D")
    else
        println("Needs improvement (> 10% error)")
        println("Consider:")
        println("  - Increase grid size (GN)")
        println("  - Expand degree range")
        println("  - Increase refinement time per point")
        println("  - Add more time points in objective")
    end
else
    println("Refinement did not converge for any of the $(refined.n_raw) raw critical points")
    println("\nPossible actions:")
    println("  - Check if degrees are appropriate for 4D")
    println("  - Increase max_time_per_point")
    println("  - Review raw critical point quality")
    println("\nTotal time: $(round((stage1_time + stage2_time) / 60, digits=2)) minutes")
end

# ==============================================================================
# Summary
# ==============================================================================

display_section("4D Test Complete!")
println("4D parameter estimation test finished")
println("Total time: $(round((stage1_time + stage2_time) / 60, digits=2)) minutes")
println()
