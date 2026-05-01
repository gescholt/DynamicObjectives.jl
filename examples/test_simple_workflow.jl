#!/usr/bin/env julia
"""
Simple test of the 2-stage workflow:

Stage 1 (Globtim): Find raw critical points
Stage 2 (globtimpostprocessing): Refine critical points with local optimization

This demonstrates the Phase 2 migration architecture where refinement
is separated from critical point finding.
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing:
    GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using Printf

# Include ExperimentCLI module to access ExperimentParams
if !isdefined(Main, :ExperimentCLI)
    globtim_path = joinpath(dirname(@__DIR__), "..", "globtim")
    include(joinpath(globtim_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

# ==============================================================================
# Main Workflow
# ==============================================================================

display_section("2-STAGE WORKFLOW TEST", subtitle = """
  This test demonstrates the complete workflow:
    1. Model Setup & Verification
    2. Experiment Configuration
    3. Stage 1: Find raw critical points (Globtim)
    4. Stage 2: Refine with local optimization (globtimpostprocessing)
    5. Parameter Recovery Verification""")

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

display_section("Step 1: Model Setup and Verification")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

display_model_summary(model)

display_parameters([Symbol("α"), Symbol("β")], hcat(p_true), labels = ["True Values"])

# Create and verify objective function
# Note: eval_timeout removed - it adds 30x overhead per call due to @async/timedwait
objective = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    [0.0, 20.0],
    120,
    L2_norm;
    return_inf_on_error = true,
)

obj_at_true = objective(p_true)
println("Objective at true parameters: $(@sprintf("%.2e", obj_at_true))")
if obj_at_true < 1e-6
    println("Model validation: PASSED")
else
    error("Objective should be near zero at true parameters, got $obj_at_true")
end

# ==============================================================================
# Step 2: Configure Experiment
# ==============================================================================

display_section("Step 2: Experiment Configuration")

config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,
    degree_range = 4:6,
    max_time = 3600.0,
    basis = :chebyshev,
)

display_results(
    [
        "Grid size (GN)" => "$(config.GN) ($(config.GN^2) points for 2D)",
        "Degree range" => "$(config.degree_range)",
        "Domain size" => config.domain_size,
        "Basis" => "$(config.basis)",
        "Max time" => "$(config.max_time)s",
    ],
    title = "Experiment Configuration",
)

output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "simple_workflow"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

display_section("Step 3: Stage 1 - Find Raw Critical Points")
display_subsection("Running Globtim to search for critical points...")

t_start = time()
result = run_standard_experiment(
    objective_function = objective,
    objective_name = "lv2d_simple_test",
    bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String,Any}("experiment_type" => "lv2d_simple_test"),
    true_params = p_true,
)
stage1_time = time() - t_start

display_results(
    [
        "Degrees processed" => result[:degrees_processed],
        "Total critical points" => result[:total_critical_points],
        "Best raw objective" =>
            minimum([r.best_objective for r in result[:degree_results]]),
        "Time elapsed (s)" => stage1_time,
    ],
    title = "Stage 1 Results",
)

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

display_section("Step 4: Stage 2 - Refine Critical Points")
display_subsection("Running globtimpostprocessing for local optimization...")

refinement_config = ode_refinement_config(max_time_per_point = 30.0, show_progress = false)

display_results(
    [
        "Method" => "BFGS local optimization",
        "Max time per point" => "$(refinement_config.max_time_per_point)s",
        "Objective" => "L2 distance from true data",
    ],
    title = "Refinement Settings",
)

t_start = time()
refined = refine_experiment_results(result[:output_dir], objective, refinement_config)
stage2_time = time() - t_start

display_results(
    [
        "Raw points" => refined.n_raw,
        "Converged points" => refined.n_converged,
        "Success rate (%)" => 100.0 * refined.n_converged / refined.n_raw,
        "Mean improvement" => refined.mean_improvement,
        "Best refined objective" => refined.best_refined_value,
        "Time elapsed (s)" => stage2_time,
    ],
    title = "Stage 2 Results",
)

# ==============================================================================
# Step 5: Verification
# ==============================================================================

display_section("Step 5: Parameter Recovery Verification")

if refined.n_converged > 0
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    # Display parameter comparison
    display_parameters(
        [Symbol("α"), Symbol("β")],
        hcat(p_true, best_params),
        labels = ["True", "Recovered"],
    )

    # Display recovery metrics
    display_results(
        [
            "Relative error" => "$(@sprintf("%.2f", 100*recovery_error))%",
            "Objective value" => refined.best_refined_value,
            "Stage 1 time" => "$(round(stage1_time, digits=2))s",
            "Stage 2 time" => "$(round(stage2_time, digits=2))s",
            "Total time" => "$(round(stage1_time + stage2_time, digits=2))s",
        ],
        title = "Recovery Metrics",
    )

    # Final assessment
    println()
    if recovery_error < 0.01
        println("RESULT: Excellent recovery (< 1% relative error)")
    elseif recovery_error < 0.05
        println("RESULT: Good recovery (< 5% relative error)")
    else
        println("RESULT: Recovery needs improvement (> 5% relative error)")
        println("Consider: increase GN, expand degree_range, adjust refinement settings")
    end
else
    println(
        "WARNING: Refinement did not converge for any of the $(refined.n_raw) raw critical points",
    )
    println("\nPossible actions:")
    println("  - Adjust refinement configuration (max_time_per_point, tolerances)")
    println("  - Review raw critical point quality from Stage 1")
    println("  - Check objective function behavior")
    println("\nTotal time: $(round(stage1_time + stage2_time, digits=1))s")
end

# ==============================================================================
# Summary
# ==============================================================================

display_section("WORKFLOW COMPLETE")
println("All stages executed successfully!")
