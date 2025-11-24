#!/usr/bin/env julia
"""
Simple test of the 2-stage workflow:

Stage 1 (globtimcore): Find raw critical points
Stage 2 (globtimpostprocessing): Refine critical points with local optimization

This demonstrates the Phase 2 migration architecture where refinement
is separated from critical point finding.
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using PrettyTables
using Printf

# Include ExperimentCLI module to access ExperimentParams
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

# ==============================================================================
# Display Helper Functions
# ==============================================================================

"""Display a workflow section header"""
function print_section(title::String)
    println("\n" * "="^80)
    println(title)
    println("="^80)
end

"""Display a subsection header"""
function print_subsection(title::String)
    println("\n" * title)
    println("-"^80)
end

"""Display stage results as a formatted table"""
function display_stage_table(stage_name::String, results::Vector{Pair{String, Any}}; color=:cyan)
    data = Matrix{Any}(undef, length(results), 2)
    for (i, (key, value)) in enumerate(results)
        data[i, 1] = key
        if value isa AbstractFloat
            if abs(value) < 1e-3 || abs(value) > 1e4
                data[i, 2] = @sprintf("%.4e", value)
            else
                data[i, 2] = round(value, digits=4)
            end
        else
            data[i, 2] = value
        end
    end

    println()
    pretty_table(
        data,
        header=["Metric", "Value"],
        header_crayon=crayon"bold $color",
        border_crayon=crayon"$color",
        alignment=[:l, :r]
    )
end

# ==============================================================================
# Main Workflow
# ==============================================================================

print_section("2-STAGE WORKFLOW TEST")
println("""
This test demonstrates the complete workflow:
  1. Model Setup & Verification
  2. Experiment Configuration
  3. Stage 1: Find raw critical points (globtimcore)
  4. Stage 2: Refine with local optimization (globtimpostprocessing)
  5. Parameter Recovery Verification
""")

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

print_section("Step 1: Model Setup and Verification")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

display_model_summary(model)

display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)

# Create and verify objective function
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    eval_timeout = 10.0
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

print_section("Step 2: Experiment Configuration")

config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,
    degree_range = 4:6,
    max_time = 3600.0,
    basis = :chebyshev
)

config_data = [
    "Grid size (GN)" "$(config.GN) ($(config.GN^2) points for 2D)";
    "Degree range" "$(config.degree_range)";
    "Domain size" "$(config.domain_size)";
    "Basis" "$(config.basis)";
    "Max time" "$(config.max_time)s"
]

println()
pretty_table(
    config_data,
    header=["Parameter", "Value"],
    header_crayon=crayon"bold blue",
    border_crayon=crayon"blue"
)

output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "simple_workflow"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

print_section("Step 3: Stage 1 - Find Raw Critical Points")
println("Running globtimcore to search for critical points...")

t_start = time()
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String, Any}("experiment_type" => "lv2d_simple_test"),
    true_params = p_true
)
stage1_time = time() - t_start

stage1_results = [
    "Degrees processed" => result[:degrees_processed],
    "Total critical points" => result[:total_critical_points],
    "Best raw objective" => minimum([r.best_objective for r in result[:degree_results]]),
    "Time elapsed (s)" => stage1_time
]

display_stage_table("Stage 1", stage1_results; color=:yellow)

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

print_section("Step 4: Stage 2 - Refine Critical Points")
println("Running globtimpostprocessing for local optimization...")

refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    show_progress = false
)

println("\nRefinement settings:")
println("  Method: BFGS local optimization")
println("  Max time per point: $(refinement_config.max_time_per_point)s")
println("  Objective: L2 distance from true data")

t_start = time()
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
stage2_time = time() - t_start

stage2_results = [
    "Raw points" => refined.n_raw,
    "Converged points" => refined.n_converged,
    "Success rate (%)" => 100.0 * refined.n_converged / refined.n_raw,
    "Mean improvement" => refined.mean_improvement,
    "Best refined objective" => refined.best_refined_value,
    "Time elapsed (s)" => stage2_time
]

display_stage_table("Stage 2", stage2_results; color=:green)

# ==============================================================================
# Step 5: Verification
# ==============================================================================

print_section("Step 5: Parameter Recovery Verification")

if refined.n_converged > 0
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    # Display parameter comparison
    display_parameters(
        [Symbol("α"), Symbol("β")],
        hcat(p_true, best_params),
        labels=["True", "Recovered"]
    )

    # Display recovery metrics
    metrics_data = [
        "Relative error" "$(@sprintf("%.2f", 100*recovery_error))%";
        "Objective value" "$(@sprintf("%.2e", refined.best_refined_value))";
        "Stage 1 time" "$(round(stage1_time, digits=2))s";
        "Stage 2 time" "$(round(stage2_time, digits=2))s";
        "Total time" "$(round(stage1_time + stage2_time, digits=2))s"
    ]

    println()
    pretty_table(
        metrics_data,
        header=["Metric", "Value"],
        header_crayon=crayon"bold blue",
        border_crayon=crayon"blue",
        alignment=[:l, :r]
    )

    # Final assessment
    println()
    if recovery_error < 0.01
        println("RESULT: Excellent recovery (< 1% relative error)")
        println("The optimization successfully recovered the true parameters with high accuracy.")
    elseif recovery_error < 0.05
        println("RESULT: Good recovery (< 5% relative error)")
        println("The optimization found parameters close to the true values.")
    else
        println("RESULT: Recovery needs improvement (> 5% relative error)")
        println("Consider adjusting experiment configuration:")
        println("  - Increase grid size (GN)")
        println("  - Expand degree range")
        println("  - Adjust refinement settings")
    end
else
    println("WARNING: Refinement did not converge for any of the $(refined.n_raw) raw critical points")
    println("\nPossible actions:")
    println("  - Adjust refinement configuration (max_time_per_point, tolerances)")
    println("  - Review raw critical point quality from Stage 1")
    println("  - Check objective function behavior")
    println("  - Consider alternative optimization methods")
    println("\nTotal time: $(round(stage1_time + stage2_time, digits=1))s")
end

# ==============================================================================
# Summary
# ==============================================================================

print_section("WORKFLOW COMPLETE")
println("All stages executed successfully!")
println()
