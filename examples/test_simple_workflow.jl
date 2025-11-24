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
using Globtim: Globtim, run_standard_experiment  # Load first to bring in ConstructionBase
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using Term
using Printf

# Include ExperimentCLI module to access ExperimentParams
# Note: Must load Globtim first since ExperimentCLI uses ConstructionBase
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

# ==============================================================================
# Helper Functions for Rich Output
# ==============================================================================

"""Display a section header with Term panel"""
function display_section_header(title::String, subtitle::String=""; style::String="bold cyan")
    content = subtitle == "" ? "" : subtitle
    panel = Panel(
        content,
        title=title,
        title_style=style,
        style=style,
        fit=true,
        padding=(0, 2, 0, 2)
    )
    println(panel)
end

"""Display stage results in a formatted panel"""
function display_stage_results(stage_name::String, results::Dict; style::String="green")
    content = ""
    for (key, value) in results
        if value isa AbstractFloat
            if abs(value) < 1e-3 || abs(value) > 1e4
                content *= "  $(rpad(key, 25)): $(@sprintf("%.4e", value))\n"
            else
                content *= "  $(rpad(key, 25)): $(round(value, digits=6))\n"
            end
        else
            content *= "  $(rpad(key, 25)): $value\n"
        end
    end

    panel = Panel(
        content,
        title="✓ $stage_name Complete",
        title_style="bold $style",
        style=style,
        fit=true,
        padding=(1, 2, 1, 2)
    )
    println(panel)
end

"""Display configuration in a formatted table"""
function display_config(config::ExperimentParams)
    content = """
  Grid size (GN)      : $(config.GN) → $(config.GN^2) points for 2D
  Degree range        : $(config.degree_range)
  Domain size         : $(config.domain_size)
  Basis               : $(config.basis)
  Max time            : $(config.max_time)s
"""

    panel = Panel(
        content,
        title="Experiment Configuration",
        title_style="bold blue",
        style="blue",
        fit=true,
        padding=(1, 2, 0, 2)
    )
    println(panel)
end

# ==============================================================================
# Main Workflow
# ==============================================================================

# Title banner
title_panel = Panel(
    "Stage 1: Find raw critical points (globtimcore)\n" *
    "Stage 2: Refine with local optimization (globtimpostprocessing)",
    title="2-Stage Workflow Test",
    title_style="bold magenta",
    style="magenta",
    fit=true,
    padding=(1, 2, 1, 2)
)
println("\n", title_panel, "\n")

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

display_section_header("Step 1: Model Setup & Verification", "Defining the Lotka-Volterra model")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Display model summary
display_model_summary(model)

# Display parameters
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)

# Create objective function
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

# Verify objective works at true parameters
obj_at_true = objective(p_true)
if obj_at_true < 1e-6
    check_panel = Panel(
        "Objective value at true parameters: $(@sprintf("%.2e", obj_at_true))\n" *
        "Model validation passed ✓",
        title="Objective Function Verified",
        title_style="bold green",
        style="green",
        fit=true,
        padding=(0, 2, 0, 2)
    )
    println(check_panel)
else
    error("Objective should be near zero at true parameters, got $obj_at_true")
end
println()

# ==============================================================================
# Step 2: Configure Experiment
# ==============================================================================

display_section_header("Step 2: Experiment Configuration", "Configuring globtimcore parameters")

# Configure globtimcore using ExperimentParams
# Note: GN = 10 gives 10^2 = 100 grid points (fast for testing)
#       For production, use GN = 50-100 for better accuracy
config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,           # Small grid for quick test (100 points for 2D)
    degree_range = 4:6, # Reduced range for speed
    max_time = 3600.0,
    basis = :chebyshev
)

display_config(config)
println()

# Create output directory
output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "simple_workflow"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

display_section_header("Step 3: Stage 1 - Raw Critical Points", "Running globtimcore to find critical points"; style="bold yellow")

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

stage1_results = Dict(
    "Degrees processed" => result[:degrees_processed],
    "Total critical points" => result[:total_critical_points],
    "Best raw value" => minimum([r.best_objective for r in result[:degree_results]]),
    "Time elapsed" => stage1_time
)

display_stage_results("Stage 1", stage1_results; style="yellow")
println()

# ==============================================================================
# Step 4: Stage 2 - Refine Critical Points
# ==============================================================================

display_section_header("Step 4: Stage 2 - Refinement", "Using globtimpostprocessing for local optimization"; style="bold green")

# Configure refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    show_progress = false
)

ref_config_content = """
  Max time per point  : $(refinement_config.max_time_per_point)s
  Method              : BFGS local optimization
  Objective           : L2 distance from true data
"""
ref_config_panel = Panel(
    ref_config_content,
    title="Refinement Configuration",
    title_style="bold cyan",
    style="cyan",
    fit=true,
    padding=(0, 2, 0, 2)
)
println(ref_config_panel)
println()

# Run refinement
t_start = time()
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
stage2_time = time() - t_start

stage2_results = Dict(
    "Raw points" => refined.n_raw,
    "Converged" => refined.n_converged,
    "Success rate" => 100*refined.n_converged/refined.n_raw,
    "Mean improvement" => refined.mean_improvement,
    "Best refined value" => refined.best_refined_value,
    "Time elapsed" => stage2_time
)

display_stage_results("Stage 2", stage2_results; style="green")
println()

# ==============================================================================
# Step 5: Verification
# ==============================================================================

display_section_header("Step 5: Parameter Recovery Verification", "Comparing recovered parameters with true values")

# Check if any points converged
if refined.n_converged > 0
    # Get best refined parameters
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    # Display comparison
    display_parameters(
        [Symbol("α"), Symbol("β")],
        hcat(p_true, best_params),
        labels=["True", "Recovered"]
    )

    # Display recovery metrics
    recovery_metrics = Dict(
        "Relative error (%)" => 100*recovery_error,
        "Objective value" => refined.best_refined_value,
        "Total time (s)" => stage1_time + stage2_time
    )

    metrics_content = ""
    for (key, value) in recovery_metrics
        if value isa AbstractFloat
            if abs(value) < 1e-3 || abs(value) > 1e4
                metrics_content *= "  $(rpad(key, 25)): $(@sprintf("%.4e", value))\n"
            else
                metrics_content *= "  $(rpad(key, 25)): $(round(value, digits=6))\n"
            end
        else
            metrics_content *= "  $(rpad(key, 25)): $value\n"
        end
    end

    metrics_panel = Panel(
        metrics_content,
        title="Recovery Metrics",
        title_style="bold blue",
        style="blue",
        fit=true,
        padding=(1, 2, 1, 2)
    )
    println(metrics_panel)
    println()

    # Final verdict
    verdict_content = ""
    verdict_style = ""
    verdict_title = ""

    if recovery_error < 0.01
        verdict_content = "Excellent recovery (< 1% relative error)\n\nThe optimization successfully recovered the true parameters\nwith high accuracy!"
        verdict_style = "green"
        verdict_title = "✅ Success - Excellent Recovery"
    elseif recovery_error < 0.05
        verdict_content = "Good recovery (< 5% relative error)\n\nThe optimization found parameters close to the true values.\nPerformance is acceptable for most applications."
        verdict_style = "green"
        verdict_title = "✅ Success - Good Recovery"
    else
        verdict_content = "Recovery needs improvement (> 5% relative error)\n\nConsider adjusting experiment configuration:\n- Increase grid size (GN)\n- Expand degree range\n- Adjust refinement settings"
        verdict_style = "yellow"
        verdict_title = "⚠️  Needs Improvement"
    end

    verdict_panel = Panel(
        verdict_content,
        title=verdict_title,
        title_style="bold $verdict_style",
        style=verdict_style,
        fit=true,
        padding=(1, 2, 1, 2)
    )
    println(verdict_panel)
else
    # No points converged - refinement needs improvement
    warning_content = """
  Refinement did not converge for any of the $(refined.n_raw) raw critical points

  Possible actions:
    • Adjust refinement configuration (max_time_per_point, tolerances)
    • Review raw critical point quality from Stage 1
    • Check objective function behavior
    • Consider alternative optimization methods

  Total time: $(round(stage1_time + stage2_time, digits=1))s
"""

    warning_panel = Panel(
        warning_content,
        title="⚠️  No Convergence",
        title_style="bold yellow",
        style="yellow",
        fit=true,
        padding=(1, 2, 1, 2)
    )
    println(warning_panel)
end
println()

# Final summary banner
completion_panel = Panel(
    "All stages completed successfully!\n\n" *
    "Stage 1: Critical point search\n" *
    "Stage 2: Local refinement\n" *
    "Stage 3: Verification & analysis",
    title="Workflow Complete",
    title_style="bold magenta",
    style="magenta",
    fit=true,
    padding=(1, 2, 1, 2)
)
println(completion_panel)
println()
