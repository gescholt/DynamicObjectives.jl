#!/usr/bin/env julia
# Quick Start: Your First Test with globtim
# Usage: julia --project=. examples/quick_start.jl
#
# Prerequisites: Run setup_dev_packages.jl once to install globtimcore

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using Printf

# Load ExperimentCLI for configuration
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

display_section("Quick Start: Testing Dynamic_objectives with globtim")

# ==============================================================================
# STEP 1: Choose a Model (starting with easiest)
# ==============================================================================

display_section("Step 1: Setting up model")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

display_results([
    "Model" => "LV 2D v3 with 2 outputs (EASIEST)",
    "Parameters" => 2,
    "Outputs" => "2 (full observability)",
    "True params" => "[$(p_true[1]), $(p_true[2])]",
    "Bounds" => "[($(bounds[1][1]), $(bounds[1][2])), ($(bounds[2][1]), $(bounds[2][2]))]"
], title="Model Configuration")

# ==============================================================================
# STEP 2: Create Objective Function
# ==============================================================================

display_section("Step 2: Creating objective function")

error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    [0.0, 20.0],  # time interval
    30,           # number of time points
    L2_norm,      # distance metric
    first,        # aggregation function
    nothing;      # no noise
    return_inf_on_error = true,
    eval_timeout = nothing
)

println("Objective function created")

# ==============================================================================
# STEP 3: Verify Setup
# ==============================================================================

display_section("Step 3: Verifying setup")

error_at_true = error_func(p_true)
@printf "Error at true parameters: %.6e\n" error_at_true

if error_at_true < 1e-6
    println("Setup verified (error near 0)")
else
    error("Setup failed - error should be near 0 at true parameters")
end

# Test a few points
display_subsection("Testing objective function at different points")
test_points = [
    [1.0, 0.5],    # true
    [1.1, 0.5],    # perturb p1
    [1.0, 0.6],    # perturb p2
    [0.5, 0.5],    # far from true
    [2.0, 1.0],    # far from true
]

for p in test_points
    err = error_func(p)
    @printf "  p = [%.1f, %.1f]  ->  error = %.6e\n" p[1] p[2] err
end

println("\nObjective function responds correctly to parameter changes")

# ==============================================================================
# STEP 4: Run globtimcore Optimization
# ==============================================================================

display_section("Step 4: Running globtimcore optimization")

# Configure experiment (fast settings for quick start)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 12,              # Grid resolution (fast)
    degree_range = 4:6,   # Polynomial degrees
    max_time = 3600.0,
    basis = :chebyshev
)

# Create output directory
output_dir = mktempdir()

display_results([
    "GN" => config.GN,
    "Degree range" => "$(config.degree_range)",
    "Output dir" => output_dir
], title="Experiment Config")

# Stage 1: Find raw critical points using polynomial approximation + HomotopyContinuation
display_subsection("Stage 1: Finding critical points via polynomial approximation")
raw_result = run_standard_experiment(
    objective_function = error_func,
    problem_params = nothing,  # 1-arg function
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String, Any}("experiment_type" => "quick_start_lv2d"),
    true_params = p_true
)

n_raw = raw_result[:total_critical_points]
best_raw = minimum([r.best_objective for r in raw_result[:degree_results]])

display_results([
    "Critical points found" => n_raw,
    "Degrees processed" => raw_result[:degrees_processed],
    "Best raw objective" => best_raw
], title="Stage 1 Results")

# Stage 2: Refine critical points using local BFGS optimization
display_subsection("Stage 2: Refining critical points with BFGS")
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    f_abstol = 1e-10,
    show_progress = false
)

refined_result = refine_experiment_results(
    raw_result[:output_dir],
    error_func,
    refinement_config
)

display_results([
    "Refined points" => "$(refined_result.n_converged)/$(refined_result.n_raw)",
    "Best refined value" => refined_result.best_refined_value,
    "Mean improvement" => "$(round(refined_result.mean_improvement, digits=1))x"
], title="Stage 2 Results")

# ==============================================================================
# STEP 5: Analyze Results
# ==============================================================================

display_section("Step 5: Results Analysis")

# Extract best parameters
if refined_result.n_converged > 0
    p_best = refined_result.refined_points[refined_result.best_refined_idx]
    error_best = refined_result.best_refined_value
    recovery_error = norm(p_best .- p_true) / norm(p_true)
    used_refined = true
else
    # Fall back to best raw point if no refinement converged
    raw_points = Vector{Float64}[]
    for dr in raw_result[:degree_results]
        append!(raw_points, dr.critical_points)
    end
    if isempty(raw_points)
        error("No critical points found in any degree result")
    end
    raw_errors = [error_func(p) for p in raw_points]
    best_raw_local_idx = argmin(raw_errors)
    p_best = raw_points[best_raw_local_idx]
    error_best = raw_errors[best_raw_local_idx]
    recovery_error = norm(p_best .- p_true) / norm(p_true)
    used_refined = false
    println("WARNING: No refinement converged, using best raw critical point")
end

display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true, p_best),
    labels=["True", "Recovered"]
)

display_results([
    "Objective at solution" => error_best,
    "Relative param error" => "$(round(100*recovery_error, digits=4))%",
    "Source" => used_refined ? "refined" : "raw",
    "Raw critical points" => n_raw,
    "Refined points" => refined_result.n_converged,
    "Stage 1 time (s)" => raw_result[:total_time]
], title="Recovery Summary")

if recovery_error < 0.01
    println("\nSUCCESS - Excellent parameter recovery (<1% error)")
elseif recovery_error < 0.05
    println("\nGood parameter recovery (<5% error)")
else
    println("\nPoor recovery - consider increasing GN or degree_range")
end
