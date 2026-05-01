#!/usr/bin/env julia
"""
Integration Template: 2-Stage Pipeline

This template demonstrates a BALANCED approach between simple and advanced:
- Stage 1: run_standard_experiment (raw critical points via globtim)
- Stage 2: refine_experiment_results (NM refinement via globtimpostprocessing)
- Good for production workflows

USE WHEN:
- You want reasonable control without too much boilerplate
- You're comfortable with config structs
- You want both raw and refined results
- This is a production workflow, not a quick test

ARCHITECTURE:
┌──────────────────────────────────────────────────────────────┐
│ Stage 1: run_standard_experiment (raw critical points)      │
│ Stage 2: refine_experiment_results (NM refinement)          │
│ Returns: (raw_result, refined_result)                       │
└──────────────────────────────────────────────────────────────┘
"""

using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))  # Activate Dynamic_objectives environment

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing:
    GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra

# Include ExperimentCLI for config creation
if !isdefined(Main, :ExperimentCLI)
    globtim_path = joinpath(dirname(dirname(@__DIR__)), "..", "globtim")
    include(joinpath(globtim_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

#===============================================================================
STEP 1: Define your model and create objective function
===============================================================================#

# Example: 2D Lotka-Volterra model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]  # True parameters to recover
ic = [1.0, 0.5]      # Initial condition
bounds = [(0.0, 3.0), (0.0, 2.0)]  # Parameter search bounds

# Create error function (1-arg function: p -> error)
objective = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    [0.0, 20.0],
    30,  # Time span and number of points
    L2_norm;
    eval_timeout = 10.0,
)

#===============================================================================
STEP 2: Create configuration
===============================================================================#

# Pipeline wrapper accepts ExperimentParams (or NamedTuple with same fields)
config = ExperimentParams(
    domain_size = 1.5,
    GN = 50,
    degree_range = 4:8,
    max_time = 3600.0,
    basis = :chebyshev,
)

# Optional: Create custom refinement config
refinement_config =
    ode_refinement_config(max_time_per_point = 30.0, f_tol = 1e-10, verbose = false)

output_dir = mkpath(joinpath(@__DIR__, "..", "..", "test_results", "template_pipeline"))

println("="^80)
println("Pipeline Integration Template")
println("="^80)

#===============================================================================
STEP 3: Run pipeline (both stages in one call)
===============================================================================#

# NOTE: run_globtim_pipeline() is being updated for Phase 2
# For now, we'll run the stages manually in pipeline style

println("\n" * "="^80)
println("Running 2-stage pipeline")
println("="^80)

# Stage 1: Find raw critical points
println("\nStage 1: Finding raw critical points...")
raw_result = run_standard_experiment(
    objective_function = objective,
    objective_name = "lv2d_pipeline_template",
    bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict("experiment_type" => "lv2d_pipeline_template"),
    true_params = p_true,
)

println(
    "  ✓ Found $(raw_result[:total_critical_points]) critical points across $(raw_result[:degrees_processed]) degrees",
)
best_raw = minimum([r.best_objective for r in raw_result[:degree_results]])
println("  Best raw value: $(round(best_raw, digits=6))")

# Stage 2: Refine critical points
println("\nStage 2: Refining critical points...")
refined_result =
    refine_experiment_results(raw_result[:output_dir], objective, refinement_config)

println("  ✓ Refined $(refined_result[:n_converged])/$(refined_result[:n_raw]) points")
println("  Best refined value: $(round(refined_result[:best_refined_value], digits=8))")
println("  Mean improvement: $(round(refined_result[:mean_improvement], digits=2))x")

#===============================================================================
STEP 4: Analyze results
===============================================================================#

println("\n" * "="^80)
println("Results Analysis")
println("="^80)

# Extract best parameters
best_params = refined_result[:refined_points][refined_result[:best_refined_idx]]
recovery_error = norm(best_params .- p_true) / norm(p_true)

println("\nParameter Recovery:")
println("  True parameters:      $p_true")
println("  Recovered parameters: $best_params")
println("  Recovery error:       $(round(100*recovery_error, digits=3))%")

println("\nComputational Summary:")
println("  Raw critical points:   $(raw_result[:total_critical_points])")
println("  Refined points:        $(refined_result[:n_converged])")
println(
    "  Convergence rate:      $(round(100*refined_result[:n_converged]/refined_result[:n_raw], digits=1))%",
)
println("  Stage 1 time:          $(round(raw_result[:total_time], digits=2))s")
println("  Best objective value:  $(round(refined_result[:best_refined_value], digits=8))")

if recovery_error < 0.05
    println("\n✅ SUCCESS - Excellent parameter recovery!")
elseif recovery_error < 0.10
    println("\n✓ Good parameter recovery")
else
    println("\n⚠️  Poor recovery - consider adjusting parameters")
end

#===============================================================================
STEP 5: Access both result structures
===============================================================================#

println("\n" * "="^80)
println("Available Results")
println("="^80)

println("\nRaw result fields:")
println("  :experiment_id       = $(raw_result[:experiment_id])")
println("  :total_critical_points = $(raw_result[:total_critical_points])")
println("  :degrees_processed   = $(raw_result[:degrees_processed])")
println("  :success_rate        = $(round(100*raw_result[:success_rate], digits=1))%")
println(
    "  :degree_results      = Array of DegreeResult ($(length(raw_result[:degree_results])) entries)",
)
println("  :output_dir          = $(raw_result[:output_dir])")

println("\nRefined result fields:")
println("  :n_raw               = $(refined_result[:n_raw])")
println("  :n_converged         = $(refined_result[:n_converged])")
println("  :best_refined_value  = $(round(refined_result[:best_refined_value], digits=8))")
println("  :best_refined_idx    = $(refined_result[:best_refined_idx])")
println(
    "  :refined_points      = Array{Vector{Float64}} ($(length(refined_result[:refined_points])) entries)",
)
println(
    "  :refined_values      = Array{Float64} ($(length(refined_result[:refined_values])) entries)",
)
println("  :mean_improvement    = $(round(refined_result[:mean_improvement], digits=2))x")

#===============================================================================
STEP 6: Save or process results further
===============================================================================#

# You can now:
# 1. Pass results to plotting functions (globtimplots)
# 2. Run further analysis (globtimpostprocessing)
# 3. Compare across multiple experiments
# 4. Export to other formats

# Example: Find all points within 10% of optimal
threshold = refined_result[:best_refined_value] * 1.10
good_indices = findall(v -> v <= threshold, refined_result[:refined_values])
println("\nPoints within 10% of optimal: $(length(good_indices))")

if length(good_indices) > 1
    println("Alternative solutions found:")
    for idx in good_indices[1:min(3, length(good_indices))]
        p = refined_result[:refined_points][idx]
        v = refined_result[:refined_values][idx]
        println("  params=$p, value=$(round(v, digits=8))")
    end
end

#===============================================================================
CUSTOMIZATION OPTIONS
===============================================================================#

# Pipeline wrapper benefits:
# ✓ Single call for both stages
# ✓ Returns both raw and refined results
# ✓ Less boilerplate than manual approach
# ✓ Still configurable via ExperimentParams
#
# Configuration options:
#
# config (ExperimentParams or NamedTuple):
#   - domain_size: Search domain size
#   - GN: Grid resolution
#   - degree_range: Polynomial degrees
#   - basis: :chebyshev or :legendre
#   - max_time: Time limit per degree
#
# refinement_config (optional):
#   - max_time_per_point: BFGS time limit
#   - f_tol, x_tol: Convergence tolerances
#   - max_iterations: BFGS iterations
#   - verbose: Detailed output
#
# refine::Bool (optional):
#   - Set to false to skip Stage 2
#   - Returns only raw_result
#
# When to use each template:
# ┌─────────────────────────────────────────────────────────────┐
# │ Simple    → Quick tests, don't care about details          │
# │ Pipeline  → Production workflows, want both results         │ ← YOU ARE HERE
# │ Advanced  → Maximum control, HPC, custom analysis           │
# └─────────────────────────────────────────────────────────────┘
