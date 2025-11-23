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

# Include ExperimentCLI module to access ExperimentParams
# Note: Must load Globtim first since ExperimentCLI uses ConstructionBase
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(dirname(@__DIR__), "..", "globtimcore")
    include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
end
using .ExperimentCLI

println("Testing simple 2-stage workflow")
println("="^80)

# Create test objective
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

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

# Create output directory
output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "simple_workflow"))

# Stage 1: Get raw critical points
println("\nStage 1: run_standard_experiment...")
result = run_standard_experiment(
    objective_function = objective,
    problem_params = nothing,
    domain_bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict{String, Any}("experiment_type" => "lv2d_simple_test"),
    true_params = p_true
)

println("✓ Found critical points across $(result[:degrees_processed]) degrees")
println("  Total critical points: $(result[:total_critical_points])")
println("  Best raw value: $(minimum([r.best_objective for r in result[:degree_results]]))")

# Configure refinement
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,
    verbose = false
)

# Stage 2: Refine in post-processing
println("\nStage 2: refine_experiment_results...")
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
println("✓ Refined $(refined[:n_converged])/$(refined[:n_raw]) points")
println("  Best refined value: $(refined[:best_refined_value])")
println("  Improvement: $(round(refined[:mean_improvement], digits=2))x")

# Verify
best_params = refined[:refined_points][refined[:best_refined_idx]]
recovery_error = norm(best_params .- p_true) / norm(p_true)
println("\nRecovery error: $(round(100*recovery_error, digits=2))%")

if recovery_error < 0.05
    println("✅ SUCCESS!")
else
    println("⚠️  Needs improvement")
end
