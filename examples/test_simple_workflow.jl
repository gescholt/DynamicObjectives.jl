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

println("="^80)
println("Testing Simple 2-Stage Workflow")
println("="^80)
println()

# ==============================================================================
# Step 1: Model Setup and Verification
# ==============================================================================

println("Step 1: Model Setup")
println("-"^80)

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Display model summary
display_model_summary(model)
println()

# Display parameters
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)
println()

# Create objective function
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

# Verify objective works at true parameters
obj_at_true = objective(p_true)
println("✓ Objective at true parameters: $(round(obj_at_true, digits=8))")
@assert obj_at_true < 1e-6 "Objective should be near zero at true parameters"
println()

# ==============================================================================
# Step 2: Configure Experiment
# ==============================================================================

println("Step 2: Experiment Configuration")
println("-"^80)

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

println("Configuration:")
println("  Grid size (GN): $(config.GN) ($(config.GN^2) points for 2D)")
println("  Degree range: $(config.degree_range)")
println("  Domain size: $(config.domain_size)")
println("  Basis: $(config.basis)")
println()

# Create output directory
output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "simple_workflow"))

# ==============================================================================
# Step 3: Stage 1 - Find Raw Critical Points
# ==============================================================================

println("Step 3: Stage 1 - Find Raw Critical Points (globtimcore)")
println("-"^80)

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

println()
println("✓ Stage 1 Complete")
println("  Degrees processed: $(result[:degrees_processed])")
println("  Total critical points: $(result[:total_critical_points])")
println("  Best raw value: $(round(minimum([r.best_objective for r in result[:degree_results]]), digits=6))")
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
    verbose = false
)

println("Refinement configuration:")
println("  Max time per point: $(refinement_config.max_time_per_point)s")
println("  Method: BFGS local optimization")
println()

# Run refinement
t_start = time()
refined = refine_experiment_results(
    result[:output_dir],
    objective,
    refinement_config
)
stage2_time = time() - t_start

println()
println("✓ Stage 2 Complete")
println("  Raw points: $(refined[:n_raw])")
println("  Converged: $(refined[:n_converged])")
println("  Success rate: $(round(100*refined[:n_converged]/refined[:n_raw], digits=1))%")
println("  Mean improvement: $(round(refined[:mean_improvement], digits=2))x")
println("  Best refined value: $(round(refined[:best_refined_value], digits=8))")
println("  Time elapsed: $(round(stage2_time, digits=1))s")
println()

# ==============================================================================
# Step 5: Verification
# ==============================================================================

println("Step 5: Parameter Recovery Verification")
println("-"^80)

# Get best refined parameters
best_params = refined[:refined_points][refined[:best_refined_idx]]
recovery_error = norm(best_params .- p_true) / norm(p_true)

# Display comparison
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true, best_params),
    labels=["True", "Recovered"]
)
println()

println("Recovery metrics:")
println("  Relative error: $(round(100*recovery_error, digits=2))%")
println("  Objective value: $(round(refined[:best_refined_value], digits=8))")
println("  Total time: $(round(stage1_time + stage2_time, digits=1))s")
println()

# Final verdict
if recovery_error < 0.01
    println("✅ SUCCESS: Excellent recovery (< 1% error)")
elseif recovery_error < 0.05
    println("✅ SUCCESS: Good recovery (< 5% error)")
else
    println("⚠️  Needs improvement (> 5% error)")
end
println()

println("="^80)
println("2-Stage Workflow Complete!")
println("="^80)
