#!/usr/bin/env julia
"""
Integration Template: Advanced Manual 2-Stage Pattern

This template demonstrates FULL CONTROL over the 2-stage pipeline:
- Stage 1 (Globtim): Polynomial approximation + HomotopyContinuation
- Stage 2 (globtimpostprocessing): Local BFGS refinement

USE WHEN:
- You need fine control over experiment parameters
- You're running on HPC or batch systems
- You want to inspect intermediate results
- You need custom experiment configs

PATTERN:
- Manual include of ExperimentCLI module
- Explicit ExperimentParams creation
- Direct calls to run_standard_experiment() and refine_experiment_results()
- Full access to all configuration options

ARCHITECTURE:
┌──────────────────────────────────────────────────────────────┐
│ Stage 1: Globtim                                              │
│   Input:  objective, bounds, config                          │
│   Output: raw_critical_points_deg_X.csv                      │
│   Method: Polynomial approx + HomotopyContinuation           │
└──────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│ Stage 2: globtimpostprocessing                               │
│   Input:  raw_critical_points_deg_X.csv, objective           │
│   Output: refined_critical_points_deg_X.csv                  │
│   Method: Local BFGS optimization                            │
└──────────────────────────────────────────────────────────────┘
"""

using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))  # Activate DynamicObjectives environment

using DynamicObjectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing:
    GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra

#===============================================================================
STEP 0: Include ExperimentCLI module for config creation
===============================================================================#

# The ExperimentCLI module provides ExperimentParams struct
# It's not exported from Globtim, so we include it manually
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
STEP 2: Configure Stage 1 (Globtim)
===============================================================================#

# Create ExperimentParams with full control over all parameters
config = ExperimentParams(
    # Domain configuration
    domain_size = 1.5,              # Centered on middle of bounds
    GN = 50,                        # Grid size (50 points per dim)

    # Polynomial approximation
    degree_range = 4:10,            # Try degrees 4 through 10
    basis = :chebyshev,             # Chebyshev or Legendre

    # Time limits
    max_time = 3600.0,              # 1 hour max per degree

    # Optimization tolerances (for grid-based initial guesses)
    optim_f_tol = 1e-6,
    optim_x_tol = 1e-6,
    max_iterations = 300,

    # Feature toggles (Phase 2: these are disabled in Globtim)
    enable_gradient_computation = false,    # No gradients in Stage 1
    enable_hessian_computation = false,     # No Hessians in Stage 1
    enable_bfgs_refinement = false,          # No BFGS in Stage 1
)

# Create output directory
output_dir = mkpath(joinpath(@__DIR__, "..", "..", "test_results", "template_advanced"))

println("="^80)
println("Advanced 2-Stage Integration Template")
println("="^80)
println("Configuration:")
println("  Domain size:    $(config.domain_size)")
println("  Grid size:      $(config.GN)")
println("  Degree range:   $(config.degree_range)")
println("  Basis:          $(config.basis)")
println("  Output dir:     $output_dir")

#===============================================================================
STEP 3: Stage 1 - Find raw critical points
===============================================================================#

println("\n" * "="^80)
println("Stage 1: Finding raw critical points with Globtim")
println("="^80)

stage1_start = time()

result = run_standard_experiment(
    objective_function = objective,
    objective_name = "lv2d_advanced_template",
    bounds = bounds,
    experiment_config = config,
    output_dir = output_dir,
    metadata = Dict(
        "experiment_type" => "lv2d_advanced_template",
        "description" => "Advanced 2-stage integration example",
    ),
    true_params = p_true,
)

stage1_time = time() - stage1_start

println("\nStage 1 Results:")
println("  Degrees processed:   $(result[:degrees_processed])")
println("  Total critical pts:  $(result[:total_critical_points])")
println("  Success rate:        $(round(100*result[:success_rate], digits=1))%")
println("  Computation time:    $(round(stage1_time, digits=2))s")

# Extract best raw value across all degrees
best_raw_value = minimum([r.best_objective for r in result[:degree_results]])
println("  Best raw value:      $(round(best_raw_value, digits=6))")

# Inspect per-degree results
println("\n  Per-degree breakdown:")
for deg_result in result[:degree_results]
    status_symbol = deg_result.status == "success" ? "✓" : "✗"
    println(
        "    $status_symbol Degree $(deg_result.degree): $(deg_result.n_critical_points) points, best=$(round(deg_result.best_objective, digits=6))",
    )
end

#===============================================================================
STEP 4: Configure Stage 2 (globtimpostprocessing)
===============================================================================#

println("\n" * "="^80)
println("Stage 2: Refining critical points with globtimpostprocessing")
println("="^80)

# Create refinement config for ODE problems
refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,      # 30 seconds per point max
    f_tol = 1e-10,                  # Tighter tolerance for refinement
    x_tol = 1e-10,
    max_iterations = 1000,
    verbose = false,                  # Set to true for detailed output
)

println("Refinement configuration:")
println("  Method:           BFGS")
println("  Max time/point:   $(refinement_config[:max_time_per_point])s")
println("  f_tol:            $(refinement_config[:f_tol])")
println("  x_tol:            $(refinement_config[:x_tol])")

#===============================================================================
STEP 5: Stage 2 - Refine critical points
===============================================================================#

stage2_start = time()

refined = refine_experiment_results(
    result[:output_dir],      # Where to find raw_critical_points_deg_X.csv
    objective,                # Same objective function
    refinement_config,
)

stage2_time = time() - stage2_start

println("\nStage 2 Results:")
println("  Raw points:          $(refined.n_raw)")
println("  Converged points:    $(refined.n_converged)")
println("  Convergence rate:    $(round(100*refined.n_converged/refined.n_raw, digits=1))%")
println("  Best refined value:  $(round(refined.best_refined_value, digits=8))")
println("  Mean improvement:    $(round(refined.mean_improvement, digits=2))x")
println("  Computation time:    $(round(stage2_time, digits=2))s")

#===============================================================================
STEP 6: Analyze results
===============================================================================#

println("\n" * "="^80)
println("Final Analysis")
println("="^80)

# Check if any points converged
if refined.n_converged > 0
    best_params = refined.refined_points[refined.best_refined_idx]
    recovery_error = norm(best_params .- p_true) / norm(p_true)

    println("True parameters:        $p_true")
    println("Best recovered params:  $best_params")
    println("Recovery error:         $(round(100*recovery_error, digits=3))%")
    println()
    println("Total pipeline time:    $(round(stage1_time + stage2_time, digits=2))s")
    println(
        "  Stage 1 (raw):        $(round(stage1_time, digits=2))s ($(round(100*stage1_time/(stage1_time+stage2_time), digits=1))%)",
    )
    println(
        "  Stage 2 (refine):     $(round(stage2_time, digits=2))s ($(round(100*stage2_time/(stage1_time+stage2_time), digits=1))%)",
    )

    if recovery_error < 0.01
        println("\n✅ Excellent recovery (< 1% error)")
    elseif recovery_error < 0.05
        println("\n✅ Good recovery (< 5% error)")
    elseif recovery_error < 0.10
        println("\n✓ Acceptable recovery (< 10% error)")
    else
        println("\n⚠️  Poor recovery - consider:")
        println("   - Increasing domain_size")
        println("   - Increasing GN (grid size)")
        println("   - Expanding degree_range")
        println("   - Checking if objective is well-behaved")
    end
else
    # No points converged - refinement needs improvement
    # TODO: Investigate postprocessing step - understand optimization methods,
    #       convergence criteria, and alternative refinement approaches
    println(
        "⚠️  Refinement did not converge for any of the $(refined.n_raw) raw critical points",
    )
    println()
    println("Total pipeline time:    $(round(stage1_time + stage2_time, digits=2))s")
    println("  Stage 1 (raw):        $(round(stage1_time, digits=2))s")
    println("  Stage 2 (refine):     $(round(stage2_time, digits=2))s")
end

#===============================================================================
STEP 7: Access saved files
===============================================================================#

println("\n" * "="^80)
println("Saved Files")
println("="^80)
println("Output directory: $output_dir")
println()
println("Raw critical points:")
for deg_result in result[:degree_results]
    if deg_result.status == "success"
        csv_path = joinpath(output_dir, "critical_points_raw_deg_$(deg_result.degree).csv")
        if isfile(csv_path)
            println("  ✓ $csv_path")
        end
    end
end

println("\nRefined critical points:")
for deg_result in result[:degree_results]
    if deg_result.status == "success"
        csv_path =
            joinpath(output_dir, "critical_points_refined_deg_$(deg_result.degree).csv")
        if isfile(csv_path)
            println("  ✓ $csv_path")
        end
    end
end

summary_path = joinpath(output_dir, "experiment_summary.json")
if isfile(summary_path)
    println("\nExperiment summary:")
    println("  ✓ $summary_path")
end

#===============================================================================
CUSTOMIZATION NOTES
===============================================================================#

# This template shows ALL the knobs you can turn:
#
# ExperimentParams fields:
#   - domain_size: How far from center of bounds to search
#   - GN: Grid size (points per dimension)
#   - degree_range: Polynomial degrees to try
#   - basis: :chebyshev or :legendre
#   - max_time: Time limit per degree
#   - optim_f_tol, optim_x_tol, max_iterations: Grid refinement params
#   - enable_*: Feature toggles (Phase 2: set to false)
#
# run_standard_experiment arguments:
#   - objective_function: Your 1-arg function p -> error
#   - bounds: [(min, max), ...] for each dimension
#   - experiment_config: ExperimentParams instance
#   - output_dir: Where to save results
#   - metadata: Dict with experiment info
#   - true_params: For computing recovery error
#
# ode_refinement_config arguments:
#   - max_time_per_point: Time limit for BFGS per point
#   - f_tol, x_tol: Convergence tolerances
#   - max_iterations: Max BFGS iterations
#   - verbose: Print detailed refinement info
#
# For HPC batch jobs:
#   1. Create config from CLI args: config = parse_experiment_args(ARGS)
#   2. Run Stage 1 on compute node
#   3. Run Stage 2 on compute node or locally
#   4. Analyze results locally with globtimpostprocessing
