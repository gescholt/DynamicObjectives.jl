#!/usr/bin/env julia
"""
Batch testing of Dynamic_objectives models with globtim optimizer

Tests all EASY and MEDIUM models to validate integration and measure performance.

Usage:
    julia --project=. examples/globtim_batch_test.jl [easy|medium|all]

Arguments:
    easy   - Test only EASY models (4 models, 2D)
    medium - Test only MEDIUM models (5 models, 3-4D)
    all    - Test both EASY and MEDIUM (default)

Output:
    - Individual results in test_results/batch_globtim/<model_name>/
    - Summary table printed to console
    - Summary CSV: test_results/batch_globtim/batch_summary.csv
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dynamic_objectives
using Globtim
using LinearAlgebra
using Printf
using Dates
using CSV
using DataFrames

# Load test configurations
include("test_configs.jl")

# Parse command line arguments
test_phase = length(ARGS) > 0 ? ARGS[1] : "all"

println("="^80)
println("Batch Testing: Dynamic_objectives + globtim")
println("="^80)
println("  Phase: $test_phase")
println("  Start time: $(now())")
println("="^80)

# ==============================================================================
# Configuration
# ==============================================================================

# globtim settings (standard configuration)
const GLOBTIM_CONFIG = (
    GN = 8,
    degree_range = 4:8,
    max_time = 3600,
    basis = :chebyshev
)

# Output directory
output_base = mkpath(joinpath(@__DIR__, "..", "test_results", "batch_globtim"))

# Success criteria
const RECOVERY_TOLERANCE_STRICT = 0.01
const RECOVERY_TOLERANCE_LOOSE = 0.1
const OBJECTIVE_TOLERANCE = 1e-6

# ==============================================================================
# Helper Functions
# ==============================================================================

"""Run a single model test and return summary statistics"""
function test_single_model(config, output_dir)
    println("\nTesting: $(config.name)")
    println("  Description: $(config.description)")
    println("  Parameters: $(length(config.p_true))D")

    start_time = time()

    try
        # Define model
        model, params, states, outputs = config.model_fn()

        # Create objective
        objective = create_globtim_objective(
            model, outputs, config.ic, config.p_true,
            config.time_interval, config.numpoints;
            distance = config.distance,
            eval_timeout = config.timeout
        )

        # Run optimization
        result = run_globtim_optimization(
            objective,
            config.bounds,
            config.p_true;
            GN = GLOBTIM_CONFIG.GN,
            degree_range = GLOBTIM_CONFIG.degree_range,
            output_dir = joinpath(output_dir, config.name),
            model_name = config.name,
            max_time = GLOBTIM_CONFIG.max_time,
            basis = GLOBTIM_CONFIG.basis
        )

        elapsed = time() - start_time

        # Extract metrics
        success = result[:success_rate] > 0.0
        recovery_error = result[:recovery_error]
        best_objective = result[:best_objective]
        n_critical_points = result[:total_critical_points]

        # Classify result
        if !success
            status = "FAILED"
        elseif recovery_error < RECOVERY_TOLERANCE_STRICT && best_objective < OBJECTIVE_TOLERANCE
            status = "EXCELLENT"
        elseif recovery_error < RECOVERY_TOLERANCE_LOOSE
            status = "GOOD"
        else
            status = "POOR"
        end

        @printf "  ✓ Status: %s\n" status
        @printf "    Time: %.2fs\n" elapsed
        @printf "    Recovery error: %.6e\n" recovery_error
        @printf "    Objective: %.6e\n" best_objective
        @printf "    Critical points: %d\n" n_critical_points

        return (
            name = config.name,
            dimension = length(config.p_true),
            status = status,
            success = success,
            recovery_error = recovery_error,
            objective_value = best_objective,
            time_seconds = elapsed,
            critical_points = n_critical_points,
            success_rate = result[:success_rate],
            best_degree = result[:best_degree],
            error_message = nothing
        )

    catch e
        elapsed = time() - start_time

        @printf "  ❌ EXCEPTION: %s\n" e
        println("    $(sprint(showerror, e))")

        return (
            name = config.name,
            dimension = length(config.p_true),
            status = "EXCEPTION",
            success = false,
            recovery_error = Inf,
            objective_value = Inf,
            time_seconds = elapsed,
            critical_points = 0,
            success_rate = 0.0,
            best_degree = nothing,
            error_message = string(e)
        )
    end
end

# ==============================================================================
# Run Tests
# ==============================================================================

results = []

# Select which configurations to test
test_configs = if test_phase == "easy"
    println("\nRunning EASY models only (4 models, 2D)...")
    easy_configs
elseif test_phase == "medium"
    println("\nRunning MEDIUM models only (5 models, 3-4D)...")
    medium_configs
else
    println("\nRunning EASY + MEDIUM models (9 total)...")
    vcat(easy_configs, medium_configs)
end

println("  Total models: $(length(test_configs))")
println("  Configuration: GN=$(GLOBTIM_CONFIG.GN), degrees=$(GLOBTIM_CONFIG.degree_range)")

# Run each test
for (i, config) in enumerate(test_configs)
    println("\n" * "="^80)
    println("[$i/$(length(test_configs))]")

    result = test_single_model(config, output_base)
    push!(results, result)
end

# ==============================================================================
# Summary Statistics
# ==============================================================================

println("\n" * "="^80)
println("BATCH TEST SUMMARY")
println("="^80)

# Convert to DataFrame for easy analysis
df = DataFrame(results)

# Overall statistics
n_total = nrow(df)
n_success = sum(df.success)
n_excellent = sum(df.status .== "EXCELLENT")
n_good = sum(df.status .== "GOOD")
n_poor = sum(df.status .== "POOR")
n_failed = sum(.!df.success)

println("\nOverall Results:")
@printf "  Total models: %d\n" n_total
@printf "  Successful: %d (%.1f%%)\n" n_success (n_success / n_total * 100)
@printf "  Failed: %d (%.1f%%)\n" n_failed (n_failed / n_total * 100)
println("\nQuality Breakdown:")
@printf "  EXCELLENT: %d (recovery <0.01, objective <1e-6)\n" n_excellent
@printf "  GOOD: %d (recovery <0.1)\n" n_good
@printf "  POOR: %d (recovery >0.1)\n" n_poor

# Statistics by dimension
println("\nBy Dimension:")
for dim in sort(unique(df.dimension))
    df_dim = filter(row -> row.dimension == dim, df)
    n_dim = nrow(df_dim)
    n_success_dim = sum(df_dim.success)
    @printf "  %dD: %d/%d success (%.1f%%)\n" dim n_success_dim n_dim (n_success_dim / n_dim * 100)
end

# Timing statistics
successful_results = filter(row -> row.success, df)
if nrow(successful_results) > 0
    println("\nTiming (successful runs only):")
    @printf "  Mean: %.2f seconds\n" mean(successful_results.time_seconds)
    @printf "  Median: %.2f seconds\n" median(successful_results.time_seconds)
    @printf "  Min: %.2f seconds\n" minimum(successful_results.time_seconds)
    @printf "  Max: %.2f seconds\n" maximum(successful_results.time_seconds)
end

# Detailed table
println("\n" * "="^80)
println("DETAILED RESULTS")
println("="^80)
println("\n", df)

# ==============================================================================
# Save Results
# ==============================================================================

summary_file = joinpath(output_base, "batch_summary.csv")
CSV.write(summary_file, df)

println("\n" * "="^80)
println("Results saved to:")
println("  Summary CSV: $summary_file")
println("  Individual results: $output_base/<model_name>/")
println("="^80)

# ==============================================================================
# Success Criteria Evaluation
# ==============================================================================

println("\n" * "="^80)
println("SUCCESS CRITERIA EVALUATION")
println("="^80)

# EASY models: >75% success
easy_results = filter(row -> row.dimension == 2, df)
if nrow(easy_results) > 0
    easy_success_rate = sum(easy_results.success) / nrow(easy_results)
    if easy_success_rate >= 0.75
        @printf "✓ EASY models: %.1f%% success (target: >75%%)\n" (easy_success_rate * 100)
    else
        @printf "✗ EASY models: %.1f%% success (target: >75%%)\n" (easy_success_rate * 100)
    end
end

# MEDIUM models: >50% success
medium_results = filter(row -> row.dimension >= 3 && row.dimension <= 4, df)
if nrow(medium_results) > 0
    medium_success_rate = sum(medium_results.success) / nrow(medium_results)
    if medium_success_rate >= 0.50
        @printf "✓ MEDIUM models: %.1f%% success (target: >50%%)\n" (medium_success_rate * 100)
    else
        @printf "✗ MEDIUM models: %.1f%% success (target: >50%%)\n" (medium_success_rate * 100)
    end
end

println("\nEnd time: $(now())")
println("="^80)
