#!/usr/bin/env julia
"""
Analyze globtim batch test results using globtimpostprocessing

Loads results from batch runs and generates comprehensive analysis:
- Statistical summaries
- Success rate trends
- Performance comparisons
- Quality diagnostics

Usage:
    julia --project=. examples/analyze_globtim_results.jl [results_dir]

Arguments:
    results_dir - Directory containing batch test results (default: test_results/batch_globtim)

Dependencies:
    - GlobtimPostProcessing package must be available in environment
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Core packages
using DataFrames
using Statistics
using CSV
using Printf
using Dates

using GlobtimPostProcessing

# ==============================================================================
# Configuration
# ==============================================================================

# Get results directory from command line or use default
results_dir = if length(ARGS) > 0
    ARGS[1]
else
    joinpath(@__DIR__, "..", "test_results", "batch_globtim")
end

println("="^80)
println("globtim Results Analysis")
println("="^80)
println("  Results directory: $results_dir")
println("  Analysis time: $(now())")
println("="^80)

# ==============================================================================
# Load Batch Summary
# ==============================================================================

summary_file = joinpath(results_dir, "batch_summary.csv")

if !isfile(summary_file)
    error(
        "Summary file not found: $summary_file\nExpected batch_summary.csv in results directory.",
    )
end

println("\n[1/5] Loading batch summary...")
df = CSV.read(summary_file, DataFrame)
println("  Loaded $(nrow(df)) model results")

# ==============================================================================
# Basic Statistics
# ==============================================================================

println("\n[2/5] Computing basic statistics...")

# Overall metrics
n_total = nrow(df)
n_success = sum(df.success)
success_rate = n_success / n_total

# Quality counts
n_excellent = sum(df.status .== "EXCELLENT")
n_good = sum(df.status .== "GOOD")
n_poor = sum(df.status .== "POOR")
n_failed = n_total - n_success

println("\n--- Overall Performance ---")
@printf "  Models tested: %d\n" n_total
@printf "  Success rate: %.1f%% (%d/%d)\n" (success_rate * 100) n_success n_total
println("\n  Quality distribution:")
@printf "    EXCELLENT: %d (%.1f%%) - recovery <0.01, objective <1e-6\n" n_excellent (
    n_excellent / n_total * 100
)
@printf "    GOOD: %d (%.1f%%) - recovery <0.1\n" n_good (n_good / n_total * 100)
@printf "    POOR: %d (%.1f%%) - recovery >0.1\n" n_poor (n_poor / n_total * 100)
@printf "    FAILED: %d (%.1f%%) - no solution found\n" n_failed (n_failed / n_total * 100)

# ==============================================================================
# Analysis by Dimension
# ==============================================================================

println("\n[3/5] Analysis by dimension...")

println("\n--- Performance by Parameter Dimension ---")
for dim in sort(unique(df.dimension))
    df_dim = filter(row -> row.dimension == dim, df)
    n_dim = nrow(df_dim)
    n_success_dim = sum(df_dim.success)
    success_rate_dim = n_success_dim / n_dim

    # Timing stats for successful runs
    df_dim_success = filter(row -> row.success, df_dim)
    time_stats = if nrow(df_dim_success) > 0
        @sprintf(
            "mean=%.2fs, median=%.2fs",
            mean(df_dim_success.time_seconds),
            median(df_dim_success.time_seconds)
        )
    else
        "N/A"
    end

    # Recovery error stats
    recovery_stats = if nrow(df_dim_success) > 0
        finite_errors = filter(!isinf, df_dim_success.recovery_error)
        if length(finite_errors) > 0
            @sprintf("mean=%.2e, median=%.2e", mean(finite_errors), median(finite_errors))
        else
            "all infinite"
        end
    else
        "N/A"
    end

    @printf "\n  %dD Models:\n" dim
    @printf "    Count: %d\n" n_dim
    @printf "    Success: %d/%d (%.1f%%)\n" n_success_dim n_dim (success_rate_dim * 100)
    @printf "    Time: %s\n" time_stats
    @printf "    Recovery error: %s\n" recovery_stats
end

# ==============================================================================
# Critical Points Analysis
# ==============================================================================

println("\n[4/5] Critical points analysis...")

successful_results = filter(row -> row.success, df)

if nrow(successful_results) > 0
    println("\n--- Critical Points Found ---")
    @printf "  Total critical points: %d\n" sum(successful_results.critical_points)
    @printf "  Mean per model: %.2f\n" mean(successful_results.critical_points)
    @printf "  Median per model: %.0f\n" median(successful_results.critical_points)
    @printf "  Range: %d - %d\n" minimum(successful_results.critical_points) maximum(
        successful_results.critical_points,
    )

    # Degree analysis
    if "best_degree" in names(df)
        degrees = filter(!isnothing, successful_results.best_degree)
        if length(degrees) > 0
            println("\n--- Best Polynomial Degree ---")
            degree_counts =
                combine(groupby(DataFrame(degree = degrees), :degree), nrow => :count)
            sort!(degree_counts, :count, rev = true)
            for row in eachrow(degree_counts)
                @printf "    Degree %d: %d times (%.1f%%)\n" row.degree row.count (
                    row.count / length(degrees) * 100
                )
            end
        end
    end
else
    println("  No successful results to analyze")
end

# ==============================================================================
# Model-Specific Analysis
# ==============================================================================

println("\n[5/5] Detailed model breakdown...")

println("\n--- Individual Model Performance ---")
println("")
@printf "%-25s %4s %10s %12s %12s %8s\n" "Model" "Dim" "Status" "Recovery Err" "Objective" "Time (s)"
println("-"^80)

for row in eachrow(df)
    recovery_str = isinf(row.recovery_error) ? "Inf" : @sprintf("%.2e", row.recovery_error)
    objective_str =
        isinf(row.objective_value) ? "Inf" : @sprintf("%.2e", row.objective_value)
    status_icon = row.success ? "✓" : "✗"

    @printf "%-25s %4d %3s %-6s %12s %12s %8.2f\n"
    row.name
    row.dimension
    status_icon
    row.status
    recovery_str
    objective_str
    row.time_seconds
end

# ==============================================================================
# Advanced Analysis (if GlobtimPostProcessing available)
# ==============================================================================

if @isdefined(GlobtimPostProcessing)
    println("\n" * "="^80)
    println("ADVANCED ANALYSIS (using GlobtimPostProcessing)")
    println("="^80)

    println("\n[Loading individual experiment results...]")

    # Try to load detailed results for each model
    for row in eachrow(df)
        if !row.success
            continue
        end

        model_dir = joinpath(results_dir, row.name)
        if !isdir(model_dir)
            continue
        end

        try
            # Load experiment results
            println("\n  Analyzing: $(row.name)")

            # Look for results_summary.json
            summary_json = joinpath(model_dir, "results_summary.json")
            if isfile(summary_json)
                # Could load and analyze JSON here
                println("    ✓ Found results_summary.json")
            end

            # Look for critical points CSV
            cp_files =
                filter(f -> startswith(f, "critical_points_deg_"), readdir(model_dir))
            if length(cp_files) > 0
                println("    ✓ Found $(length(cp_files)) critical point files")

                # Could analyze distribution of critical points here
                # - Spatial distribution
                # - Objective values
                # - Distance to true parameters
            end

        catch e
            @warn "  Failed to analyze $(row.name): $e"
        end
    end
else
    println("\n" * "="^80)
    println("Note: Install GlobtimPostProcessing for advanced analysis features:")
    println("  - Critical point clustering")
    println("  - Convergence diagnostics")
    println("  - Landscape fidelity metrics")
    println("  - Comparative visualizations")
    println("="^80)
end

# ==============================================================================
# Success Criteria Check
# ==============================================================================

println("\n" * "="^80)
println("SUCCESS CRITERIA ASSESSMENT")
println("="^80)

all_criteria_met = true

# Criterion 1: EASY models >75% success
easy_results = filter(row -> row.dimension == 2, df)
if nrow(easy_results) > 0
    easy_success_rate = sum(easy_results.success) / nrow(easy_results)
    if easy_success_rate >= 0.75
        @printf "\n✓ PASS: EASY models (2D)\n"
        @printf "  Success rate: %.1f%% (target: ≥75%%)\n" (easy_success_rate * 100)
    else
        @printf "\n✗ FAIL: EASY models (2D)\n"
        @printf "  Success rate: %.1f%% (target: ≥75%%)\n" (easy_success_rate * 100)
        all_criteria_met = false
    end
end

# Criterion 2: MEDIUM models >50% success
medium_results = filter(row -> row.dimension >= 3 && row.dimension <= 4, df)
if nrow(medium_results) > 0
    medium_success_rate = sum(medium_results.success) / nrow(medium_results)
    if medium_success_rate >= 0.50
        @printf "\n✓ PASS: MEDIUM models (3-4D)\n"
        @printf "  Success rate: %.1f%% (target: ≥50%%)\n" (medium_success_rate * 100)
    else
        @printf "\n✗ FAIL: MEDIUM models (3-4D)\n"
        @printf "  Success rate: %.1f%% (target: ≥50%%)\n" (medium_success_rate * 100)
        all_criteria_met = false
    end
end

# Criterion 3: No crashes (all runs completed)
if "error_message" in names(df)
    n_exceptions = sum(.!ismissing.(df.error_message))
    if n_exceptions == 0
        @printf "\n✓ PASS: Stability\n"
        @printf "  No crashes or exceptions\n"
    else
        @printf "\n⚠ WARNING: Stability\n"
        @printf "  %d models threw exceptions\n" n_exceptions
    end
end

# Criterion 4: Parameter recovery accuracy (for EXCELLENT models)
if n_excellent > 0
    excellent_models = filter(row -> row.status == "EXCELLENT", df)
    max_recovery = maximum(excellent_models.recovery_error)
    @printf "\n✓ PASS: Recovery accuracy\n"
    @printf "  %d models with recovery <0.01 (max: %.2e)\n" n_excellent max_recovery
end

# Final assessment
println("\n" * "="^80)
if all_criteria_met
    println("✅ OVERALL: All success criteria met!")
else
    println("❌ OVERALL: Some criteria not met - see details above")
end
println("="^80)

println("\nAnalysis complete: $(now())")
