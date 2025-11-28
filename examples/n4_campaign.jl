#!/usr/bin/env julia
"""
Testing Campaign for n=4 Dimension Models

This script runs a comprehensive testing campaign for all 4-parameter models
in the Dynamic_objectives benchmark suite.

Features:
- Timing information saved for each model and stage
- Results saved in Markdown format (preserving rich formatting)
- CSV summary for data analysis
- Non-interactive output suitable for batch jobs
- Configurable via environment variables

Usage:
    julia --project=. examples/n4_hpc_campaign.jl

    # Or use the runner script:
    ./scripts/run_n4_campaign.sh --gn 10 --max-time 1800

Environment Variables:
    CAMPAIGN_OUTPUT_DIR   - Output directory (default: test_results/n4_campaign)
    CAMPAIGN_GN           - Grid parameter (default: 8)
    CAMPAIGN_MAX_TIME     - Max time per model in seconds (default: 3600)
    CAMPAIGN_MODELS       - Comma-separated model names (default: all n=4 models)

Output Files:
    - campaign_report.md       - Full Markdown report with timing and results
    - campaign_summary.csv     - CSV summary for data analysis
    - timing_breakdown.csv     - Detailed timing information
    - <model_name>/            - Per-model result directories
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
using Statistics

# Load test configurations
include("test_configs.jl")

# ==============================================================================
# Configuration from Environment Variables
# ==============================================================================

const OUTPUT_DIR = get(ENV, "CAMPAIGN_OUTPUT_DIR",
    joinpath(@__DIR__, "..", "test_results", "n4_campaign"))
const GN = parse(Int, get(ENV, "CAMPAIGN_GN", "8"))
const MAX_TIME = parse(Int, get(ENV, "CAMPAIGN_MAX_TIME", "3600"))
const DEGREE_RANGE = 4:8
const BASIS = :chebyshev

# Model selection - default to all n=4 models
const DEFAULT_N4_MODELS = [
    "DAISY_Ex3_with_input",
    "DAISY_Ex3_no_input",
    "LV_4D_Constrained"
]

function get_selected_models()
    model_env = get(ENV, "CAMPAIGN_MODELS", "")
    if isempty(model_env)
        return DEFAULT_N4_MODELS
    else
        return split(model_env, ",") .|> strip .|> String
    end
end

# Success criteria
const RECOVERY_TOLERANCE_STRICT = 0.01
const RECOVERY_TOLERANCE_LOOSE = 0.1
const OBJECTIVE_TOLERANCE = 1e-6

# ==============================================================================
# Markdown Output Functions
# ==============================================================================

"""Generate a Markdown table from a DataFrame"""
function df_to_markdown(df::DataFrame; precision::Int=6)
    io = IOBuffer()

    # Header
    cols = names(df)
    println(io, "| ", join(cols, " | "), " |")
    println(io, "| ", join(fill("---", length(cols)), " | "), " |")

    # Rows
    for row in eachrow(df)
        row_strs = String[]
        for col in cols
            val = row[col]
            if val isa AbstractFloat
                if abs(val) < 1e-3 || abs(val) > 1e4
                    push!(row_strs, @sprintf("%.4e", val))
                else
                    push!(row_strs, string(round(val, digits=precision)))
                end
            elseif val === nothing
                push!(row_strs, "-")
            else
                push!(row_strs, string(val))
            end
        end
        println(io, "| ", join(row_strs, " | "), " |")
    end

    return String(take!(io))
end

"""Generate timing table in Markdown format"""
function timing_to_markdown(timings::Vector{NamedTuple})
    io = IOBuffer()

    println(io, "| Model | Stage 1 (s) | Stage 2 (s) | Total (s) | Status |")
    println(io, "| --- | ---: | ---: | ---: | --- |")

    for t in timings
        s1 = isnothing(t.stage1_time) ? "-" : @sprintf("%.2f", t.stage1_time)
        s2 = isnothing(t.stage2_time) ? "-" : @sprintf("%.2f", t.stage2_time)
        total = @sprintf("%.2f", t.total_time)
        println(io, "| $(t.model) | $s1 | $s2 | $total | $(t.status) |")
    end

    return String(take!(io))
end

# ==============================================================================
# Test Execution Functions
# ==============================================================================

"""Run a single model test and return detailed results with timing"""
function test_model_with_timing(config, output_dir)
    model_output_dir = mkpath(joinpath(output_dir, config.name))

    timing = (
        model = config.name,
        stage1_time = nothing,
        stage2_time = nothing,
        total_time = 0.0,
        status = "PENDING"
    )

    total_start = time()

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

        # Stage 1: Core optimization
        stage1_start = time()
        result = run_globtim_optimization(
            objective,
            config.bounds,
            config.p_true;
            GN = GN,
            degree_range = DEGREE_RANGE,
            output_dir = model_output_dir,
            model_name = config.name,
            max_time = MAX_TIME,
            basis = BASIS
        )
        stage1_time = time() - stage1_start

        # Stage 2 timing is included in the result if available
        stage2_time = get(result, :stage2_time, nothing)

        total_time = time() - total_start

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

        timing = (
            model = config.name,
            stage1_time = stage1_time,
            stage2_time = stage2_time,
            total_time = total_time,
            status = status
        )

        return (
            name = config.name,
            dimension = length(config.p_true),
            status = status,
            success = success,
            recovery_error = recovery_error,
            objective_value = best_objective,
            time_seconds = total_time,
            stage1_time = stage1_time,
            stage2_time = stage2_time,
            critical_points = n_critical_points,
            success_rate = result[:success_rate],
            best_degree = result[:best_degree],
            error_message = nothing
        ), timing

    catch e
        total_time = time() - total_start

        timing = (
            model = config.name,
            stage1_time = nothing,
            stage2_time = nothing,
            total_time = total_time,
            status = "EXCEPTION"
        )

        return (
            name = config.name,
            dimension = length(config.p_true),
            status = "EXCEPTION",
            success = false,
            recovery_error = Inf,
            objective_value = Inf,
            time_seconds = total_time,
            stage1_time = nothing,
            stage2_time = nothing,
            critical_points = 0,
            success_rate = 0.0,
            best_degree = nothing,
            error_message = string(e)
        ), timing
    end
end

# ==============================================================================
# Report Generation
# ==============================================================================

"""Generate full Markdown report"""
function generate_markdown_report(results_df, timings, campaign_start, campaign_end, selected_models)
    io = IOBuffer()

    # Header
    println(io, "# N=4 Testing Campaign Report")
    println(io)
    println(io, "**Generated:** $(Dates.format(campaign_end, "yyyy-mm-dd HH:MM:SS"))")
    println(io)
    println(io, "## Campaign Configuration")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "| --- | --- |")
    println(io, "| Start Time | $(Dates.format(campaign_start, "yyyy-mm-dd HH:MM:SS")) |")
    println(io, "| End Time | $(Dates.format(campaign_end, "yyyy-mm-dd HH:MM:SS")) |")
    println(io, "| Total Duration | $(round((campaign_end - campaign_start).value / 1000 / 60, digits=2)) minutes |")
    println(io, "| Grid Parameter (GN) | $GN |")
    println(io, "| Degree Range | $(DEGREE_RANGE) |")
    println(io, "| Max Time per Model | $(MAX_TIME)s |")
    println(io, "| Models Tested | $(length(selected_models)) |")
    println(io)

    # Summary Statistics
    println(io, "## Summary Statistics")
    println(io)

    n_total = nrow(results_df)
    n_success = sum(results_df.success)
    n_excellent = sum(results_df.status .== "EXCELLENT")
    n_good = sum(results_df.status .== "GOOD")
    n_poor = sum(results_df.status .== "POOR")
    n_failed = sum(.!results_df.success)

    println(io, "| Metric | Count | Percentage |")
    println(io, "| --- | ---: | ---: |")
    println(io, "| Total Models | $n_total | 100% |")
    println(io, "| Successful | $n_success | $(round(n_success/n_total*100, digits=1))% |")
    println(io, "| Failed | $n_failed | $(round(n_failed/n_total*100, digits=1))% |")
    println(io)

    println(io, "### Quality Breakdown")
    println(io)
    println(io, "| Quality | Count | Criteria |")
    println(io, "| --- | ---: | --- |")
    println(io, "| EXCELLENT | $n_excellent | recovery < 0.01, objective < 1e-6 |")
    println(io, "| GOOD | $n_good | recovery < 0.1 |")
    println(io, "| POOR | $n_poor | recovery >= 0.1 |")
    println(io)

    # Timing Summary
    println(io, "## Timing Analysis")
    println(io)
    println(io, "### Per-Model Timing")
    println(io)
    println(io, timing_to_markdown(timings))
    println(io)

    # Overall timing statistics
    successful_times = [t.total_time for t in timings if t.status in ["EXCELLENT", "GOOD", "POOR"]]
    if !isempty(successful_times)
        println(io, "### Timing Statistics (Successful Runs)")
        println(io)
        println(io, "| Statistic | Value (seconds) |")
        println(io, "| --- | ---: |")
        println(io, "| Mean | $(round(mean(successful_times), digits=2)) |")
        println(io, "| Median | $(round(median(successful_times), digits=2)) |")
        println(io, "| Min | $(round(minimum(successful_times), digits=2)) |")
        println(io, "| Max | $(round(maximum(successful_times), digits=2)) |")
        println(io, "| Std Dev | $(round(std(successful_times), digits=2)) |")
        println(io)
    end

    # Detailed Results
    println(io, "## Detailed Results")
    println(io)

    # Select columns for display
    display_cols = [:name, :dimension, :status, :recovery_error, :objective_value, :critical_points, :best_degree]
    display_df = select(results_df, intersect(display_cols, names(results_df)))

    println(io, df_to_markdown(display_df))
    println(io)

    # Per-model details
    println(io, "## Per-Model Analysis")
    println(io)

    for row in eachrow(results_df)
        println(io, "### $(row.name)")
        println(io)
        println(io, "| Property | Value |")
        println(io, "| --- | --- |")
        println(io, "| Dimension | $(row.dimension) |")
        println(io, "| Status | **$(row.status)** |")
        println(io, "| Success | $(row.success) |")

        if row.recovery_error !== nothing && isfinite(row.recovery_error)
            println(io, "| Recovery Error | $(@sprintf("%.6e", row.recovery_error)) |")
        end

        if row.objective_value !== nothing && isfinite(row.objective_value)
            println(io, "| Objective Value | $(@sprintf("%.6e", row.objective_value)) |")
        end

        println(io, "| Critical Points | $(row.critical_points) |")
        println(io, "| Best Degree | $(something(row.best_degree, "-")) |")
        println(io, "| Total Time | $(round(row.time_seconds, digits=2))s |")

        if row.error_message !== nothing
            println(io, "| Error | $(row.error_message) |")
        end

        println(io)
    end

    # Footer
    println(io, "---")
    println(io)
    println(io, "*Report generated by Dynamic_objectives n4_hpc_campaign.jl*")

    return String(take!(io))
end

# ==============================================================================
# Main Execution
# ==============================================================================

function main()
    campaign_start = now()

    # Setup output directory
    output_dir = mkpath(OUTPUT_DIR)

    # Get selected models
    selected_models = get_selected_models()

    # Header
    println("=" ^ 80)
    println("N=4 HPC Testing Campaign")
    println("=" ^ 80)
    println("  Start time: $campaign_start")
    println("  Output directory: $output_dir")
    println("  Configuration: GN=$GN, degrees=$DEGREE_RANGE, max_time=$(MAX_TIME)s")
    println("  Models: $(join(selected_models, ", "))")
    println("=" ^ 80)

    # Get configurations for selected models
    test_configs_selected = [configs_by_name[name] for name in selected_models if haskey(configs_by_name, name)]

    if isempty(test_configs_selected)
        error("No valid model configurations found for: $selected_models")
    end

    println("\nRunning $(length(test_configs_selected)) model(s)...")

    # Run tests
    results = []
    timings = NamedTuple[]

    for (i, config) in enumerate(test_configs_selected)
        println("\n" * "-" ^ 80)
        println("[$i/$(length(test_configs_selected))] Testing: $(config.name)")
        println("  Description: $(config.description)")
        println("  Parameters: $(length(config.p_true))D")
        println("-" ^ 80)

        result, timing = test_model_with_timing(config, output_dir)
        push!(results, result)
        push!(timings, timing)

        # Print immediate feedback
        @printf "  Status: %s\n" result.status
        @printf "  Total time: %.2fs\n" result.time_seconds
        if result.success
            @printf "  Recovery error: %.6e\n" result.recovery_error
            @printf "  Objective: %.6e\n" result.objective_value
        end
    end

    campaign_end = now()

    # Convert results to DataFrame
    results_df = DataFrame(results)

    # Generate Markdown report
    println("\n" * "=" ^ 80)
    println("Generating Reports...")
    println("=" ^ 80)

    markdown_report = generate_markdown_report(results_df, timings, campaign_start, campaign_end, selected_models)

    # Save outputs
    report_file = joinpath(output_dir, "campaign_report.md")
    open(report_file, "w") do f
        write(f, markdown_report)
    end
    println("  Markdown report: $report_file")

    summary_file = joinpath(output_dir, "campaign_summary.csv")
    CSV.write(summary_file, results_df)
    println("  Summary CSV: $summary_file")

    timing_df = DataFrame(timings)
    timing_file = joinpath(output_dir, "timing_breakdown.csv")
    CSV.write(timing_file, timing_df)
    println("  Timing CSV: $timing_file")

    # Print summary to console
    println("\n" * "=" ^ 80)
    println("CAMPAIGN COMPLETE")
    println("=" ^ 80)

    n_total = nrow(results_df)
    n_success = sum(results_df.success)

    println("\nResults:")
    @printf "  Total: %d models\n" n_total
    @printf "  Successful: %d (%.1f%%)\n" n_success (n_success / n_total * 100)
    @printf "  Failed: %d (%.1f%%)\n" (n_total - n_success) ((n_total - n_success) / n_total * 100)

    println("\nTiming:")
    total_campaign_time = (campaign_end - campaign_start).value / 1000
    @printf "  Total campaign time: %.2f seconds (%.2f minutes)\n" total_campaign_time (total_campaign_time / 60)

    println("\nOutput files:")
    println("  $report_file")
    println("  $summary_file")
    println("  $timing_file")

    println("\n" * "=" ^ 80)

    return results_df, timings
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
