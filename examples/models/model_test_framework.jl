# Model Testing Framework
# Utilities for running and analyzing per-model tests
#
# Include this file in your model test scripts after activating the project:
#   using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))
#   include("model_test_framework.jl")

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing:
    GlobtimPostProcessing,
    refine_experiment_results,
    ode_refinement_config,
    check_l2_quality,
    detect_stagnation,
    check_objective_distribution_quality,
    has_ground_truth,
    compute_parameter_recovery_stats
using LinearAlgebra
using ForwardDiff
using Dates
using Printf

# Load ExperimentCLI
globtim_path = joinpath(@__DIR__, "..", "..", "..", "globtim")
isdir(globtim_path) ||
    error("globtim not found at $globtim_path — run from the monorepo root")
include(joinpath(globtim_path, "src", "ExperimentCLI.jl"))

"""
    ModelTestConfig

Configuration for a model test run.
"""
Base.@kwdef struct ModelTestConfig
    # Model identification
    name::String
    description::String

    # Model setup function (returns model, params, states, outputs)
    model_fn::Function

    # Problem parameters
    p_true::Vector{Float64}
    ic::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    time_interval::Vector{Float64} = [0.0, 20.0]
    numpoints::Int = 100

    # Optimization parameters
    GN::Int = 10
    degree_range::UnitRange{Int} = 6:10
    max_time::Float64 = 1800.0
    basis::Symbol = :chebyshev
    domain_size::Float64 = 1.5

    # Refinement parameters
    max_time_per_point::Float64 = 30.0

    # Distance function
    distance_function::Function = L2_norm

    # Output
    output_dir::String = ""
end

"""
    run_model_test(config::ModelTestConfig; verbose=true)

Run a complete model test with both stages and generate reports.

Returns a NamedTuple with all results.
"""
function run_model_test(config::ModelTestConfig; verbose::Bool = true)
    test_start = now()
    n_params = length(config.p_true)
    param_names = [Symbol("p$i") for i in 1:n_params]

    # Set up output directory
    output_dir = if isempty(config.output_dir)
        mkpath(joinpath(@__DIR__, "..", "..", "test_results", "models", config.name))
    else
        mkpath(config.output_dir)
    end

    if verbose
        display_section("Model Test: $(config.name)", subtitle = config.description)
    end

    # ===========================================================================
    # Step 1: Model Setup
    # ===========================================================================

    if verbose
        display_section("Step 1: Model Setup")
    end

    model, params, states, outputs = config.model_fn()

    if verbose
        display_model_summary(model)
        display_parameters(param_names, hcat(config.p_true), labels = ["True Values"])
    end

    # Create objective function
    objective = make_error_distance(
        model,
        outputs,
        config.ic,
        config.p_true,
        config.time_interval,
        config.numpoints,
        config.distance_function;
        return_inf_on_error = true,
    )

    # Verify objective
    obj_at_true = objective(config.p_true)
    if verbose
        println("Objective at true params: $(round(obj_at_true, digits=8))")
    end

    # ===========================================================================
    # Step 2: Stage 1 - Find Raw Critical Points
    # ===========================================================================

    if verbose
        display_section("Step 2: Stage 1 - Find Raw Critical Points")
        display_subsection("Running Globtim...")
    end

    # Configure experiment
    experiment_config = ExperimentCLI.ExperimentParams(
        domain_size = config.domain_size,
        GN = config.GN,
        degree_range = config.degree_range,
        max_time = config.max_time,
        basis = config.basis,
    )

    if verbose
        display_results(
            [
                "Grid size (GN)" => "$(config.GN) ($(config.GN^n_params) points)",
                "Degree range" => "$(config.degree_range)",
                "Max time" => "$(config.max_time)s",
            ],
            title = "Stage 1 Configuration",
        )
    end

    # Run Stage 1
    stage1_start = time()
    result = run_standard_experiment(
        objective_function = objective,
        objective_name = config.name,
        bounds = config.bounds,
        experiment_config = experiment_config,
        output_dir = output_dir,
        metadata = Dict{String,Any}(),
        true_params = config.p_true,
    )
    stage1_time = time() - stage1_start

    best_raw_value = minimum([dr.best_objective for dr in result[:degree_results]])

    if verbose
        display_results(
            [
                "Total critical points" => result[:total_critical_points],
                "Best raw value" => best_raw_value,
                "Time elapsed" => "$(round(stage1_time, digits=2))s",
            ],
            title = "Stage 1 Complete",
        )

        # Degree comparison
        degree_comparison = [
            (
                degree = dr.degree,
                n_critical_points = dr.n_critical_points,
                best_objective = dr.best_objective,
            ) for dr in result[:degree_results]
        ]
        display_degree_comparison(degree_comparison)
    end

    # ===========================================================================
    # Step 3: Stage 2 - Refine Critical Points
    # ===========================================================================

    if verbose
        display_section("Step 3: Stage 2 - Refine Critical Points")
        display_subsection("Running GlobtimPostProcessing...")
    end

    refinement_config = ode_refinement_config(
        max_time_per_point = config.max_time_per_point,
        show_progress = verbose,
    )

    stage2_start = time()
    refined = refine_experiment_results(result[:output_dir], objective, refinement_config)
    stage2_time = time() - stage2_start

    if verbose
        display_results(
            [
                "Raw points" => refined.n_raw,
                "Converged" => refined.n_converged,
                "Success rate" => "$(round(100*refined.n_converged/max(refined.n_raw,1), digits=1))%",
                "Mean improvement" => "$(round(refined.mean_improvement, digits=2))x",
                "Best refined value" => refined.best_refined_value,
                "Time elapsed" => "$(round(stage2_time, digits=2))s",
            ],
            title = "Stage 2 Complete",
        )

        display_quality_summary(refined, title = "Critical Point Quality")
    end

    # ===========================================================================
    # Step 4: Quality Diagnostics
    # ===========================================================================

    l2_result = nothing
    stagnation_result = nothing
    dist_result = nothing

    if verbose
        display_section("Step 4: Quality Diagnostics")
    end

    # L2 approximation quality
    try
        l2_result = check_l2_quality(result[:output_dir])
        if verbose
            grade_str = string(l2_result.grade)
            grade_color =
                l2_result.grade == :excellent ? "✓" :
                l2_result.grade == :good ? "○" : l2_result.grade == :acceptable ? "△" : "✗"
            display_results(
                [
                    "L2 Grade" => "$(grade_color) $(uppercase(grade_str))",
                    "L2 Error" => @sprintf("%.4e", l2_result.l2_error),
                ],
                title = "L2 Approximation Quality",
            )
        end
    catch e
        if verbose
            println("  L2 quality check skipped: $(e)")
        end
    end

    # Stagnation detection
    try
        stagnation_result = detect_stagnation(result[:output_dir])
        if verbose
            if stagnation_result.detected
                display_results(
                    [
                        "Status" => "⚠ Stagnation detected",
                        "Stagnation degree" => stagnation_result.stagnation_degree,
                    ],
                    title = "Convergence Analysis",
                )
            else
                display_results(
                    ["Status" => "✓ No stagnation detected"],
                    title = "Convergence Analysis",
                )
            end
        end
    catch e
        if verbose
            println("  Stagnation detection skipped: $(e)")
        end
    end

    # Objective distribution quality
    try
        dist_result = check_objective_distribution_quality(result[:output_dir])
        if verbose
            display_results(
                [
                    "Outliers" => "$(dist_result.n_outliers)/$(dist_result.n_points)",
                    "Distribution" =>
                        dist_result.is_healthy ? "✓ Healthy" : "⚠ Issues detected",
                ],
                title = "Objective Distribution",
            )
        end
    catch e
        if verbose
            println("  Distribution check skipped: $(e)")
        end
    end

    # ===========================================================================
    # Step 5: Display Top Critical Points
    # ===========================================================================

    if verbose && refined.n_converged > 0
        display_section("Step 5: Top Critical Points")

        display_top_critical_points(
            refined.refined_points[1:refined.n_converged],
            refined.refined_values[1:refined.n_converged],
            param_names;
            n = 5,
            true_params = config.p_true,
            title = "Top 5 Critical Points",
        )
    end

    # ===========================================================================
    # Step 6: Gradient Validation
    # ===========================================================================

    grad_norms = Float64[]
    if refined.n_converged > 0
        if verbose
            display_section("Step 6: Gradient Validation")
            pb = progress_bar(
                refined.n_converged;
                description = "Computing gradients",
                width = 30,
            )
        end

        for i in 1:refined.n_converged
            params_i = refined.refined_points[i]
            grad = ForwardDiff.gradient(objective, params_i)
            push!(grad_norms, norm(grad))
            if verbose
                update_progress!(pb)
            end
        end

        if verbose
            finish_progress!(pb; description = "Gradients computed")
            display_gradient_analysis(grad_norms, tolerance = 1e-6)
        end
    end

    # ===========================================================================
    # Step 7: Parameter Recovery
    # ===========================================================================

    recovery_error = Inf
    best_params = Float64[]

    if refined.n_converged > 0
        best_params = refined.refined_points[refined.best_refined_idx]
        recovery_error = norm(best_params .- config.p_true) / norm(config.p_true)

        if verbose
            display_section("Step 7: Parameter Recovery")

            display_parameters(
                param_names,
                hcat(config.p_true, best_params),
                labels = ["True", "Recovered"],
            )

            param_errors =
                abs.(best_params .- config.p_true) ./ max.(abs.(config.p_true), 1e-10)

            error_pairs = Pair{String,Any}["Overall relative error"=>@sprintf(
                "%.2f%%",
                recovery_error*100
            )]
            for i in 1:n_params
                push!(
                    error_pairs,
                    "Parameter $i error" => @sprintf("%.2f%%", param_errors[i] * 100)
                )
            end
            push!(error_pairs, "Best objective" => refined.best_refined_value)

            display_results(error_pairs, title = "Recovery Metrics")

            # Success message
            println()
            if recovery_error < 0.01
                println("SUCCESS: Excellent recovery (< 1% error)")
            elseif recovery_error < 0.05
                println("SUCCESS: Good recovery (< 5% error)")
            elseif recovery_error < 0.10
                println("ACCEPTABLE: Reasonable recovery (< 10% error)")
            else
                println("Needs improvement (> 10% error)")
            end
        end
    end

    test_end = now()
    total_time = stage1_time + stage2_time

    # ===========================================================================
    # Generate Report
    # ===========================================================================

    if verbose
        display_section("Test Complete")
        display_results(
            [
                "Total time" => "$(round(total_time, digits=2))s ($(round(total_time/60, digits=2)) min)",
                "Stage 1" => "$(round(stage1_time, digits=2))s",
                "Stage 2" => "$(round(stage2_time, digits=2))s",
                "Output directory" => output_dir,
            ],
            title = "Summary",
        )
    end

    # Generate markdown report
    report_path = generate_test_report(
        config,
        result,
        refined,
        grad_norms,
        stage1_time,
        stage2_time,
        output_dir,
        param_names;
        l2_result = l2_result,
        stagnation_result = stagnation_result,
        dist_result = dist_result,
    )

    if verbose
        println("\nReport saved to: $report_path")
    end

    return (
        config = config,
        stage1_result = result,
        refined = refined,
        grad_norms = grad_norms,
        stage1_time = stage1_time,
        stage2_time = stage2_time,
        total_time = total_time,
        recovery_error = recovery_error,
        best_params = best_params,
        output_dir = output_dir,
        report_path = report_path,
        # Quality diagnostics
        l2_result = l2_result,
        stagnation_result = stagnation_result,
        dist_result = dist_result,
    )
end

"""
    generate_test_report(...)

Generate a detailed markdown report for the model test.
"""
function generate_test_report(
    config,
    result,
    refined,
    grad_norms,
    stage1_time,
    stage2_time,
    output_dir,
    param_names;
    l2_result = nothing,
    stagnation_result = nothing,
    dist_result = nothing,
)
    io = IOBuffer()
    n_params = length(config.p_true)
    total_time = stage1_time + stage2_time

    # Compute recovery metrics
    recovery_error = Inf
    best_params = Float64[]
    if refined.n_converged > 0
        best_params = refined.refined_points[refined.best_refined_idx]
        recovery_error = norm(best_params .- config.p_true) / norm(config.p_true)
    end

    # Determine overall status
    status_emoji, status_text = if recovery_error < 0.01
        ("✅", "EXCELLENT")
    elseif recovery_error < 0.05
        ("✅", "SUCCESS")
    elseif recovery_error < 0.10
        ("⚠️", "ACCEPTABLE")
    else
        ("❌", "NEEDS IMPROVEMENT")
    end

    # ===========================================================================
    # Title and Executive Summary
    # ===========================================================================

    println(io, "# Model Test Report: $(config.name)")
    println(io)
    println(io, "> **Generated:** $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io)

    # Executive Summary Box
    println(io, "## $(status_emoji) Executive Summary")
    println(io)
    println(io, "| Metric | Value | Status |")
    println(io, "| :--- | :---: | :---: |")

    # Recovery error row
    recovery_status =
        recovery_error < 0.01 ? "🟢" :
        recovery_error < 0.05 ? "🟢" : recovery_error < 0.10 ? "🟡" : "🔴"
    println(
        io,
        "| **Parameter Recovery** | $(@sprintf("%.2f%%", recovery_error * 100)) error | $(recovery_status) |",
    )

    # Best objective row
    obj_status =
        refined.best_refined_value < 1e-6 ? "🟢" :
        refined.best_refined_value < 1e-3 ? "🟡" : "🔴"
    println(
        io,
        "| **Best Objective** | $(@sprintf("%.4e", refined.best_refined_value)) | $(obj_status) |",
    )

    # Convergence rate row
    conv_rate = refined.n_converged / max(refined.n_raw, 1)
    conv_status = conv_rate > 0.8 ? "🟢" : conv_rate > 0.5 ? "🟡" : "🔴"
    println(
        io,
        "| **Refinement Convergence** | $(@sprintf("%.0f%%", conv_rate * 100)) ($(refined.n_converged)/$(refined.n_raw)) | $(conv_status) |",
    )

    # Total time
    println(io, "| **Total Time** | $(round(total_time, digits=1))s | ⏱️ |")
    println(io)

    # Visual progress bar for convergence
    conv_bar = _make_progress_bar(conv_rate, 20)
    println(io, "**Convergence:** `$(conv_bar)` $(round(conv_rate * 100, digits=0))%")
    println(io)

    # ===========================================================================
    # Model Information
    # ===========================================================================

    println(io, "## Model Information")
    println(io)
    println(io, "**$(config.description)**")
    println(io)
    println(io, "| Property | Value |")
    println(io, "| :--- | :--- |")
    println(io, "| Dimensions | $(n_params) parameters |")
    println(
        io,
        "| Time Interval | `[$(config.time_interval[1]), $(config.time_interval[2])]` |",
    )
    println(io, "| Sample Points | $(config.numpoints) |")
    println(io)

    # Parameters table with bounds
    println(io, "### Parameters and Bounds")
    println(io)
    println(io, "| Parameter | True Value | Lower | Upper | Range |")
    println(io, "| :---: | ---: | ---: | ---: | :--- |")
    for (i, name) in enumerate(param_names)
        lower, upper = config.bounds[i]
        range_width = upper - lower
        true_pos = (config.p_true[i] - lower) / range_width
        pos_bar = _make_position_bar(true_pos, 10)
        println(
            io,
            "| **$(name)** | $(round(config.p_true[i], digits=4)) | $(lower) | $(upper) | $(pos_bar) |",
        )
    end
    println(io)

    # ===========================================================================
    # Optimization Configuration
    # ===========================================================================

    println(io, "## Optimization Configuration")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "| :--- | ---: |")
    println(io, "| Grid Size (GN) | $(config.GN) |")
    println(io, "| Total Grid Points | $(config.GN^n_params) |")
    println(
        io,
        "| Polynomial Degrees | $(first(config.degree_range)) to $(last(config.degree_range)) |",
    )
    println(io, "| Stage 1 Max Time | $(config.max_time)s |")
    println(io, "| Refinement Time/Point | $(config.max_time_per_point)s |")
    println(io)

    # ===========================================================================
    # Stage 1 Results
    # ===========================================================================

    println(io, "## Stage 1: Global Search")
    println(io)

    best_raw = minimum([dr.best_objective for dr in result[:degree_results]])

    println(io, "| Metric | Value |")
    println(io, "| :--- | ---: |")
    println(io, "| Critical Points Found | $(result[:total_critical_points]) |")
    println(io, "| Best Raw Objective | $(@sprintf("%.4e", best_raw)) |")
    println(io, "| Time | $(round(stage1_time, digits=2))s |")
    println(io)

    # Convergence trend with sparkline
    objectives = [dr.best_objective for dr in result[:degree_results]]
    sparkline = _make_sparkline(objectives)
    println(io, "**Convergence Trend:** $(sparkline)")
    println(io)

    # Per-degree results in collapsible section
    println(io, "<details>")
    println(io, "<summary><strong>Per-Degree Results</strong> (click to expand)</summary>")
    println(io)
    println(io, "| Degree | Critical Points | Best Objective | Trend |")
    println(io, "| ---: | ---: | ---: | :---: |")

    prev_obj = Inf
    for dr in result[:degree_results]
        trend =
            dr.best_objective < prev_obj * 0.99 ? "↗" :
            dr.best_objective > prev_obj * 1.01 ? "↘" : "→"
        println(
            io,
            "| $(dr.degree) | $(dr.n_critical_points) | $(@sprintf("%.4e", dr.best_objective)) | $(trend) |",
        )
        prev_obj = dr.best_objective
    end
    println(io)
    println(io, "</details>")
    println(io)

    # ===========================================================================
    # Stage 2 Results
    # ===========================================================================

    println(io, "## Stage 2: Local Refinement")
    println(io)
    println(io, "| Metric | Value |")
    println(io, "| :--- | ---: |")
    println(io, "| Input Points | $(refined.n_raw) |")
    println(io, "| Converged | $(refined.n_converged) |")
    println(io, "| Success Rate | $(round(100*conv_rate, digits=1))% |")
    println(io, "| Mean Improvement | $(round(refined.mean_improvement, digits=1))× |")
    println(io, "| Best Refined Value | $(@sprintf("%.4e", refined.best_refined_value)) |")
    println(io, "| Time | $(round(stage2_time, digits=2))s |")
    println(io)

    # Improvement visualization
    if best_raw > 0 && refined.best_refined_value > 0
        improvement = best_raw / refined.best_refined_value
        println(
            io,
            "**Improvement:** Raw → Refined = $(@sprintf("%.4e", best_raw)) → $(@sprintf("%.4e", refined.best_refined_value)) ($(round(improvement, digits=1))× better)",
        )
        println(io)
    end

    # ===========================================================================
    # Quality Diagnostics
    # ===========================================================================

    has_diagnostics =
        !isnothing(l2_result) || !isnothing(stagnation_result) || !isnothing(dist_result)
    if has_diagnostics
        println(io, "## Quality Diagnostics")
        println(io)

        # L2 Quality with star rating
        if !isnothing(l2_result)
            stars =
                l2_result.grade == :excellent ? "★★★★" :
                l2_result.grade == :good ? "★★★☆" :
                l2_result.grade == :acceptable ? "★★☆☆" : "★☆☆☆"
            grade_emoji =
                l2_result.grade == :excellent ? "🟢" :
                l2_result.grade == :good ? "🟢" : l2_result.grade == :acceptable ? "🟡" : "🔴"
            println(io, "### L2 Approximation Quality")
            println(io)
            println(io, "| Grade | Rating | Error |")
            println(io, "| :--- | :---: | ---: |")
            println(
                io,
                "| $(grade_emoji) **$(uppercase(string(l2_result.grade)))** | $(stars) | $(@sprintf("%.4e", l2_result.l2_error)) |",
            )
            println(io)
        end

        # Stagnation Analysis
        if !isnothing(stagnation_result)
            println(io, "### Convergence Analysis")
            println(io)
            if stagnation_result.detected
                println(
                    io,
                    "⚠️ **Stagnation detected** at degree $(stagnation_result.stagnation_degree)",
                )
                println(io)
                println(
                    io,
                    "> Consider increasing the polynomial degree range or adjusting grid density.",
                )
            else
                println(io, "✅ **No stagnation detected** - Convergence is healthy")
            end
            println(io)
        end

        # Distribution Quality
        if !isnothing(dist_result)
            println(io, "### Objective Distribution")
            println(io)
            outlier_pct = dist_result.n_outliers / max(dist_result.n_points, 1) * 100
            if dist_result.is_healthy
                println(
                    io,
                    "✅ **Healthy distribution** ($(dist_result.n_outliers) outliers, $(@sprintf("%.1f%%", outlier_pct)))",
                )
            else
                println(
                    io,
                    "⚠️ **$(dist_result.n_outliers)/$(dist_result.n_points) outliers** ($(@sprintf("%.1f%%", outlier_pct)))",
                )
                println(io)
                println(
                    io,
                    "> High outlier count may indicate numerical instability or difficult landscape.",
                )
            end
            println(io)
        end
    end

    # ===========================================================================
    # Parameter Recovery (Enhanced)
    # ===========================================================================

    if refined.n_converged > 0
        println(io, "## Parameter Recovery")
        println(io)

        # Overall status
        println(
            io,
            "**Status:** $(status_emoji) $(status_text) ($(@sprintf("%.2f%%", recovery_error * 100)) relative error)",
        )
        println(io)

        # Detailed parameter table
        println(io, "| Parameter | True | Recovered | Error | Position in Bounds |")
        println(io, "| :---: | ---: | ---: | ---: | :--- |")
        for (i, name) in enumerate(param_names)
            err = abs(best_params[i] - config.p_true[i]) / max(abs(config.p_true[i]), 1e-10)
            err_emoji = err < 0.01 ? "🟢" : err < 0.05 ? "🟡" : "🔴"

            lower, upper = config.bounds[i]
            range_width = upper - lower
            rec_pos = clamp((best_params[i] - lower) / range_width, 0.0, 1.0)
            true_pos = (config.p_true[i] - lower) / range_width
            pos_bar = _make_comparison_bar(true_pos, rec_pos, 12)

            println(
                io,
                "| **$(name)** | $(round(config.p_true[i], digits=4)) | $(round(best_params[i], digits=4)) | $(err_emoji) $(@sprintf("%.2f%%", err*100)) | $(pos_bar) |",
            )
        end
        println(io)
        println(io, "*Position bar: `○` = true value, `●` = recovered value*")
        println(io)

        # ===========================================================================
        # Top Critical Points
        # ===========================================================================

        println(io, "## Top Critical Points")
        println(io)
        println(
            io,
            format_critical_points_markdown(
                refined.refined_points[1:refined.n_converged],
                refined.refined_values[1:refined.n_converged],
                param_names;
                n = 5,
                true_params = config.p_true,
            ),
        )
        println(io)

        # ===========================================================================
        # Gradient Validation
        # ===========================================================================

        if !isempty(grad_norms)
            println(io, "## Gradient Validation")
            println(io)

            n_valid = count(n -> n < 1e-6, grad_norms)
            valid_pct = n_valid / length(grad_norms) * 100
            valid_status = valid_pct > 90 ? "🟢" : valid_pct > 70 ? "🟡" : "🔴"

            println(io, "| Metric | Value | Status |")
            println(io, "| :--- | ---: | :---: |")
            println(
                io,
                "| Valid points (‖∇f‖ < 1e-6) | $(n_valid)/$(length(grad_norms)) ($(@sprintf("%.0f%%", valid_pct))) | $(valid_status) |",
            )
            println(io, "| Min ‖∇f‖ | $(@sprintf("%.4e", minimum(grad_norms))) | |")
            println(
                io,
                "| Mean ‖∇f‖ | $(@sprintf("%.4e", sum(grad_norms)/length(grad_norms))) | |",
            )
            println(io, "| Max ‖∇f‖ | $(@sprintf("%.4e", maximum(grad_norms))) | |")
            println(io)
        end
    end

    # ===========================================================================
    # Timing Summary (Visual)
    # ===========================================================================

    println(io, "## Timing Summary")
    println(io)

    stage1_pct = stage1_time / total_time
    stage2_pct = stage2_time / total_time

    println(io, "| Stage | Time | Share |")
    println(io, "| :--- | ---: | :--- |")
    println(
        io,
        "| Stage 1 (Global) | $(round(stage1_time, digits=2))s | $(_make_progress_bar(stage1_pct, 15)) $(@sprintf("%.0f%%", stage1_pct*100)) |",
    )
    println(
        io,
        "| Stage 2 (Local) | $(round(stage2_time, digits=2))s | $(_make_progress_bar(stage2_pct, 15)) $(@sprintf("%.0f%%", stage2_pct*100)) |",
    )
    println(io, "| **Total** | **$(round(total_time, digits=2))s** | |")
    println(io)

    # Performance metrics
    points_per_sec = result[:total_critical_points] / max(stage1_time, 0.1)
    println(
        io,
        "**Performance:** $(@sprintf("%.1f", points_per_sec)) critical points/second",
    )
    println(io)

    # ===========================================================================
    # Output Files
    # ===========================================================================

    println(io, "## Output Files")
    println(io)
    println(io, "| File | Description |")
    println(io, "| :--- | :--- |")
    println(io, "| `report.md` | This report |")

    # List CSV files
    csv_files = filter(f -> endswith(f, ".csv"), readdir(output_dir))
    for f in csv_files
        desc = if contains(f, "refined")
            "Refined critical points"
        elseif contains(f, "critical")
            "Raw critical points"
        elseif contains(f, "comparison")
            "Raw vs refined comparison"
        else
            "Data file"
        end
        println(io, "| `$(f)` | $(desc) |")
    end

    # List JSON files
    json_files = filter(f -> endswith(f, ".json"), readdir(output_dir))
    for f in json_files
        desc = if contains(f, "summary")
            "Refinement statistics"
        elseif contains(f, "metadata")
            "Experiment metadata"
        else
            "Configuration/results"
        end
        println(io, "| `$(f)` | $(desc) |")
    end
    println(io)

    # ===========================================================================
    # Recommendations
    # ===========================================================================

    recommendations = String[]

    if recovery_error > 0.10
        push!(
            recommendations,
            "❗ **High recovery error** - Consider increasing `GN` (grid density) or extending the `degree_range`",
        )
    end

    if !isnothing(stagnation_result) && stagnation_result.detected
        push!(
            recommendations,
            "⚠️ **Stagnation detected** - Try higher polynomial degrees or adjust domain bounds",
        )
    end

    if !isnothing(dist_result) && !dist_result.is_healthy
        push!(
            recommendations,
            "⚠️ **Outliers detected** - Check for numerical issues or narrow the search bounds",
        )
    end

    if conv_rate < 0.5
        push!(
            recommendations,
            "⚠️ **Low convergence rate** - Increase `max_time_per_point` for refinement",
        )
    end

    if !isempty(grad_norms)
        valid_pct = count(n -> n < 1e-6, grad_norms) / length(grad_norms)
        if valid_pct < 0.7
            push!(
                recommendations,
                "⚠️ **Gradient validation issues** - Some critical points may not be true minima",
            )
        end
    end

    if !isempty(recommendations)
        println(io, "## Recommendations")
        println(io)
        for rec in recommendations
            println(io, "- $(rec)")
        end
        println(io)
    end

    # ===========================================================================
    # Footer
    # ===========================================================================

    println(io, "---")
    println(io)
    println(io, "*Generated by Dynamic_objectives model test framework*")
    println(io)
    println(io, "**Output directory:** `$(output_dir)`")

    # Save report
    report_path = joinpath(output_dir, "report.md")
    open(report_path, "w") do f
        write(f, String(take!(io)))
    end

    return report_path
end

# ===========================================================================
# Helper functions for visual elements
# ===========================================================================

"""
    _make_progress_bar(fraction, width)

Create a Unicode progress bar.
"""
function _make_progress_bar(fraction::Real, width::Int)
    filled = round(Int, fraction * width)
    empty = width - filled
    return "█" ^ filled * "░" ^ empty
end

"""
    _make_position_bar(position, width)

Create a position indicator bar showing where a value lies in a range.
"""
function _make_position_bar(position::Real, width::Int)
    pos = clamp(round(Int, position * (width - 1)) + 1, 1, width)
    bar = collect("─" ^ width)
    bar[pos] = '◆'
    return "`" * join(bar) * "`"
end

"""
    _make_comparison_bar(true_pos, recovered_pos, width)

Create a bar comparing true and recovered positions.
"""
function _make_comparison_bar(true_pos::Real, recovered_pos::Real, width::Int)
    bar = collect("─" ^ width)

    true_idx = clamp(round(Int, true_pos * (width - 1)) + 1, 1, width)
    rec_idx = clamp(round(Int, recovered_pos * (width - 1)) + 1, 1, width)

    if true_idx == rec_idx
        bar[true_idx] = '◉'  # Overlapping
    else
        bar[true_idx] = '○'   # True value (hollow)
        bar[rec_idx] = '●'    # Recovered value (filled)
    end

    return "`" * join(bar) * "`"
end

"""
    _make_sparkline(values)

Create a Unicode sparkline from a series of values.
"""
function _make_sparkline(values::Vector)
    if isempty(values)
        return ""
    end

    # Use log scale for better visualization of objective values
    log_vals = [v > 0 ? log10(v) : -20 for v in values]

    min_val = minimum(log_vals)
    max_val = maximum(log_vals)
    range_val = max_val - min_val

    if range_val ≈ 0
        return "▅" ^ length(values)
    end

    blocks = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█']

    sparkline = ""
    for v in log_vals
        # Invert: lower objective = higher bar
        normalized = 1.0 - (v - min_val) / range_val
        idx = clamp(round(Int, normalized * 7) + 1, 1, 8)
        sparkline *= blocks[idx]
    end

    return sparkline
end
