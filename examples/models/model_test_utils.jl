# Shared utilities for per-model test scripts
# Include this file to get common functionality

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dynamic_objectives
using Globtim: Globtim, run_standard_experiment
using GlobtimPostProcessing: GlobtimPostProcessing, refine_experiment_results, ode_refinement_config
using LinearAlgebra
using ForwardDiff
using Dates
using Printf

# Include ExperimentCLI module for config
if !isdefined(Main, :ExperimentCLI)
    globtimcore_path = joinpath(@__DIR__, "..", "..", "..", "globtimcore")
    if isdir(globtimcore_path)
        include(joinpath(globtimcore_path, "src", "ExperimentCLI.jl"))
    end
end

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
    bounds::Vector{Tuple{Float64, Float64}}
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
    distance::Function = L2_norm

    # Output
    output_dir::String = ""
end

"""
    run_model_test(config::ModelTestConfig; verbose=true)

Run a complete model test with both stages and generate reports.

Returns a NamedTuple with all results.
"""
function run_model_test(config::ModelTestConfig; verbose::Bool=true)
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
        display_section("Model Test: $(config.name)", subtitle=config.description)
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
        display_parameters(param_names, hcat(config.p_true), labels=["True Values"])
    end

    # Create objective function
    objective = make_error_distance(
        model, outputs, config.ic, config.p_true,
        config.time_interval, config.numpoints,
        config.distance, first, nothing;
        return_inf_on_error = true
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
        display_subsection("Running globtimcore...")
    end

    # Configure experiment
    experiment_config = if @isdefined(ExperimentCLI)
        ExperimentCLI.ExperimentParams(
            domain_size = config.domain_size,
            GN = config.GN,
            degree_range = config.degree_range,
            max_time = config.max_time,
            basis = config.basis
        )
    else
        (
            domain_size = config.domain_size,
            GN = config.GN,
            degree_range = config.degree_range,
            max_time = config.max_time,
            basis = config.basis
        )
    end

    if verbose
        display_results([
            "Grid size (GN)" => "$(config.GN) ($(config.GN^n_params) points)",
            "Degree range" => "$(config.degree_range)",
            "Max time" => "$(config.max_time)s",
        ], title="Stage 1 Configuration")
    end

    # Run Stage 1
    stage1_start = time()
    result = run_standard_experiment(
        objective_function = objective,
        problem_params = nothing,
        domain_bounds = config.bounds,
        experiment_config = experiment_config,
        output_dir = output_dir,
        metadata = Dict{String, Any}("model" => config.name),
        true_params = config.p_true
    )
    stage1_time = time() - stage1_start

    best_raw_value = minimum([dr.best_objective for dr in result[:degree_results]])

    if verbose
        display_results([
            "Total critical points" => result[:total_critical_points],
            "Best raw value" => best_raw_value,
            "Time elapsed" => "$(round(stage1_time, digits=2))s",
        ], title="Stage 1 Complete")

        # Degree comparison
        degree_comparison = [(
            degree = dr.degree,
            n_critical_points = dr.n_critical_points,
            best_objective = dr.best_objective
        ) for dr in result[:degree_results]]
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
        show_progress = verbose
    )

    stage2_start = time()
    refined = refine_experiment_results(
        result[:output_dir],
        objective,
        refinement_config
    )
    stage2_time = time() - stage2_start

    if verbose
        display_results([
            "Raw points" => refined.n_raw,
            "Converged" => refined.n_converged,
            "Success rate" => "$(round(100*refined.n_converged/max(refined.n_raw,1), digits=1))%",
            "Mean improvement" => "$(round(refined.mean_improvement, digits=2))x",
            "Best refined value" => refined.best_refined_value,
            "Time elapsed" => "$(round(stage2_time, digits=2))s",
        ], title="Stage 2 Complete")

        display_quality_summary(refined, title="Critical Point Quality")
    end

    # ===========================================================================
    # Step 4: Display Top Critical Points
    # ===========================================================================

    if verbose && refined.n_converged > 0
        display_section("Step 4: Top Critical Points")

        display_top_critical_points(
            refined.refined_points[1:refined.n_converged],
            refined.refined_values[1:refined.n_converged],
            param_names;
            n=5,
            true_params=config.p_true,
            title="Top 5 Critical Points"
        )
    end

    # ===========================================================================
    # Step 5: Gradient Validation
    # ===========================================================================

    grad_norms = Float64[]
    if refined.n_converged > 0
        if verbose
            display_section("Step 5: Gradient Validation")
            pb = progress_bar(refined.n_converged; description="Computing gradients", width=30)
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
            finish_progress!(pb; description="Gradients computed")
            display_gradient_analysis(grad_norms, tolerance=1e-6)
        end
    end

    # ===========================================================================
    # Step 6: Parameter Recovery
    # ===========================================================================

    recovery_error = Inf
    best_params = Float64[]

    if refined.n_converged > 0
        best_params = refined.refined_points[refined.best_refined_idx]
        recovery_error = norm(best_params .- config.p_true) / norm(config.p_true)

        if verbose
            display_section("Step 6: Parameter Recovery")

            display_parameters(
                param_names,
                hcat(config.p_true, best_params),
                labels=["True", "Recovered"]
            )

            param_errors = abs.(best_params .- config.p_true) ./ max.(abs.(config.p_true), 1e-10)

            error_pairs = Pair{String, Any}[
                "Overall relative error" => @sprintf("%.2f%%", recovery_error * 100)
            ]
            for i in 1:n_params
                push!(error_pairs, "Parameter $i error" => @sprintf("%.2f%%", param_errors[i] * 100))
            end
            push!(error_pairs, "Best objective" => refined.best_refined_value)

            display_results(error_pairs, title="Recovery Metrics")

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
        display_results([
            "Total time" => "$(round(total_time, digits=2))s ($(round(total_time/60, digits=2)) min)",
            "Stage 1" => "$(round(stage1_time, digits=2))s",
            "Stage 2" => "$(round(stage2_time, digits=2))s",
            "Output directory" => output_dir,
        ], title="Summary")
    end

    # Generate markdown report
    report_path = generate_test_report(
        config, result, refined, grad_norms,
        stage1_time, stage2_time, output_dir, param_names
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
        report_path = report_path
    )
end

"""
    generate_test_report(...)

Generate a detailed markdown report for the model test.
"""
function generate_test_report(
    config, result, refined, grad_norms,
    stage1_time, stage2_time, output_dir, param_names
)
    io = IOBuffer()
    n_params = length(config.p_true)

    println(io, "# Model Test Report: $(config.name)")
    println(io)
    println(io, "**Generated:** $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io)
    println(io, "## Model Information")
    println(io)
    println(io, "| Property | Value |")
    println(io, "| --- | --- |")
    println(io, "| Name | $(config.name) |")
    println(io, "| Description | $(config.description) |")
    println(io, "| Parameters | $(n_params) |")
    println(io, "| Time Interval | $(config.time_interval) |")
    println(io, "| Sample Points | $(config.numpoints) |")
    println(io)

    println(io, "### True Parameters")
    println(io)
    println(io, "```")
    println(io, "p_true = $(config.p_true)")
    println(io, "```")
    println(io)

    println(io, "### Bounds")
    println(io)
    println(io, "```")
    for (i, b) in enumerate(config.bounds)
        println(io, "p[$i]: $(b[1]) to $(b[2])")
    end
    println(io, "```")
    println(io)

    println(io, "## Optimization Configuration")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "| --- | --- |")
    println(io, "| GN | $(config.GN) |")
    println(io, "| Grid Points | $(config.GN^n_params) |")
    println(io, "| Degree Range | $(config.degree_range) |")
    println(io, "| Max Time | $(config.max_time)s |")
    println(io, "| Refinement Time/Point | $(config.max_time_per_point)s |")
    println(io)

    println(io, "## Stage 1: Raw Critical Points")
    println(io)
    println(io, "| Metric | Value |")
    println(io, "| --- | --- |")
    println(io, "| Total Critical Points | $(result[:total_critical_points]) |")

    best_raw = minimum([dr.best_objective for dr in result[:degree_results]])
    println(io, "| Best Raw Objective | $(@sprintf("%.4e", best_raw)) |")
    println(io, "| Time | $(round(stage1_time, digits=2))s |")
    println(io)

    println(io, "### Per-Degree Results")
    println(io)
    println(io, "| Degree | Critical Points | Best Objective |")
    println(io, "| ---: | ---: | ---: |")
    for dr in result[:degree_results]
        println(io, "| $(dr.degree) | $(dr.n_critical_points) | $(@sprintf("%.4e", dr.best_objective)) |")
    end
    println(io)

    println(io, "## Stage 2: Refined Critical Points")
    println(io)
    println(io, "| Metric | Value |")
    println(io, "| --- | --- |")
    println(io, "| Raw Points | $(refined.n_raw) |")
    println(io, "| Converged | $(refined.n_converged) |")
    println(io, "| Success Rate | $(round(100*refined.n_converged/max(refined.n_raw,1), digits=1))% |")
    println(io, "| Mean Improvement | $(round(refined.mean_improvement, digits=2))x |")
    println(io, "| Best Refined Value | $(@sprintf("%.4e", refined.best_refined_value)) |")
    println(io, "| Time | $(round(stage2_time, digits=2))s |")
    println(io)

    if refined.n_converged > 0
        println(io, "## Top 5 Critical Points")
        println(io)
        println(io, format_critical_points_markdown(
            refined.refined_points[1:refined.n_converged],
            refined.refined_values[1:refined.n_converged],
            param_names;
            n=5,
            true_params=config.p_true
        ))
        println(io)

        println(io, "## Parameter Recovery")
        println(io)
        best_params = refined.refined_points[refined.best_refined_idx]
        recovery_error = norm(best_params .- config.p_true) / norm(config.p_true)

        println(io, "| Parameter | True | Recovered | Error |")
        println(io, "| --- | ---: | ---: | ---: |")
        for (i, name) in enumerate(param_names)
            err = abs(best_params[i] - config.p_true[i]) / max(abs(config.p_true[i]), 1e-10)
            println(io, "| $name | $(round(config.p_true[i], digits=4)) | $(round(best_params[i], digits=4)) | $(@sprintf("%.2f%%", err*100)) |")
        end
        println(io)
        println(io, "**Overall Relative Error:** $(@sprintf("%.2f%%", recovery_error * 100))")
        println(io)

        if !isempty(grad_norms)
            println(io, "## Gradient Validation")
            println(io)
            println(io, "| Metric | Value |")
            println(io, "| --- | ---: |")
            println(io, "| Valid points (||∇f|| < 1e-6) | $(count(n -> n < 1e-6, grad_norms))/$(length(grad_norms)) |")
            println(io, "| Min ||∇f|| | $(@sprintf("%.4e", minimum(grad_norms))) |")
            println(io, "| Mean ||∇f|| | $(@sprintf("%.4e", sum(grad_norms)/length(grad_norms))) |")
            println(io, "| Max ||∇f|| | $(@sprintf("%.4e", maximum(grad_norms))) |")
            println(io)
        end
    end

    println(io, "## Timing Summary")
    println(io)
    println(io, "| Stage | Time |")
    println(io, "| --- | ---: |")
    println(io, "| Stage 1 | $(round(stage1_time, digits=2))s |")
    println(io, "| Stage 2 | $(round(stage2_time, digits=2))s |")
    println(io, "| **Total** | **$(round(stage1_time + stage2_time, digits=2))s** |")
    println(io)

    println(io, "---")
    println(io, "*Generated by Dynamic_objectives model test framework*")

    # Save report
    report_path = joinpath(output_dir, "report.md")
    open(report_path, "w") do f
        write(f, String(take!(io)))
    end

    return report_path
end
