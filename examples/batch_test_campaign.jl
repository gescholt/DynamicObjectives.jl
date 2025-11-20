#!/usr/bin/env julia
# Batch Testing Campaign
# Usage: julia --project=. examples/batch_test_campaign.jl

using Dynamic_objectives
using Printf
using Dates

println("="^80)
println("Dynamic_objectives - Batch Testing Campaign")
println("="^80)

# ==============================================================================
# TEST CONFIGURATIONS
# ==============================================================================

# Define test cases (from TESTING_GUIDE.md)
test_cases = [
    # PHASE 1: EASY (2 parameters)
    (
        name = "LV_2D_v3_2outputs",
        phase = "EASY",
        model_fn = define_lotka_volterra_2D_model_v3_two_outputs,
        p_true = [1.0, 0.5],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance = L2_norm,
        timeout = nothing
    ),
    (
        name = "LV_2D_v1",
        phase = "EASY",
        model_fn = define_lotka_volterra_2D_model,
        p_true = [0.5, -0.3],
        ic = [1.0, 0.5],
        bounds = [(-1.0, 2.0), (-1.0, 0.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance = L2_norm,
        timeout = nothing
    ),

    # PHASE 2: MEDIUM (3-4 parameters)
    (
        name = "LV_3D_v2",
        phase = "MEDIUM",
        model_fn = define_lotka_volterra_3D_model_v2,
        p_true = [1.0, 0.5, 0.3],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)],
        time_interval = [0.0, 20.0],
        numpoints = 40,
        distance = L2_norm,
        timeout = nothing
    ),
    (
        name = "DAISY_Ex3_no_input",
        phase = "MEDIUM",
        model_fn = define_daisy_ex3_model_4D_no_input,
        p_true = [0.5, 0.3, 0.4, 0.2],
        ic = [1.0, 0.5, 0.3],
        bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)],
        time_interval = [0.0, 20.0],
        numpoints = 50,
        distance = L2_norm,
        timeout = nothing
    ),

    # PHASE 3: HARD
    (
        name = "FitzHugh_Nagumo",
        phase = "HARD",
        model_fn = define_fitzhugh_nagumo_3D_model,
        p_true = [0.8, 0.7, 0.8],
        ic = [0.0, 0.0],
        bounds = [(0.1, 2.0), (0.0, 2.0), (0.0, 2.0)],
        time_interval = [0.0, 50.0],
        numpoints = 100,
        distance = L2_norm,
        timeout = 10.0
    ),
    (
        name = "Simple_2D_square",
        phase = "HARD",
        model_fn = define_simple_2D_model_locally_identifiable_square,
        p_true = [0.5, 2.0],
        ic = [1.0],
        bounds = [(-1.0, 2.0), (-3.0, 3.0)],
        time_interval = [0.0, 5.0],
        numpoints = 20,
        distance = L2_norm,
        timeout = nothing
    ),
]

# ==============================================================================
# OPTIMIZER PLACEHOLDER
# ==============================================================================

"""
Replace this with your actual globtim optimizer call
"""
function run_optimizer(error_func, bounds, config)
    # PLACEHOLDER: Replace with your globtim call
    # Example:
    # result = globtim_optimize(
    #     objective = error_func,
    #     bounds = bounds,
    #     max_evals = config.max_evals,
    #     # ... your options ...
    # )
    # return result

    # For now, just do a simple test
    println("    [PLACEHOLDER] Would call globtim optimizer here")
    println("    [PLACEHOLDER] Returning dummy result")

    return (
        minimizer = config.p_true .+ randn(length(config.p_true)) * 0.1,
        minimum = rand() * 0.01,
        converged = true,
        iterations = 100,
        f_calls = 500
    )
end

# ==============================================================================
# RUN CAMPAIGN
# ==============================================================================

results = []
start_time = now()

println("\nStarting test campaign at $(Dates.format(start_time, "HH:MM:SS"))")
println("Testing $(length(test_cases)) models\n")

for (i, config) in enumerate(test_cases)
    println("="^80)
    println("[$i/$(length(test_cases))] $(config.name) ($(config.phase))")
    println("="^80)

    test_start = now()

    try
        # Setup
        println("  Setting up...")
        model, params, states, outputs = config.model_fn()

        error_func = make_error_distance(
            model, outputs, config.ic, config.p_true,
            config.time_interval, config.numpoints,
            config.distance, first, nothing;
            return_inf_on_error = true,
            eval_timeout = config.timeout
        )

        # Verify
        error_at_true = error_func(config.p_true)
        @printf "  ✓ Setup complete (error @ true: %.6e)\n" error_at_true

        if error_at_true > 1e-6
            @warn "  ⚠ Error at true parameters is not near zero!"
        end

        # Run optimizer
        println("  Running optimizer...")
        result = run_optimizer(error_func, config.bounds, config)

        # Evaluate results
        p_best = result.minimizer
        error_best = result.minimum
        param_error = norm(p_best - config.p_true)

        test_time = (now() - test_start).value / 1000  # seconds

        # Store results
        push!(results, (
            name = config.name,
            phase = config.phase,
            success = result.converged,
            error_best = error_best,
            param_error = param_error,
            f_calls = result.f_calls,
            time_sec = test_time,
            p_true = config.p_true,
            p_best = p_best
        ))

        # Print summary
        println("\n  Results:")
        @printf "    Best error: %.6e\n" error_best
        @printf "    Parameter error: %.6e\n" param_error
        @printf "    Function calls: %d\n" result.f_calls
        @printf "    Time: %.2f seconds\n" test_time
        println("    Converged: $(result.converged ? "✓" : "✗")")

        if error_best < 1e-4
            println("\n  ✓ SUCCESS: Found good solution!")
        else
            println("\n  ⚠ WARNING: Did not converge to optimal solution")
        end

    catch e
        println("  ✗ FAILED: $e")
        push!(results, (
            name = config.name,
            phase = config.phase,
            success = false,
            error_best = Inf,
            param_error = Inf,
            f_calls = 0,
            time_sec = 0.0,
            p_true = config.p_true,
            p_best = nothing
        ))
    end

    println()
end

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

total_time = (now() - start_time).value / 1000  # seconds

println("="^80)
println("CAMPAIGN SUMMARY")
println("="^80)

# Overall stats
n_total = length(results)
n_success = count(r -> r.success && r.error_best < 1e-4, results)
n_converged = count(r -> r.success, results)

println("\nOverall Performance:")
@printf "  Total tests: %d\n" n_total
@printf "  Successful: %d (%.1f%%)\n" n_success (n_success/n_total*100)
@printf "  Converged: %d (%.1f%%)\n" n_converged (n_converged/n_total*100)
@printf "  Total time: %.2f seconds\n" total_time

# By phase
println("\nBy Difficulty Phase:")
for phase in ["EASY", "MEDIUM", "HARD"]
    phase_results = filter(r -> r.phase == phase, results)
    if !isempty(phase_results)
        n = length(phase_results)
        n_succ = count(r -> r.success && r.error_best < 1e-4, phase_results)
        @printf "  %s: %d/%d succeeded (%.1f%%)\n" phase n_succ n (n_succ/n*100)
    end
end

# Detailed results table
println("\nDetailed Results:")
println("-"^80)
@printf "%-25s %-8s %-12s %-12s %-8s %-8s\n" "Model" "Phase" "Error" "Param Err" "Calls" "Time(s)"
println("-"^80)

for r in results
    status = r.success && r.error_best < 1e-4 ? "✓" : "✗"
    @printf "%s %-23s %-8s " status r.name r.phase

    if isfinite(r.error_best)
        @printf "%.4e  %.4e  " r.error_best r.param_error
    else
        @printf "%-12s %-12s " "FAILED" "FAILED"
    end

    @printf "%-8d %.2f\n" r.f_calls r.time_sec
end

println("-"^80)

# Save results to file
results_file = "test_results_$(Dates.format(start_time, "yyyymmdd_HHMMSS")).txt"
open(results_file, "w") do io
    println(io, "Dynamic_objectives Test Campaign Results")
    println(io, "Started: $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
    println(io, "="^80)
    println(io)

    for r in results
        println(io, "Model: $(r.name)")
        println(io, "  Phase: $(r.phase)")
        println(io, "  Success: $(r.success)")
        println(io, "  Best error: $(r.error_best)")
        println(io, "  Parameter error: $(r.param_error)")
        println(io, "  Function calls: $(r.f_calls)")
        println(io, "  Time: $(r.time_sec) seconds")
        println(io, "  True params: $(r.p_true)")
        println(io, "  Best params: $(r.p_best)")
        println(io)
    end
end

println("\n✓ Results saved to: $results_file")

println("\n" * "="^80)
println("Campaign completed!")
println("="^80)
