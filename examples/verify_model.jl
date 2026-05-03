#!/usr/bin/env julia
"""
Quick verification script with rich terminal display

This script:
1. Visualizes the ODE model response using terminal plots
2. Tests the objective function on a few sample points
3. Shows timing information for grid evaluation estimates
4. Benchmarks ODE tolerance sweep — speedup vs accuracy
5. Benchmarks ODE solver comparison — Tsit5 vs Vern7 vs Vern9
6. Demonstrates TolerantObjective mutable tolerance switching

Run this BEFORE the full integration test to verify everything works!
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using DynamicObjectives
using LinearAlgebra
using Printf
using ModelingToolkit: complete, ODEProblem
using DynamicObjectives: Vern9, Vern7, Tsit5

display_section("Model Verification with Rich Display")

# ==============================================================================
# Step 1: Create model and display summary
# ==============================================================================

display_section("Step 1: Model Setup")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Display model summary with rich formatting
display_model_summary(model)

# Display parameters
display_parameters([Symbol("α"), Symbol("β")], hcat(p_true), labels = ["True Values"])

# ==============================================================================
# Step 2: Generate and visualize time series
# ==============================================================================

display_section("Step 2: Visualizing Model Response")

# Generate time series data
time_interval = [0.0, 20.0]
numpoints = 30

# Create ODEProblem for sample_data
problem = ODEProblem(
    complete(model),
    merge(Dict(states .=> ic), Dict(params .=> p_true)),
    time_interval,
)

data = sample_data(problem, model, outputs, time_interval, p_true, ic, numpoints)

# Display time series with terminal plots
display_subsection("Time series for true parameters")
display_time_series(data)

# ==============================================================================
# Step 3: Test objective function
# ==============================================================================

display_section("Step 3: Testing Objective Function")

# Create objective function
# Note: eval_timeout removed - it adds 30x overhead per call due to @async/timedwait
objective = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints,
    L2_norm;
    return_inf_on_error = true,
)

# Test at true parameters (should be ~0)
display_subsection("Evaluating objective at different parameter values")
test_points = [
    (p_true, "True parameters"),
    ([1.1, 0.6], "Slightly perturbed"),
    ([1.5, 0.75], "Moderately perturbed"),
    ([2.0, 1.0], "Heavily perturbed"),
]

results = []
for (p, label) in test_points
    t_start = time()
    obj_val = objective(p)
    elapsed = time() - t_start
    push!(results, (label, p, obj_val, elapsed))

    status = obj_val < 1e-6 ? "PASS" : obj_val == Inf ? "FAIL" : "OK"
    println(
        "  [$status] $label: $(round(obj_val, digits=8)) ($(round(1000*elapsed, digits=1))ms)",
    )
end

# ==============================================================================
# Step 4: Performance estimate for grid search
# ==============================================================================

display_section("Step 4: Grid Search Performance Estimates")

# Use average timing from test evaluations
avg_time_per_eval = sum(r[4] for r in results) / length(results)

# Different grid sizes
grid_configs = [
    (5, "Very fast test"),
    (10, "Fast test"),
    (20, "Medium test"),
    (50, "Full grid (slow!)"),
    (100, "Production (very slow!)"),
]

display_subsection("Estimated time for different grid sizes (2D problem)")

# Build results for display
grid_results = Pair{String,Any}[]
for (GN, desc) in grid_configs
    n_points = GN^2
    est_time = avg_time_per_eval * n_points

    time_str = if est_time < 60
        "$(round(est_time, digits=1))s"
    elseif est_time < 3600
        "$(round(est_time/60, digits=1))min"
    else
        "$(round(est_time/3600, digits=1))hr"
    end

    push!(grid_results, "GN=$GN ($n_points pts)" => time_str)
end

display_results(grid_results, title = "Grid Size Estimates")

# ==============================================================================
# Step 5: Recommendations
# ==============================================================================

display_section("Recommendations")

display_results(
    [
        "Avg time per ODE eval" => "$(round(1000*avg_time_per_eval, digits=1))ms",
        "Quick test (GN=10)" => "$(round(avg_time_per_eval*100, digits=1))s",
        "Accurate (GN=50)" => "$(round(avg_time_per_eval*2500/60, digits=1))min",
    ],
    title = "Performance Summary",
)

# Calculate multi-degree estimate
single_degree_time = avg_time_per_eval * 100  # GN=10
total_time = single_degree_time * 3
println("\nFor degree_range = 4:6 (3 degrees) with GN=10:")
println("  Total time estimate: $(round(total_time, digits=1))s")

# ==============================================================================
# Step 6: ODE Tolerance Sweep Benchmark
# ==============================================================================

display_section("Step 6: ODE Tolerance Sweep Benchmark")
display_subsection("Comparing objective accuracy and speed across ODE solver tolerances")

tolerance_levels = [1e-4, 1e-6, 1e-8, 1e-10]
reference_value = objective(p_true)  # tight default (1e-10)

tol_results = []
for tol in tolerance_levels
    # Create objective at this tolerance
    obj_tol = make_error_distance(
        model,
        outputs,
        ic,
        p_true,
        time_interval,
        numpoints,
        L2_norm;
        return_inf_on_error = true,
        abstol = tol,
        reltol = tol,
    )

    # Benchmark: evaluate at several test points
    n_bench = 20
    test_p = [1.1, 0.6]  # fixed test point for consistency
    t_start = time()
    vals = [obj_tol(test_p) for _ in 1:n_bench]
    elapsed = (time() - t_start) / n_bench

    # Compare accuracy at true params
    val_at_true = obj_tol(p_true)

    push!(
        tol_results,
        (
            tol = tol,
            val_at_test = vals[1],
            val_at_true = val_at_true,
            avg_ms = 1000 * elapsed,
        ),
    )

    @printf(
        "  tol=%.0e: val@test=%.8e  val@true=%.2e  avg=%.1fms\n",
        tol,
        vals[1],
        val_at_true,
        1000 * elapsed
    )
end

# Show speedup relative to tightest tolerance
if length(tol_results) >= 2
    t_tight = tol_results[end].avg_ms
    println("\n  Speedup relative to tol=1e-10:")
    for r in tol_results
        speedup = t_tight / max(r.avg_ms, 1e-6)
        @printf("    tol=%.0e: %.1fx\n", r.tol, speedup)
    end
end

# ==============================================================================
# Step 7: ODE Solver Comparison Benchmark
# ==============================================================================

display_section("Step 7: ODE Solver Comparison Benchmark")
display_subsection("Comparing ODE solvers at different tolerance levels")

solvers =
    [(Tsit5(), "Tsit5 (RK5/4)"), (Vern7(), "Vern7 (RK7/6)"), (Vern9(), "Vern9 (RK9/8)")]

# Test at two tolerance regimes
for tol in [1e-6, 1e-10]
    println("\n  --- Tolerance = $(tol) ---")
    for (slvr, name) in solvers
        obj_s = make_error_distance(
            model,
            outputs,
            ic,
            p_true,
            time_interval,
            numpoints,
            L2_norm;
            return_inf_on_error = true,
            solver = slvr,
            abstol = tol,
            reltol = tol,
        )

        # Benchmark
        n_bench = 20
        test_p = [1.1, 0.6]
        t_start = time()
        vals = [obj_s(test_p) for _ in 1:n_bench]
        elapsed = (time() - t_start) / n_bench

        val_at_true = obj_s(p_true)

        @printf(
            "    %-16s: val@test=%.8e  val@true=%.2e  avg=%.1fms\n",
            name,
            vals[1],
            val_at_true,
            1000 * elapsed
        )
    end
end

# ==============================================================================
# Step 8: TolerantObjective Demo
# ==============================================================================

display_section("Step 8: TolerantObjective Demo")
display_subsection("Demonstrating mutable tolerance switching")

tol_obj = TolerantObjective(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints,
    L2_norm;
    return_inf_on_error = true,
    abstol = 1e-10,
    reltol = 1e-10,
)

test_p = [1.1, 0.6]

# Evaluate at tight tolerance
val_tight = tol_obj(test_p)
println("  Tight (1e-10): $(round(val_tight, digits=10))")

# Switch to coarse
set_tolerance!(tol_obj, 1e-4)
val_coarse = tol_obj(test_p)
println("  Coarse (1e-4): $(round(val_coarse, digits=10))")

# Switch back to tight
set_tolerance!(tol_obj, 1e-10)
val_tight2 = tol_obj(test_p)
println("  Tight again:   $(round(val_tight2, digits=10))")

# Verify consistency
diff = abs(val_tight - val_tight2)
println("  Consistency check (tight vs tight again): diff = $(diff)")
@assert diff < 1e-12 "TolerantObjective should produce identical results after round-trip"

# Switch solver
set_solver!(tol_obj, Tsit5())
val_tsit5 = tol_obj(test_p)
println("  Tsit5 (1e-10): $(round(val_tsit5, digits=10))")

display_section("Verification Complete")
println("Model is working correctly.")
println("Next: julia --project=. examples/test_simple_workflow.jl")
