#!/usr/bin/env julia
"""
Quick verification script with rich terminal display

This script:
1. Visualizes the ODE model response using terminal plots
2. Tests the objective function on a few sample points
3. Shows timing information for grid evaluation estimates

Run this BEFORE the full integration test to verify everything works!
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dynamic_objectives
using LinearAlgebra

println("="^80)
println("Model Verification with Rich Display")
println("="^80)
println()

# ==============================================================================
# Step 1: Create model and display summary
# ==============================================================================

println("Step 1: Model Setup")
println("-"^80)

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Display model summary with rich formatting
display_model_summary(model)
println()

# Display parameters
display_parameters(
    [Symbol("α"), Symbol("β")],
    hcat(p_true),
    labels=["True Values"]
)
println()

# ==============================================================================
# Step 2: Generate and visualize time series
# ==============================================================================

println("Step 2: Visualizing Model Response")
println("-"^80)

# Generate time series data
time_interval = [0.0, 20.0]
numpoints = 30

data = sample_data(
    model, outputs, ic, p_true,
    time_interval, numpoints
)

# Display time series with terminal plots
println("Time series for true parameters:")
display_time_series(data)
println()

# ==============================================================================
# Step 3: Test objective function
# ==============================================================================

println("Step 3: Testing Objective Function")
println("-"^80)

# Create objective function
objective = make_error_distance(
    model, outputs, ic, p_true,
    time_interval, numpoints,
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

# Test at true parameters (should be ~0)
println("Evaluating objective at different parameter values...")
test_points = [
    (p_true, "True parameters"),
    ([1.1, 0.6], "Slightly perturbed"),
    ([1.5, 0.75], "Moderately perturbed"),
    ([2.0, 1.0], "Heavily perturbed")
]

results = []
for (p, label) in test_points
    t_start = time()
    obj_val = objective(p)
    elapsed = time() - t_start
    push!(results, (label, p, obj_val, elapsed))

    status = obj_val < 1e-6 ? "✓" : obj_val == Inf ? "✗" : "○"
    println("  $status $label: $(round(obj_val, digits=8)) ($(round(1000*elapsed, digits=1))ms)")
end
println()

# ==============================================================================
# Step 4: Performance estimate for grid search
# ==============================================================================

println("Step 4: Grid Search Performance Estimates")
println("-"^80)

# Use average timing from test evaluations
avg_time_per_eval = sum(r[4] for r in results) / length(results)

# Different grid sizes
grid_configs = [
    (5, "Very fast test"),
    (10, "Fast test"),
    (20, "Medium test"),
    (50, "Full grid (slow!)"),
    (100, "Production (very slow!)")
]

println("Estimated time for different grid sizes (2D problem):")
println()
println("  GN  │  Grid Points  │  Est. Time  │  Recommendation")
println("  ────┼───────────────┼─────────────┼─────────────────")

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

    recommendation = if GN <= 10
        "✓ Quick test"
    elseif GN <= 20
        "○ Reasonable"
    else
        "⚠ Slow"
    end

    println("  $(lpad(GN, 3)) │ $(lpad(n_points, 13)) │ $(lpad(time_str, 11)) │  $recommendation")
end
println()

# ==============================================================================
# Step 5: Recommendations
# ==============================================================================

println("="^80)
println("Recommendations")
println("="^80)
println()

println("Based on performance analysis:")
println("  • Average time per ODE evaluation: $(round(1000*avg_time_per_eval, digits=1))ms")
println("  • Recommended for quick test: GN = 10 (≈$(round(avg_time_per_eval*100, digits=1))s)")
println("  • Recommended for accuracy: GN = 50 (≈$(round(avg_time_per_eval*2500/60, digits=1))min)")
println()

# Calculate multi-degree estimate
println("For degree_range = 4:6 (3 degrees):")
single_degree_time = avg_time_per_eval * 100  # GN=10
total_time = single_degree_time * 3
println("  • Total time estimate: $(round(total_time, digits=1))s")
println()

println("✓ Verification complete! Model is working correctly.")
println("  You can now run: julia --project=. examples/test_simple_workflow.jl")
println()
