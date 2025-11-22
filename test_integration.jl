#!/usr/bin/env julia
"""
Quick test of run_globtim_optimization with display_optimization_progress integration

This script tests the newly added run_globtim_optimization() function with
real-time progress monitoring.
"""

using Pkg
Pkg.activate(@__DIR__)

println("Loading Dynamic_objectives...")
using Dynamic_objectives
using LinearAlgebra

println("\n" * "="^80)
println("Testing run_globtim_optimization with Real-Time Progress")
println("="^80)

# ==============================================================================
# Setup: Create a simple 2D model
# ==============================================================================

println("\n[1/3] Setting up model...")

# Use the simplest model: Lotka-Volterra 2D v3 with 2 outputs
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Define true parameters and bounds
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

println("  ✓ Model: Lotka-Volterra 2D v3 (2 outputs)")
println("  ✓ True parameters: $p_true")
println("  ✓ Bounds: $bounds")

# ==============================================================================
# Create objective function
# ==============================================================================

println("\n[2/3] Creating objective function...")

objective = create_globtim_objective(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30;
    distance = L2_norm,
    eval_timeout = nothing
)

# Verify objective works
error_at_true = objective(p_true, nothing)
println("  ✓ Objective created")
println("  ✓ Error at true params: $(round(error_at_true, digits=8))")

if error_at_true > 1e-8
    @warn "Error at true parameters is unexpectedly high!"
end

# ==============================================================================
# Run optimization with progress display
# ==============================================================================

println("\n[3/3] Running optimization with real-time progress monitoring...")

# NOTE: This will show real-time progress including:
# - Grid evaluation progress
# - display_optimization_progress() during critical point evaluation
# - Summary with timing and recovery metrics

result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 6,  # Smaller grid for quick test (6^2 = 36 points)
    degree_range = 4:5,  # Just 2 degrees for quick test
    output_dir = joinpath(@__DIR__, "test_output"),
    model_name = "LV_2D_test",
    show_progress = true,  # Enable real-time progress display
    save_results = false   # Don't save files for this test
)

# ==============================================================================
# Validate results
# ==============================================================================

println("\n" * "="^80)
println("Validation")
println("="^80)

# Check that we got results
if result[:best_params] !== nothing
    println("  ✓ Found solution")
    println("    Best parameters: $(result[:best_params])")
    println("    True parameters: $p_true")
    println("    Recovery error: $(round(100 * result[:recovery_error], digits=2))%")
    println("    Objective value: $(round(result[:best_objective], digits=6))")

    # Check quality
    if result[:recovery_error] < 0.01
        println("\n  ✅ EXCELLENT: Recovery error < 1%")
    elseif result[:recovery_error] < 0.1
        println("\n  ✅ GOOD: Recovery error < 10%")
    else
        println("\n  ⚠️  ACCEPTABLE: Found solution but recovery error > 10%")
    end
else
    println("  ❌ No solution found")
    println("    This might happen with a small grid (GN=6)")
    println("    Try increasing GN or degree_range for better results")
end

println("\n" * "="^80)
println("Test Complete!")
println("="^80)
println("\nKey features demonstrated:")
println("  ✓ create_globtim_objective() - Creates optimization-ready objective")
println("  ✓ run_globtim_optimization() - Runs global optimization")
println("  ✓ Real-time progress display with display_optimization_progress()")
println("  ✓ Automatic result collection and metrics computation")
println("  ✓ No file output (save_results=false)")
println("="^80)
