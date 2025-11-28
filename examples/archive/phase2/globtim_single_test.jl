#!/usr/bin/env julia
"""
Example: Single model test with globtim optimizer

Tests the simplest 2D model (Lotka-Volterra v3 with 2 outputs) to verify
the globtim integration works correctly.

This script demonstrates:
1. Creating a globtim-compatible objective from Dynamic_objectives
2. Running optimization with StandardExperiment
3. Validating parameter recovery

Usage:
    julia --project=. examples/globtim_single_test.jl

Expected output:
- Parameter recovery error < 0.01
- Objective value at solution < 1e-6
- Completion in < 60 seconds (2D with GN=8)
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dynamic_objectives
using Globtim
using LinearAlgebra
using Printf

println("="^80)
println("globtim Integration Test - Single Model (2D LV v3)")
println("="^80)

# ==============================================================================
# Model Configuration: EASIEST model (2D, full observability)
# ==============================================================================

println("\n[1/5] Defining model...")

# Define Lotka-Volterra 2D v3 with 2 outputs (fully observable)
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# True parameters to recover
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Time series configuration
time_interval = [0.0, 20.0]
numpoints = 30

println("  Model: Lotka-Volterra 2D v3 (2 outputs)")
println("  True parameters: $p_true")
println("  Parameter bounds: $bounds")
println("  Time interval: $time_interval")
println("  Sample points: $numpoints")

# ==============================================================================
# Create globtim Objective
# ==============================================================================

println("\n[2/5] Creating objective function...")

# Use the integration wrapper
objective = create_globtim_objective(
    model, outputs, ic, p_true,
    time_interval, numpoints;
    distance = L2_norm,
    eval_timeout = nothing  # No timeout needed for this easy model
)

println("  ✓ Objective function created")

# Verify objective evaluates correctly
error_at_true = objective(p_true, nothing)
@printf "  Error at true params: %.6e (should be ≈0)\n" error_at_true

if error_at_true > 1e-10
    @warn "Unexpected error at true parameters! Check ODE integration."
end

# ==============================================================================
# Run globtim Optimization
# ==============================================================================

println("\n[3/5] Running globtim optimization...")
println("  Configuration:")
println("    Grid points (GN): 8")
println("    Grid size: 8^2 = 64 evaluations")
println("    Polynomial degrees: 4:8")
println("    Basis: Chebyshev")

# Create output directory
output_dir = mkpath(joinpath(@__DIR__, "..", "test_results", "lv2d_single_test"))

# Run optimization using the helper function
result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 8,
    degree_range = 4:8,
    output_dir = output_dir,
    model_name = "LV_2D_v3_two_outputs",
    max_time = 3600
)

println("  ✓ Optimization completed")

# ==============================================================================
# Analyze Results
# ==============================================================================

println("\n[4/5] Analyzing results...")

@printf "  Total time: %.2f seconds\n" result[:total_time]
@printf "  Success rate: %.1f%% (%d/%d degrees)\n" (result[:success_rate] * 100) result[:total_time] length(result[:degree_results])
println("  Total critical points found: $(result[:total_critical_points])")

# Find best solution
if result[:best_params] !== nothing
    println("\n  Best solution:")
    @printf "    Degree: %d\n" result[:best_degree]
    @printf "    Objective value: %.6e\n" result[:best_objective]
    println("    Best parameters: $(result[:best_params])")
    println("    True parameters: $p_true")
    @printf "    Recovery error: %.6e\n" result[:recovery_error]
else
    println("\n  ❌ No solution found!")
end

# ==============================================================================
# Validate Success Criteria
# ==============================================================================

println("\n[5/5] Validation...")

success = true

# Criterion 1: At least one degree succeeded
if result[:success_rate] == 0.0
    println("  ❌ FAIL: No degrees produced valid solutions")
    success = false
else
    println("  ✓ PASS: At least one degree succeeded")
end

# Criterion 2: Parameter recovery accuracy
if result[:best_params] !== nothing
    if result[:recovery_error] < 0.01
        println("  ✓ PASS: Parameter recovery error < 0.01")
    elseif result[:recovery_error] < 0.1
        println("  ⚠  PARTIAL: Recovery error < 0.1 (acceptable, could improve)")
    else
        println("  ❌ FAIL: Recovery error > 0.1")
        success = false
    end
else
    println("  ❌ FAIL: No solution to validate")
    success = false
end

# Criterion 3: Objective value at solution
if result[:best_objective] !== nothing && result[:best_objective] < 1e-6
    println("  ✓ PASS: Objective value < 1e-6 (near-perfect fit)")
elseif result[:best_objective] !== nothing && result[:best_objective] < 1e-3
    println("  ⚠  PARTIAL: Objective value < 1e-3 (good fit)")
elseif result[:best_objective] !== nothing
    println("  ❌ FAIL: Objective value > 1e-3 (poor fit)")
    success = false
end

# Criterion 4: Computational efficiency (2D should be fast)
if result[:total_time] < 60
    println("  ✓ PASS: Completed in < 60 seconds")
elseif result[:total_time] < 300
    println("  ⚠  PARTIAL: Completed in < 5 minutes (acceptable)")
else
    println("  ⚠  SLOW: Took > 5 minutes (unexpected for 2D)")
end

# ==============================================================================
# Summary
# ==============================================================================

println("\n" * "="^80)
if success
    println("✅ SUCCESS: Integration test passed!")
    println("\nFiles created:")
    println("  - $output_dir/results_summary.json")
    println("  - $output_dir/critical_points_deg_*.csv")
else
    println("❌ FAILED: Integration test did not meet criteria")
    println("\nDebugging:")
    println("  - Check output directory: $output_dir")
    println("  - Review results_summary.json for details")
    println("  - Verify globtimcore configuration is appropriate")
end
println("="^80)
