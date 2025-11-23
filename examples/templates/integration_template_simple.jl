#!/usr/bin/env julia
"""
Integration Template: Simple Wrapper Pattern

This template demonstrates the SIMPLEST way to integrate Dynamic_objectives
with globtimcore + globtimpostprocessing for critical point finding.

USE WHEN:
- You want a quick test or prototype
- You don't need fine control over experiment parameters
- You're okay with sensible defaults

PATTERN:
- Uses high-level wrapper from Dynamic_objectives
- No manual config structs needed
- Clean keyword arguments
- Built-in progress monitoring

ARCHITECTURE:
┌─────────────────────────────────────────────────────────────┐
│ Dynamic_objectives wrapper                                  │
│   └─> globtimcore (find critical points)                   │
│   └─> globtimpostprocessing (refine with BFGS)             │
└─────────────────────────────────────────────────────────────┘

NOTE: This wrapper uses a simplified grid-based approach. For full polynomial
approximation + HomotopyContinuation, use the advanced or pipeline templates.
"""

using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))  # Activate Dynamic_objectives environment

using Dynamic_objectives
using LinearAlgebra

#===============================================================================
STEP 1: Define your model and create objective function
===============================================================================#

# Example: 2D Lotka-Volterra model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
p_true = [1.0, 0.5]  # True parameters to recover
ic = [1.0, 0.5]      # Initial condition
bounds = [(0.0, 3.0), (0.0, 2.0)]  # Parameter search bounds

# Create error function (1-arg function: p -> error)
objective = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,  # Time span and number of points
    L2_norm, first, nothing;
    eval_timeout = 10.0
)

#===============================================================================
STEP 2: Run optimization with simple wrapper
===============================================================================#

result = run_globtim_optimization(
    objective,
    bounds,
    p_true;
    GN = 50,                    # Grid size (50 points per dimension)
    degree_range = 4:8,         # Polynomial degrees to try
    output_dir = joinpath(@__DIR__, "..", "..", "test_results", "template_simple"),
    model_name = "lv2d_example",
    show_progress = true,
    save_results = true
)

#===============================================================================
STEP 3: Check results
===============================================================================#

println("\n" * "="^80)
println("Results Summary")
println("="^80)
println("Best parameters found: $(result[:best_params])")
println("True parameters:       $p_true")
println("Recovery error:        $(round(100*result[:recovery_error], digits=2))%")
println("Best objective value:  $(round(result[:best_value], digits=6))")
println("Total computation time: $(round(result[:total_time], digits=2))s")

if result[:recovery_error] < 0.05
    println("\n✅ Excellent recovery (< 5% error)")
elseif result[:recovery_error] < 0.10
    println("\n✓ Good recovery (< 10% error)")
else
    println("\n⚠️  Poor recovery - consider increasing GN or degree_range")
end

#===============================================================================
CUSTOMIZATION OPTIONS
===============================================================================#

# You can customize the wrapper with these keyword arguments:
#
# GN::Int = 8
#   Grid size (number of points per dimension)
#   Larger = more accurate but slower
#   Typical values: 8-20 for quick tests, 50-100 for production
#
# degree_range = 4:8
#   Range of polynomial degrees to try
#   Higher degrees can capture more complex landscapes
#   But also more expensive
#
# output_dir::String = "."
#   Where to save results
#
# model_name::String = "model"
#   Name for output files
#
# max_time::Real = 3600
#   Maximum time in seconds (default: 1 hour)
#
# basis::Symbol = :chebyshev
#   Polynomial basis (:chebyshev or :legendre)
#
# show_progress::Bool = true
#   Show progress messages
#
# save_results::Bool = true
#   Save CSV files with critical points
