#!/usr/bin/env julia
"""
Test: Constrained Lotka-Volterra 4D

A 4-parameter model with skew-symmetric interaction structure.
This model tests global optimization on a moderately difficult nonlinear system.

Run:
    julia --project=../.. test_lv_4d_constrained.jl
"""

include("model_test_utils.jl")

# Define model configuration
config = ModelTestConfig(
    name = "LV_4D_Constrained",
    description = "4D Lotka-Volterra with skew-symmetric perturbations",

    model_fn = define_constrained_lotka_volterra_4D,

    p_true = [1.0, 0.5, 1.0, 0.5],
    ic = [1.0, 0.5, 0.5, 0.5],
    bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 3.0), (0.0, 2.0)],

    time_interval = [0.0, 20.0],
    numpoints = 120,

    # Optimization settings
    GN = 10,
    degree_range = 6:10,
    max_time = 1800.0,
    max_time_per_point = 30.0,

    distance = L2_norm
)

# Run the test
result = run_model_test(config)

println("\nTest completed successfully!")
println("Results saved to: $(result.output_dir)")
