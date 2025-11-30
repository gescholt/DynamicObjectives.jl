#!/usr/bin/env julia
"""
Test: DAISY Ex3 without Input

A 4-parameter DAISY benchmark model without time-dependent input.
Simplified variant of the DAISY Ex3 compartmental model.

Run:
    julia --project=../.. test_daisy_ex3_no_input.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("model_test_framework.jl")

# Define model configuration
config = ModelTestConfig(
    name = "DAISY_Ex3_no_input",
    description = "DAISY benchmark without input",

    model_fn = define_daisy_ex3_model_4D_no_input,

    p_true = [0.5, 0.3, 0.4, 0.2],
    ic = [1.0, 0.5, 0.3],
    bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)],

    time_interval = [0.0, 20.0],
    numpoints = 100,

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
