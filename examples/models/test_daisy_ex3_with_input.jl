#!/usr/bin/env julia
"""
Test: DAISY Ex3 with Time-Dependent Input

A 4-parameter DAISY benchmark model with time-dependent input.
This is a linear compartmental model from the DAISY database.

Run:
    julia --project=../.. test_daisy_ex3_with_input.jl
"""

include("model_test_utils.jl")

# Define model configuration
config = ModelTestConfig(
    name = "DAISY_Ex3_with_input",
    description = "DAISY benchmark with time-dependent input",

    model_fn = define_daisy_ex3_model_4D,

    p_true = [0.5, 0.3, 0.4, 0.2],
    ic = [1.0, 0.5, 0.3, 0.0],
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
