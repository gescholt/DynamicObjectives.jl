#!/usr/bin/env julia
# Validate all 13 models in Dynamic_objectives
# Run: julia --project=. test/validate_all_models.jl

using Dynamic_objectives
using Test

println("="^80)
println("Dynamic_objectives - Model Validation Suite")
println("="^80)

# Define all model tests with basic configurations
models = [
    # Easy 2D models
    (name = "LV 2D v1", fn = define_lotka_volterra_2D_model,
     p_true = [0.5, -0.3], ic = [1.0, 0.5]),
    (name = "LV 2D v2", fn = define_lotka_volterra_2D_model_v2,
     p_true = [0.5, -0.3], ic = [1.0, 0.5]),
    (name = "LV 2D v3", fn = define_lotka_volterra_2D_model_v3,
     p_true = [1.0, 0.5], ic = [1.0, 0.5]),
    (name = "LV 2D v3 (2 outputs)", fn = define_lotka_volterra_2D_model_v3_two_outputs,
     p_true = [1.0, 0.5], ic = [1.0, 0.5]),

    # Medium 3D models
    (name = "LV 3D v1", fn = define_lotka_volterra_3D_model,
     p_true = [0.5, -0.3, 0.2], ic = [1.0, 0.5]),
    (name = "LV 3D v2", fn = define_lotka_volterra_3D_model_v2,
     p_true = [1.0, 0.5, 0.3], ic = [1.0, 0.5]),

    # Medium 4D models
    (name = "DAISY Ex3 (with input)", fn = define_daisy_ex3_model_4D,
     p_true = [0.5, 0.3, 0.4, 0.2], ic = [1.0, 0.5, 0.3, 0.0]),
    (name = "DAISY Ex3 (no input)", fn = define_daisy_ex3_model_4D_no_input,
     p_true = [0.5, 0.3, 0.4, 0.2], ic = [1.0, 0.5, 0.3]),
    (name = "LV 4D Constrained", fn = define_constrained_lotka_volterra_4D,
     p_true = [0.01, 0.02, -0.01, 0.03], ic = [0.8, 1.2, 0.8, 1.2]),

    # Hard models
    (name = "FitzHugh-Nagumo", fn = define_fitzhugh_nagumo_3D_model,
     p_true = [0.8, 0.7, 0.8], ic = [0.0, 0.0]),
    (name = "Simple 1D (a²)", fn = define_simple_1D_model_locally_identifiable,
     p_true = [1.5], ic = [1.0]),
    (name = "Simple 2D (a*b)", fn = define_simple_2D_model_locally_identifiable,
     p_true = [2.0, 3.0], ic = [1.0]),
    (name = "Simple 2D (b²)", fn = define_simple_2D_model_locally_identifiable_square,
     p_true = [0.5, 2.0], ic = [1.0]),
]

# Test basic configuration (without 20D model - too slow for quick validation)
passed = 0
failed = 0

for (i, spec) in enumerate(models)
    print("[$i/13] Testing $(spec.name)... ")
    try
        # Define model
        model, params, states, outputs = spec.fn()

        # Create objective function
        time_interval = [0.0, 10.0]
        numpoints = 20
        error_func = make_error_distance(
            model, outputs, spec.ic, spec.p_true,
            time_interval, numpoints,
            L2_norm, first, nothing;
            return_inf_on_error = true,
            eval_timeout = 10.0
        )

        # Test at true parameters
        error_at_true = error_func(spec.p_true)
        @test error_at_true < 1e-6

        # Test at perturbed parameters
        p_test = spec.p_true .+ 0.1
        error_at_test = error_func(p_test)
        @test error_at_test > 0

        println("✓ PASS (error @ true: $(round(error_at_true, sigdigits=3)), @ perturbed: $(round(error_at_test, sigdigits=3)))")
        passed += 1
    catch e
        println("✗ FAIL")
        println("  Error: ", sprint(showerror, e))
        failed += 1
    end
end

# Test 20D model separately (may be slow)
println("\n" * "="^80)
print("[14/14] Testing LV 4D Generalized (20 params)... ")
try
    model, params, states, outputs = define_generalized_lotka_volterra_4D()
    p_true = [
        1.0, 1.0, 1.0, 1.0,  # Growth rates
        -0.5, -0.1, -0.1, -0.1,
        -0.1, -0.5, -0.1, -0.1,
        -0.1, -0.1, -0.5, -0.1,
        -0.1, -0.1, -0.1, -0.5
    ]
    ic = [0.5, 0.5, 0.5, 0.5]

    error_func = make_error_distance(
        model, outputs, ic, p_true,
        [0.0, 10.0], 20,
        L2_norm, first, nothing;
        return_inf_on_error = true,
        eval_timeout = 15.0  # Longer timeout for 20D
    )

    error_at_true = error_func(p_true)
    @test error_at_true < 1e-6

    println("✓ PASS (error @ true: $(round(error_at_true, sigdigits=3)))")
    passed += 1
catch e
    println("✗ FAIL")
    println("  Error: ", sprint(showerror, e))
    failed += 1
end

# Summary
println("="^80)
println("Validation Summary:")
println("  PASSED: $passed/14")
println("  FAILED: $failed/14")
println("="^80)

if failed > 0
    println("\n⚠️  Some models failed validation. Check errors above.")
    exit(1)
else
    println("\n✅ All 14 models validated successfully!")
    println("\nReady for globtim testing campaign!")
    println("See TESTING_GUIDE.md for recommended test configurations.")
    exit(0)
end
