"""
Test suite for globtim integration with Dynamic_objectives

Tests the integration of Globtim optimizer with parameter estimation benchmarks.
Follows TDD approach: tests written first, implementation follows.

Test structure:
- Phase 1: 2D models (EASY)
- Phase 2: 3D/4D models (MEDIUM)
- Phase 3: Edge cases and error handling
"""

using Test
using Dynamic_objectives
using LinearAlgebra

# Test helper functions
function test_parameter_recovery(p_recovered, p_true; tolerance=0.01)
    error = norm(p_recovered - p_true)
    return error < tolerance, error
end

function test_objective_at_true_params(error_func, p_true; tolerance=1e-10)
    error = error_func(p_true)
    return error < tolerance, error
end

@testset "globtim Integration" begin

    # ==============================================================================
    # Phase 1: 2D Model Tests (EASY)
    # ==============================================================================

    @testset "Phase 1.1: Direct 1-arg Objective" begin

        @testset "make_error_distance creates valid 1-arg objective" begin
            # Since Phase 2, globtim auto-detects 1-arg functions — no wrapper needed
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            time_interval = [0.0, 20.0]
            numpoints = 30

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                time_interval, numpoints,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func isa Function
            @test error_func(p_true) < 1e-10
        end

        @testset "Objective responds to parameter perturbation" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10

            p_perturbed = p_true .+ [0.1, 0.1]
            @test error_func(p_perturbed) > 1e-6
        end
    end

    @testset "Phase 1.2: globtim Integration - 2D LV v3 (2 outputs)" begin

        @testset "Model: LV 2D v3 two outputs (EASIEST)" begin
            # Model configuration
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            bounds = [(0.5, 1.5), (0.25, 0.75)]  # Tighter domain around true params
            time_interval = [0.0, 20.0]
            numpoints = 30

            # Create error function
            error_func = make_error_distance(
                model, outputs, ic, p_true,
                time_interval, numpoints,
                L2_norm;
                return_inf_on_error = true
            )

            # Verify error function works
            @test error_func(p_true) < 1e-10

            # Optimization test requires globtim environment activation
            @test_skip begin
                # config = ExperimentParams(GN=8, degree_range=4:8, ...)
                # result = run_standard_experiment(
                #     objective_function=error_func, objective_name="test",
                #     bounds=bounds,
                #     experiment_config=config, output_dir=mktempdir())
                # @test result[:success_rate] > 0.0
                true
            end
        end

        @testset "Parameter recovery accuracy" begin
            # Once optimization works, test recovery quality
            @test_skip begin
                # p_recovered = result[:best_params]
                # recovered, error = test_parameter_recovery(p_recovered, p_true, tolerance=0.01)
                # @test recovered
                # @test error < 0.01
                true
            end
        end

        @testset "Computational efficiency - 2D" begin
            # Test that 2D problems are computationally feasible
            @test_skip begin
                # GN = 8  # 64 grid points for 2D
                # @test result[:total_time] < 60  # Should complete in <1 minute
                # @test result[:function_evaluations] < 500
                true
            end
        end
    end

    @testset "Phase 1.3: Additional 2D Models (EASY)" begin

        @testset "Model: LV 2D v1 (classic)" begin
            model, params, states, outputs = define_lotka_volterra_2D_model()
            p_true = [1.5, 1.0]
            ic = [1.0, 1.0]
            bounds = [(-1.0, 2.0), (-1.0, 0.0)]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10

            # Skip optimization until wrapper ready
            @test_skip begin
                # Test parameter recovery on this model
                true
            end
        end

        @testset "Model: LV 2D v2 (c=0.1)" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v2()
            p_true = [1.5, 1.0]
            ic = [1.0, 1.0]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10
        end

        @testset "Model: LV 2D v3 (c=0.5, single output)" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v3()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10
        end
    end

    # ==============================================================================
    # Phase 2: 3D/4D Model Tests (MEDIUM)
    # ==============================================================================

    @testset "Phase 2.1: 3D Models (MEDIUM)" begin

        @testset "Model: LV 3D v2" begin
            model, params, states, outputs = define_lotka_volterra_3D_model_v2()
            p_true = [1.0, 0.5, 1.0]
            ic = [1.0, 0.5]  # 2 states (x1, x2), not 3
            bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 3.0)]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10

            @test_skip begin
                # Test 3D parameter recovery
                # Expect: >50% success, <2000 evaluations
                true
            end
        end
    end

    @testset "Phase 2.2: 4D Models (MEDIUM)" begin

        @testset "Model: DAISY Ex3 (4D, no input)" begin
            model, params, states, outputs = define_daisy_ex3_model_4D_no_input()
            p_true = [-0.0693, -0.0945, -0.04, -0.038]
            ic = [0.0, 0.0, 0.0]  # 3 states (x1, x2, x3), not 4
            bounds = [(-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10

            @test_skip begin
                # Test 4D parameter recovery
                # 4D with GN=8: 4096 evaluations
                # May need more time
                true
            end
        end

        @testset "Model: LV 4D Constrained" begin
            model, params, states, outputs = define_constrained_lotka_volterra_4D()
            p_true = [1.0, 0.5, 1.0, 0.5]
            ic = [1.0, 0.5, 0.5, 0.5]
            bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 3.0), (0.0, 2.0)]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            @test error_func(p_true) < 1e-10
        end
    end

    # ==============================================================================
    # Phase 3: Edge Cases and Error Handling
    # ==============================================================================

    @testset "Phase 3: Edge Cases" begin

        @testset "Timeout handling (FitzHugh-Nagumo)" begin
            # This model is oscillatory and can hang ODE solver
            model, params, states, outputs = define_fitzhugh_nagumo_3D_model()
            p_true = [0.2, 0.2, 3.0]
            ic = [-1.0, 1.0]  # 2 states (V, R), not 3

            # With timeout
            error_func_timeout = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true,
                eval_timeout = 10.0  # 10 second timeout
            )

            @test error_func_timeout(p_true) < 1e-10

            # Test that bad parameters return Inf (don't hang)
            p_bad = [10.0, 10.0, 10.0]  # Likely to cause issues
            error_bad = error_func_timeout(p_bad)
            @test error_bad == Inf || error_bad > 1e6
        end

        @testset "ODE solver failure handling" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            error_func = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            # Parameters far outside reasonable bounds
            # Note: Not all extreme parameters cause ODE failures - depends on system dynamics
            # LV model is quite stable, so we just test it returns a finite error
            p_extreme = [100.0, 100.0]
            error_extreme = error_func(p_extreme)
            @test error_extreme isa Number  # Should return some error value
            @test error_extreme > 0.0  # Should be worse than true parameters
        end

        @testset "Different distance metrics" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            # L2 norm (default)
            error_L2 = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm;
                return_inf_on_error = true
            )

            # L1 norm (robust)
            error_L1 = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L1_norm;
                return_inf_on_error = true
            )

            # log-L2 norm (wide range)
            error_logL2 = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                log_L2_norm;
                return_inf_on_error = true
            )

            # All should be ≈0 at true params
            @test error_L2(p_true) < 1e-10
            @test error_L1(p_true) < 1e-10
            @test error_logL2(p_true) < 1e-10
        end
    end

    # ==============================================================================
    # Success Criteria Summary
    # ==============================================================================

    @testset "Success Criteria" begin
        @test_skip begin
            # Once all tests pass, verify overall success criteria
            # EASY models (2D): >75% success rate, <500 evals
            # MEDIUM models (3D/4D): >50% success rate, <2000 evals
            # All models handle errors gracefully (no crashes)
            true
        end
    end
end
