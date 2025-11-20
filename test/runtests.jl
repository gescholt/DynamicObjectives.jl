using Dynamic_objectives
using Test

@testset "Dynamic_objectives.jl" begin

    # Core functionality tests
    @testset "DAISY LV4D Definition" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        @test length(params) == 4
        @test length(states) == 4
        @test length(outputs) == 2
    end

    @testset "Generalized LV4D Definition" begin
        model, params, states, outputs = define_generalized_lotka_volterra_4D()
        @test length(params) == 20  # 4 growth rates + 16 interaction matrix elements
        @test length(states) == 4
        @test length(outputs) == 4
    end

    @testset "Constrained LV4D Definition" begin
        model, params, states, outputs = define_constrained_lotka_volterra_4D()
        @test length(params) == 4  # Only epsilon perturbations
        @test length(states) == 4
        @test length(outputs) == 4
    end

    @testset "Error Function Creation" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        p_true = [0.1, 0.2, 0.3, 0.4]
        ic = [1.0, 2.0, 1.0, 1.0]
        time_interval = [0.0, 10.0]

        error_func = make_error_distance(
            model, outputs, ic, p_true, time_interval, 25
        )

        # Test at true parameters (should be ~0)
        @test error_func(p_true) < 1e-6

        # Test at perturbed parameters (should be > 0)
        p_test = p_true .+ 0.1
        @test error_func(p_test) > 0
    end

    @testset "Distance Functions" begin
        y_true = [1.0, 2.0, 3.0]
        y_test = [1.1, 2.1, 3.1]

        @test L1_norm(y_true, y_test) > 0
        @test L2_norm(y_true, y_test) > 0
        @test log_L2_norm(y_true, y_test) < 0  # log of small number
    end

end

# globtim integration tests (comprehensive test suite)
include("test_globtim_integration.jl")
