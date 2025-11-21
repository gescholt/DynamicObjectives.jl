using Dynamic_objectives
using Test

@testset "Dynamic_objectives.jl" begin
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

    @testset "Display Infrastructure" begin
        # Test DisplayConfig creation
        config = DisplayConfig()
        @test config.use_color == true
        @test config.table_backend == :unicode
        @test config.precision == 4

        # Test custom config
        custom_config = DisplayConfig(use_color=false, precision=6)
        @test custom_config.use_color == false
        @test custom_config.precision == 6

        # Test display functions don't error (output is visual, so we just check they run)
        model, params, states, outputs = define_daisy_ex3_model_4D()

        # Suppress output during tests
        original_stdout = stdout
        redirect_stdout(devnull)

        try
            # Test model summary
            @test_nowarn display_model_summary(model, config=DisplayConfig(use_color=false))

            # Test parameter display
            param_names = [:α, :β, :γ, :δ]
            p_values = [0.1, 0.2, 0.3, 0.4]
            @test_nowarn display_parameters(param_names, p_values, config=DisplayConfig(use_color=false))

            # Test error metrics display
            errors = Dict("L1" => 0.123, "L2" => 0.045)
            @test_nowarn display_error_metrics(errors, config=DisplayConfig(use_color=false))

            # Test optimization result display
            result = (params=p_values, error=0.001, iterations=10, converged=true, time_elapsed=1.0)
            @test_nowarn display_optimization_result(result, config=DisplayConfig(use_color=false))

            # Test time series display
            ic = [1.0, 2.0, 1.0, 1.0]
            time_interval = [0.0, 10.0]
            p_true = [0.1, 0.2, 0.3, 0.4]

            using ModelingToolkit: complete
            problem = ModelingToolkit.ODEProblem(
                complete(model),
                merge(
                    Dict(states .=> ic),
                    Dict(params .=> p_true)
                ),
                time_interval
            )

            data = sample_data(problem, model, outputs, time_interval, p_true, ic, 10)
            @test_nowarn display_time_series(data, show_plot=false, config=DisplayConfig(use_color=false))
        finally
            redirect_stdout(original_stdout)
        end
    end
end
