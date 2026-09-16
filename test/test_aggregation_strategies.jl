using Test
using DynamicObjectives
using DynamicObjectives: AGGREGATION_STRATEGIES, resolve_aggregation
import DynamicObjectives: Tsit5

# Tests for the named aggregation registry.
# The `aggregate_distances` knob on make_error_distance is already a
# positional arg accepting any callable — this suite asserts the named
# registry is present and that distinct strategies produce distinguishable
# scalar objectives on a fixed fixture.

@testset "Aggregation strategies" begin
    @testset "Registry contents" begin
        @test AGGREGATION_STRATEGIES[:sum] === sum
        @test AGGREGATION_STRATEGIES[:maximum] === maximum
        @test AGGREGATION_STRATEGIES[:minimum] === minimum
        @test AGGREGATION_STRATEGIES[:first] === first
        @test haskey(AGGREGATION_STRATEGIES, :mean)
        @test haskey(AGGREGATION_STRATEGIES, :rms)
    end

    @testset "resolve_aggregation" begin
        @test resolve_aggregation(:sum) === sum
        @test resolve_aggregation(:maximum) === maximum
        @test resolve_aggregation(:rms)(Float64[3.0, 4.0]) ≈ sqrt(12.5)
        @test_throws ArgumentError resolve_aggregation(:nonexistent_strategy)
    end

    @testset "Strategies produce distinguishable objectives (LV4D constrained)" begin
        model, params, states, outputs = define_constrained_lotka_volterra_4D()
        # Constrained LV4D parameters are epsilon perturbations around a
        # fixed-point regime; zeros give a stable reference trajectory.
        p_true = [0.0, 0.0, 0.0, 0.0]
        ic = [1.0, 2.0, 1.0, 2.0]
        time_interval = [0.0, 10.0]

        # Keep CI fast: Tsit5 at loose tolerance with 10 sample points.
        build(strategy) = make_error_distance(
            model,
            outputs,
            ic,
            p_true,
            time_interval,
            10,
            L2_norm,
            strategy,
            nothing;
            return_inf_on_error = false,
            solver = Tsit5(),
            abstol = 1e-6,
            reltol = 1e-6,
        )

        err_sum = build(AGGREGATION_STRATEGIES[:sum])
        err_max = build(AGGREGATION_STRATEGIES[:maximum])
        err_mean = build(AGGREGATION_STRATEGIES[:mean])
        err_first = build(AGGREGATION_STRATEGIES[:first])

        # Different per-output sensitivities → aggregation must differ
        p_test = p_true .+ [0.1, 0.2, 0.3, 0.4]
        v_sum = err_sum(p_test)
        v_max = err_max(p_test)
        v_mean = err_mean(p_test)
        v_first = err_first(p_test)

        @test all(isfinite, (v_sum, v_max, v_mean, v_first))
        # At least three distinct scalar values across the four strategies.
        @test length(unique((v_sum, v_max, v_mean, v_first))) >= 3

        # Max is >= mean >= any single component, all non-negative.
        @test v_max >= v_mean
        @test v_mean >= 0
        @test v_first >= 0
    end
end
