using Test
using DynamicObjectives
import DynamicObjectives: Tsit5

# Tests for the partial-observability helper.
# The underlying capability (passing a subset of outputs to
# make_error_distance) already worked — this suite asserts the helper
# validates indices and that restricting outputs changes the landscape.

@testset "Partial observability" begin
    model, params, states, outputs = define_constrained_lotka_volterra_4D()
    @assert length(outputs) == 4  # fixture precondition

    p_true = [0.0, 0.0, 0.0, 0.0]
    ic = [1.0, 2.0, 1.0, 2.0]
    time_interval = [0.0, 10.0]

    solver_kwargs =
        (return_inf_on_error = false, solver = Tsit5(), abstol = 1e-6, reltol = 1e-6)

    err_full = make_error_distance(
        model,
        outputs,
        ic,
        p_true,
        time_interval,
        10,
        L2_norm,
        sum,
        nothing;
        solver_kwargs...,
    )
    err_partial = make_partial_observability_distance(
        model,
        outputs,
        [1, 3],
        ic,
        p_true,
        time_interval,
        10,
        L2_norm,
        sum,
        nothing;
        solver_kwargs...,
    )

    @testset "Both match reference at p_true" begin
        @test err_full(p_true) < 1e-4
        @test err_partial(p_true) < 1e-4
    end

    @testset "Restricting outputs changes the landscape" begin
        p_test = p_true .+ [0.1, 0.2, 0.3, 0.4]
        v_full = err_full(p_test)
        v_partial = err_partial(p_test)

        @test isfinite(v_full)
        @test isfinite(v_partial)
        @test v_full > 0
        @test v_partial > 0
        # Distinct landscapes: a 2-output objective should not evaluate
        # identically to a 4-output objective at a perturbed p_test.
        @test v_full != v_partial
        # With L2_norm + sum, full-output error upper-bounds partial.
        @test v_full >= v_partial
    end

    @testset "Index validation" begin
        @test_throws ArgumentError make_partial_observability_distance(
            model,
            outputs,
            Int[],
            ic,
            p_true,
            time_interval,
        )
        @test_throws ArgumentError make_partial_observability_distance(
            model,
            outputs,
            [0, 1],
            ic,
            p_true,
            time_interval,
        )
        @test_throws ArgumentError make_partial_observability_distance(
            model,
            outputs,
            [1, 5],
            ic,
            p_true,
            time_interval,
        )
    end
end
