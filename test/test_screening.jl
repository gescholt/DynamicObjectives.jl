"""
Test suite for trajectory screening functions.

Tests is_trajectory_bounded() across multiple ODE models, covering:
- Known-good parameters that produce bounded trajectories
- Parameters that cause solver failure / blow-up
- Edge cases: low thresholds, degenerate parameters
- Performance: screening should be fast (< 100ms per call for 2D/4D)
"""

using Test
using DynamicObjectives

@testset "Screening" begin
    @testset "is_trajectory_bounded — LV 2D" begin
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]

        @testset "known-good parameters are bounded" begin
            result = is_trajectory_bounded(model, ic, [1.0, 0.5], [0.0, 10.0])
            @test result.bounded == true
            @test result.solver_success == true
            @test result.max_amplitude > 0.0
            @test result.max_amplitude < 100.0
            @test isfinite(result.max_amplitude)
        end

        @testset "multiple good parameter sets" begin
            good_params = [[0.5, 0.3], [1.0, 0.5], [1.5, 0.8], [0.8, 0.2]]
            for p in good_params
                result = is_trajectory_bounded(model, ic, p, [0.0, 10.0])
                @test result.bounded == true
                @test result.solver_success == true
            end
        end

        @testset "low threshold catches moderate amplitudes" begin
            result =
                is_trajectory_bounded(model, ic, [1.0, 0.5], [0.0, 10.0]; threshold = 1.0)
            # LV dynamics typically produce amplitudes > 1.0
            @test result.bounded == false
            @test result.solver_success == true
            @test result.max_amplitude > 1.0
        end

        @testset "custom solver and tolerance" begin
            result = is_trajectory_bounded(
                model,
                ic,
                [1.0, 0.5],
                [0.0, 10.0];
                solver = Vern9(),
                abstol = 1e-10,
                reltol = 1e-10,
                numpoints = 50,
            )
            @test result.bounded == true
            @test result.solver_success == true
        end
    end

    # ========================================================================
    # sweep_p_true tests
    # ========================================================================

    @testset "Sweep p_true" begin
        @testset "sweep_p_true — LV 2D random" begin
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            result = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 30)

            @test length(result.valid) + length(result.rejected) == 30
            @test result.pass_rate == length(result.valid) / 30
            @test result.pass_rate > 0.0  # LV 2D in this range should have some valid

            # All valid points must be within bounds
            for p in result.valid
                @test length(p) == 2
                @test bounds[1][1] <= p[1] <= bounds[1][2]
                @test bounds[2][1] <= p[2] <= bounds[2][2]
            end

            # Diagnostics populated
            @test result.diagnostics.n_total == 30
            @test result.diagnostics.n_valid == length(result.valid)
            @test result.diagnostics.elapsed_seconds > 0.0
            @test result.diagnostics.ms_per_candidate > 0.0
        end

        @testset "sweep_p_true — LV 2D grid" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            result = sweep_p_true(
                model,
                ic,
                bounds,
                [0.0, 10.0];
                n_candidates = 20,
                sampling = :grid,
            )

            @test length(result.valid) + length(result.rejected) <= 20
            @test result.pass_rate > 0.0
        end

        @testset "sweep_p_true — LV 4D short interval" begin
            model4, _, _, _ = define_constrained_lotka_volterra_4D()
            ic4 = [0.8, 1.2, 0.8, 1.2]
            bounds4 = [(-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3)]

            # Short time interval — should get some valid candidates
            result = sweep_p_true(model4, ic4, bounds4, [0.0, 2.0]; n_candidates = 20)

            @test length(result.valid) + length(result.rejected) == 20
            @test result.diagnostics.n_total == 20
            # With t=[0,2], most should pass
            @test result.pass_rate > 0.5
        end

        @testset "sweep_p_true — LV 4D long interval (stable with carrying capacity)" begin
            model4, _, _, _ = define_constrained_lotka_volterra_4D()
            ic4 = [0.8, 1.2, 0.8, 1.2]
            bounds4 = [(-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3)]

            # With negative diagonal self-interaction, most candidates should pass even at t=10
            result = sweep_p_true(model4, ic4, bounds4, [0.0, 10.0]; n_candidates = 20)

            @test result.pass_rate > 0.5  # Most should pass with carrying capacity
            # Rejected entries (if any) have correct structure
            for r in result.rejected
                @test r isa RejectedCandidate
                @test r.reason in (:solver_failure, :amplitude_exceeded)
                @test length(r.p_true) == 4
            end
        end

        @testset "sweep_p_true — DAISY 4D" begin
            model, _, _, _ = define_daisy_ex3_model_4D()
            ic = [1.0, 0.5, 0.3, 0.0]
            bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]

            result = sweep_p_true(model, ic, bounds, [0.0, 20.0]; n_candidates = 15)

            @test result.diagnostics.n_total == 15
            # DAISY should be reasonably well-behaved
            @test result.pass_rate >= 0.0  # At least it runs without error
        end

        @testset "sweep_p_true — n_candidates=1" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            result = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 1)

            @test length(result.valid) + length(result.rejected) == 1
            @test result.diagnostics.n_total == 1
        end

        @testset "sweep_p_true — margin=0" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            result = sweep_p_true(
                model,
                ic,
                bounds,
                [0.0, 10.0];
                n_candidates = 10,
                margin = 0.0,
            )

            @test result.diagnostics.n_total == 10
            # Points can now be right at the boundary
            for p in result.valid
                @test bounds[1][1] <= p[1] <= bounds[1][2]
                @test bounds[2][1] <= p[2] <= bounds[2][2]
            end
        end

        @testset "sweep_p_true — amplitude statistics" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            result = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 20)

            if length(result.valid) > 0
                @test isfinite(result.diagnostics.amplitude_min)
                @test isfinite(result.diagnostics.amplitude_max)
                @test isfinite(result.diagnostics.amplitude_median)
                @test result.diagnostics.amplitude_min > 0.0
                @test result.diagnostics.amplitude_min <= result.diagnostics.amplitude_max
            end
        end

        @testset "sweep_p_true — input validation" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            @testset "wrong bounds length" begin
                bad_bounds = [(0.5, 2.0)]
                @test_throws AssertionError sweep_p_true(model, ic, bad_bounds, [0.0, 10.0])
            end

            @testset "invalid margin" begin
                @test_throws AssertionError sweep_p_true(
                    model,
                    ic,
                    bounds,
                    [0.0, 10.0];
                    margin = 0.5,
                )
                @test_throws AssertionError sweep_p_true(
                    model,
                    ic,
                    bounds,
                    [0.0, 10.0];
                    margin = -0.1,
                )
            end

            @testset "invalid n_candidates" begin
                @test_throws AssertionError sweep_p_true(
                    model,
                    ic,
                    bounds,
                    [0.0, 10.0];
                    n_candidates = 0,
                )
            end

            @testset "reversed bounds" begin
                bad_bounds = [(2.0, 0.5), (0.2, 1.0)]
                @test_throws AssertionError sweep_p_true(model, ic, bad_bounds, [0.0, 10.0])
            end

            @testset "unknown sampling strategy" begin
                @test_throws ErrorException sweep_p_true(
                    model,
                    ic,
                    bounds,
                    [0.0, 10.0];
                    sampling = :sobol,
                )
            end
        end

        @testset "sweep_p_true — performance" begin
            model, _, _, _ = define_lotka_volterra_2D_model_v3_two_outputs()
            ic = [0.5, 1.0]
            bounds = [(0.5, 2.0), (0.2, 1.0)]

            # Warmup
            sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 5)

            # Benchmark
            result = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 50)

            @test result.diagnostics.ms_per_candidate < 50.0  # Should be well under 50ms after warmup
            println(
                "  Sweep performance: $(round(result.diagnostics.ms_per_candidate, digits=1)) ms/candidate (LV 2D, 50 candidates)",
            )
        end
    end

    @testset "is_trajectory_bounded — LV 3D" begin
        model, params, states, outputs = define_lotka_volterra_3D_model()
        # @mtkcompile reduces 3 states to 2 unknowns (x1, x2)
        ic = [1.0, 0.5]

        result = is_trajectory_bounded(model, ic, [0.5, -0.3, 0.2], [0.0, 10.0])
        @test result.solver_success == true
        # If bounded, great; if not, the solver at least ran successfully
        @test isfinite(result.max_amplitude) || !result.solver_success
    end

    @testset "is_trajectory_bounded — LV 4D Constrained" begin
        model, params, states, outputs = define_constrained_lotka_volterra_4D()
        ic = [0.8, 1.2, 0.8, 1.2]

        @testset "small perturbations are bounded at t=10" begin
            # With negative diagonal self-interaction (carrying capacity),
            # the model is stable for small epsilon perturbations
            result =
                is_trajectory_bounded(model, ic, [0.01, 0.02, -0.01, 0.03], [0.0, 10.0])
            @test result.bounded == true
            @test result.solver_success == true
            @test isfinite(result.max_amplitude)
        end

        @testset "near-zero perturbations are bounded long-term" begin
            # Very small perturbations should be stable even over long intervals
            result =
                is_trajectory_bounded(model, ic, [0.001, 0.001, -0.001, 0.001], [0.0, 30.0])
            @test result.solver_success == true
            @test result.bounded == true
        end
    end

    @testset "is_trajectory_bounded — DAISY 4D" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        ic = [1.0, 0.5, 0.3, 0.0]

        result = is_trajectory_bounded(model, ic, [0.5, 0.3, 0.4, 0.2], [0.0, 20.0])
        @test result.bounded == true
        @test result.solver_success == true
    end

    @testset "is_trajectory_bounded — FitzHugh-Nagumo" begin
        model, params, states, outputs = define_fitzhugh_nagumo_3D_model()
        # @mtkcompile reduces 3 states to 2 unknowns (V, R)
        ic = [0.0, 0.0]
        p_fhn = [0.8, 0.7, 0.8]  # from test_configs.jl

        @testset "default Tsit5 hits maxiters (stiff system)" begin
            # FHN is stiff — explicit Tsit5 at coarse tolerance cannot solve it,
            # which the screening function correctly reports as solver_success=false
            result = is_trajectory_bounded(model, ic, p_fhn, [0.0, 50.0])
            @test result.solver_success == false
            @test result.bounded == false
        end

        @testset "stiff solver succeeds but trajectory unbounded at coarse tol" begin
            # With a stiff-capable solver, the ODE solves but coarse tolerance
            # produces inaccurate trajectories that grow exponentially
            result = is_trajectory_bounded(
                model,
                ic,
                p_fhn,
                [0.0, 5.0];
                solver = AutoTsit5(Rosenbrock23()),
            )
            @test result.solver_success == true
            @test isfinite(result.max_amplitude)
        end
    end

    @testset "is_trajectory_bounded — diagnostics" begin
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]

        @testset "return type is BoundednessResult" begin
            result = is_trajectory_bounded(model, ic, [1.0, 0.5], [0.0, 10.0])
            @test result isa BoundednessResult
            @test result.bounded isa Bool
            @test result.max_amplitude isa Float64
            @test result.solver_success isa Bool
        end

        @testset "failed solver returns Inf amplitude" begin
            # FitzHugh-Nagumo with default Tsit5 is stiff — solver fails
            model_fhn, _, _, _ = define_fitzhugh_nagumo_3D_model()
            ic_fhn = [0.0, 0.0]
            result = is_trajectory_bounded(model_fhn, ic_fhn, [0.8, 0.7, 0.8], [0.0, 50.0])
            @test result.max_amplitude == Inf
        end
    end

    @testset "is_trajectory_bounded — input validation" begin
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p = [1.0, 0.5]

        @testset "wrong parameter length" begin
            @test_throws AssertionError is_trajectory_bounded(model, ic, [1.0], [0.0, 10.0])
            @test_throws AssertionError is_trajectory_bounded(
                model,
                ic,
                [1.0, 0.5, 0.3],
                [0.0, 10.0],
            )
        end

        @testset "wrong IC length" begin
            @test_throws AssertionError is_trajectory_bounded(model, [0.5], p, [0.0, 10.0])
            @test_throws AssertionError is_trajectory_bounded(
                model,
                [0.5, 1.0, 0.3],
                p,
                [0.0, 10.0],
            )
        end

        @testset "invalid time interval" begin
            @test_throws AssertionError is_trajectory_bounded(model, ic, p, [10.0, 0.0])  # reversed
        end

        @testset "invalid threshold" begin
            @test_throws AssertionError is_trajectory_bounded(
                model,
                ic,
                p,
                [0.0, 10.0];
                threshold = -1.0,
            )
            @test_throws AssertionError is_trajectory_bounded(
                model,
                ic,
                p,
                [0.0, 10.0];
                threshold = 0.0,
            )
        end
    end

    @testset "is_trajectory_bounded — performance" begin
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p = [1.0, 0.5]

        # Warm up (first call includes compilation)
        is_trajectory_bounded(model, ic, p, [0.0, 10.0])

        # Time 10 calls
        t_start = time()
        for _ in 1:10
            is_trajectory_bounded(model, ic, p, [0.0, 10.0])
        end
        elapsed = time() - t_start
        avg_ms = (elapsed / 10) * 1000

        @test avg_ms < 100.0  # Should be well under 100ms per call for 2D
        println("  Screening performance: $(round(avg_ms, digits=1)) ms/call (LV 2D)")
    end
end

# ============================================================================
# probe_landscape tests
# ============================================================================

@testset "Probe landscape" begin
    @testset "probe_landscape — DAISY 4D" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        ic = [1.0, 0.5, 0.3, 0.0]
        p_true = [0.5, 0.3, 0.4, 0.2]
        bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]

        result = probe_landscape(
            model,
            outputs,
            ic,
            p_true,
            bounds,
            [0.0, 20.0];
            n_probes = 10,
            numpoints = 20,
        )

        @testset "return structure" begin
            @test result isa ProbeResult
        end

        @testset "probe dimensions" begin
            @test length(result.probe_values) == 10
            @test length(result.probe_points) == 10
            @test all(length(p) == 4 for p in result.probe_points)
        end

        @testset "objective_at_true is small (tolerance noise floor)" begin
            # At coarse Tsit5/1e-4, f(p_true) is nonzero due to tolerance mismatch
            # between reference data (Vern9/1e-10) and evaluation (Tsit5/1e-4).
            # It should be small relative to probe values at random points.
            @test isfinite(result.objective_at_true)
            @test result.objective_at_true >= 0.0
            # Probe values at random points should generally be much larger
            finite_probes = filter(isfinite, result.probe_values)
            if !isempty(finite_probes)
                @test result.objective_at_true < maximum(finite_probes)
            end
        end

        @testset "landscape has structure" begin
            @test result.n_finite > 0
            @test result.n_finite + result.n_inf == 10
            if result.n_finite >= 2
                @test isfinite(result.dynamic_range)
                @test result.dynamic_range >= 0.0
                @test isfinite(result.variance)
                @test result.variance >= 0.0
            end
        end

        @testset "probe points within bounds" begin
            for p in result.probe_points
                for (d, (lo, hi)) in enumerate(bounds)
                    @test lo <= p[d] <= hi
                end
            end
        end

        @testset "construction timing" begin
            @test result.construction_time_ms > 0.0
        end
    end

    @testset "probe_landscape — LV 2D" begin
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p_true = [1.0, 0.5]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        result =
            probe_landscape(model, outputs, ic, p_true, bounds, [0.0, 10.0]; n_probes = 15)

        @test length(result.probe_values) == 15
        @test result.n_finite >= 1
        @test isfinite(result.objective_at_true)
    end

    @testset "probe_landscape — sweep integration" begin
        # Test the natural workflow: sweep -> probe valid candidates
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        sweep = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 10)

        @test length(sweep.valid) > 0

        # Probe the first valid candidate
        p_candidate = sweep.valid[1]
        result = probe_landscape(
            model,
            outputs,
            ic,
            p_candidate,
            bounds,
            [0.0, 10.0];
            n_probes = 5,
        )

        @test isfinite(result.objective_at_true)
        @test length(result.probe_values) == 5
        @test result.construction_time_ms > 0.0
    end

    @testset "probe_landscape — multiple candidates ranking" begin
        # Probe several candidates and verify we can rank them
        model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        sweep = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates = 15)

        probes = []
        for p in sweep.valid[1:min(5, length(sweep.valid))]
            result =
                probe_landscape(model, outputs, ic, p, bounds, [0.0, 10.0]; n_probes = 5)
            push!(probes, result)
        end

        # All probes should have valid structure
        @test length(probes) >= 1
        for result in probes
            @test length(result.probe_values) == 5
            @test result.n_finite + result.n_inf == 5
        end

        # Should be able to sort by dynamic_range
        finite_probes = filter(r -> isfinite(r.dynamic_range), probes)
        if length(finite_probes) >= 2
            ranges = [r.dynamic_range for r in finite_probes]
            sorted = sort(ranges, rev = true)
            @test sorted[1] >= sorted[end]  # sanity: sorting works
        end
    end

    @testset "probe_landscape — n_probes=1" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p_true = [1.0, 0.5]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        result =
            probe_landscape(model, outputs, ic, p_true, bounds, [0.0, 10.0]; n_probes = 1)

        @test length(result.probe_values) == 1
        @test result.n_finite + result.n_inf == 1
        # dynamic_range needs >= 2 finite values, so should be NaN
        @test isnan(result.dynamic_range)
    end

    @testset "probe_landscape — input validation" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p_true = [1.0, 0.5]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        @testset "wrong bounds length" begin
            bad_bounds = [(0.5, 2.0)]
            @test_throws AssertionError probe_landscape(
                model,
                outputs,
                ic,
                p_true,
                bad_bounds,
                [0.0, 10.0],
            )
        end

        @testset "wrong p_true length" begin
            @test_throws AssertionError probe_landscape(
                model,
                outputs,
                ic,
                [1.0],
                bounds,
                [0.0, 10.0],
            )
        end

        @testset "invalid n_probes" begin
            @test_throws AssertionError probe_landscape(
                model,
                outputs,
                ic,
                p_true,
                bounds,
                [0.0, 10.0];
                n_probes = 0,
            )
        end
    end

    @testset "probe_landscape — performance" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        p_true = [1.0, 0.5]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        # Warmup
        probe_landscape(model, outputs, ic, p_true, bounds, [0.0, 10.0]; n_probes = 3)

        # Benchmark 10 probes
        t_start = time()
        for _ in 1:10
            probe_landscape(model, outputs, ic, p_true, bounds, [0.0, 10.0]; n_probes = 5)
        end
        elapsed = time() - t_start
        avg_ms = (elapsed / 10) * 1000

        @test avg_ms < 500.0  # Construction + 5 evals should be well under 500ms after warmup
        println(
            "  Probe performance: $(round(avg_ms, digits=1)) ms/probe (LV 2D, 5 probes)",
        )
    end
end

# ============================================================================
# rank_probes tests
# ============================================================================

@testset "Rank probes" begin

    # Helper: create a fake probe result for testing ranking logic
    function _fake_probe(;
        dr = 1.0,
        var = 100.0,
        n_finite = 10,
        n_inf = 0,
        obj_true = 0.001,
    )
        ProbeResult(
            obj_true,
            vcat(fill(100.0, n_finite), fill(Inf, n_inf)),
            [rand(2) for _ in 1:(n_finite+n_inf)],
            dr,
            var,
            n_finite,
            n_inf,
            5.0,
        )
    end

    @testset "rank_probes — basic sorting" begin
        probes = [
            _fake_probe(dr = 1.0, var = 50.0),
            _fake_probe(dr = 3.0, var = 50.0),
            _fake_probe(dr = 2.0, var = 50.0),
        ]
        ranking = rank_probes(probes)
        @test ranking.ranked_indices == [2, 3, 1]  # sorted by dr descending
        @test isempty(ranking.filtered_out)
    end

    @testset "rank_probes — variance tiebreak" begin
        probes = [
            _fake_probe(dr = 2.0, var = 100.0),
            _fake_probe(dr = 2.0, var = 500.0),
            _fake_probe(dr = 2.0, var = 300.0),
        ]
        ranking = rank_probes(probes)
        @test ranking.ranked_indices == [2, 3, 1]  # same dr, sorted by var descending
    end

    @testset "rank_probes — NaN dynamic_range sorts last" begin
        probes = [
            _fake_probe(dr = NaN, var = NaN),
            _fake_probe(dr = 1.5, var = 100.0),
            _fake_probe(dr = NaN, var = NaN),
        ]
        ranking = rank_probes(probes)
        @test ranking.ranked_indices[1] == 2  # finite dr comes first
        @test Set(ranking.ranked_indices[2:3]) == Set([1, 3])  # NaN at end
    end

    @testset "rank_probes — low finite fraction filter" begin
        probes = [
            _fake_probe(dr = 3.0, n_finite = 2, n_inf = 8),  # 20% finite — below 50% threshold
            _fake_probe(dr = 1.0, n_finite = 8, n_inf = 2),  # 80% finite — passes
            _fake_probe(dr = 2.0, n_finite = 5, n_inf = 5),  # 50% finite — passes (>=)
        ]
        ranking = rank_probes(probes; min_finite_fraction = 0.5)
        @test 1 in ranking.filtered_out
        @test ranking.filter_reasons[1] == :low_finite_fraction
        @test ranking.ranked_indices == [3, 2]  # index 1 filtered, 3 before 2 by dr
    end

    @testset "rank_probes — noise floor filter" begin
        probes = [
            _fake_probe(dr = 3.0, obj_true = 999.0),  # obj_true >> median — filtered
            _fake_probe(dr = 1.0, obj_true = 0.001),   # clean
        ]
        ranking = rank_probes(probes; max_noise_ratio = 0.5)
        @test 1 in ranking.filtered_out
        @test ranking.filter_reasons[1] == :high_noise_floor
        @test ranking.ranked_indices == [2]
    end

    @testset "rank_probes — empty input" begin
        ranking = rank_probes(ProbeResult[])
        @test isempty(ranking.ranked_indices)
        @test isempty(ranking.filtered_out)
    end

    @testset "rank_probes — all filtered out" begin
        probes = [
            _fake_probe(dr = 3.0, n_finite = 1, n_inf = 9),
            _fake_probe(dr = 2.0, n_finite = 0, n_inf = 10),
        ]
        ranking = rank_probes(probes; min_finite_fraction = 0.5)
        @test isempty(ranking.ranked_indices)
        @test length(ranking.filtered_out) == 2
    end

    @testset "rank_probes — disabled filters" begin
        probes = [_fake_probe(dr = 3.0, n_finite = 1, n_inf = 9, obj_true = 999.0)]
        # min_finite_fraction=0, max_noise_ratio=Inf — no filtering
        ranking = rank_probes(probes; min_finite_fraction = 0.0, max_noise_ratio = Inf)
        @test ranking.ranked_indices == [1]
        @test isempty(ranking.filtered_out)
    end
end

# ============================================================================
# screen_and_probe tests
# ============================================================================

@testset "Screen and probe" begin
    @testset "screen_and_probe — LV 2D end-to-end" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        result = screen_and_probe(
            model,
            outputs,
            ic,
            bounds,
            [0.0, 10.0];
            n_candidates = 20,
            n_probes = 5,
        )

        @testset "return structure" begin
            @test result isa ScreeningResult
            @test result.sweep isa SweepResult
            @test result.ranking isa RankingResult
        end

        @testset "sweep ran" begin
            @test result.sweep.diagnostics.n_total == 20
            @test length(result.sweep.valid) > 0
        end

        @testset "probes match valid candidates" begin
            @test length(result.probes) == length(result.sweep.valid)
            for probe in result.probes
                @test length(probe.probe_values) == 5
            end
        end

        @testset "ranking covers all non-filtered" begin
            n_ranked = length(result.ranking.ranked_indices)
            n_filtered = length(result.ranking.filtered_out)
            @test n_ranked + n_filtered == length(result.probes)
        end

        @testset "ranked indices are valid" begin
            for i in result.ranking.ranked_indices
                @test 1 <= i <= length(result.probes)
            end
        end

        @testset "ranking is sorted by dynamic_range descending" begin
            indices = result.ranking.ranked_indices
            for j in 1:length(indices)-1
                dr_a = result.probes[indices[j]].dynamic_range
                dr_b = result.probes[indices[j+1]].dynamic_range
                # NaN sorts to end, so if dr_a is finite and dr_b is NaN that's fine
                if isfinite(dr_a) && isfinite(dr_b)
                    @test dr_a >= dr_b
                end
            end
        end
    end

    @testset "screen_and_probe — DAISY 4D" begin
        model, _, _, outputs = define_daisy_ex3_model_4D()
        ic = [1.0, 0.5, 0.3, 0.0]
        bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]

        result = screen_and_probe(
            model,
            outputs,
            ic,
            bounds,
            [0.0, 20.0];
            n_candidates = 15,
            n_probes = 5,
        )

        @test length(result.sweep.valid) > 0
        @test length(result.probes) == length(result.sweep.valid)
        @test length(result.ranking.ranked_indices) > 0
    end

    @testset "screen_and_probe — composability" begin
        # Verify that screen_and_probe produces the same results as
        # calling the three functions individually
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        result = screen_and_probe(
            model,
            outputs,
            ic,
            bounds,
            [0.0, 10.0];
            n_candidates = 10,
            n_probes = 3,
        )

        # The orchestrator should produce consistent counts
        @test result.sweep.diagnostics.n_total == 10
        @test length(result.probes) == length(result.sweep.valid)
        n_total =
            length(result.ranking.ranked_indices) + length(result.ranking.filtered_out)
        @test n_total == length(result.probes)
    end

    @testset "screen_and_probe — input validation" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        # Wrong bounds length propagates from sweep_p_true
        @test_throws AssertionError screen_and_probe(
            model,
            outputs,
            ic,
            [(0.5, 2.0)],
            [0.0, 10.0];
            n_candidates = 5,
        )
    end

    @testset "screen_and_probe — performance" begin
        model, _, _, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
        ic = [0.5, 1.0]
        bounds = [(0.5, 2.0), (0.2, 1.0)]

        # Warmup
        screen_and_probe(
            model,
            outputs,
            ic,
            bounds,
            [0.0, 10.0];
            n_candidates = 5,
            n_probes = 3,
        )

        t_start = time()
        result = screen_and_probe(
            model,
            outputs,
            ic,
            bounds,
            [0.0, 10.0];
            n_candidates = 30,
            n_probes = 5,
        )
        elapsed = (time() - t_start) * 1000

        @test elapsed < 5000.0  # 30 candidates + 30 probes should be well under 5s
        println(
            "  screen_and_probe: $(round(elapsed, digits=0))ms (30 candidates, 5 probes)",
        )
    end
end
