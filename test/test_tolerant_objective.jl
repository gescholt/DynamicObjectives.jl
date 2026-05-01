# Tests for TolerantObjective — mutable solver/tolerance wrapper

@testset "TolerantObjective" begin

    # Use a small fast model for all tests
    model, params, states, outputs = define_lotka_volterra_2D_model_v3()
    p_true = [1.0, 0.5]
    ic = [1.0, 0.5]
    time_interval = [0.0, 20.0]
    numpoints = 15

    # ── Construction ─────────────────────────────────────────────────────
    @testset "Construction" begin
        @testset "Default solver/tolerance" begin
            tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints)
            @test get_tolerance(tol_obj) == (1e-10, 1e-10)
            @test get_solver(tol_obj) isa Vern9
        end

        @testset "Custom solver/tolerance" begin
            tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
                solver=Tsit5(), abstol=1e-4, reltol=1e-6)
            abstol, reltol = get_tolerance(tol_obj)
            @test abstol == 1e-4
            @test reltol == 1e-6
            @test get_solver(tol_obj) isa Tsit5
        end

        @testset "With non-default distance and aggregation" begin
            tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints,
                L1_norm, maximum)
            # Should construct without error and be callable
            @test tol_obj(p_true) isa Float64
        end
    end

    # ── Callable interface ───────────────────────────────────────────────
    @testset "Callable Interface" begin
        tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints)

        @testset "Returns Float64" begin
            val = tol_obj(p_true)
            @test val isa Float64
        end

        @testset "At true params ~0 (tight tolerance)" begin
            val = tol_obj(p_true)
            @test val < 1e-6
        end

        @testset "Perturbed params > 0" begin
            p_perturbed = p_true .+ 0.1
            val = tol_obj(p_perturbed)
            @test val > 0.0
        end
    end

    # ── set_tolerance! ───────────────────────────────────────────────────
    @testset "set_tolerance!" begin
        tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
            solver=Vern9(), abstol=1e-10, reltol=1e-10)

        @testset "Single-arg form" begin
            set_tolerance!(tol_obj, 1e-4)
            abstol, reltol = get_tolerance(tol_obj)
            @test abstol == 1e-4
            @test reltol == 1e-4
        end

        @testset "Two-arg form" begin
            set_tolerance!(tol_obj, 1e-6, 1e-8)
            abstol, reltol = get_tolerance(tol_obj)
            @test abstol == 1e-6
            @test reltol == 1e-8
        end

        @testset "No-op when unchanged" begin
            set_tolerance!(tol_obj, 1e-6, 1e-8)  # already at these values
            # Should return quickly without rebuilding
            abstol, reltol = get_tolerance(tol_obj)
            @test abstol == 1e-6
            @test reltol == 1e-8
        end

        @testset "Objective still callable after switch" begin
            set_tolerance!(tol_obj, 1e-4)
            val = tol_obj(p_true)
            @test val isa Float64
            @test isfinite(val)
        end
    end

    # ── set_solver! ──────────────────────────────────────────────────────
    @testset "set_solver!" begin
        tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints)
        @test get_solver(tol_obj) isa Vern9

        set_solver!(tol_obj, Tsit5())
        @test get_solver(tol_obj) isa Tsit5

        # Still callable
        val = tol_obj(p_true)
        @test val isa Float64
        @test isfinite(val)

        # Switch back
        set_solver!(tol_obj, Vern9())
        @test get_solver(tol_obj) isa Vern9
    end

    # ── Tolerance switching preserves correctness ────────────────────────
    @testset "Coarse vs Tight Tolerance" begin
        # Build two objectives: one tight, one coarse
        tight_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
            solver=Vern9(), abstol=1e-10, reltol=1e-10)
        coarse_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
            solver=Tsit5(), abstol=1e-4, reltol=1e-4)

        # At true params: tight should be ~0, coarse has noise floor
        tight_val = tight_obj(p_true)
        coarse_val = coarse_obj(p_true)
        @test tight_val < 1e-6
        # Coarse tolerance means the reference and evaluation data differ slightly
        # so objective at p_true is nonzero (typically 1-10 range)
        @test coarse_val >= 0.0  # non-negative by construction

        # At perturbed params: both should be positive
        p_perturbed = [1.2, 0.6]
        tight_perturbed = tight_obj(p_perturbed)
        coarse_perturbed = coarse_obj(p_perturbed)
        @test tight_perturbed > 0.0
        @test coarse_perturbed > 0.0

        # Ordering should be preserved: same perturbed point should be "worse" than true
        # (for tight tolerance)
        @test tight_perturbed > tight_val
    end

    @testset "Switch preserves tight accuracy" begin
        # Start coarse, switch to tight — should match a fresh tight construction
        switchable = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
            solver=Tsit5(), abstol=1e-4, reltol=1e-4)

        # Evaluate at coarse
        coarse_val = switchable(p_true)
        @test isfinite(coarse_val)

        # Switch to tight
        set_tolerance!(switchable, 1e-10)
        set_solver!(switchable, Vern9())

        # Now should match tight accuracy
        tight_val = switchable(p_true)
        @test tight_val < 1e-6

        # Compare with a fresh tight objective
        fresh_tight = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
            solver=Vern9(), abstol=1e-10, reltol=1e-10)
        fresh_val = fresh_tight(p_true)

        # They should agree closely (both use same reference data from construction)
        @test abs(tight_val - fresh_val) < 1e-12
    end

    # ── Re-exported solvers ──────────────────────────────────────────────
    @testset "Re-exported Solvers" begin
        @test Tsit5() isa Tsit5
        @test Vern9() isa Vern9
    end

end
