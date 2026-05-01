# Tests for gradient compatibility of ODE objectives
#
# Verifies that ODE-based objective functions from make_error_distance / create_objective
# produce correct gradients via both FiniteDiff and ForwardDiff.
#
# ForwardDiff compatibility was enabled by replacing the in-place mutation
# `problem.p.tunable .= p_true` with `SciMLBase.remake(problem; p=...)` in
# data_generation.jl, which preserves ForwardDiff.Dual types through the ODE solve.

using ForwardDiff
using FiniteDiff

@testset "Gradient Compatibility" begin

    # Use a small fast 2D model
    model, params, states, outputs = define_lotka_volterra_2D_model_v3()
    p_true = [1.0, 0.5]
    ic = [1.0, 0.5]
    time_interval = [0.0, 20.0]
    numpoints = 15

    objective = make_error_distance(model, outputs, ic, p_true, time_interval, numpoints)

    # ── Objective sanity checks ──────────────────────────────────────────

    @testset "Objective returns near-zero at p_true" begin
        @test objective(p_true) < 1e-6
    end

    @testset "Objective returns positive away from p_true" begin
        p_away = [1.5, 0.3]
        @test objective(p_away) > 0.01
    end

    # ── Finite-difference gradient ───────────────────────────────────────

    @testset "Finite-difference gradient is non-zero at non-stationary point" begin
        # Pick a point clearly away from any minimum
        p_test = [1.5, 0.3]
        f0 = objective(p_test)

        # Manual central finite differences
        ε = 1e-6
        grad = zeros(2)
        for i in 1:2
            p_plus  = copy(p_test); p_plus[i]  += ε
            p_minus = copy(p_test); p_minus[i] -= ε
            grad[i] = (objective(p_plus) - objective(p_minus)) / (2ε)
        end

        grad_norm = sqrt(sum(grad .^ 2))

        # The gradient must be non-zero at this non-stationary point
        @test grad_norm > 1e-4
        # Each component should be finite
        @test all(isfinite, grad)
    end

    @testset "Gradient near p_true is small" begin
        # At the global minimum, gradient should be approximately zero
        ε = 1e-6
        grad = zeros(2)
        for i in 1:2
            p_plus  = copy(p_true); p_plus[i]  += ε
            p_minus = copy(p_true); p_minus[i] -= ε
            grad[i] = (objective(p_plus) - objective(p_minus)) / (2ε)
        end

        grad_norm = sqrt(sum(grad .^ 2))

        # Should be very small (numerical noise level)
        @test grad_norm < 1e-2
    end

    # ── ForwardDiff gradient ─────────────────────────────────────────────

    @testset "ForwardDiff gradient is non-zero at non-stationary point" begin
        p_test = [1.5, 0.3]
        grad_ad = ForwardDiff.gradient(objective, p_test)

        @test all(isfinite, grad_ad)
        @test !all(grad_ad .== 0)  # Must NOT be all-zero (the old bug)
        @test sqrt(sum(grad_ad .^ 2)) > 1e-4
    end

    @testset "ForwardDiff gradient near p_true" begin
        # Note: at exactly p_true, the L2_norm objective is ||Y_true - Y_test||₂ = sqrt(v'v)
        # where v → 0. The gradient of sqrt(v'v) at v=0 is undefined (subgradient),
        # so ForwardDiff may return a non-small gradient at the exact minimum.
        # This is mathematically correct — the function is non-smooth at its minimum.
        # We test slightly offset from p_true to avoid this non-smoothness.
        p_near = p_true .+ 1e-4
        grad_ad = ForwardDiff.gradient(objective, p_near)

        @test all(isfinite, grad_ad)
        # Near (not at) the minimum, the ForwardDiff gradient should be moderate
        @test sqrt(sum(grad_ad .^ 2)) < 1e3
    end

    @testset "ForwardDiff matches FiniteDiff gradient" begin
        p_test = [1.1, 0.6]
        grad_ad = ForwardDiff.gradient(objective, p_test)
        grad_fd = FiniteDiff.finite_difference_gradient(objective, p_test)

        # ForwardDiff should agree with FiniteDiff to ~6 digits
        # (FiniteDiff is O(h^2) accurate, typically ~1e-8 relative error)
        for i in 1:2
            rel_err = abs(grad_ad[i] - grad_fd[i]) / max(abs(grad_ad[i]), 1e-10)
            @test rel_err < 1e-4
        end
    end

    # ── ForwardDiff Hessian ──────────────────────────────────────────────

    @testset "ForwardDiff Hessian is non-zero at non-stationary point" begin
        p_test = [1.1, 0.6]
        H_ad = ForwardDiff.hessian(objective, p_test)

        @test all(isfinite, H_ad)
        @test !all(H_ad .== 0)
        # Hessian should be symmetric
        @test H_ad ≈ H_ad' atol=1e-10
    end

    # ── create_objective pathway ─────────────────────────────────────────

    @testset "create_objective gradient is non-zero at non-stationary point" begin
        entry = CatalogueEntry(
            name = "test_grad",
            description = "gradient test",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = p_true,
            ic = ic,
            bounds = [(0.5, 2.0), (0.1, 1.0)],
            time_interval = time_interval,
            numpoints = numpoints,
        )
        obj = create_objective(entry)

        p_test = [1.5, 0.3]

        ε = 1e-6
        grad = zeros(2)
        for i in 1:2
            p_plus  = copy(p_test); p_plus[i]  += ε
            p_minus = copy(p_test); p_minus[i] -= ε
            grad[i] = (obj(p_plus) - obj(p_minus)) / (2ε)
        end

        grad_norm = sqrt(sum(grad .^ 2))
        @test grad_norm > 1e-4
        @test all(isfinite, grad)
    end

    @testset "create_objective ForwardDiff gradient works" begin
        entry = CatalogueEntry(
            name = "test_grad_ad",
            description = "ForwardDiff gradient test",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = p_true,
            ic = ic,
            bounds = [(0.5, 2.0), (0.1, 1.0)],
            time_interval = time_interval,
            numpoints = numpoints,
        )
        obj = create_objective(entry)

        p_test = [1.5, 0.3]
        grad_ad = ForwardDiff.gradient(obj, p_test)

        @test all(isfinite, grad_ad)
        @test !all(grad_ad .== 0)
        @test sqrt(sum(grad_ad .^ 2)) > 1e-4
    end

    # ── Objective value consistency ──────────────────────────────────────

    @testset "Objective is deterministic (same input → same output)" begin
        p_test = [1.3, 0.7]
        v1 = objective(p_test)
        v2 = objective(p_test)
        @test v1 == v2
    end

    # ── log_L2_norm ForwardDiff ──────────────────────────────────────────

    @testset "ForwardDiff works with log_L2_norm objective" begin
        obj_log = make_error_distance(model, outputs, ic, p_true, time_interval, numpoints,
                                       log_L2_norm)
        p_test = [1.1, 0.6]
        grad_ad = ForwardDiff.gradient(obj_log, p_test)

        @test all(isfinite, grad_ad)
        @test !all(grad_ad .== 0)
    end
end
