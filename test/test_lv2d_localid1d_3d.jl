"""
    test_lv2d_localid1d_3d.jl

Tests for the hybrid LV2D + LocalID-1D 3D objective function.

Verifies:
1. Global minima at (1.0, 1.5, +/-1.5) have F ~ 0
2. Non-minimum points have F > 0
3. FiniteDiff gradient at global minima ~ 0
4. Symmetry-breaking: Hessians differ at +1.5 vs -1.5

Run with:
    julia --project=profiles/dev -e 'include("pkg/Dynamic_objectives/test/test_lv2d_localid1d_3d.jl")'
"""

using Test
using Dynamic_objectives
using FiniteDiff: finite_difference_gradient, finite_difference_hessian

@testset "LV2D + LocalID-1D 3D Hybrid Objective" begin

    @testset "Construction" begin
        hybrid = create_lv2d_localid1d_3d_objective()
        @test hybrid isa Function

        # alpha must be < 1
        @test_throws ErrorException create_lv2d_localid1d_3d_objective(alpha=1.0)
        @test_throws ErrorException create_lv2d_localid1d_3d_objective(alpha=1.5)
    end

    @testset "Dimension check" begin
        hybrid = create_lv2d_localid1d_3d_objective()
        @test_throws ErrorException hybrid([1.0, 2.0])       # 2D
        @test_throws ErrorException hybrid([1.0, 2.0, 3.0, 4.0])  # 4D
    end

    hybrid = create_lv2d_localid1d_3d_objective(alpha=0.3, beta=0.068)

    @testset "Global minima" begin
        # Both (1.0, 1.5, +1.5) and (1.0, 1.5, -1.5) should give F ~ 0
        f_pos = hybrid([1.0, 1.5, 1.5])
        f_neg = hybrid([1.0, 1.5, -1.5])
        @test f_pos < 1e-6
        @test f_neg < 1e-6
    end

    @testset "Non-minimum points are positive" begin
        @test hybrid([2.0, 3.0, 0.0]) > 0.0
        @test hybrid([5.0, 5.0, 1.0]) > 0.0
        @test hybrid([1.0, 1.5, 0.0]) > 0.0  # true 2D but wrong 1D
        @test hybrid([0.5, 0.5, 1.5]) > 0.0  # true 1D but wrong 2D
    end

    @testset "Gradient at global minima" begin
        grad_pos = finite_difference_gradient(hybrid, [1.0, 1.5, 1.5])
        grad_neg = finite_difference_gradient(hybrid, [1.0, 1.5, -1.5])
        @test maximum(abs, grad_pos) < 1e-2
        @test maximum(abs, grad_neg) < 1e-2
    end

    @testset "Symmetry-breaking via modulation" begin
        # The sin(pi*a3) modulation makes the basins at +1.5 and -1.5 different
        # At a3 = +1.5: sin(3pi/2) = -1, scale = (1 - alpha) = 0.7
        # At a3 = -1.5: sin(-3pi/2) = +1, scale = (1 + alpha) = 1.3
        # So the Hessians should differ in the a1/a2 block
        H_pos = finite_difference_hessian(hybrid, [1.0, 1.5, 1.5])
        H_neg = finite_difference_hessian(hybrid, [1.0, 1.5, -1.5])
        # The 2D Hessian block (a1, a2) should be scaled differently
        @test !isapprox(H_pos[1:2, 1:2], H_neg[1:2, 1:2], rtol=0.1)
    end

    @testset "Non-separability" begin
        # If the objective were separable, F(a1, a2, a3) = g(a1,a2) + h(a3),
        # then the mixed partial d^2F / da1 da3 would be zero everywhere.
        # With modulation, it's nonzero at points where E_LV2D > 0.
        # Use a3=0.25 where cos(pi*0.25) = cos(pi/4) ~ 0.707 (nonzero)
        # At a3=0.5, cos(pi/2)=0 so mixed partials vanish — avoid that
        H_off = finite_difference_hessian(hybrid, [2.0, 3.0, 0.25])
        # Off-diagonal block (a1/a2 vs a3) should be nonzero
        @test abs(H_off[1, 3]) > 1e-8 || abs(H_off[2, 3]) > 1e-8
    end
end

println("\nAll LV2D + LocalID-1D 3D hybrid tests passed!")
