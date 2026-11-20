# The k-brick kappa=0 product family. Its value is that the
# critical set of the 2k-parameter misfit is the product of k planar censuses,
# which only holds if brick i of the product IS the planar brick with the same
# rates. These tests pin that: trajectories of the product blocks equal the
# planar brick's, k = 2 equals the hand-written coupled model at kappa = 0, and
# the registered k = 3 model has the dimensions the catalogue entries assume.

using Test
using DynamicObjectives
using DynamicObjectives: Tsit5, ODEProblem
using DynamicObjectives.ModelingToolkit: unknowns, parameters
using DynamicObjectives.SciMLBase: solve

function _outputs(fn, p, ic; T = 20.0, n = 60)
    model, params, states, outs = fn()
    prob = ODEProblem(model, merge(Dict(unknowns(model) .=> ic), Dict(params .=> p)), (0.0, T))
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12, saveat = range(0.0, T; length = n))
    return reduce(hcat, [collect(sol[eq.rhs]) for eq in outs])
end

@testset "lv2d product family (5618.1)" begin
    brick = DynamicObjectives.define_lotka_volterra_2D_sciml_benchmark
    YA = _outputs(brick, [1.0, 1.5], ones(2))
    YB = _outputs(brick, [1.2, 1.3], ones(2))
    YC = _outputs(brick, [1.4, 1.1], ones(2))

    @testset "k = 1 is the brick" begin
        Y1 = _outputs(() -> DynamicObjectives._define_lv2d_sciml_product(1), [1.0, 1.5], ones(2))
        @test maximum(abs.(Y1 .- YA)) < 1e-9
    end

    @testset "k = 2 is the coupled model at kappa = 0" begin
        p4 = [1.0, 1.5, 1.2, 1.3]
        Y2 = _outputs(() -> DynamicObjectives._define_lv2d_sciml_product(2), p4, ones(4))
        Yc = _outputs(DynamicObjectives.define_lv2d_sciml_coupled_k000, p4, ones(4))
        @test maximum(abs.(Y2 .- Yc)) < 1e-9
        model, params, states, outs = DynamicObjectives._define_lv2d_sciml_product(2)
        @test string.(params) == ["a1", "a2", "a3", "a4"]
    end

    @testset "k = 3 is registered, six-dimensional, block-diagonal" begin
        @test haskey(DynamicObjectives.MODEL_REGISTRY, "define_lv2d_sciml_product3_k000")
        model, params, states, outs = DynamicObjectives.define_lv2d_sciml_product3_k000()
        @test length(params) == 6 && length(states) == 6 && length(outs) == 6
        @test length(parameters(model)) == 6 && length(unknowns(model)) == 6
        Y3 = _outputs(DynamicObjectives.define_lv2d_sciml_product3_k000,
            [1.0, 1.5, 1.2, 1.3, 1.4, 1.1], ones(6))
        @test maximum(abs.(Y3[:, 1:2] .- YA)) < 1e-9
        @test maximum(abs.(Y3[:, 3:4] .- YB)) < 1e-9
        @test maximum(abs.(Y3[:, 5:6] .- YC)) < 1e-9
        # Changing brick C's rates leaves bricks A and B untouched (no coupling).
        Y3b = _outputs(DynamicObjectives.define_lv2d_sciml_product3_k000,
            [1.0, 1.5, 1.2, 1.3, 0.9, 1.9], ones(6))
        @test maximum(abs.(Y3b[:, 1:4] .- Y3[:, 1:4])) < 1e-9
        @test maximum(abs.(Y3b[:, 5:6] .- YC)) > 1e-3
    end

    @test_throws ArgumentError DynamicObjectives._define_lv2d_sciml_product(0)
end
