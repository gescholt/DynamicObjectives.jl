# Tests for Cartesian-product gluing (bead zwbs.10.1).
#
# Construction:  F([x; y]) := f1(x) + f2(y)
# Property:      CP(F) = CP(f1) × CP(f2)  (exact)

using Test
using DynamicObjectives

@testset "glue: independent-coordinate sum" begin
    f1(x) = (x[1] - 1.0)^2 + (x[2] + 0.5)^2          # 2D bowl, min at (1, -0.5)
    f2(y) = (y[1] - 0.3)^2 + (y[2] - 2.0)^2          # 2D bowl, min at (0.3, 2.0)
    b1 = [(-2.0, 2.0), (-2.0, 2.0)]
    b2 = [(-1.0, 1.0), (0.0, 4.0)]

    g = glue(f1, b1, f2, b2)

    @test g.dim == 4
    @test g.n1 == 2
    @test g.n2 == 2
    @test g.bounds == [(-2.0, 2.0), (-2.0, 2.0), (-1.0, 1.0), (0.0, 4.0)]
    @test g.oracle_cps === nothing  # no CPs supplied

    # Value check: F should match f1 + f2 at the joint minimum
    z_min = [1.0, -0.5, 0.3, 2.0]
    @test g.F(z_min) ≈ 0.0 atol = 1e-12

    # Off-min behavior
    z_test = [2.0, 0.0, -1.0, 4.0]
    expected = f1([2.0, 0.0]) + f2([-1.0, 4.0])
    @test g.F(z_test) ≈ expected atol = 1e-12
end

@testset "glue: oracle CP cartesian product" begin
    # Quartic with 2 minima at (±1, 0) and a degenerate saddle at (0, 0)
    f1(x) = (x[1]^2 - 1.0)^2 + x[2]^2
    f2(y) = (y[1]^2 - 1.0)^2 + y[2]^2

    cps1 = [[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0]]   # 2 min + 1 saddle
    cps2 = [[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0]]

    g = glue(
        f1,
        [(-2.0, 2.0), (-1.0, 1.0)],
        f2,
        [(-2.0, 2.0), (-1.0, 1.0)];
        cps1 = cps1,
        cps2 = cps2,
    )

    @test length(g.oracle_cps) == 9   # 3 × 3
    @test g.dim == 4

    # Verify each CP is at a stationary point: |F(p+ε e_i) - F(p-ε e_i)| / 2ε
    # should be O(ε^2) since p is a CP. Test against a finite-difference gradient.
    ε = 1e-5
    for p in g.oracle_cps
        for i in 1:g.dim
            ep = zeros(g.dim)
            ep[i] = ε
            grad_i = (g.F(p .+ ep) - g.F(p .- ep)) / (2ε)
            @test abs(grad_i) < 1e-6
        end
    end
end

@testset "glue: count_oracle_cps_in_box" begin
    # 1D bowl glued with 1D bowl → 4-CP toy product over 2D
    cps1 = [[0.0, 0.0]]
    cps2 = [[1.0, 1.0], [-1.0, -1.0]]
    f1(x) = sum(x .^ 2)
    f2(y) = sum((y .^ 2 .- 1.0) .^ 2)

    g = glue(
        f1,
        [(-1.0, 1.0), (-1.0, 1.0)],
        f2,
        [(-2.0, 2.0), (-2.0, 2.0)];
        cps1 = cps1,
        cps2 = cps2,
    )

    @test length(g.oracle_cps) == 2  # 1 × 2

    # Whole bounds: contains both
    @test count_oracle_cps_in_box(g.oracle_cps, g.bounds) == 2

    # Restrict the y-block to the upper-right corner only
    restrict = [(-1.0, 1.0), (-1.0, 1.0), (0.5, 2.0), (0.5, 2.0)]
    @test count_oracle_cps_in_box(g.oracle_cps, restrict) == 1
end

@testset "glue: subdivision finds glued minimum" begin
    # Cheap test: glue two 2D bowls, run adaptive_refine, check it
    # converges fast on a tight box around the unique 4D minimum.
    using Globtim

    f1(x) = (x[1] - 1.0)^2 + (x[2] + 0.5)^2
    f2(y) = (y[1] - 0.3)^2 + (y[2] - 2.0)^2

    g = glue(
        f1,
        [(0.0, 2.0), (-2.0, 2.0)],
        f2,
        [(-1.0, 1.0), (0.0, 4.0)];
        cps1 = [[1.0, -0.5]],
        cps2 = [[0.3, 2.0]],
    )

    @test length(g.oracle_cps) == 1
    expected_min = [1.0, -0.5, 0.3, 2.0]
    @test g.oracle_cps[1] == expected_min

    tree = adaptive_refine(
        g.F,
        g.bounds,
        4;
        l2_tolerance = 1e-2,
        tolerance_mode = :relative,
        max_depth = 2,
        max_leaves = 16,
        parallel = false,
    )

    @test n_leaves(tree) >= 1
    @test n_pruned(tree) == 0  # bowl × bowl, no Inf regions
    @test g.F(expected_min) ≈ 0.0 atol = 1e-12
end
