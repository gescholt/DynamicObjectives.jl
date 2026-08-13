# check_test_deps: skip — runtests.jl includes this file only when
# HomotopyContinuation is resolvable (`Base.identify_package(...) !== nothing`),
# so its deps are intentionally absent from the standalone mirror's test env.
# Recovery regression for Cartesian-product gluing (bead zwbs.10.1).
#
# The other glue tests check that the analytic oracle is stationary and that
# adaptive_refine runs; NONE of them run globtim's critical-point solver and
# verify it RECOVERS the full oracle set. This closes that gap with a cheap
# 2D case: glue two 1D double-wells → a 2D objective with exactly 9 known CPs
# (4 min, 4 saddle, 1 max), fit at exact degree, enumerate CPs, and assert
# recall = precision = 1.0 broken down by Morse index.
#
# Cheap by construction: 2D, degree 4 ⇒ 5×5 grid, HC gradient Bezout 3²=9,
# sub-second. Guarded on HomotopyContinuation availability (workspace runs
# have it; the standalone mirror does not) — same pattern as test_ext_toml_pipeline.jl.

using Test
using DynamicObjectives
using Globtim
import HomotopyContinuation              # triggers the HC solver extension
using DynamicPolynomials
using DataFrames
using LinearAlgebra

@testset "glue: globtim recovers the full oracle CP set" begin
    # 1D double-well w(t) = (t²-1)²: CPs at t=0 (max) and t=±1 (min).
    w(v) = (v[1]^2 - 1.0)^2
    cps1d = [[-1.0], [0.0], [1.0]]
    idx1d = [0, 1, 0]                    # min, max(1D index 1), min

    g = glue(w, [(-2.0, 2.0)], w, [(-2.0, 2.0)]; cps1 = cps1d, cps2 = cps1d)
    @test g.dim == 2
    @test length(g.oracle_cps) == 9

    # aligned oracle points + Morse indices (same nested order glue() uses)
    oracle_pts = Vector{Float64}[]
    oracle_idx = Int[]
    for (p, ip) in zip(cps1d, idx1d), (q, iq) in zip(cps1d, idx1d)
        push!(oracle_pts, vcat(p, q))
        push!(oracle_idx, ip + iq)
    end
    @test count(==(0), oracle_idx) == 4      # minima
    @test count(==(1), oracle_idx) == 4      # saddles
    @test count(==(2), oracle_idx) == 1      # maximum

    # exact-degree fit (double-well is degree 4 ⇒ p.nrm ≈ 0)
    T = TestInput(g.F; dim = 2, center = [0.0, 0.0], sample_range = [2.0, 2.0],
                  GN = 4, tolerance = nothing)
    p = Constructor(T, 4)
    @test p.nrm < 1e-8                       # fit is exact ⇒ we test the solver

    @polyvar x[1:2]
    sols = solve_polynomial_system(x, p; basis = p.basis, normalized = p.normalized)
    df = process_crit_pts(sols, g.F, T)
    found = [[df.x1[i], df.x2[i]] for i in 1:nrow(df)]
    @test !isempty(found)

    # match recovered ↔ oracle (greedy NN, physical coords)
    morse_index(H; tol = 1e-8) = count(λ -> λ < -tol, eigvals(Symmetric(H)))
    tau = 0.05
    used = falses(length(oracle_pts))
    matched_idx = Int[]                      # oracle indices that got matched
    n_spurious = 0
    for a in found
        best, bj = tau, 0
        for (j, b) in enumerate(oracle_pts)
            used[j] && continue
            d = norm(a - b)
            (d < best) && ((best, bj) = (d, j))
        end
        if bj != 0
            used[bj] = true
            push!(matched_idx, oracle_idx[bj])
        else
            n_spurious += 1
        end
    end

    # recall = 1.0 (every oracle CP recovered), precision = 1.0 (no spurious)
    @test length(matched_idx) == 9           # all 9 oracle CPs matched
    @test n_spurious == 0
    @test count(==(0), matched_idx) == 4     # all 4 minima recovered
    @test count(==(1), matched_idx) == 4     # all 4 saddles recovered
    @test count(==(2), matched_idx) == 1     # the maximum recovered

    # recovered points are numerically on the oracle
    for a in found
        @test minimum(norm(a - b) for b in oracle_pts) < 1e-6
    end
end
