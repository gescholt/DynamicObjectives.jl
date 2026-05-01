"""
Test suite for grid-based interestingness scoring.

Tests use synthetic grid data (no ODE solving) to verify:
- Local minima detection on structured grids
- Dynamic range computation
- Basin fraction computation
- Deceptive minima detection
- Composite interestingness score
- GridScoreResult construction
"""

using Test
using Dynamic_objectives
using Dynamic_objectives: _detect_local_minima, _compute_dynamic_range,
    _compute_basin_fraction, _detect_deceptive_minima,
    _index_to_point, _evaluate_grid, _default_points_per_dim
using LinearAlgebra: norm

@testset "Grid Scoring" begin

    # ========================================================================
    # _default_points_per_dim
    # ========================================================================

    @testset "default points per dim" begin
        @test _default_points_per_dim(1) == 20
        @test _default_points_per_dim(2) == 20
        @test _default_points_per_dim(3) == 12
        @test _default_points_per_dim(4) == 8
        @test _default_points_per_dim(5) == 6
        @test _default_points_per_dim(10) == 6
    end

    # ========================================================================
    # _detect_local_minima — 1D
    # ========================================================================

    @testset "local minima detection — 1D" begin
        @testset "single interior minimum" begin
            # V-shaped: [5, 3, 1, 3, 5]
            grid = Float64[5, 3, 1, 3, 5]
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test minima[1] == CartesianIndex(3)
        end

        @testset "two minima" begin
            # Two valleys: [5, 1, 5, 2, 5]
            grid = Float64[5, 1, 5, 2, 5]
            minima = _detect_local_minima(grid)
            @test length(minima) == 2
            vals = sort([grid[m] for m in minima])
            @test vals == [1.0, 2.0]
        end

        @testset "flat region — no minima" begin
            # All equal values: no point is strictly less than neighbors
            grid = Float64[3, 3, 3, 3, 3]
            minima = _detect_local_minima(grid)
            @test length(minima) == 0
        end

        @testset "monotone decreasing — boundary minimum" begin
            # [5, 4, 3, 2, 1] — boundary point has only one neighbor
            grid = Float64[5, 4, 3, 2, 1]
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test minima[1] == CartesianIndex(5)
        end

        @testset "Inf values are not minima" begin
            grid = Float64[5, Inf, 1, Inf, 5]
            minima = _detect_local_minima(grid)
            # Inf at indices 2,4 are skipped. Indices 1,5 (val 5) have Inf neighbors
            # which don't block minimum status, so they are minima too.
            # Index 3 (val 1) is the global minimum.
            @test length(minima) == 3
            @test any(grid[m] == 1.0 for m in minima)
            @test !any(grid[m] == Inf for m in minima)
        end

        @testset "Inf neighbors don't block minima" begin
            # [Inf, 3, Inf] — 3 has only Inf neighbors, should be a minimum
            grid = Float64[Inf, 3, Inf]
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test grid[minima[1]] == 3.0
        end
    end

    # ========================================================================
    # _detect_local_minima — 2D
    # ========================================================================

    @testset "local minima detection — 2D" begin
        @testset "single central minimum" begin
            # Bowl: minimum at center
            grid = [
                5.0 4.0 5.0;
                4.0 1.0 4.0;
                5.0 4.0 5.0
            ]
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test minima[1] == CartesianIndex(2, 2)
        end

        @testset "two minima in 2D" begin
            # Two bowls
            grid = [
                5.0 4.0 5.0 4.0 5.0;
                4.0 1.0 4.0 2.0 4.0;
                5.0 4.0 5.0 4.0 5.0
            ]
            minima = _detect_local_minima(grid)
            @test length(minima) == 2
            vals = sort([grid[m] for m in minima])
            @test vals == [1.0, 2.0]
        end

        @testset "corner minimum" begin
            # Minimum at corner (1,1) — only 2 face-adjacent neighbors
            grid = [
                1.0 3.0 5.0;
                3.0 4.0 5.0;
                5.0 5.0 5.0
            ]
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test minima[1] == CartesianIndex(1, 1)
        end

        @testset "face-adjacent only (not diagonal)" begin
            # Center = 3, face-adjacent = 4, diagonal = 1
            # Center should be a minimum because only face-adjacent neighbors matter
            grid = [
                1.0 4.0 1.0;
                4.0 3.0 4.0;
                1.0 4.0 1.0
            ]
            minima = _detect_local_minima(grid)
            # Center at (2,2) is a local min (face neighbors all 4 > 3)
            # Corners at (1,1),(1,3),(3,1),(3,3) have value 1 with face neighbors 4
            @test length(minima) == 5
        end

        @testset "saddle point is not a minimum" begin
            # Center = 3 with lower values on two sides
            grid = [
                5.0 2.0 5.0;
                5.0 3.0 5.0;
                5.0 2.0 5.0
            ]
            minima = _detect_local_minima(grid)
            # (2,2) has neighbor (1,2)=2 < 3, so NOT a minimum
            @test !(CartesianIndex(2, 2) in minima)
        end
    end

    # ========================================================================
    # _detect_local_minima — 3D
    # ========================================================================

    @testset "local minima detection — 3D" begin
        @testset "single minimum in 3x3x3" begin
            grid = fill(5.0, 3, 3, 3)
            grid[2, 2, 2] = 1.0
            minima = _detect_local_minima(grid)
            @test length(minima) == 1
            @test minima[1] == CartesianIndex(2, 2, 2)
        end
    end

    # ========================================================================
    # _compute_dynamic_range
    # ========================================================================

    @testset "dynamic range" begin
        @testset "known range" begin
            grid = Float64[1.0, 10.0, 100.0, 1000.0]
            @test _compute_dynamic_range(grid) ≈ 3.0  # log10(1000/1)
        end

        @testset "all equal" begin
            grid = Float64[5.0, 5.0, 5.0]
            @test _compute_dynamic_range(grid) == 0.0
        end

        @testset "with Inf values" begin
            grid = Float64[1.0, Inf, 100.0, Inf]
            @test _compute_dynamic_range(grid) ≈ 2.0  # log10(100/1)
        end

        @testset "with zeros and negatives" begin
            # Only positive finite values count
            grid = Float64[-1.0, 0.0, 1.0, 100.0]
            @test _compute_dynamic_range(grid) ≈ 2.0  # log10(100/1)
        end

        @testset "single positive value" begin
            grid = Float64[0.0, 5.0, -1.0]
            @test _compute_dynamic_range(grid) == 0.0
        end
    end

    # ========================================================================
    # _index_to_point
    # ========================================================================

    @testset "index to point" begin
        axes = [
            [0.0, 0.5, 1.0],  # dim 1
            [10.0, 20.0, 30.0],  # dim 2
        ]
        @test _index_to_point(CartesianIndex(1, 1), axes) == [0.0, 10.0]
        @test _index_to_point(CartesianIndex(2, 3), axes) == [0.5, 30.0]
        @test _index_to_point(CartesianIndex(3, 2), axes) == [1.0, 20.0]
    end

    # ========================================================================
    # _compute_basin_fraction
    # ========================================================================

    @testset "basin fraction" begin
        @testset "single minimum — all belong to it" begin
            grid = [
                5.0 4.0 5.0;
                4.0 1.0 4.0;
                5.0 4.0 5.0
            ]
            axes = [[0.0, 1.0, 2.0], [0.0, 1.0, 2.0]]
            minima_idx = [CartesianIndex(2, 2)]
            minima_pts = [[1.0, 1.0]]
            p_true = [1.0, 1.0]

            frac = _compute_basin_fraction(grid, axes, minima_idx, minima_pts, p_true)
            @test frac == 1.0
        end

        @testset "two minima — split basin" begin
            # Left half is minimum at (1,1), right half at (1,5)
            # 5x1 grid: [1, 3, 5, 3, 1]
            grid = Float64[1, 3, 5, 3, 1]
            axes = [[0.0, 1.0, 2.0, 3.0, 4.0]]
            minima_idx = [CartesianIndex(1), CartesianIndex(5)]
            minima_pts = [[0.0], [4.0]]
            p_true = [0.0]  # closest to first minimum

            frac = _compute_basin_fraction(grid, axes, minima_idx, minima_pts, p_true)
            # Points 1,2 are closer to [0.0], point 3 equidistant (argmin picks first), points 4,5 closer to [4.0]
            # So basin of first minimum: points at x=0,1,2 → 3/5 = 0.6
            @test frac ≈ 0.6
        end

        @testset "no minima — returns 0" begin
            grid = Float64[1.0, 2.0, 3.0]
            axes = [[0.0, 1.0, 2.0]]
            frac = _compute_basin_fraction(grid, axes, CartesianIndex[], Vector{Float64}[], [1.0])
            @test frac == 0.0
        end
    end

    # ========================================================================
    # _detect_deceptive_minima
    # ========================================================================

    @testset "deceptive minima" begin
        bounds = [(0.0, 10.0), (0.0, 10.0)]
        # Domain diagonal = sqrt(100 + 100) ≈ 14.14

        @testset "far minimum with similar value is deceptive" begin
            p_true = [1.0, 1.0]
            minima_pts = [[1.0, 1.0], [9.0, 9.0]]
            minima_vals = [0.5, 0.8]  # both low, within 2x

            deceptive = _detect_deceptive_minima(
                minima_pts, minima_vals, p_true, bounds, 2.0, 0.1,
            )
            # [9,9] is far from p_true (dist ≈ 11.3, normalized ≈ 0.80 > 0.1)
            # and has val 0.8 < 2.0 * 0.5 = 1.0
            @test length(deceptive) == 1
            @test deceptive[1] == [9.0, 9.0]
        end

        @testset "nearby minimum is not deceptive" begin
            p_true = [5.0, 5.0]
            minima_pts = [[5.0, 5.0], [5.5, 5.5]]
            minima_vals = [0.5, 0.6]

            deceptive = _detect_deceptive_minima(
                minima_pts, minima_vals, p_true, bounds, 2.0, 0.1,
            )
            # [5.5, 5.5] is close to p_true (dist ≈ 0.71, normalized ≈ 0.05 < 0.1)
            @test length(deceptive) == 0
        end

        @testset "high-value minimum is not deceptive" begin
            p_true = [1.0, 1.0]
            minima_pts = [[1.0, 1.0], [9.0, 9.0]]
            minima_vals = [0.5, 50.0]  # 50 >> 2 * 0.5

            deceptive = _detect_deceptive_minima(
                minima_pts, minima_vals, p_true, bounds, 2.0, 0.1,
            )
            @test length(deceptive) == 0
        end

        @testset "negative objective values (log-scale)" begin
            p_true = [1.0, 1.0]
            minima_pts = [[1.0, 1.0], [9.0, 9.0]]
            # log-scale: more negative = better
            # global_min = -100, threshold = 2.0
            # "low" means val < -100 / 2.0 = -50
            minima_vals = [-100.0, -80.0]

            deceptive = _detect_deceptive_minima(
                minima_pts, minima_vals, p_true, bounds, 2.0, 0.1,
            )
            # -80 < -50? Yes. Far from p_true? Yes.
            @test length(deceptive) == 1
        end

        @testset "empty minima" begin
            deceptive = _detect_deceptive_minima(
                Vector{Float64}[], Float64[], [1.0, 1.0], bounds, 2.0, 0.1,
            )
            @test isempty(deceptive)
        end
    end

    # ========================================================================
    # _evaluate_grid
    # ========================================================================

    @testset "evaluate grid" begin
        @testset "1D quadratic" begin
            obj(p) = (p[1] - 2.0)^2
            axes = [[0.0, 1.0, 2.0, 3.0, 4.0]]
            vals = _evaluate_grid(obj, axes)
            @test size(vals) == (5,)
            @test vals[1] ≈ 4.0
            @test vals[3] ≈ 0.0
            @test vals[5] ≈ 4.0
        end

        @testset "2D separable" begin
            obj(p) = p[1]^2 + p[2]^2
            axes = [[0.0, 1.0], [0.0, 1.0]]
            vals = _evaluate_grid(obj, axes)
            @test size(vals) == (2, 2)
            @test vals[1, 1] ≈ 0.0
            @test vals[2, 2] ≈ 2.0
        end
    end

    # ========================================================================
    # interestingness_score
    # ========================================================================

    @testset "interestingness score" begin
        @testset "boring landscape (single minimum, full basin)" begin
            gs = GridScoreResult(
                zeros(5, 5),                     # grid_values
                [[0.0, 1.0, 2.0, 3.0, 4.0],
                 [0.0, 1.0, 2.0, 3.0, 4.0]],   # grid_axes
                [5, 5],                           # grid_size
                1,                                # n_local_minima
                [CartesianIndex(3, 3)],           # local_minima_indices
                [0.0],                            # local_minima_values
                [[2.0, 2.0]],                     # local_minima_points
                0.5,                              # dynamic_range
                1.0,                              # basin_fraction (100% in p_true basin)
                0,                                # n_deceptive
                Vector{Float64}[],                # deceptive_points
                25,                               # n_finite
                25,                               # n_total
                10.0,                             # evaluation_time_ms
                zeros(2),                         # directional_variation
                0.0,                              # curvature_score (flat)
                0.0,                              # basin_depth_ratio (no competitor)
                0.0,                              # conditioning
                0.0,                              # ruggedness
                0.0,                              # plateau_fraction
            )
            score = interestingness_score(gs)
            # Only dyn contributes: 0.10 * (0.5/6) ≈ 0.008
            @test score < 0.1
        end

        @testset "interesting landscape (many minima, small basin, deceptive)" begin
            gs = GridScoreResult(
                zeros(5, 5),
                [[0.0, 1.0, 2.0, 3.0, 4.0],
                 [0.0, 1.0, 2.0, 3.0, 4.0]],
                [5, 5],
                15,                               # 15 local minima
                CartesianIndex{2}[],              # (details not needed for score)
                Float64[],
                Vector{Float64}[],
                8.0,                              # 8 orders of magnitude
                0.05,                             # only 5% basin for p_true
                5,                                # 5 deceptive minima
                Vector{Float64}[],
                25, 25, 10.0,
                ones(2),                          # directional_variation
                1.0,                              # curvature_score (max)
                1.0,                              # basin_depth_ratio (equally deep competitor)
                4.0,                              # conditioning (saturates at 4)
                1.0,                              # ruggedness (saturates at 1)
                0.7,                              # plateau_fraction (saturates at 0.7)
            )
            score = interestingness_score(gs)
            # After rebalance: basin=max(0,1-4*0.05)=0.8 (weight 0.12 ⇒ 0.096)
            # All other components saturate ⇒ score ≈ 0.976
            @test score > 0.9
        end

        @testset "basin narrowness rebalance (slope 4)" begin
            # basin_fraction=0.25 is the cutoff — basin term contributes 0.
            make_gs = (bf) -> GridScoreResult(
                zeros(5, 5), [[0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 2.0, 3.0, 4.0]], [5, 5],
                1, CartesianIndex{2}[], Float64[], Vector{Float64}[],
                0.0,               # dynamic_range (zero so only basin contributes)
                bf, 0, Vector{Float64}[],
                25, 25, 1.0,
                zeros(2), 0.0, 0.0, 0.0, 0.0, 0.0,
            )
            @test interestingness_score(make_gs(0.0)) ≈ 0.12 atol=1e-9
            @test interestingness_score(make_gs(0.1)) ≈ 0.12 * 0.6 atol=1e-9
            @test interestingness_score(make_gs(0.25)) ≈ 0.0 atol=1e-9
            @test interestingness_score(make_gs(0.5)) ≈ 0.0 atol=1e-9
        end

        @testset "basin depth ratio discounted by multimodal context" begin
            # depth_ratio=1.0 against a single competitor is degenerate → 0 credit.
            make_gs = (n_min, dr) -> GridScoreResult(
                zeros(5, 5), [[0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 2.0, 3.0, 4.0]], [5, 5],
                n_min, CartesianIndex{2}[], Float64[], Vector{Float64}[],
                0.0, 1.0, 0, Vector{Float64}[],   # basin_fraction=1 ⇒ basin term=0
                25, 25, 1.0,
                zeros(2), 0.0, dr, 0.0, 0.0, 0.0,
            )
            # n_min=1 ⇒ depth_context=0 regardless of dr
            @test interestingness_score(make_gs(1, 1.0)) ≈ 0.0 atol=1e-9
            # n_min=5 ⇒ depth_context=1.0 ⇒ full depth credit + multi contribution
            @test interestingness_score(make_gs(5, 1.0)) ≈ 0.12 + 0.20 * (4 / 14) atol=1e-9
        end

        @testset "score is bounded [0, 1]" begin
            # Extreme case
            gs = GridScoreResult(
                zeros(1), [[0.0]], [1],
                100, CartesianIndex{1}[], Float64[], Vector{Float64}[],
                100.0, 0.0, 100, Vector{Float64}[],
                1, 1, 1.0,
                ones(1),                          # directional_variation
                1.0, 1.0, 100.0, 100.0, 1.0,     # curv, depth, cond, rugged, plateau (all clamp to 1)
            )
            score = interestingness_score(gs)
            @test 0.0 <= score <= 1.0
        end
    end

    # ========================================================================
    # print_grid_score (smoke test)
    # ========================================================================

    @testset "print_grid_score" begin
        gs = GridScoreResult(
            [1.0 2.0; 3.0 4.0],
            [[0.0, 1.0], [0.0, 1.0]],
            [2, 2],
            1, [CartesianIndex(1, 1)], [1.0], [[0.0, 0.0]],
            0.6, 0.5, 0, Vector{Float64}[],
            4, 4, 5.0,
            [0.5, 0.5],                          # directional_variation
            0.5, 0.5, 2.0, 0.5, 0.3,             # curv, depth, cond, rugged, plateau
        )
        buf = IOBuffer()
        print_grid_score(gs; io=buf)
        output = String(take!(buf))
        @test occursin("Grid-Based Landscape Score", output)
        @test occursin("Local minima:", output)
        @test occursin("Dynamic range:", output)
        @test occursin("Interestingness:", output)
    end

end  # @testset "Grid Scoring"
