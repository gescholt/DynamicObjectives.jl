using Test
using DynamicObjectives
import DynamicObjectives: Tsit5

# End-to-end regression test for grid-based interestingness scoring of the
# top p_true candidates.
#
# The implementation — screen_and_probe → score_top_candidates →
# interestingness_score — already exists in src/screening.jl and
# src/grid_scoring.jl. This suite asserts the three layers compose and
# that interestingness_score stays within its documented [0, 1] range
# on the pipeline output.
#
# Runtime budget: LV 2D fixture, tiny candidate counts, points_per_dim=4
# → ~32 ODE solves total at Tsit5/1e-4.

@testset "Interestingness end-to-end" begin
    model, params, states, outputs = define_lotka_volterra_2D_model_v3()
    ic = [1.0, 1.0]
    bounds = [(0.5, 2.0), (0.5, 2.0)]
    time_interval = [0.0, 5.0]

    screening = screen_and_probe(
        model,
        outputs,
        ic,
        bounds,
        time_interval;
        n_candidates = 20,
        n_probes = 6,
        numpoints_screen = 10,
        numpoints_probe = 15,
        solver = Tsit5(),
        abstol = 1e-4,
        reltol = 1e-4,
        verbose = false,
    )

    @testset "screen_and_probe produced a ranked result" begin
        @test isa(screening, DynamicObjectives.ScreeningResult)
        @test !isempty(screening.sweep.valid)
        # Ranking length bounded by number of candidates that passed probe.
        @test length(screening.ranking.ranked_indices) <= length(screening.probes)
    end

    # Skip the grid-scoring leg if the pre-screen rejected every candidate
    # (possible on a degenerate fixture — don't want a flaky failure here).
    if isempty(screening.ranking.ranked_indices)
        @warn "Screening produced no ranked candidates; skipping grid-scoring leg"
    else
        top_n = min(2, length(screening.ranking.ranked_indices))
        grid_results = DynamicObjectives.score_top_candidates(
            screening,
            model,
            outputs,
            ic,
            bounds,
            time_interval;
            top_n = top_n,
            points_per_dim = 4,
            numpoints = 10,
            solver = Tsit5(),
            abstol = 1e-4,
            reltol = 1e-4,
        )

        @testset "score_top_candidates returns top_n entries" begin
            @test length(grid_results) == top_n
        end

        successes = filter(!isnothing, grid_results)
        @test !isempty(successes)

        @testset "interestingness_score lies in [0, 1]" begin
            for gs in successes
                s = DynamicObjectives.interestingness_score(gs)
                @test isfinite(s)
                @test 0.0 <= s <= 1.0
            end
        end
    end
end
