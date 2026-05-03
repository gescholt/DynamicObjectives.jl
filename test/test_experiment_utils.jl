# Tests for experiment_utils.jl — build_bounds, CandidateResult, print_candidate_summary_table

@testset "Experiment Utilities" begin

    # ── build_bounds (scalar radius) ─────────────────────────────────────

    @testset "build_bounds — scalar radius" begin
        @testset "basic 2D" begin
            bounds = build_bounds([1.0, 2.0], 0.5)
            @test length(bounds) == 2
            @test bounds[1] == (0.5, 1.5)
            @test bounds[2] == (1.5, 2.5)
        end

        @testset "1D" begin
            bounds = build_bounds([3.0], 1.0)
            @test bounds == [(2.0, 4.0)]
        end

        @testset "4D" begin
            bounds = build_bounds([0.1, 0.2, 0.3, 0.4], 0.01)
            @test length(bounds) == 4
            for (i, (lo, hi)) in enumerate(bounds)
                @test hi - lo ≈ 0.02
                @test (lo + hi) / 2 ≈ 0.1 * i
            end
        end

        @testset "zero radius" begin
            bounds = build_bounds([1.0, 2.0], 0.0)
            @test bounds[1] == (1.0, 1.0)
            @test bounds[2] == (2.0, 2.0)
        end

        @testset "negative center" begin
            bounds = build_bounds([-1.0, -2.0], 0.5)
            @test bounds[1] == (-1.5, -0.5)
            @test bounds[2] == (-2.5, -1.5)
        end
    end

    # ── build_bounds (per-dimension radii) ───────────────────────────────

    @testset "build_bounds — per-dimension radii" begin
        @testset "basic anisotropic" begin
            bounds = build_bounds([1.0, 2.0], [0.1, 0.5])
            @test bounds[1] == (0.9, 1.1)
            @test bounds[2] == (1.5, 2.5)
        end

        @testset "length mismatch errors" begin
            @test_throws ErrorException build_bounds([1.0, 2.0], [0.1])
            @test_throws ErrorException build_bounds([1.0], [0.1, 0.2])
        end
    end

    # ── CandidateResult ──────────────────────────────────────────────────

    @testset "CandidateResult construction" begin
        r = CandidateResult("test", [1.0, 2.0], 6, 3, 0.001, 0.01, 0.05, 0.1, 0.8, 0.95)
        @test r.name == "test"
        @test r.p_true == [1.0, 2.0]
        @test r.degree == 6
        @test r.n_cps == 3
        @test r.raw_obj == 0.001
        @test r.raw_recovery == 0.01
        @test r.l2_error == 0.05
        @test r.rel_l2 == 0.1
        @test r.capture_1pct == 0.8
        @test r.capture_5pct == 0.95
    end

    # ── print_candidate_summary_table ────────────────────────────────────

    @testset "print_candidate_summary_table" begin
        results = [
            CandidateResult("Model_A", [1.0, 0.5], 4, 2, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8),
            CandidateResult("Model_A", [1.0, 0.5], 6, 3, 0.005, 0.02, 0.08, 0.15, 0.7, 0.9),
            CandidateResult("Model_B", [0.3, 0.7], 4, 1, 0.1, 0.2, 0.3, 0.4, 0.0, 0.5),
        ]
        entries = [(name = "Model_A",), (name = "Model_B",)]

        # Capture output — just verify it doesn't error
        original_stdout = stdout
        redirect_stdout(devnull)
        try
            best = print_candidate_summary_table(
                results,
                entries;
                title = "Test",
                radius = 0.1,
                degree_range = 4:2:6,
            )
            @test length(best) == 2
            @test best[1].name == "Model_A"
            @test best[1].degree == 6  # lower raw_recovery
            @test best[2].name == "Model_B"
        finally
            redirect_stdout(original_stdout)
        end
    end

    @testset "print_candidate_summary_table — NaN values" begin
        results = [CandidateResult("NaN_model", [1.0], 4, 0, NaN, NaN, 0.0, NaN, NaN, NaN)]
        entries = [(name = "NaN_model",)]

        original_stdout = stdout
        redirect_stdout(devnull)
        try
            best = print_candidate_summary_table(
                results,
                entries;
                title = "NaN Test",
                radius = 0.1,
                degree_range = 4:2:4,
            )
            # All NaN recovery — no valid results
            @test isempty(best)
        finally
            redirect_stdout(original_stdout)
        end
    end
end
