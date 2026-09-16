using Test
using DynamicObjectives

# Candidate-level parallelism — assert that serial and parallel
# `run_catalogue_experiments` return identical result vectors for a
# deterministic runner. Uses `CatalogueEntry` fixtures so the code path
# actually touches the type signature consumers will hit in practice.

@testset "run_catalogue_experiments" begin

    # Build a handful of distinct entries. The model_fn / bounds values are
    # irrelevant here — the deterministic runner only reads `entry.name`.
    function _make_entry(name::String)
        CatalogueEntry(
            name = name,
            description = "fixture entry $name",
            model_fn = () -> nothing,
            p_true = [0.0],
            ic = [0.0],
            bounds = [(0.0, 1.0)],
        )
    end

    entries = [_make_entry("a"), _make_entry("b"), _make_entry("c"), _make_entry("d")]
    runner = entry -> hash(entry.name)

    @testset "serial path" begin
        results = run_catalogue_experiments(entries, runner)
        @test length(results) == 4
        @test results == [hash("a"), hash("b"), hash("c"), hash("d")]
    end

    @testset "parallel path matches serial" begin
        serial = run_catalogue_experiments(entries, runner; parallel = false)
        parallel = run_catalogue_experiments(entries, runner; parallel = true)
        @test parallel == serial
    end

    @testset "parallel preserves order with varied chunking" begin
        # max_workers kwarg controls chunk granularity; results must be
        # index-aligned regardless of how @spawn schedules the chunks.
        for max_workers in (1, 2, 4, 8)
            results = run_catalogue_experiments(
                entries,
                runner;
                parallel = true,
                max_workers = max_workers,
            )
            @test results == [hash("a"), hash("b"), hash("c"), hash("d")]
        end
    end

    @testset "degenerate inputs" begin
        @test run_catalogue_experiments(CatalogueEntry[], runner) == UInt64[]
        single = [_make_entry("solo")]
        @test run_catalogue_experiments(single, runner; parallel = true) == [hash("solo")]
    end
end
